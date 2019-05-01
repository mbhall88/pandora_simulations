from pathlib import Path
from contextlib import ExitStack
import subprocess
import pysam
import json

REF_PANEL_FLANK_WIDTH = 100
MIN_QUERY_PROBE_LEN = 70


class Sequence(str):
    def get_left_flank(self, pos, flank_width):
        start = pos - flank_width

        if start < 0:
            return self[start:] + self[:pos]
        else:
            return self[start:pos]

    def get_right_flank(self, pos, flank_width):
        end = pos + flank_width + 1

        if end > len(self):
            return self[pos + 1 :] + self[: end - len(self)]
        else:
            return self[pos + 1 : end]

    def get_probe(self, pos, flank_width):
        return (
            self.get_left_flank(pos[0], flank_width)
            + self[slice(*pos)]
            + self.get_right_flank(pos[1] - 1, flank_width)
        )


# =================================================
# test Sequence methods
s = Sequence("0123456789")
assert s.get_left_flank(0, 4) == "6789"
assert s.get_left_flank(2, 4) == "8901"
assert s.get_left_flank(9, 4) == "5678"
assert s.get_right_flank(3, 4) == "4567"
assert s.get_right_flank(8, 4) == "9012"
assert s.get_right_flank(0, 4) == "1234"
assert s.get_probe((5, 6), 4) == "123456789"
assert s.get_probe((2, 3), 4) == "890123456"
assert s.get_probe((7, 8), 4) == "345678901"
assert s.get_probe((0, 2), 4) == "6789012345"
assert s.get_probe((7, 9), 4) == "3456789012"
# =================================================


class BWA:
    def __init__(self, threads: int):
        self.threads = threads

    def index(self, reference):
        self.reference = reference

        completed_process = subprocess.run(
            ["bwa", "index", self.reference],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        completed_process.check_returncode()

    def align(self, query: str) -> tuple:
        options = self.get_options()
        self.alignment = subprocess.run(
            ["bwa", "mem", *options, str(self.reference), "-"],
            input=query.encode(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.alignment.check_returncode()

        return self._get_samfile()

    def get_options(self):
        options = []
        options.extend(["-t", str(self.threads)])

        return options

    def _get_samfile(self):
        header = ""
        sam_lines = []
        for line in self.alignment.stdout.decode().split("\n"):
            if line.startswith("@"):
                header += line + "\n"
            else:
                sam_lines.append(line)

        header = pysam.AlignmentHeader.from_text(header)

        return (
            header,
            [pysam.AlignedSegment.fromstring(sam, header) for sam in sam_lines if sam],
        )


class Reference:
    def __init__(self, vcf: Path, sequence: Path):
        self.vcf = vcf
        self.sequence = sequence

    def make_panel(self, output_path: Path, flank_width=REF_PANEL_FLANK_WIDTH):
        with ExitStack() as stack:
            vcf = stack.enter_context(pysam.VariantFile(self.vcf))
            ref = stack.enter_context(pysam.FastxFile(str(self.sequence)))
            panel = stack.enter_context(output_path.open("w"))

            for i, record in enumerate(ref):
                ref_seq = Sequence(record.sequence)

            assert i == 0

            for entry in vcf:
                assert (
                    str(ref_seq[entry.start]) == entry.alts[0]
                ), f"Pos: {entry.pos}\nGot {ref_seq[entry.start]}. Expecting {entry.alts[0]}"

                probe_name = (
                    f">{entry.ref}{entry.pos}{entry.alts[0]} flank_width={flank_width}"
                )
                probe = ref_seq.get_probe((entry.start, entry.stop), flank_width)

                panel.write(f"{probe_name}\n{probe}\n")


def is_invalid_vcf_entry(entry: pysam.VariantRecord) -> bool:
    return any(gt is None for gt in entry.samples["sample"]["GT"])


class Query:
    def __init__(self, vcf: Path, genes: Path):
        self.vcf = Path(
            pysam.tabix_index(str(vcf), preset="vcf", keep_original=True, force=True)
        )
        self.genes = genes
        self._probe_names = set()
        self._min_probe_length = MIN_QUERY_PROBE_LEN
        self._entry_number = 0

    def make_probes(self) -> str:
        query_probes = ""

        with ExitStack() as stack:
            vcf = stack.enter_context(pysam.VariantFile(self.vcf))
            genes_fasta = stack.enter_context(pysam.FastxFile(str(self.genes)))

            for gene in genes_fasta:
                try:
                    entries = vcf.fetch(contig=gene.name)
                except ValueError as error:
                    if str(error).startswith("invalid contig"):
                        continue
                    else:
                        raise

                probes_for_gene = self.create_probes_for_gene_variants(gene, entries)
                query_probes += probes_for_gene

        return query_probes

    def create_probes_for_gene_variants(
        self, gene: pysam.FastxRecord, variants: pysam.tabix_iterator
    ) -> str:
        probes = ""

        for variant in variants:
            if is_invalid_vcf_entry(variant):
                continue

            start, end = self.calculate_probe_boundaries_for_entry(variant)
            probe = self.create_probe_for_variant(variant)
            probe.set_sequence(gene.sequence[start:end])
            probes += str(probe) + "\n"

        return probes

    def create_probe_for_variant(
        self, variant: pysam.VariantRecord
    ) -> pysam.FastxRecord:
        probe = pysam.FastxRecord()
        probe.set_name(f"{variant.chrom}_pos{variant.pos}")

        if probe.name in self._probe_names:
            self._entry_number += 1
        else:
            self._probe_names.add(probe.name)
            self._entry_number = 0

        probe.set_name(probe.name + f"_entry{self._entry_number}")

        return probe

    def calculate_probe_boundaries_for_entry(self, entry: pysam.VariantRecord) -> tuple:
        start = entry.start
        end = entry.stop

        if entry.rlen < self._min_probe_length:
            diff = self._min_probe_length - entry.rlen
            start = max(0, start - diff // 2)
            end += diff // 2

        return start, end


def is_mapping_invalid(record: pysam.AlignedSegment) -> bool:
    return any([record.is_unmapped, record.is_secondary, record.is_supplementary])


def is_snp_called_correctly(record: pysam.AlignedSegment) -> bool:
    expected_base = record.reference_name[-1]
    snp_idx = 100 - record.reference_start
    actual_base = record.query_alignment_sequence[snp_idx]
    return expected_base == actual_base


def map_probes_to_panel(probes: str, reference_panel: Path, threads: int) -> dict:
    bwa = BWA(threads)
    bwa.index(reference_panel)
    header, sam = bwa.align(probes)

    results = {"snps_called_correctly": [], "mismatches": []}

    valid_pandora_calls = 0
    sites_seen = set()

    for record in sam:
        if is_mapping_invalid(record):
            continue
        elif not record.reference_start <= 100 < record.reference_end:
            continue
        valid_pandora_calls += 1

        if record.reference_name not in sites_seen:
            sites_seen.add(record.reference_name)

        results["snps_called_correctly"].append(is_snp_called_correctly(record))
        results["mismatches"].append(record.get_tag("NM"))

    results["total_pandora_calls"] = len(sam)
    results["pandora_calls_crossing_ref_site"] = valid_pandora_calls
    results["reference_sites_called"] = len(sites_seen)

    return results


def write_results(results: dict, output: Path):
    with output.open("w") as fh:
        json.dump(results, fh, indent=4)


def main():
    reference_vcf = Path("combined_random_paths_mutated.vcf")
    reference_seq = Path("combined_random_paths_mutated_1.fasta")
    reference_panel = Path("panel.fa")
    query_vcf = Path("pandora_genotyped.vcf")
    pandora_consensus = Path("pandora.consensus.fq.gz")
    threads = 1
    num_snps = 250
    output = Path("output.json")

    reference = Reference(reference_vcf, reference_seq)

    reference.make_panel(reference_panel)

    query = Query(query_vcf, pandora_consensus)

    query_probes = query.make_probes()

    results = map_probes_to_panel(query_probes, reference_panel, threads)
    results["total_reference_sites"] = num_snps

    write_results(results, output)


if __name__ == "__main__":
    main()
