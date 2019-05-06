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


class BWA:
    def __init__(self, threads=1):
        self.threads = threads
        self.reference = ""

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

        if self.alignment.returncode != 0:
            if b"fail to locate the index" in self.alignment.stderr:
                raise IndexError("Reference must be indexed by BWA before alignment.")
            else:
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

            left_idxs, right_idxs = self.calculate_probe_boundaries_for_entry(variant)
            left_flank = gene.sequence[slice(*left_idxs)]
            right_flank = gene.sequence[slice(*right_idxs)]

            probe = self.create_probe_for_variant(variant)
            probe.set_sequence(left_flank + get_variant_sequence(variant) + right_flank)
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

        probe.set_name(
            probe.name
            + f"_entry{self._entry_number}_CONF{get_genotype_confidence(variant)}"
        )

        return probe

    def calculate_probe_boundaries_for_entry(self, entry: pysam.VariantRecord) -> tuple:
        variant_len = get_variant_length(entry)
        delta_len = self._min_probe_length - variant_len
        left = [entry.start, entry.start]
        right = [entry.stop, entry.stop]
        variant_shorter_than_min_probe_len = delta_len > 0

        if variant_shorter_than_min_probe_len:
            flank_width = delta_len // 2
            left[0] = max(0, entry.start - flank_width)
            right[1] = entry.stop + flank_width

        return left, right


def is_invalid_vcf_entry(entry: pysam.VariantRecord) -> bool:
    genotype = get_genotype(entry)

    return genotype is None


def get_genotype_confidence(variant: pysam.VariantRecord) -> float:
    return float(variant.samples["sample"].get("GT_CONF", 0))


def get_genotype(variant: pysam.VariantRecord) -> int:
    return variant.samples["sample"]["GT"][0]


def get_variant_sequence(variant: pysam.VariantRecord) -> str:
    genotype = get_genotype(variant)

    if genotype is None:
        return variant.ref
    else:
        return variant.alleles[genotype]


def get_variant_length(variant: pysam.VariantRecord) -> int:
    return len(get_variant_sequence(variant))


def is_mapping_invalid(record: pysam.AlignedSegment) -> bool:
    return any([record.is_unmapped, record.is_secondary, record.is_supplementary])


def is_snp_called_correctly(record: pysam.AlignedSegment) -> bool:
    expected_base = record.reference_name[-1]
    snp_idx = REF_PANEL_FLANK_WIDTH - record.reference_start
    actual_base = record.query_alignment_sequence[snp_idx]
    return expected_base == actual_base


def map_probes_to_panel(probes: str, reference_panel: Path, threads=1) -> dict:
    bwa = BWA(threads)
    bwa.index(reference_panel)
    header, sam = bwa.align(probes)

    results = {"snps_called_correctly": [], "mismatches": [], "ids": []}

    valid_pandora_calls = 0
    sites_seen = set()

    for record in sam:
        if is_mapping_invalid(record):
            continue
        elif not (record.reference_start <= 100 < record.reference_end):
            continue
        valid_pandora_calls += 1

        if record.reference_name not in sites_seen:
            sites_seen.add(record.reference_name)

        results["snps_called_correctly"].append(is_snp_called_correctly(record))
        results["mismatches"].append(record.get_tag("NM"))
        results["ids"].append(record.query_name)

    results["total_pandora_calls"] = len(sam)
    results["pandora_calls_crossing_ref_site"] = valid_pandora_calls
    results["reference_sites_called"] = len(sites_seen)

    return results


def write_results(results: dict, output: Path):
    with output.open("w") as fh:
        json.dump(results, fh, indent=4)


def main():
    reference_panel = Path(snakemake.output.reference_panel)
    reference = Reference(
        Path(snakemake.input.reference_vcf), Path(snakemake.input.reference_seq)
    )
    reference.make_panel(reference_panel)

    query = Query(
        Path(snakemake.input.query_vcf), Path(snakemake.input.pandora_consensus)
    )
    query_probes = query.make_probes()

    results = map_probes_to_panel(query_probes, reference_panel, snakemake.threads)
    results["total_reference_sites"] = snakemake.wildcards.num_snps

    write_results(results, Path(snakemake.output.results))


if __name__ == "__main__":
    main()
