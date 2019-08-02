from pathlib import Path
from contextlib import ExitStack
import subprocess
import pysam
import sys
import json
from typing import List, Tuple, Dict

REF_PANEL_FLANK_WIDTH = 25
QUERY_PROBE_FLANK_WIDTH = 50


class OverlappingRecordsError(Exception):
    pass


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

    def index(self, reference: str):
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
    def __init__(
        self, vcf: Path, genes: Path, min_probe_length: int = QUERY_PROBE_FLANK_WIDTH
    ):
        self.vcf = Path(
            pysam.tabix_index(str(vcf), preset="vcf", keep_original=True, force=True)
        )
        self.genes = genes
        self._probe_names = set()
        self._min_probe_length = min_probe_length
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
                        raise error

                probes_for_gene = self._create_probes_for_gene_variants(gene, entries)
                query_probes += probes_for_gene

        return query_probes

    def _create_probes_for_gene_variants(
        self, gene: pysam.FastxRecord, variants: pysam.tabix_iterator
    ) -> str:
        """Note: An assumption is made with this function that the variants you pass in
        are from the gene passed with them."""
        probes = ""
        variants = [entry for entry in variants if not is_invalid_vcf_entry(entry)]

        # if any_entries_overlap(variants):
        #     raise OverlappingRecordsError(
        #         f"Found overlapping variant records in {str(self.vcf)}"
        #     )

        # intervals = merge_overlap_intervals(
        #     [self.calculate_probe_boundaries_for_entry(variant) for variant in variants]
        # )
        intervals = [
            self.calculate_probe_boundaries_for_entry(variant) for variant in variants
        ]
        intervals_to_variants = assign_variants_to_intervals(variants, intervals)

        intervals_to_probes = dict()

        for variant in variants:
            interval = self.calculate_probe_boundaries_for_entry(variant)
            if interval in intervals_to_probes and float(
                intervals_to_probes[interval].name.split("=")[-1]
            ) > get_genotype_confidence(variant):
                continue

            mutated_consensus = ""
            consensus = gene.sequence[slice(*interval)]
            last_idx = 0

            start_idx_of_variant_on_consensus = variant.start - interval[0]
            mutated_consensus += consensus[last_idx:start_idx_of_variant_on_consensus]
            mutated_consensus += get_variant_sequence(variant)
            last_idx = start_idx_of_variant_on_consensus + variant.rlen
            mutated_consensus += consensus[last_idx:]
            probe = pysam.FastxRecord()
            probe.set_name(
                f"{variant.chrom}_POS={variant.pos}_interval={interval}_GT_CONF={get_genotype_confidence(variant)}".replace(" ", "")
            )
            probe.set_sequence(mutated_consensus)
            intervals_to_probes[interval] = probe

        for probe in intervals_to_probes.values():
            probes += str(probe) + "\n"

        return probes

    def calculate_probe_boundaries_for_entry(
        self, entry: pysam.VariantRecord
    ) -> Tuple[int, int]:
        variant_len = get_variant_length(entry)
        delta_len = self._min_probe_length - variant_len
        left = [entry.start, entry.start]
        right = [entry.stop, entry.stop]
        variant_shorter_than_min_probe_len = delta_len > 0

        if variant_shorter_than_min_probe_len:
            flank_width = delta_len // 2
            left[0] = max(0, entry.start - flank_width)
            right[1] = entry.stop + flank_width

        return left[0], right[-1]


def merge_overlap_intervals(intervals: List[List[int]]) -> List[Tuple[int, int]]:
    """Checks consecutive intervals and if they overlap it merges them into a
    single interval.
    Args:
        intervals: A list of intervals where each interval is a List with two
        elements corresponding to the start and end of the interval
        respectively.
    Returns:
        A new intervals list where any intervals that overlapped have been
        merged into a single interval.
    Example:
        >>> intervals = [[1, 4], [3, 7], [10, 14]]
        >>> merge_overlap_intervals(intervals)
        [(1, 7), (10, 14)]
    """
    merged_intervals = []
    cached_interval = None

    for interval in intervals:
        if cached_interval is None:
            cached_interval = interval
            continue

        if outside_interval(cached_interval, interval):
            merged_intervals.append(tuple(cached_interval))
            cached_interval = interval
        else:
            cached_interval = extend_interval(cached_interval, interval)

    if cached_interval is not None:
        merged_intervals.append(tuple(cached_interval))

    return merged_intervals


def outside_interval(first_interval: List[int], second_interval: List[int]) -> bool:
    """Determines whether two intervals overlap.
    Args:
        first_interval: The interval with the lower start index.
        second_interval: The interval with the higher start index.
    Returns:
        Whether the start of the second interval is less than the end of the
        first interval. i.e do they overlap?
    Notes:
        If the end index of the first interval is equal to the start of the
        second interval than they are deemed to NOT be overlapping.
    Example:
        >>> first_interval = [0, 4]
        >>> second_interval = [3, 7]
        >>> outside_interval(first_interval, second_interval)
        False
    """
    return second_interval[0] >= first_interval[1]


def extend_interval(interval_to_extend: List[int], interval: List[int]) -> List[int]:
    """Extends an interval to encompass another.
    Args:
        interval_to_extend: The interval to extend.
        interval: The interval to extend by.
    Returns:
        A new interval with the same start as interval_to_extend and the same
        end as interval.
    """
    interval_to_extend[1] = interval[1]

    return interval_to_extend


def find_index_in_intervals(intervals: List[Tuple[int, int]], query: int) -> int:
    """Return the index of the interval that the query lies within"""
    for i, (start, end) in enumerate(intervals):
        if start <= query <= end:
            return i
    return -1


def is_invalid_vcf_entry(entry: pysam.VariantRecord) -> bool:
    genotype = get_genotype(entry)

    return genotype is None


def get_genotype_confidence(variant: pysam.VariantRecord) -> float:
    return float(variant.samples["sample"].get("GT_CONF", 0))


def get_genotype(variant: pysam.VariantRecord) -> int:
    samples = list(variant.samples.keys())
    assert len(samples) == 1
    return variant.samples[samples[0]]["GT"][0]


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
    for query_pos, ref_pos, ref_base in record.get_aligned_pairs(with_seq=True):
        if ref_pos == 100:
            if ref_base.islower():
                return False
            else:
                return True


def map_probes_to_panel(probes: str, reference_panel: Path, threads=1) -> dict:
    bwa = BWA(threads)
    bwa.index(str(reference_panel))
    header, sam = bwa.align(probes)

    results = {"snps_called_correctly": [], "mismatches": [], "ids": [], "ref_ids": []}

    valid_pandora_calls = 0
    sites_seen = set()

    for record in sam:
        if is_mapping_invalid(record):
            continue
        elif not (
            record.reference_start <= REF_PANEL_FLANK_WIDTH < record.reference_end
        ):
            continue
        valid_pandora_calls += 1

        if record.reference_name not in sites_seen:
            sites_seen.add(record.reference_name)

        results["snps_called_correctly"].append(is_snp_called_correctly(record))
        results["mismatches"].append(record.get_tag("NM"))
        results["ids"].append(record.query_name)
        results["ref_ids"].append(record.reference_name)

    results["total_pandora_calls"] = len(sam)
    results["pandora_calls_crossing_ref_site"] = valid_pandora_calls
    results["reference_sites_called"] = len(sites_seen)

    return results


def map_panel_to_probes(
    panel: Path, probes: Path, threads: int = 1
) -> List[pysam.AlignedSegment]:
    bwa = BWA(threads)
    bwa.index(str(probes))
    header, sam = bwa.align(panel.read_text())

    return [record for record in sam if not is_mapping_invalid(record)]


def get_results_from_alignment(alignment: List[pysam.AlignedSegment]) -> dict:
    probes_mapped = dict()

    for entry in alignment:
        if entry.is_unmapped:
            continue
        if entry.query_name in probes_mapped:
            raise KeyError(f"{entry.query_name} has already been mapped.")
        probes_mapped[entry.query_name] = {
            "correct": record_contains_expected_snp(entry),
            "id": entry.reference_name,
        }

    return probes_mapped


def write_results(results: dict, output: Path):
    with output.open("w") as fh:
        json.dump(results, fh, indent=4)


def record_contains_expected_snp(record: pysam.AlignedSegment) -> bool:
    expected_base = record.query_name[-1]

    for query_pos, ref_pos, ref_base in record.get_aligned_pairs(with_seq=True):
        if query_pos == REF_PANEL_FLANK_WIDTH:
            return expected_base == ref_base

    return False


def evaluate_candidates(candidate_files: List[Path], panel: str, threads: int) -> dict:
    probes_mapped_to_candidate = dict()
    slices_containing_mutation = set()

    for candidate in candidate_files:
        bwa = BWA(threads)
        bwa.index(str(candidate))
        header, sam = bwa.align(panel)

        for entry in sam:
            if entry.is_unmapped or (
                entry.query_name in probes_mapped_to_candidate
                and probes_mapped_to_candidate[entry.query_name]
            ):
                continue

            if record_contains_expected_snp(entry):
                probes_mapped_to_candidate[entry.query_name] = True
                slices_containing_mutation.add(candidate)
            else:
                probes_mapped_to_candidate[entry.query_name] = False

    results = dict(
        slices_containing_mutation=len(slices_containing_mutation),
        total_slices=len(candidate_files),
        variant_sites_denovo_ran_on=len(probes_mapped_to_candidate),
        variant_sites_denovo_correctly_discovered=sum(
            probes_mapped_to_candidate.values()
        ),
    )

    return results


def any_entries_overlap(entries: List[pysam.VariantRecord]) -> bool:
    """Checks whether any VCF entry's alleles overlap"""
    if len(entries) > 1:
        for i in range(len(entries) - 1):
            this_entry = entries[i]
            next_entry = entries[i + 1]
            if (this_entry.start + get_variant_length(this_entry)) > next_entry.start:
                print(f"{this_entry.pos} overlaps {next_entry.pos}")
                if get_genotype(this_entry) != get_genotype(next_entry):
                    print(this_entry.chrom)
                    print(
                        f"Genotypes are different: {get_genotype(this_entry)} {get_genotype(next_entry)}"
                    )
                # return True

    return False


def assign_variants_to_intervals(
    variants: List[pysam.VariantRecord], intervals: List[Tuple[int, int]]
) -> Dict[Tuple[int, int], pysam.VariantRecord]:
    """Assigns each variant to an interval based on the variant's position.
    If a variant doesn't fall within an interval it is not added."""
    interval_to_variants = dict()
    for variant in variants:
        i = find_index_in_intervals(intervals, variant.start)

        if i == -1:
            print(
                "WARNING: variant does not fall within any intervals.", file=sys.stderr
            )
            continue

        if intervals[i] not in interval_to_variants:
            interval_to_variants[intervals[i]] = variant
        else:
            if get_genotype_confidence(variant) > get_genotype_confidence(
                interval_to_variants[intervals[i]]
            ):
                interval_to_variants[intervals[i]] = variant

    return interval_to_variants


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

    query_probes_path = Path(snakemake.output.query_probes)
    
    with query_probes_path.open("w") as fh:
        fh.write(query_probes)

    alignment = map_panel_to_probes(reference_panel, query_probes_path, snakemake.threads)
    probe_results = dict(pandora_calls=get_results_from_alignment(alignment))
    # todo: convert probe_results into expected result format
    probe_results["total_reference_sites"] = snakemake.wildcards.num_snps

    panel = reference_panel.read_text()
    candidate_files = list(Path(snakemake.input.denovo_dir).rglob("*.fa"))
    candidate_results = evaluate_candidates(candidate_files, panel, snakemake.threads)

    results = {**probe_results, **candidate_results}

    write_results(results, Path(snakemake.output.results))


if __name__ == "__main__":
    main()
