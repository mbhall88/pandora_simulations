import json
import re
from pathlib import Path
from collections import defaultdict
from typing import List, Tuple


class RegexError(Exception):
    pass


class VariantCall:
    def __init__(self, name, correct, ref_id):
        self.name = name.strip()
        self.correct = correct
        self.ref_pos = int(ref_id[1:-1])

    def __eq__(self, other):
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)

    def gene(self) -> str:
        regex = re.compile(r"(.+)_POS=")
        match = regex.search(self.name)
        if not match:
            raise RegexError(f"Gene name could not be parsed from {self.name}")

        return match.group(1)

    def pos(self) -> int:
        regex = re.compile(r"POS=(\d+)_interval")
        match = regex.search(self.name)
        if not match:
            raise RegexError(f"POS could not be parsed from {self.name}")

        return int(match.group(1))

    def interval(self) -> Tuple[int, int]:
        regex = re.compile(r"interval=\((\d+),(\d+)\)_GT_CONF")
        match = regex.search(self.name)
        if not match:
            raise RegexError(f"POS could not be parsed from {self.name}")

        return int(match.group(1)), int(match.group(2))

    def gt_conf(self) -> float:
        regex = re.compile(r"GT_CONF=([0-9]+\.?[0-9]*)")
        match = regex.search(self.name)
        if not match:
            raise RegexError(f"GT_CONF could not be parsed from {self.name}")

        return float(match.group(1))


class Result:
    def __init__(
        self, denovo_kmer_size=0, coverage=0, read_quality="", num_snps=0, max_nesting=0
    ):
        self.denovo_kmer_size = denovo_kmer_size
        self.coverage = coverage
        self.read_quality = read_quality
        self.num_snps = num_snps
        self.max_nesting = max_nesting
        self.data = dict()
        self.variant_calls = []

    def __eq__(self, other: "Result"):
        return all(
            [
                self.denovo_kmer_size == other.denovo_kmer_size,
                self.coverage == other.coverage,
                self.read_quality == other.read_quality,
                self.num_snps == other.num_snps,
                self.max_nesting == other.max_nesting,
                self.data == other.data,
            ]
        )

    def __str__(self):
        return f"Denovo kmer size:\t{self.denovo_kmer_size}\nCoverage:\t{self.coverage}\nRead quality:\t{self.read_quality}\nNumber of SNPs:\t{self.num_snps}\nMax nesting:\t{self.max_nesting}\n"

    def load_data_from_json(self, path: Path):
        with path.open() as json_file:
            self.data = json.load(json_file)
        self.variant_calls = self._get_variant_calls()

    def unique_snps_called(self) -> int:
        return len(self.data.get("pandora_calls", {}))

    def total_denovo_slices(self) -> int:
        return self.data.get("total_slices", 0)

    def denovo_slices_containing_variant_site(self) -> int:
        return self.data.get("slices_containing_mutation", 0)

    def variant_sites_denovo_ran_on(self) -> int:
        return self.data.get("variant_sites_denovo_ran_on", 0)

    def variant_sites_denovo_correctly_discovered(self) -> int:
        return self.data.get("variant_sites_denovo_correctly_discovered", 0)

    def _get_variant_calls(self) -> List[VariantCall]:
        calls = self.data.get("pandora_calls", {})
        variant_calls = []
        for ref_id, call_data in calls.items():
            name = call_data["id"]
            correct = call_data["correct"]
            variant_calls.append(VariantCall(name, correct, ref_id))

        return variant_calls

    def denovo_recall(self) -> float:
        """% of simulated mutations where the mutant allele was included in the
        candidates.
        """
        try:
            return self.variant_sites_denovo_correctly_discovered() / self.num_snps
        except ZeroDivisionError:
            return 0.0

    def denovo_precision(self) -> float:
        """% of slices where we perform de novo, that include a simulated-mutation
        within
        """
        try:
            return (
                self.variant_sites_denovo_correctly_discovered()
                / self.total_denovo_slices()
            )
        except ZeroDivisionError:
            return 0.0

    def _calculate_metrics(self, conf_threshold=0):
        position_calls = defaultdict(list)
        for call in self.variant_calls:
            position_calls[call.ref_pos - 1].append((call.correct, call.gt_conf()))

        positive_calls = 0
        true_positives = 0

        for pos, calls in position_calls.items():
            confident_calls = [
                i
                for i, (is_correct, conf) in enumerate(calls)
                if conf >= conf_threshold
            ]
            if confident_calls:
                positive_calls += 1
                if any(calls[i][0] for i in confident_calls):
                    true_positives += 1

        false_positives = positive_calls - true_positives

        return true_positives, false_positives

    def overall_true_positives(self, conf_threshold=0):
        return self._calculate_metrics(conf_threshold)[0]

    def overall_false_negatives(self, conf_threshold=0):
        return self.num_snps - self.overall_true_positives(conf_threshold)

    def overall_false_positives(self, conf_threshold=0):
        return self._calculate_metrics(conf_threshold)[-1]

    def overall_recall(self, conf_threshold=0) -> float:
        """True variants called relative to all variants
        Calculation: TP/TP+FN
        """
        true_positives = self.overall_true_positives(conf_threshold)

        try:
            return true_positives / self.num_snps
        except ZeroDivisionError:
            return 0.0

    def overall_precision(self, conf_threshold=0) -> float:
        """True variants called relative to total calls
        Calculation: TP/TP+FP
        """
        true_positives = self.overall_true_positives(conf_threshold)
        false_positives = self.overall_false_positives(conf_threshold)

        try:
            return true_positives / (true_positives + false_positives)
        except ZeroDivisionError:
            return 0.0

    def overall_accuracy(self, conf_threshold=0):
        """Ratio of correct calls to total calls and variants
        Calculation: TP+TN/TP+FP+TN+FN
        Note: we are not using TN so cancels out of equation
        """
        true_positives = self.overall_true_positives(conf_threshold)
        false_positive = self.overall_false_positives(conf_threshold)
        false_negatives = self.overall_false_negatives(conf_threshold)

        try:
            return true_positives / (true_positives + false_negatives + false_positive)
        except ZeroDivisionError:
            return 0.0

    def overall_error_rate(self, conf_threshold=0):
        return 1 - self.overall_accuracy(conf_threshold)

    def as_dict(self, conf_threshold=0) -> dict:
        """Returns a dict representation.
        The contents of this dict were chosen for the purposes of creating a pandas
        dataframe to use for plotting results.
        """
        return dict(
            denovo_precision=self.denovo_precision(),
            denovo_recall=self.denovo_recall(),
            overall_recall=self.overall_recall(conf_threshold),
            overall_precision=self.overall_precision(conf_threshold),
            overall_accuracy=self.overall_accuracy(conf_threshold),
            overall_error_rate=self.overall_error_rate(conf_threshold),
            denovo_kmer_size=self.denovo_kmer_size,
            num_snps=self.num_snps,
            coverage=self.coverage,
            max_nesting=self.max_nesting,
            read_quality=self.read_quality,
            conf_threshold=conf_threshold,
        )

    @staticmethod
    def extract_parameters_from_path(path: Path) -> dict:
        """Extracts the simulation parameters from a given path.

        Note: Path must end in the following format:
        <rest of path>/<max_nesting>/<num_snps>/<read_quality>/<coverage>/<denovo_kmer_size>/<can be anything>/<filename>

        :param path: Path to extract parameters from.
        :returns: dictionary containing the extracted paramters.
        """
        return dict(
            denovo_kmer_size=int(path.parts[-3]),
            coverage=int(path.parts[-4]),
            read_quality=path.parts[-5],
            num_snps=int(path.parts[-6]),
            max_nesting=int(path.parts[-7]),
        )

    @staticmethod
    def from_dict(data: dict) -> "Result":
        return Result(
            data.get("denovo_kmer_size", 0),
            data.get("coverage", 0),
            data.get("read_quality", ""),
            data.get("num_snps", 0),
            data.get("max_nesting", 0),
        )

    @staticmethod
    def from_json(path: Path) -> "Result":
        data = Result.extract_parameters_from_path(path)
        result = Result.from_dict(data)
        result.load_data_from_json(path)

        return result
