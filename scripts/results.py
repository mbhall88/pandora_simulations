import json
from pathlib import Path
from typing import List


class VariantCall:
    def __init__(self, name, correct, mismatches):
        self.name = name
        parts = name.split("_")
        self.confidence = float(parts[-1].replace("CONF", ""))
        self.entry = int(parts[-2].replace("entry", ""))
        self.pos = int(parts[-3].replace("pos", ""))
        self.gene = name.split("_pos")[0]
        self.correct = bool(correct)
        self.mismatches = int(mismatches)

    def __eq__(self, other):
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)


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

    def load_data_from_json(self, path: Path):
        with path.open() as json_file:
            self.data = json.load(json_file)

    def unique_snps_called(self) -> int:
        return self.data.get("reference_sites_called", 0)

    def get_variant_calls(self) -> List[VariantCall]:
        try:
            variant_calls = [
                VariantCall(name, correct, mismatch)
                for name, correct, mismatch in zip(
                    self.data["ids"],
                    self.data["snps_called_correctly"],
                    self.data["mismatches"],
                )
            ]
        except KeyError:
            variant_calls = []

        return variant_calls

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
