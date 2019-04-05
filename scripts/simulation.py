from pathlib import Path


class Simulation:
    """A class to hold all parameters for a simulation run and provide methods
    for generating filepaths for each parameter combination.
    """

    def __init__(
        self,
        gene: str,
        num_genes: int,
        max_nesting_lvl: int,
        num_snps: int,
        is_imperfect: bool,
        coverage: int,
        denovo_kmer_size: int,
    ) -> "Simulation":
        self.gene = gene
        self.num_genes = num_genes
        self.max_nesting_lvl = max_nesting_lvl
        self.num_snps = num_snps
        self.read_quality = "imperfect" if is_imperfect else "perfect"
        self.coverage = coverage
        self.denovo_kmer_size = denovo_kmer_size

    def get_directory(self) -> Path:
        """Generates a directory path describing the paramters of the simulation."""
        path = Path(str(self.num_genes))
        path /= str(self.max_nesting_lvl)
        path /= self.gene
        path /= str(self.num_snps)
        path /= self.read_quality
        path /= str(self.coverage)
        path /= str(self.denovo_kmer_size)
        return path
