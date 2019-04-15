import unittest
from pathlib import Path
from simulation import Simulation


class TestSimulation(unittest.TestCase):
    def setUp(self):
        gene = "test"
        num_genes = 3
        max_nesting_lvl = 1
        num_snps = 2
        is_imperfect = True
        coverage = 30
        denovo_kmer_size = 11
        self.simulation = Simulation(
            gene,
            num_genes,
            max_nesting_lvl,
            num_snps,
            is_imperfect,
            coverage,
            denovo_kmer_size,
        )

    def test_get_directory_withImperfectReads(self):
        actual = self.simulation.get_directory()
        expected = Path("3") / "1" / "test" / "6" / "imperfect" / "30" / "11"

        self.assertEqual(actual, expected)

    def test_get_directory_withPerfectReads(self):
        self.simulation.read_quality = "perfect"

        actual = self.simulation.get_directory()
        expected = Path("3") / "1" / "test" / "6" / "perfect" / "30" / "11"

        self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
