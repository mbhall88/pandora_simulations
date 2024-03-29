import itertools
import random
from pathlib import Path
from scripts.simulation import Simulation


def extract_gene_name(string: str) -> str:
    """Extract gene name from panX filenames"""
    return string.split("_na_")[0]


def pick_genes_for_simulation(num: int, exclude: set = None) -> list:
    if exclude is None:
        exclude = set()
    all_genes = sorted(
        [
            p
            for p in Path("data/all_gene_alignments").glob("*.fa")
            if extract_gene_name(p.name) not in exclude
        ]
    )
    random.seed(88)
    return random.sample(all_genes, num)


def generate_all_simulations(config: dict, genes: list) -> list:

    parameter_combinations = list(
        itertools.product(
            genes,
            [config["num_genes"]],
            config["prg_nesting_lvls"],
            config["snps_per_gene"],
            [True, False],
            config["coverages"],
            config["denovo_kmer_sizes"],
        )
    )

    return [Simulation(*params) for params in parameter_combinations]


# ======================================================
# Config files
# ======================================================
configfile: "config.yaml"


CONTAINERS = config["containers"]
SCRIPTS = Path("scripts").resolve()

# ======================================================
# Rules
# ======================================================
exclude_genes = {"GC00000393_5"}
genes_for_simulation = pick_genes_for_simulation(
    config["num_genes"], exclude=exclude_genes
)
GENE_NAMES = [extract_gene_name(gene.name) for gene in genes_for_simulation]
all_simulations = generate_all_simulations(
    config, [extract_gene_name(gene.name) for gene in genes_for_simulation]
)

files = set()
for sim in all_simulations:
    sim_path = Path("analysis") / sim.get_directory()

    files.add(f"{sim_path}/evaluation/results.json")
    files.add(f"{sim_path}/evaluate_no_denovo/results.json")
    files.add(f"{sim_path}/evaluation/happy/results.summary.csv")


rule all:
    input:
        files,


rules_dir = Path("rules/")


include: str(rules_dir / "build_prg.smk")
include: str(rules_dir / "random_path.smk")
include: str(rules_dir / "mutate.smk")
include: str(rules_dir / "simulate.smk")
include: str(rules_dir / "discover.smk")
include: str(rules_dir / "map.smk")
include: str(rules_dir / "evaluate.smk")
