rule simulate_reads:
    input:
        "analysis/{max_nesting_lvl}/{num_snps}/combined_random_paths_mutated_1.fasta"
    output:
        reads = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/reads.simulated.fa",
        log = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/reads.simulated.log",
        errors = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/reads.simulated.errors.txt"
    params:
        profile = "ecoli_R9_1D",
        perfect_reads = lambda wildcards: wildcards.read_quality == "perfect",
        num_reads = 500 * config["num_genes"],
        extra = "--circular --rnf --seed 88 --unalign-rate 0 ",
    singularity: CONDA_IMG
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/simulate_reads.log"
    wrapper:
        "https://bitbucket.org/mbhall88/snakemake-wrappers/raw/f18bb46319967126a4671f2e6f22c5ce0f2702a0/bio/nanosim-h/wrapper.py"
        # "0.32.0/bio/nanosim-h"
