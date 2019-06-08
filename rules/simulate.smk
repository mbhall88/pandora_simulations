rule simulate_reads:
    input:
        "analysis/{max_nesting_lvl}/{num_snps}/combined_random_paths_mutated_1.fasta"
    output:
        reads = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/reads.simulated.fa",
        log = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/reads.simulated.log",
        errors = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/reads.simulated.errors.txt"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 1000
    params:
        profile = "/hps/nobackup/research/zi/mbhall/nanosim_training/ecoli/profile",
        perfect_reads = lambda wildcards: wildcards.read_quality == "perfect",
        num_reads = 500 * config["num_genes"],
        extra = "--circular --rnf --seed 88 --unalign-rate 0 ",
    singularity: CONDA_IMG
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/simulate_reads.log"
    wrapper:
        "0.33.0/bio/nanosim-h"
