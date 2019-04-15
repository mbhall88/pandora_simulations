rule mutate_random_path:
    input:
        "analysis/{max_nesting_lvl}/{gene}_random_path.fa"
    output:
        sequences = "analysis/{max_nesting_lvl}/{num_snps}/{gene}_random_path_mutated_1.fasta",
        vcf = "analysis/{max_nesting_lvl}/{num_snps}/{gene}_random_path_mutated.vcf"
    params:
        num_simulations = 1,
        extra = (
            " ".join([
                "--num-substitutions {num_snps}",
                "--num-insertions 0",
                "--num-deletions 0",
                "--random-seed 88",  # make analysis reproducible
            ]),
        )
    singularity: CONDA_IMG
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{gene}_mutate_random_path.log"
    wrapper:
        "0.32.0/bio/snp-mutator"
