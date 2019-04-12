rule mutate_random_path:
    input:
        "analysis/{max_nesting_lvl}/{gene}/random_path.fa"
    output:
        sequences = "analysis/{max_nesting_lvl}/{gene}/{num_snps}/random_path_mutated_1.fasta",
        vcf = "analysis/{max_nesting_lvl}/{gene}/{num_snps}/random_path_mutated.vcf"
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
        "logs/{max_nesting_lvl}/{gene}/{num_snps}/mutate_random_path.log"
    wrapper:
        "COMMIT/bio/snp-mutator"  # TODO: add commit when PR is merged for wrapper
