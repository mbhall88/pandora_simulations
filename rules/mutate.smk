rule mutate_random_path:
    input:
        "analysis/{max_nesting_lvl}/combined_random_paths.fa"
    output:
        sequences = "analysis/{max_nesting_lvl}/{num_snps}/combined_random_paths_mutated_1.fasta",
        vcf = "analysis/{max_nesting_lvl}/{num_snps}/combined_random_paths_mutated.vcf"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 1000
    params:
        num_simulations = 1,
        extra = lambda wildcards: (
            " ".join([
                "--num-substitutions {}".format(wildcards.num_snps),
                "--num-insertions 0",
                "--num-deletions 0",
                "--random-seed 88",  # make analysis reproducible
            ]),
        )
    singularity: CONDA_IMG
    conda:
        "../envs/mutate.yaml"
    log:
        "logs/{max_nesting_lvl}/{num_snps}/mutate_random_paths.log"
    shell:
        """
        input=$(realpath {input})
        log=$(realpath {log})
        vcf=$(realpath {output.vcf})
        cd $(dirname {output.sequences}) || exit 1
        snpmutator --num-simulations {params.num_simulations} \
            {params.extra} \
            --vcf $vcf \
            $input 2> $log
        """
