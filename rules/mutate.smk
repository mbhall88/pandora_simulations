rule mutate_random_path:
    input:
        rules.join_random_paths_into_single_reference_sequence.output[0],
    output:
        sequences="analysis/{max_nesting_lvl}/{num_snps}/combined_random_paths_mutated_1.fasta",
        vcf="analysis/{max_nesting_lvl}/{num_snps}/combined_random_paths_mutated.vcf",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
    params:
        num_simulations=1,
        extra=lambda wildcards: (
            " ".join(
                [
                    "--num-substitutions {}".format(wildcards.num_snps),
                    "--num-insertions 0",
                    "--num-deletions 0",
                    "--random-seed 88",  # make analysis reproducible
                ]
            ),
        ),
    conda:
        "../envs/mutate.yaml"
    log:
        "logs/{max_nesting_lvl}/{num_snps}/mutate_random_paths.log",
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

rule index_mutated_path:
    input:
        rules.mutate_random_path.output.sequences,
    output:
        "analysis/{max_nesting_lvl}/{num_snps}/combined_random_paths_mutated_1.fasta.fai",
    params:
        "-f" # optional params string
    wrapper:
        "0.77.0/bio/samtools/faidx"