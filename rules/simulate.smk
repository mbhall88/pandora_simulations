# generates the reference to simulate reads from
rule join_mutated_random_paths:
    input:
        expand(
                "analysis/{{max_nesting_lvl}}/{gene}/{{num_snps}}/random_path_mutated_1.fasta",
                gene=[extract_gene_name(gene.name) for gene in genes_for_simulation],
        )
    output:
        "analysis/{max_nesting_lvl}/combined_mutated_random_paths_with_{num_snps}_snps.fa"
    log:
        "logs/{max_nesting_lvl}/join_mutated_random_paths_{num_snps}_snps.log"
    run:
        from pathlib import Path

        with open(output[0], "w") as output_fh:
            header = ">combined_genome_of_genes "
            sequence = ""

            for input_path in input:
                current_gene = Path(input_path).parts[2]
                header += f"{current_gene} "

                with open(input_path) as input_fh:
                    for line in input_fh:
                        if not line.startswith(">"):
                            sequence += line.rstrip()

            output_fh.write(f"{header}\n{sequence}\n")




# rule simulate_reads:
#     input: "analysis/{sample}/{sample}_random_path_mutated.fa"
#     output: "analysis/{sample}/simulate/simulated.fa"
#     params:
#         profile = "ecoli_R9_1D",
#         perfect = "--perfect" if config["perfect"] else "",
#         num_reads = 500,
#         min_len = 200,
#         prefix = "analysis/{sample}/simulate/simulated",
#     log: "logs/simulate_reads/{sample}.log"
#     shell:
#         """
#         max_len=$(grep -v '>' {input} | wc | awk '{{print $3-$1-10}}')
#         nanosim-h --number {params.num_reads} \
#             --rnf \
#             {params.perfect} \
#             --profile {params.profile} \
#             --max-len $max_len \
#             --min-len {params.min_len} \
#             --out-pref {params.prefix} \
#             {input} 2> {log}
#         """
