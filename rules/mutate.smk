rule mutate_random_path:
    input:
        "analysis/{max_nesting_lvl}/{gene}_random_path.fa"
    output:
        fasta = "analysis/{sample}/{sample}_random_path_mutated.fa",
        variants = "analysis/{sample}/random_path_snps.tsv"
    log: "logs/random_path/{sample}.log"
    params:
        num_snps = config["num_snps"]
    shell:
        """
        snpmutator
        """

# rule join_random_paths:
#     input:
#         expand(
#                 "analysis/{{max_nesting_lvl}}/{gene}_random_path.fa",
#                 gene=[extract_gene_name(gene.name) for gene in genes_for_simulation],
#         )
#     output:
#         "analysis/{}"
