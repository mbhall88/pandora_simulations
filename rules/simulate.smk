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


rule simulate_reads:
    input:
        "analysis/{max_nesting_lvl}/combined_mutated_random_paths_with_{num_snps}_snps.fa"
    output:
        reads = "analysis/{max_nesting_lvl}/{gene}/{num_snps}/{read_quality}/reads.simulated.fa",
        log = "analysis/{max_nesting_lvl}/{gene}/{num_snps}/{read_quality}/reads.simulated.log",
        errors = "analysis/{max_nesting_lvl}/{gene}/{num_snps}/{read_quality}/reads.simulated.errors.txt"
    params:
        profile = "ecoli_R9_1D",
        perfect_reads = lambda wildcards: wildcards.read_quality == "perfect",
        num_reads = 500 * config["num_genes"],
        extra = "--circular --rnf --seed 88 ",
    singularity: CONDA_IMG
    log:
        "logs/{max_nesting_lvl}/{gene}/{num_snps}/{read_quality}/simulate_reads.log"
    wrapper:
        "0.32.0/bio/nanosim-h"
