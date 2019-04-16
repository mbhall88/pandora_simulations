import itertools
import random
from pathlib import Path
from scripts.simulation import Simulation


CONDA_IMG = "docker://continuumio/miniconda3:4.5.12"


def extract_gene_name(string: str) -> str:
    """Extract gene name from panX filenames"""
    return string.split("_na_")[0]


def pick_genes_for_simulation(num: int) -> list:
    all_genes = list(Path("data/all_gene_alignments").rglob("*.fa.gz"))
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

# ======================================================
# Rules
# ======================================================
genes_for_simulation = pick_genes_for_simulation(config["num_genes"])
GENE_NAMES = [extract_gene_name(gene.name) for gene in genes_for_simulation]
all_simulations = generate_all_simulations(
    config, [extract_gene_name(gene.name) for gene in genes_for_simulation]
)

files = set()
for sim in all_simulations:
    sim_path = Path("analysis") / sim.get_directory()

    files.add(
        f"{sim_path}/map_without_discovery/pandora_genotyped.vcf"
    )


rule all:
    input: files


rules_dir = Path("rules/")
include: str(rules_dir / "multiple_sequence_alignment.smk")
include: str(rules_dir / "build_prg.smk")
include: str(rules_dir / "random_path.smk")
include: str(rules_dir / "mutate.smk")
include: str(rules_dir / "simulate.smk")
include: str(rules_dir / "map.smk")


#
# """
# Do multiple sequence alignment to align new paths with where they should be
# """
# rule run_msa:
#     input: "analysis/{sample}/simulate/{sample}_with_denovo.fa"
#     output: "analysis/{sample}/simulate/{sample}_with_denovo.msa.fa"
#     log: "logs/run_msa/{sample}.log"
#     shell:
#         """
#         clustalo --infile {input} --outfile {output} -v --force 2> {log}
#         """
#
#
# rule make_denovo_prg:
#     input: "analysis/{sample}/simulate/{sample}_with_denovo.msa.fa"
#     output:
#         prg = "analysis/{sample}/simulate/{sample}_with_denovo.prg.fa",
#     params:
#         script = "scripts/make_prg_from_msa.py",
#         prefix = "analysis/{sample}/simulate/{sample}_with_denovo",
#     log: "logs/make_denovo_prg/{sample}.log"
#     shell:
#         """
#         python3 {params.script} --prefix {params.prefix} {input} 2> {log}
#         mv {params.prefix}.max_nest10.min_match7.prg {output.prg} 2>> {log}
#         echo '>{wildcards.sample}' | cat - {output.prg} > temp && mv temp {output.prg} 2>> {log}
#         echo '' >> {output.prg} 2>> {log}
#         rm summary.tsv 2>> {log}
#         """
#
#
# rule index_denovo_prg:
#     input: "analysis/{sample}/simulate/{sample}_with_denovo.prg.fa"
#     output: "analysis/{sample}/simulate/{sample}_with_denovo.prg.fa.k15.w14.idx"
#     params:
#         pandora = config["pandora"],
#     log: "logs/index_denovo_prg/{sample}.log"
#     shell:
#         """
#         {params.pandora} index {input} &> {log}
#         """
#
#
# rule map_with_denovo:
#     input:
#         prg = "analysis/{sample}/simulate/{sample}_with_denovo.prg.fa",
#         index = "analysis/{sample}/simulate/{sample}_with_denovo.prg.fa.k15.w14.idx",
#         reads = "analysis/{sample}/simulate/simulated.shuffle.fa",
#         ref = "analysis/{sample}/{sample}_random_path_mutated.fa",
#     output:
#         consensus = "analysis/{sample}/simulate/pandora_map_with_denovo/pandora.consensus.fq.gz",
#         genotype_vcf = "analysis/{sample}/simulate/pandora_map_with_denovo/pandora_genotyped.vcf",
#     params:
#         pandora = config["pandora"],
#         out_prefix = "analysis/{sample}/simulate/pandora_map_with_denovo/",
#     log: "logs/map_with_denovo/{sample}.log"
#     shell:
#         """
#         genome_size=$(grep -v '>' {input.ref} | wc | awk '{{print $3-$1-10}}')
#         {params.pandora} map -p {input.prg} \
#             -r {input.reads} \
#             -o {params.out_prefix} \
#             --output_kg \
#             --output_covgs \
#             --output_vcf \
#             --genotype \
#             --genome_size $genome_size &> {log}
#         """
#
# rule consensus_fastq_to_fasta:
#     input: "analysis/{sample}/simulate/pandora_map_with_denovo/pandora.consensus.fq.gz"
#     output: "analysis/{sample}/simulate/pandora_map_with_denovo/pandora.consensus.fa"
#     log: "logs/consensus_fastq_to_fasta/{sample}.log"
#     shell: "fastaq to_fasta -l 0 {input} {output} 2> {log}"
#
#
# """
# Before applying the VCF changes to the fasta file, need to filter the VCF as it
# contains some invalid entries. And run by first indexing the gzipped version of
# the file with `tabix` and then copying the name of this file to equal that of
# the non-bgziped version of the VCF (for `pysam`)
# """
# rule filter_genotype_vcf:
#     input:
#         vcf = "analysis/{sample}/simulate/pandora_map_with_denovo/pandora_genotyped.vcf"
#     output:
#         vcf = "analysis/{sample}/simulate/pandora_map_with_denovo/pandora_genotyped.filtered.vcf.gz",
#         index = "analysis/{sample}/simulate/pandora_map_with_denovo/pandora_genotyped.filtered.vcf.gz.tbi"
#     params:
#         script = "scripts/filter_invalid_vcf_lines.py"
#     log: "logs/filter_genotype_vcf/{sample}.log"
#     shell:
#         """
#         bgzip -c {input.vcf} > {input.vcf}.gz 2> {log}
#         tabix -p vcf {input.vcf}.gz 2>> {log}
#         cp {input.vcf}.gz.tbi {input.vcf}.tbi 2>> {log}
#         python3 {params.script} {input.vcf} | bgzip -c > {output.vcf} 2>> {log}
#         tabix -p vcf {output.vcf} 2>> {log}
#         """
#
#
# """
# In this step we apply the genotyped VCF variants to the pandora consensus
# fasta file - as this is the reference that the VCF variants is with respect to.
# """
# rule consensus_to_genotype:
#     input:
#         vcf = "analysis/{sample}/simulate/pandora_map_with_denovo/pandora_genotyped.filtered.vcf.gz",
#         ref = "analysis/{sample}/simulate/pandora_map_with_denovo/pandora.consensus.fa",
#     output: "analysis/{sample}/simulate/pandora_map_with_denovo/pandora.genotyped.fa"
#     params:
#         script = "scripts/apply_vcf.py"
#     log: "logs/consensus_to_genotype/{sample}.log"
#     shell:
#         """
#         python3 {params.script} --fasta {input.ref} \
#             --vcf {input.vcf} --output {output} 2> {log}
#         """
#
# rule dnadiff:
#     input:
#         final_sequence = "analysis/{sample}/simulate/pandora_map_with_denovo/pandora.genotyped.fa",
#         expected_sequence = "analysis/{sample}/{sample}_random_path_mutated.fa",
#     output:
#         report = "analysis/{sample}/simulate/dnadiff/out.report",
#         snps_file = "analysis/{sample}/simulate/dnadiff/out.snps"
#     params:
#         prefix = "analysis/{sample}/simulate/dnadiff/out"
#     log: "logs/dnadiff/{sample}.log"
#     shell:
#         """
#         dnadiff -p {params.prefix} {input.expected_sequence} {input.final_sequence} 2> {log}
#         """


# rule consensus_fastq_to_fasta:
#     input:
#     output:
#     params:
#     log:
#     shell:
