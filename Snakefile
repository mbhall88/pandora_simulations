from pathlib import Path

RULES_DIRECTORY = Path("rules/")

# ======================================================
# Config files
# ======================================================
configfile: "config.yaml"


# ======================================================
# Functions and Classes
# ======================================================
convert = {"A": "C", "C": "T", "G": "A", "T": "G"}
def mutate(base):
    base = base.upper()
    return convert[base]

rule all:
    input:
        "data/all_gene_alignments/organise_panx_data.complete",
        # expand(
        #     "analysis/{sample}/simulate/pandora_map_with_discovery/pandora_genotyped.filtered.vcf.gz",
        #     sample=config["sample"]
        # )

include: str(RULES_DIRECTORY / "panx.smk")

#
# rule make_msa:
#     output:
#         msa = "analysis/{sample}/{sample}.msa.fa",
#     params:
#         num_variants = 3,
#     run:
#         import requests
#         import random
#
#         gene = wildcards.sample
#         url = f"http://pangenome.tuebingen.mpg.de/dataset/Escherichia_coli/geneCluster/{gene}_na_aln_reduced.fa"
#         r = requests.get(url)
#         fasta = r.content.decode("utf-8").split("\n")
#         seq = ""
#         # get first sequence from MSA
#         for line in fasta[1:]:
#             if line.startswith(">"):
#                 break
#             seq += line.strip()
#
#         # remove gaps
#         seq = seq.replace("-", "")
#
#         # generate random variant positions
#         random.seed(1)
#         variant_positions = sorted(
#             random.sample(range(len(seq)), k=params.num_variants)
#             )
#         sequences = []
#         for pos in variant_positions:
#             new_base = mutate(seq[pos])
#             new_seq = list(seq)
#             new_seq[pos] = new_base
#             sequences.append("".join(new_seq))
#
#         with open(output.msa, "w") as fout:
#             print(f">consensus\n{seq}", file=fout)
#             for i, s in enumerate(sequences):
#                 variant_pos = variant_positions[i]
#                 original_base = seq[variant_pos]
#                 new_base = mutate(original_base)
#                 print(f">alt{i} {original_base}{variant_pos}{new_base}\n{s}",
#                       file=fout)
#
#
# rule make_prg:
#     input: "analysis/{sample}/{sample}.msa.fa"
#     output:
#         prg = "analysis/{sample}/{sample}.prg.fa",
#     params:
#         script = "scripts/make_prg_from_msa.py",
#         prefix = "analysis/{sample}/{sample}",
#     log: "logs/make_prg/{sample}.log"
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
# rule index_prg:
#     input: "analysis/{sample}/{sample}.prg.fa"
#     output: "analysis/{sample}/{sample}.prg.fa.k15.w14.idx"
#     params:
#         pandora = config["pandora"],
#     log: "logs/index_prg/{sample}.log"
#     shell:
#         """
#         {params.pandora} index {input} &> {log}
#         """
#
#
# rule random_path:
#     input:
#         prg = "analysis/{sample}/{sample}.prg.fa",
#         index = "analysis/{sample}/{sample}.prg.fa.k15.w14.idx"
#     output: "analysis/{sample}/{sample}_random_path.fa"
#     params:
#         pandora = config["pandora"],
#         num_paths = 1,
#     log: "logs/random_path/{sample}.log"
#     shell:
#         """
#         {params.pandora} random_path {input.prg} {params.num_paths} &> {log}
#         gzip -d random_paths.fa.gz && mv random_paths.fa {output} 2>> {log}
#         """
#
#
# rule mutate_random_path:
#     input: "analysis/{sample}/{sample}_random_path.fa"
#     output:
#         fasta = "analysis/{sample}/{sample}_random_path_mutated.fa",
#         variants = "analysis/{sample}/random_path_snps.tsv"
#     log: "logs/random_path/{sample}.log"
#     params:
#         num_snps = config["num_snps"]
#     run:
#         import random
#         snp_info_fh = open(output.variants, "w")
#         print("pos\tref\talt", file=snp_info_fh)
#         seq = ""
#         with open(input[0]) as fasta:
#             # get first sequence from MSA
#             for line in fasta:
#                 if line.startswith(">"):
#                     header = line.strip()
#                     continue
#                 seq += line.strip()
#
#         # generate random snps
#         snps = random.sample(range(len(seq)), k=params.num_snps)
#         new_seq = list(seq)
#         for pos in snps:
#             new_base = mutate(seq[pos])
#             new_seq[pos] = new_base
#             print("{}\t{}\t{}".format(pos, seq[pos], new_base), file=snp_info_fh)
#
#         with open(output.fasta, "w") as fout:
#             print("{}\n{}\n".format(header, "".join(new_seq)), file=fout)
#
#         snp_info_fh.close()
#
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
#
# rule shuffle:
#     input: "analysis/{sample}/simulate/simulated.fa"
#     output: "analysis/{sample}/simulate/simulated.shuffle.fa"
#     log: "logs/shuffle/{sample}.log"
#     shell:
#         """
#         fastaq to_fasta -l 0 {input} - | paste - - | sort -R | awk '{{print $1; print $2}}' > {output} 2> {log}
#         """
#
# rule map_with_discovery:
#     input:
#         prg = "analysis/{sample}/{sample}.prg.fa",
#         index = "analysis/{sample}/{sample}.prg.fa.k15.w14.idx",
#         reads = "analysis/{sample}/simulate/simulated.shuffle.fa",
#         ref = "analysis/{sample}/{sample}_random_path_mutated.fa",
#     output:
#         directory("analysis/{sample}/simulate/pandora_map_with_discovery/denovo_paths"),
#         consensus = "analysis/{sample}/simulate/pandora_map_with_discovery/pandora.consensus.fq.gz",
#         genotype_vcf = "analysis/{sample}/simulate/pandora_map_with_discovery/pandora_genotyped.vcf",
#     params:
#         pandora = config["pandora"],
#         out_prefix = "analysis/{sample}/simulate/pandora_map_with_discovery/",
#     log: "logs/map_with_discovery/{sample}.log"
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
#             --genome_size $genome_size \
#             --discover \
#             --log_level debug &> {log}
#         """
#
#
# """
# Add the de novo paths to the original multiple sequence alignment fasta file.
# """
# rule add_denovo_to_msa:
#     input:
#         denovo_dir = "analysis/{sample}/simulate/pandora_map_with_discovery/denovo_paths",
#         msa = "analysis/{sample}/{sample}.msa.fa",
#     output: "analysis/{sample}/simulate/{sample}_with_denovo.fa"
#     log: "logs/add_denovo_to_msa/{sample}.log"
#     run:
#         denovo_paths = list(Path(input.denovo_dir).glob("*.fa"))
#         # append contents of all denovo paths to the msa file
#         with open(output[0], "w") as fout:
#             fout.write(Path(input.msa).read_text())
#             for p in denovo_paths:
#                 fout.write(p.read_text())
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
