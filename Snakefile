from pathlib import Path

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
        expand(
            "analysis/{sample}/simulate_perfect/pandora_map_with_discovery/pandora.genotyped.fa",
            sample=config["sample"]
            ),
        # expand(
        #     "analysis/{sample}/simulate_perfect/pandora_map_with_discovery/pandora_genotyped.filtered.vcf.gz",
        #     sample=config["sample"]
        # )


rule make_msa:
    output:
        msa = "analysis/{sample}/{sample}.msa.fa",
    params:
        num_variants = 3,
    run:
        import requests
        import random

        gene = wildcards.sample
        url = f"http://pangenome.tuebingen.mpg.de/dataset/Escherichia_coli/geneCluster/{gene}_na_aln_reduced.fa"
        r = requests.get(url)
        fasta = r.content.decode("utf-8").split("\n")
        seq = ""
        # get first sequence from MSA
        for line in fasta[1:]:
            if line.startswith(">"):
                break
            seq += line.strip()

        # remove gaps
        seq = seq.replace("-", "")

        # generate random variant positions
        random.seed(1)
        variant_positions = sorted(
            random.sample(range(len(seq)), k=params.num_variants)
            )
        sequences = []
        for pos in variant_positions:
            new_base = mutate(seq[pos])
            new_seq = list(seq)
            new_seq[pos] = new_base
            sequences.append("".join(new_seq))

        with open(output.msa, "w") as fout:
            print(f">consensus\n{seq}", file=fout)
            for i, s in enumerate(sequences):
                variant_pos = variant_positions[i]
                original_base = seq[variant_pos]
                new_base = mutate(original_base)
                print(f">alt{i} {original_base}{variant_pos}{new_base}\n{s}",
                      file=fout)


rule make_prg:
    input: "analysis/{sample}/{sample}.msa.fa"
    output:
        prg = "analysis/{sample}/{sample}.prg.fa",
    params:
        script = "scripts/make_prg_from_msa.py",
        prefix = "analysis/{sample}/{sample}",
    log: "logs/make_prg/{sample}.log"
    shell:
        """
        python3 {params.script} --prefix {params.prefix} {input} 2> {log}
        mv {params.prefix}.max_nest10.min_match7.prg {output.prg} 2>> {log}
        echo '>{wildcards.sample}' | cat - {output.prg} > temp && mv temp {output.prg} 2>> {log}
        echo '' >> {output.prg} 2>> {log}
        rm summary.tsv 2>> {log}
        """


rule index_prg:
    input: "analysis/{sample}/{sample}.prg.fa"
    output: "analysis/{sample}/{sample}.prg.fa.k15.w14.idx"
    params:
        pandora = config["pandora"],
    log: "logs/index_prg/{sample}.log"
    shell:
        """
        {params.pandora} index {input} &> {log}
        """


rule random_path:
    input:
        prg = "analysis/{sample}/{sample}.prg.fa",
        index = "analysis/{sample}/{sample}.prg.fa.k15.w14.idx"
    output: "analysis/{sample}/{sample}_random_path.fa"
    params:
        pandora = config["pandora"],
        num_paths = 1,
    log: "logs/random_path/{sample}.log"
    shell:
        """
        {params.pandora} random_path {input.prg} {params.num_paths} &> {log}
        gzip -d random_paths.fa.gz && mv random_paths.fa {output} 2>> {log}
        """


rule mutate_random_path:
    input: "analysis/{sample}/{sample}_random_path.fa"
    output:
        fasta = "analysis/{sample}/{sample}_random_path_mutated.fa",
        vcf = "analysis/{sample}/snpmutator.vcf"
    params:
        snps = 3,
        simulations = 1,
        insertions = 0,
        deletions = 0,
        random_seed = 1,
    log: "logs/random_path/{sample}.log"
    shell:
        """
        snpmutator --num-simulations {params.simulations} \
            --num-substitutions {params.snps} \
            --num-insertions {params.insertions} \
            --num-deletions {params.deletions} \
            --random-seed {params.random_seed} \
            --vcf {output.vcf} \
            {input} 2> {log}
        mv {wildcards.sample}_random_path_mutated_1.fasta {output.fasta} 2>> {log}
        """


rule simulate_reads:
    input: "analysis/{sample}/{sample}_random_path_mutated.fa"
    output: "analysis/{sample}/simulate_perfect/simulated.fa"
    params:
        profile = "ecoli_R9_1D",
        num_reads = 500,
        min_len = 200,
        prefix = "analysis/{sample}/simulate_perfect/simulated",
    log: "logs/simulate_reads/{sample}.log"
    shell:
        """
        max_len=$(grep -v '>' {input} | wc | awk '{{print $3-$1-10}}')
        nanosim-h --number {params.num_reads} \
            --rnf \
            --perfect \
            --profile {params.profile} \
            --max-len $max_len \
            --min-len {params.min_len} \
            --out-pref {params.prefix} \
            {input} 2> {log}
        """


rule map_with_discovery:
    input:
        prg = "analysis/{sample}/{sample}.prg.fa",
        index = "analysis/{sample}/{sample}.prg.fa.k15.w14.idx",
        reads = "analysis/{sample}/simulate_perfect/simulated.fa",
        ref = "analysis/{sample}/{sample}_random_path_mutated.fa",
    output:
        directory("analysis/{sample}/simulate_perfect/pandora_map_with_discovery/denovo_paths"),
        consensus = "analysis/{sample}/simulate_perfect/pandora_map_with_discovery/pandora.consensus.fq.gz",
        genotype_vcf = "analysis/{sample}/simulate_perfect/pandora_map_with_discovery/pandora_genotyped.vcf",
    params:
        pandora = config["pandora"],
        out_prefix = "analysis/{sample}/simulate_perfect/pandora_map_with_discovery/",
    log: "logs/map_with_discovery/{sample}.log"
    shell:
        """
        genome_size=$(grep -v '>' {input.ref} | wc | awk '{{print $3-$1-10}}')
        {params.pandora} map -p {input.prg} \
            -r {input.reads} \
            -o {params.out_prefix} \
            --output_kg \
            --output_covgs \
            --output_vcf \
            --genotype \
            --genome_size $genome_size \
            --discover &> {log}
        """


rule consensus_fastq_to_fasta:
    input: "analysis/{sample}/simulate_perfect/pandora_map_with_discovery/pandora.consensus.fq.gz"
    output: "analysis/{sample}/simulate_perfect/pandora_map_with_discovery/pandora.consensus.fa"
    log: "logs/consensus_fastq_to_fasta/{sample}.log"
    shell: "fastaq to_fasta -l 0 {input} {output} 2> {log}"


"""
Before applying the VCF changes to the fasta file, need to filter the VCF as it
contains some invalid entries. And run by first indexing the gzipped version of
the file with `tabix` and then copying the name of this file to equal that of
the non-bgziped version of the VCF (for `pysam`)
"""
rule filter_genotype_vcf:
    input:
        vcf = "analysis/{sample}/simulate_perfect/pandora_map_with_discovery/pandora_genotyped.vcf"
    output:
        vcf = "analysis/{sample}/simulate_perfect/pandora_map_with_discovery/pandora_genotyped.filtered.vcf.gz",
        index = "analysis/{sample}/simulate_perfect/pandora_map_with_discovery/pandora_genotyped.filtered.vcf.gz.tbi"
    params:
        script = "scripts/filter_invalid_vcf_lines.py"
    log: "logs/filter_genotype_vcf/{sample}.log"
    shell:
        """
        bgzip -c {input.vcf} > {input.vcf}.gz 2> {log}
        tabix -p vcf {input.vcf}.gz 2>> {log}
        cp {input.vcf}.gz.tbi {input.vcf}.tbi 2>> {log}
        python3 {params.script} {input.vcf} | bgzip -c > {output.vcf} 2>> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """



"""
In this step we apply the genotyped VCF variants to the pandora consensus
fasta file - as this is the reference that the VCF variants is with respect to.
"""
rule consensus_to_genotype:
    input:
        vcf = "analysis/{sample}/simulate_perfect/pandora_map_with_discovery/pandora_genotyped.filtered.vcf.gz",
        ref = "analysis/{sample}/simulate_perfect/pandora_map_with_discovery/pandora.consensus.fa",
    output: "analysis/{sample}/simulate_perfect/pandora_map_with_discovery/pandora.genotyped.fa"
    params:
        script = "scripts/apply_vcf.py"
    log: "logs/consensus_to_genotype/{sample}.log"
    shell:
        """
        python3 {params.script} --fasta {input.ref} \
            --vcf {input.vcf} --output {output} 2> {log}
        """
# rule consensus_fastq_to_fasta:
#     input:
#     output:
#     params:
#     log:
#     shell:
