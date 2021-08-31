from pathlib import Path


rule pandora_discover:
    input:
        prg=rules.index_initial_combined_prg.input[0],
        index=rules.index_initial_combined_prg.output[0],
        reads=rules.subsample_reads.output.reads,
        ref=rules.mutate_random_path.output.sequences,
    output:
        denovo_dir=directory(
            "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/discover/denovo_paths"
        ),
        consensus="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/discover/pandora.consensus.fq.gz",
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000,
    params:
        outdir=lambda wildcards, output: Path(output.denovo_dir).parent,
        opts=" ".join(["-v", "--discover-k", "{denovo_kmer_size}"]),
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery.log",
    container:
        CONTAINERS["pandora"]
    shadow:
        "shallow"
    shell:
        """
        genome_size=$(grep -v '>' {input.ref} | wc | awk '{{print $3-$1-10}}')

        pandora discover {params.opts} \
            -o {params.outdir} \
            -t {threads} \
            -g $genome_size \
            {input.prg} {input.reads} &> {log}
        """


rule map_without_discovery:
    input:
        prg=rules.index_initial_combined_prg.input[0],
        index=rules.index_initial_combined_prg.output[0],
        reads=rules.subsample_reads.output.reads,
        ref=rules.mutate_random_path.output.sequences,
    output:
        consensus="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_without_discovery/pandora.consensus.fq.gz",
        genotype_vcf="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_without_discovery/pandora_genotyped.vcf",
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000,
    params:
        outdir=lambda wildcards, output: Path(output.genotype_vcf).parent,
        opts=" ".join(["-v", "--genotype"]),
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_without_discovery.log",
    container:
        CONTAINERS["pandora"]
    shell:
        """
        genome_size=$(grep -v '>' {input.ref} | wc | awk '{{print $3-$1-10}}')

        pandora map {params.opts} \
            -o {params.outdir} \
            -t {threads} \
            -g $genome_size \
            {input.prg} {input.reads} &> {log}
        """


rule map_with_discovery:
    input:
        prg=rules.index_combined_prg_after_adding_denovo_paths.input[0],
        index=rules.index_combined_prg_after_adding_denovo_paths.output[0],
        reads=rules.subsample_reads.output.reads,
        ref=rules.mutate_random_path.output.sequences,
    output:
        consensus="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/pandora.consensus.fq.gz",
        genotype_vcf="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/pandora_genotyped.vcf",
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000,
    params:
        outdir=lambda wildcards, output: Path(output.genotype_vcf).parent,
        opts=" ".join(["-v", "--genotype"]),
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery.log",
    container:
        CONTAINERS["pandora"]
    shell:
        """
        genome_size=$(grep -v '>' {input.ref} | wc | awk '{{print $3-$1-10}}')
        pandora map {params.opts} \
            -o {params.outdir} \
            -t {threads} \
            -g $genome_size \
            {input.prg} {input.reads} &> {log}
        """
