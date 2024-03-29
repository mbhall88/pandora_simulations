from pathlib import Path



rule map_without_discovery:
    input:
        prg=rules.index_initial_combined_prg.input[0],
        index=rules.index_initial_combined_prg.output[0],
        reads=rules.subsample_reads.output.reads,
        ref=rules.mutate_random_path.output.sequences,
        vcf_ref=rules.get_random_paths_from_prg.output[0],
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
            --vcf-refs {input.vcf_ref} \
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
        vcf_ref=rules.get_random_paths_from_prg.output[0],
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
            --vcf-refs {input.vcf_ref} \
            {input.prg} {input.reads} &> {log}
        """
