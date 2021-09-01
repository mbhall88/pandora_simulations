
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


rule add_denovo_paths_to_msa:
    input:
        denovo_dir=rules.pandora_discover.output.denovo_dir,
        msa_dir="data/all_gene_alignments/",
    output:
        msa_dir=directory(
            "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/updated_msas/"
        ),
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: 16_000 * attempt,
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/add_denovo_to_msa.log",
    conda:
        "../envs/update_msas.yaml"
    params:
        script="../scripts/add_denovo_paths_to_msa.py",
        options="",
    shell:
        """
        python {params.script} {params.options} -o {output.msa_dir} \
            -j {threads} -M {input.msa_dir} {input.denovo_dir} 2> {log}
        """


rule build_prg_after_adding_denovo_paths:
    input:
        rules.add_denovo_paths_to_msa.output.msa_dir,
    output:
        outdir=directory(
            "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/updated_prgs/"
        ),
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: (16_000 + (attempt * 1000)) * attempt,
    params:
        opts=["-O", "p", "-N", "{max_nesting_lvl}"],
    container:
        CONTAINERS["make_prg"]
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/build_prg_after_adding_denovo_paths.log",
    script:
        "../scripts/update_prgs.py"


rule combine_prgs_after_adding_denovo_paths:
    input:
        prg_dir=rules.build_prg_after_adding_denovo_paths.output.outdir,
    output:
        "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/updated.prg.fa",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 500,
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/combine_prgs_after_adding_denovo_paths.log",
    conda:
        "../envs/fd.yaml"
    shell:
        "fd -e prg -d 1 -X awk 1 \; . {input.prg_dir} > {output[0]} 2> {log}"


rule index_combined_prg_after_adding_denovo_paths:
    input:
        rules.combine_prgs_after_adding_denovo_paths.output[0],
    output:
        rules.combine_prgs_after_adding_denovo_paths.output[0] + ".k15.w14.idx",
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000,
    container:
        CONTAINERS["pandora"]
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/index_combined_prg_after_adding_denovo_paths.log",
    shell:
        "pandora index -t {threads} -v {input} > {log} 2>&1"
