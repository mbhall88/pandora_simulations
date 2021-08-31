
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
        msa="data/all_gene_alignments/{gene}_na_aln.fa.gz",
    output:
        "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/updated_msas/{gene}/msa_with_denovo_paths.fa",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 500,
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/{gene}/add_denovo_to_msa.log",
    run:
        gene_paths = list(Path(input.denovo_dir).glob(f"{wildcards.gene}.*.fa"))
        new_msa_path = Path(output[0])

        if not new_msa_path.parent.is_dir():
            new_msa_path.parent.mkdir(parents=True, exist_ok=True)

        with new_msa_path.open("w") as fh_out:
            original_msa = Path(input.msa)
            fh_out.write(original_msa.read_text())

            for p in gene_paths:
                read_counter = 1

                with p.open() as fasta:
                    for line in fasta:
                        if line.startswith(">"):
                            fh_out.write(
                                line.rstrip() + p.stem + f"_path{read_counter}\n"
                            )
                            read_counter += 1
                        else:
                            fh_out.write(line)


rule run_msa_after_adding_denovo_paths:
    input:
        rules.add_denovo_paths_to_msa.output[0],
    output:
        msa="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/updated_msas/{gene}/msa_with_denovo_paths.msa.fa",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
    container:
        CONTAINERS["mafft"]
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/{gene}/run_msa_after_adding_denovo_paths.log",
    shell:
        "mafft --auto --thread {threads} {input} > {output.msa} 2> {log}"


rule build_prg_after_adding_denovo_paths:
    input:
        rules.run_msa_after_adding_denovo_paths.output.msa,
    output:
        prg="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/updated_msas/{gene}.prg",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (1000 + (attempt * 1000)) * attempt,
    params:
        prg_name=lambda wildcards, output: Path(output.prg).with_suffix(""),
        opts="-O p -N {max_nesting_lvl} -S {gene}",
    container:
        CONTAINERS["make_prg"]
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/{gene}/build_prg_after_adding_denovo_paths.log",
    shell:
        "make_prg from_msa {params.opts} -n {params.prg_name} {input} 2> {log}"


rule combine_prgs_after_adding_denovo_paths:
    input:
        expand(
            "analysis/{{max_nesting_lvl}}/{{num_snps}}/{{read_quality}}/{{coverage}}/{{denovo_kmer_size}}/map_with_discovery/updated_msas/{gene}.prg",
            gene=GENE_NAMES,
        ),
    output:
        "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/updated_msas/combined.prg.fa",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 500,
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/combine_prgs_after_adding_denovo_paths.log",
    shell:
        """
        awk 1 {input} > {output} 2> {log}
        """


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
