rule add_denovo_paths_to_msa:
    input:
        denovo_dir="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/denovo_paths",
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
        "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/updated_msas/{gene}/msa_with_denovo_paths.fa",
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
