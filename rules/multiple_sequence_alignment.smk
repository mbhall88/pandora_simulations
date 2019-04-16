rule dealign_original_msa:
    input:
        "data/all_gene_alignments/{gene}_na_aln.fa.gz"
    output:
        "data/realignments/{gene}.clustalo.fa"
    threads: 2
    params:
        extra = "--dealign "
    singularity: CONDA_IMG
    log:
        "logs/dealign_original_msa/{gene}.log"
    wrapper:
        "0.32.0/bio/clustalo"


rule add_denovo_paths_to_msa:
    input:
        denovo_dir = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/denovo_paths",
        msa = "data/realignments/{gene}.clustalo.fa"
    output:
        "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/updated_msas/{gene}/msa_with_denovo_paths.fa"
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/{gene}/add_denovo_to_msa.log"
    run:
        gene_paths = list(Path(input.denovo_dir).glob(f"{wildcards.gene}*.fa"))
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
                            fh_out.write(line.rstrip() + p.stem + f"_path{read_counter}\n")
                            read_counter += 1
                        else:
                            fh_out.write(line)


rule run_msa_after_adding_denovo_paths:
    input:
        "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/updated_msas/{gene}/msa_with_denovo_paths.fa"
    output:
        "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/updated_msas/{gene}/msa_with_denovo_paths.clustalo.fa"
    threads: 2
    singularity: CONDA_IMG
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/{gene}/run_msa_after_adding_denovo_paths.log"
    wrapper:
        "0.32.0/bio/clustalo"
