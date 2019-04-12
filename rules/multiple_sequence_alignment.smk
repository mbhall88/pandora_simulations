rule dealign_original_msa:
    input:
        "data/all_gene_alignments/{gene}_na_aln.fa.gz"
    output:
        "data/realignments/{gene}.clustalo.fa"
    threads: 2
    params:
        extra = "--dealign"
    singularity: CONDA_IMG
    log:
        "logs/dealign_original_msa/{gene}.log"
    wrapper:
        "0.31.1-12-g7fc736b/bio/clustalo"
