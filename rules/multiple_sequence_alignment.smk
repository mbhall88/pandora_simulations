rule dealign_original_msa:
    input:
        "data/all_gene_alignments/{subdir}/{gene}_na_aln.fa.gz"
    output:
        "data/all_gene_alignments/{subdir}/{gene}.clustalo.fa"
    threads: 2
    params:
        extra = "--dealign"
    log:
        "logs/dealign_original_msa/{subdir}/{gene}.log"
    wrapper:
        "0.31.1-12-g7fc736b/bio/clustalo"
