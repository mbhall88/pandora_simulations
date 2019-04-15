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


# rule add_denovo_paths_to_msa:
#     input:
#         denovo_dir = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/denovo_paths",
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
