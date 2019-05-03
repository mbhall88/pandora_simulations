rule evaluate:
    input:
        reference_vcf = "analysis/{max_nesting_lvl}/{num_snps}/combined_random_paths_mutated.vcf",
        reference_seq = "analysis/{max_nesting_lvl}/{num_snps}/combined_random_paths_mutated_1.fasta",
        query_vcf = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_without_discovery/pandora_genotyped.vcf",
        pandora_consensus = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_without_discovery/pandora.consensus.fq.gz",
    output:
        reference_panel = Path("analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/evaluation/reference_panel.fa"),
        results = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/evaluation/results.json"
    singularity: CONDA_IMG
    conda:
        "../envs/evaluate.yaml"
    script:
        "../scripts/evaluate.py"