rule evaluate:
    input:
        reference_vcf = "analysis/{max_nesting_lvl}/{num_snps}/combined_random_paths_mutated.vcf",
        reference_seq = "analysis/{max_nesting_lvl}/{num_snps}/combined_random_paths_mutated_1.fasta",
        query_vcf = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_without_discovery/pandora_genotyped.vcf",
        pandora_consensus = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_without_discovery/pandora.consensus.fq.gz",
        denovo_dir = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/denovo_paths",
    output:
        reference_panel = Path("analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/evaluation/reference_panel.fa"),
        query_probes = Path("analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/evaluation/query_probes.fa"),
        results = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/evaluation/results.json"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000
    singularity: CONDA_IMG
    conda:
        "../envs/evaluate.yaml"
    script:
        "../scripts/evaluate.py"


rule evaluate_no_denovo:
    input:
        reference_vcf = "analysis/{max_nesting_lvl}/{num_snps}/combined_random_paths_mutated.vcf",
        reference_seq = "analysis/{max_nesting_lvl}/{num_snps}/combined_random_paths_mutated_1.fasta",
        query_vcf = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/pandora_genotyped.vcf",
        pandora_consensus = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/pandora.consensus.fq.gz",
    output:
        reference_panel = Path("analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/evaluate_no_denovo/reference_panel.fa"),
        query_probes = Path("analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/evaluate_no_denovo/query_probes.fa"),
        results = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/evaluate_no_denovo/results.json"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000
    singularity: CONDA_IMG
    conda:
        "../envs/evaluate.yaml"
    script:
        "../scripts/evaluate.py"
