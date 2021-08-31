rule evaluate:
    input:
        reference_vcf=rules.mutate_random_path.output.vcf,
        reference_seq=rules.mutate_random_path.output.sequences,
        query_vcf=rules.map_with_discovery.output.genotype_vcf,
        pandora_consensus=rules.map_with_discovery.output.consensus,
        denovo_dir=rules.pandora_discover.output.denovo_dir,
    output:
        reference_panel="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/evaluation/reference_panel.fa",
        query_probes="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/evaluation/query_probes.fa",
        results="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/evaluation/results.json",
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000,
    conda:
        "../envs/evaluate.yaml"
    script:
        "../scripts/evaluate.py"


rule evaluate_no_denovo:
    input:
        reference_vcf=rules.mutate_random_path.output.vcf,
        reference_seq=rules.mutate_random_path.output.sequences,
        query_vcf=rules.map_without_discovery.output.genotype_vcf,
        pandora_consensus=rules.map_without_discovery.output.consensus,
    output:
        reference_panel="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/evaluate_no_denovo/reference_panel.fa",
        query_probes="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/evaluate_no_denovo/query_probes.fa",
        results="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/evaluate_no_denovo/results.json",
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000,
    conda:
        "../envs/evaluate.yaml"
    script:
        "../scripts/evaluate.py"
