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


rule norm_truth_vcf:
    input:
        truth_vcf=rules.evaluate.input.reference_vcf,
        truth_idx=rules.index_mutated_path.output[0],
        truth_ref=rules.join_random_paths_into_single_reference_sequence.output[0],
        ref=rules.get_random_paths_from_prg.output[0],
    output:
        truth_vcf="analysis/{max_nesting_lvl}/{num_snps}/truth.vcf",
    log:
        "logs/{max_nesting_lvl}/{num_snps}/fix_truth_vcf_chrom.log",
    script:
        "../scripts/norm_truth_vcf.py"


rule norm_pandora_vcf:
    input:
        vcf=rules.map_with_discovery.output.genotype_vcf,
        ref=rules.index_random_path.input[0],
        ref_idx=rules.index_random_path.output[0],
    output:
        vcf="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/pandora_genotyped.norm.vcf",
    resources:
        mem_mb=int(0.5 * GB),
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/norm_pandora_vcf.log",
    container:
        CONTAINERS["bcftools"]
    shell:
        """
        ( bcftools view -a -f .,PASS {input.vcf} \
          | bcftools norm -o {output.vcf} -f {input.ref} ) 2> {log}
        """


rule happy_eval:
    input:
        truth_vcf=rules.norm_truth_vcf.output.truth_vcf,
        query_vcf=rules.norm_pandora_vcf.output.vcf,
        ref=rules.index_mutated_path.input[0],
        ref_idx=rules.index_mutated_path.output[0],
    output:
        summary=(
            "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/evaluation/happy/results.summary.csv"
        ),
    resources:
        mem_mb=lambda wildcards, attempt: int(1_024) * attempt,
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/happy_eval.log",
    container:
        CONTAINERS["happy"]
    params:
        opts=" ".join(
            (
                "--set-gt hom",
                "--pass-only",
                "--write-vcf",
                "--leftshift",
            )
        ),
        prefix=lambda wc, output: output.summary.split(".")[0],
        bcftools_opts="-H -f PASS,.",
    shell:
        """
        truth_count=$(bcftools view {params.bcftools_opts} {input.truth_vcf} | wc -l)
        query_count=$(bcftools view {params.bcftools_opts} {input.query_vcf} | wc -l)

        if [ "$truth_count" -eq 0 ] || [ "$query_count" -eq 0 ]; then
          printf 'TP,FN,FP\n0,%d,%d\n' "$truth_count" "$query_count" > {output.summary} 2> {log}
        else
          hap.py {params.opts} -o {params.prefix} -r {input.ref} \
            {input.truth_vcf} {input.query_vcf} 2> {log}
        fi
        """
