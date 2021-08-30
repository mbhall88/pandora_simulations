rule build_initial_prg:
    input:
        "data/realignments/{gene}.clustalo.fa",
    output:
        prg="data/prgs/max_nesting_lvl_{max_nesting_lvl}/{gene}/prg.fa",
    params:
        prefix="data/prgs/max_nesting_lvl_{max_nesting_lvl}/{gene}/prg",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (1000 + (attempt * 1000)) * attempt,
    container:
        CONTAINERS["make_prg"]
    log:
        "logs/{max_nesting_lvl}/{gene}_build_initial_prg.log",
    shell:
        """
        python3 {params.script} -v --max_nesting {wildcards.max_nesting_lvl} \
            --prefix {params.prefix} {input} 2> {log}
        mv {params.prefix}.max_nest{wildcards.max_nesting_lvl}.min_match7.prg {output.prg} 2>> {log}
        tmp_fname=$(mktemp)
        echo '>{wildcards.gene}' | cat - {output.prg} > "$tmp_fname" && mv "$tmp_fname" {output.prg} 2>> {log}
        echo '' >> {output.prg} 2>> {log}
        rm -f summary.tsv 2>> {log}
        """


rule combine_prgs:
    input:
        expand(
            "data/prgs/max_nesting_lvl_{{max_nesting_lvl}}/{gene}/prg.fa",
            gene=GENE_NAMES,
        ),
    output:
        "data/prgs/max_nesting_lvl_{max_nesting_lvl}/combined.prg.fa",
    threads: 1
    resources:
        mem_mb=500,
    log:
        "logs/{max_nesting_lvl}/combine_prgs.log",
    shell:
        """
        awk 1 {input} > {output} 2> {log}
        """


rule index_initial_combined_prg:
    input:
        "data/prgs/max_nesting_lvl_{max_nesting_lvl}/combined.prg.fa",
    output:
        "data/prgs/max_nesting_lvl_{max_nesting_lvl}/combined.prg.fa.k15.w14.idx",
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000,
    container:
        CONTAINERS["pandora"]
    log:
        "logs/{max_nesting_lvl}/index_initial_combined_prg.log",
    shell:
        """
        pandora index -t {threads} --log_level debug {input} &> {log}
        """


rule build_prg_after_adding_denovo_paths:
    input:
        "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/updated_msas/{gene}/msa_with_denovo_paths.clustalo.fa",
    output:
        prg="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/updated_msas/{gene}/prg.fa",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (1000 + (attempt * 1000)) * attempt,
    params:
        script="scripts/make_prg_from_msa.py",
        prefix="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/updated_msas/{gene}/prg",
    container:
        CONTAINERS["make_prg"]
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/{gene}/build_prg_after_adding_denovo_paths.log",
    shell:
        """
        python3 {params.script} -v --max_nesting {wildcards.max_nesting_lvl} \
            --prefix {params.prefix} {input} 2> {log}
        mv {params.prefix}.max_nest{wildcards.max_nesting_lvl}.min_match7.prg {output.prg} 2>> {log}
        tmp_fname=$(mktemp)
        echo '>{wildcards.gene}' | cat - {output.prg} > "$tmp_fname" && mv "$tmp_fname" {output.prg} 2>> {log}
        echo '' >> {output.prg} 2>> {log}
        rm -f summary.tsv 2>> {log}
        """


rule combine_prgs_after_adding_denovo_paths:
    input:
        expand(
            "analysis/{{max_nesting_lvl}}/{{num_snps}}/{{read_quality}}/{{coverage}}/{{denovo_kmer_size}}/map_with_discovery/updated_msas/{gene}/prg.fa",
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
        "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/updated_msas/combined.prg.fa",
    output:
        "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/updated_msas/combined.prg.fa.k15.w14.idx",
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000,
    container:
        CONTAINERS["pandora"]
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/index_combined_prg_after_adding_denovo_paths.log",
    shell:
        """
        pandora index -t {threads} --log_level debug {input} > {log} 2>&1
        """
