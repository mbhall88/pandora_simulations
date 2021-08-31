from pathlib import Path


rule build_initial_prg:
    input:
        "data/all_gene_alignments/{gene}_na_aln.fa.gz",
    output:
        prg="data/prgs/max_nesting_lvl_{max_nesting_lvl}/{gene}.prg",
    params:
        prg_name=lambda wildcards, output: Path(output.prg).with_suffix(""),
        opts="-O p -N {max_nesting_lvl} -S {gene}",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (1000 + (attempt * 1000)) * attempt,
    container:
        CONTAINERS["make_prg"]
    log:
        "logs/{max_nesting_lvl}/{gene}_build_initial_prg.log",
    shadow:
        "shallow"
    shell:
        "make_prg from_msa {params.opts} -n {params.prg_name} {input} 2> {log}"


rule combine_prgs:
    input:
        expand(
            "data/prgs/max_nesting_lvl_{{max_nesting_lvl}}/{gene}.prg",
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
        "awk 1 {input} > {output} 2> {log}"


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
        "pandora index -t {threads} -v {input} &> {log}"
