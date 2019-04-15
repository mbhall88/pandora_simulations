rule build_initial_prg:
    input:
        "data/realignments/{gene}.clustalo.fa"
    output:
        prg = "data/prgs/max_nesting_lvl_{max_nesting_lvl}/{gene}/prg.fa",
    params:
        script = "scripts/make_prg_from_msa.py",
        prefix = "data/prgs/max_nesting_lvl_{max_nesting_lvl}/{gene}/prg",
        max_nesting_lvl = "{max_nesting_lvl}",
    singularity: CONDA_IMG
    conda:
        "../envs/make_prg.yaml"
    log:
        "logs/{max_nesting_lvl}/{gene}_build_initial_prg.log"
    shell:
        """
        python3 {params.script} --max_nesting {params.max_nesting_lvl} \
            --prefix {params.prefix} {input} 2> {log}
        mv {params.prefix}.max_nest{params.max_nesting_lvl}.min_match7.prg {output.prg} 2>> {log}
        echo '>{wildcards.gene}' | cat - {output.prg} > temp && mv temp {output.prg} 2>> {log}
        echo '' >> {output.prg} 2>> {log}
        rm summary.tsv 2>> {log}
        """

rule combine_prgs:
    input:
        expand(
            "data/prgs/max_nesting_lvl_{{max_nesting_lvl}}/{gene}/prg.fa",
            gene=[extract_gene_name(gene.name) for gene in genes_for_simulation]
        )
    output:
        "data/prgs/max_nesting_lvl_{max_nesting_lvl}/combined.prg.fa",
    log:
        "logs/{max_nesting_lvl}/combine_prgs.log"
    shell:
        """
        cat {input} > {output} 2> {log}
        """

rule index_initial_combined_prg:
    input:
        "data/prgs/max_nesting_lvl_{max_nesting_lvl}/combined.prg.fa",
    output:
        "data/prgs/max_nesting_lvl_{max_nesting_lvl}/combined.prg.fa.k15.w14.idx"
    # singularity:
    #     "shub://mbhall88/Singularity_recipes:pandora@ac594f67db8a2f66e1c5cc049cfe1968"
    log:
        "logs/{max_nesting_lvl}/index_initial_combined_prg.log"
    shell:
        """
        pandora index {input} &> {log}
        """
