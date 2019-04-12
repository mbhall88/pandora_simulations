rule build_initial_prg:
    input:
        "data/realignments/{gene}.clustalo.fa"
    output:
        prg = "data/prgs/max_nesting_lvl_{max_nesting_lvl}/{gene}.prg.fa",
    params:
        script = "scripts/make_prg_from_msa.py",
        prefix = "data/prgs/max_nesting_lvl_{max_nesting_lvl}/{gene}",
        max_nesting_lvl = "{max_nesting_lvl}",
    singularity:
        "docker://continuumio/miniconda3:4.5.12"
    conda:
        "../envs/make_prg.yaml"
    log:
        "logs/make_prg/max_nesting_lvl_{max_nesting_lvl}/{gene}.log"
    shell:
        """
        python3 {params.script} --max_nesting {params.max_nesting_lvl} \
            --prefix {params.prefix} {input} 2> {log}
        mv {params.prefix}.max_nest{params.max_nesting_lvl}.min_match7.prg {output.prg} 2>> {log}
        echo '>{wildcards.gene}' | cat - {output.prg} > temp && mv temp {output.prg} 2>> {log}
        echo '' >> {output.prg} 2>> {log}
        rm summary.tsv 2>> {log}
        """

rule index_initial_prg:
    input:
        "data/prgs/max_nesting_lvl_{max_nesting_lvl}/{gene}.prg.fa"
    output:
        "data/prgs/max_nesting_lvl_{max_nesting_lvl}/{gene}.prg.fa.k15.w14.idx"
    # singularity:
    #     "shub://mbhall88/Singularity_recipes:pandora@ac594f67db8a2f66e1c5cc049cfe1968"
    log:
        "logs/index_initial_prg/max_nesting_lvl_{max_nesting_lvl}/{gene}.log"
    shell:
        """
        pandora index {input} &> {log}
        """
