rule get_random_path_from_prg:
    input:
        prg = "data/prgs/max_nesting_lvl_{max_nesting_lvl}/{gene}/prg.fa",
        index = "data/prgs/max_nesting_lvl_{max_nesting_lvl}/{gene}/prg.fa.k15.w14.idx"
    output:
        "analysis/{max_nesting_lvl}/{gene}_random_path.fa"
    params:
        num_paths = 1,
    # singularity:
    #     "shub://mbhall88/Singularity_recipes:pandora@ac594f67db8a2f66e1c5cc049cfe1968"
    log:
        "logs/{max_nesting_lvl}/{gene}_get_random_path_from_prg.log"
    shell:
        """
        pandora random_path {input.prg} {params.num_paths} &> {log}
        gzip -d random_paths.fa.gz && mv random_paths.fa {output} 2>> {log}
        """
