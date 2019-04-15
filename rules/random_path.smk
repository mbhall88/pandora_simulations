rule get_random_paths_from_prg:
    input:
        prg = "data/prgs/max_nesting_lvl_{max_nesting_lvl}/combined.prg.fa",
        index = "data/prgs/max_nesting_lvl_{max_nesting_lvl}/combined.prg.fa.k15.w14.idx"
    output:
        "analysis/{max_nesting_lvl}/random_paths.fa"
    params:
        num_paths = 1,
    # singularity:
    #     "shub://mbhall88/Singularity_recipes:pandora@ac594f67db8a2f66e1c5cc049cfe1968"
    log:
        "logs/{max_nesting_lvl}/get_random_paths_from_prg.log"
    shell:
        """
        pandora random_path {input.prg} {params.num_paths} &> {log}
        gzip -d random_paths.fa.gz && mv random_paths.fa {output} 2>> {log}
        """


rule join_random_paths_into_single_reference_sequence:
    input:
        "analysis/{max_nesting_lvl}/random_paths.fa"
    output:
        "analysis/{max_nesting_lvl}/combined_random_paths.fa"
    log:
        "logs/{max_nesting_lvl}/join_random_paths_into_single_reference_sequence.log"
    run:
        from pathlib import Path

        with open(output[0], "w") as output_fh:
            header = ">combined_genome_of_genes "
            sequence = ""

            for input_path in input:
                with open(input_path) as input_fh:
                    for line in input_fh:
                        if not line.startswith(">"):
                            sequence += line.rstrip()

            output_fh.write(f"{header}\n{sequence}\n")
