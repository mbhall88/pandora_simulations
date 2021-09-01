rule get_random_paths_from_prg:
    input:
        prg=rules.index_initial_combined_prg.input[0],
        index=rules.index_initial_combined_prg.output[0],
    output:
        "analysis/{max_nesting_lvl}/random_paths.fa",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 500,
    params:
        num_paths=1,
    container:
        CONTAINERS["pandora"]
    log:
        "logs/{max_nesting_lvl}/get_random_paths_from_prg.log",
    shadow:
        "shallow"
    shell:
        """
        pandora random {input.prg} -n {params.num_paths} &> {log}
        gzip -d random_paths.fa.gz && mv random_paths.fa {output} 2>> {log}
        """


rule join_random_paths_into_single_reference_sequence:
    input:
        rules.get_random_paths_from_prg.output[0],
    output:
        "analysis/{max_nesting_lvl}/combined_random_paths.fa",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 500,
    log:
        "logs/{max_nesting_lvl}/join_random_paths_into_single_reference_sequence.log",
    run:
        with open(output[0], "w") as output_fh:
            header = ">combined_genome_of_genes "
            sequence = ""

            for input_path in input:
                with open(input_path) as input_fh:
                    for line in input_fh:
                        if not line.startswith(">"):
                            sequence += line.rstrip()

            output_fh.write(f"{header}\n{sequence}\n")
