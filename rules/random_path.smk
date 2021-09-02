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
        r"""
        pandora random {input.prg} -n {params.num_paths} &> {log}
        ( gzip -d -c random_paths.fa.gz | sed -r "s/^(>.+)_0$/\1/" ) > {output} 2>> {log}
        """

rule index_random_path:
    input:
        rules.get_random_paths_from_prg.output[0],
    output:
        "analysis/{max_nesting_lvl}/random_paths.fa.fai",
    params:
        "-f" # optional params string
    wrapper:
        "0.77.0/bio/samtools/faidx"

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
                    length = 0
                    name = next(input_fh)[1:].split()[0]
                    for line in input_fh:
                        if not line.startswith(">"):
                            sequence += line.rstrip()
                            length += len(line.rstrip())
                        else:
                            header += f"name={name};len={length} "
                            name = line[1:].split()[0]
                            length = 0
                    header += f"name={name};len={length} "

            output_fh.write(f"{header}\n{sequence}\n")
