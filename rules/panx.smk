rule download_panx_data:
    output:
        temp("data/all_gene_alignments.tar.gz")
    params:
        url = "http://pangenome.tuebingen.mpg.de/dataset/Escherichia_coli/all_gene_alignments.tar.gz"
    log:
        "logs/download_panx_data.log"
    shell:
        """
            wget {params.url} -O {output} 2> {log}
        """


rule extract_panx_data:
    input:
        rules.download_panx_data.output
    output:
        directory("data/all_gene_alignments")
    log:
        "logs/extract_panx_data.log"
    shell:
        """
        tar -xzf {input} 2> {log}
        """


rule organise_panx_data:
    input:
        rules.extract_panx_data.output
    output:
        "data/all_gene_alignments/organise_panx_data.complete"
    params:
        files_per_dir = 4000
    log:
        "logs/organise_panx_data.log"
    shell:
        """
        bash scripts/organise_panx_data.sh {input} {output} {params.files_per_dir} 2> {log}
        """
