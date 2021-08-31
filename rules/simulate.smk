rule simulate_reads:
    input:
        "analysis/{max_nesting_lvl}/{num_snps}/combined_random_paths_mutated_1.fasta",
    output:
        reads="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/reads.simulated.fa",
        log="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/reads.simulated.log",
        errors="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/reads.simulated.errors.txt",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
    params:
        profile="/hps/nobackup/research/zi/mbhall/nanosim_training/ecoli/profile",
        perfect_reads=lambda wildcards: wildcards.read_quality == "perfect",
        num_reads=config["num_reads_to_simulate"],
        extra="--circular --rnf --seed 88 --unalign-rate 0 ",
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/simulate_reads.log",
    wrapper:
        "0.33.0/bio/nanosim-h"


rule subsample_reads:
    input:
        reads=rules.simulate_reads.output.reads,
        ref=rules.mutate_random_path.output.sequences,
    output:
        reads="analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/reads.simulated.{coverage}x.fa",
    container:
        CONTAINERS["rasusa"]
    log:
        "logs/subsample_reads/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}x.log",
    resources:
        mem_mb=int(1_000 * 0.5),
    shell:
        """
        genome_size=$(grep -v '>' {input.ref} | wc | awk '{{print $3-$1-10}}') 2> {log}
        rasusa -c {wildcards.coverage} -g $genome_size -i {input.reads} -o {output.reads} 2>> {log}
        """
