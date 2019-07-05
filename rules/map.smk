rule map_with_discovery:
    input:
        prg = "data/prgs/max_nesting_lvl_{max_nesting_lvl}/combined.prg.fa",
        index = "data/prgs/max_nesting_lvl_{max_nesting_lvl}/combined.prg.fa.k15.w14.idx",
        reads = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/reads.simulated.fa",
        ref = "analysis/{max_nesting_lvl}/{num_snps}/combined_random_paths_mutated_1.fasta",
    output:
        directory("analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/denovo_paths"),
        consensus = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/pandora.consensus.fq.gz",
        genotype_vcf = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/pandora_genotyped.vcf",
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000
    params:
        outdir = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/",
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery.log"
    shell:
        """
        genome_size=$(grep -v '>' {input.ref} | wc | awk '{{print $3-$1-10}}')
        read_file=$(realpath {input.reads})
        prg_file=$(realpath {input.prg})
        log_file=$(realpath {log})
        mkdir -p {params.outdir} 
        cd {params.outdir} || exit 1
        
        pandora map --prg_file $prg_file \
            --read_file $read_file \
            --outdir $(pwd) \
            --output_kg \
            -t {threads} \
            --output_covgs \
            --max_covg {wildcards.coverage} \
            --output_vcf \
            --genotype \
            --genome_size $genome_size \
            --discover \
            --denovo_kmer_size {wildcards.denovo_kmer_size} \
            --max_num_kmers_to_average 100000\
            --log_level debug &> $log_file
        """


rule map_without_discovery:
    input:
        prg = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/updated_msas/combined.prg.fa",
        index = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/updated_msas/combined.prg.fa.k15.w14.idx",
        reads = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/reads.simulated.fa",
        ref = "analysis/{max_nesting_lvl}/{num_snps}/combined_random_paths_mutated_1.fasta",
    output:
        consensus = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_without_discovery/pandora.consensus.fq.gz",
        genotype_vcf = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_without_discovery/pandora_genotyped.vcf",
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000
    params:
        outdir = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_without_discovery/",
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_without_discovery.log"
    shell:
        """
        genome_size=$(grep -v '>' {input.ref} | wc | awk '{{print $3-$1-10}}')
        pandora map --prg_file {input.prg} \
            --read_file {input.reads} \
            --outdir {params.outdir} \
            --output_kg \
            --output_covgs \
            -t {threads} \
            --max_covg {wildcards.coverage} \
            --output_vcf \
            --max_num_kmers_to_average 100000 \
            --genotype \
            --genome_size $genome_size \
            --log_level debug &> {log}
        """
