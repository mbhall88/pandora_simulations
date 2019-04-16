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
    params:
        outdir = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/",
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery.log"
    shell:
        """
        genome_size=$(grep -v '>' {input.ref} | wc | awk '{{print $3-$1-10}}')
        pandora map --prg_file {input.prg} \
            --read_file {input.reads} \
            --outdir {params.outdir} \
            --output_kg \
            --output_covgs \
            --max_covg {wildcards.coverage} \
            --output_vcf \
            --genotype \
            --genome_size $genome_size \
            --discover \
            --denovo_kmer_size {wildcards.denovo_kmer_size} \
            --log_level debug &> {log}
        """
