rule map_with_discovery:
    input:
        prg = "data/prgs/max_nesting_lvl_{max_nesting_lvl}/combined.prg.fa",
        index = "data/prgs/max_nesting_lvl_{max_nesting_lvl}/combined.prg.fa.k15.w14.idx",
        reads = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/reads.simulated.fa",
        ref = "analysis/{max_nesting_lvl}/{num_snps}/combined_mutated_random_paths.fa"
    output:
        directory("analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/denovo_paths"),
        consensus = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/pandora.consensus.fq.gz",
        genotype_vcf = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/pandora_genotyped.vcf",
    params:
        out_prefix = "analysis/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery/",
    log:
        "logs/{max_nesting_lvl}/{num_snps}/{read_quality}/{coverage}/{denovo_kmer_size}/map_with_discovery.log"
    shell:
        """
        genome_size=$(grep -v '>' {input.ref} | wc | awk '{{print $3-$1-10}}')
        pandora map -p {input.prg} \
            -r {input.reads} \
            -o {params.out_prefix} \
            --output_kg \
            --output_covgs \
            --output_vcf \
            --genotype \
            --genome_size $genome_size \
            --discover \
            --log_level debug &> {log}
        """
