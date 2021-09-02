import sys

sys.stderr = open(snakemake.log[0], "w")

in_truth_vcf = snakemake.input.truth_vcf
out_truth_vcf = snakemake.output.truth_vcf
ref = snakemake.input.ref
truth_faidx = snakemake.input.truth_idx
truth_ref = snakemake.input.truth_ref

with open(truth_faidx) as f:
    CHROM, length = next(f).split("\t")[:2]
    vcf_contig_header = f"##contig=<ID={CHROM},length={length}>"

pos_idx = dict()
with open(truth_ref) as fa:
    start = 0
    comments = next(fa).split()[1:]
    for com in comments:
        name = com.split(";")[0].split("=")[-1]
        l = int(com.split(";")[1].split("=")[-1])
        end = start + l
        pos_idx[name] = (start, end)
        start = end

assert end == length, f"Reference length {length} does not equal sum of genes {end}"

with open(in_truth_vcf) as fin, open(out_truth_vcf, "w") as fout:
    written_extra = False
    for line in fin:
        if not written_extra and line[:2].count("#") == 1:
            for gene, iv in pos_idx.items():
                print(f"##contig=<ID={gene},length={iv[-1]-iv[0]}>", file=fout)
            written_extra = True
        elif line[0] != "#":
            cols = line.split("\t")
            pos = int(cols[1]) - 1
            flag = False
            for gene, (s, e) in pos_idx.items():
                if s <= pos < e:
                    cols[0] = gene
                    cols[1] = str(pos - s + 1)
                    flag = True
            if not flag:
                raise ValueError(f"Couldnt find gene for {cols}")
            fout.write("\t".join(cols))
        else:
            fout.write(line)
