import pysam
import sys
from pathlib import Path

def main():
    vcf_path = Path(sys.argv[1])

    # print out the VCF header as pysam adds a bunch of shit to it
    with vcf_path.open() as fh:
        for line in fh:
            if not line.startswith('#'):
                break
            print(line.strip())

    with pysam.VariantFile(vcf_path) as vcf:
        for record in vcf:
            assert len(record.samples.keys()) == 1
            
            # if genotype is ref or empty (.) then skip
            if any(gt in record.samples['sample']['GT'] for gt in (0, None)):
                continue
            elif record.alts is None or record.ref == '.' or '.' in record.alts:
                continue

            sys.stdout.write(str(record))



if __name__ == "__main__":
    main()
