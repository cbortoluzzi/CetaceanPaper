#!/usr/bin/env python



# Author : Chiara Bortoluzzi


import gzip
import argparse
from pathlib import Path



parser = argparse.ArgumentParser(description = 'Filter called variants in VCF')
parser.add_argument('--vcf', help = 'Called variants in VCF format')
parser.add_argument('--minGQ', help = 'Minimum genotype quality to include a genotype [default: 15]', type = int, default = 15)
parser.add_argument('--minDP', help = 'Minimum depth to include a genotype [default: 4]', type = int, default = 4)
parser.add_argument('--maxDP', help = 'Maximum depth to include a genotype [default: 100]', type = int, default = 100)
parser.add_argument('--minQ', help = 'Minimum quality score to include a genotype [default: 20]', type = int, default = 20)




class FilterVCF:

        def filter_deepvariant_vcf(self, vcf_f, minGQ, minDP, maxDP, minQ):
                self.vcf_write = Path(vcf_f.replace('.vcf.gz', '.postfiltering.vcf'))
                with gzip.open(vcf_f, 'rt') as vcf_reader, open(self.vcf_write, 'w') as vcf_writer:
                        for line in vcf_reader:
                                if line.startswith('#'):
                                        line = line.strip()
                                        vcf_writer.write('{}\n'.format(line))
                                else:
                                     	chrom, pos, id, ref, alt, qual, filter, info, format, sample = line.strip().split()
                                        # we are filtering out all sites except bi-allelic sites
                                        if len(ref) == len(alt) == 1:
                                                # exclude sites with a quality score below minQ
                                                if float(qual) >= minQ:
                                                        # filter bi-allelic SNPs based on minDP, maxDP and GQ
                                                        gt, gq, dp, ad, vaf, pl = sample.split(':')
                                                        if (int(dp) >= minDP and int(dp) <= maxDP) and int(gq) >= minGQ:
                                                                vcf_writer.write('{}'.format(line))




if __name__ == "__main__":
        args = parser.parse_args()
        filter_vcf = FilterVCF()
        filter_vcf.filter_deepvariant_vcf(args.vcf, args.minGQ, args.minDP, args.maxDP, args.minQ)

