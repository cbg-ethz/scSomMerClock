#!/usr/bin/env python3

import argparse
import gzip

BASES = ['A', 'C', 'G', 'T']


def vcf_to_pileup(vcf_file, out_pileup, out_samples='', out_sample_types=''):
    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

    pileup = ''
    with file_stream as f_in:
        for line in f_in:
            # Skip VCF header lines
            if line.startswith('#'):
                # Safe cell/sample names
                if line.startswith('#CHROM'):
                    sample_names = line.strip().split('\t')[9:]
                continue
            # VCF records
            cols = line.strip().split('\t')
            ref = cols[3]

            new_pileup = ''
            for s_rec in cols[9:]:
                s_rec_format = s_rec.split(':')
                try:
                    DP = int(s_rec_format[1])
                except ValueError:
                    DP = 0
                    import pdb; pdb.set_trace()
                else:
                    if DP == 0:
                        new_pileup += '0\t*\t*\t'
                        continue
                s_bases = ''
                for base_i, base_no in enumerate(s_rec_format[2].split(',')):
                    if BASES[base_i] == ref:
                        s_bases += '.' * int(base_no)
                    else:
                        s_bases += BASES[base_i] * int(base_no)

                new_pileup += '{}\t{}\t{}\t'.format(DP, s_bases, '~' * DP)
            pileup += '{}\t{}\t{}\t{}\n' \
                .format(cols[0], cols[1], ref, new_pileup.rstrip())
    
    with open(out_pileup, 'w') as f_pileup:
        f_pileup.write(pileup.rstrip())

    if not out_samples:
        out_samples = '{}.SampleNames.txt'.format(vcf_file)
    with open(out_samples, 'w') as f_names:
        f_names.write('\n'.join(sample_names))

    if not out_sample_types:
        out_sample_types = '{}.SampleTypes.txt'.format(vcf_file)
    with open(out_sample_types, 'w') as f_types:
        f_types.write('\n'.join(
            ['{}\tCT'.format(i) if not i.startswith('healthy') \
                else '{}\tCN'.format(i) for i in sample_names])
        )


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='Input file')
    parser.add_argument('-o', '--output', type=str, help='Output file.')
    parser.add_argument('-os', '--out_samples', type=str,  defaul='',
        help='Output file for sample names. Default=<INPUT>.SampleNames.txt')
    parser.add_argument('-ot', '--out_types', type=str, default='',
        help='Output file for sample types. Default=<INPUT>.SampleTypes.txt')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        vcf_to_pileup(snakemake.input[0], snakemake.output.pileup,
            snakemake.output.samples, snakemake.output.s_types)
    else:
        args = parse_args()
        vcf_to_pileup(args.input, args.output, args.out_samples, args.out_types)