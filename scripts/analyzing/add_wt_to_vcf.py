#!/usr/bin/env python3

import argparse
import gzip
import os


def add_outg_sample(vcf_file, out_file, outg_name):
    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

    fmt_map = {'GT': '0|0', 'DP': '100', 'AD': '100,0', 'RC': '100,0,0,0',
        'GQ': '999'}

    with open(out_file, 'w') as f_out:
        with file_stream as f_in:
            for line in f_in:
                try:
                    line = line.decode()
                except:
                    pass
                # Skip VCF header lines
                if line.startswith('#'):
                    # Safe column headers
                    if line.startswith('#CHROM'):
                        f_out.write(line.strip() + f'\t{outg_name}\n')
                    else:
                        f_out.write(line)
                    continue

                line_cols = line.strip().split('\t')
                outg_entry = ''
                for fmt_id in line_cols[8].split(':'):
                    try:
                        outg_entry += f'{fmt_map[fmt_id]}:'
                    except KeyError:
                        outg_entry += '.:'
                f_out.write(line.strip() + '\t' + outg_entry.rstrip(':') + '\n')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_file', type=str,
        help='VCF file.')
    parser.add_argument('-o', '--out_file', type=str,
        default='', help='Output file.')
    parser.add_argument('-n', '--outg_name', type=str,
        default='healthycell', help='Outgroup name. Default = healthycell.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    if not (args.out_file):
        if args.in_file.endswith('.vcf'):
            args.out_file = args.in_file[:-4] + '_outg.vcf'
        elif args.in_file.endswith('.vcf.gz'):
            args.out_file = args.in_file[:-7]  + '_outg.vcf'
        else:
            args.out_file = args.in_file + '_outg.vcf'

    add_outg_sample(args.in_file, args.out_file, args.outg_name)