#!/usr/bin/env python3

import argparse
import os

from SigProfilerExtractor import sigpro as sig


def main(args):
    try:
        sig.sigProfilerExtractor(
            'vcf',
            args.out_dir,
            args.input,
            minimum_signatures=args.min_signatures,
            maximum_signatures=args.max_signatures,
            context_type='96',
            exome=args.exome,
            cpu=args.cpu,
        )
    except KeyError:
        out_file = os.path.join(args.out_dir, 'SBS96', 'Suggested_Solution',
            'COSMIC_SBS96_Decomposed_Solution',
            'De_Novo_map_to_COSMIC_SBS96.csv')
        if not os.path.exists(out_file):
            raise RunTimeError('\n!ERROR! sigProfilerExtractor failed\n')
        else:
            print('success')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str,  help='Input file')
    parser.add_argument('-o', '--out_dir', type=str, default='',
        help='Output directory. Default = <INPUT>/poissonTests_all')
    parser.add_argument('-mis', '--min_signatures', type=int, default=1,
        help='Minimum of signatures.')
    parser.add_argument('-mas', '--max_signatures', type=int, default=10,
        help='Maximum of signatures.')
    parser.add_argument('-n', '--cpu', type=int, default=1,
        help='Number of cpus.')
    parser.add_argument('-e', '--exome', action='store_true', default=False,
        help='Extract exome.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    main(args)