#!/usr/bin/env python3

import argparse
import os

from SigProfilerExtractor import sigpro as sig


def main(args):
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