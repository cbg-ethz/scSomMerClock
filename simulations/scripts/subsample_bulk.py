#!/usr/bin/env python3

import argparse
import gzip
import numpy as np

import re


def subsample_vcf(vcf_file, prefix, subsamples, reps, skip=[], outg_id=-1):
    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

    header = ''
    ss_ids = []
    body = [''] * len(subsamples) * reps

    with file_stream as f_in:
        for line in f_in:
            try:
                line = line.decode()
            except AttributeError:
                pass
            # Skip VCF header lines
            if line.startswith('#'):
                # Safe column headers
                if line.startswith('#CHROM'):
                    line_cols = line.strip().split('\t')
                    samples_all = [i.strip() for i in line_cols[9:]]
                    sample_no_all = len(samples_all)
                    skip = [i if i >= 0 else sample_no_all + i for i in skip]

                    samples = np.array([j for i, j in enumerate(samples_all) \
                        if i not in skip])
                    if outg_id < 0:
                        outg_id = len(samples) + outg_id
                    samples = np.delete(samples, outg_id)

                    idx = 0
                    for subsample in subsamples:
                        for rep in range(reps):
                            ss_id = np.random.choice(np.arange(len(samples)),
                                    size=subsample, replace=False)
                            ss_id = np.sort(ss_id)

                            ss_name = np.append(samples[ss_id], ['healthycell'])
                            ss_id = np.append(ss_id, outg_id)

                            body[idx] += '\t'.join(line_cols[:9] + list(ss_name)) + '\n'
                            ss_ids.append(ss_id)
                            idx += 1
                else:
                    header += line
                continue
            elif line.strip() == '':
                break

            # VCF records
            line_cols = line.strip().split('\t')
            snv_cols = '\t'.join(line_cols[:9])
            for i, ss_id in enumerate(ss_ids):
                sample_cols = '\t'.join(
                    [k for j, k in enumerate(line_cols[9:]) if j in ss_id])
                body[i] += f'{snv_cols}\t{sample_cols}\n'

    idx = 0
    for sample in subsamples:
        for rep in range(reps):
            out_file = f'{prefix}.ss{sample}.{rep}.gz'
            with gzip.open(out_file, 'wb') as f_out:
                f_out.write(f'{header}{body[idx]}'.encode())
            idx += 1


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str,
        help='Absolute or relative path(s) to input file(s)')
    parser.add_argument('-p', '--prefix', type=str, default='',
        help='Prefix to output files. Default = <INPUT>.ss<NO>.<REP>.gz')
    parser.add_argument('-n', '--no', nargs='+', type=int, required=True,
        help='Size(s) of subsample(s).')
    parser.add_argument('-r', '--reps', type=int, default=1,
        help='Number of replicates.')
    parser.add_argument('-s', '--skip', nargs='+', type=int,
        help='Samples to skip. Default = [] (CellCoal user tree artefact).')
    parser.add_argument('-outg', '--outgroup', type=int, default=-1,
        help='Outgroup sample number AFTER skipping. Default = -1')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        subsample_vcf(
            vcf_file=snakemake.input[0],
            prefix=snakemake.params.prefix,
            subsamples=snakemake.params.subsamples,
            reps=snakemake.params.reps,
            outg_id=-1
        )
    else:
        args = parse_args()
        if not args.prefix:
            args.prefix = args.input
        if not args.skip:
            args.skip = []

        subsample_vcf(
            vcf_file=args.input,
            prefix=args.prefix,
            subsamples=args.no,
            reps=args.reps,
            outg_id=args.outgroup,
            skip=args.skip,
        )