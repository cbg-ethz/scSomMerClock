#!/usr/bin/env python3

import argparse
import gzip
import numpy as np

import re



def subsample_vcf(vcf_file, out_files, no, reps, skip=[], outg_id=-1):
    if vcf_file.endswith('gz'):
        file_stream = gzip.open(vcf_file, 'rb')
    else:
        file_stream = open(vcf_file, 'r')

    header = ''
    body = [''] * reps

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
                    samples_all = line_cols[9:]
                    sample_no_all = len(samples_all)
                    skip = [i if i >= 0 else sample_no_all + i for i in skip]

                    samples = np.array([j for i, j in enumerate(samples_all) \
                        if i not in skip])
                    if outg_id < 0:
                        outg_id = len(samples) + outg_id
                    samples = np.delete(samples, outg_id)

                    ss_ids = np.zeros((reps, no,))
                    for i in range(reps):
                        ss_ids[i] = np.random.choice(np.arange(len(samples)),
                            size=no, replace=False)
                    ss_ids = np.sort(ss_ids)

                    subsamples = np.append(samples[ss_ids],
                        np.full((reps, 1),'healthycell') , axis=1)
                    ss_ids = np.append(ss_ids, np.full((reps, 1), outg_id), axis=1)
                    sample_no = no + 1

                    for i, subsample in enumerate(subsamples):

                        body[i] += '\t'.join(line_cols[:9] + list(subsample)) + '\n'
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

            # for s_i, s_rec_raw in enumerate(line_cols[9:]):
            #     for i, ss_id in enumerate(ss_ids):
            #         if s_i in ss_id:
            #             body[i] += '\t' + s_rec_raw
            #         body += '\n'

    for i, out_file in enumerate(out_files):
        with gzip.open(out_file, 'wb') as f_out:
            f_out.write(f'{header}{body[i]}'.encode())


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str,
        help='Absolute or relative path(s) to input file(s)')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Path to the output file. Default = <INPUT_DIR>.subsample')
    parser.add_argument('-n', '--no', type=int, required=True,
        help='Size of subsample.')
    parser.add_argument('-r', '--reps', type=int, default=1,
        help='Number of replicates.')
    parser.add_argument('-s', '--skip', nargs='+', type=int, default=[-1],
        help='Samples to skip. Default = [-1] (CellCoal user tree artefact).')
    parser.add_argument('-outg', '--outgroup', type=int, default=-1,
        help='Outgroup sample number AFTER skipping. Default = -1')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        subsample_vcf(snakemake.input[0], snakemake.output,
            int(snakemake.wildcards.sample), snakemake.params.reps,
            skip=[-1], outg_id=-1)
    else:
        args = parse_args()
        if not args.output:
            args.output = args.input + f'.subsample{args.no}'
        subsample_vcf(args.input, args.output, args.no, args.reps,
            skip=args.skip, outg_id=args.outgroup)