#!/usr/bin/env python3

import os
from statistics import mean, stdev


def merge_LRT(in_files, out_file):
    out_str = ''
    avg = [[], [], [], [], 0, 0]
    for i, in_file in enumerate(in_files):
        with open(in_file, 'r') as f_score:
            score_raw = f_score.read().strip().split('\n')

        # Copy header line
        if i == 0:
            out_str += score_raw[0] + '\n'
        out_str += score_raw[1] + '\n'

        details = score_raw[1].split('\t')
        for i, j in [(0, 1), (1, 2), (2, 3), (3, 5),]:
            avg[i].append(float(details[j]))

        if details[-1] == 'H0':
            avg[4] += 1
        else:
            avg[5] += 1

    avg_line = f'\nAvg.\t{mean(avg[0])}\t{mean(avg[1])}\t{mean(avg[2])}\t\t' \
        f'{mean(avg[3])}\tH0:{avg[4]};H1:{avg[5]}\n'

    with open(out_file, 'w') as f_out:
        f_out.write(out_str + avg_line)

    print(avg_line)
    import pdb; pdb.set_trace()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, nargs='+', help='Input files')
    parser.add_argument('-o', '--output', type=str, help='Output file.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        merge_LRT(snakemake.input, snakemake.output[0])
    else:
        import argparse
        args = parse_args()
        merge_LRT(args.input, args.output)