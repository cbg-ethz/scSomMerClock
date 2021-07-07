#!/usr/bin/env python3

import os
import re
from statistics import mean, stdev
from scipy.stats.distributions import chi2


def get_LRT(in_files, out_file, cell_no, alpha=0.05):
    scores = {}
    for in_file in in_files:
        _, run, _, model, _, _ = os.path.basename(in_file).split('.')
        with open(in_file, 'r') as f_score:
            score_raw = f_score.read().strip().split('\n')
        try:
            score = -float(score_raw[-1].split('\t')[1])
        except (ValueError, IndexError):
            print(f'Cannot read: {in_file}')
            continue

        if not run in scores:
            scores[run] ={'clock': -1, 'noClock': -1}
        scores[run][model] = score

    clock = re.search('clock(\d+.?\d*)_', out_file).group(1) == '0'

    out_str = ''
    avg = [[], [], [], [], [], 0]
    for run, run_info in sorted(scores.items()):
        h0 = run_info['clock']
        h1 = run_info['noClock']
        if h0 == -1 or h1 == -1:
            continue

        LR = -2 * (h0 - h1)
        dof = cell_no - 1
        p_val = chi2.sf(LR, dof)

        for i, j in [(0, h0), (1, h1), (2, LR), (3, dof), (4, p_val)]:
            avg[i].append(j)

        if p_val < alpha:
            hyp = 'H1'
            if not clock:
                avg[5] += 1
        else:
            hyp = 'H0'
            if clock:
                avg[5] += 1

        out_str += f'{run}\t{h0:0>5.2f}\t{h1:0>5.2f}\t{LR:0>5.2f}\t{dof}\t' \
            f'{p_val:.2E}\t{hyp}\n'

    avg_line = f'\n-1\t{mean(avg[0]):0>5.2f}\t{mean(avg[1]):0>5.2f}\t' \
        f'{mean(avg[2]):0>5.2f}\t{mean(avg[3]):.2f}\t{mean(avg[4]):.2E}\t' \
        f'{avg[5]}/{len(scores)}\n'

    with open(out_file, 'w') as f_out:
        f_out.write('run\tH0\tH1\t-2logLR\tdof\tp-value\thypothesis\n')
        f_out.write(out_str + avg_line)

    print(avg_line)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, nargs='+', help='Input files')
    parser.add_argument('-o', '--output', type=str, help='Output file.')
    parser.add_argument('-nc', '--no_cells', type=int, help='Number of cells.')
    parser.add_argument('-a', '--alpha', type=float, default=0.05,
        help='Significance threshold. Default = 0.05.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        get_LRT(snakemake.input, snakemake.output[0],
            snakemake.params.no_cells)
    else:
        import argparse
        args = parse_args()
        get_LRT(args.input, args.output, args.no_cells , args.alpha)