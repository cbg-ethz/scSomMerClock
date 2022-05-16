#!/usr/bin/env python3

import os
import re
import pandas as pd


def merge_PAUP(in_files, out_file, alpha=0.05):
    if '_clock0_' in out_file:
        true_hyp = 'H0'
    else:
        true_hyp = 'H1'

    df = pd.DataFrame(columns=['run', '-2logLR', 'p-value', 'hypothesis'])
    for in_file in sorted(in_files):

        run = int(re.search('paup.(\d+)', in_file).group(1))
        with open(in_file, 'r') as f_score:
            for i, line in enumerate(f_score):
                print(line)
                if not '-lnL' in line:
                    continue
                line_vals = re.split('\s+', line.strip())
                lnl_idx = line_vals.index('-lnL')
                pVal_idx = lnl_idx + 4

                try:
                    pVal = float(line_vals[pVal_idx])
                    LR = float(line_vals[pVal_idx - 1])
                except:
                    import pdb; pdb.set_trace()

                if pVal <= alpha:
                    hyp = 'H0'
                else:
                    hyp = 'H1'

                break
        df.loc[df.shape[0]] = [run, LR, pVal, hyp]

    hyp_avg = f'{(df["hypothesis"] == true_hyp).sum()}/{df.shape[0]}'
    df.loc[df.shape[0]] = [-1, df['-2logLR'].mean(), df['p-value'].mean(), hyp_avg]
    df.set_index('run', inplace=True)

    df.round(4).to_csv(out_file, sep='\t', index=True)
    print(df.iloc[-1])


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, nargs='+', help='Input files')
    parser.add_argument('-o', '--output', type=str, help='Output file.')
    parser.add_argument('-a', '--alpha', type=float, default=0.05,
        help='Significance threshold. Default = 0.05.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        merge_PAUP(snakemake.input, snakemake.output[0])
    else:
        import argparse
        args = parse_args()
        merge_PAUP(args.input, args.output, args.alpha)