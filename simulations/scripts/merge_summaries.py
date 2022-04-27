#!/usr/bin/env python3

import gzip
import os
import re
import pandas as pd


def merge_summaries(in_files, out_file):
    if not in_files:
        raise IOError('!ERROR: No summary files to merge!')

    df = pd.DataFrame()
    sim_dir = in_files[0][:re.search('minGQ\d+', in_files[0]).end()]
    subsampling = False

    for i, in_file in enumerate(sorted(in_files)):
        new_df = pd.read_csv(in_file, sep='\t')
        if 'subsample_size' in new_df.columns:
            new_df.set_index(['run', 'subsample_size', 'subsample_rep'],
                inplace=True)
            subsampling = True
        else:
            new_df.set_index('run', inplace=True, drop=True)
            new_df.rename(index={'Avg.': '-1'}, inplace=True)
            new_df.index = new_df.index.astype(int)
        new_df.dropna(axis=1, inplace=True)

        base_name = os.path.basename(in_file)
        if base_name.split('.')[0] == 'poissonDisp':
            new_df.columns = [f'{i}_poissonDisp' for i in new_df.columns]
        elif base_name.startswith('poissonTree'):
            tree = in_file.split(os.sep)[-2].split('_')[-1]
            new_df.columns = [f'{i}_{tree}' for i in new_df.columns]

        df = pd.concat([df, new_df], axis=1)

    df.sort_index(inplace=True)
    clock = '_clock0_' in out_file
    if not clock:
        cellcoal_log = os.sep.join(
            in_file.split(os.path.sep)[:-3] + ['cellcoal.out'])
        with open(cellcoal_log, 'r') as f:
            for line in f:
                if line.startswith('Data set'):
                    run_no = int(line.strip().split(' ')[2])
                elif line.startswith('  node '):
                    cells_aff = {int(i) for i in \
                        re.search('(?<=\( )[\d ]+(?= \))', line).group().split()}
                    cell_no = len(cells_aff)
                    df.loc[run_no, 'aff. cells'] = cell_no
                    if subsampling:
                        for subs, subs_data in df.loc[run_no].groupby('subsample_size'):
                            for rep_no, _ in subs_data.groupby('subsample_rep'):
                                vcf_file = os.path.join(sim_dir, 'vcf_dir',
                                    f'vcf.{run_no:0>4}.ss{subs}.{rep_no}.gz')

                                cells_sampled = get_samples(vcf_file)
                                ss_id = (run_no, subs, rep_no)
                                df.loc[ss_id, 'aff. cells sampled'] \
                                    = len(cells_sampled & cells_aff)
        if subsampling:
            df.loc[-1, 'aff. cells sampled'] = df['aff. cells sampled'].mean()
        df.loc[-1, 'aff. cells'] = df['aff. cells'].mean()

    snv_cols = [i for i in df.columns if i.startswith('SNVs_')]
    dof_cols = [i for i in df.columns if i.startswith('dof_')]
    H_cols = [i for i in df.columns \
        if i.startswith('H0_') or i.startswith('H1_')]

    df['dof'] = df[dof_cols[0]]
    df.drop(snv_cols + H_cols + dof_cols, inplace=True, axis=1)
    if not clock:
        if subsampling:
            df = df[['aff. cells', 'aff. cells sampled', 'dof'] \
                + list(df.columns[:-3])]
        else:
            df = df[['aff. cells', 'dof'] + list(df.columns[:-2])]
    else:
        df = df[['dof'] + list(df.columns[:-1])]

    import pdb; pdb.set_trace()
    idx = df.index.tolist()
    df.reindex(idx[1:] + [idx[0]]).to_csv(out_file, sep='\t', index=True)


def get_samples(vcf_file):
    with gzip.open(vcf_file, 'rb') as f:
        for line in f:
            line = line.decode()
            if line.startswith('#CHROM'):
                sampled_raw = line.strip().split('\t')[9:-1]
                sampled = {int(i[4:]) for i in sampled_raw}
                break
    return sampled


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, nargs='+', help='Input files')
    parser.add_argument('-o', '--output', type=str, help='Output file.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        merge_summaries(snakemake.input, snakemake.output[0])
    else:
        import argparse
        args = parse_args()
        merge_summaries(args.input, args.output)