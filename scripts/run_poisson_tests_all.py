#!/usr/bin/env python3

import argparse
import numpy as np
import os
import pandas as pd
import subprocess


def run_bash(cmd_raw, bsub=True):
    if bsub:
        cmd = f"sbatch -t 30 -p amd-shared --qos amd-shared --mem 2G " \
            f"--wrap '{cmd_raw}'"
    else:
        cmd = cmd_raw

    subp = subprocess.Popen(cmd,
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = subp.communicate()
    subp.wait()

    print(f'Running: {cmd}')
    if not bsub:
        print(str(stdout), str(stderr))
    print('\n')


def run_poisson_disp(vcf_files, exe, out_dir, replace):
    out_file = os.path.join(out_dir, 'Poisson_dispersion_all.tsv')
    if not replace and os.path.exists(out_file):
        return
    cmd = f'python {exe} {" ".join(vcf_files)} -o {out_file} -b'
    run_bash(cmd)


def run_poisson_tree(tree, vcf_file, args, replace, bsub=True, only_name=False):
    path_strs = vcf_file.split(os.path.sep)
    try:
        clock_dir_no = path_strs.index('ClockTest')
    except ValueError:
        dataset = 'unknown'
        subset = 'unknown'
        filters = 'unknown'
    else:
        subset = path_strs[clock_dir_no + 1]
        file_ids = path_strs[-1].split('.')
        dataset = file_ids[0]
        filters = file_ids[1]

    out_file = os.path.join(args.out_dir,
        f'Poisson_tree_{tree}_{dataset}_{subset}_{filters}.tsv')
    if only_name:
        return out_file

    if not replace and os.path.exists(out_file):
        return

    if tree == 'cellphy':
        tree_file = vcf_file + '.raxml.bestTree'
    elif tree == 'scite':
        vcf_dir = os.path.dirname(vcf_file)
        tree_file = os.path.join(vcf_dir, 'scite_dir',
            f'{dataset}.{filters}_ml0.newick')
    else:
        raise RuntimeError(f'Unknown tree file: {tree}')

    if not os.path.exists(tree_file):
        print(f'!WARNING! Missing {tree: >7} tree file: {tree_file}')
        return

    cmd = f'python {args.exe_tree} {vcf_file} {tree_file} -o {out_file} ' \
        f'-e {args.exe_paup} -b'

    run_bash(cmd, bsub)


def merge_datasets(disp_file, tree_files, out_dir):
    if disp_file:
        df_disp = pd.read_csv(disp_file, sep='\t', index_col=[0, 1, 2])
        df_disp.drop(['H0', 'H1', 'hypothesis'], axis=1, inplace=True)
        df_disp.columns = [f'{i}.dispersion' for i in df_disp.columns]

    if tree_files:
        df_trees = pd.DataFrame()
        for tree_file in tree_files:
            if not tree_file or not os.path.exists(tree_file):
                print(f'!WARNING! Missing tree file: {tree_file}')
                continue
            tree = os.path.basename(tree_file).split('_')[2]
            df_tree = pd.read_csv(tree_file, sep='\t', index_col=[0, 1, 2])
            df_tree.drop(['H0', 'H1', 'hypothesis'], axis=1, inplace=True)
            df_tree.columns = [f'{i}.tree.{tree}' for i in df_tree.columns]
            dataset = df_tree.iloc[0].name
            if dataset in df_trees.index:
                df_trees = df_trees.merge(df_tree, how='outer', left_index=True,
                    right_index=True)
            else:
                df_trees.loc[dataset] = np.nan
                df_trees.loc[dataset, df_tree.columns] = df_tree.iloc[0]

    if disp_file and tree_files:
        df = df_disp.merge(df_trees, left_index=True, right_index=True)
    elif disp_file:
        df = df_disp
    else:
        df = df_trees
    dof_cols = [i for i in df.columns if 'dof' in i]
    df.rename({dof_cols[0]: 'dof'}, axis=1, inplace=True)
    df.drop(dof_cols[1:], axis=1, inplace=True)

    out_file = os.path.join(out_dir, 'Summary_biological_data.tsv')
    df.to_csv(out_file, sep='\t')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str,  help='vcf master file')
    parser.add_argument('-o', '--out_dir', type=str,
        default='poisson_tests_all', help='Output file.')
    parser.add_argument('-m', '--mode', type=str, choices=['run', 'merge'],
        default='run', help='Which task to do: run (on hpc)|merge')
    parser.add_argument('-et', '--exe_tree', type=str,
        default='simulations/scripts/get_poisson_tree_LRT.py',
        help='Poisson Tree exe.')
    parser.add_argument('-ed', '--exe_disp', type=str,
        default='simulations/scripts/get_poisson_LRT.py',
        help='Poisson Dispersion exe.')
    parser.add_argument('-p', '--exe_paup', type=str,
        default='../paup4a168_ubuntu64', help='PAUP*  exe.')
    parser.add_argument('-t', '--tests', choices=['both', 'tree', 'dispersion'],
        default='both', help='Which tests to perform.')
    parser.add_argument('-l', '--local', action='store_false',
        help='Run locally instead of HPC.')
    parser.add_argument('-r', '--replace', action='store_true',
        help='Overwrite already existing files.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    with open(args.input, 'r') as f:
        vcf_files = f.read().strip().split('\n')

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    if args.mode == 'run':
        if args.tests == 'both' or args.tests == 'dispersion':
            run_poisson_disp(vcf_files, args.exe_disp, args.out_dir, args.replace)

        if args.tests == 'both' or args.tests == 'tree':
            for vcf_file in vcf_files:
                run_poisson_tree('cellphy', vcf_file, args, args.replace,
                    bsub=args.local)
                run_poisson_tree('scite', vcf_file, args, args.replace,
                    bsub=args.local)
    else:
        if args.tests == 'both' or args.tests == 'dispersion':
            disp_file = os.path.join(args.out_dir, 'Poisson_dispersion_all.tsv')
        else:
            disp_file = None

        tree_files = []
        if args.tests == 'both' or args.tests == 'tree':
            for vcf_file in vcf_files:
                tree_files.append(run_poisson_tree('cellphy', vcf_file, args,
                    args.replace, only_name=True))
                tree_files.append(run_poisson_tree('scite', vcf_file, args,
                    args.replace, only_name=True))
        merge_datasets(disp_file, tree_files, args.out_dir)


