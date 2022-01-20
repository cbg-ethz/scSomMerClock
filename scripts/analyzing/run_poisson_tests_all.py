#!/usr/bin/env python3

import argparse
from datetime import datetime
import numpy as np
import os
import pandas as pd
import shutil
import subprocess
import tarfile


MODULE_STR = 'module load ete;'


def run_bash(cmd_raw, bsub=True):
    if bsub:
        cmd = f"sbatch -t 30 -p amd-shared --qos amd-shared --mem 2G " \
            f"--wrap '{MODULE_STR} {cmd_raw}'"
    else:
        cmd = f'{MODULE_STR} {cmd_raw}'

    print(f'Running:\n{cmd}')
    subp = subprocess.Popen(cmd,
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = subp.communicate()
    subp.wait()

    if not bsub:
        print(str(stdout), str(stderr))
    print()


def run_poisson_disp(vcf_files, args):
    out_file = os.path.join(args.out_dir, 'Poisson_dispersion_all.tsv')
    if not args.replace and os.path.exists(out_file):
        return
    cmd = f'python {args.exe_disp} {" ".join(vcf_files)} -o {out_file} -b'
    run_bash(cmd, args.local)


def run_poisson_tree(tree, vcf_file, args):
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

    w_max = np.arange(100, 1001, 100)
    w_max_str = ' '.join([str(i) for i in w_max])

    out_file = os.path.join(args.out_dir,
        f'Poisson_tree_{tree}_{dataset}_{subset}_{filters}.tsv')

    if os.path.exists(out_file) and not args.replace:
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
        f'-w {w_max_str} -b'

    run_bash(cmd, args.local)


def run_plotting(vcf_files, args):
    for vcf_file in vcf_files:
        path_strs = vcf_file.split(os.path.sep)
        clock_dir_no = path_strs.index('ClockTest')
        subset = path_strs[clock_dir_no + 1]
        file_ids = path_strs[-1].split('.')
        dataset = file_ids[0]
        filters = file_ids[1]

        w_max = 500

        cellphy_tree = vcf_file + '.raxml.bestTree'
        cellphy_plot = cellphy_tree + f'_w{w_max}_mapped.png'
        if os.path.exists(cellphy_plot) and not args.replace:
            print(f'\tTree file exists: {cellphy_plot}')
        elif not os.path.exists(cellphy_tree):
            print(f'\tMissing tree file: {cellphy_tree}')
        else:
            cmd = f'python {args.exe_tree} {vcf_file} {cellphy_tree} -w {w_max} ' \
                f'-b -p'
            run_bash(cmd, args.local)

        vcf_dir = os.path.dirname(vcf_file)
        scite_tree = os.path.join(vcf_dir, 'scite_dir',
            f'{dataset}.{filters}_ml0.newick')
        scite_plot = scite_tree + f'_w{w_max}_mapped.png'
        if os.path.exists(scite_plot) and not args.replace:
            print(f'\tTree file exists: {scite_plot}')
        elif not os.path.exists(scite_tree):
            print(f'\tMissing tree file: {scite_tree}')
        else:
            cmd = f'python {args.exe_tree} {vcf_file} {scite_tree} -w {w_max} ' \
                f'-b -p'
            run_bash(cmd, args.local)



def merge_datasets(disp_file, tree_files, out_dir):
    if disp_file:
        df_disp = pd.read_csv(disp_file, sep='\t', index_col=[0, 1, 2])
        df_disp.drop(['H0', 'H1', 'hypothesis', '-2logLR'], axis=1, inplace=True)
        df_disp.columns = [f'{i}.dispersion' for i in df_disp.columns]

    if tree_files:
        df_trees = pd.DataFrame()
        for tree_file in tree_files:
            if not tree_file or not os.path.exists(tree_file):
                print(f'!WARNING! Missing tree file: {tree_file}')
                continue
            tree = os.path.basename(tree_file).split('_')[2]
            weight = os.path.basename(tree_file).split('_')[3].lstrip('w')
            df_tree = pd.read_csv(tree_file, sep='\t', index_col=[0, 1, 2])
            drop_cols = [i for i in df_tree.columns \
                if 'hypothesis' in i or '-2logLR' in i]
            df_tree.drop(drop_cols, axis=1, inplace=True)
            df_tree.columns = [f'{i}.tree.{tree}.{weight}' for i in df_tree.columns]
            dataset = df_tree.iloc[0].name
            if dataset in df_trees.index:
                try:
                    df_trees.loc[dataset, df_tree.columns] = df_tree.iloc[0]
                except KeyError:
                    df_trees = pd.concat([df_trees, df_tree], axis=1)
            else:
                try:
                    df_trees.loc[dataset, df_tree.columns] = df_tree.iloc[0]
                except KeyError:
                    df_trees = df_tree

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


def get_poisson_tree_files(vcf_files):
    files = []
    for vcf_file in vcf_files:
        path_strs = vcf_file.split(os.path.sep)
        clock_dir_no = path_strs.index('ClockTest')
        subset = path_strs[clock_dir_no + 1]
        file_ids = path_strs[-1].split('.')
        dataset = file_ids[0]
        filters = file_ids[1]

        files.append(os.path.join(args.out_dir,
            f'Poisson_tree_cellphy_{dataset}_{subset}_{filters}.tsv'))
        files.append(os.path.join(args.out_dir,
            f'Poisson_tree_scite_{dataset}_{subset}_{filters}.tsv'))
    return files


def get_plot_files(vcf_files):
    files = []
    for vcf_file in vcf_files:
        path_strs = vcf_file.split(os.path.sep)
        file_ids = path_strs[-1].split('.')
        dataset = file_ids[0]
        filters = file_ids[1]

        files.append(vcf_file + '.raxml.bestTree')
        files.append(os.path.join(os.path.dirname(vcf_file), 'scite_dir',
            f'{dataset}.{filters}_ml0.newick'))
    return files


def get_summary_files(vcf_files):
    files = []
    for vcf_file in vcf_files:
        path_strs = vcf_file.split(os.path.sep)
        clock_dir_no = path_strs.index('ClockTest')
        subset = path_strs[clock_dir_no + 1]
        file_ids = path_strs[-1].split('.')
        dataset = file_ids[0]
        filters = file_ids[1]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str,  help='vcf master file')
    parser.add_argument('-o', '--out_dir', type=str,
        default='poisson_tests_all', help='Output file.')
    parser.add_argument('-m', '--mode', type=str,
        choices=['run', 'merge', 'compress', 'plot'],
        default='run', help='Which task to do: run|merge|compress|plot')
    parser.add_argument('-et', '--exe_tree', type=str,
        default='simulations/scripts/get_poisson_tree_LRT.py',
        help='Poisson Tree exe.')
    parser.add_argument('-ed', '--exe_disp', type=str,
        default='simulations/scripts/get_poisson_LRT.py',
        help='Poisson Dispersion exe.')
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
            run_poisson_disp(vcf_files, args)

        if args.tests == 'both' or args.tests == 'tree':
            for vcf_file in vcf_files:
                run_poisson_tree('cellphy', vcf_file, args)
                run_poisson_tree('scite', vcf_file, args)
    elif args.mode == 'merge':
        if args.tests == 'both' or args.tests == 'dispersion':
            disp_file = os.path.join(args.out_dir, 'Poisson_dispersion_all.tsv')
        else:
            disp_file = None

        poisson_tree_files = []
        if args.tests == 'both' or args.tests == 'tree':
            poisson_tree_files.extend(get_poisson_tree_files(vcf_files))
        merge_datasets(disp_file, poisson_tree_files, args.out_dir)
    elif args.mode == 'plot':
        run_plotting(vcf_files, args)
    else:
        comp_dir = os.path.join(args.out_dir,
            f'{datetime.now():%Y%m%d_%H:%M:%S}_compressed')
        print(f'Writing files to: {comp_dir}.tar.gz')
        os.mkdir(comp_dir)

        print('Copying: '\
            f'{os.path.join(args.out_dir, "Summary_biological_data.tsv")}')
        shutil.copyfile(os.path.join(args.out_dir, 'Summary_biological_data.tsv'),
            os.path.join(comp_dir, 'Summary_biological_data.tsv'))

        plot_files = []
        if args.tests == 'both' or args.tests == 'tree':
            plot_files.extend(get_plot_files(vcf_files))

        for plot_file in plot_files:
            print(f'Copying: {plot_file + "_w500_mapped.png"}')
            try:
                shutil.copyfile(plot_file + '_w500_mapped.png',
                    os.path.join(comp_dir, plot_file + '_w500_mapped.png'))
            except FileNotFoundError:
                print('\t!WARNING! - Missing')
            except shutil.SameFileError:
                print(f'\t!WARNING! - Double: {plot_file + "_w500_mapped.png"} -' \
                    f' {os.path.join(comp_dir, plot_file + "_w500_mapped.png")}')

        tar = tarfile.open(comp_dir + '.tar.gz', 'w:gz')
        tar.add(comp_dir)
        tar.close()




