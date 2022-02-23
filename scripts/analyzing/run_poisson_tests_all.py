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


def run_bash(cmd_raw, bsub=True, module_str=MODULE_STR):
    if bsub:
        cmd = f"sbatch -t 30 -p amd-shared --qos amd-shared --mem 2G " \
            f"--wrap '{module_str} {cmd_raw}'"
    else:
        cmd = f'{module_str} {cmd_raw}'

    print(f'Running:\n{cmd}')
    subp = subprocess.Popen(cmd,
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = subp.communicate()
    subp.wait()

    if not bsub:
        print(str(stdout), str(stderr))
    print()



def run_tests(vcf_files, args):
    if args.tests == 'both' or args.tests == 'dispersion':
        run_poisson_disp(vcf_files, args)
    if args.tests == 'both' or args.tests == 'tree':
        run_poissonTree(vcf_files, args)


def run_poisson_disp(vcf_files, args):
    out_file = os.path.join(args.out_dir, 'Poisson_dispersion_all.tsv')
    if not args.replace and os.path.exists(out_file):
        return
    cmd = f'python {args.exe_disp} {" ".join(vcf_files)} -o {out_file} -b'
    run_bash(cmd, args.local, '')


def run_poissonTree(vcf_files, args):
    for vcf_file in vcf_files:
        for tree in ['cellphy', 'scite']:
            run_poissonTree_single(vcf_file, tree, args)


def run_poissonTree_single(vcf_file, tree, args):
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


def run_plotting(vcf_files, args, gather_only=False):
    w_max = 500

    if gather_only:
        phyl_dir = os.path.join(args.out_dir, 'phylogeny')
        if not os.path.exists(phyl_dir):
            os.makedirs(phyl_dir)

    for vcf_file in vcf_files:
        path_strs = vcf_file.split(os.path.sep)
        clock_dir_no = path_strs.index('ClockTest')
        subset = path_strs[clock_dir_no + 1]
        if subset == 'all':
            continue

        file_ids = path_strs[-1].split('.')
        dataset = file_ids[0]
        filters = file_ids[1]

        for tree in ['cellphy', 'scite']:
            if tree == 'cellphy':
                tree_file = vcf_file + '.raxml.bestTree'
                log_file = vcf_file + '.raxml.log'
            else:
                tree_file = os.path.join(os.path.dirname(vcf_file), 'scite_dir',
                    f'{dataset}.{filters}_ml0.newick')
                log_file = os.path.join(os.path.dirname(vcf_file), 'scite_dir',
                    f'{dataset}.{filters}.log')

            base_name = f'{dataset}_{subset}_{filters.replace("_outg", "")}_{tree}'
            plot_file_raw = os.path.join(args.out_dir, base_name)
            plot_file = plot_file_raw + f'_w{w_max:.0f}_mapped.png'

            if gather_only:
                if not os.path.exists(tree_file):
                    print(f'\tMissing tree file: {tree_file}')
                else:
                    shutil.copyfile(tree_file,
                        os.path.join(phyl_dir, f'{base_name}.newick'))
                    shutil.copyfile(log_file,
                        os.path.join(phyl_dir, f'{base_name}.log'))

                    shutil.copyfile(vcf_file,
                        os.path.join(phyl_dir, f'{base_name}.vcf.gz'))
                continue

            if os.path.exists(plot_file) and not args.replace:
                print(f'\tExisting tree plot: {plot_file}')
            elif not os.path.exists(tree_file):
                print(f'\tMissing tree file: {tree_file}')
            else:
                cmd = f'python {args.exe_tree} {vcf_file} {tree_file} ' \
                    f'-o {plot_file_raw} -w {w_max} -b -p'
                run_bash(cmd, args.local)

    if gather_only:
        tar = tarfile.open(phyl_dir + '.tar.gz', 'w:gz')
        tar.add(phyl_dir)
        tar.close()


def merge_datasets(vcf_files, args):
    if args.tests == 'both' or args.tests == 'dispersion':
        disp_file = os.path.join(args.out_dir, 'Poisson_dispersion_all.tsv')
    else:
        disp_file = None

    tree_files = []
    if args.tests == 'both' or args.tests == 'tree':
        tree_files.extend(get_poisson_tree_files(vcf_files))

    if disp_file:
        df_disp = pd.read_csv(disp_file, sep='\t', index_col=[0, 1, 2])
        df_disp.drop(['H0', 'H1', 'hypothesis', '-2logLR'], axis=1, inplace=True)
        df_disp.columns = [f'{i}.dispersion' for i in df_disp.columns]

    if tree_files:
        df_trees = pd.DataFrame()
        for tree_file in tree_files:
            if not os.path.exists(tree_file):
                print(f'!WARNING! Missing tree file: {tree_file}')
                continue
            tree = os.path.basename(tree_file).split('_')[2]
            df_tree = pd.read_csv(tree_file, sep='\t', index_col=[0, 1, 2])
            drop_cols = [i for i in df_tree.columns \
                if 'hypothesis' in i or '-2logLR' in i]
            df_tree.drop(drop_cols, axis=1, inplace=True)
            df_tree.columns = [f'{i}.tree.{tree}' for i in df_tree.columns]
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
                    df_trees = df_trees.append(df_tree)

    if disp_file and tree_files:
        df = df_disp.merge(df_trees, left_index=True, right_index=True)
    elif disp_file:
        df = df_disp
    else:
        df = df_trees
    dof_cols = [i for i in df.columns if 'dof' in i]
    df.rename({dof_cols[0]: 'dof'}, axis=1, inplace=True)
    df.drop(dof_cols[1:], axis=1, inplace=True)
    df.sort_index(inplace=True)

    out_file = os.path.join(args.out_dir, 'Summary_biological_data.tsv')
    df.to_csv(out_file, sep='\t')


def compress_results(vcf_files, args):
    comp_dir = os.path.join(args.out_dir,
        f'{datetime.now():%Y%m%d_%H:%M:%S}_compressed')
    print(f'Writing files to: {comp_dir}.tar.gz')
    os.mkdir(comp_dir)

    old_file = os.path.join(args.out_dir, 'Summary_biological_data.tsv')
    new_file = os.path.join(comp_dir, os.path.basename(old_file))
    run_bash(f'cp {old_file} {new_file}', False, '')

    for old_file in get_plot_files(vcf_files):
        path_parts = old_file.split(os.path.sep)
        old_name = path_parts[-1]
        if 'scite_dir' in path_parts:
            tree = 'scite'
            subset = path_parts[-3]
            new_name = old_name \
                .replace('_outg_ml0.newick_w500_mapped.png', '')
        else:
            tree = 'cellphy'
            subset = path_parts[-2]
            new_name = old_name \
                .replace('_outg.vcf.gz.raxml.bestTree_w500_mapped.png', '')

        dataset, filters = new_name.split('.')
        new_name = f'{dataset}_{subset}_{filters}_{tree}.png'
        new_file = os.path.join(comp_dir, new_name)
        if not os.path.exists(old_file):
            print(f'\tMissing file: {old_file}')
            continue
        run_bash(f'cp {old_file} {new_file}', False, '')

    tar = tarfile.open(comp_dir + '.tar.gz', 'w:gz')
    tar.add(comp_dir)
    tar.close()


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
        if 'all' in path_strs[-2]:
            continue

        files.append(vcf_file + '.raxml.bestTree_w500_mapped.png')
        files.append(os.path.join(os.path.dirname(vcf_file), 'scite_dir',
            f'{dataset}.{filters}_ml0.newick_w500_mapped.png'))
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
        choices=['run', 'merge', 'compress', 'plot', 'gather_plot'],
        default='run', help='Which task to do. Default = run.')
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
        run_tests(vcf_files, args)
    elif args.mode == 'merge':
        merge_datasets(vcf_files, args)
    elif args.mode == 'plot' or args.mode == 'gather_plot':
        run_plotting(vcf_files, args, args.mode == 'gather_plot')
    else:
        compress_results(vcf_files, args)