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
DRIVER_FILE = '../data/resources/2020-02-02_IntOGen-Drivers-20200213/Compendium_Cancer_Genes.tsv'
DATA_DIRS = ['H65', 'Li55', 'Lo-P1', 'Lo-P2', 'Lo-P3', 'Ni8', 'S21_P1',
    'S21_P2', 'W32', 'W55', 'Wu61', 'Wu63', 'X25']


def run_bash(cmd_raw, bsub=True, module_str=MODULE_STR):
    if bsub:
        cmd = f"sbatch -t 30 -p amd-shared --qos amd-shared --mem 2G " \
            f"--wrap '{module_str} {cmd_raw}'"
    else:
        cmd = f'{cmd_raw}'

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
    if 'dispersion' in args.tests:
        run_poisson_disp(vcf_files, args)
    if 'poissonTree' in args.tests:
        run_poissonTree(vcf_files, args)


def run_poisson_disp(vcf_files, args):
    out_file = os.path.join(args.out_dir, 'Poisson_dispersion_all.tsv')
    if not os.path.exists(out_file) and args.check:
        print(f'!Missing! Poisson Dispersion file: {out_file}')
        return

    if os.path.exists(out_file) and not args.replace:
        return
    cmd = f'python {args.exe_disp} {" ".join(vcf_files)} -o {out_file} -b'
    run_bash(cmd, args.local, '')


def run_poissonTree(vcf_files, args):
    for vcf_file in vcf_files:
        for tree in ['cellphy', 'scite']:
            run_poissonTree_single(vcf_file, tree, args)


def run_poissonTree_single(vcf_file, tree, args):
    file_ids = os.path.basename(vcf_file).replace('.vcf.gz', '').split('_')

    dataset = file_ids[0]

    if dataset == 'S21':
        dataset = '_'.join(file_ids[:2])
        subset = file_ids[2]
        filters = file_ids[3]
        if filters == 'all':
            return
    else:
        subset = file_ids[1]
        filters = file_ids[2]

    if dataset == 'Wu61':
        dataset = file_ids[0]
        subset = '_'.join(file_ids[1:-1])
        filters = file_ids[-1]

    if dataset not in args.dataset:
        return

    out_base = f'{dataset}_{subset}_{filters}_{tree}'
    tree_file = os.path.join(args.input,
        f'{dataset}_{subset}_{filters}.{tree}.newick')

    # Run poisson tree test
    out_file = os.path.join(args.out_dir, f'{out_base}.poissonTree.tsv')
    if os.path.exists(out_file) and not args.replace:
        pass
    else:
        if args.check:
            print(f'!Missing! Poisson tree file: {out_file}')
            return

        if not os.path.exists(tree_file):
            print(f'!WARNING! Missing file: {tree_file}')
            return
        cmd = f'python {args.exe_tree} {vcf_file} {tree_file} -o {out_file} -b'
        run_bash(cmd, args.local)

    # Plot phylogenic tree
    fig_file = os.path.join(args.out_dir,
        f'{out_base}.w{args.plotting_wmax:.0f}_mapped.pdf')
    if os.path.exists(fig_file) and not args.replace:
        pass
    else:
        if args.check:
            print(f'!Missing! Phylogenetic tree file: {fig_file}')
            return

        if not os.path.exists(tree_file):
            print(f'!WARNING! Missing file: {tree_file}')
            return
        prefix = os.path.join(args.out_dir, out_base)
        cmd_plt = f'python {args.exe_tree} {vcf_file} {tree_file} -o {prefix} ' \
            f'-b --plotting --w_max {args.plotting_wmax}'
        if args.drivers:
            cmd_plt += f' --drivers {args.drivers}'
        run_bash(cmd_plt, args.local)

    # Generate mutation-branch mapping file
    mapping_dir = os.path.join(args.out_dir, 'Mappings')
    if not os.path.exists(mapping_dir):
        os.mkdir(mapping_dir)
    map_file = os.path.join(mapping_dir,
        f'{out_base}.w{args.plotting_wmax:.0f}_branchMuts.tsv')
    if os.path.exists(map_file) and not args.replace:
        pass
    else:
        if args.check:
            print(f'!Missing! Branch-mutation mapping file: {map_file}')
            return

        if not os.path.exists(tree_file):
            print(f'!WARNING! Missing file: {tree_file}')
            return
        prefix = os.path.join(args.out_dir, out_base)
        cmd_map = f'python {args.exe_tree} {vcf_file} {tree_file} -b ' \
            f'--w_max {args.plotting_wmax} -mbm {map_file}'
        run_bash(cmd_map, args.local)



def merge_datasets(vcf_files, args):
    if 'dispersion' in args.tests:
        disp_file = os.path.join(args.out_dir, 'Poisson_dispersion_all.tsv')
    else:
        disp_file = None

    tree_files = []
    if 'poissonTree' in args.tests:
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
        if 'all' in path_strs[-2] and not dataset.startswith(('Lo-P', 'S21_P1')):
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
    parser.add_argument('input', type=str,  help='Input directory')
    parser.add_argument('-o', '--out_dir', type=str, default='',
        help='Output directory. Default = <INPUT>/poissonTests_all')
    parser.add_argument('-m', '--mode', type=str, default='run',
        choices=['run', 'merge', 'compress'],
        help='Which task to do. Default = run.')
    parser.add_argument('-et', '--exe_tree', type=str,
        default='simulations/scripts/get_poisson_tree_LRT.py',
        help='Poisson Tree exe.')
    parser.add_argument('-ed', '--exe_disp', type=str,
        default='simulations/scripts/get_poisson_LRT.py',
        help='Poisson Dispersion exe.')
    parser.add_argument('-dr', '--drivers', type=str,
        default=DRIVER_FILE,
        help=f'Path to IntOGen driver file. Default = {DRIVER_FILE}.')
    parser.add_argument('-plt_w', '--plotting_wmax', type=int, default=400,
        help='W_max value used for coloring braches in Phylogentic tree. ' \
            'Default = 400.')
    parser.add_argument('-t', '--tests', default=['poissonTree', 'dispersion'],
        choices=['poissonTree', 'dispersion'], help='Tests to perform.')
    parser.add_argument('-da', '--dataset', type=str, nargs='+',
        choices=DATA_DIRS, default=DATA_DIRS,
        help='Datasets to process. Default = all.')
    parser.add_argument('-l', '--local', action='store_false',
        help='Run locally instead of HPC.')
    parser.add_argument('-r', '--replace', action='store_true',
        help='Overwrite already existing files.')
    parser.add_argument('-c', '--check', action='store_true',
        help='Check only if files exist, do not run anything.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    vcf_files = []
    for file in sorted(os.listdir(args.input)):
        if not file.endswith(('.vcf.gz', '.vcf')):
            continue
        vcf_files.append(os.path.join(args.input, file))

    if not args.out_dir:
        args.out_dir = os.path.join(args.input, 'poissonTests_all')

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    if args.mode == 'run':
        run_tests(vcf_files, args)
    elif args.mode == 'merge':
        merge_datasets(vcf_files, args)
    else:
        compress_results(vcf_files, args)