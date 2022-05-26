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
DATA_DIRS = ['CRC08', 'CRC09', 'H65', 'Li55', 'Lo-P1', 'Lo-P2', 'Lo-P3', 'Ni8', 'S21_P1',
    'S21_P2', 'W32', 'W55', 'Wu61', 'Wu63', 'X25']
DATA_FILTERS = ['all', '33nanFilter', '50nanFilter', '99nanFilter']

SLURM = True


def run_bash(cmd_raw, bsub=True, module_str=MODULE_STR):
    if bsub:
        if SLURM:
            cmd = f"sbatch -t 30 --mem 2G -p amd-shared --qos amd-shared " \
                f"--wrap '{module_str} {cmd_raw}'"
        else:
            cmd = f'bsub -W 30 -R "rusage[mem=2048]" "{cmd_raw}"'
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
        for tree in args.method:
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

    if dataset not in args.dataset or filters not in args.filters:
        return

    out_base = f'{dataset}_{subset}_{filters}_{tree}'
    tree_file = os.path.join(args.input,
        f'{dataset}_{subset}_{filters}.{tree}.newick')

    # Run poisson tree test
    if 'poissonTree' in args.mode:
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

    if 'phylogeny' in args.mode:
        # Plot phylogenic tree
        fig_file = os.path.join(args.out_dir,
            f'{out_base}.w{args.plotting_wmax:.0f}.png')
        if os.path.exists(fig_file) and not args.replace:
            pass
        else:
            if args.check:
                print(f'!Missing! Phylogenetic tree file: {fig_file}')
                return

            if not os.path.exists(tree_file):
                print(f'!WARNING! Missing file: {tree_file}')
                return
            cmd_plt = f'python {args.exe_tree} {vcf_file} {tree_file} -o {fig_file} ' \
                f'-b --plotting --w_max {args.plotting_wmax}'
            if args.drivers:
                cmd_plt += f' --drivers {args.drivers}'
            run_bash(cmd_plt, args.local)

    if 'mutMapping' in args.mode:
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


def merge_datasets(tsv_files, args):
    for file in tsv_files:
        file_name = os.path.basename(file)

        if file_name.startswith('Poisson_dispersion'):
            continue
        elif file_name.startswith('S21_P1'):
            dataset = 'S21_P1'
            _, _, subset, filter_str, tree_raw = file_name.split('_')
        elif file_name.startswith('S21_P2'):
            dataset = 'S21_P2'
            _, _, subset, filter_str, tree_raw = file_name.split('_')
        elif file_name.startswith('Wu61'):
            dataset = 'Wu61'
            tree_raw = file_name.split('_')[-1]
            filter_str = file_name.split('_')[-2]
            subset = '_'.join(file_name.split('_')[1:-2])
        else:
            dataset, subset, filter_str, tree_raw = file_name.split('_')
        tree = tree_raw.split('.')[0]
        if tree not in args.method:
            continue
        df_new = pd.read_csv(file, sep='\t')
        df_new.loc[0, 'subset'] = subset
        df_new.loc[0, 'filters'] = filter_str

        df_new.insert(3, 'tree', tree)

        drop_cols = [i for i in df_new.columns \
            if 'hypothesis' in i or '-2logLR' in i]
        df_new.drop(drop_cols, axis=1, inplace=True)

        try:
            df = pd.concat([df, df_new], axis=0, ignore_index=True)
        except NameError:
            df = df_new

    dof_cols = [i for i in df.columns if 'dof' in i]
    df.rename({dof_cols[0]: 'dof'}, axis=1, inplace=True)
    df.drop(dof_cols[1:], axis=1, inplace=True)
    df.sort_index(inplace=True)

    out_file = os.path.join(args.out_dir, '00_summary_biological_data.tsv')
    df.to_csv(out_file, sep='\t', index=False)


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


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str,  help='Input directory')
    parser.add_argument('-o', '--out_dir', type=str, default='',
        help='Output directory. Default = <INPUT>/poissonTests_all')
    parser.add_argument('-m', '--mode', nargs='+', type=str,
        choices=['poissonTree', 'phylogeny', 'mutMapping'],
        default=['poissonTree', 'phylogeny', 'mutMapping'],
        help='Which task to do. Default = run.')
    parser.add_argument('-et', '--exe_tree', type=str,
        default='simulations/scripts/get_poisson_tree_LRT.py',
        help='Poisson Tree exe.')
    parser.add_argument('-ed', '--exe_disp', type=str,
        default='simulations/scripts/get_poisson_LRT.py',
        help='Poisson Dispersion exe.')
    parser.add_argument('-dr', '--drivers', type=str, default=DRIVER_FILE,
        help=f'Path to IntOGen driver file. Default = {DRIVER_FILE}.')
    parser.add_argument('-plt_w', '--plotting_wmax', type=int, default=1000,
        help='W_max value used for coloring braches in Phylogentic tree. ' \
            'Default = 1000.')
    parser.add_argument('-t', '--tests', default=['poissonTree'],
        choices=['poissonTree', 'dispersion'], help='Tests to perform.')
    parser.add_argument('-da', '--dataset', type=str, nargs='+',
        choices=DATA_DIRS, default=DATA_DIRS,
        help='Datasets to process. Default = all.')
    parser.add_argument('-fl', '--filters', type=str, nargs='+',
        choices=DATA_FILTERS, default=DATA_FILTERS ,
        help=f'Datafilters to process. Default = {DATA_FILTERS}.')
    parser.add_argument('-me', '--method', nargs='+', type=str,
        choices=['cellphy', 'scite'], default=['cellphy'],
        help=f'Tree inference method. Default = cellphy.')
    parser.add_argument('-l', '--local', action='store_false',
        help='Run locally instead of HPC.')
    parser.add_argument('-r', '--replace', action='store_true',
        help='Overwrite already existing files.')
    parser.add_argument('-c', '--check', action='store_true',
        help='Check only if files exist, do not run anything.')
    parser.add_argument('--merge', action='store_true',
        help='Merge individual clock test to final table.')
    parser.add_argument('--lsf', action='store_true',
        help='Run lsf queue submit command (bsub) instead of slurm (sbatch).')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    if args.lsf:
        SLURM = False
        MODULE_STR = 'conda activate ete3;'

    vcf_files = []
    tsv_files = []
    for file in sorted(os.listdir(args.input)):
        if file.endswith(('.vcf.gz', '.vcf')):
            vcf_files.append(os.path.join(args.input, file))
        elif file.endswith('tsv') and not 'summary' in file:
            tsv_files.append(os.path.join(args.input, file))

    if args.merge:
        if not tsv_files:
            raise IOError(f'\nNo tsv files found in: {args.input}\n')
        if not args.out_dir:
            args.out_dir = args.input
        merge_datasets(tsv_files, args)
    else:
        if not vcf_files:
            raise IOError(f'\nNo vcf files found in: {args.input}\n')

        if not args.out_dir:
                args.out_dir = os.path.join(args.input, 'poissonTests_all')
        if not os.path.exists(args.out_dir):
            os.mkdir(args.out_dir)

        run_tests(vcf_files, args)