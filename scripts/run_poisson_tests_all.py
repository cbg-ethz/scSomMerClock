#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import subprocess


def run_poisson_disp(vcf_files, exe, out_dir):
    out_file = os.path.join(out_dir, 'Poisson_dispersion_all.tsv')
    cmmd = ' '.join(
        ['python', exe, ' '.join(vcf_files), '-o', out_file, '-b'])
    print('\nShell command:\n{}\n'.format(cmmd))

    poissonDisp = subprocess.Popen(
        cmmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = poissonDisp.communicate()
    poissonDisp.wait()
    stdout = str(stdout)

    if stderr:
        for i in stdout.split('\\n'):
            print(i)
        raise RuntimeError('Error in Poisson Dispersion: {}'.format(stderr))
    return out_file


def run_poisson_tree(tree, vcf_file, args):
    path_strs = vcf_file.split(os.path.sep)
    try:
        clock_dir_no = path_strs.index('ClockTest')
    except ValueError:
        dataset = 'unknown'
        subset = 'unknown'
    else:
        dataset = path_strs[clock_dir_no - 1]
        subset = path_strs[clock_dir_no + 1]

    out_file = os.path.join(args.out_dir, f'Poisson_tree_{tree}_{dataset}_{subset}.tsv')
    if tree == 'cellphy':
        tree_file = vcf_file + '.raxml.bestTree'
    elif tree == 'scite':
        vcf_dir = os.path.dirname(vcf_file)
        tree_file = os.path.join(vcf_dir, 'scite_dir', 'scite_tree_ml0.newick')
    else:
        raise RuntimeError(f'Unknown tree file: {tree}')

    if not os.path.exists(tree_file):
        print(f'!WARNING! Missing {tree: >7} tree file: {tree_file}')
        return

    cmmd = ' '.join(
        ['python', args.exe_tree, vcf_file, tree_file, '-o', out_file,
            '-e', args.exe_paup, '-b'])
    print('\nShell command:\n{}\n'.format(cmmd))

    poissonTree = subprocess.Popen(
        cmmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = poissonTree.communicate()
    poissonTree.wait()
    stdout = str(stdout)

    if not 'success' in stdout and stderr:
        for i in stdout.split('\\n'):
            print(i)
        raise RuntimeError('Error in Poisson Tree: {}'.format(stderr))
    return out_file


def merge_datasets(disp_file, tree_files, out_dir):
    if disp_file:
        df_disp = pd.read_csv(disp_file, sep='\t', index_col=[0, 1])
        df_disp.drop(['H0', 'H1', 'hypothesis'], axis=1, inplace=True)
        df_disp.columns = [f'{i}.dispersion' for i in df_disp.columns]

    if tree_files:
        df_trees = pd.DataFrame()
        for tree_file in tree_files:
            if not tree_file:
                continue
            tree = os.path.basename(tree_file).split('_')[2]
            df_tree = pd.read_csv(tree_file, sep='\t', index_col=[0, 1])
            df_tree.drop(['H0', 'H1', 'hypothesis'], axis=1, inplace=True)
            df_tree.columns = [f'{i}.tree.{tree}' for i in df_tree.columns]
            dataset = df_tree.iloc[0].name
            if dataset in df_trees.index:
                df_trees.loc[dataset, df_tree.columns] = df_tree.iloc[0]
            else:
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

    out_file = os.path.join(out_dir, 'Summary_biological_data.tsv')
    df.to_csv(out_file, sep='\t')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str,  help='vcf master file')
    parser.add_argument('-o', '--out_dir', type=str,
        default='poisson_tests_all', help='Output file.')
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
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    with open(args.input, 'r') as f:
        vcf_files = f.read().strip().split('\n')

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    if args.tests == 'both' or args.tests == 'dispersion':
        disp_file = run_poisson_disp(vcf_files, args.exe_disp, args.out_dir)
    else:
        disp_file = None

    # # <TODO> remove test case
    # disp_file = os.path.join(args.out_dir, 'Poisson_dispersion_all.tsv')
    # tree_files = [os.path.join(args.out_dir, i) for i in os.listdir(args.out_dir) if 'tree' in i]
    # merge_datasets(disp_file, tree_files, args.out_dir)
    # import pdb; pdb.set_trace()
    # # -----------------------

    tree_files = []
    if args.tests == 'both' or args.tests == 'tree':
        for vcf_file in vcf_files:
            tree_files.append(run_poisson_tree('cellphy', vcf_file, args))
            tree_files.append(run_poisson_tree('scite', vcf_file, args))

    merge_datasets(disp_file, tree_files, args.out_dir)


