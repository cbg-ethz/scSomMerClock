#!/usr/bin/env python3

import argparse
import os
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


def run_poisson_tree(tree, vcf_file, exe, paup_exe, out_dir):
    path_strs = vcf_file.split(os.path.sep)
    try:
        clock_dir_no = path_strs.index('ClockTest')
    except ValueError:
        dataset = 'unknown'
        subset = 'unknown'
    else:
        dataset = path_strs[clock_dir_no - 1]
        subset = path_strs[clock_dir_no + 1]

    out_file = os.path.join(out_dir, f'Poisson_tree_{tree}_{dataset}_{subset}.tsv')
    if tree == 'cellphy':
        tree_file = vcf_file + '.raxml.bestTree'
    elif tree == 'scite':
        vcf_dir = os.path.dirname(vcf_file)
        tree_file = os.path.join(vcf_dir, 'scite_dir', 'scite_tree_ml0.newick')
    else:
        raise RuntimeError(f'Unknown tree file: {tree}')

    cmmd = ' '.join(
        ['python', exe, vcf_file, tree_file, '-o', out_file, '-e', paup_exe, '-b'])
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
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    with open(args.input, 'r') as f:
        vcf_files = f.read().strip().split('\n')

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)
    # <TODO> uncomment
    run_poisson_disp(vcf_files, args.exe_disp, args.out_dir)

    # # <TODO> remove test case
    # test_vcf = '../data/Ni9/Ni8_cancer.vcf.gz'
    # run_poisson_tree('cellphy', test_vcf, args.exe_tree, args.exe_paup, args.out_dir)
    # run_poisson_tree('scite', test_vcf, args.exe_tree, args.exe_paup, args.out_dir)
    # -----------------------

    for vcf_file in vcf_files:
        run_poisson_tree_cellphy(vcf_file, args.exe_tree, args.exe_paup, args.out_dir)
        run_poisson_tree_scite(vcf_file, args.exe_tree, args.exe_paup, args.out_dir)

