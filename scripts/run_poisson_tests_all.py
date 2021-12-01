#!/usr/bin/env python3

import argparse
import os
import subprocess


def run_poisson_disp(vcf_files, exe, out_dir):
    out_file = os.path.join(out_dir, 'Poisson_dispersion_all.tsv')
    cmmd = ' '.join(
        ['python', exe, '-i', ' '.join(vcf_files), '-o', out_file, '-b'])
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
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    with open(args.input, 'r') as f:
        vcf_files = f.read().strip().split('\n')

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)
    run_poisson_disp(vcf_files, args.exe_disp, args.out_dir)
    exit()
    for vcf_file in vcf_files:

        import pdb; pdb.set_trace()
        run_poisson_tree_cellphy(vcf_file, args.exe_tree)
        run_poisson_tree_scite(vcf_file, args.exe_tree)

