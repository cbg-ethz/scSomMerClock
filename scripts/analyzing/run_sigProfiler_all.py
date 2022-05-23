#!/usr/bin/env python3

import argparse
import os
import subprocess


DATASETS = ['CRC08', 'CRC09', 'H65', 'Li55', 'Lo-P1', 'Lo-P2', 'Lo-P3', 'Ni8', 'S21_P1',
    'S21_P2', 'W32', 'W55', 'Wu61', 'Wu63', 'X25']

MODULE_STR = ''
SLURM = False


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


def run_sigProfiler(vcf_files, args):
    for vcf_file in vcf_files:
        file_ids = os.path.basename(vcf_file).replace('.vcf.gz', '').split('_')
        dataset = file_ids[0]
        if dataset == 'S21':
            dataset = '_'.join(file_ids[:2])

        if dataset not in args.dataset:
            return

        out_dir = vcf_file.replace('.vcf.gz', '') + '_signatures'
        out_file = f'{vcf_file}'
        if os.path.exists(out_file) and not args.replace:
            pass
        else:
            if args.check:
                print(f'!Missing! Signature file: {out_file}')
                return

            cmd = f'python {args.exe_profiler} {vcf_file} -o {out_file}'
            run_bash(cmd, args.local)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str,  help='Input directory')
    parser.add_argument('-o', '--out_dir', type=str, default='',
        help='Output directory. Default = <INPUT>/poissonTests_all')
    parser.add_argument('-ep', '--exe_profiler', type=str,
        default='./get_signatures.py', help='Signature Profiler exe.')
    parser.add_argument('-da', '--dataset', type=str, nargs='+',
        choices=DATASETS, default=DATASETS,
        help='Datasets to process. Default = all.')
    parser.add_argument('-l', '--local', action='store_false',
        help='Run locally instead of HPC.')
    parser.add_argument('-r', '--replace', action='store_true',
        help='Overwrite already existing files.')
    parser.add_argument('-c', '--check', action='store_true',
        help='Check only if files exist, do not run anything.')
    parser.add_argument('--slurm', action='store_true',
        help='Run slurm (sbatch) submit command instead of lsf queue (bsub).')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    if args.slurm:
        SLURM = True

    vcf_files = []
    for file in sorted(os.listdir(args.input)):
        if file.endswith(('.vcf.gz', '.vcf')):
            vcf_files.append(os.path.join(args.input, file))

    if not vcf_files:
        raise IOError(f'\nNo vcf files found in: {args.input}\n')

    if args.out_dir and not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    run_sigProfiler(vcf_files, args)