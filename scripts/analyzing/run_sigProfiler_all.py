#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess


DATASETS = ['CRC08', 'CRC09', 'H65', 'Li55', 'Lo-P1', 'Lo-P2', 'Lo-P3', 'Ni8', 'S21_P1',
    'S21_P2', 'W32', 'W55', 'Wu61', 'Wu63', 'X25']

MODULE_STR = 'module load bcftools; '
SLURM = False


def run_bash(cmd_raw, bsub=True, module_str=MODULE_STR):
    if bsub:
        if SLURM:
            cmd = f"sbatch -n 2 -t 60 --mem 4G -p amd-shared --qos amd-shared " \
                f"--wrap '{module_str} {cmd_raw}'"
        else:
            cmd = f'bsub -n 2 -W 60 -R "rusage[mem=4096]" "{cmd_raw}"'
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


def run_sigProfiler(args):
    for file_name in sorted(os.listdir(args.input)):
        if not file_name.endswith(('.vcf.gz', '.vcf')):
            continue
        in_file = os.path.join(args.input, file_name)
        file_name_raw = file_name.replace('.vcf.gz', '').replace('.vcf', '')

        dataset = file_name_raw.split('_')[0]
        if dataset == 'S21':
            dataset = '_'.join(file_name_raw.split('_')[:2])

        if dataset not in args.dataset:
            continue

        out_dir = os.path.join(args.out_dir, file_name_raw + '.mutSigs')
        final_file = os.path.join(out_dir, 'SBS96', 'Suggested_Solution',
                'COSMIC_SBS96_Decomposed_Solution',
                'De_Novo_map_to_COSMIC_SBS96.csv')
        if args.check:
            if not os.path.exists(final_file):
                print(f'!Missing! Signature file: {final_file}')
            continue

        if os.path.exists(out_dir):
            if args.replace:
                shutil.rmtree(out_dir)
            else:
                print(f'Signature directory already present: {out_dir}')
                if len(os.listdir(out_dir)) > 1:
                    print('\t and continues more than vcf file: Skipping')
                    continue

        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        vcf_file_uncomp = os.path.join(out_dir, file_name_raw + '.vcf')
        if file_name.endswith('.vcf'):
            shutil.copy(in_file, vcf_file_uncomp)
        else:
            # Unzip vcf file
            uncomp_cmd = f'bcftools view {in_file} -O v -o {vcf_file_uncomp}'
            run_bash(uncomp_cmd, False)

        out_dir_temp = os.path.join(args.out_dir, file_name_raw + '.temp')
        if os.path.exists(out_dir_temp):
            shutil.rmtree(out_dir_temp)
        os.mkdir(out_dir_temp)
        cmd = f'python {args.exe_profiler} {out_dir} -o {out_dir_temp} -n 2 ' \
            f'&& mv {out_dir_temp}/SBS96 {out_dir} && rm -r {out_dir_temp}'
        run_bash(cmd, args.local)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str,  help='Input directory')
    parser.add_argument('-o', '--out_dir', type=str, default='',
        help='Output directory. Default = <INPUT>/sigProfiles')
    parser.add_argument('-ep', '--exe_profiler', type=str,
        default='scripts/analyzing/run_sigProfiler.py',
        help='Signature Profiler exe.')
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

    if not args.out_dir:
        args.out_dir = os.path.join(args.input, 'sigProfiles')

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    run_sigProfiler(args)