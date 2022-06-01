#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess


DATASETS = ['CRC08', 'CRC09', 'H65', 'Li55', 'Lo-P1', 'Lo-P2', 'Lo-P3', 'Ni8', 'S21_P1',
    'S21_P2', 'W32', 'W55', 'Wu61', 'Wu63', 'X25']
DATA_FILTERS = ['all', '33nanFilter', '50nanFilter', '99nanFilter']

MODULE_STR = 'module load bcftools; '
SLURM = False


def run_bash(cmd_raw, bsub=True, module_str=MODULE_STR):
    if bsub:
        if SLURM:
            cmd = f"sbatch -n 2 -t 90 --mem 4G -p amd-shared --qos amd-shared " \
                f"--wrap '{module_str} {cmd_raw}'"
        else:
            cmd = f'bsub -n 2 -W 90 -R "rusage[mem=4096]" "{cmd_raw}"'
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
    return str(stdout), str(stderr)


def run_sigProfiler(args):
    all_sigs = True

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

        filters = file_name_raw.split('_')[-1]
        if filters not in args.filters:
            continue

        out_dir = os.path.join(args.out_dir, file_name_raw + '.mutSigs')
        final_file = os.path.join(out_dir, 'SBS96', 'Suggested_Solution',
                'COSMIC_SBS96_Decomposed_Solution',
                'De_Novo_map_to_COSMIC_SBS96.csv')
        if args.check:
            if not os.path.exists(final_file):
                print(f'!Missing! Signature file: {final_file}')
                all_sigs = False
            continue

        if os.path.exists(out_dir):
            if args.replace:
                shutil.rmtree(out_dir)
            else:
                if os.path.exists(final_file):
                    continue
                else:
                    print(f'Signature directory BUT NOT signature present: {out_dir}')
                    if len(os.listdir(out_dir)) > 1:
                        print('\t and continues more than vcf file: Skipping')
                        continue

        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        vcf_file_uncomp = os.path.join(out_dir, file_name_raw + '.vcf')
        cell_names_raw, _ = run_bash(f'bcftools query -l {in_file}', False)
        cell_names = cell_names_raw[2:-3].split('\\n')
        try:
            cell_names.remove('healthycell')
        except ValueError:
            pass

        for cell_name in cell_names:
            uncomp_cmd = f'bcftools view -s {cell_name} {in_file} -O v ' \
                f'-o {out_dir}/{file_name_raw}.{cell_name}.vcf'
            run_bash(uncomp_cmd, False)

        out_dir_temp = os.path.join(args.out_dir, file_name_raw + '.temp')
        if os.path.exists(out_dir_temp):
            shutil.rmtree(out_dir_temp)
        os.mkdir(out_dir_temp)
        cmd = f'python {args.exe_profiler} {out_dir} -o {out_dir_temp} -n 2 ' \
            f'&& mv {out_dir_temp}/SBS96 {out_dir} && rm -r {out_dir_temp}'
        run_bash(cmd, args.local)

    if args.check and all_sigs:
        print('All signature files exist!')


def merge_sigProfiler(args):
    str_out = 'dataset\tsubset\tfilters\tsignatures\n'
    for dataset_dir in sorted(os.listdir(args.input)):
        sig_file = os.path.join(args.input, dataset_dir, 'SBS96',
            'Suggested_Solution', 'COSMIC_SBS96_Decomposed_Solution',
            'De_Novo_map_to_COSMIC_SBS96.csv')
        if not os.path.exists(sig_file):
            print(f'\n!WARNING! Missing file: {sig_file}')
            continue

        file_ids = dataset_dir.replace('.mutSigs', '').split('_')
        dataset = '_'.join(file_ids[:-2])
        subset = file_ids[-2]
        filters = file_ids[-1]

        with open(sig_file, 'r') as f:
            sig_raw = f.read()
        sig = sig_raw.split('\n')[1].split(',')[1]
        str_out += f'{dataset}\t{subset}\t{filters}\t{sig}\n'

    with open(args.merge, 'w') as f:
        f.write(str_out)



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
        help=f'Datasets to process. Default = {DATASETS}.')
    parser.add_argument('-fl', '--filters', type=str, nargs='+',
        choices=DATA_FILTERS, default=DATA_FILTERS,
        help=f'Datafilters to process. Default = {DATA_FILTERS}.')
    parser.add_argument('-l', '--local', action='store_false',
        help='Run locally instead of HPC.')
    parser.add_argument('-r', '--replace', action='store_true',
        help='Overwrite already existing files.')
    parser.add_argument('-c', '--check', action='store_true',
        help='Check only if files exist, do not run anything.')
    parser.add_argument('--merge', type=str, default='',
        help='Merge individual signature to summary file.')
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

    if args.merge:
        merge_sigProfiler(args)
    else:
        run_sigProfiler(args)