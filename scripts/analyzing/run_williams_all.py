#!/usr/bin/env python3

import argparse
from datetime import datetime
import os
import re
import shutil
import subprocess
import tarfile

import numpy as np
import pandas as pd


MODULE_STR = ''

DEPTH = {
    'CRC08_ND_WGS': {'CRC08_ND_WGS': 41.7},
    'CRC08_NP_WGS': {'CRC08_NP_WGS': 42.8},
    'CRC08_TD_WGS': {'CRC08_TD_WGS': 46.5},
    'CRC09_TI_WGS_MG': {'CRC09_TI_WGS_MG': 42.9},
    'CRC09_TM_WGS_MG': {'CRC09_TM_WGS_MG': 39.8},
    'Li55': {'BGI_BC-T': 176.1, 'BGI_BN-T': 24.3},
    'Ni8_P01M01E': {'SRR975206': 35.5, 'SRR975210': 51.0},
    'Ni8_P01P01E': {'SRR975206': 35.5, 'SRR975212': 65.9},
    'W32': {'SRR1163508': 86.5, 'SRR1298936': 59.1},
    'W55': {'SRR1153400': 71.6, 'SRR1153401': 18.1},
    'Wu61_Cancer': {'SRR3086496': 57.7, 'SRR3086497': 51.4},
    'Wu61_Polyps': {'SRR3086496': 57.7, 'SRR3086498': 48.3},
    'Wu63_CRC0827-Ca-1': {'CRC0827-Ca-1': 25, 'CRC0827-Normal': 53.6},
    'Wu63_CRC0827-Ca-2': {'CRC0827-Ca-2': 29.9, 'CRC0827-Normal': 53.6},
    'Wu63_Polyps': {'CRC0827-Adenoma_Polyps': 58.8, 'CRC0827-Normal': 53.6},
    'X25': {'BGI_RC-T': 126.8, 'BGI_RN-T': 37.3},
    'H65_BGI_LC-T1': {'BGI_LC-T1': 54.2, 'BGI_YH-Control': 24.1},
    'H65_BGI_LN-T1': {'BGI_LN-T1': 35.9, 'BGI_YH-Control': 24.1}
}

CELLULARITY = {
    'CRC08_NP_WGS': {'CRC08_NP_WGS': 0.99},
    'CRC08_TD_WGS': {'CRC08_TD_WGS': 0.99},
    'CRC09_TI_WGS_MG': {'CRC09_TI_WGS_MG': 0.92},
    'CRC09_TM_WGS_MG': {'CRC09_TM_WGS_MG': 0.87},
    'Li55': {'BGI_BC-T': 0.91},
    'Ni8_P01M01E': {'SRR975210': 0.98},
    'Ni8_P01P01E': {'SRR975212': 0.19},
    'W32': {'SRR1298936': 0.97},
    'W55': {'SRR1153400': 0.34},
    'Wu61_Cancer': {'SRR3086497': 0.33},
    'Wu61_Polyps': {'SRR3086498': 0.28},
    'Wu63_CRC0827-Ca-1': {'CRC0827-Ca-1': 0.1},
    'Wu63_CRC0827-Ca-2': {'CRC0827-Ca-2': 0.1},
    'Wu63_Polyps': {'CRC0827-Adenoma_Polyps': 0.44},
    'X25': {'BGI_RC-T': 0.4},
    'H65_BGI_LC-T1': {'BGI_LC-T1': 1.0},
    'H65_BGI_LN-T1': {'BGI_LN-T1': 0.99}
}


def run_bash(cmd_raw, bsub=True, module_str=MODULE_STR):
    if bsub:
        cmd = f"sbatch -t 60 -p amd-shared --qos amd-shared --mem 16G " \
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


def convert_vcf(vcf_file, args):
    basename = os.path.splitext(os.path.basename(vcf_file))[0]
    out_file = os.path.join(args.out_dir, basename)

    if not args.replace:
        existing = []
        for file in os.listdir(args.out_dir):
            if file.startswith(basename) and file.endswith('.VAF'):
                existing.append(os.path.join(args.out_dir, file))
        if len(existing) > 0:
            print('Existing files:\n\t' + '\n\t'.join(existing))
            return existing
    cmd = f'python {args.exe_VAF} {vcf_file} -o {out_file} -b'
    run_bash(cmd, False)

    VAFs = []
    for file in os.listdir(args.out_dir):
        if file.startswith(basename) and file.endswith('.VAF'):
            VAFs.append(os.path.join(args.out_dir, file))
    return VAFs


def run_williams(VAF_file, args):
    basename = os.path.splitext(os.path.basename(VAF_file))[0]
    out_file = os.path.join(args.out_dir, basename + '.mobster')

    if not args.replace and os.path.exists(out_file):
        print(f'Existing files:\n\t{out_file}')
        return

    dataset = basename.split('.')[0]
    sample = re.search('PASS_(.*)\.VAF', VAF_file).group(1)
    if 'tetra' in VAF_file:
        ploidy = 4
    else:
        ploidy = 2

    try:
        depth = DEPTH[dataset][sample]
        cellularity = CELLULARITY[dataset][sample]
    except KeyError:
        raise KeyError(f'Unknown key {sample} in: {dataset}\t({VAF_file})')

    # From https://cran.r-project.org/web/packages/neutralitytestr/vignettes/neutraltytestr.html
    fmax = cellularity / ploidy - 2 * np.sqrt(1 / depth)
    if fmax < 0.1:
        fmax = 0.25
    print(f'Neutralitytest fmax: {fmax:.4f}')

    cmd = f'{args.exe_will} {VAF_file} {out_file} --fmax {fmax} -K 3 --plot'
    run_bash(cmd, args.local)


def merge_datasets(args):
    cols = ['dataset', 'R^2_pVal', 'area_pVal', 's_Bayes', 'clones_Bayes']
    df = pd.DataFrame([], columns=cols)

    i = 0
    for in_file in os.listdir(args.input):
        if not in_file.endswith('.mobster'):
            continue

        df.loc[i, 'dataset'] = os.path.basename(in_file).split('.')[0]

        with open(os.path.join(args.input, in_file), 'r') as f:
            log_lines = f.read().strip().split('\n')[1:]

            for j, line_raw in enumerate(log_lines):
                line = line_raw.strip()
                if line.startswith('mu'):
                    header = re.split('\s+', line)
                    if len(header) == 2: # only mu & exponent: no clone
                        df.loc[i, 's_Bayes'] = 0
                        df.loc[i, 'clones_Bayes'] = 0
                    elif 's' in header:
                        s_idx = header.index('s') + 1
                        s1 = float(re.split('\s+', log_lines[j + 1])[s_idx])
                        if log_lines[j + 2].strip().startswith('2'):
                            freq_idx = header.index('subclonefrequency') + 1
                            freq1 = float(
                                re.split('\s+', log_lines[j + 1])[freq_idx])
                            freq2 = float(
                                re.split('\s+', log_lines[j + 2])[freq_idx])
                            s2 = float(re.split('\s+', log_lines[j + 2])[s_idx])
                            df.loc[i, 's_Bayes'] = (s1 * freq1 + s2 * freq2) \
                                / (freq1 + freq2)
                            df.loc[i, 'clones_Bayes'] = 2
                        else:
                            df.loc[i, 's_Bayes'] = s1
                            df.loc[i, 'clones_Bayes'] = 1

                elif line.startswith('s'):
                    if line.endswith('subclone'):
                        freq1 = float(re.split('\s+', log_lines[j - 2])[5])
                        freq2 = float(re.split('\s+', log_lines[j - 1])[5])
                        s1 = float(re.split('\s+', log_lines[j + 1])[1])
                        s2 = float(re.split('\s+', log_lines[j + 2])[1])
                        # weighted mean
                        df.loc[i, 's_Bayes'] = (s1 * freq1 + s2 * freq2) \
                            / (freq1 + freq2)
                        df.loc[i, 'clones_Bayes'] = 2
                    else:
                        df.loc[i, 's_Bayes'] = float(
                            log_lines[j + 1].strip().split(' ')[1])
                        df.loc[i, 'clones_Bayes'] = 1
                # Best statistic, accrording to Williams et al. 2018, p. 11
                elif line == 'Area:':
                    df.loc[i, 'area_pVal'] = float(
                        log_lines[j + 1].strip().split(' ')[-1])
                elif line == 'R^2:':
                    df.loc[i, 'R^2_pVal'] = float(
                        log_lines[j + 1].strip().split(' ')[-1])
                    # Last relevant line
                    break
        i += 1

    out_file = os.path.join(args.out_dir, 'Summary_williams.tsv')
    df.sort_values('dataset').to_csv(out_file, sep='\t', index=False)



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str,
        default='/mnt/lustre/scratch/home/uvi/be/mva/singlecell/Projects/mol_clock/VC_diploid/',
        help='Input dir to files with pattern: <DATASET>.mutect2filtered.PASS.vcf')
    parser.add_argument('-o', '--out_dir', type=str,
        default='poisson_tests_all', help='Output file.')
    parser.add_argument('--merge', action='store_true',
        help='Merge individual clock test to final table.')
    parser.add_argument('-evaf', '--exe_VAF', type=str,
        default='simulations/scripts/VCF_to_VAF.py',
        help='VCF to VAF converter exe.')
    parser.add_argument('-ewill', '--exe_will', type=str,
        default='simulations/scripts/run_williams2016.R',
        help='Williams 2016 R script exe.')
    parser.add_argument('-l', '--local', action='store_false',
        help='Run locally instead of HPC.')
    parser.add_argument('-r', '--replace', action='store_true',
        help='Overwrite already existing files.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    if args.merge:
        merge_datasets(args)
    else:
        for in_file in sorted(os.listdir(args.input)):
            if not in_file.endswith('.mutect2filtered.PASS.vcf'):
                continue
            if in_file == 'W32.mutect2filtered.PASS.vcf':
                in_file = 'W32.mutect2filtered.tetra.PASS.vcf'
            VAF_files = convert_vcf(os.path.join(args.input, in_file), args)
            for VAF_file in VAF_files:
                run_williams(VAF_file, args)




