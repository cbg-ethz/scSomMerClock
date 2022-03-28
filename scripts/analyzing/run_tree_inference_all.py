#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess


# monica_dir = '/mnt/lustre/scratch/home/uvi/be/mva/singlecell/Projects/mol_clock/VariantCallsApril/dec21/'
monica_dir = '/mnt/lustre/scratch/home/uvi/be/mva/singlecell/Projects/mol_clock/VariantCallsApril/filter2'
base_dir = '/home/uvi/be/nbo/data/data/'

data_dirs = {
    'H65_Monica': ['all', 'cancer'],
    'Li55': ['all', 'cancer', 'normal'],
    'Lo-P1': ['all'],
    'Lo-P2': ['all'],
    'Lo-P3': ['all'],
    'Ni8_Monica': ['all', 'cancer'],
    'S21_P1': ['all'],
    'S21_P2': ['all', 'cancer'],
    'W32_Monica': ['all', 'aneuploid', 'haploid', 'normal'],
    'W55': ['all', 'cancer', 'normal'],
    'Wu61': ['all', 'cancer_C', 'cancer_CA', 'normal', 'polyps'],
    'Wu63': ['all', 'cancer', 'normal', 'polyps'],
    'X25': ['all', 'cancer']
}
data_filters = ['all', '33nanFilter'] #, '50nanFilter']

WT_col_script = '/home/uvi/be/nbo/MolClockAnalysis/scripts/analyzing/add_wt_to_vcf.py'
cellphy_exe = '/home/uvi/be/nbo/cellphy/cellphy.sh'
scite_script = '/home/uvi/be/nbo/MolClockAnalysis/simulations/scripts/run_scite.py'
scite_exe = '/home/uvi/be/nbo/infSCITE/infSCITE'

scite_time = 1000
scite_mem = 10
cellphy_time = 1000
cellphy_mem = 2

KEEP_GOING = False
CLOCK_DIR = ''


def run_bash(cmd_raw, bsub=True, time=30, mem=2):
    if bsub:
        cmd = f"sbatch -t {time} -p amd-shared --qos amd-shared --mem {mem}G " \
            f"--wrap '{cmd_raw}'"
    else:
        cmd = cmd_raw

    print(f'Running: {cmd}')
    subp = subprocess.Popen(cmd,
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = subp.communicate()
    subp.wait()

    if not cmd.startswith('sbatch'):
        print(str(stdout), str(stderr))
    if stderr != b'' and not KEEP_GOING:
        exit()
    print('\n')


def print_masterlist(out_file='vcf_masterlist.txt'):
    out_str = ''
    for data_dir, sub_dirs in data_dirs.items():
        data_set = data_dir.replace('_Monica', '')
        for sub_dir in sub_dirs:
            if sub_dir == 'all' and len(sub_dirs) > 1:
                continue
            vcf_dir = os.path.join(base_dir, data_dir, CLOCK_DIR, sub_dir)
            for data_filter in data_filters:
                vcf_name = f'{data_set}.{data_filter}_outg.vcf.gz'
                out_str += os.path.join(vcf_dir, vcf_name) + '\n'

    with open(out_file, 'w') as f:
        f.write(out_str)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--local', action='store_true',
        help='Run locally instead of HPC.')
    parser.add_argument('-r', '--replace', action='store_true',
        help='Overwrite already existing files.')
    parser.add_argument('-k', '--keep_going', action='store_true',
        help='Dont exit on errors.')
    parser.add_argument('-cd', '--clock_dir', type='str', default='ClockTest'
        help='Clock directory name.')
    parser.add_argument('-c', '--check', action='store_true',
        help='Check only if files exist, do not run anything.')
    parser.add_argument('-m', '--master', default='', type=str,
        help='Print master file list to file.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    KEEP_GOING = args.keep_going
    CLOCK_DIR = args.clock_dir

    if args.master:
        print_masterlist(args.master)
    if args.check:
        vcf_exist = True
        tree_exist = True

    # Iterate data sets
    for data_dir, sub_dirs in data_dirs.items():
        data_set = data_dir.replace('_Monica', '')
        # Iterate sub sets
        for sub_dir in sub_dirs:
            vcf_dir = os.path.join(base_dir, data_dir, CLOCK_DIR, sub_dir)
            # Iterate nan filters
            for data_filter in data_filters:
                vcf_raw_name = f'{data_set}.{data_filter}.vcf'
                vcf_raw_file = os.path.join(vcf_dir, vcf_raw_name)
                vcf_name = f'{data_set}.{data_filter}_outg.vcf.gz'
                vcf_file = os.path.join(vcf_dir, vcf_name)

                cellphy_out = vcf_file + '.raxml.bestTree'
                scite_out = os.path.join(vcf_dir, 'scite_dir',
                    f'{data_set}.{data_filter}_outg_ml0.newick')

                if args.check:
                    if not os.path.exists(vcf_file):
                        print(f'\tMissing vcf:\t\t{vcf_file}')
                        all_exist = False
                    if not os.path.exists(cellphy_out) and (sub_dir != 'all' \
                            or (sub_dir == 'all' and len(sub_dirs) == 1)):
                        print(f'\tMissing cellphy tree:\t{cellphy_out}')
                        tree_exist = False
                    if not os.path.exists(scite_out) and (sub_dir != 'all' \
                            or (sub_dir == 'all' and len(sub_dirs) == 1)):
                        print(f'\tMissing scite tree:\t{scite_out}')
                        tree_exist = False
                    continue

                # Copy file from monicas dir
                if sub_dir == 'all':
                    if 'dec21' in monica_dir:
                        # monica_file = os.path.join(monica_dir, f'{data_set}_dec21',
                        #     'all.all_chr.filtered.vcf')
                        monica_file = os.path.join(monica_dir, f'{data_set}_dec21',
                            f'{data_set}_ado_filtered.vcf')
                    else:
                        monica_file = os.path.join(monica_dir, f'{data_set}.vcf')
                    # Zip, add WT column, and index
                    if data_filter == 'all':
                        if not os.path.exists(vcf_file) or args.replace:
                            shutil.copyfile(monica_file, vcf_raw_file)
                            wt_col_cmd = f'python {WT_col_script} -i {vcf_raw_file}'
                            run_bash(wt_col_cmd, False)

                            unzip_file = vcf_file.replace('.gz', '')
                            zip_cmd = f'bgzip -f {unzip_file} && tabix {vcf_file}' \
                                f'&& chmod 755 {vcf_file}'
                            run_bash(zip_cmd, False)
                    # Filter
                    else:
                        if not os.path.exists(vcf_file) or args.replace:
                            base_file = os.path.join(vcf_dir, f'{data_set}.all_outg.vcf.gz')
                            flt_val = float(data_filter[:2]) / 100
                            flt_cmd = f'bcftools filter -i \'F_PASS(GT!="mis") ' \
                                f'> {flt_val}\' -O z -o {vcf_file} {base_file} ' \
                                f'&& chmod 755 {vcf_file}'
                            run_bash(flt_cmd, False)
                    if len(sub_dirs) > 1:
                        continue
                # Copy file from 'all' dir
                else:
                    if not os.path.exists(vcf_file) or args.replace:
                        base_file = os.path.join(base_dir, data_dir, CLOCK_DIR,
                            'all', vcf_name)
                        sample_file = os.path.join(vcf_dir, 'samples.txt')
                        cp_cmd = f'bcftools view --samples-file {sample_file} ' \
                            f'{vcf_file} {base_file} ' \
                            f'&& chmod 755 {vcf_file}'
                        run_bash(cp_cmd, False)


                tree_cmds = []
                if not os.path.exists(cellphy_out) or args.replace:
                    tree_cmds.append(
                        (f"{cellphy_exe} SEARCH -r -t 1 -z -l {vcf_file}",
                            cellphy_time, cellphy_mem)
                    )

                if not os.path.exists(scite_out) or args.replace:
                    tree_cmds.append(
                        (f"python3 {scite_script} -e {scite_exe} -s 1000000 " \
                            f"--verbose -p {data_set}.{data_filter}_outg {vcf_file}",
                        scite_time, scite_mem)
                    )

                for cmd, time, mem in tree_cmds:
                    if args.local:
                        run_bash(cmd, False)
                    else:
                        run_bash(cmd, True, time, mem)

    if args.check and vcf_exist:
        print('All vcf files exist')
    if args.check and tree_exist:
        print('All tree files exist')