#!/usr/bin/env python3

import argparse
import os
import re
import shutil
import subprocess
import tarfile


# MONICA_DIR = '/mnt/lustre/scratch/home/uvi/be/mva/singlecell/Projects/mol_clock/VariantCallsApril/filter2'
# MONICA_DIR = '/mnt/lustre/scratch/home/uvi/be/mva/singlecell/Projects/mol_clock/vcfs2022'

DATA_DIRS = {
    'CRC08': ['all', 'normal', 'TC', 'TD', 'TP'],
    'CRC09': ['all', 'TI', 'TM'],
    'H65': ['all'],
    'Li55': ['all', 'cancer', 'normal'],
    'Lo-P1': ['all'],
    'Lo-P2': ['all'],
    'Lo-P3': ['all'],
    'Ni8': ['all'],
    'S21_P1': ['all'],
    'S21_P2': ['all'],
    'W32': ['all', 'aneuploid', 'haploid', 'normal'],
    'W55': ['all', 'cancer', 'normal'],
    'Wu61': ['all', 'cancer_C', 'cancer_CA', 'polyps'],
    'Wu63': ['all', 'cancer', 'normal', 'polyps'],
    'X25': ['all', 'cancer']
}
WGS = ['S21_P1', 'S21_P2'] # 'CRC08', 'CRC09', 'Lo-P1', 'Lo-P2', 'Lo-P3'
DATA_FILTERS = ['all', '33nanFilter', '50nanFilter']

WT_col_script = '/home/uvi/be/nbo/MolClockAnalysis/scripts/analyzing/add_wt_to_vcf.py'
cellphy_exe = '/home/uvi/be/nbo/cellphy/cellphy.sh'
scite_script = '/home/uvi/be/nbo/MolClockAnalysis/simulations/scripts/run_scite.py'
scite_exe = '/home/uvi/be/nbo/infSCITE/infSCITE'

scite_time = 2880
scite_time_WGS = 4320
scite_mem = 12
cellphy_cores = 8
cellphy_time = 2880
cellphy_mem = 6

KEEP_GOING = False
SLURM = True


def run_bash(cmd_raw, bsub=True, cores=1, time=30, mem=2):
    if bsub:
        if SLURM:
            cmd = f"sbatch --cpus-per-task {cores} -t {time} -p shared " \
                f"--qos shared --mem {mem}G --wrap '{cmd_raw}'"
        else:
            cmd = f'bsub -n {cores} -W {time} -R "rusage[mem={mem*1000}]" ' \
                f'"{cmd_raw}"'
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
    if not KEEP_GOING and (stderr != b'' and not stderr.startswith(b'Warn:')):
        exit()
    print('\n')


def run_inference(args):
    if args.check:
        vcf_exist = True
        tree_exist = True

    # Iterate data sets
    for data_set, sub_dirs in DATA_DIRS.items():
        if data_set not in args.dataset:
            continue

        # Iterate sub sets
        for sub_dir in sub_dirs:
            vcf_dir = os.path.join(args.out_dir, data_set, sub_dir)
            if not os.path.exists(vcf_dir):
                os.makedirs(vcf_dir)

            if data_set.startswith('S21_'):
                data_filters = DATA_FILTERS + ['99nanFilter']
            else:
                data_filters = DATA_FILTERS

            # Iterate nan filters
            for data_filter in data_filters:
                if data_filter not in args.filters:
                    continue

                vcf_raw_name = f'{data_set}.{data_filter}.vcf.gz'
                vcf_raw_file = os.path.join(vcf_dir, vcf_raw_name)
                vcf_name = f'{data_set}.{data_filter}_outg.vcf.gz'
                vcf_file = os.path.join(vcf_dir, vcf_name)

                sample_file = os.path.join(vcf_dir, 'samples.txt')

                cellphy_out = vcf_file + '.raxml.bestTree'
                scite_out = os.path.join(vcf_dir, 'scite_dir',
                    f'{data_set}.{data_filter}_outg_ml0.newick')

                if args.check:
                    if not os.path.exists(vcf_file):
                        print(f'\tMissing vcf:\t\t{vcf_file}')
                        vcf_exist = False
                    if args.files_only:
                        continue
                    if sub_dir == 'all' and len(sub_dirs) > 1:
                        continue

                    if 'cellphy' in args.method \
                            and not os.path.exists(cellphy_out):
                        print(f'\tMissing cellphy tree:\t{cellphy_out}')
                        tree_exist = False
                    if 'scite' in args.method and not os.path.exists(scite_out) :
                        print(f'\tMissing scite tree:\t{scite_out}')
                        tree_exist = False
                    continue

                if not args.tree_only:
                    # Arrange files for processing
                    if sub_dir == 'all':
                        # Zip, add WT column, and index
                        if data_filter == 'all':
                            if not os.path.exists(vcf_file) or args.replace:
                                in_file = os.path.join(args.input, f'{data_set}.vcf')
                                # Compress and index if not done
                                if not os.path.exists(in_file + '.gz'):
                                    idx_cmd = f'bgzip -f {in_file} && tabix ' \
                                        f'{in_file}.gz && chmod 755 {in_file}.gz'
                                    run_bash(idx_cmd, False)
                                # Copy and filter cells that did not pass QC
                                cp_cmd = f'bcftools view --samples-file {sample_file} ' \
                                    f'--force-samples {in_file}.gz -O z ' \
                                    f'-o {vcf_raw_file}'
                                run_bash(cp_cmd, False)

                                # Add wild type column
                                wt_col_cmd = f'python {WT_col_script} -i {vcf_raw_file}'
                                run_bash(wt_col_cmd, False)

                                # Zip and tabify
                                unzip_file = vcf_file.replace('.gz', '')
                                zip_cmd = f'bgzip -f {unzip_file} && tabix {vcf_file} ' \
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
                        # Skip 'all' subdir if other subsets exist
                        if len(sub_dirs) > 1:
                            continue
                    # Copy file from 'all' dir
                    else:
                        if not os.path.exists(vcf_file) or args.replace:
                            base_file = os.path.join(args.out_dir, data_set,
                                    'all', vcf_name)
                            cp_cmd = f'bcftools view --samples-file {sample_file} ' \
                                f'{base_file} | bcftools filter -i ' \
                                f'\'N_PASS(GT="alt") != 0\' -O z -o {vcf_file} - ' \
                                f'&& chmod 755 {vcf_file}'
                            run_bash(cp_cmd, False)

                if args.files_only:
                    continue

                tree_cmds = []
                if 'cellphy' in args.method:
                    if not os.path.exists(cellphy_out) or args.replace:
                        tree_cmds.append(
                            (f'{cellphy_exe} FULL -o healthycell -r -z -l -y ' \
                                f'-t {cellphy_cores} {vcf_file}',
                            cellphy_time, cellphy_mem, cellphy_cores)
                        )
                if 'scite' in args.method:
                    if not os.path.exists(scite_out) or args.replace:
                        tree_cmds.append(
                            (f'python3 {scite_script} -e {scite_exe} -s 1000000 ' \
                                f'--verbose -p {data_set}.{data_filter}_outg ' \
                                f'{vcf_file}',
                            scite_time, scite_mem, 1)
                        )

                for cmd, time, mem, cores in tree_cmds:
                    if args.local:
                        run_bash(cmd, False)
                    else:
                        run_bash(cmd, True, cores, time, mem)

    if args.check and vcf_exist:
        print('All vcf files exist')
    if args.check and tree_exist:
        print('All tree files exist')


def compress_files(args):
    if not os.path.exists(args.compress):
        os.makedirs(args.compress)

    for data_set, sub_dirs in DATA_DIRS.items():
        if data_set not in args.dataset:
            continue

        maf_file = os.path.join(args.out_dir, f'{data_set}.maf')
        if os.path.exists(maf_file):
            shutil.copyfile(maf_file,
                os.path.join(args.compress, f'{data_set}.maf'))

        for sub_dir in sub_dirs:
            if sub_dir == 'all' and len(sub_dirs) > 1:
                continue

            vcf_dir = os.path.join(args.out_dir, data_set, sub_dir)
            if data_set.startswith('S21_'):
                data_filters = DATA_FILTERS + ['99nanFilter']
            else:
                data_filters = DATA_FILTERS

            for data_filter in data_filters:
                if data_filter not in args.filters:
                    continue

                vcf_name = f'{data_set}.{data_filter}_outg.vcf.gz'
                vcf_file = os.path.join(vcf_dir, vcf_name)

                base_name = f'{data_set}_{sub_dir}_{data_filter}'

                if not os.path.exists(vcf_file):
                    print(f'\tMissing vcf file: {vcf_file}')
                    continue

                shutil.copyfile(vcf_file,
                    os.path.join(args.compress, f'{base_name}.vcf.gz'))

                if 'cellphy' in args.method:
                    tree_file = os.path.join(vcf_dir, f'{vcf_name}.raxml.supportFBP')
                    log_file = os.path.join(vcf_dir, f'{vcf_name}.raxml.log')
                    if not os.path.exists(tree_file):
                        print(f'\tMissing cellphy tree file: {tree_file}')
                    else:
                        shutil.copyfile(tree_file,
                            os.path.join(args.compress, f'{base_name}.cellphy.newick'))
                        shutil.copyfile(log_file,
                            os.path.join(args.compress, f'{base_name}.cellphy.log'))

                if 'scite' in args.method:
                    tree_file = os.path.join(vcf_dir, 'scite_dir',
                        f'{data_set}.{data_filter}_outg_ml0.newick')
                    log_file = os.path.join(vcf_dir, 'scite_dir',
                        f'{data_set}.{data_filterz}_outg.log')
                    if not os.path.exists(tree_file):
                        print(f'\tMissing scite tree file: {tree_file}')
                    else:
                        shutil.copyfile(tree_file,
                            os.path.join(args.compress, f'{base_name}.scite.newick'))
                        shutil.copyfile(log_file,
                            os.path.join(args.compress, f'{base_name}.scite.log'))

    tar = tarfile.open(args.compress + '.tar.gz', 'w:gz')
    tar.add(args.compress)
    tar.close()
    # Remove uncompressed folder
    shutil.rmtree(args.compress)


def check_errors(args):
    for data_set, sub_dirs in DATA_DIRS.items():
        if data_set not in args.dataset:
            continue
        for sub_dir in sub_dirs:
            if sub_dir == 'all' and len(sub_dirs) > 1:
                continue
            vcf_dir = os.path.join(args.input, data_set, sub_dir)
            if data_set.startswith('S21_'):
                data_filters = DATA_FILTERS + ['99nanFilter']
            else:
                data_filters = DATA_FILTERS

            for data_filter in data_filters:
                if data_filter not in args.filters:
                    continue

                data_id = f'{data_set}_{sub_dir}_{data_filter}'
                if 'cellphy' in args.method:
                    log_file1 = os.path.join(vcf_dir,
                        f'{data_set}.{data_filter}_outg.vcf.gz.raxml.log')
                    log_file2 = os.path.join(args.input,
                        f'{data_set}_{sub_dir}_{data_filter}.cellphy.log')
                    if not os.path.exists(log_file1) \
                            and not os.path.exists(log_file2):
                        print(f'\tMissing cellphy log file: {data_id}')
                    else:
                        if os.path.exists(log_file1):
                            log_file = log_file1
                        else:
                            log_file = log_file2

                        with open(log_file, 'r') as f:
                            log = f.read().strip()
                        FN = float(re.search('ADO_RATE: (0.\d+(e-\d+)?)', log) \
                            .group(1))
                        if FN == 0:
                            print(f'\tCELLPHY - No errors detected in: {log_file}')
                if 'scite' in args.method:
                    log_file1 = os.path.join(vcf_dir, 'scite_dir',
                        f'{data_set}.{data_filter}_outg.log')
                    log_file2 = os.path.join(args.input,
                        f'{data_set}_{sub_dir}_{data_filter}.scite.log')
                    if not os.path.exists(log_file1) \
                            and not os.path.exists(log_file2):
                        print(f'\tMissing cellphy log file: {log_file}')
                    else:
                        if os.path.exists(log_file1):
                            log_file = log_file1
                        else:
                            log_file = log_file2

                        with open(log_file, 'r') as f:
                            log = f.read()
                        FN = float(re.search(
                            'best value for beta:\\\\t(\d.\d+(e-\d+)?)', log) \
                                .group(1))
                        if FN == 0:
                            print(f'\tSCITE - No errors detected in: {log_file}')

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help=f'Input directory.')
    parser.add_argument('-o', '--out_dir', type=str, default='',
        help='Output directory. Default = BASEDIR.')
    parser.add_argument('-l', '--local', action='store_true',
        help='Run locally instead of HPC.')
    parser.add_argument('-r', '--replace', action='store_true',
        help='Overwrite already existing files.')
    parser.add_argument('-k', '--keep_going', action='store_true',
        help='Dont exit on errors.')
    parser.add_argument('-f', '--files_only', action='store_true',
        help='Create only subset files, do not run tree inference.')
    parser.add_argument('-t', '--tree_only', action='store_true',
        help='Rerun only tree inference, do not re-create subset files.')
    parser.add_argument('-c', '--check', action='store_true',
        help='Check only if files exist, do not run anything.')
    parser.add_argument('-m', '--method', nargs='+', type=str,
        choices=['cellphy', 'scite'], default=['cellphy'],
        help=f'Tree inference method. Default = cellphy.')
    parser.add_argument('-comp', '--compress', default='', type=str,
        help='Compress files to folder and exit.')
    parser.add_argument('-d', '--dataset', type=str, nargs='+',
        choices=DATA_DIRS.keys(), default=DATA_DIRS.keys(),
        help=f'Datasets to process. Default = {DATA_DIRS.keys()}.')
    parser.add_argument('-fl', '--filters', type=str, nargs='+',
        choices=DATA_FILTERS + ['99nanFilter'],
        default=DATA_FILTERS + ['99nanFilter'],
        help=f'Datafilters to process. ' \
            f'Default = {[DATA_FILTERS + ["99nanFilter"]]}.')
    parser.add_argument('-cerr', '--check_errors', action='store_true',
        help='Check for runs where inferred error is 0.')
    parser.add_argument('--lsf', action='store_true',
        help='Run lsf queue submit command (bsub) instead of slurm (sbatch).')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    KEEP_GOING = args.keep_going

    if not args.out_dir:
        args.out_dir = args.input

    if args.lsf:
        SLURM = False

    if args.compress:
        compress_files(args)
    elif args.check_errors:
        check_errors(args)
    else:
        run_inference(args)