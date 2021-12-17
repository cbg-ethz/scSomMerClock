#!/usr/bin/env python3

import os
import shutil
import subprocess


monica_dir = '/mnt/lustre/scratch/home/uvi/be/mva/singlecell/Projects/mol_clock/VariantCallsApril/dec21/'
base_dir = '/home/uvi/be/nbo/data/data/'

data_dirs = {
    'H65_Monica': ['all', 'cancer', 'normal'],
    'Li55': ['all', 'cancer', 'normal'],
    'Lo-P1': ['all'],
    'Lo-P2': ['all'],
    'Lo-P3': ['all'],
    'Ni8_Monica': ['all', 'cancer'],
    'S21_P1': ['all'],
    'S21_P2': ['all', 'cancer', 'left'],
    'W32_Monica': ['all', 'aneuploid', 'cancer', 'haploid', 'normal'],
    'W55': ['all', 'cancer', 'normal'],
    'Wu61': ['all', 'cancer', 'cancer_C', 'cancer_CA', 'cancer_CA', 'normal',
        'polyps'],
    'Wu63': ['all', 'cancer', 'cancer_polyps', 'normal', 'polyps'],
    'X25': ['all', 'cancer', 'normal']
}
data_filters = ['all', '33nanFilter', '50nanFilter']

cellphy_exe = '/home/uvi/be/nbo/cellphy/cellphy.sh'
scite_script = '/home/uvi/be/nbo/MolClockAnalysis/simulations/scripts/run_scite.py'
scite_exe = '/home/uvi/be/nbo/infSCITE/infSCITE'

scite_time = 600
cellphy_time = 300


def run_bash(cmd):
    subp = subprocess.Popen(cmd,
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = subp.communicate()
    subp.wait()
    print(f'Running: {cmd}')
    print(str(stdout), str(stderr))
    print('\n')


if __name__ == '__main__':
    # Iterate data sets
    for data_dir, sub_dirs in data_dirs.items():
        data_set = data_dir.replace('_Monica', '')
        # Iterate sub sets
        for sub_dir in sub_dirs:
            vcf_dir = os.path.join(base_dir, data_dir, 'ClockTest', sub_dir)
            # Iterate nan filters
            for data_filter in data_filters:
                vcf_name = f'{data_set}.{data_filter}.vcf.gz'
                vcf_file = os.path.join(vcf_dir, vcf_name)

                # Copy file from monicas dir
                if sub_dir == 'all':
                    monica_file = os.path.join(monica_dir, f'{data_set}_dec21',
                        'all.all_chr.filtered.vcf')
                    # Zip and index
                    if data_filter == 'all':
                        if not os.path.exists(vcf_file):
                            unzip_file = vcf_file.replace('.gz', '')
                            shutil.copyfile(monica_file, unzip_file)
                            zip_cmd = f'bgzip {unzip_file} && tabix {vcf_file}'
                            run_bash(zip_cmd)
                    # Filter
                    else:
                        if not os.path.exists(vcf_file):
                            base_file = os.path.join(vcf_dir, f'{data_set}.all.vcf.gz')
                            flt_val = float(data_filter[:2]) / 100
                            flt_cmd = f'bcftools filter -i \'F_PASS(GT!="mis") ' \
                                f'> {flt_val}\' -O z -o {vcf_file} {base_file}'
                            run_bash(flt_cmd)
                # Copy file from 'all' dir
                else:
                    if not os.path.exists(vcf_file):
                        base_file = os.path.join(base_dir, data_dir, 'ClockTest',
                            'all', vcf_file)
                        sample_file = os.path.join(vcf_dir, 'samples.txt')
                        cp_cmd = f"bcftools view --samples-file {sample_file} -O z " \
                            "-o {vcf_file} {base_file}"
                        run_bash(cp_cmd)

                tree_cmds = []

                cellphy_out = vcf_file + '.raxml.bestTree'
                if not os.path.exists(cellphy_out):
                    tree_cmds.append(
                        f"sbatch -t {cellphy_time} -p amd-shared --qos amd-shared " \
                        f"--mem 2G --wrap '{cellphy_exe} SEARCH -r -t 1 -z -l " \
                        f"{vcf_file}'"
                    )

                scite_out = os.path.join(vcf_dir, 'scite_dir',
                    f'{data_set}.{data_filter}_ml0.newick')
                if not os.path.exists(scite_out):
                    tree_cmds.append(
                        f"sbatch -t {scite_time} -p amd-shared --qos amd-shared " \
                        f"--mem 10G --wrap 'python3 {scite_script} -e {scite_exe} " \
                        f"-s 1000000 --verbose -p {data_set}.{data_filter} {vcf_file}'"
                    )

                for cmd in tree_cmds:
                    run_bash(cmd)
