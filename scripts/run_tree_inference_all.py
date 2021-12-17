#!/usr/bin/env python3

import os
import subprocess


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
    'Wu63': ['all', 'cancer', 'cancer_polyp', 'normal', 'polyps'],
    'X25': ['all', 'cancer', 'normal']
}
data_filters = ['all', '33nanFilter', '50nanFilter']

cellphy_exe = '/home/uvi/be/nbo/cellphy/cellphy.sh'
scite_script = '/home/uvi/be/nbo/MolClockAnalysis/simulations/scripts/run_scite.py'
scite_exe = '/home/uvi/be/nbo/infSCITE/infSCITE'


if __name__ == '__main__':
    for data_dir, sub_dirs in data_dirs.items():
        data_set = data_dir.replace('_Monica', '')
        for sub_dir in sub_dirs:
            vcf_dir = os.path.join(base_dir, data_dir, 'ClockTest', sub_dir)
            for data_filter in data_filters:
                vcf_name = f'{data_set}.{data_filter}.vcf.gz'
                vcf_file = os.path.join(vcf_dir, vcf_name)

                cellphy_cmd = f"sbatch -t 720 -p amd-shared --qos amd-shared " \
                    f"--mem 2G --wrap '{cellphy_exe} SEARCH -r -t 1 -z -l " \
                    f"{vcf_file}'"
                scite_cmd = f"sbatch -t 1440 -p amd-shared --qos amd-shared " \
                    f"--mem 10G --wrap 'python3 {scite_script} -e {scite_exe} " \
                    f"-s 1000000 --verbose -p {data_set}.{data_filter} {vcf_file}'"


                for cmd in [cellphy_cmd, scite_cmd]:
                    tree_inf = subprocess.Popen(cmd,
                        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                    )
                    stdout, stderr = tree_inf.communicate()
                    tree_inf.wait()
                    print(str(stdout), str(stderr))
