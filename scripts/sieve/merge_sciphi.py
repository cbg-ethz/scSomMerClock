#!/usr/bin/env python3

import os


def merge_standard(in_files, out_file):
    # file_endings = ['_mut2Sample.tsv', '.probs', '.gv', '.params.txt', '.vcf']

    header = ''
    body = ''
    sorted_files = sorted(
        [(i.split('.')[-2].replace('X', '23').replace('Y', '24'), i) \
                for i in in_files],
        key=lambda x: x[0])
    for i, (chrom, chrom_file) in enumerate(sorted_files):
        with open(chrom_file , 'r') as f:
            for line in f:
                if line.startswith('#'):
                    if i == 0:
                        header += line
                    continue
                body += line
   
    if not out_file:
        out_dir = os.path.sep.join(in_files[0].split(os.path.sep)[:-3])
        base_name = '.'.join(os.path.basename(in_files[0]).split('.')[:-2])
        out_file = os.path.join(out_dir, '{}.all.vcf'.format(base_name))
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

    with open(out_file, 'w') as f_out:
        f_out.write('{}{}'.format(header, body))


def merge_readCounts(in_files, out_file):
    bg = [{}, {}, {}, {}, {}]
    
    cand_sites = 0
    bg_sites = 0

    def set_chr(x):
        try:
            return int(x)
        except ValueError:
            if x == 'X':
                return 23
            elif x == 'Y':
                return 24
            else:
                raise IndexError

    try:
        in_files_sorted = sorted(in_files,
            key=lambda x: set_chr(x.split(os.path.sep)[-3].split('.')[-1]))
    except IndexError:
        in_files_sorted = sorted(in_files,
            key=lambda x: set_chr(x.split(os.path.sep)[-1].split('.')[-2]))

    for i, in_file in enumerate(in_files_sorted):
        with open(in_file , 'r') as f:
            file_raw = f.read().strip().split('\n')

        if i == 0:
            sample_str = '\n'.join(file_raw[:2])
            mut_str = file_raw[6]

        cand_sites += int(file_raw[3])
        bg_sites += int(file_raw[5])

        bg_row = -1
        for line_no, line in enumerate(file_raw[7:]):
            if line  == '=background=':
                bg_row = 0
            elif bg_row >= 0:
                bg_elements = line.split('\t')
                for bg_j, bg_element in enumerate(bg_elements):
                    try:
                        i1, i2 = bg_element.split(',')
                    except:
                        print(in_file, line_no, bg_elements[bg_j-1], bg_element)
                        exit()

                    try:
                        bg[bg_row][i1] += int(i2)
                    except KeyError:
                        bg[bg_row][i1] = int(i2)

                bg_row += 1
            else:
                mut_str += '\n' + line

    bg_str = '=background='
    for bg_line_out in bg:
        bg_str_new = '\t'.join(['{},{}'.format(*i) \
            for i in sorted(bg_line_out.items(), key=lambda x: int(x[0]))])
        bg_str += '\n' + bg_str_new

    par_str = '=numCandidateMutatedSites=\n{}\n=numBackgroundSites=\n{}' \
        .format(cand_sites, bg_sites)

    if not out_file:
        out_dir = os.path.dirname(in_files[0])
        out_file = os.path.join(out_dir, 'SciPhi_merged.tsv')
        print('Out_file : {}'.format(out_file))
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

    out_str = '{}\n{}\n{}\n{}'.format(sample_str, par_str, mut_str, bg_str)
    with open(out_file, 'w') as f_out:
        f_out.write(out_str)
    

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str,  nargs='*',
        help='Path to input tsv files')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Output file. Default = <INPUT_DIR>/SciPhi_merged.tsv.')
    parser.add_argument('-st', '--standard', action='store_true',
        help='Merge standard SciPhi output instead of read counts.')
    
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        merge_readCounts(snakemake.input, snakemake.output[0])
    else:
        import argparse
        args = parse_args()
        if args.standard:
            merge_standard(args.input, args.output)
        else:
            merge_readCounts(args.input, args.output)
