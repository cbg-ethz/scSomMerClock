#!/usr/bin/env python3

import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser(prog='merge_sciphi',
        description='*** Merge SciPhi candidate sites identified per chromosome. ***')
    parser.add_argument('input', type=str,  nargs='*',
        help='Absolute or relative path(s) to input tsv files')
    parser.add_argument('-m', '--mode', type=int, default=0,
        help='SCIPHI output: 0 for Readcount, 1 for standard. Default = 0.')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Path to the output directory. Default = <INPUT_DIR>.')
    args = parser.parse_args()
    return args


def merge_standard(args):
    bg = [{}, {}, {}]
    for i, in_file in enumerate(args.input):
        with open(in_file , 'r') as f:
            file_raw = f.read().strip().split('\n')

        if i == 0:
            sample_str = '\n'.join(file_raw[:2])
            par_str = '\n'.join(file_raw[2:8])
            mut_str = file_raw[8]

        for mut_line in file_raw[9:]:
            if mut_line == '=background=':
                break
            else:
                mut_str += '\n' + mut_line

        for bg_i, bg_line in enumerate(file_raw[-3:]):
            bg_raw = bg_line.split('\t')

            for bg_j in range(len(bg_raw) // 2):
                try:
                    bg[bg_i][bg_raw[2 * bg_j]] += int(bg_raw[2 * bg_j + 1])
                except KeyError:
                    bg[bg_i][bg_raw[2 * bg_j]] = int(bg_raw[2 * bg_j + 1])

    bg_str = '=background='
    for bg_line_out in bg:
        bg_str_new = '\t'.join(['{}\t{}'.format(*i) for i in bg_line_out.items()])
        bg_str += '\n' + bg_str_new

    return sample_str, par_str, mut_str, bg_str


def merge_readCounts(args):
    bg = [{}, {}, {}, {}, {}] 
    cand_sites = 0
    bg_sites = 0

    for i, in_file in enumerate(args.input):
        with open(in_file , 'r') as f:
            file_raw = f.read().strip().split('\n')

        if i == 0:
            sample_str = '\n'.join(file_raw[:2])
            mut_str = file_raw[8]

        cand_sites += int(file_raw[3])
        bg_sites += int(file_raw[5])

        bg_flag = False
        for line_no, line in enumerate(file_raw[9:]):
            if bg_flag or line == '=background=':
                bg_flag = True
                bg_elements = line.split('\t')
                for bg_j, bg_element in enumerate(bg_elements):
                    try:
                        i1, i2 = bg_element.split(',')
                    except:
                        print(in_file, line_no, bg_elements[bg_j-1], bg_element)
                        exit()

                    try:
                        bg[bg_i][i1] += int(i2)
                    except KeyError:
                        bg[bg_i][i1] = int(i2)
                    except ValueError:
                        print(bg_j)
            else:
                mut_str += '\n' + line

    bg_str = '=background='
    for bg_line_out in bg:
        bg_str_new = '\t'.join(['{},{}'.format(*i) for i in bg_line_out.items()])
        bg_str += '\n' + bg_str_new

    par_str = '=numCandidateMutatedSites={}\n=numBackgroundSites={}' \
        .format(cand_sites, bg_sites)

    return sample_str, par_str, mut_str, bg_str


def main(args):
    if not args.output:
        args.output = os.path.sep.join(args.input[0].split(os.path.sep)[:-3])
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    if args.mode == 0:
        sample_str, par_str, mut_str, bg_str = merge_readCounts(args)
    else:
        sample_str, par_str, mut_str, bg_str = merge_standard(args)


    out_str = '{}\n{}\n{}\n{}'.format(sample_str, par_str, mut_str, bg_str)
    out_file = os.path.join(args.output, 'SciPhi_merged.tsv')
    with open(out_file, 'w') as f_out:
        f_out.write(out_str)
    

if __name__ == '__main__':
    args = parse_args()
    main(args)