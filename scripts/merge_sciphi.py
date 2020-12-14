#!/usr/bin/env python3

import argparse
import os


def parse_args():
    parser = argparse.ArgumentParser(prog='merge_sciphi',
        description='*** Merge SciPhi candidate sites identified per chromosome. ***')
    parser.add_argument('input', type=str,  nargs='*',
        help='Absolute or relative path(s) to input tsv files')
    parser.add_argument('-o', '--output', type=str, default='',
        help='Path to the output directory. Default = <INPUT_DIR>.')
    args = parser.parse_args()
    return args


def main(args):
    if not args.output:
        args.output = os.path.sep.join(args.input[0].split(os.path.sep)[:-3])
    if not os.path.exists(args.output):
        os.makedirs(args.output)

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

    out_str = '{}\n{}\n{}\n{}'.format(sample_str, par_str, mut_str, bg_str)
    out_file = os.path.join(args.output, 'SciPhi_merged.tsv')
    with open(out_file, 'w') as f_out:
        f_out.write(out_str)
    

if __name__ == '__main__':
    args = parse_args()
    main(args)