#!/usr/bin/env python3

import os
import re
import math
from statistics import mean, stdev
from utils import tail
    

def get_marginal_ll_str(h0, h1):
    diff = h1[0] - h0[0]
    if diff > 50:
        B01 = math.inf
    elif diff < -50:
        B01 = math.inf
    else:
        B01 = math.exp(diff)

    logB_01 = 2 * diff

    if logB_01 < 2:
        evidence = 'None'
        evidence_flag = 0
    elif logB_01 < 6:
        evidence = 'Positive'
        evidence_flag = 1
    elif logB_01 < 10:
        evidence = 'Strong'
        evidence_flag = 1
    else:
        evidence = 'Very Strong'
        evidence_flag = 1

    out_str = f'{h0[0]:.1f}\t{h0[1]:.0f}\t{h1[0]:.1f}\t{h1[1]:.0f}\t' \
        f'{logB_01:.1f}\t{B01:.0f}\t{evidence}'
    return out_str, logB_01, evidence_flag


def get_Bayes_factor(in_files, out_file, ss=False):
    scores = {}
    runtimes = {}
    for in_file in in_files:
        _, run, steps, model, _ = os.path.basename(in_file).split('.')
        steps = int(steps)
        with open(in_file, 'r') as f_score:
            score_raw = f_score.read().strip().split('\n')
        score = float(score_raw[-1].split('\t')[2])
        score_diff = abs(float(score_raw[-2].split('\t')[2]) - \
            float(score_raw[-3].split('\t')[2]))
        
        try:
            log_tail = tail(open(in_file.replace('.lstat', '.log'), 'r'), 1000)
        except FileNotFoundError:
            runtime = (0, -1)
        else:
            runtime_start = log_tail.index('Analysis completed in')
            runtime_end = log_tail[runtime_start:].index('seconds') + runtime_start
            runtime = re.findall('\d+', log_tail[runtime_start: runtime_end])

        if ss:
            if not log_tail:
                raise IOError('Log file with SS score not found: {}' \
                    .format(in_file.replace('.lstat', '.log')))
            ss_start = log_tail.index('Marginal likelihood (ln)')
            ss_end = log_tail[ss_start:].index('More statistics on stepping-stone')
            ss_raw = [i.strip() for i in \
                    log_tail[ss_start:ss_start + ss_end].strip().split('\n') \
                if i.strip() and not i.strip().startswith('-')]

            ss_mean = float(re.search('Mean: +(-\d+.\d+)',
                log_tail[ss_start:ss_start + ss_end]).group(1))

            ss_chain1 = float(re.search('1 +(-\d+.\d+)',
                log_tail[ss_start:ss_start + ss_end]).group(1))
            ss_chain2 = float(re.search('2 +(-\d+.\d+)',
                log_tail[ss_start:ss_start + ss_end]).group(1))
            ss_diff = abs(ss_chain1 - ss_chain2)
        else:
            ss_mean = None
            ss_diff = None

        # HH:MM:SS
        if len(runtime) == 3:
            runtime = 60 * 60 * int(runtime[0]) + 60 * int(runtime[1]) \
                + int(runtime[2])
        # MM:SS
        elif len(runtime) == 2:
            runtime = 60 * int(runtime[0]) + int(runtime[1])
        # SS
        else:
            runtime = int(runtime[0])
            
        if not steps in scores:
            scores[steps] = {run: {'clock': [-1] * 4, 'noClock': [-1] * 4}}
            runtimes[steps] = []
        if not run in scores[steps]:
            scores[steps][run] = {'clock': [-1] * 4, 'noClock': [-1] * 4}
        scores[steps][run][model] = (score, score_diff, ss_mean, ss_diff)
        runtimes[steps].append(runtime)

    out_str = ''
    summary_str = ''
    for step, step_info in sorted(scores.items()):
        summary_df = {'harmonic': {'h0': [], 'h0_diff': [], 'h1': [],
                'h1_diff': [], 'logB_01': [], 'evidence': []},
            'ss': {'h0': [], 'h0_diff': [], 'h1': [], 'h1_diff': [],
                'logB_01': [], 'evidence': []}}
        
        for run, run_info in step_info.items():
            h0 = run_info['clock']
            h1 = run_info['noClock']
            
            if h0[0] == -1 or h1[0] == -1:
                continue

            summary_df['harmonic']['h0'].append(h0[0])
            summary_df['harmonic']['h0_diff'].append(h0[1])
            summary_df['harmonic']['h1'].append(h1[0])
            summary_df['harmonic']['h1_diff'].append(h1[1])

            new_out_str, logB_01, evidence_flag = \
                get_marginal_ll_str(h0[:2], h1[:2])

            summary_df['harmonic']['logB_01'].append(logB_01)
            summary_df['harmonic']['evidence'].append(evidence_flag)

            out_str += f'{step}\t{run}\t' + new_out_str

            if ss:
                summary_df['ss']['h0'].append(h0[0])
                summary_df['ss']['h0_diff'].append(h0[1])
                summary_df['ss']['h1'].append(h1[0])
                summary_df['ss']['h1_diff'].append(h1[1])

                new_out_str_ss, logB_01_ss, evidence_flag_ss = \
                    get_marginal_ll_str(h0[2:], h1[2:])

                summary_df['ss']['logB_01'].append(logB_01_ss)
                summary_df['ss']['evidence'].append(evidence_flag_ss)

                out_str += '\t' + new_out_str_ss

            out_str += '\n'

        if len(step_info) < 2:
            continue

        summary_str += f"{step:.0E}\t{max(runtimes[step]):.2f}\t" \
            f"{mean(summary_df['harmonic']['h0']):.1f}\t" \
            f"{stdev(summary_df['harmonic']['h0']):.1f}\t" \
            f"{mean(summary_df['harmonic']['h0_diff']):.0f}\t" \
            f"{stdev(summary_df['harmonic']['h0_diff']):.0f}\t" \
            f"{mean(summary_df['harmonic']['h1']):.1f}\t" \
            f"{stdev(summary_df['harmonic']['h1']):.1f}\t" \
            f"{mean(summary_df['harmonic']['h1_diff']):.0f}\t" \
            f"{stdev(summary_df['harmonic']['h1_diff']):.0f}\t" \
            f"{mean(summary_df['harmonic']['logB_01']):.1f}\t" \
            f"{stdev(summary_df['harmonic']['logB_01']):.1f}\t" \
            f"{sum(summary_df['harmonic']['evidence'])}/" \
            f"{len(summary_df['harmonic']['evidence'])}"

        if ss:
            summary_str += f"\t{mean(summary_df['ss']['h0']):.1f}\t" \
            f"{stdev(summary_df['ss']['h0']):.1f}\t" \
            f"{mean(summary_df['ss']['h0_diff']):.0f}\t" \
            f"{stdev(summary_df['ss']['h0_diff']):.0f}\t" \
            f"{mean(summary_df['ss']['h1']):.1f}\t" \
            f"{stdev(summary_df['ss']['h1']):.1f}\t" \
            f"{mean(summary_df['ss']['h1_diff']):.0f}\t" \
            f"{stdev(summary_df['ss']['h1_diff']):.0f}\t" \
            f"{mean(summary_df['ss']['logB_01']):.1f}\t" \
            f"{stdev(summary_df['ss']['logB_01']):.1f}\t" \
            f"{sum(summary_df['ss']['evidence'])}/" \
            f"{len(summary_df['ss']['evidence'])}"
        summary_str += '\n'

    with open(out_file, 'w') as f_out:
        if ss:
            f_out.write('steps\trun\t' \
                'H0:clock (harmonic mean)\tH0 LL diff (harmonic mean)\t'
                'H1:noClock (harmonic mean)\tH1 LL diff (harmonic mean)\t' \
                '2log_e(B_01) (harmonic mean)\tB_01 (harmonic mean)\t' \
                'Evidence (harmonic mean)\t' \
                'H0:clock (ss)\tH0 LL diff (ss)\tH1:noClock (ss)\t'
                'H1 LL diff (ss)\t2log_e(B_01) (ss)\tB_01 (ss)\tEvidence (ss)\n')
        else:
            f_out.write('steps\trun\t' \
                'H0:clock (harmonic mean)\tH0 LL diff (harmonic mean)\t'
                'H1:noClock (harmonic mean)\tH1 LL diff (harmonic mean)\t' \
                '2log_e(B_01) (harmonic mean)\tB_01 (harmonic mean)\t' \
                'Evidence (harmonic mean)\n')
        f_out.write(out_str.strip())

    if summary_str:
        with open(out_file.replace('.tsv', '.short.tsv'), 'w') as f_out:
            if ss:
                f_out.write('steps\tMax. runtime [secs]\t'
                    'Avg. H0:clock (harmonic mean)\tStd. H0:clock (harmonic mean)\t'
                    'Avg. H0 LL diff (harmonic mean)\tStd. H0 LL diff (harmonic mean)\t'
                    'Avg. H1:noClock (harmonic mean)\tStd. H1:noClock (harmonic mean)\t'
                    'Avg. H1 LL diff (harmonic mean)\tStd. H1 LL diff (harmonic mean)\t'
                    'Avg. 2log_e(B_01) (harmonic mean)\t'
                    'Std. 2log_e(B_01) (harmonic mean)\tEvidence (harmonic mean)\t'
                    'Avg. H0:clock (ss)\tStd. H0:clock (ss)\tAvg. H0 LL diff (ss)\t'
                    'Std. H0 LL diff (ss)\tAvg. H1:noClock (ss)\tStd. H1:noClock (ss)\t'
                    'Avg. H1 LL diff (ss)\tStd. H1 LL diff (ss)\t'
                    'Avg. 2log_e(B_01) (ss)\tStd. 2log_e(B_01) (ss)\t'
                    'Evidence (ss)\n')
            else:
                f_out.write('steps\tMax. runtime [secs]\t'
                    'Avg. H0:clock (harmonic mean)\tStd. H0:clock (harmonic mean)\t'
                    'Avg. H0 LL diff (harmonic mean)\tStd. H0 LL diff (harmonic mean)\t'
                    'Avg. H1:noClock (harmonic mean)\tStd. H1:noClock (harmonic mean)\t'
                    'Avg. H1 LL diff (harmonic mean)\tStd. H1 LL diff (harmonic mean)\t'
                    'Avg. 2log_e(B_01)\tStd. 2log_e(B_01)\tEvidence\n')
            f_out.write(summary_str.strip())

    print(summary_str)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, nargs='+', help='Input files')
    parser.add_argument('-o', '--output', type=str, help='Output file.')
    parser.add_argument('-s', '--ss', action='store_true',
        help='Use stepping stone sampling instead of MCMC.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    if 'snakemake' in globals():
        get_Bayes_factor(snakemake.input, snakemake.output[0],
            snakemake.params.ss)
    else:
        import argparse
        args = parse_args()
        get_Bayes_factor(args.input, args.output, args.ss)