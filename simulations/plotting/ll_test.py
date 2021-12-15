#!/usr/bin/env python3

import argparse
import numpy as np

def calc_G10N_likelihood(reads, eps=None, delta=None, gamma=None):
    # Line 4647 in cellcoal
    # CellCoal order: AA AC AG AT CC CG CT GG GT TT
    ll_GT = {}
    for ia1, A1 in enumerate(['A', 'C', 'G', 'T']):
        for ia2, A2 in [(0, 'A'), (1, 'C'), (2, 'G'), (3, 'T')][ia1:]:
            ll = 0

            g0 = 0
            g1 = 0
            g2 = 0
            for ib, read in enumerate(reads):
                if gamma:
                    p_bA1 = four_temp_ampl(ib == ia1, eps, gamma) 
                    p_bA2 = four_temp_ampl(ib == ia2, eps, gamma)
                else:
                    p_bA1 = GATK(ib == ia1, eps) 
                    p_bA2 = GATK(ib == ia2, eps)

                if delta:
                    # Lines 4688f
                    g0 += read * np.log10(0.5 * p_bA1 + 0.5 * p_bA2)
                    g1 += read * np.log10(p_bA1)
                    g2 += read * np.log10(p_bA2)
                else:
                    ll += read * np.log10(0.5 * p_bA1 + 0.5 * p_bA2)

            if delta:
                g0 = np.log10(1 - delta) + g0
                g1 = np.log10(delta/2) + g1
                g2 = np.log10(delta/2) + g2

                # t0 = (1 - delta) * g0
                # t1 = (delta/2) * g1
                # t2 = (delta/2) * g2

                t0 = g0
                t1 = g1
                t2 = g2

                ll = np.logaddexp(t0, np.logaddexp(t1, t2))
                
                ll_GT[A1 + A2] = round(ll, 2)

    x = np.array(list(ll_GT.values()))
    ll_norm = _normalize_log(x)
    ll_norm_0 = ll_norm - ll_norm.max()
    return ll_norm_0.round(2), x


def _normalize_log(probs):
        max_i = np.argmax(probs, axis=0)
        try:
            log_probs_norm = probs - probs[max_i] - np.log1p(np.sum(
                np.exp(probs[np.arange(probs.size) != max_i] - probs[max_i])
            ))
        except FloatingPointError:
            if probs[0] > probs[1]:
                return np.array([0, log_EPSILON])
            else:
                return np.array([log_EPSILON, 0])
        else:
            return log_probs_norm


def GATK(is_same, eps):
    if is_same:
        return 1 - eps
    else:
        return eps / 3


# Line 4675
def four_temp_ampl(is_same, eps, gamma):
    if is_same:
        return (1 - gamma) * (1 - eps) + gamma * eps / 3
    else:
        return (1 - gamma) * eps / 3 + (gamma / 3) * (1 - eps / 3)



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('read_counts', type=int, nargs=4,
        help='Read counts in the order: A C G T')
    parser.add_argument('-e', '--epsilon', type=float, default=1e-6,
        help='Sequencing error. Default = 1e-6.')
    parser.add_argument('-d', '--delta', type=float, default=0.2,
        help='ADO rate. Default = 0.2.')
    parser.add_argument('-g', '--gamma', type=float, default=1e-2,
        help='Amplification error. Default = 1e-2.')
    args = parser.parse_args()
    return args



if __name__ == '__main__':
    args = parse_args()
    ll, ll_raw = calc_G10N_likelihood(args.read_counts,
        eps=args.epsilon, delta=args.delta, gamma=args.gamma)
    print(np.exp(ll_raw * 1/ np.sum(args.read_counts)).sum())
    A = ['A', 'C', 'G', 'T']
    i = 0
    for i1, A1 in enumerate(A):
        for i2, A2 in enumerate(A[i1:]):
            print(f'{A1}|{A2}:\t{ll[i]: >8.2f}\t({ll_raw[i]: >8.2f})')
            i += 1