#!/usr/bin/env python3

import argparse
import numpy as np
from scipy.special import comb

np.seterr(all='raise')


def get_combo_rec(l, t, r):
    # base case 1: less than 2 ref
    if l[-1][0] < 2:
        return
    # base case 2: negative alt
    elif l[-1][1] < 0:
        return
    # base case 3: true combination
    elif len(l) == t:
        r.append(l[1:])
        return
    # recursive case
    else:
        # substract blue\ref
        get_combo_rec(l + [(l[-1][0] - 1, l[-1][1])], t, r)
        # substract red\alt
        get_combo_rec(l + [(l[-1][0], l[-1][1] - 1)], t, r)
    return r


def get_prob_rec(k, n, a, s, p, inv=-1):
    # base case 1: negative alt
    if k < 0:
        pass
    # base case 2: less than 2 ref
    elif n == 2 and k != 0:
        pass
    elif n - k < 1:
        pass
    else:
        if inv == 0:
            s_plus = np.log((1-a) * k + a * (n - k)) - np.log(n)
        elif inv == 1:
            s_plus = np.log((1-a) * (n - k) + a * k) - np.log(n)
        else:
            s_plus = 0
        # base case 0: too unlikely, prob nearly 0
        if s + s_plus < -20:
            p.append(s + s_plus)
        # base case 3: true combination
        elif n == 2 and k == 0:
            p.append(s + s_plus)
        # recursive case
        else:
            # substract blue\ref
            get_prob_rec(k, n - 1, a, s + s_plus, p, inv=1)
            # substract red\alt
            get_prob_rec(k - 1, n - 1, a, s + s_plus, p, inv=0)
    return p


def log_sum(x):
    m = np.max(x)
    s = 0
    for x_i in x:
        if x_i != m:
            s += np.exp(x_i - m)
    return m + np.log1p(s)


def log_norm(x):
    max_i = np.argmax(x)
    x_exp = np.exp(x[np.arange(x.size) != max_i] - x[max_i])
    x_norm = x - x[max_i] - np.log1p(np.sum(x_exp))
    return np.exp(np.clip(x_norm, None, 0))


def get_prob(k, n, a):
    # start = np.log((1-a) * (n - k) + a * k) - np.log(n)
    x = get_prob_rec(k, n, a, 0, [])
    if len(x) > 0:
        return np.exp(log_sum(x))
    else:
        return 0


def get_prob_short(k, n, a):
    no_comb = np.log(comb(n - 2, k))
    try:
        prob_no_error = (n - k - 2) * np.log(1 - a) \
            + np.log(np.arange(2, n-k)).sum()
        prob_error = np.log(np.arange(0, k) * (1 - a) + (n - k) * a).sum()
    except:
        import pdb; pdb.set_trace()
    norm = np.log(np.arange(2, n)).sum()
    return np.exp(no_comb + prob_no_error + prob_error - norm)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', type=int, default='', help='# alt reads')
    parser.add_argument('-n', type=int, help='# total reads.')
    parser.add_argument('-a', type=float, default=0.001, help='Ampl. error rate.')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    print(f'Exact: {get_prob(args.k, args.n, args.a):6f}')
    print(f'Upper bound: {get_prob_short(args.k, args.n, args.a):6f}')
    exit()
    p_all = np.zeros((args.n-1, 2))
    for k in range(0, args.n-1, 1):
        p_all[k, 0] = get_prob(k, args.n, args.a)
        p_all[k, 1] = get_prob_short(k, args.n, args.a)
        print(f'p({k: >2}|{args.n}) = {p_all[k, 0]:.6f}\t'
            f'(upper bound: {p_all[k, 1]:.6f})')
    print(f'\nsum: {np.sum(p_all, axis=0)[0]:6f}')