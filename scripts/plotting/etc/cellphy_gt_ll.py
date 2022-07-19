#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt


TICK_FONTSIZE = 16
LABEL_FONTSIZE = 20


def plot_error_split(error, splits=20):
    fig, ax = plt.subplots(figsize=(16, 12))

    x = np.zeros(splits)
    y = np.zeros(splits)
    labels = []
    for i, split in enumerate(np.linspace(0, 0.2, splits)):
        epsilon = split
        gamma = error - epsilon
        x[i] = (1 - gamma) * epsilon / 3 + gamma * (1 - epsilon / 3) / 3
        y[i] = (1 - gamma) * (1 - epsilon) + gamma * epsilon / 3

        labels.append(r'$\epsilon$' + f': {epsilon:.3f}, ' + r'$\gamma$' \
            + f': {gamma:.3f}')

    plt.plot(x, y, marker='x', lw=0)
    x_offset = 0.000025
    y_offset = 0.0001

    for j, label in enumerate(labels):
        if j < splits / 2:
            plt.text(x[j] + x_offset, y[j] + y_offset, label)
        else:
            plt.text(x[j] + x_offset, y[j] - y_offset, label)

    ax.set_ylabel(r'$b = A$', fontsize=LABEL_FONTSIZE)
    ax.set_xlabel(r'$b \neq A$', fontsize=LABEL_FONTSIZE)
    plt.title(r'$\gamma + \epsilon$ = ' + str(error))

    plt.show()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--error', type=float, default=0.2,
        help='Sum of sequencing and amplification error. Default = 0.2.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    plot_error_split(args.error)