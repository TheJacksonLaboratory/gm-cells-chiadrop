#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 qizai <jianhao2@illinois.edu>
#
# Distributed under terms of the MIT license.

"""
This is a script version for jupyter notebook: 
    generate_random_loop_annotation_for_NULL_hypothesis.ipynb

It takes the real annotation file, random selected convergent loop with 
length >= [x]kbs. Repeatedly sample across the same chromosome [k] times. 
And add the annotation with lopp_ID: 'pseduo_lp[i]'
"""

import pandas as pd
import numpy as np
import time
import os
import ipdb
import tqdm
import argparse


def random_choose_from_list(low, high, sample_size=1000, random_seed=None, replace=True):
    if not random_seed is None:
        np.random.seed(random_seed)
    if replace:
        sampled_list = np.random.randint(low=low, high=high, size=sample_size)
    else:
        sampled_list = np.random.choice(high, sample_size, replace=False)
    return sampled_list


def repeat_sample_for_each_loop(loop, chrom_size, repeat=5, random_seed=42, replace=False):
    '''
    every loop from the sampled loop is a convergent loop.
    '''
    df_fake_loop_annot = pd.DataFrame(columns=loop.keys(), index=range(5))

    # sample different locations in chromosome for new starting sites.
    start_index_list = random_choose_from_list(1, chrom_size, repeat,
                                               random_seed=random_seed, replace=replace)

    columns_to_be_shifted = ['left_start', 'left_end', 'right_start', 'right_end',
                             'left_max_index', 'right_max_index']

    motif_pos_columns = ['left_motif_start', 'left_motif_end',
                         'right_motif_start', 'right_motif_end']

    left_motif_length = loop['left_motif_end'] - loop['left_motif_start']
    right_motif_length = loop['right_motif_end'] - loop['right_motif_start']

    for idx in range(repeat):
        # first repeat the loop.
        df_fake_loop_annot.loc[idx] = loop
        fake_left_start = start_index_list[idx]

        # modify the start/end.
        shifted_distance = loop['left_start'] - fake_left_start
        df_fake_loop_annot.loc[idx][columns_to_be_shifted] = loop[columns_to_be_shifted] - shifted_distance

        # change the [unextended] motif_location to 8kbs region around start and end.
        fake_left_motif_start = df_fake_loop_annot.loc[idx]['left_start'] + 4000
        fake_left_motif_end = fake_left_motif_start + left_motif_length

        fake_right_motif_end = df_fake_loop_annot.loc[idx]['right_end'] - 4000
        fake_right_motif_start = fake_right_motif_end - right_motif_length

        df_fake_loop_annot.loc[idx][motif_pos_columns] = [fake_left_motif_start, fake_left_motif_end,
                                                          fake_right_motif_start, fake_right_motif_end]

    return df_fake_loop_annot


def combine_all_sample_loop(df_sampled_loop, df_chrom_size, repeat=5,
                            random_seed=42, replace=True):
    num_sampled_loop = df_sampled_loop.shape[0]
    total_sample_size = repeat * num_sampled_loop
    df_fake_loop_annot = pd.DataFrame(columns=df_sampled_loop.columns,
                                      index=range(total_sample_size))

    np.random.seed(random_seed)
    for idx in range(num_sampled_loop):
        cur_loop = df_sampled_loop.iloc[idx]
        cur_chrom = cur_loop['left_chr']
        chrom_size = df_chrom_size.loc[cur_chrom]['size']

        df_tmp = repeat_sample_for_each_loop(cur_loop, chrom_size, repeat,
                                             replace=replace,
                                             random_seed=None)

        df_fake_loop_annot.loc[repeat * idx: repeat *
                               (idx + 1) - 1, :] = df_tmp.values

    # Reset loop_ID after genereation
    df_fake_loop_annot['loop_ID'] = [
        'lp{}'.format(x) for x in range(total_sample_size)]

    return df_fake_loop_annot


if __name__ == '__main__':
    #     p2real_annot = '../../../ChIA_PET/results/LHG0052H.e500.clusters.cis.bothanchint_G250.PETcnt_G9.motifannot_added_labels'
    # p2chrom_size = '../../../ChIA_PET/hg38.chrom.sizes'
    parser = argparse.ArgumentParser()
    parser.add_argument('--p2real_annot', type = str,
            default = '../../../ChIA_PET/LHG0052H.e500.clusters.cis.bothanchint_G250.PETcnt_G9.motifannot_added_labels')
    parser.add_argument('--p2chrom_size', type = str,
            default = '../../../ChIA_PET/hg38.chrom.sizes')
    parser.add_argument('--nloops', type = int,
            default = 1000)
    parser.add_argument('--repeat', type = int,
            default = 30)
    parser.add_argument('--p2pseudo_annot', type = str,
            default = 'pseudo_loop_annot')
    args = parser.parse_args()

    p2real_annot = args.p2real_annot
    p2chrom_size = args.p2chrom_size
    sample_size = args.nloops
    repeat = args.repeat

    df_real_loop_annot = pd.read_csv(p2real_annot, sep='\t',
                                     names=['left_chr', 'left_start', 'left_end', 'right_chr', 'right_start', 'right_end',
                                            'PET count', 'left_max_intensity', 'right_max_intensity',
                                            'left_max_index', 'right_max_index', 'loop_ID',
                                            'left_motif_chr', 'left_motif_start', 'left_motif_end', 'left_motif_strand', 'left_distance',
                                            'right_motif_chr', 'right_motif_start', 'right_motif_end', 'right_motif_strand', 'right_distance'])
    df_real_loop_annot = pd.read_csv(p2real_annot, sep='\t', index_col=0)

    df_chrom_size = pd.read_csv(p2chrom_size, sep='\t',
                                names=['size'],
                                index_col=0)

    df_real_convergent_loop = df_real_loop_annot[df_real_loop_annot['convergence'] == 'convergence']

    real_200kbs_loop_idx = df_real_convergent_loop.apply(lambda x: (x.right_end - x.left_start) > 2e5,
                                                         axis=1)
    df_real_conv_200kbs_loop = df_real_convergent_loop[real_200kbs_loop_idx]

    random_seed = 42

    selected_loop_index = random_choose_from_list(0, len(df_real_conv_200kbs_loop.index), sample_size=sample_size,
                                                  random_seed=random_seed, replace=False)

    df_sampled_loop = df_real_conv_200kbs_loop.iloc[selected_loop_index]

    df_sampled_loop_dedup = df_sampled_loop.drop_duplicates()

    df_tmp = combine_all_sample_loop(df_sampled_loop_dedup, df_chrom_size,
                                     repeat=repeat, random_seed=random_seed,
                                     replace=True)

    df_tmp = df_tmp.drop(['bias', 'convergence', 'NULL motif'], axis=1)

    df_tmp.to_csv(args.p2pseudo_annot,
                  sep='\t', header=False, index=False)
