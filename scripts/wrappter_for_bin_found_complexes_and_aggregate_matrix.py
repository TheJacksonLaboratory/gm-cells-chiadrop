#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 qizai <jianhao2@illinois.edu>
#
# Distributed under terms of the MIT license.

"""
This script will search for all lp{0, 1, ..., N}_{right, left, middle}_100_region.csv
file. Each line in the csv file is a complex, with at least two fragment intersected with the 
loop region. 'complex_direction' indicate the side of motif this complex has at least 
one fragment overlap with.

NOTE [TODO]:
    annotation file need to contain 'convergence' info. Original ChIA-PET annot has a extended 
    annotated file which has it. 
    pseudo_annot needs to be change accordingly.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
import ipdb
import os
import functools
import tqdm
import argparse


def mainfn(args):
    p2loop_annot = args.p2loop_annot
    p2loop_tag = args.p2loop_tag
    p2root = args.p2intersected_complex_folder
    p2res_root = args.p2binning_results_saved_folder
    expr_name = args.expr_name
    nbins = args.nbins
    pseudo_loop = args.pseudo
    
    
    annot_col_name = ['left_chr', 'left_start', 'left_end', 
                      'right_chr', 'right_start', 'right_end',
                      'PET count', 'left_max_intensity', 'right_max_intensity',
                      'left_max_index', 'right_max_index', 'loop_ID',
                      'left_motif_chr', 'left_motif_start', 'left_motif_end', 
                      'left_motif_strand', 'left_distance',
                      'right_motif_chr', 'right_motif_start', 'right_motif_end', 
                      'right_motif_strand', 'right_distance']
    if pseudo_loop:
        df_loop_annot = pd.read_csv(p2loop_annot, sep='\t',
                                    names = annot_col_name)
        df_loop_annot = df_loop_annot.assign(convergence = 'convergence')
    else:
        # annot_col_name = ['bias', 'convergence', 'NULL motif']
        df_loop_annot_raw = pd.read_csv(p2loop_annot, sep='\t',
                                    names = annot_col_name)
        df_loop_tag = pd.read_csv(p2loop_tag, sep = '\t', index_col = 0)
        loop_conv_tag = df_loop_tag.loc[df_loop_annot_raw.index]['convergence'].values
        df_loop_annot = df_loop_annot_raw.assign(convergence = loop_conv_tag)

    p2region = os.path.join(p2root, expr_name)

    left_region_set = set([x.split('_')[0]
                           for x in os.listdir(p2region) if ('right' in x)])
    right_region_set = set([x.split('_')[0]
                            for x in os.listdir(p2region) if ('left' in x)])

    middle_region_set = set([x.split('_')[0]
                            for x in os.listdir(p2region) if ('middle' in x)])

    common_regions = list(right_region_set.intersection(left_region_set).intersection(middle_region_set))

    # ---- filter out the larger loops.
    com_idxs = df_loop_annot['loop_ID'].isin(common_regions)
    df_coms = df_loop_annot[com_idxs]
    com_convergence_idxs = (df_coms['convergence'] == 'convergence')
    df_coms_conv = df_coms[com_convergence_idxs]
    coms_conv_2kb_idx = df_coms_conv.apply(lambda x: (x.right_end - x.left_start) > 2e5,
                                           axis=1)
    df_coms_conv_2kb = df_coms_conv[coms_conv_2kb_idx]

    # ----- main loop ----
    df_combine_mid_right_bin_vectors = df_combine_mid_left_bin_vectors = pd.DataFrame()
    t1 = time.time()
    qbar = tqdm.tqdm(total=100)
    r_count = 0
    test_regions = common_regions
    total_regions = len(test_regions)
    one_tenth_total = int(total_regions / 100)
    for region in test_regions:
        r_count += 1
        if r_count % one_tenth_total == 0:
            qbar.update(1)
        region_loop_idx = (df_loop_annot['loop_ID'] == region)
        region_start, region_end = df_loop_annot[region_loop_idx][[
            'left_start', 'right_end']].values[0]
        left_motif_region = df_loop_annot[region_loop_idx][[
            'left_motif_start', 'left_motif_end']].values[0]
        right_motif_region = df_loop_annot[region_loop_idx][[
            'right_motif_start', 'right_motif_end']].values[0]

        middle_file_name = '{}_middle_100_region.csv'.format(region)
        p2df_middle = os.path.join(p2region, middle_file_name)
        df_middle = pd.read_csv(p2df_middle, index_col=0).drop_duplicates()

        left_file_name = '{}_left_100_region.csv'.format(region)
        p2df_left = os.path.join(p2region, left_file_name)
        df_left = pd.read_csv(p2df_left, index_col=0).drop_duplicates()

        right_file_name = '{}_right_100_region.csv'.format(region)
        p2df_right = os.path.join(p2region, right_file_name)
        df_right = pd.read_csv(p2df_right, index_col=0).drop_duplicates()

        if ((region_end - region_start) <= 2 * 1e5) or (df_loop_annot[region_loop_idx]['convergence'] != 'convergence').all():
            continue

        # Gather middle and right
        inline_snippet = functools.partial(combine_two_part, nbins=nbins, region = region,
                                           region_start=region_start, region_end=region_end,
                                           left_motif_region=left_motif_region,
                                           right_motif_region=right_motif_region)

        df_mid_overlap_right, df_right_tmp = inline_snippet(df_middle.copy(), df_right.copy(),
                                                            mid_direction='Right',
                                                            sid_direction='Left')
        df_combine_mid_right_bin_vectors = pd.concat(
            (df_combine_mid_right_bin_vectors, df_mid_overlap_right, df_right_tmp))

        df_mid_overlap_left, df_left_tmp = inline_snippet(df_middle.copy(), df_left.copy(),
                                                          mid_direction='Left',
                                                          sid_direction='Right')
        df_combine_mid_left_bin_vectors = pd.concat(
            (df_combine_mid_left_bin_vectors, df_mid_overlap_left, df_left_tmp))

    t2 = time.time()
    print(t2 - t1)
    qbar.close()

    # ---- save the output files.
    file_name_template = '{}bin_intersected_complex_in_both_region_{}_100_expansion.csv'
    p2res_root = os.path.join(p2res_root, expr_name)
    p2mid_left_bin = os.path.join(p2res_root, 
                                  file_name_template.format(nbins, 'left'))
    p2mid_right_bin = os.path.join(p2res_root,
                                   file_name_template.format(nbins, 'right'))

    if not os.path.isdir(p2res_root):
        os.makedirs(p2res_root)

    df_combine_mid_left_bin_vectors.to_csv(p2mid_left_bin)
    df_combine_mid_right_bin_vectors.to_csv(p2mid_right_bin)
    return df_combine_mid_left_bin_vectors, df_combine_mid_right_bin_vectors


def get_chrom_name_and_coord_array(List_of_frag_cord):
    list_of_frag_coord = List_of_frag_cord.split(';')

    # get rid of chr<>, split number into a tuple.
    tmp_list = [x.split(':')[1] for x in list_of_frag_coord]
    chrom_name = list_of_frag_coord[0].split(':')[0]

    tmp_list = [x.split('-') for x in tmp_list]

    # order the fragment by starting site.
    tmp_list.sort(key=lambda x: int(x[0]))

    # convert fragment to integer.
    # Note that end site has '(x)' at the end.
    tmp_list = [[int(x[0]), int(x[1].split('(')[0])] for x in tmp_list]
    frag_coord_filtered_array = np.array(tmp_list)
    num_intersected_frag = len(frag_coord_filtered_array)

    return chrom_name, frag_coord_filtered_array, num_intersected_frag


def bin_frag_array_by_loop_and_midpoint(frag_array, nbins, region_start, region_end):
    region_span = region_end - region_start

    # use midpoint of fragment as indicator. Normalized it in region span.
    frag_midpoint_array = np.mean(frag_array, axis=1)
    norm_frag_midpoint_array = (
        frag_midpoint_array - region_start) / region_span

    # scale up the normalization to nbins, and use largest integer as bin index.
    # |.bin0.|.bin1.|.bin2.|
    #   a     b   c   d
    # [a, b, c ,d ] -> [0, 1, 1, 2]
    digitize_frag_array = np.floor(
        norm_frag_midpoint_array * nbins).astype(int)
    digitize_frag_len = max(digitize_frag_array) - min(digitize_frag_array)

    return digitize_frag_array, digitize_frag_len


def List_of_frag_coord_to_all(List_of_frag_cord, nbins, region_start, region_end):
    chrom_name, frag_coord_filtered_array, num_intersect_frag = get_chrom_name_and_coord_array(
        List_of_frag_cord)
    digitize_frag_array, digitize_frag_len = bin_frag_array_by_loop_and_midpoint(frag_coord_filtered_array,
                                                                                 nbins,
                                                                                 region_start,
                                                                                 region_end)
    return digitize_frag_array, num_intersect_frag, digitize_frag_len, chrom_name


def get_motif_region_bin_idx(motif_region, nbins, region_start, region_end):
    '''
    motif_region [m_start, m_end]
    return range(m_start_bin, m_end_bin)
    '''
    ext_1kbs = np.array([-1000, 1000])
    expand_motif_region = motif_region + 4 * ext_1kbs

    region_span = region_end - region_start
    norm_expand_motif_region = (
        expand_motif_region - region_start) / region_span
    digitize_norm_exp_m_region = np.floor(
        norm_expand_motif_region * nbins).astype(int)
    clipped_digitize_norm_exp_m_region = np.clip(
        digitize_norm_exp_m_region, 0, nbins - 1).astype(int)

    return clipped_digitize_norm_exp_m_region


def df_inplace_binning_modify(df_middle, region, nbins, region_start, region_end,
                              left_motif_region, right_motif_region):

    return_col = ['{}_bin_index_vectors'.format(nbins), 'gem_intersected_span',
                  'num_intersected_frag',
                  'left_motif_{}_bin_start'.format(
                      nbins), 'left_motif_{}_bin_end'.format(nbins),
                  'right_motif_{}_bin_start'.format(
                      nbins), 'right_motif_{}_bin_end'.format(nbins),
                  'chrom', 'GEM_ID', 'loop_ID']
    if df_middle.shape[0] == 0:
        return pd.DataFrame(columns=return_col)

    def test_fn(x): return List_of_frag_coord_to_all(
        x, nbins, region_start, region_end)
    complex_bin_idx, complex_num_intersect_frag, complex_bin_len, complex_chrom = zip(
        *(df_middle['Intersected_list_of_frag_coordinates'].apply(test_fn)))
    complex_bin_idx = list(complex_bin_idx)
    def fn(x): return ';'.join([str(i) for i in x])
    complex_bin_idx = [fn(x) for x in complex_bin_idx]

    df_middle = df_middle.assign(**{'{}_bin_index_vectors'.format(nbins): complex_bin_idx,
                                    'num_intersected_frag': complex_num_intersect_frag,
                                    'gem_intersected_span': complex_bin_len,
                                    'chrom': complex_chrom})
    df_middle['loop_ID'] = region

    left_motif_bin_region = get_motif_region_bin_idx(left_motif_region, nbins,
                                                     region_start, region_end)
    df_middle['left_motif_{}_bin_start'.format(
        nbins)] = left_motif_bin_region[0]
    df_middle['left_motif_{}_bin_end'.format(nbins)] = left_motif_bin_region[1]

    right_motif_bin_region = get_motif_region_bin_idx(right_motif_region, nbins,
                                                      region_start, region_end)
    df_middle['right_motif_{}_bin_start'.format(
        nbins)] = right_motif_bin_region[0]
    df_middle['right_motif_{}_bin_end'.format(
        nbins)] = right_motif_bin_region[1]

    return df_middle[return_col]


def combine_two_part(df_mid, df_side, region, mid_direction, sid_direction,
                     nbins, region_start, region_end,
                     left_motif_region, right_motif_region):
    df_middle_anchor_motif_overlap = df_mid[df_mid['complex_direction']
                                            == mid_direction]
    df_side_anchor_motif_overlap = df_side[df_side['complex_direction']
                                           == sid_direction]

    df_middle_tmp = df_inplace_binning_modify(df_middle_anchor_motif_overlap, region,
                                              nbins, region_start, region_end,
                                              left_motif_region, right_motif_region)
    df_middle_tmp['position'] = 'mid_region_overlap_{}_motif'.format(
        mid_direction)

    df_side_tmp = df_inplace_binning_modify(df_side_anchor_motif_overlap, region,
                                            nbins, region_start, region_end,
                                            left_motif_region, right_motif_region)
    df_side_tmp['position'] = '{}_expand_region'.format(mid_direction)
    return df_middle_tmp, df_side_tmp


if __name__ == '__main__':
    # p2loop_annot = '../../results/pseudo_loop_annotataion_1000_loops_5_random_sample_per_loop.csv'
    # p2root = '/Users/qizai/projects/ChIA_Drop/results/intersection_files/pseudo_loops/'
    # expr_type = 'SHG0180-181-182NR_hg38_cohesin_FDR_0.1_PASS_pseudo'
    # p2res_root = '../../results/bin_complex_for_left_middle_right_aggregate/CDH0002NR_hg38_CTCF_FDR_0.1_PASS_ALL/'

    parser = argparse.ArgumentParser()
    parser.add_argument('--p2loop_annot', type = str, 
            default = '../../results/pseudo_loop_annotataion_1000_loops_5_random_sample_per_loop.csv')
    parser.add_argument('--p2loop_tag', type = str,
            default = '../../../ChIA_PET/loop_tags',
            help = 'loop tags file from ChIA_PET result. Contains convergence tag for each loop \
                    by default the pseudo loop are all convergent, when not set, should used \
                    the .motifannot_added_labels file as loop annotation file.')

    parser.add_argument('--p2intersected_complex_folder', type = str,
            default = '../../results/intersection_files/')
    parser.add_argument('--p2binning_results_saved_folder', type = str,
            default = '../../results/bin_complex_for_left_middle_right_aggregate/')
    parser.add_argument('--expr_name', type = str,
            help = 'e.g. CDH002NF_hg38_CTCF_FDR_0.1_ALL/. It is the subfolder inside \
                    p2intersected_complex_folder and p2binning_results_saved_folder \
                    to separate different experiments')
    parser.add_argument('--pseudo', type = int,
            choices = [0, 1], default = 0, help = 'Has to be 1 when use pseudo loop')
    parser.add_argument('--nbins', type = int,
            default = 100)
    args = parser.parse_args()
    print(args)
    tmp_left, tmp_right = mainfn(args)
