#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 qizai <jianhao2@illinois.edu>
#
# Distributed under terms of the MIT license.

"""
This script will take the bedgraph file as input, and process it to create the output 
binning intensity.
"""
import os
from pyBedGraph import BedGraph

import numpy as np
import pandas as pd
import scipy
from scipy.stats import binom_test
import ipdb
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl


def get_max_intensity_in_same_len_bins(bedGraph, nbins, left_start, chrom_left, right_end,
                                       chrom_right=None, chrom_size = np.infty, flank_per=5):
    '''
    if chrom_right != None, then check if chrom_left == chrom_right.
        pyBedGraph can only query [chr, start, end] tuple.
    ----
    left_start: left anchor starting site
    right_end: right anchor ending site
    nbins: number of bins in the loop
    flank_per: percent of loop length to extend on both side.
    '''
    if chrom_right != None:
        if chrom_left != chrom_right:
            raise ValueError('row has anchors in different chromosome {}, {}'.format(chrom_left,
                                                                                     left_start))
    loop_length = right_end - left_start
    assert loop_length > 0
    flank_length = int(loop_length * flank_per / 100)


    start_idx = max(left_start - flank_length, 0)
    # ipdb.set_trace()
    end_idx = min(right_end + flank_length, chrom_size.values[0] - 1)

    if start_idx < 0 or start_idx > chrom_size.values[0] - 1:
        ipdb.set_trace()
    

    nbins_edges = np.linspace(start_idx, end_idx, nbins + 1, dtype=np.int32)

    start_list = nbins_edges[:-1]
    end_list = nbins_edges[1:]

    try:
        bin_values = bedGraph.stats(start_list=start_list,
                                    end_list=end_list,
                                    chrom_name=chrom_left,
                                    stat='max')
    except:
        print(chrom_left)
        print(end_list)
        print(start_idx)
        ipdb.set_trace()
    return bin_values


def get_aggregated_inten_for_each_class(df_binned_intensity_per_loop, nbins, catag):
    '''
    nbins \in {100, 500, 1000}
    catag \in {'bias', 'convergence', 'NULL motif'}
    '''
    bin_name = '{} binned intensity'.format(nbins)

    set_of_label = set(df_binned_intensity_per_loop[catag])
    label_list = list([x for x in set_of_label if x != 'na'])
    label_list.sort()

    total_num_loops_in_catag = (
        df_binned_intensity_per_loop[catag] != 'na').sum()

    chrom_list = list(set(df_binned_intensity_per_loop['chrom']))
    chrom_list.sort(key=lambda x: int(x[3:]) if x != 'chrX' else 24)
    chrom_list.append('whole genome')

    df_aggregate_sum = pd.DataFrame(columns=label_list, index=chrom_list)
    df_aggregate_mean = pd.DataFrame(columns=label_list, index=chrom_list)
    df_aggregate_var = pd.DataFrame(columns=label_list, index=chrom_list)

    for label in label_list:
        label_loop_idx = (df_binned_intensity_per_loop[catag] == label)
        for chrom in chrom_list[:-1]:
            # avoid whole genome.
            chrom_loop_idx = (df_binned_intensity_per_loop['chrom'] == chrom)

            tmp_df = df_binned_intensity_per_loop.loc[chrom_loop_idx &
                                                      label_loop_idx]
            sum_of_intensity = tmp_df[bin_name].sum()
            mean_of_intensity = tmp_df[bin_name].mean()
            var_of_intensity = np.stack(tmp_df[bin_name]).var(axis=0)

            df_aggregate_sum.loc[chrom][label] = sum_of_intensity
            df_aggregate_mean.loc[chrom][label] = mean_of_intensity
            df_aggregate_var.loc[chrom][label] = var_of_intensity

        df_aggregate_mean.loc['whole genome'][label] = df_binned_intensity_per_loop.loc[label_loop_idx][bin_name].mean()
        df_aggregate_sum.loc['whole genome'][label] = df_binned_intensity_per_loop.loc[label_loop_idx][bin_name].sum()
        total_var_per_label = np.stack(
            df_binned_intensity_per_loop.loc[label_loop_idx][bin_name]).var(axis=0)
        df_aggregate_var.loc['whole genome'][label] = total_var_per_label

    return df_aggregate_sum, df_aggregate_mean, df_aggregate_var


def aggre_by_sum(df_agg_sum, label, chrom='whole genome', scilent=False, p2f=None):
    y = df_agg_sum.loc[chrom][label]
    x = np.arange(len(y))

    fig_name = '{} - {} - aggregated by sum of max intensity'.format(
        chrom, label)
    plt.plot(x, y)
    # plt.grid()
    plt.title(fig_name)
    plt.xlabel('bins')
    plt.ylabel('sum of intensity')
    if p2f != None:
        # if not os.path.isdir(p2f):
        #     os.makedirs(p2f)
        # p2savefig = os.path.join(p2f, fig_name)
        p2savefig = p2f + '.png'
    else:
        p2savefig = 'results/sum_agg_plots/{}'.format(fig_name)
    plt.savefig(p2savefig, dpi=150)
    if not scilent:
        plt.show()
    plt.close()


def aggre_by_mean_var(df_agg_mean, df_agg_var, label, chrom='whole genome',
                      scilent=False, p2f=None):
    y = df_agg_mean.loc[chrom][label]
    x = np.arange(len(y))
    y_err = np.sqrt(df_agg_var.loc[chrom][label])

    fig_name = '{} - {} - aggregated by mean and std intensity'.format(
        chrom, label)
    plt.plot(x, y)
    plt.grid()
    plt.title(fig_name)
    plt.xlabel('bins')
    plt.ylabel('mean and std of intensity')

    plt.fill_between(x, y - y_err, y + y_err, alpha=0.5, color='grey')

    if p2f != None:
        # if not os.path.isdir(p2f):
        #     os.makedirs(p2f)
        # p2savefig = os.path.join(p2f, fig_name)
        p2savefig = p2f + '.png'
    else:
        p2savefig = 'results/sum_agg_plots/{}'.format(fig_name)
    plt.savefig(p2savefig, dpi=150)

    if not scilent:
        plt.show()
    plt.close()


def mainfn(args):
    bedgraph = args.p2bedgraph
    expr_name = args.expr_name
    chrom_size = args.p2chrom
    p2annot_loop = args.p2loop_annot
    p2loop_tag = args.p2loop_tag
    nbins = args.nbins
    p2save_dir = args.p2save_dir
    pseudo = args.pseudo

    p2bedgraph = os.path.join(bedgraph, expr_name)
    bg = BedGraph(chrom_size, p2bedgraph)

    annot_columns = ['left_chr', 'left_start', 'left_end', 'right_chr', 'right_start', 'right_end',
                     'PET count', 'left_max_intensity', 'right_max_intensity',
                     'left_max_index', 'right_max_index', 'loop_ID',
                     'left_motif_chr', 'left_motif_start', 'left_motif_end', 'left_motif_strand', 'left_distance',
                     'right_motif_chr', 'right_motif_start', 'right_motif_end', 'right_motif_strand', 'right_distance']
    df_loop = pd.read_csv(p2annot_loop, names=annot_columns, sep='\t')


    chromfile = pd.read_table(chrom_size, names=['chrom', 'size'])
    for row in chromfile.iterrows():
        chrom_name = row[1]['chrom']
        if chrom_name not in df_loop['left_chr'].values:
            continue
        bg.load_chrom_data(chrom_name)


    if pseudo == 0:
        loop_tag = pd.read_csv(
            p2loop_tag, sep='\t', index_col=0)
    elif pseudo == 1:
        # def drop_exceed_loop(df_loop_row, chromfile):
        #     chrom_name = df_loop_row['left_chr']
        #     chrom_len = chromfile[chromfile['chrom'] == chrom_name]['size'].values
        #     loop_end = df_loop_row['right_end']
        #     return loop_end > chrom_len
        # tmp_fn = lambda x: drop_exceed_loop(x, chromfile)
        # df_loop['exceed_chorm'] = df_loop.apply(tmp_fn, axis = 1)
        # df_loop = df_loop[df_loop['exceed_chorm'] == False]
        # df_loop = df_loop.reset_index()

        loop_tag = pd.DataFrame(index = df_loop.index, 
                columns = ['bias', 'convergence', 'NULL motif'])
        loop_tag.bias = 'balance'
        loop_tag.convergence = 'convergent'
        loop_tag['NULL motif'] = 'na'

    df_binned_intensity_per_loop = pd.DataFrame(
        index=df_loop.index, columns=['bias', '{} binned intensity'.format(nbins)])
    df_binned_intensity_per_loop['bias'] = loop_tag['bias']

    tmp_df = df_loop.apply(lambda x:
                           get_max_intensity_in_same_len_bins(bg,
                                                              nbins,
                                                              x.left_start,
                                                              x.left_chr,
                                                              x.right_end,
                                                              x.right_chr,
                                                              chrom_size = chromfile[chromfile['chrom'] == x.left_chr]['size']),
                           axis=1)
    df_binned_intensity_per_loop['{} binned intensity'.format(nbins)] = tmp_df

    df_binned_intensity_per_loop['convergence'] = loop_tag['convergence']
    df_binned_intensity_per_loop['NULL motif'] = loop_tag['NULL motif']
    df_binned_intensity_per_loop['chrom'] = df_loop['left_chr']

    binned_intensity_per_loop_name = 'binned_results_{}'.format(
        expr_name)
    if pseudo == 1:
        binned_intensity_per_loop_name = 'pseudo_' + binned_intensity_per_loop_name
    if not os.path.isdir(p2save_dir):
        os.makedirs(p2save_dir)
    p2binned_intensity_per_loop = os.path.join(p2save_dir, binned_intensity_per_loop_name)
    df_binned_intensity_per_loop.to_pickle(p2binned_intensity_per_loop)

    norm_df_binned_intensity_per_loop = df_binned_intensity_per_loop.copy()
    binned_name_list = []
    for name in df_binned_intensity_per_loop.columns:
        if 'binned intensity' in name:
            binned_name_list.append(name)
            norm_fn = lambda x: x / max(x)
            norm_df_binned_intensity_per_loop[name] = norm_df_binned_intensity_per_loop[name].apply(norm_fn)
    if args.norm == 0:
        df_agg_sum, df_agg_mean, df_agg_var = get_aggregated_inten_for_each_class(df_binned_intensity_per_loop,
                                                                                  nbins=nbins,
                                                                                  catag='bias')
    else:
        df_agg_sum, df_agg_mean, df_agg_var = get_aggregated_inten_for_each_class(norm_df_binned_intensity_per_loop,
                                                                                  nbins=nbins,
                                                                                  catag='bias')


    for label in df_agg_mean.columns:
        fig_name = 'norm_sum_agg_plot_{}_{}'.format(label.replace(' ', '_'), expr_name)
        if pseudo == 1:
            fig_name = 'pseudo_' + fig_name
        p2avg_fig = os.path.join(p2save_dir, 'aggregated_plots', fig_name)
        # aggre_by_mean_var(df_agg_mean, df_agg_var, label=label,
        #                   chrom='whole genome', scilent=True,
        #                   p2f=p2avg_fig)
        aggre_by_sum(df_agg_sum, label=label,
                          chrom='whole genome', scilent=True,
                          p2f=p2avg_fig)


# ----------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--p2bedgraph', type=str,
                        default='../data/intensity_bedgraph/')
    parser.add_argument('--expr_name', type=str,
                        default='SHG0180-181-182NR_hg38_cohesin.bedgraph')

    parser.add_argument('--p2chrom', type=str,
                        help='path to chromosome size file',
                        default='../../ChIA_PET/hg38.chrom.sizes')

    parser.add_argument('--p2loop_annot', type=str,
                        default='../data/LHG0052H.e500.clusters.cis.bothanchint_G250.PETcnt_G9.motifannot')
    parser.add_argument('--p2loop_tag', type=str,
                        default='../../ChIA_PET/results/loop_tags')

    parser.add_argument('--p2save_dir', type = str,
                        default = '../results/binned_intensity_per_loop/')
    parser.add_argument('--nbins', type=int,
                        default=1000)

    parser.add_argument('--pseudo', type = int,
                        choices = [0, 1], default = 0)
    parser.add_argument('--norm', type = int,
            choices = [0, 1], default = 1)

#     chrom_size = '../../ChIA_PET/hg38.chrom.sizes'
#     bedgraph = '../data/intensity_bedgraph/LHG0052H.for.BROWSER.sorted.bedgraph'

    mpl.rcParams['figure.dpi'] = 100
    args = parser.parse_args()

    mainfn(args)
