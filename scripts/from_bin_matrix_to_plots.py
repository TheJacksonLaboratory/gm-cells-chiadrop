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


def aggre_by_avg_with_pseudo(df_agg_sum, df_agg_sum_pseudo,
        label, chrom='whole genome', scilent=False, p2f=None, color = 'b', linestyle = '-', 
        y_min = 0, y_max = 1):
    y = df_agg_sum.loc[chrom][label]
    x = np.arange(len(y))

    y_pseudo = df_agg_sum_pseudo.loc[chrom]['balance']
    x_pseudo = np.arange(len(y_pseudo))

    len_of_y = len(y)
    left_peak = max(y[:int(len_of_y * 0.2)])
    right_peak = max(y[int(len_of_y * 0.8):])
    # p_val = binom_test([left_peak, right_peak])

    fig_name = '{} - {} - aggregated by sum of max intensity'.format(
        chrom, label)
    plt.plot(x, y, color = color, linestyle = linestyle, lw = 4)

    plt.plot(x_pseudo, y_pseudo, color = '#606060', lw = 4)
    # plt.grid()
    plt.title(fig_name)
    plt.xlabel('bins')
    plt.ylabel('average intensity')
    plt.ylim(y_min, y_max)
    if p2f != None:
        # if not os.path.isdir(p2f):
        #     os.makedirs(p2f)
        # p2savefig = os.path.join(p2f, fig_name)
        p2savefig = p2f + '_left_peak_{:.2f}_right_peak_{:.2f}'.format(left_peak, right_peak) +'.png'
    else:
        p2savefig = 'results/sum_agg_plots/{}'.format(fig_name)
    plt.savefig(p2savefig, dpi=150)
    if not scilent:
        plt.show()
    plt.close()

def aggre_by_avg(df_agg_sum, label, chrom='whole genome', scilent=False, p2f=None, color = 'b', linestyle = '-', 
        y_min = 0, y_max = 1):
    y = df_agg_sum.loc[chrom][label]
    x = np.arange(len(y))

    len_of_y = len(y)
    left_peak = max(y[:int(len_of_y * 0.2)])
    right_peak = max(y[int(len_of_y * 0.8):])
    # p_val = binom_test([left_peak, right_peak])

    fig_name = '{} - {} - aggregated by sum of max intensity'.format(
        chrom, label)
    plt.plot(x, y, color = color, linestyle = linestyle)
    # plt.grid()
    plt.title(fig_name)
    plt.xlabel('bins')
    plt.ylabel('average intensity')
    plt.ylim(y_min, y_max)
    if p2f != None:
        # if not os.path.isdir(p2f):
        #     os.makedirs(p2f)
        # p2savefig = os.path.join(p2f, fig_name)
        p2savefig = p2f + '_left_peak_{:.2f}_right_peak_{:.2f}'.format(left_peak, right_peak) +'.png'
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
    expr_name = args.expr_name
    nbins = args.nbins
    p2save_dir = args.p2save_dir

    color_dict = {
            'CTCF_ChIA_PET': '#0000FF',
            'CTCF_ChIA_Drop': '#0000B2',
            'SMC1A_ChIA_PET': '#008000', # 'Cohesin_ChIA_PET': '#008000', 
                'Cohesin_ChIP_seq': '#008000',
            'Cohesin_ChIA_Drop': '#005900',
                'SMC1A_ChIA_Drop': '#005900',
                'RAD21_ChIA_Drop': '#005900',

            'RNAPII_ChIA_PET': '#800080', # 'RNAPII_ChIA_Drop': '#590059',
            'NIPBL_ChIP_seq': '#FFA500',
            'RAD21_ChIP_seq': '#008000',
            'MED1_ChIP_seq': '#A52A2A', 
                'unknown_1': '#0000B2',
                'unknown_2': '#0000B2',

            'CTCF_ChIA_Drop_new': '#0000B2',
            'Cohesin_ChIA_Drop_new': '#005900',
            'RNAPII_ChIA_Drop_new' : '#590059',
            'RNAPII_ChIA_PET_new' : '#800080',
            'RAD21_ChIA_PET_new' : '#008000'
            
            }
    expr_dict = {
            'CDH0002NR_hg38_CTCF.bedgraph': 'CTCF_ChIA_Drop',
            'LHG0035V.for.BROWSER.sorted.bedgraph': 'RNAPII_ChIA_PET',
            'LHG0051H.for.BROWSER.sorted.bedgraph': 'SMC1A_ChIA_PET',
            'LHG0052H.for.BROWSER.sorted.bedgraph': 'CTCF_ChIA_PET',
            # 'SHG0180-181-182NR_hg38_cohesin.bedgraph': 'unknown',
            'selected_ChIP-Seq.Rad21.ENCFF916YVI.hg38.sorted.bedgraph': 'RAD21_ChIP_seq',
            'selected_ChIP-Seq.MED1.rep1.GSM2443457.SRR5139372.rmdup.q30.hg38.sorted.bedgraph': 'MED1_ChIP_seq',
            'selected_ChIP-Seq.NIPBL.rep1.GSM2443453.SRR5139368.rmdup.q30.hg38.sorted.bedgraph': 'NIPBL_ChIP_seq',
                'selected_ChIP-Seq.SMC1A.rep1.GSM2443455.SRR5139370.rmdup.q30.hg38.sorted.bedgraph': 'Cohesin_ChIP_seq',
            'CDH0002NR9L_hg38_CTCF_filt_comp.bedgraph': 'CTCF_ChIA_Drop',
            'CDH0004NR_hg38_cohesin_filt_comp.bedgraph': 'Cohesin_ChIA_Drop',
                'SHG0180-181-195-196NR_hg38_RAD21_filt_comp.bedgraph': 'RAD21_ChIA_Drop',
                'SHG0182-186-197NR_hg38_SMC1A_filt_comp.bedgraph': 'SMC1A_ChIA_Drop',
                'SHG8113-8132NR_hg38_CTCF_filt_comp.bedgraph': 'unknown_1',
                'SHG8117-8133NR_hg38_CTCF_filt_comp.bedgraph': 'unknown_2',

            'GM12878-CTCF-pooled.bedgraph': 'CTCF_ChIA_Drop_new',
            'GM12878-cohesin-pooled.bedgraph': 'Cohesin_ChIA_Drop_new',
            'GM12878-RNAPII-pooled.bedgraph' : 'RNAPII_ChIA_Drop_new',
            'LHG0035N_0035V_0045V.for.BROWSER.sorted.bedgraph' : 'RNAPII_ChIA_PET_new',
            'LHG0104V.for.BROWSER.sorted.bedgraph' : 'RAD21_ChIA_PET_new'
            }
    if expr_name not in expr_dict:
        print(expr_name, ' is not specified in the slides with type of expr')
        return None
    else:
        expr_tag = expr_dict[expr_name]
        if expr_tag not in color_dict:
            print('{} belongs to {}, but the color of it is not specified, use yellow dash line'.format(expr_name, expr_tag))
            color = 'y'
            linestyle = '-.'
        else:
            color = color_dict[expr_tag]
            linestyle = '-'


    binned_intensity_per_loop_name = 'binned_results_{}'.format(
        expr_name)
    p2binned_intensity_per_loop = os.path.join(p2save_dir, binned_intensity_per_loop_name)
    p2binned_intensity_per_loop_pseudo = os.path.join(p2save_dir, 'pseudo_' + binned_intensity_per_loop_name)


    df_binned_intensity_per_loop = pd.read_pickle(p2binned_intensity_per_loop)
    df_binned_intensity_per_loop_pseudo = pd.read_pickle(p2binned_intensity_per_loop_pseudo)

    # Select only the convergent loops.
    df_binned_intensity_per_loop = df_binned_intensity_per_loop[df_binned_intensity_per_loop.convergence == 'convergence']
    df_binned_intensity_per_loop_pseudo = df_binned_intensity_per_loop_pseudo[df_binned_intensity_per_loop_pseudo.convergence == 'convergent']
    
    norm_df_binned_intensity_per_loop = df_binned_intensity_per_loop.copy()
    norm_df_binned_intensity_per_loop_pseudo = df_binned_intensity_per_loop_pseudo.copy()

    binned_matrix_col = '{} binned intensity'.format(nbins)
    max_of_each_row_in_observed_loop = df_binned_intensity_per_loop[binned_matrix_col].apply(max).mean()
    max_of_each_row_in_pseudo_loop = df_binned_intensity_per_loop_pseudo[binned_matrix_col].apply(max).mean()
    norm_factor = max(max_of_each_row_in_pseudo_loop, max_of_each_row_in_observed_loop)
    binned_name_list = []

    for name in df_binned_intensity_per_loop.columns:
        if 'binned intensity' in name:
            binned_name_list.append(name)
            norm_fn = lambda x: x / norm_factor
            norm_df_binned_intensity_per_loop[name] = norm_df_binned_intensity_per_loop[name].apply(norm_fn)
            norm_df_binned_intensity_per_loop_pseudo[name] = norm_df_binned_intensity_per_loop_pseudo[name].apply(norm_fn)

    if args.norm == 0:
        df_agg_sum, df_agg_mean, df_agg_var = get_aggregated_inten_for_each_class(df_binned_intensity_per_loop,
                                                                                  nbins=nbins,
                                                                                  catag='bias')
    else:
        df_agg_sum, df_agg_mean, df_agg_var = get_aggregated_inten_for_each_class(norm_df_binned_intensity_per_loop,
                                                                                  nbins=nbins,
                                                                                  catag='bias')
        df_agg_sum_pseudo, df_agg_mean_pseudo, df_agg_var_pseudo = get_aggregated_inten_for_each_class(norm_df_binned_intensity_per_loop_pseudo,
                                                                                  nbins=nbins,
                                                                                  catag='bias')

    number_of_region_per_label = {}
    y_max = max(df_agg_mean.loc['whole genome']['left biased'].max(),
            df_agg_mean.loc['whole genome']['right biased'].max(),
            df_agg_mean.loc['whole genome']['balance'].max(),
            df_agg_mean_pseudo.loc['whole genome']['balance'].max())
    y_min = min(df_agg_mean.loc['whole genome']['left biased'].min(),
            df_agg_mean.loc['whole genome']['right biased'].min(),
            df_agg_mean.loc['whole genome']['balance'].min(),
            df_agg_mean_pseudo.loc['whole genome']['balance'].min())

    for label in df_agg_mean.columns:
        num_region = (df_binned_intensity_per_loop.bias == label).sum()
        number_of_region_per_label[label] = num_region

        fig_name = '{}_norm_mean_agg_plot_{}_{}_num_of_region_{}'.format(expr_tag, label.replace(' ', '_'), expr_name, num_region)
        p2avg_fig = os.path.join(p2save_dir, 'colored_aggregated_plots', fig_name)
        # aggre_by_mean_var(df_agg_mean, df_agg_var, label=label,
        #                   chrom='whole genome', scilent=True,
        #                   p2f=p2avg_fig)
        aggre_by_avg_with_pseudo(df_agg_mean, df_agg_mean_pseudo,
                          label=label,
                          chrom='whole genome', scilent=True,
                          p2f=p2avg_fig, color = color, linestyle = linestyle,
                          y_max = y_max * 1.1, y_min = y_min * 0.9)

        if label == 'balance':
            num_region_pseudo = (df_binned_intensity_per_loop_pseudo.bias == label).sum()
            fig_name_pseudo = '{}_pseudo_norm_mean_agg_plot_{}_{}_num_of_region_{}'.format(expr_tag, label.replace(' ', '_'), 
                    expr_name, num_region_pseudo)

            p2avg_fig = os.path.join(p2save_dir, 'colored_aggregated_plots', fig_name_pseudo)
            aggre_by_avg(df_agg_mean_pseudo, label = label,
                    chrom = 'whole genome', scilent = True,
                    p2f = p2avg_fig, color = color, linestyle = linestyle,
                    y_max = y_max * 1.1, y_min = y_min * 0.9)


# ----------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--expr_name', type=str,
                        default='SHG0180-181-182NR_hg38_cohesin.bedgraph')

    parser.add_argument('--p2save_dir', type = str,
                        default = '../results/binned_intensity_per_loop/')
    parser.add_argument('--nbins', type=int,
                        default=1000)

    parser.add_argument('--norm', type = int,
            choices = [0, 1], default = 1)

#     chrom_size = '../../ChIA_PET/hg38.chrom.sizes'
#     bedgraph = '../data/intensity_bedgraph/LHG0052H.for.BROWSER.sorted.bedgraph'

    mpl.rcParams['figure.dpi'] = 100
    args = parser.parse_args()

    mainfn(args)
