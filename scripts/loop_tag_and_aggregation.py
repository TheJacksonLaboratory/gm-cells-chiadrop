#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 qizai <jianhao2@illinois.edu>
#
# Distributed under terms of the MIT license.

"""
This is a python script file for ChIA-PET annotated region/loop aggregation.
"""

import numpy as np
import pandas as pd
import scipy
from scipy.stats import binom_test
import ipdb
import argparse
import matplotlib.pyplot as plt
from pyBedGraph import BedGraph 
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 100
import os


def main_fn(args):
    p2loop_file = args.p2loop_file
    p2bedgraph = args.p2bedgraph
    p2save_loop_tag = args.p2save_loop_tag
    nbins = args.nbins
    p2chrom_size = args.p2chrom_size
    p2binned_intensity_per_loop = args.p2binned_intensity_per_loop
    p2agg_stats = args.p2agg_stats

    annot_col_names = ['left_chr', 'left_start', 'left_end', 'right_chr', 'right_start', 'right_end',
                       'PET count', 'left_max_intensity', 'right_max_intensity',
                       'left_max_index', 'right_max_index', 'loop_ID',
                       'left_motif_chr',
                       'left_motif_start',
                       'left_motif_end',
                       'left_motif_strand',
                       'left_distance',
                       'right_motif_chr',
                       'right_motif_start',
                       'right_motif_end',
                       'right_motif_strand',
                       'right_distance']

    conv_dict = {
            '+-': 'convergence',
            '-+': 'divergence',
            '++': 'right tandem',
            '--': 'left tandem'
        }

    null_dict = {
            '.+': 'NULL-right',
            '.-': 'NULL-left',
            '-.': 'left-NULL',
            '+.': 'right-NULL',
            '..': 'NULL'
        }

    df_loop = pd.read_table(p2loop_file, 
                      names = annot_col_names)
    
    loop_tag = pd.DataFrame(columns = ['bias', 'convergence', 'NULL motif'], index = df_loop.index)
    loop_tag['bias'] = df_loop.apply(lambda x: binomial_test_fn(x.left_max_intensity, 
                                                       x.right_max_intensity), 
                            axis = 1)
    loop_tag['convergence'] = df_loop.apply(lambda x: motif_convergence_fn(x.left_motif_strand, 
                                                       x.right_motif_strand, conv_dict), 
                            axis = 1)
    loop_tag['NULL motif'] = df_loop.apply(lambda x: find_NULL_motif(x.left_motif_strand, 
                                                       x.right_motif_strand, null_dict), 
                            axis = 1)
    
    # save loop tag and added label loop annotation file.
    df_loop_new = df_loop.copy()
    df_loop_new[['bias', 'convergence', 'NULL motif']] = loop_tag[['bias', 'convergence', 'NULL motif']]
    
    loop_tag.to_csv(p2save_loop_tag, sep='\t')
    
    p2labeled_loop = p2loop_file + '_added_labels'
    df_loop_new.to_csv(p2labeled_loop, sep = '\t')
    
    whole_genome_balance_count = (loop_tag['bias'] == 'balance').sum()
    whole_genome_left_biased_count = (loop_tag['bias'] == 'left biased').sum()
    whole_genome_right_biased_count = (loop_tag['bias'] == 'right biased').sum()
    
    # aggregate bias
    chrom_list = list(set(df_loop['left_chr']).union(set(df_loop['right_chr'])))
    chrom_list.sort(key = lambda x: int(x[3:]) if x != 'chrX' else 24)
    chrom_list.append('whole genome')
    df_bias_count = pd.DataFrame(columns = ['balance_loop_count', 'balance_PET_count',
                                            'left_biased_loop_count', 'left_biased_PET_count',
                                            'right_biased_loop_count', 'right_biased_PET_count',],
                                index = chrom_list)
    
    for chrom in chrom_list[:-1]:
        chrom_loop_idx = (df_loop['left_chr'] == chrom)
        
        balance_tag_idx = (loop_tag['bias'] == 'balance')
        left_bias_tag_idx = (loop_tag['bias'] == 'left biased')
        right_bias_tag_idx = (loop_tag['bias'] == 'right biased')
    
        chrom_balance_idx = (balance_tag_idx & chrom_loop_idx)
        chrom_left_biased_idx = (left_bias_tag_idx & chrom_loop_idx)
        chrom_right_biased_idx = (right_bias_tag_idx & chrom_loop_idx)
        
        chrom_balance_count = chrom_balance_idx.sum()
        chrom_left_biased_count = chrom_left_biased_idx.sum()
        chrom_right_biased_count = chrom_right_biased_idx.sum()
        
        chrom_balance_PET = df_loop.loc[chrom_balance_idx]['PET count'].sum()
        chrom_left_biased_PET = df_loop.loc[chrom_left_biased_idx]['PET count'].sum()
        chrom_right_biased_PET = df_loop.loc[chrom_right_biased_idx]['PET count'].sum()
        
        df_bias_count.loc[chrom] = {'balance_loop_count': chrom_balance_count,
                                   'balance_PET_count': chrom_balance_PET,
                                   'left_biased_loop_count': chrom_left_biased_count,
                                   'left_biased_PET_count': chrom_left_biased_PET,
                                   'right_biased_loop_count': chrom_right_biased_count,
                                   'right_biased_PET_count': chrom_right_biased_PET}
    
    df_bias_count.loc['whole genome'] = df_bias_count.loc[chrom_list[:-1]].sum(axis = 0)
    
    df_bias_count['loop_count_proportion_blr'] = df_bias_count.apply(lambda x: count_proportion_fn(x,
                                                                                             'balance_loop_count',
                                                                                             'left_biased_loop_count',
                                                                                             'right_biased_loop_count'),
                                                               axis = 1)
    df_bias_count['PET_count_proportion_blr'] = df_bias_count.apply(lambda x: count_proportion_fn(x,
                                                                                             'balance_PET_count',
                                                                                             'left_biased_PET_count',
                                                                                             'right_biased_PET_count'),
                                                               axis = 1)
    
    
    p2df_bias_count = p2agg_stats + '_bias_count.csv'
    df_bias_count.to_csv(p2df_bias_count)
    
    # aggregate convergence results.
    conv_column_list = ['convergence_loop_count', 'convergence_PET_count',
                        'divergence_loop_count', 'divergence_PET_count',
                        'left_tandem_loop_count', 'left_tandem_PET_count',
                        'right_tandem_loop_count', 'right_tandem_PET_count']
    
    df_convergence_count = pd.DataFrame(columns = conv_column_list, index = chrom_list)
    
    for chrom in chrom_list[:-1]:
        chrom_loop_idx = (df_loop['left_chr'] == chrom)
        
        convergence_tag_idx = (loop_tag['convergence'] == 'convergence')
        divergence_tag_idx = (loop_tag['convergence'] == 'divergence')
        left_tendem_tag_idx = (loop_tag['convergence'] == 'left tandem')
        right_tendem_tag_idx = (loop_tag['convergence'] == 'right tandem')
    
        chrom_convergence_idx = (convergence_tag_idx & chrom_loop_idx)
        chrom_divergence_idx = (divergence_tag_idx & chrom_loop_idx)
        chrom_left_tendem_idx = (left_tendem_tag_idx & chrom_loop_idx)
        chrom_right_tendem_idx = (right_tendem_tag_idx & chrom_loop_idx)
        
        chrom_convergence_count = chrom_convergence_idx.sum()
        chrom_divergence_count = chrom_divergence_idx.sum()
        chrom_left_tendem_count = chrom_left_tendem_idx.sum()
        chrom_right_tendem_count = chrom_right_tendem_idx.sum()
        
        chrom_convergence_PET = df_loop.loc[chrom_convergence_idx]['PET count'].sum()
        chrom_divergence_PET = df_loop.loc[chrom_divergence_idx]['PET count'].sum()
        chrom_left_tendem_PET = df_loop.loc[chrom_left_tendem_idx]['PET count'].sum()
        chrom_right_tendem_PET = df_loop.loc[chrom_right_tendem_idx]['PET count'].sum()
        
        count_list = [chrom_convergence_count, chrom_convergence_PET,
                      chrom_divergence_count, chrom_divergence_PET,
                      chrom_left_tendem_count, chrom_left_tendem_PET,
                      chrom_right_tendem_count, chrom_right_tendem_PET]
    
        df_convergence_count.loc[chrom] = dict(zip(conv_column_list, count_list))
        
    
    df_convergence_count.loc['whole genome'] = df_convergence_count.loc[chrom_list[:-1]].sum(axis = 0)
    
    df_convergence_count['PET_count_proportion_cdlr'] = df_convergence_count.apply(
        lambda x: convergence_proportion_fn(x, 
                                            'convergence_PET_count',
                                           'divergence_PET_count',
                                           'left_tandem_PET_count',
                                           'right_tandem_PET_count'),
        axis = 1)
    
    p2df_convergence_count = p2agg_stats + '_convergence_count.csv'
    df_convergence_count.to_csv(p2df_convergence_count)
    
    # aggregate NULL motif.
    NULL_name_list = list(set(loop_tag['NULL motif']))
    NULL_name_list.sort()
    
    NULL_column_list = []
    for n in NULL_name_list:
        if n == 'na':
            continue
        NULL_column_list.append('{}_loop_count'.format(n))
        NULL_column_list.append('{}_PET_count'.format(n))
    
    df_NULL_count = pd.DataFrame(columns = NULL_column_list, index = chrom_list)
    
    for chrom in chrom_list[:-1]:
        chrom_loop_idx = (df_loop['left_chr'] == chrom)
        
        NULL_val_list = []
        for n in NULL_column_list:
            cur_type = n.split('_')[0]
            cur_tag_idx = (loop_tag['NULL motif'] == cur_type)
            
            chrom_cur_tag_idx = (cur_tag_idx & chrom_loop_idx)
            
            if n.split('_')[1] == 'loop':
                chrom_cur_count = chrom_cur_tag_idx.sum()
            elif n.split('_')[1] == 'PET':
                chrom_cur_count = df_loop.loc[chrom_cur_tag_idx]['PET count'].sum()
            
            NULL_val_list.append(chrom_cur_count)
        
        df_NULL_count.loc[chrom] = dict(zip(NULL_column_list, NULL_val_list))
        
    
    df_NULL_count.loc['whole genome'] = df_NULL_count.loc[chrom_list[:-1]].sum()
    
    loop_count_name_list = [x for x in NULL_column_list if 'loop' in x]
    df_NULL_count['loop_nn_nl_nr_ln_rn'] = df_NULL_count.apply(
        lambda x: NULL_proportion_fn(x, loop_count_name_list),
        axis = 1)
    
    PET_count_name_list = [x for x in NULL_column_list if 'PET' in x]
    df_NULL_count['PET_nn_nl_nr_ln_rn'] = df_NULL_count.apply(
        lambda x: NULL_proportion_fn(x, PET_count_name_list),
        axis = 1)
    
    p2df_NULL_count = p2agg_stats + '_NULL_motif_count.csv'
    df_NULL_count.to_csv(p2df_NULL_count)
    
    
    # READ bedgraph file and get intensity
    ipdb.set_trace()
    bg = BedGraph(p2chrom_size, p2bedgraph)
    chromfile = pd.read_table(p2chrom_size, names = ['chrom', 'size'])
    for row in chromfile.iterrows():
        bg.load_chrom_data(row[1]['chrom'])
    
    bin_name = '{} binned intensity'.format(nbins)
    df_binned_intensity_per_loop = pd.DataFrame(index = df_loop.index, columns = ['bias', bin_name])
    df_binned_intensity_per_loop['bias'] = loop_tag['bias']
    
    my_bg = bg
    tmp_df = df_loop.apply(lambda x: 
                           get_max_intensity_in_same_len_bins(my_bg, 
                                                              nbins,
                                                             x.left_start,
                                                             x.left_chr,
                                                             x.right_end,
                                                             x.right_chr),
                          axis = 1)
    df_binned_intensity_per_loop[bin_name] = tmp_df
    df_binned_intensity_per_loop['convergence'] = loop_tag['convergence']
    df_binned_intensity_per_loop['NULL motif'] = loop_tag['NULL motif']
    df_binned_intensity_per_loop['chrom'] = df_loop['left_chr']
    
    df_binned_intensity_per_loop.to_pickle(p2binned_intensity_per_loop)
    
    # aggregate intensity for each class
    # USE ChIA_Drop/2.4 from_binned_matrix_to_plot.py for aggregation plots
    # --- end of main_fn -----


def aggre_by_sum(df_agg_sum, label, chrom = 'whole genome', scilent = False, p2f = None):
    y = df_agg_sum.loc[chrom][label]
    x = np.arange(len(y))
    
    fig_name = '{} - {} - aggregated by sum of max intensity'.format(chrom, label)
    plt.plot(x, y)
    plt.grid()
    plt.title(fig_name)
    plt.xlabel('bins')
    plt.ylabel('sum of intensity')
    if p2f != None:
        if not os.path.isdir(p2f):
            os.makedirs(p2f)
        p2savefig = os.path.join(p2f, fig_name)
    else:
        p2savefig = 'results/sum_agg_plots/{}'.format(fig_name)
    plt.savefig(p2savefig, dpi = 150)
    if not scilent: 
        plt.show()
    plt.close()

def aggre_by_mean_var(df_agg_mean, df_agg_var, label, chrom = 'whole genome', 
                      scilent = False, p2f = None):
    y = df_agg_mean.loc[chrom][label]
    x = np.arange(len(y))
    y_err = np.sqrt(df_agg_var.loc[chrom][label])
    
    fig_name = '{} - {} - aggregated by mean and std intensity'.format(chrom, label)
    plt.plot(x, y)
    plt.grid()
    plt.title(fig_name)
    plt.xlabel('bins')
    plt.ylabel('mean and std of intensity')
    
    plt.fill_between(x, y - y_err, y + y_err, alpha = 0.5, color = 'grey')
    
    if p2f != None:
        if not os.path.isdir(p2f):
            os.makedirs(p2f)
        p2savefig = os.path.join(p2f, fig_name)
    else:
        p2savefig = 'results/sum_agg_plots/{}'.format(fig_name)
    plt.savefig(p2savefig, dpi = 150)

    if not scilent:
        plt.show()
    plt.close()


def get_aggregated_inten_for_each_class(df_binned_intensity_per_loop, nbins, catag):
    '''
    nbins \in {100, 500, 1000}
    catag \in {'bias', 'convergence', 'NULL motif'}
    '''
    bin_name = '{} binned intensity'.format(nbins)
    
    set_of_label = set(df_binned_intensity_per_loop[catag])
    label_list = list([x for x in set_of_label if x != 'na'])
    label_list.sort()
    
    total_num_loops_in_catag = (df_binned_intensity_per_loop[catag] != 'na').sum()
    
    chrom_list = list(set(df_binned_intensity_per_loop['chrom']))
    chrom_list.sort(key = lambda x: int(x[3:]) if x != 'chrX' else 24)
    chrom_list.append('whole genome')
    
    df_aggregate_sum = pd.DataFrame(columns = label_list, index = chrom_list)
    df_aggregate_mean = pd.DataFrame(columns = label_list, index = chrom_list)
    df_aggregate_var = pd.DataFrame(columns = label_list, index = chrom_list)
    

    for label in label_list:
        label_loop_idx = (df_binned_intensity_per_loop[catag] == label)
        for chrom in chrom_list[:-1]:
            # avoid whole genome.
            chrom_loop_idx = (df_binned_intensity_per_loop['chrom'] == chrom)
            
            tmp_df = df_binned_intensity_per_loop.loc[chrom_loop_idx & label_loop_idx]
            sum_of_intensity = tmp_df[bin_name].sum()
            mean_of_intensity = tmp_df[bin_name].mean()
            var_of_intensity = np.stack(tmp_df[bin_name]).var(axis = 0)
            
            df_aggregate_sum.loc[chrom][label] = sum_of_intensity
            df_aggregate_mean.loc[chrom][label] = mean_of_intensity
            df_aggregate_var.loc[chrom][label] = var_of_intensity
            
        df_aggregate_mean.loc['whole genome'][label] = df_binned_intensity_per_loop.loc[label_loop_idx][bin_name].mean()
        df_aggregate_sum.loc['whole genome'][label] = df_binned_intensity_per_loop.loc[label_loop_idx][bin_name].sum()
        total_var_per_label = np.stack(df_binned_intensity_per_loop.loc[label_loop_idx][bin_name]).var(axis = 0)
        df_aggregate_var.loc['whole genome'][label] = total_var_per_label



    return df_aggregate_sum, df_aggregate_mean, df_aggregate_var




def get_max_intensity_in_same_len_bins(bedGraph, nbins, left_start, chrom_left, right_end,  
                                       chrom_right = None, flank_per = 5):
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
    end_idx = right_end + flank_length
    
    nbins_edges = np.linspace(start_idx, end_idx, nbins + 1, dtype = np.int32)
    
    start_list = nbins_edges[:-1]
    end_list = nbins_edges[1:]
    
    bin_values = bedGraph.stats(start_list = start_list,
                               end_list = end_list,
                               chrom_name = chrom_left,
                               stat = 'max')
    return bin_values



def NULL_proportion_fn(NULL_count_row, list_of_column_name):
    count_list = np.array([NULL_count_row.get(x) for x in list_of_column_name]) 
    total_count = sum(count_list)
    percent_list = count_list / total_count * 100
    
    template = ['{:.0f}'] * len(percent_list)
    template = ' : '.join(template)
    
    return template.format(*list(percent_list))


def convergence_proportion_fn(convergence_count_row, conv, div, right, left):
    ccount = convergence_count_row.get(conv)
    dcount = convergence_count_row.get(div)
    lcount = convergence_count_row.get(left)
    rcount = convergence_count_row.get(right)
    template = '{:.0f} : {:.0f} : {:.0f} : {:.0f}'
    
    total_count = sum((ccount, dcount, lcount, rcount))
    cprop = ccount / total_count * 100
    dprop = dcount / total_count * 100
    lprop = lcount / total_count * 100
    rprop = rcount / total_count * 100
    
    return template.format(cprop, dprop, lprop, rprop)


def count_proportion_fn(bias_count_row, balance, left, right):
    bcount = bias_count_row.get(balance)
    lcount = bias_count_row.get(left)
    rcount = bias_count_row.get(right)
    template = '{:.0f} : {:.0f} : {:.0f}'
    
    total_count = sum((bcount, lcount, rcount))
    bprop = bcount / total_count * 100
    lprop = lcount / total_count * 100
    rprop = rcount / total_count * 100
    
    return template.format(bprop, lprop, rprop)

def find_NULL_motif(left_motif, right_motif, null_dict):
    '''
    - means <, + means >, . means NULL(*)
    *>: .+ => NULL-right
    *<: .- => NULL-left
    <*: -. => left-NULL
    >*: +. => right-NULL
    **: .. => NULL
    '''
    if left_motif == '.' or right_motif == '.':
        pattern = left_motif + right_motif
        return null_dict.get(pattern)
    return 'na'

def motif_convergence_fn(left_motif, right_motif, conv_dict):
    '''
    - means <, + means >.
    ><: +- convergence
    <>: -+ divergence
    >>: ++ right tandem
    <<: -- left tandem

     conv_dict = {
         '+-': 'convergence',
         '-+': 'divergence',
         '++': 'right tandem',
         '--': 'left tandem'
     }

    '''
    if left_motif != '.' and right_motif != '.':
        pattern = left_motif + right_motif
        return conv_dict.get(pattern)
    return 'na'


def binomial_test_fn(left_intensity, right_intensity, p = 0.5, sig = 5e-2):
    total_intensity = left_intensity + right_intensity
    p_val = binom_test(left_intensity, total_intensity, 
                       p, alternative = 'two-sided')
    if p_val > sig:
        # not significant
        bias = 'balance'
    else:
        # reject Null hypo
        if left_intensity > right_intensity:
            bias = 'left biased'
        else:
            bias = 'right biased'
    return bias


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--p2loop_file', type=str,
                        default='../../ChIA_PET/LHG0052H.e500.clusters.cis.bothanchint_G250.PETcnt_G9.motifannot')
    parser.add_argument('--p2bedgraph', type=str,
                        default='../../data/LHG0052H.for.BROWSER.sorted.bedgraph')
    parser.add_argument('--p2save_loop_tag', type=str,
                        default='loop_tag')
    parser.add_argument('--nbins', type = int,
            default = 1000)
    parser.add_argument('--p2chrom_size', type = str,
            default = '../../ChIA_PET/hg38.chrom.sizes')

    parser.add_argument('--p2agg_stats', type = str,
            default = 'ChIA_PET_region_file_aggregation_stats')
    parser.add_argument('--p2binned_intensity_per_loop', type = str,
            default = 'binned_intensity_per_loop_chia_pet')


    args = parser.parse_args()

    main_fn(args)
