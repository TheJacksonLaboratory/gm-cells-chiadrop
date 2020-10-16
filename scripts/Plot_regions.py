#!/usr/bin/env python
# coding: utf-8

# In[40]:


# The main contribution of this code is from Jianhao Peng
# This code modify the code from Jianhao

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from multiprocessing import Pool
# from tqdm import tqdm
import functools
import warnings
warnings.filterwarnings('ignore')

def get_list_of_frags(df_region,i):

    frag_left_start = df_region.loc[i]['left_start']
    frag_right_end = df_region.loc[i]['right_end']
    
    Mcount = df_region.loc[i]['Fragment_number']
    Num_mid_frags = Mcount - 2
    
    Lst_of_frags = [[df_region.loc[i]['left_start'],df_region.loc[i]['left_end']]]
    if Num_mid_frags > 0:
        list_of_frag_coord = list(map(int, df_region.loc[i]['mid_frags'].split(',')))
        for k in range(Num_mid_frags):
            Lst_of_frags.append([list_of_frag_coord[2*k],list_of_frag_coord[2*k+1]])
    Lst_of_frags.append([df_region.loc[i]['right_start'],df_region.loc[i]['right_end']])
    
    # order the fragment by starting site.
    Lst_of_frags.sort(key = lambda x: int(x[0]))
    
    Lst_of_frags_array = np.array(Lst_of_frags)
    
    return Lst_of_frags_array

def plot_parallel_lines(Mdata,list_of_gems, left_most_site = 0, right_most_site = 100,
                        left_anchor_start= 0, right_anchor_end= 100,
                        left_anchor_end = 5, right_anchor_start = 95, sorted_by_len = False,
                       title = 'ordered GEM', silent = False , saveplotpath):
    '''
    plot a list of gems in parallel,
    where each row is a gem: 
        [[frag1_start, frag1_end], [frag2_start, frag2_end], ..., [fragk_start, fragk_end]]
    so the list if 
    1. plot a straight line. thin, from start to end.
    2. plot each frag. thick, from its start to end.
    3. plot vertically left most and right most site.
    '''
    if not sorted_by_len:
        if not isinstance(list_of_gems, list):
            sorted_list_of_gems = list_of_gems.tolist()
        sorted_list_of_gems.sort(key = lambda x: x[-1][-1] - x[0][0])
    else:
        if not isinstance(list_of_gems, list):
            sorted_list_of_gems = list_of_gems.tolist()
    
    y_idx = np.arange(len(list_of_gems))[::-1]
    y_buffer = min(max(0.05, (y_idx.max() - y_idx.min()) * 0.05), 2)
    x_buffer = (right_most_site - left_most_site) * 0.05
    
    total_length = right_most_site - left_most_site
#     left_anchor_percent = (left_anchor_end - left_most_site)/total_length * 100
#     right_anchor_percent = (right_most_site - right_anchor_start)/total_length * 100
    
    height = max(1, len(y_idx) * 0.05)
    Pos = np.array(Mdata['start']+Mdata['end'])/2
    Strand = np.array(Mdata['orientation'])
    CTCF_Pet_int = np.array(Mdata['CTCF PET intensity'])
    CTCF_Drop_int = np.array(Mdata['CTCF Drop intensity'])
    Coh_Pet_int = np.array(Mdata['Coh PET intensity'])
    Coh_Drop_int = np.array(Mdata['Coh Drop intensity'])
    Ypos = np.ones(len(Mdata))
    
    if len(Pos)>0:
        if len(Pos) > 1:
            mBw = min(min(abs(Pos[0:-1]-Pos[1:])),total_length/20)
        else:
            mBw = total_length/20
        Bw = 0.9*mBw*Ypos
        fig, (ax0,axCTCF_P,axCTCF_D,axCoh_P,axCoh_D, ax1) = plt.subplots(6, 1, sharex= True,figsize = (6, height+0.2+4*0.4), dpi = 150 ,gridspec_kw={'height_ratios': [0.2,0.4,0.4,0.4,0.4, height]})
        fig.subplots_adjust(hspace=0)

    #     ax0.scatter(Pos,Ypos)
        axCTCF_P.bar(Pos,CTCF_Pet_int,Bw,color = 'b')
        axCTCF_D.bar(Pos,CTCF_Drop_int,Bw,color = 'g')
        axCoh_P.bar(Pos,Coh_Pet_int,Bw,color = 'darkorange')
        axCoh_D.bar(Pos,Coh_Drop_int,Bw,color = 'k')

        for i in range(len(Mdata)):
            txt = Mdata.iloc[i,4]
            if Strand[i] == '+':
                ax0.scatter(Pos[i],Ypos[i], marker= '>', color = 'r' )
                if txt[0:3] == 'smt':
                    NUM = txt.split('-')[-1]
                    ax0.annotate(NUM, (Pos[i], Ypos[i]),xytext = (Pos[i], Ypos[i]*0.96),ha='center', size=4)    
            elif Strand[i] == '-':
                ax0.scatter(Pos[i],Ypos[i], marker= '<', color = 'r' )
                if txt[0:3] == 'smt':
                    NUM = txt.split('-')[-1]
                    ax0.annotate(NUM, (Pos[i], Ypos[i]),xytext = (Pos[i], Ypos[i]*0.96),ha='center', size=4)

        ax0.set_yticks([])
        axCTCF_P.set_yticks([max(CTCF_Pet_int)])
        axCTCF_P.set_yticklabels([max(CTCF_Pet_int)],fontsize=4)
        axCTCF_P.set_ylabel('CTCF ChIA_PET',fontsize=4,rotation=0)
        axCTCF_D.set_yticks([max(CTCF_Drop_int)])
        axCTCF_D.set_yticklabels([max(CTCF_Drop_int)],fontsize=4)
        axCTCF_D.set_ylabel('CTCF ChIA_Drop',fontsize=4,rotation=0)
        axCoh_P.set_yticks([max(Coh_Pet_int)])
        axCoh_P.set_yticklabels([max(Coh_Pet_int)],fontsize=4)
        axCoh_P.set_ylabel('Cohesin ChIA_PET  ',fontsize=4,rotation=0)
        axCoh_D.set_yticks([max(Coh_Drop_int)])
        axCoh_D.set_yticklabels([max(Coh_Drop_int)],fontsize=4)
        axCoh_D.set_ylabel('Cohesin ChIA_Drop   ',fontsize=4,rotation=0)
    else:
        fig, (ax0,axCTCF_P,axCTCF_D,axCoh_P,axCoh_D, ax1) = plt.subplots(6, 1, sharex= True,figsize = (6, height+0.2+4*0.4), dpi = 150 ,gridspec_kw={'height_ratios': [0.2,0.4,0.4,0.4,0.4, height]})
        fig.subplots_adjust(hspace=0)

#     fig = plt.figure()
    for idx, line in enumerate(sorted_list_of_gems):
        
        left_end = line[0][0]
        right_end = line[-1][-1]
        x_thin = [left_end, right_end]
        y_thin = np.ones_like(x_thin) * y_idx[idx]
        ax1.plot(x_thin, y_thin, c = 'grey', alpha = 0.8, lw = 0.4)
        for frag in line:
            x_thick = frag
            if (x_thick[-1] - x_thick[0]) < (0.001 * total_length):
                x_thick[-1] = x_thick[0] + 0.001 * total_length 
            y_thick = np.ones_like(x_thick) * y_idx[idx]
            ax1.plot(x_thick, y_thick, c = 'b', lw = 2)
    
    # 3. add vertical lines.
    if left_anchor_start>-0.5:
        ax1.axvspan(left_anchor_start, left_anchor_end, alpha = 0.3, color = 'r')
    if right_anchor_start>-0.5:
        ax1.axvspan(right_anchor_start, right_anchor_end, alpha = 0.3, color = 'r')
    
    # set limits on x & y axis
    ax1.set_ylim( y_idx.min() - y_buffer, y_idx.max() + y_buffer)    
    ax1.set_xlim(left_most_site - x_buffer, right_most_site + x_buffer)
    
    # change x & y axis ticks and labels.
#     xlabels = ['left anchor\n({:.2f}%)'.format(left_anchor_percent),
#               'right anchor\n({:.2f}%)'.format(right_anchor_percent)]
    if left_anchor_start>-0.5 and right_anchor_start>-0.5:
        ax1.set_xticks([(left_anchor_start + left_anchor_end)/2,
                        (right_anchor_start + right_anchor_end)/2])
        ax1.set_xticklabels([(left_anchor_start + left_anchor_end)/2,
                             (right_anchor_start + right_anchor_end)/2],fontsize=4)
    elif left_anchor_start<0 and right_anchor_start>-0.5:
        ax1.set_xticks([left_most_site,
                        (right_anchor_start + right_anchor_end)/2])
        ax1.set_xticklabels([left_most_site,
                             (right_anchor_start + right_anchor_end)/2],fontsize=4)
    elif left_anchor_start>-0.5 and right_anchor_start<-0.5:
        ax1.set_xticks([(left_anchor_start + left_anchor_end)/2,
                        right_most_site])
        ax1.set_xticklabels([(left_anchor_start + left_anchor_end)/2,
                             right_most_site],fontsize=4)
    else:
        ax1.set_xticks([left_most_site,
                        right_most_site])
        ax1.set_xticklabels([left_most_site,
                             right_most_site],fontsize=4)
        
    ax1.tick_params(axis = 'y', which = 'both', left = False, top = False, 
                   labelleft = False)
    
    ax0.set_title(title, fontsize = 5)
#     saveplotpath = 'Minji_data/Cohesin_results/01ALL/4kbext_dm/Plots/'
    plt.savefig(saveplotpath+'{}.png'.format(title), 
                dpi = 600, bbox_inches='tight')
    
    if not silent:
        plt.show()
    plt.close()
        
def region_plot(Mdata,df_region_boundary, region, direction, silent = False, saveplotpath,bedfilepath,lib_name):
    '''
    GEM file is large, avoid overloading it
    it contains a GEM and its fragments in each line.
    df_region_boundary also didn't change when switching region. 
    Avoid overloading.
    '''
#     bedfilepath = 'Minji_data/Cohesin_results/01ALL/4kbext_dm/Bedfiles/'
    test_file = bedfilepath+'{}_{}.bed'.format(region, direction)
    df_region = pd.read_csv(test_file, sep = '\t', 
                           names = ['name', 'left_start', 'left_end', 'GEM_ID', 
                                    'Fragment_number','???','mid_frags', 'right_start', 'right_end'])
    if len(df_region) == 0:
#         print('region {} does not have any complex in direction: {}'.format(
#         region, direction))
        return None
    
    df_region['frag_coord_filtered_array'] = None
    
    for i in range(len(df_region)):
#         print(df_region)
        df_region.at[i,'frag_coord_filtered_array'] = get_list_of_frags(df_region,i)
    
    left_most_site, right_most_site, _ = df_region_boundary[df_region_boundary['loop_ID'] == region][['left_start', 'right_end', 'loop_ID']].values[0]
#     left_most_site -= 4000
#     right_most_site += 4000
    left_anchor_start, right_anchor_end = df_region_boundary[df_region_boundary['loop_ID'] == region][['left_motif_start', 'right_motif_end']].values[0]
    left_anchor_end, right_anchor_start = df_region_boundary[df_region_boundary['loop_ID'] == region][['left_motif_end', 'right_motif_start']].values[0]
    
    left_anchor_start -= 4000
    right_anchor_end += 4000
    left_anchor_end += 4000
    right_anchor_start -= 4000
#     print((left_most_site + left_anchor_end)/2,(right_anchor_start + right_most_site)/2)
#     print(right_anchor_start)
    Total = len(df_region['Fragment_number'])
#     lib_name = 'SHG0180-181-182NR_hg38_cohesin_FDR_0.1_ALL'
    title = lib_name+'_'+region+'_'+direction+'_'+str(Mdata.iloc[0]['chr'])+'_'+str(left_most_site)+'_'+str(right_most_site)+'_ComplexNum_'+str(Total)

#     Find corresponding motifs
    Mdata = Mdata[Mdata['end']>=left_most_site]
    Mdata = Mdata[Mdata['start']<=right_most_site]
    
#     print(Mdata)
    
    plot_parallel_lines(Mdata,df_region['frag_coord_filtered_array'], left_most_site, right_most_site, 
                        left_anchor_start= left_anchor_start, right_anchor_end= right_anchor_end,
                        left_anchor_end= left_anchor_end, right_anchor_start= right_anchor_start,
                        sorted_by_len = False, title = title, silent = silent, saveplotpath)



# t1 = time.time()
def fn(x,Mdata,df_region_boundary, direction, saveplotpath,lib_name):
    Mdata_new = Mdata[Mdata['chr'] ==  df_region_boundary['left_chr'][x]]
    region_plot(Mdata_new, df_region_boundary, 
                df_region_boundary.iloc[x,11], direction, 
                 silent = True, saveplotpath,bedfilepath,lib_name)
    return None



# Thread = 0
# Chuck = int(len(df_region_boundary)/10)
# # for direction in ['Left', 'Right', 'Both', 'None']:
# #     with Pool(8) as p:
# # #         r = list(tqdm.tqdm(p.imap(_foo, range(30)), total=30))
# #         r = list(tqdm(p.imap(fn, range(100)), total = 100))
# for direction in ['Left']:
#     test_fn = functools.partial(fn,Mdata = Mdata, df_region_boundary = df_region_boundary,direction = direction)
#     for i in tqdm(range(0,1)):
#         test_fn(i)
# t2 = time.time()
# print('pass: {:2f}s'.format(t2 - t1))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--PoolNum',type = int)
    parser.add_argument('--p2region',type = str)
    parser.add_argument('--Motifpath',type = str)
    parser.add_argument('--saveplotpath',type = str)
    parser.add_argument('--bedfilepath',type = str)
    parser.add_argument('--lib_name',type = str)

    args = parser.parse_args()
    
    PoolNum = args.PoolNum
    p2region = args.p2region # path of files for loops/domains
    Motifpath = args.Motifpath # path of Motif file
    saveplotpath = args.saveplotpath # path for saving plots (i.e. ./Plots/)
    bedfilepath = args.bedfilepath # path for reading bedfiles, a directory (i.e. ./Bedfiles/)
    lib_name = args.lib_name # name for the library (i.e. SHG0180-181-182NR_hg38_cohesin_FDR_0.1_ALL)
#     p2region = 'Minji_data/LHG0052H.e500.clusters.cis.bothanchint_G250.PETcnt_G9.motifannot.sorted.domains'
    df_region_boundary = pd.read_csv(p2region, sep = '\t', 
                                    names = ['left_chr', 'left_start', 'left_end', 'right_chr', 
                                             'right_start', 'right_end','PET count', 'left_max_intensity',
                                             'right_max_intensity', 'left_max_index', 'right_max_index', 'loop_ID',
                                             'left_motif_chr', 'left_motif_start', 'left_motif_end', 'left_motif_strand',
                                             'left_distance','right_motif_chr', 'right_motif_start', 'right_motif_end',
                                             'right_motif_strand', 'right_distance'])
    # Motifpath = 'Minji_data/CTCF_motifs_STORM_hg38_Ext4kbBoth_chiapet_chiadrop_intensity_id.bed'
#     Motifpath = 'Minji_data/CTCF_motifs_STORM_hg38_Ext4kbBoth_with_supermotif_id.bed'
    Mdata = pd.read_csv(Motifpath, sep = '\t',
                        names = ['chr','start','end','orientation','Mt_ID',
                             'dmID','Side',
                             'CTCF PET intensity','CTCF Drop intensity',
                             'Coh PET intensity','Coh Drop intensity'])
    
    for direction in ['Left', 'Right']:
        test_fn = functools.partial(fn,Mdata = Mdata,
                                    df_region_boundary = df_region_boundary,
                                    direction = direction,
                                    saveplotpath=saveplotpath,bedfilepath=bedfilepath,lib_name=lib_name)
        with Pool(PoolNum) as p:
            p.map(test_fn, range(len(df_region_boundary)))
# t2 = time.time()
# print('pass: {:2f}s'.format(t2 - t1))

