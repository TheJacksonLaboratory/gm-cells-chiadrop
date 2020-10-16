# The main contribution of this code is from Jianhao Peng
# This code modify the code from Jianhao

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from multiprocessing import Pool
from tqdm import tqdm
import functools
import warnings
import glob
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

def plot_each_block(list_of_gems,fig,ax,total_length,Colorcode,sorted_by_len = False):
    if not sorted_by_len:
        if not isinstance(list_of_gems, list):
            sorted_list_of_gems = list_of_gems.tolist()
        sorted_list_of_gems.sort(key = lambda x: x[-1][-1] - x[0][0])
    else:
        if not isinstance(list_of_gems, list):
            sorted_list_of_gems = list_of_gems.tolist()
    
    y_idx = np.arange(len(list_of_gems))[::-1]
    num_frags = len(y_idx)
    height = len(y_idx) * 0.05
    
#     ipdb.set_trace()
    for idx, line in enumerate(sorted_list_of_gems):
        left_end = line[0][0]
        right_end = line[-1][-1]
        x_thin = [left_end, right_end]
        y_thin = np.ones_like(x_thin) * y_idx[idx]
        ax.plot(x_thin, y_thin, c = 'grey', alpha = 0.8, lw = 0.4)
        for frag in line:
            x_thick = frag
            if (x_thick[-1] - x_thick[0]) < (0.001 * total_length):
                x_thick[-1] = x_thick[0] + 0.001 * total_length 
            y_thick = np.ones_like(x_thick) * y_idx[idx]
            # If CTCF, c = #0000B2. If Cohesin, c = #005900. If RNAPII, c = #590059
            ax.plot(x_thick, y_thick, c = Colorcode, lw = 2)
    
    ax.tick_params(axis = 'y', which = 'both', left = False, top = False,labelleft = False)
    ax.set_ylabel('\n\n C#:'+str(num_frags),fontsize=6,rotation='horizontal', ha='right')
    ax.set_ylim(ymin=-1.25, ymax=len(y_idx)+0.75)
    return num_frags
        
def region_plot(region_boundary,bedfilepath,LoopName,Loop, saveplotpath,LoopInfo,RAIDLen,Colorcode, silent = True):
    '''
    GEM file is large, avoid overloading it
    it contains a GEM and its fragments in each line.
    df_region_boundary also didn't change when switching region. 
    Avoid overloading.
    '''
#     ipdb.set_trace()
    Typelist = ['Left','Right','Both','LoL','RoR']
    total_length = region_boundary[-1]-region_boundary[0]
    df_region = [[] for dummy in range(6)]
    Heights = np.zeros([6])
    
    for kk in range(1,6):
        Type = Typelist[kk-1]
        df_region_path = bedfilepath+'_'+Type+'.bed'
#         ipdb.set_trace()
        df_region[kk] = pd.read_csv(df_region_path, sep = '\t', 
                                            names = ['name', 'left_start', 'left_end', 'GEM_ID', 
                                                    'Fragment_number','mid_frags', 'right_start',
                                                     'right_end'])
        
        df_region[kk]['frag_coord_filtered_array'] = None
        for i in range(len(df_region[kk])):
            df_region[kk].at[i,'frag_coord_filtered_array'] = get_list_of_frags(df_region[kk],i)
        
        Heights[kk] = max(len(np.arange(len(df_region[kk]['frag_coord_filtered_array']))[::-1])* 0.05+0.1,0.25)        
    
    Heights[0] = 0.4
    
    
    fig, ax = plt.subplots(6, 1, sharex= True,
                           figsize = (6, np.sum(Heights)),
                           gridspec_kw={'height_ratios': Heights})
    
    
    ax[0].get_yaxis().set_ticks([])
#     Add loop info here
    ax[0].set_xlim(xmin=region_boundary[0], xmax=region_boundary[-1])
    ax[0].set_xticks([region_boundary[0],region_boundary[-1]])
    ax[0].set_xticklabels([region_boundary[0],region_boundary[-1]],fontsize=4)
    
    Num_frags = np.zeros([6])
    for kk in range(1,6):
        Num_frags[kk] = plot_each_block(df_region[kk]['frag_coord_filtered_array'],fig,ax[kk],total_length,Colorcode,sorted_by_len = False)
    
    
#     ipdb.set_trace()
    if Loop.iloc[6][0:5] == 'CRNAP':
        marker = 's'
        color = 'r'
    elif Loop.iloc[6][0:5] == 'CRNEN':
        marker = 'o'
        color = 'yellow'
    elif Loop.iloc[6][0:5] == 'CRNWP':
        marker = 's'
        color = 'darkorange'
    elif Loop.iloc[6][0:5] == 'CRNOT':
        marker = 'hexagon2'
        color = 'b'
    elif Loop.iloc[6][0:5] == 'CRNPE':
        marker = '^'
        color = 'g' 
    ax[0].scatter((Loop.iloc[1]+Loop.iloc[2])/2,0.2, marker= marker, color = color )
    ax[0].annotate(LoopInfo[0], ((Loop.iloc[1]+Loop.iloc[2])/2,0.2),xytext = ((Loop.iloc[1]+Loop.iloc[2])/2, 0.275),ha='center', size=4)
    
    if Loop.iloc[7][0:5] == 'CRNAP':
        marker = 's'
        color = 'r'
    elif Loop.iloc[7][0:5] == 'CRNEN':
        marker = 'o'
        color = 'yellow'
    elif Loop.iloc[7][0:5] == 'CRNWP':
        marker = 's'
        color = 'darkorange'
    elif Loop.iloc[7][0:5] == 'CRNOT':
        marker = 'hexagon2'
        color = 'b'
    elif Loop.iloc[7][0:5] == 'CRNPE':
        marker = '^'
        color = 'g' 
    ax[0].scatter((Loop.iloc[4]+Loop.iloc[5])/2,0.2, marker= marker, color = color )
    ax[0].annotate(LoopInfo[1], ((Loop.iloc[4]+Loop.iloc[5])/2,0.2),xytext = ((Loop.iloc[4]+Loop.iloc[5])/2, 0.075),ha='center', size=4)
    
    
    title = LoopName
    ax[0].set_title(title+'('+str(RAIDLen)+'bp)', fontsize = 10)
    fig.subplots_adjust(hspace=0)
    
    tmplist=[]
    for num in Num_frags[1:]:
        tmplist.append(str(int(num)))
    
    plt.savefig(saveplotpath+'GM12878-cohesin-pooled_comp_FDR_0.2_PASS_'+title+'_C_'+'_'.join(tmplist)+'.pdf', bbox_inches='tight')
    
    if not silent:
        plt.show()
    plt.close()
         
    return Num_frags[1:]

def mainfunc(bedfiledir,saveplotpath,regionfilepath,RAIDfilepath,LoopInfopath,
                 savecsvpath,Colorcode):
# bedfiledir: path for sorted bedfiles
# saveplotpath: path to save plots
# regionfilepath: path for the region file (i.e. 'CRN-loops.G100kb.20200705.bedpe')
# RAIDfilepath: path for RAID file (i.e. 'CRN-loops.G100kb.RAID-annot.20200705.bed')
# LoopInfopath: path for LoopInfo file (i.e. 'GM12878_ChIA-PET_no-CTCF_yes-cohesin_yes-RNAPII_yes-NIPBL_peaks_chromHMM_RNAseq_rep1_transcriptquant_ENCFF879KFK_annot_20200708.bed')
# savecsvpath: path for saving stats.csv
# Colorcode: colorcode for fragments
    header_list = ['First_chr', 'First_start','First_end',
               'Second_chr','Second_start','Second_end',
               'First_Name','Second_Name','RAID_Name']
    regionfile = pd.read_csv(regionfilepath,sep = '\t',names=header_list)
    RAIDfile = pd.read_csv(RAIDfilepath,sep = '\t',names = ['chr','start','end','Name','??'])
    LoopInfofile = pd.read_csv(LoopInfopath,sep = '\t',names = ['chr','start','end','Name','+',
                                                                'chr+','chr+start','chr+end','gene+',
                                                                'TPM+','-',
                                                                'chr-','chr-start','chr-end','gene-',
                                                                'TPM-'])

    Record = []
    for i in tqdm(range(len(regionfile))):
    #     ipdb.set_trace()
        Loop = regionfile.iloc[i]
        RAID = RAIDfile[RAIDfile['Name'] == Loop.iloc[8]]
        Lmostsite = Loop['First_start']-Loop['Second_start']+Loop['First_end']
        Rmostsite = Loop['Second_end']+Loop['Second_start']-Loop['First_end']
        Dummy = np.array([RAID['start'].iloc[0],RAID['end'].iloc[0]]).astype(int)
        RAIDLen = Dummy[-1]-Dummy[0]
        region_boundary = np.array([Lmostsite,Rmostsite]).astype(int)
    #     ipdb.set_trace()

        LoopInfo = []
        Temp = LoopInfofile[LoopInfofile['Name']==Loop.iloc[6]][['gene-','TPM-','Name','gene+','TPM+']]
        if len(Temp)>0:
            TempList = []
            for key in ['gene-','TPM-','Name','gene+','TPM+']:
                TempList.append(str(Temp[key].iloc[0]))

            LoopInfo.append(';'.join(TempList))
        else:
            LoopInfo.append('.;0;'+Loop.iloc[6]+';0;.')

        Temp = LoopInfofile[LoopInfofile['Name']==Loop.iloc[7]][['gene-','TPM-','Name','gene+','TPM+']]
        if len(Temp)>0:
            TempList = []
            for key in ['gene-','TPM-','Name','gene+','TPM+']:
                TempList.append(str(Temp[key].iloc[0]))
            LoopInfo.append(';'.join(TempList))
        else:
            LoopInfo.append('.;0;'+Loop.iloc[7]+';0;.')


    #     ipdb.set_trace()
        NumFrags = fn(Loop,region_boundary,LoopInfo,RAIDLen,Colorcode)
        Temp = []
        Temp.append(RAID['Name'].iloc[0])
        Temp.append(RAID['chr'].iloc[0]+':'+str(RAID['start'].iloc[0])+'-'+str(RAID['end'].iloc[0]))
        Temp.append(Loop.iloc[6])
        Temp.append(Loop.iloc[7])
        for num in NumFrags:
            Temp.append(int(num))

        Record.append(Temp)

    DF = pd.DataFrame(Record,columns = ['RAID ID','RAID coordinate','Left peak','Right peak',
                                '# of complexes for plot1','# of complexes for plot2',
                               '# of complexes for plot3','# of complexes for plot4',
                               '# of complexes for plot5'])
    DF.to_csv(savecsvpath+'Record.csv',index=False)    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bedfiledir',type = str)
    parser.add_argument('--saveplotpath',type = str)
    parser.add_argument('--regionfilepath',type = str)
    parser.add_argument('--RAIDfilepath',type = str)
    parser.add_argument('--LoopInfopath',type = str)
    parser.add_argument('--savecsvpath',type = str)
    parser.add_argument('--Colorcode',type = str)
    
    args = parser.parse_args()
    
    bedfiledir = args.bedfiledir
    saveplotpath = args.saveplotpath
    regionfilepath = args.regionfilepath
    savecsvpath = args.savecsvpath
    RAIDfilepath = args.RAIDfilepath
    LoopInfopath = args.LoopInfopath
    Colorcode = args.Colorcode

    mainfunc(bedfiledir,saveplotpath,regionfilepath,RAIDfilepath,LoopInfopath,
                 savecsvpath,Colorcode)