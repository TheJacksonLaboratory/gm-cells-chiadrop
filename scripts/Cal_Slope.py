#!/usr/bin/env python
# coding: utf-8

# In[28]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
from pybedtools import BedTool
import ipdb
import time
import argparse
from tqdm import tqdm

def get_list_of_fragsmp(df_region,i):

    frag_left_start = df_region.loc[i]['Lm_frag_start']
    frag_right_end = df_region.loc[i]['Rm_frag_end']
    
    Mcount = df_region.loc[i]['FragNum']
    Num_mid_frags = Mcount - 2
    
    Lst_of_frags = [(int(df_region.loc[i]['Lm_frag_start'])+int(df_region.loc[i]['Lm_frag_end']))/2]
    if Num_mid_frags > 0:
        list_of_frag_coord = list(map(int, df_region.loc[i]['Mid_frags'].split(',')))
        for k in range(Num_mid_frags):
            Lst_of_frags.append((int(list_of_frag_coord[2*k])+int(list_of_frag_coord[2*k+1]))/2)
    Lst_of_frags.append((int(df_region.loc[i]['Rm_frag_start'])+int(df_region.loc[i]['Rm_frag_end']))/2)
    return list(Lst_of_frags)

def Mt_inInterval(Mfile,start,end,option = 0):
#     Use mp as criterion
    if option == 0:
        Mfile = Mfile[(start<=(Mfile['M_start']+Mfile['M_end'])/2)&(end>=(Mfile['M_start']+Mfile['M_end'])/2)]
        return Mfile
    else:
#         ipdb.set_trace()
        Mfile = Mfile[(np.minimum(Mfile['M_end'],end)-np.maximum(Mfile['M_start'],start))>0]
        return Mfile

# print(get_list_of_fragsmp(Bfile,0))
# # Bfile.insert(9,'Fragmp',None)
# Bfile.at[0,'Fragmp'] = get_list_of_fragsmp(Bfile,0)
# Bfile.head()


# In[ ]:


def mainfunc(pathR,pathM,pathM2,RegInterval,pathB,saveBGratepath,saveMTratepath):
# pathR: path for region file (i.e. Minji_data/LHG0052H.e500.clusters.cis.bothanchint_G250.PETcnt_G9.motifannot)
# pathM: path for Motif file (i.e. Minji_data/CTCF_motifs_STORM_hg38_Ext4kbBoth_with_supermotif_domain_id_v3.bed)
# RegInterval: range(Start_pos,End_pos)
# pathB: path for directory of bedfiles extracted from data (i.e. Minji_data/Cohesin_results/01ALL/4kbext_dm/Bedfiles/)

# saveBGratepath: path for saving the background_rate table (i.e. Minji_data/Cohesin_results/01ALL/4kbext_dm/Slope/Background_rate.csv)
# saveMTratepath: path for saving the Motif_rate table (i.e. Minji_data/Cohesin_results/01ALL/4kbext_dm/Slope/Motif_rate.csv)
# Note that if you want to parallel this process, give different name for .csv files of each Thread.
# Then use MergeCSV.py to merge them (each directory should contain only 1 type of .csv file)

    # pathR = 'Minji_data/LHG0052H.e500.clusters.cis.bothanchint_G250.PETcnt_G9.motifannot'
    Region = BedTool(pathR)
    # Remove unnecessary entries
    Region_short = Region.groupby(g=[1,2,6,12,14,20,8,9,16,21], c=[12], o=['count'])

    # pathM = 'Minji_data/CTCF_motifs_STORM_hg38_Ext4kbBoth_with_supermotif_domain_id_v3.bed'
    Mfile = pd.read_csv(pathM, sep = '\t',names = ['chr','M_start','M_end','Sign','M_name',
                                                  'dmID','Side',
                                                  'CTCF_Pet_int','CTCF_Drop_int','Coh_Pet_int','Coh_Drop_int'])
    
#     pathM2 = 'Minji_data/CTCF_motifs_STORM_hg38_Ext4kbBoth.sorted.id.bed'
    M2file = pd.read_csv(pathM2, sep = '\t',names = ['chr','M_start','M_end','Sign','M_name'])
    
    
    List = []
    List_mt = []
    window_size = 8000

    for i in RegInterval:
    # for i in tqdm(range(7,8)):
        # print(Region_short[i])
        total_length = int(Region_short[i][2])-int(Region_short[i][1])
        # print(total_length)
        # Count how many slicing should be
        Num_window = int(np.ceil(total_length/window_size))
        lpID = Region_short[i][3]
    #     print(lpID)

        R_start = int(Region_short[i][1])
        R_end = int(Region_short[i][2])

        for direction in ['Left','Right']:
            
#             pathB = 'Minji_data/Cohesin_results/01ALL/4kbext_dm/Bedfiles/'
            pathB_all = pathB+'{}_{}.bed'.format(lpID,direction)
            Bfile = pd.read_csv(pathB_all, sep = '\t',names = ['chr','Lm_frag_start','Lm_frag_end','GEMID','FragNum',
                                                          '???','Mid_frags','Rm_frag_start','Rm_frag_end'])
            # Bfile preprocessing
            Bfile.insert(9,'Fragmp',None)
            for k in range(len(Bfile)):
                Bfile.at[k,'Fragmp'] = get_list_of_fragsmp(Bfile,k)
        #     ipdb.set_trace()

    #         #         Do for Motif_rate
            Lmp_start = int(Region_short[i][1])
            Rmp_end = int(Region_short[i][2])
        #         Find the corresponding mt
            RegMfile = Mfile[Mfile['chr'] == Region_short[i][0]]
            RegM2file = M2file[M2file['chr'] == Region_short[i][0]]
    #         ipdb.set_trace()
            TMP = RegMfile[RegMfile['dmID']==Region_short[i][3]]
            Lmt = TMP[TMP['Side']=='L']
    #         Rmt = RegMfile[RegMfile['dmID']==Region_short[i][3]]
            Rmt = TMP[TMP['Side']=='R']
    #         ipdb.set_trace()
            if len(Lmt)>0:
                if len(Rmt)>0:
                    TempMt = Mt_inInterval(RegMfile,int(Lmt['M_start']),int(Rmt['M_start'])).sort_values(by = ['M_start'])
                else:
                    TempMt = Mt_inInterval(RegMfile,int(Lmt['M_start']),Rmp_end).sort_values(by = ['M_start'])
            else:
                if len(Rmt)>0:
                    TempMt = Mt_inInterval(RegMfile,Lmp_start,int(Rmt['M_start'])).sort_values(by = ['M_start'])
                else:
                    TempMt = Mt_inInterval(RegMfile,Lmp_start,Rmp_end).sort_values(by = ['M_start'])

            if direction == 'Left' and len(Lmt)>0:
            #     Do for back_ground, use binning
                Count = 0
                Table = np.zeros(Num_window)
                for k in range(len(Bfile)):
                    Fragmp = np.array(Bfile.loc[k,'Fragmp'])
                    Idx = np.minimum((np.maximum(np.array((Fragmp - R_start)/window_size),0)),Num_window-1).astype(int)
                    Table[Idx] += 1
                for j in range(Num_window):
                    W_start = R_start + window_size*j
                    W_end = R_start + window_size*(j+1)-1
                    if j == Num_window:
                        W_end = R_end
                    TempMfile = Mt_inInterval(RegM2file,W_start,W_end,1)
                    Count = int(Table[j])

                    if len(TempMfile)>0:
                        for l in range(len(TempMfile)):
                            List.append([Region_short[i][3],Region_short[i][0],W_start,W_end,
                                         TempMfile.iloc[l,0]+':'+str(TempMfile.iloc[l,1])
                                         +'-'+str(TempMfile.iloc[l,2])+','+TempMfile.iloc[l,3]+','
                                         +TempMfile.iloc[l,4]
                                         ,Count,'L2R'])
                    else:
                        List.append([Region_short[i][3],Region_short[i][0],W_start,W_end,'.',Count,'L2R'])

        #         loop for all mt within Lmt_id, Rmt_id
                for j in range(1,len(TempMt)):
        #             mt = TempMt.iloc[j]
                    mt_start = TempMt.iloc[j]['M_start']
                    mt_end = TempMt.iloc[j]['M_end']
                    Count = 0
                    for k in range(len(Bfile)):
        #                 ipdb.set_trace()
                        Fragmp = np.array(Bfile.loc[k,'Fragmp'])
                        if any((Fragmp >= mt_start)&(Fragmp<=mt_end)):
                            Count += 1
                    List_mt.append([Region_short[i][3],
                                    Lmt.iloc[0,0]+':'+str(Lmt.iloc[0,1])+'-'+str(Lmt.iloc[0,2])+
                                    ','+Lmt.iloc[0,3],
                                    Lmt.iloc[0,4],Lmt.iloc[0,7],Lmt.iloc[0,8],
                                    Lmt.iloc[0,9],Lmt.iloc[0,10],
                                    TempMt.iloc[j,0]+':'+str(TempMt.iloc[j,1])+'-'+str(TempMt.iloc[j,2])
                                     +','+TempMt.iloc[j,3],TempMt.iloc[j,4],TempMt.iloc[j,7],
                                    TempMt.iloc[j,8],TempMt.iloc[j,9],
                                    TempMt.iloc[j,10],Count,'L2R'])

            elif direction == 'Right' and len(Rmt)>0:
                Count = 0
                Table = np.zeros(Num_window)
                for k in range(len(Bfile)):
                    Fragmp = np.array(Bfile.loc[k,'Fragmp'])
                    Idx = np.minimum(np.maximum(np.array((R_end-Fragmp)/window_size),0),Num_window-1).astype(int)
                    Table[Idx] += 1
                for j in range(Num_window):
                    W_start = R_end - window_size*(j+1)+1
                    W_end = R_end - window_size*j
                    if j == Num_window:
                        W_start = R_start
                    TempMfile = Mt_inInterval(RegM2file,W_start,W_end,1)
                    Count = int(Table[j])

                    if len(TempMfile)>0:
                        for l in range(len(TempMfile)):
                            List.append([Region_short[i][3],Region_short[i][0],W_start,W_end,
                                         TempMfile.iloc[l,0]+':'+str(TempMfile.iloc[l,1])
                                         +'-'+str(TempMfile.iloc[l,2])+','+TempMfile.iloc[l,3]+','
                                         +TempMfile.iloc[l,4]
                                         ,Count,'R2L'])
                    else:
                        List.append([Region_short[i][3],Region_short[i][0],W_start,W_end,'.',Count,'R2L'])
        # loop for all mt within Lmt_id, Rmt_id
                for j in range(0,len(TempMt)-1):
                    mt = TempMt.iloc[j]
                    mt_start = mt['M_start']
                    mt_end = mt['M_end']
                    Count = 0
                    for k in range(len(Bfile)):
        #                 ipdb.set_trace()
                        Fragmp = np.array(Bfile.loc[k,'Fragmp'])
                        if any((Fragmp >= mt_start)&(Fragmp<=mt_end)):
                            Count += 1
                    List_mt.append([Region_short[i][3],
                                    TempMt.iloc[j,0]+':'+str(TempMt.iloc[j,1])+'-'+str(TempMt.iloc[j,2])
                                     +','+TempMt.iloc[j,3],TempMt.iloc[j,4],TempMt.iloc[j,7],
                                    TempMt.iloc[j,8],TempMt.iloc[j,9],
                                    TempMt.iloc[j,10],
                                    Rmt.iloc[0,0]+':'+str(Rmt.iloc[0,1])+'-'+str(Rmt.iloc[0,2])+
                                    ','+Rmt.iloc[0,3],
                                    Rmt.iloc[0,4],Rmt.iloc[0,7],Rmt.iloc[0,8],
                                    Rmt.iloc[0,9],Rmt.iloc[0,10],
                                    Count,'R2L'])

    Background_rate = pd.DataFrame(List,columns = ['LoopID','Chr','W_start','W_end','Motif overlap','GEM count', 'Direction'])
    saveBGratepath = 'Minji_data/Cohesin_results/01ALL/4kbext_dm/Slope/Background_rate.csv'
    Background_rate.to_csv(saveBGratepath,index=False)
    Motif_rate = pd.DataFrame(List_mt,columns = ['LoopID',
                                                 'Left motif','LM ID','LM CTCF ChIA_PET intensity','LM CTCF ChIA_Drop intensity',
                                                 'LM Cohesin ChIA_PET intensity','LM Cohesin ChIA_Drop intensity',
                                                 'Right motif','RM ID','RM CTCF ChIA_PET intensity','RM CTCF ChIA_Drop intensity',
                                                 'RM Cohesin ChIA_PET intensity','RM Cohesin ChIA_Drop intensity',
                                                 'GEM count', 'Direction'])
    saveMTratepath = 'Minji_data/Cohesin_results/01ALL/4kbext_dm/Slope/Motif_rate.csv'
    Motif_rate.to_csv(saveMTratepath,index=False)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pathR',type = str)
    parser.add_argument('--pathM',type = str)
    parser.add_argument('--pathM2',type = str)
    parser.add_argument('--pathB',type = str)
    parser.add_argument('--saveBGratepath',type = str)
    parser.add_argument('--saveMTratepath',type = str)
    
    parser.add_argument('--Start_pos',type = int)
    parser.add_argument('--End_pos',type = int)
    args = parser.parse_args()
    
    pathR = args.pathR
    pathM = args.pathM
    pathM2 = args.pathM2
    pathB = args.pathB
    saveBGratepath = args.saveBGratepath
    saveMTratepath = args.saveMTratepath
    Start_pos = args.Start_pos
    End_pos = args.End_pos
# pathR: path for region file (i.e. Minji_data/LHG0052H.e500.clusters.cis.bothanchint_G250.PETcnt_G9.motifannot)
# pathM: path for Motif file (motif_rate) (i.e. Minji_data/CTCF_motifs_STORM_hg38_Ext4kbBoth_with_supermotif_domain_id_v3.bed)
# pathM2: path for Motif file (background_rate) (i.e. Minji_data/CTCF_motifs_STORM_hg38_Ext4kbBoth_with_supermotif_domain_id_v3.bed)
# RegInterval: range(Start_pos,End_pos)
# pathB: path for directory of bedfiles extracted from data (i.e. Minji_data/Cohesin_results/01ALL/4kbext_dm/Bedfiles/)

# saveBGratepath: path for saving the background_rate table (i.e. Minji_data/Cohesin_results/01ALL/4kbext_dm/Slope/Background_rate.csv)
# saveMTratepath: path for saving the Motif_rate table (i.e. Minji_data/Cohesin_results/01ALL/4kbext_dm/Slope/Motif_rate.csv)
# Note that if you want to parallel this process, give different name for .csv files of each Thread.
# Then use MergeCSV.py to merge them (each directory should contain only 1 type of .csv file)    
    RegInterval = range(Start_pos,End_pos)
    
    mainfunc(pathR,pathM,pathM2,RegInterval,pathB,saveBGratepath,saveMTratepath)

