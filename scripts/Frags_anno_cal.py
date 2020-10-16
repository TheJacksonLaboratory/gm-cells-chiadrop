#!/usr/bin/env python
# coding: utf-8

# In[13]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
from pybedtools import BedTool
import ipdb
import argparse
import time
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

def Fg_inInterval(Fg_start,s2,Fg_end,e2):
    return (np.minimum(Fg_end,e2)-np.maximum(Fg_start,s2))>0

def get_list_of_frags(df_region,i):

    frag_left_start = df_region.loc[i]['Lm_frag_start']
    frag_right_end = df_region.loc[i]['Rm_frag_end']
    
    Mcount = df_region.loc[i]['FragNum']
    Num_mid_frags = Mcount - 2
    
    Lst_of_frags = [[int(df_region.loc[i]['Lm_frag_start']),int(df_region.loc[i]['Lm_frag_end'])]]
    if Num_mid_frags > 0:
        list_of_frag_coord = list(map(int, df_region.loc[i]['Mid_frags'].split(',')))
        for k in range(Num_mid_frags):
            Lst_of_frags.append([int(list_of_frag_coord[2*k]),int(list_of_frag_coord[2*k+1])])
    Lst_of_frags.append([int(df_region.loc[i]['Rm_frag_start']),int(df_region.loc[i]['Rm_frag_end'])])
    return list(Lst_of_frags)

# print(get_list_of_fragsmp(Bfile,0))
# # Bfile.insert(9,'Fragmp',None)
# Bfile.at[0,'Fragmp'] = get_list_of_fragsmp(Bfile,0)
# Bfile.head()


# In[14]:


def mainfunc(pathR,pathM,pathM2,RegInterval,pathB,saveFragannopath):
    # pathR: path for region file (i.e. Minji_data/LHG0052H.e500.clusters.cis.bothanchint_G250.PETcnt_G9.motifannot)
    # pathM: path for Motif file (i.e. Minji_data/CTCF_motifs_STORM_hg38_Ext4kbBoth_with_supermotif_domain_id_v3.bed)
    # RegInterval: range(Start_pos,End_pos)
    # pathB: path for directory of bedfiles extracted from data (i.e. Minji_data/Cohesin_results/01ALL/4kbext_dm/Bedfiles/)

    # saveFragannopath: path for saving the Frags_anno table (i.e. Minji_data/Cohesin_results/01ALL/4kbext_dm/Tables/Frags_anno.csv)
    # Note that if you want to parallel this process, give different name for .csv files of each Thread.
    # Then use MergeCSV.py to merge them (each directory should contain only 1 type of .csv file)
    
    
    # Region_short = BedTool('Region_short_dm.bed')

    # pathR = 'Minji_data/LHG0052H.e500.clusters.cis.bothanchint_G250.PETcnt_G9.motifannot'
    Region = BedTool(pathR)
    # Remove unnecessary entries
    Region_short = Region.groupby(g=[1,2,6,12,14,20,8,9,16,21], c=[12], o=['count'])

    # path = 'Minji_data/CTCF_motifs_STORM_hg38_Ext4kbBoth_chiapet_chiadrop_intensity_id.bed'
    # Mfile = pd.read_csv(path, sep = '\t',names = ['chr','M_start','M_end','Sign','Pet_int','Drop_int','M_name'])

#     path = 'Minji_data/CTCF_motifs_STORM_hg38_Ext4kbBoth_with_supermotif_domain_id_v7.sorted.bed'
    Mfile = pd.read_csv(pathM, sep = '\t',names = ['chr','M_start','M_end','Sign','M_name',
                                                  'dmID','Side',
                                                  'CTCF_Pet_int','CTCF_Drop_int','Coh_Pet_int','Coh_Drop_int'])

#     path = 'Minji_data/CTCF_motifs_STORM_hg38_Ext4kbBoth.sorted.id.bed'
    M2file = pd.read_csv(pathM2, sep = '\t',names = ['chr','M_start','M_end','Sign','M_name'])

#     Thread = '0'
    Max_iter = len(Region_short)
#     Chuck = int(Max_iter/20)
    # LPlist = [8476,3647,4936]
    # List = []
    # List_mt = []
    List_new = []
    window_size = 8000

    # for i in tqdm(range(Chuck*int(Thread),Chuck*(int(Thread)+1))):
    for i in RegInterval:
        # print(Region_short[i])
        total_length = int(Region_short[i][2])-int(Region_short[i][1])
        # print(total_length)
        # Count how many slicing should be
        Num_window = int(np.ceil(total_length/window_size))
        lpID = Region_short[i][3]
    #     print(lpID)

        R_start = int(Region_short[i][1])
        R_end = int(Region_short[i][2])
    # ===========================================================================================================    
        for direction in ['Left','Right','Both','None']:
#             path = 'Minji_data/Final_data_results/CTCF_NR_results/01PASS_dm/Bedfiles/{}_{}.bed'.format(lpID,direction)
            pathB_all = pathB+'{}_{}.bed'.format(lpID,direction)
            try:
                Bfile = pd.read_csv(pathB_all, sep = '\t',names = ['chr','Lm_frag_start','Lm_frag_end','GEMID','FragNum',
                                                              '???','Mid_frags','Rm_frag_start','Rm_frag_end'])
            except:
                continue
            # Bfile preprocessing
            Bfile.insert(9,'Frags',None)
            Bfile.insert(10,'Len',Bfile['Rm_frag_end']-Bfile['Lm_frag_start'])
    #         ipdb.set_trace()
            Bfile = Bfile.sort_values(by = ['Len']).reset_index(drop=True)
            for k in range(len(Bfile)):
                Bfile.at[k,'Frags'] = get_list_of_frags(Bfile,k)
        #     ipdb.set_trace()
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
                    TempMt = Mt_inInterval(RegMfile,int(Lmt['M_start']),int(Rmt['M_end']),1).sort_values(by = ['M_start'])
                else:
                    TempMt = Mt_inInterval(RegMfile,int(Lmt['M_start']),Rmp_end,1).sort_values(by = ['M_start'])
            else:
                if len(Rmt)>0:
                    TempMt = Mt_inInterval(RegMfile,Lmp_start,int(Rmt['M_end']),1).sort_values(by = ['M_start'])
                else:
                    TempMt = Mt_inInterval(RegMfile,Lmp_start,Rmp_end,1).sort_values(by = ['M_start'])



            for k in range(len(Bfile)):
                Frags = np.array(Bfile.loc[k,'Frags'])
    #             ipdb.set_trace()
                sublist = [str() for c in 'c' * len(Frags)]

                for j in range(len(TempMt)):
                    mt_start = TempMt.iloc[j]['M_start']
                    mt_end = TempMt.iloc[j]['M_end']
    #                 ipdb.set_trace()
                    Condi = Fg_inInterval(Frags[:,0],mt_start,Frags[:,1],mt_end)
                    Res = [KK for KK, val in enumerate(Condi) if val]
                    for _,KK in enumerate(Res):
    #                     ipdb.set_trace()
                        sublist[KK] += TempMt.iloc[j]['M_name']


                SubStr = ''
                for j in range(len(sublist)-1):
                    if sublist[j]:
                        SubStr += sublist[j]+','
                    else:
                        SubStr += '0'+','
                if sublist[len(sublist)-1]:
                    SubStr += sublist[len(sublist)-1]
                else:
                    SubStr += '0'     

                List_new.append([Region_short[i][3],direction,SubStr,Bfile.loc[k,'GEMID']])
    Frags_anno = pd.DataFrame(List_new,columns = ['Domain ID','Category','Fragment_annotation','GEM_ID'])
#     saveFragannopath = 'Minji_data/Final_data_results/CTCF_NR_results/01PASS_dm/Tables/Frags_anno.csv'
    Frags_anno.to_csv(saveFragannopath,index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pathR',type = str)
    parser.add_argument('--pathM',type = str)
    parser.add_argument('--pathM2',type = str)
    parser.add_argument('--pathB',type = str)
    parser.add_argument('--saveFragannopath',type = str)
    
    parser.add_argument('--Start_pos',type = int)
    parser.add_argument('--End_pos',type = int)
    args = parser.parse_args()
    
    pathR = args.pathR
    pathM = args.pathM
    pathM2 = args.pathM2
    pathB = args.pathB
    saveFragannopath = args.saveFragannopath
    Start_pos = args.Start_pos
    End_pos = args.End_pos
# pathR: path for region file (i.e. Minji_data/LHG0052H.e500.clusters.cis.bothanchint_G250.PETcnt_G9.motifannot)
# pathM: path for Motif file (motif_rate) (i.e. Minji_data/CTCF_motifs_STORM_hg38_Ext4kbBoth_with_supermotif_domain_id_v3.bed)
# pathM2: path for Motif file (background_rate) (i.e. Minji_data/CTCF_motifs_STORM_hg38_Ext4kbBoth_with_supermotif_domain_id_v3.bed)
# RegInterval: range(Start_pos,End_pos)
# pathB: path for directory of bedfiles extracted from data (i.e. Minji_data/Cohesin_results/01ALL/4kbext_dm/Bedfiles/)

# saveFragannopath: path for saving the Frags_anno table (i.e. Minji_data/Cohesin_results/01ALL/4kbext_dm/Tables/Frags_anno.csv)
# Note that if you want to parallel this process, give different name for .csv files of each Thread.
# Then use MergeCSV.py to merge them (each directory should contain only 1 type of .csv file)    
    RegInterval = range(Start_pos,End_pos)
    
    mainfunc(pathR,pathM,pathM2,RegInterval,pathB,saveFragannopath)


# In[3]:


Test = '12'
print(int(Test)==12)

