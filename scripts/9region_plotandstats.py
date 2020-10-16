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
            c = Colorcode
            # If CTCF, c = #0000B2. If Cohesin, c = #005900. If RNAPII, c = #590059
#             if Name1 == 'CTCF':
#                 c = '#0000B2'
#             elif Name1 == 'RNAPII':
#                 c = '#590059'
#             elif Name1 == 'Cohesin':
#                 c = '#005900'
                
            ax.plot(x_thick, y_thick, c = c, lw = 2)
    
    ax.tick_params(axis = 'y', which = 'both', left = False, top = False,labelleft = False)
    ax.set_ylabel('\n\n C#:'+str(num_frags),fontsize=6,rotation='horizontal', ha='right')
    ax.set_ylim(ymin=-1.25, ymax=len(y_idx)+0.75)
    return num_frags
        
def region_plot(region_boundary,bedfilepath,RegionID, saveplotpath,TTname,FFname,Colorcode, silent = True):
    '''
    GEM file is large, avoid overloading it
    it contains a GEM and its fragments in each line.
    df_region_boundary also didn't change when switching region. 
    Avoid overloading.
    '''
#     ipdb.set_trace()
    Typelist = ['Left','Right','Both','M2L','MaL','M2R','MaR']
    total_length = region_boundary[-1]-region_boundary[0]
    df_region = [[] for dummy in range(7)]
    Heights = np.zeros([7])
    
    for kk in range(7):
        Type = Typelist[kk]
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
    
    
    fig, ax = plt.subplots(7, 1, sharex= True,
                           figsize = (7, np.sum(Heights)),
                           gridspec_kw={'height_ratios': Heights})
    
    
    ax[0].get_yaxis().set_ticks([])
#     Add loop info here
    ax[0].set_xlim(xmin=region_boundary[0], xmax=region_boundary[-1])
    ax[0].set_xticks([region_boundary[0],region_boundary[-1]])
    ax[0].set_xticklabels([region_boundary[0],region_boundary[-1]],fontsize=4)
    
    Num_frags = np.zeros([7])
    for kk in range(7):
        Num_frags[kk] = plot_each_block(df_region[kk]['frag_coord_filtered_array'],fig,ax[kk],total_length,Colorcode,sorted_by_len = False)
    
#     title = 'Region'+RegionID
    title = TTname
    ax[0].set_title(title, fontsize = 10)
    fig.subplots_adjust(hspace=0)
    
    tmplist=[]
    for num in Num_frags:
        tmplist.append(str(int(num)))
    
    plt.savefig(saveplotpath+FFname+'.pdf', bbox_inches='tight')
    
    if not silent:
        plt.show()
    plt.close()
         
    return Num_frags






def fn(bedfiledir,saveplotpath,region_boundary,RegionID,TTname,FFname,Colorcode):
#     ipdb.set_trace()
#     Lmostsite = Loop.iloc[1]-Loop.iloc[4]+Loop.iloc[2]
#     Rmostsite = Loop.iloc[5]+Loop.iloc[4]-Loop.iloc[2]
#     region_boundary = np.array([Lmostsite,Rmostsite]).astype(int)
    bedfilepath = bedfiledir+'Region'+RegionID
    NumFrags = region_plot(region_boundary,bedfilepath,RegionID, saveplotpath,TTname,FFname,Colorcode, silent = True)
    return NumFrags

def mainfunc(bedfiledir,saveplotpath,regionfilepath,savecsvpath,Fname,Tname,Colorcode):
# bedfiledir: savebedpath in 9region_sort.py
# saveplotpath: path to save plots
# regionfilepath: 9 region file
# savecsvpath: path to save tables
# Fname: save plots/table name (i.e. CTCF + one of the 9 region name)
# Tname: title of plots (i.e. CTCF + one of the 9 region name)
# Colorcode: colorcode for ploting fragments


    header_list = ['chr', 'S','E','Name']
    Loops = pd.read_csv(regionfilepath,sep = '\t',names=header_list)
    Record = []
    for i in tqdm(range(0,len(Loops),3)):
    #     ipdb.set_trace()
        Loop = Loops.iloc[i:(i+3)]
        RegionID = Loop.iloc[0]['Name'][5:]

        region_boundary = np.array([Loop.iloc[0]['S'],Loop.iloc[1]['E']]).astype(int)
    #     ipdb.set_trace()
        TTname = Tname+'\n'+Loop.iloc[0]['chr']+':'+str(Loop.iloc[0]['S'])+'-'+str(Loop.iloc[1]['E'])
#         if Name1 == 'CTCF':
#             Fname = 'GM12878-CTCF-pooled_comp_FDR_0.2_PASS.RNAPII-peaks-'+Name2
#         elif Name1 == 'RNAPII':
#             Fname = 'GM12878-RNAPII-pooledv2_comp_FDR_0.2_PASS.RNAPII-peaks-'+Name2
#         elif Name1 == 'Cohesin':
#             Fname = 'GM12878-cohesin-pooled_comp_FDR_0.2_PASS.RNAPII-peaks-'+Name2
        
        LEN = region_boundary[-1]-region_boundary[0]
        FFname = Fname +'.'+Loop.iloc[0]['chr']+':'+str(Loop.iloc[0]['S'])+'-'+str(Loop.iloc[1]['E'])+'('+str(LEN)+'bp)'
        
        
        NumFrags = fn(bedfiledir,saveplotpath,region_boundary,RegionID,TTname,FFname,Colorcode)
        Temp = []
        Temp.append(RegionID)
        Temp.append(Loop.iloc[0]['chr']+':'+str(Loop.iloc[0]['S'])+'-'+str(Loop.iloc[1]['E']))
        for num in NumFrags:
            Temp.append(int(num))

        Record.append(Temp)

    DF = pd.DataFrame(Record,columns = ['Region ID','Region coordinate',
                               '# of complexes for Left','# of complexes for Right',
                               '# of complexes for Both','# of complexes for M2L',
                               '# of complexes for MaL','# of complexes for M2R',
                                '# of complexes for MaR'])
    DF.to_csv(savecsvpath+Fname+'.stats.csv',index=False)
    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bedfiledir',type = str)
    parser.add_argument('--saveplotpath',type = str)
    parser.add_argument('--regionfilepath',type = str)
    parser.add_argument('--savecsvpath',type = str)
    parser.add_argument('--Fname',type = str)
    parser.add_argument('--Tname',type = str)
    parser.add_argument('--Colorcode',type = str)
    
    args = parser.parse_args()
    
    bedfiledir = args.bedfiledir
    saveplotpath = args.saveplotpath
    regionfilepath = args.regionfilepath
    savecsvpath = args.savecsvpath
    Fname = args.Fname
    Tname = args.Tname
    Colorcode = args.Colorcode

    mainfunc(bedfiledir,saveplotpath,regionfilepath,savecsvpath,Fname,Tname,Colorcode)