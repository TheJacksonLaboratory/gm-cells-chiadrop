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
            # If CTCF, c = #0000B2. If Cohesin, c = #005900.
            ax.plot(x_thick, y_thick, c = Colorcode, lw = 2)
    
    ax.tick_params(axis = 'y', which = 'both', left = False, top = False,labelleft = False)
    ax.set_ylabel('\n\n C#:'+str(num_frags),fontsize=6,rotation='horizontal', ha='right')
    ax.set_ylim(ymin=-1.25, ymax=len(y_idx)+0.75)
    return num_frags
        
def region_plot(region_boundary,bedfilepath,bedfiledir2, saveplotpath,title,Colorcode, silent = False):
    '''
    GEM file is large, avoid overloading it
    it contains a GEM and its fragments in each line.
    df_region_boundary also didn't change when switching region. 
    Avoid overloading.
    '''
#     ipdb.set_trace()
    extension = 'bed'
    all_filenames = [i for i in glob.glob(bedfilepath+'*.{}'.format(extension))]
    all_filenames.sort(reverse=True)
    NumFiles = len(all_filenames)
    myorder = [1,0,2]
    Maxnum_pair = int((NumFiles-3)/2)
    if Maxnum_pair>0:
        for pair_idx in range(Maxnum_pair):
            myorder.extend([3+pair_idx*2+1,3+pair_idx*2])

    all_filenames[:] = [all_filenames[i] for i in myorder]
    
    total_length = region_boundary[-1]-region_boundary[0]
    
    region_idx = 0
    Mid_Idx = 0
    df_region = [[] for dummy in range(len(all_filenames))]
    df_region_2side = [[] for dummy in range(Maxnum_pair)]
    Heights = np.zeros(len(all_filenames))
    Heights_2side = np.zeros(Maxnum_pair)
    for df_region_path in all_filenames:
        df_region[region_idx] = pd.read_csv(df_region_path, sep = '\t', 
                                            names = ['name', 'left_start', 'left_end', 'GEM_ID', 
                                                    'Fragment_number','???','mid_frags', 'right_start',
                                                     'right_end'])
        
        
        df_region[region_idx]['frag_coord_filtered_array'] = None
        for i in range(len(df_region[region_idx])):
#             ipdb.set_trace()
            df_region[region_idx].at[i,'frag_coord_filtered_array'] = get_list_of_frags(df_region[region_idx],i)
        
        Heights[region_idx] = max(len(np.arange(len(df_region[region_idx]['frag_coord_filtered_array']))[::-1])* 0.05+0.1,0.25)        
        
        if region_idx > 2 and np.mod(region_idx,2)==0:
            tempname = df_region_path.split('/')[-1].split('_')[0]+'_'+df_region_path.split('/')[-1].split('_')[1][:-1]
            df_region_2side_path = bedfiledir2+tempname+'.bed'
            df_region_2side[Mid_Idx] = pd.read_csv(df_region_2side_path, sep = '\t', 
                                                    names = ['name', 'left_start', 'left_end', 'GEM_ID', 
                                                            'Fragment_number','???','mid_frags', 'right_start',
                                                             'right_end'])
            df_region_2side[Mid_Idx]['frag_coord_filtered_array'] = None
            for i in range(len(df_region_2side[Mid_Idx])):
#             ipdb.set_trace()
                df_region_2side[Mid_Idx].at[i,'frag_coord_filtered_array'] = get_list_of_frags(df_region_2side[Mid_Idx],i)
            
            Heights_2side[Mid_Idx] = max(len(np.arange(len(df_region_2side[Mid_Idx]['frag_coord_filtered_array']))[::-1])* 0.05+0.1,0.25)
            Mid_Idx += 1
            
        region_idx += 1
    
# #     This is for single plot
#     fig, ax = plt.subplots(NumFiles, 1, sharex= True,
#                            figsize = (NumFiles, np.sum(Heights)),
#                            gridspec_kw={'height_ratios': Heights})
    
#     Num_frags = np.zeros(len(all_filenames))
#     for region_idx in range(len(all_filenames)):
#         Num_frags[region_idx] = plot_each_block(df_region[region_idx]['frag_coord_filtered_array'],fig,ax[region_idx],total_length,sorted_by_len = False)

    
    
#     ax[0].set_xticks([region_boundary[0],region_boundary[-1]])
#     ax[0].set_xticklabels([region_boundary[0],region_boundary[-1]],fontsize=4)
    
#     title = title+'_anchcomp_'+str(int(np.sum(Num_frags[0:3])))+'_M-'+str(int((len(Num_frags)-3)/2))+'_'+str(int(np.sum(Num_frags[3:])))
#     ax[0].set_title(title, fontsize = 10)
#     fig.subplots_adjust(hspace=0)
# #     plt.savefig('Minji_data/Cohesin_results/01ALL/4kbext_dm/Plots/{}.png'.format(title), 
# #                 dpi = 600, bbox_inches='tight')
#     plt.savefig(saveplotpath+'{}.png'.format(title),dpi = 600, bbox_inches='tight')
    
#     if not silent:
#         plt.show()
#     plt.close()

# # # # This is for the split plots
    Num_frags = np.zeros(len(all_filenames))
    Num_frags_2side = np.zeros(Maxnum_pair)
    
    for figlistid in range(Maxnum_pair+1):
        if figlistid == 0:
            fig, ax = plt.subplots(3, 1, sharex= True,
                                   figsize = (3, np.sum(Heights[0:3])),
                                   gridspec_kw={'height_ratios': Heights[0:3]})
            for region_idx in range(3):
                Num_frags[region_idx] = plot_each_block(df_region[region_idx]['frag_coord_filtered_array'],fig,ax[region_idx],total_length,Colorcode,sorted_by_len = False)
            
            MetaX = ax[0].get_xlim()
            ax[0].set_xticks([region_boundary[0],region_boundary[-1]])
            ax[0].set_xticklabels([region_boundary[0],region_boundary[-1]],fontsize=4)
            Title = title+'_anchcomp_'+str(int(np.sum(Num_frags[0:3])))
            ax[0].set_title(Title, fontsize = 10)
            fig.subplots_adjust(hspace=0)
            try:
                plt.savefig(saveplotpath+'{}.pdf'.format(Title), bbox_inches='tight')
                if not silent:
                    plt.show()
                plt.close()
            except:
                plt.savefig(saveplotpath+'{}.pdf'.format(Title),dpi = 100, bbox_inches='tight')
                if not silent:
                    plt.show()
                plt.close()
        else:
            
            fig, ax = plt.subplots(3, 1, sharex= True,
                                   figsize = (3, np.sum(Heights[3+2*(figlistid-1):3+2*(figlistid)])+Heights_2side[figlistid-1]),
                                   gridspec_kw={'height_ratios': np.concatenate((Heights[3+2*(figlistid-1):3+2*(figlistid)],[Heights_2side[figlistid-1]]))})
            for region_idx in range(3+2*(figlistid-1),3+2*(figlistid)):
                Num_frags[region_idx] = plot_each_block(df_region[region_idx]['frag_coord_filtered_array'],fig,ax[np.mod(region_idx-3,2)],total_length,Colorcode,sorted_by_len = False)
            
            Num_frags_2side[figlistid-1] = plot_each_block(df_region_2side[figlistid-1]['frag_coord_filtered_array'],fig,ax[2],total_length,Colorcode,sorted_by_len = False)
#             ipdb.set_trace()
            MidName = all_filenames[3+2*(figlistid-1)].split('/')[-1].split('_')[1][2:-1]
            ax[0].set_xlim(xmin=MetaX[0], xmax=MetaX[-1])
            ax[0].set_xticks([region_boundary[0],region_boundary[-1]])
            ax[0].set_xticklabels([region_boundary[0],region_boundary[-1]],fontsize=4)
            Title = title+'_M-'+str(MidName)+'_1side:'+str(int(np.sum(Num_frags[3+2*(figlistid-1):3+2*(figlistid)])))+'_2side:'+str(int(Num_frags_2side[figlistid-1]))
            ax[0].set_title(Title, fontsize = 10)
            fig.subplots_adjust(hspace=0)
            try:
                plt.savefig(saveplotpath+'{}.pdf'.format(Title), bbox_inches='tight')
                if not silent:
                    plt.show()
                plt.close()
            except:
                print(Title,' : This plot use low dpi since too long.')
                plt.savefig(saveplotpath+'{}.pdf'.format(Title),dpi = 100, bbox_inches='tight')
                if not silent:
                    plt.show()
                plt.close()
    return None

def region_plot_2sides(region_boundary,bedfilepath, saveplotpath,title,Colorcode, silent = False):
    '''
    GEM file is large, avoid overloading it
    it contains a GEM and its fragments in each line.
    df_region_boundary also didn't change when switching region. 
    Avoid overloading.
    '''
#     ipdb.set_trace()
    extension = 'bed'
    all_filenames = [i for i in glob.glob(bedfilepath+'*.{}'.format(extension))]
    all_filenames.sort(reverse=True)
    NumFiles = len(all_filenames)
    
    total_length = region_boundary[-1]-region_boundary[0]
    
    region_idx = 0
    df_region = [[] for dummy in range(len(all_filenames))]
    Heights = np.zeros(NumFiles)
    for df_region_path in all_filenames:
        df_region[region_idx] = pd.read_csv(df_region_path, sep = '\t', 
                                            names = ['name', 'left_start', 'left_end', 'GEM_ID', 
                                                    'Fragment_number','???','mid_frags', 'right_start',
                                                     'right_end'])
        
        df_region[region_idx]['frag_coord_filtered_array'] = None
        for i in range(len(df_region[region_idx])):
#             ipdb.set_trace()
            df_region[region_idx].at[i,'frag_coord_filtered_array'] = get_list_of_frags(df_region[region_idx],i)
        
        Heights[region_idx] = max(len(np.arange(len(df_region[region_idx]['frag_coord_filtered_array']))[::-1])* 0.05+0.1,0.25)        
        
        fig, ax = plt.subplots(1, 1,figsize = (1, Heights[region_idx]))
        Num_frags = plot_each_block(df_region[region_idx]['frag_coord_filtered_array'],fig,ax,total_length,Colorcode,sorted_by_len = False)
        ax.set_xlim(xmin=region_boundary[0], xmax=region_boundary[-1])
        ax.set_xticks([region_boundary[0],region_boundary[-1]])
        ax.set_xticklabels([region_boundary[0],region_boundary[-1]],fontsize=4)
        
        MidName = all_filenames[region_idx].split('/')[-1].split('_')[-1][2:-4]
        Title = title+'_M-'+str(MidName)+'_'+str(int(Num_frags))
        ax.set_title(Title, fontsize = 10)
        plt.savefig(saveplotpath+'{}.pdf'.format(Title), bbox_inches='tight')
        
        region_idx += 1
        try:
            plt.savefig(saveplotpath+'{}.pdf'.format(Title), bbox_inches='tight')
            if not silent:
                plt.show()
            plt.close()
        except:
            plt.savefig(saveplotpath+'{}.pdf'.format(Title), bbox_inches='tight')
            if not silent:
                plt.show()
            plt.close()
# #     This is for single plot
#     fig, ax = plt.subplots(NumFiles, 1, sharex= True,
#                            figsize = (NumFiles, np.sum(Heights)),
#                            gridspec_kw={'height_ratios': Heights})
    
#     Num_frags = np.zeros(len(all_filenames))
#     for region_idx in range(len(all_filenames)):
#         Num_frags[region_idx] = plot_each_block(df_region[region_idx]['frag_coord_filtered_array'],fig,ax[region_idx],total_length,sorted_by_len = False)

    
    
#     ax[0].set_xticks([region_boundary[0],region_boundary[-1]])
#     ax[0].set_xticklabels([region_boundary[0],region_boundary[-1]],fontsize=4)
    
#     title = title+'_anchcomp_'+str(int(np.sum(Num_frags[0:3])))+'_M-'+str(int((len(Num_frags)-3)/2))+'_'+str(int(np.sum(Num_frags[3:])))
#     ax[0].set_title(title, fontsize = 10)
#     fig.subplots_adjust(hspace=0)
# #     plt.savefig('Minji_data/Cohesin_results/01ALL/4kbext_dm/Plots/{}.png'.format(title), 
# #                 dpi = 600, bbox_inches='tight')
#     plt.savefig(saveplotpath+'{}.png'.format(title),dpi = 600, bbox_inches='tight')
    
#     if not silent:
#         plt.show()
#     plt.close()

        
    return None

def fn(regionfilepath):
    regionfile = pd.read_csv(regionfilepath,sep = '\t',header=None)
    region_boundary = np.array(regionfile.iloc[0][1:3]).astype(int)
    cr_id = regionfile.iloc[0][3]
    cr_id = cr_id[:-3]
    title = cr_id
    bedfilepath = bedfiledir+cr_id+'_'
    region_plot(region_boundary,bedfilepath,bedfiledir2, saveplotpath,title, silent = True)
#     region_plot_2sides(region_boundary,bedfilepath, saveplotpath,title, silent = True)
    
    return None

def mainplotfunc(bedfiledir,bedfiledir2,saveplotpath,regionfiledir,Multithread=20,MaxCr=3015):
#     bedfiledir: path for 1 side sorted GEMs files. i.e. './bedfiles/'
#     bedfiledir: path for 2 side sorted GEMs files. i.e. './2sides_bedfiles/'
#     saveplotpath: path to save plots. i.e. './Plots/'
#     regionfiledir: path of processed region_cr files.
#     Multithread: number of multithreading
#     MaxCr: maximum number of cr_id

    regionfilelist = [i for i in glob.glob(regionfiledir+'Region_cr*')]
    with Pool(Multithread) as p:
        r = list(tqdm(p.imap(fn, regionfilelist), total = MaxCr-1))
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bedfiledir',type = str)
    parser.add_argument('--bedfiledir2',type = str)
    parser.add_argument('--saveplotpath',type = str)
    parser.add_argument('--regionfiledir',type = str)
    
    parser.add_argument('--Multithread',type = int)
    parser.add_argument('--MaxCr',type = int)
    args = parser.parse_args()
    
    bedfiledir = args.bedfiledir
    bedfiledir2 = args.bedfiledir2
    saveplotpath = args.saveplotpath
    regionfiledir = args.regionfiledir
    Multithread = args.Multithread
    MaxCr = args.MaxCr
    
    mainplotfunc(bedfiledir,bedfiledir2,saveplotpath,regionfiledir,Multithread,MaxCr)