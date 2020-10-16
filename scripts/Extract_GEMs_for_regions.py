#!/usr/bin/env python
# coding: utf-8

# In[58]:


import pybedtools
from pybedtools import BedTool
import numpy as np
# from tqdm import tqdm
import pandas as pd
# import ipdb
import argparse

# These two lines let the pybedtool know where to put the temp files.
# cleanup() will remove all the temp file create from this session in temp file folder 
# Thread = '0'

# Process the given region, retrive all valid GEMs
def ProcessRegion(RawGEMs):
    Temp = RawGEMs.groupby(g=[1,2,3,5], c=[5,6,7,8], o=['count','collapse','collapse','collapse'])
    RefineGEMs = Temp.filter(lambda F: int(F[4]) > 1)
    # Need this to keep the result of filter! All these files can be manually removed after the job 
    Test = BedTool(RefineGEMs).saveas()
    Tempstr = ''
    for i in range(Test.count()):
        Start = np.fromstring(Test[i][6], dtype = np.int,  sep =', ' )
        End = np.fromstring(Test[i][7], dtype = np.int,  sep =', ' )
        Mcount = Test[i][5].count('P')
        Start.sort()
        End.sort()
        # chrom, start_min, end_min, GEM ID, #(fragments)  start_max, end_max,
        for j in range(len(Start)):
            if j == 0:
                Tempstr += Test[i][0]+' '+ str(Start[j]) + ' ' + str(End[j])+' '+ Test[i][3] + ' ' + str(len(Start)) + ' ' + str(Mcount)+' '
            elif len(Start)!=2 and j != (len(Start)-1):
                Tempstr += str(Start[j])+','+str(End[j])+','
            elif j == (len(Start)-1):
                Tempstr += str('-1,-1')+' ' + str(Start[j]) + ' ' + str(End[j]) + '\n'
                
    FinalGEMs = BedTool(Tempstr,from_string=True)
    return FinalGEMs

# Check for the left/right/both/none condition
def Checkintersect(s1,s2,e1,e2):
    return (min(e1, e2) - max(s1, s2)) > 0

def inInterval(FF,Temp,Type,Length,CHR):
#     ipdb.set_trace()
    if CHR == FF[0]:
#         print(CHR,FF[0],'True')
        interval = [0,0,0,0]
        interval[0] = Temp[0]-Length
        interval[1] = Temp[0]+Length+19
        interval[2] = Temp[1]-Length-19
        interval[3] = Temp[1]+Length

        NumFrag = FF[4]
        Start = list(map(int, FF[1:3]))
        End = list(map(int, FF[-2:]))
        if Type == 'Left':
            return (Checkintersect(interval[0],Start[0],interval[1],Start[1])) and not (Checkintersect(interval[2],End[0],interval[3],End[1]))
        elif Type == 'Right':
            return not (Checkintersect(interval[0],Start[0],interval[1],Start[1])) and (Checkintersect(interval[2],End[0],interval[3],End[1]))
        elif Type == 'Both':
            return (Checkintersect(interval[0],Start[0],interval[1],Start[1])) and (Checkintersect(interval[2],End[0],interval[3],End[1]))
        else:
            return not (Checkintersect(interval[0],Start[0],interval[1],Start[1])) and not (Checkintersect(interval[2],End[0],interval[3],End[1]))
    else:
        return False
# Classify FinalGEMs based on Type (left/right/both) in Region. 
# left: start_min to start_min+Length; right: end_max-Length to end_max 
# Also save the result automatically. Can change path.
def SortGEM(FinalGEMs, Region,Type,Length,savebedpath):
    Temp = list(map(int, Region[4:6]))
    CHR = Region[0]
    TypeGEMs = FinalGEMs.filter(inInterval, Temp,Type,Length,CHR).sort()#.saveas()
    # I use loop id to indicate region
#     savebedpath = 'Minji_data/Cohesin_results/01ALL/4kbext_dm/Bedfiles/'
    TypeGEMs = TypeGEMs.moveto(savebedpath+str(Region[3])+'_'+Type+'.bed')
    if Type == 'Both':
        Flag = 2
    elif Type == 'None':
        Flag = 0
    else:
        Flag = 1
    
    Count = 0
    Tot = TypeGEMs.count()
#     Check if any fragments intersect with middle motif
    for i in range(Tot):
        Mcount = int(TypeGEMs[i][5])
        if Mcount> Flag:
            Count = Count+1
    return Tot-Count, Count

def mainfunc(path1,path2,savebedpath,savecsvpath,tmpfilepath,RegInterval,Thread,Length = 4000):
# path1: path for GEMs (i.e. ___ALL.region.PEanno)
# path2: path for Region (i.e. ____PETcnt_G9.motifannot)
# savebedpath: path for saving extracted GEMs in .bed
# savecsvpath: path for saving summary table in .csv
# tmpfilepath: path for saving tmpfiles produced by pybedtool, a directory
# Thread: for naming the csv file. (i.e. '0')
# Length: Length of extension. Default = 4000 (int)
    pybedtools.helpers.cleanup()
    pybedtools.set_tempdir(tmpfilepath)
    # Specify for the path of ___ALL.region.PEanno and import it (GEMs)
#     path1 = 'Minji_data/SHG0180-181-182NR_hg38_cohesin_FDR_0.1_ALL_motifext4kbboth.region.PEanno'
    ChIA_Drop = BedTool(path1)
    
    # Specify for the path of ____PETcnt_G9.motifannot and import it (anchors, regions)
#     path2 = 'Minji_data/LHG0052H.e500.clusters.cis.bothanchint_G250.PETcnt_G9.motifannot.sorted.domains'
    Region = BedTool(path2)

    # Remove unnecessary entries
    Region_short = Region.groupby(g=[1,2,6,12,14,20,8,9,16,21], c=[12], o=['count'])
#     Region_short.moveto('Region_short.bed')
#     Region_short = BedTool('Region_short.bed')
    Max_iter = Region_short.count()
#     Length = 4000

    Dict = {}
# Dict = {'Type/loopID': ['Left_0','Left_1','Right_0','Right_1','Both_0','Both_1','None_0','None_1','Total','Left Intensity', 'Right Intensity','Left motif strand', 'Right motif strand']}
    for i in RegInterval:
        # NowRegion: chrom, start_min, end_max, loop id, ...
        # This line can be improved...
    #     NowRegion = NowRegion.saveas('NowRegion.bed')
        NowRegion = BedTool(Region_short[i:i+1]).saveas()
        # Find all fragments that intersect with Nowregion
        Intersection = ChIA_Drop.intersect(NowRegion,wa=True)
        # Append original start/and. Technical purpose for using groupby...
        results = [(f[0],'0','0',f[3],f[4],f[5],f[1],f[2]) for f in Intersection]
        Intersection = BedTool(results)

        # Sort the grouping key!!!! Otherwise the later groupby doesn't work as intended...
        Intersection = Intersection.sort(chrThenScoreA = True)
        # Extract the valid GEMs
        FinalGEMs = ProcessRegion(Intersection)
        # Classify+sort+save 
        Count_L0,Count_L1 = SortGEM(FinalGEMs, NowRegion[0],'Left',Length,savebedpath)
        Count_R0,Count_R1 = SortGEM(FinalGEMs, NowRegion[0],'Right',Length,savebedpath)
        Count_B0,Count_B1 = SortGEM(FinalGEMs, NowRegion[0],'Both',Length,savebedpath)
        Count_N0,Count_N1 = SortGEM(FinalGEMs, NowRegion[0],'None',Length,savebedpath)
        Total = Count_L0+Count_L1+Count_R0+Count_R1+Count_B0+Count_B1+Count_N0+Count_N1

        # Write into dictionary
        Dict[NowRegion[0][3]] = [NowRegion[0][3],Count_L0,Count_L1,Count_L0+Count_L1,(Count_L0+Count_L1)/Total*100,
                                 Count_R0,Count_R1,Count_R0+Count_R1,(Count_R0+Count_R1)/Total*100,
                                 Count_B0,Count_B1,Count_B0+Count_B1,(Count_B0+Count_B1)/Total*100,
                                 Count_N0,Count_N1,Count_N0+Count_N1,(Count_N0+Count_N1)/Total*100,
                                 Total,Total-(Count_N0+Count_N1),(Total-(Count_N0+Count_N1))/Total*100,
                                 NowRegion[0][6],NowRegion[0][7],NowRegion[0][8],NowRegion[0][9],
                                 NowRegion[0][0]+':'+str(NowRegion[0][1])+'-'+str(NowRegion[0][2])]
        # Clear all temp files for this session
        pybedtools.helpers.cleanup()

    RenameCol = {}
    NewCol = ['LoopID','Left_0','Left_1','Left_Tol','Left_Tol %','Right_0','Right_1','Right_Tol','Right_Tol %',
              'Both_0','Both_1','Both_Tol','Both_Tol %',
              'None_0','None_1','None_Tol','None_Tol %','Total','Total-None','Total-None %',
              'Left Intensity', 'Right Intensity','Left motif strand', 'Right motif strand',
              'Region']
    for i, name in enumerate(NewCol):
        RenameCol[i] = NewCol[i]

    DF = pd.DataFrame.from_dict(Dict,orient = 'index').rename(columns = RenameCol)
    # savecsvpath = 'Minji_data/Cohesin_results/01ALL/4kbext_dm/'
    DF.to_csv(savecsvpath+'LRBNstats_'+Thread+'.csv',index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--path1',type = str)
    parser.add_argument('--path2',type = str)
    parser.add_argument('--savebedpath',type = str)
    parser.add_argument('--savecsvpath',type = str)
    parser.add_argument('--tmpfilepath',type = str)
    
    parser.add_argument('--Start_pos',type = int)
    parser.add_argument('--End_pos',type = int)
    parser.add_argument('--Length',type = int)
    parser.add_argument('--Thread',type = str)
    args = parser.parse_args()
    
    path1 = args.path1
    path2 = args.path2
    savebedpath = args.savebedpath
    savecsvpath = args.savecsvpath
    tmpfilepath = args.tmpfilepath
    Start_pos = args.Start_pos
    End_pos = args.End_pos
    Thread = args.Thread
    Length = args.Length
    
    RegInterval = range(Start_pos,End_pos) 
    # Note that it is 0-based. 
    # It will process through exactly Start_pos to End_pos-1 (i.e. range(Start_pos,End_pos))
# path1: path for GEMs (i.e. ___ALL.region.PEanno)
# path2: path for Region (i.e. ____PETcnt_G9.motifannot)
# savebedpath: path for saving extracted GEMs in .bed
# savecsvpath: path for saving summary table in .csv
# tmpfilepath: path for saving tmpfiles produced by pybedtool, a directory
# Thread: for naming the csv file. (i.e. '0')
# Length: Length of extension. Default = 4000 (int)
    mainfunc(path1,path2,savebedpath,savecsvpath,tmpfilepath,RegInterval,Thread,Length)

