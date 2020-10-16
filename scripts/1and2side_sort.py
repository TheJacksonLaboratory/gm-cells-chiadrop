import pybedtools
from pybedtools import BedTool
import numpy as np
from tqdm import tqdm
import pandas as pd
import multiprocessing
from multiprocessing import Pool
import functools
import ipdb
import argparse

# These two lines let the pybedtool know where to put the temp files.
# cleanup() will remove all the temp file create from this session in temp file folder 
# Thread = '0'


def Find2side(RawGEMs,Region,savebedpath,Thread = '',MidRegion = '',MidName = ''):
#     ipdb.set_trace()
    Temp = RawGEMs.groupby(g=[1,2,3,5], c=[5,6,7,8], o=['count','collapse','collapse','collapse'])
    # Need this to keep the result of filter! All these files can be manually removed after the job 
    Test = BedTool(Temp).saveas('dummyfiles/Temp'+Thread)
#     ipdb.set_trace()
    Test = pd.read_csv('dummyfiles/Temp'+Thread,sep = '\t',header=None)
    
    TEMP = list(map(int, [Region[1],Region[4],Region[5],Region[2]]))
    Tempstr = ''
    interval = [0,0,0,0]
    interval[0] = TEMP[0]
    interval[1] = TEMP[1]
    interval[2] = TEMP[2]
    interval[3] = TEMP[3]
    StrList = [[] for i in range(len(Test))]

    for i in range(len(Test)):
#             ipdb.set_trace()
        Start = np.fromstring(str(Test.iloc[i][6]), dtype = np.int,  sep =', ' )
        End = np.fromstring(str(Test.iloc[i][7]), dtype = np.int,  sep =', ' )
#         Mcount = 2#Test[i][5].count('P')
        Start.sort()
        End.sort()
        try:
            GEMLen = End[-1]-Start[0]
        except:
            print(Test.iloc[i])
            ipdb.set_trace()
        # chrom, start_min, end_min, GEM ID, #(fragments)  start_max, end_max,
        subStrList = ['' for j in range(len(Start)+1)]
        
        subStrList[0] = Test.iloc[i][0]+' '+ str(Start[0]) + ' ' + str(End[0])+' '+ Test.iloc[i][3] + ' ' + str(len(Start)) + ' ' + str(GEMLen)+' '
        subStrList[-1] =' ' + str(Start[-1]) + ' ' + str(End[-1]) + '\n'
        FlagMid = False
        for j in range(1,len(Start)-1):
            FlagMid = FlagMid or Checkintersect(Start[j],MidRegion[0],End[j],MidRegion[1])
            subStrList[j] = str(Start[j])+','+str(End[j])
        if (not FlagMid) and (End[0]<MidRegion[0]) and (Start[-1]>MidRegion[1]) and (not Checkintersect(Start[0],interval[0],End[0],interval[1])) and (not Checkintersect(Start[-1],interval[2],End[-1],interval[3])):
#                     ipdb.set_trace()
            subStrList[-2] = str('-1,-1')
            StrList[i] = subStrList[0]+','.join(subStrList[1:len(Start)])+subStrList[-1]
#         ipdb.set_trace()
    StrList[:] = [X for X in StrList if X]
#     ipdb.set_trace()
    Tempstr = ''.join(StrList)
    FinalGEMs = BedTool(Tempstr,from_string=True)
    FinalGEMs = FinalGEMs.moveto(savebedpath+MidName+'.bed')
    
    if len(FinalGEMs)>0:
    #     Further filtering
        FinalGEMs = pd.read_csv(savebedpath+MidName+'.bed',sep = '\t',header=None).sort_values(by = [5])
        Rec_start = FinalGEMs.iloc[0][1]
        Rec_end = FinalGEMs.iloc[0][8]
        for row_idx in range(1,len(FinalGEMs)):
            if FinalGEMs.iloc[row_idx][1]<Rec_start and FinalGEMs.iloc[row_idx][8]>Rec_end:
                Rec_start = FinalGEMs.iloc[row_idx][1]
                Rec_end = FinalGEMs.iloc[row_idx][8]
            else:
                FinalGEMs.iat[row_idx,5] = -1

        FinalGEMs = FinalGEMs[FinalGEMs.iloc[:][5] != -1]
        Len = len(FinalGEMs)
        FinalGEMs.to_csv(savebedpath+MidName+'.bed',sep = '\t',header = None, index = None)
    else:
        Len = len(FinalGEMs)

    return Len

# Process the given region, retrive all valid GEMs
def ProcessRegion(RawGEMs,Thread):
    Temp = RawGEMs.groupby(g=[1,2,3,5], c=[5,6,7,8], o=['count','collapse','collapse','collapse'])
    RefineGEMs = Temp.filter(lambda F: int(F[4]) > 1)
    # Need this to keep the result of filter! All these files can be manually removed after the job 
    Test = BedTool(RefineGEMs).saveas('dummyfiles/RefineGEMs_'+Thread)
    if len(Test) == 0:
        return Test
    Test = pd.read_csv('dummyfiles/RefineGEMs_'+Thread,sep = '\t',header=None)
    
    Tempstr = ''
    for i in range(len(Test)):
        Start = np.fromstring(str(Test.iloc[i][6]), dtype = np.int,  sep =', ' )
        End = np.fromstring(str(Test.iloc[i][7]), dtype = np.int,  sep =', ' )
        
#         ipdb.set_trace()
        UNIQUE, COUNTS = np.unique(Test.iloc[i][5].split(','), return_counts=True)
        Dict = dict(zip(UNIQUE, COUNTS))
        if 'P' in Dict:
            Mcount = Dict['P']
        else:
            Mcount = 0
        
#         Mcount = Test[i][5].count('P')
        Start.sort()
        End.sort()
        # chrom, start_min, end_min, GEM ID, #(fragments)  start_max, end_max,
        for j in range(len(Start)):
            if j == 0:
                Tempstr += Test.iloc[i][0]+' '+ str(Start[j]) + ' ' + str(End[j])+' '+ Test.iloc[i][3] + ' ' + str(len(Start)) + ' ' + str(Mcount)+' '
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
        interval[0] = Temp[0]
        interval[1] = Temp[1]
        interval[2] = Temp[2]
        interval[3] = Temp[3]

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
    Temp = list(map(int, [Region[1],Region[4],Region[5],Region[2]]))
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

def mainfunc(path1,path2,savebedpath,savecsvpath,tmpfilepath,RegInterval,cr_id,Thread,Length = 4000):
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
    Region_short = BedTool(path2)

#     # Remove unnecessary entries
#     Region_short = Region.groupby(g=[1,2,6,12,14,20,8,9,16,21], c=[12], o=['count'])
#     Region_short.moveto('Region_short.bed')
#     Region_short = BedTool('Region_short.bed')
    Max_iter = Region_short.count()
    if RegInterval == 'All':
        RegInterval = range(Max_iter)
#     Length = 4000

    List1 = []
    List2 = []
#     ipdb.set_trace()
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
        FinalGEMs = ProcessRegion(Intersection,Thread)
#         ipdb.set_trace()
        # Classify+sort+save 
        if NowRegion[0][3][-2:] == 'SE':
            Count_L0,Count_L1 = SortGEM(FinalGEMs, NowRegion[0],'Left',Length,savebedpath)
            Count_R0,Count_R1 = SortGEM(FinalGEMs, NowRegion[0],'Right',Length,savebedpath)
            Count_B0,Count_B1 = SortGEM(FinalGEMs, NowRegion[0],'Both',Length,savebedpath)
            
            Count_L = Count_L0+Count_L1
            Count_R = Count_R0+Count_R1
            Count_B = Count_B0+Count_B1
            
            CRID = NowRegion[0][3][:-3]
            List1.append([NowRegion[0][3][:-3],'S_to_E',Count_L])
            List1.append([NowRegion[0][3][:-3],'E_to_S',Count_R])
            List1.append([NowRegion[0][3][:-3],'S_and_E',Count_B])
            
            TempList = [NowRegion[0][3][:-3],Count_L+Count_R+Count_B]
        elif NowRegion[0][3][-1] == 'S':
            Count_R0,Count_R1 = SortGEM(FinalGEMs, NowRegion[0],'Right',Length,savebedpath)
            Count_R = Count_R0+Count_R1
            
            MID = NowRegion[0][3][len(CRID)+1:-1]
            List1.append([CRID,MID+'_to_S',Count_R])
            Count_mid = Count_R
        elif NowRegion[0][3][-1] == 'E':
            Count_L0,Count_L1 = SortGEM(FinalGEMs, NowRegion[0],'Left',Length,savebedpath)
            Count_L = Count_L0+Count_L1
            
            MID = NowRegion[0][3][len(CRID)+1:-1]
            List1.append([CRID,MID+'_to_E',Count_L])
            Count_mid += Count_L
#             ipdb.set_trace()
            NowList = TempList.copy()
            NowList.extend([MID,Count_mid])
            List2.append(NowList)
#         Count_N0,Count_N1 = SortGEM(FinalGEMs, NowRegion[0],'None',Length,savebedpath)
#         Total = Count_L0+Count_L1+Count_R0+Count_R1+Count_B0+Count_B1+Count_N0+Count_N1

#         # Write into dictionary
#         Dict[NowRegion[0][3]] = [NowRegion[0][3],Count_L0,Count_L1,Count_L0+Count_L1,(Count_L0+Count_L1)/Total*100,
#                                  Count_R0,Count_R1,Count_R0+Count_R1,(Count_R0+Count_R1)/Total*100,
#                                  Count_B0,Count_B1,Count_B0+Count_B1,(Count_B0+Count_B1)/Total*100,
#                                  Count_N0,Count_N1,Count_N0+Count_N1,(Count_N0+Count_N1)/Total*100,
#                                  Total,Total-(Count_N0+Count_N1),(Total-(Count_N0+Count_N1))/Total*100,
#                                  NowRegion[0][0]+':'+str(NowRegion[0][1])+'-'+str(NowRegion[0][2])]
#         # Clear all temp files for this session
#         pybedtools.helpers.cleanup()

#     RenameCol = {}
#     NewCol = ['LoopID','Left_0','Left_1','Left_Tol','Left_Tol %','Right_0','Right_1','Right_Tol','Right_Tol %',
#               'Both_0','Both_1','Both_Tol','Both_Tol %',
#               'None_0','None_1','None_Tol','None_Tol %','Total','Total-None','Total-None %',
#               'Region']
#     for i, name in enumerate(NewCol):
#         RenameCol[i] = NewCol[i]
    DF1 = pd.DataFrame(List1,columns = ['crID','orientation','# of complexes'])
    DF2 = pd.DataFrame(List2,columns = ['crID','anchorcomp','middleID','loadcomp'])
    
    DF1.to_csv(savecsvpath+'List1_'+cr_id+'.csv',index=False)
    DF2.to_csv(savecsvpath+'List2_'+cr_id+'.csv',index=False)
#     DF = pd.DataFrame.from_dict(Dict,orient = 'index').rename(columns = RenameCol)
#     # savecsvpath = 'Minji_data/Cohesin_results/01ALL/4kbext_dm/'
#     DF.to_csv(savecsvpath+'LRBNstats_'+Thread+'.csv',index=False)


def mainfunc2(path1,path2,savebedpath,savecsvpath,tmpfilepath,RegInterval,cr_id,Thread,Length = 4000):
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
    Region_short = BedTool(path2)

#     # Remove unnecessary entries
#     Region_short = Region.groupby(g=[1,2,6,12,14,20,8,9,16,21], c=[12], o=['count'])
#     Region_short.moveto('Region_short.bed')
#     Region_short = BedTool('Region_short.bed')
    Max_iter = Region_short.count()
    if RegInterval == 'All':
        RegInterval = range(1,Max_iter)
#     Length = 4000
    
    List1 = []
    NowRegion = BedTool(Region_short[0:1]).saveas()
    # Find all fragments that intersect with Nowregion
    Intersection = ChIA_Drop.intersect(NowRegion,wa=True)
    # Append original start/and. Technical purpose for using groupby...
    results = [(f[0],'0','0',f[3],f[4],f[5],f[1],f[2]) for f in Intersection]
    Intersection = BedTool(results)

    # Sort the grouping key!!!! Otherwise the later groupby doesn't work as intended...
    Intersection = Intersection.sort(chrThenScoreA = True).saveas('dummyfiles/Intersection'+Thread)
    
#     ipdb.set_trace()
# Dict = {'Type/loopID': ['Left_0','Left_1','Right_0','Right_1','Both_0','Both_1','None_0','None_1','Total','Left Intensity', 'Right Intensity','Left motif strand', 'Right motif strand']}
    for i in RegInterval:
        TempRegion = BedTool(Region_short[i:i+1]).saveas()
        GEMid = TempRegion[0][3]
        if GEMid[-1] == 'S':
            MidRegion = np.array([TempRegion[0][5],TempRegion[0][2]]).astype(int)
        else:
            continue
        
        Len = Find2side(Intersection,NowRegion[0],savebedpath,Thread,MidRegion,GEMid[:-1])
        List1.append([GEMid[:-1],Len])
#         ipdb.set_trace()
        
    DF1 = pd.DataFrame(List1,columns = ['crID_M:x','# of complexes'])
    
    DF1.to_csv(savecsvpath+'List3_'+cr_id+'.csv',index=False)

def bedfn1(cr_id,path2_dir,path1,savebedpath,savecsvpath,tmpfilepath,RegInterval):
    path2 = path2_dir+'Region_cr'+str(cr_id)
    Thread = str(multiprocessing.current_process().pid)
    mainfunc(path1,path2,savebedpath,savecsvpath,tmpfilepath,RegInterval,str(cr_id),Thread)
#     mainfunc2(path1,path2,savebedpath,savecsvpath,tmpfilepath,RegInterval,str(cr_id),Thread)
    return None

def bedfn2(cr_id,path2_dir,path1,savebedpath,savecsvpath,tmpfilepath,RegInterval):
    path2 = path2_dir+'Region_cr'+str(cr_id)
    Thread = str(multiprocessing.current_process().pid)
#     mainfunc(path1,path2,savebedpath,savecsvpath,tmpfilepath,RegInterval,str(cr_id),Thread)
    mainfunc2(path1,path2,savebedpath,savecsvpath,tmpfilepath,RegInterval,str(cr_id),Thread)
    return None
    
def Metafunc(path1,path2_dir,savebedpath,savecsvpath,tmpfilepath,RegInterval,Side,
             Multithread=20,MaxCr=3015):
# path1: PEanno file. i.e.: 'GM12878-RNAPII-pooledv2_comp_FDR_0.2_ALL_promext1kbboth.region.PEanno'
# path2_dir: path of region_cr files. i.e.: './Multiregion_exp/Regions/'
# savebedpath: path to save sorted bedfiles. i.e.: './2sides_bedfiles/'
# savecsvpath: path to save tables. i.e.: './tables/'
# tmpfilepath: path to put temp files. i.e.: './tempfiles'
# Side: 1 side plot or 2 side plot (should be 1 or 2, integer)
# Multithread: number of multithreading
# MaxCr: maximum number of cr_id 
    if Side == 1:
        test_bedfn = functools.partial(bedfn1,path2_dir = path2_dir,path1 = path1,
                                       savebedpath = savebedpath,savecsvpath = savecsvpath,
                                       tmpfilepath = tmpfilepath,RegInterval = RegInterval)
    else:
        test_bedfn = functools.partial(bedfn2,path2_dir = path2_dir,path1 = path1,
                                       savebedpath = savebedpath,savecsvpath = savecsvpath,
                                       tmpfilepath = tmpfilepath,RegInterval = RegInterval)
    with Pool(Multithread) as p:
        r = list(tqdm(p.imap(test_bedfn, range(1,MaxCr)), total = MaxCr-1))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--path1',type = str)
    parser.add_argument('--path2_dir',type = str)
    parser.add_argument('--savebedpath',type = str)
    parser.add_argument('--savecsvpath',type = str)
    parser.add_argument('--tmpfilepath',type = str)
    
    parser.add_argument('--Start_pos',type = int)
    parser.add_argument('--End_pos',type = int)
    parser.add_argument('--Side',type = int)
    parser.add_argument('--Multithread',type = int)
    parser.add_argument('--MaxCr',type = int)
    args = parser.parse_args()
    
    path1 = args.path1
    path2_dir = args.path2_dir
    savebedpath = args.savebedpath
    savecsvpath = args.savecsvpath
    tmpfilepath = args.tmpfilepath
    Start_pos = args.Start_pos
    End_pos = args.End_pos
    Side = args.Side
    Multithread = args.Multithread
    MaxCr = args.MaxCr
    
    if End_pos <0:
        RegInterval = 'All'
    else:
        RegInterval = range(Start_pos,End_pos) 
#     # Note that it is 0-based. 
#     # It will process through exactly Start_pos to End_pos-1 (i.e. range(Start_pos,End_pos))


    Metafunc(path1,path2_dir,savebedpath,savecsvpath,tmpfilepath,RegInterval,Side,
             Multithread,MaxCr)