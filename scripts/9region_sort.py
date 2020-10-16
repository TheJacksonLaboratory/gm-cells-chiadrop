import pybedtools
from pybedtools import BedTool
import numpy as np
from tqdm import tqdm
import pandas as pd
import multiprocessing
from multiprocessing import Pool
import functools
import argparse



# Process the given region, retrive all valid GEMs
def ProcessRegion(ChIA_Drop,Loop,savebedpath):
    
#     Loop: three lines, they are
# chr, start, end, startX
# chr, start, end, endX
# chr, start, end, middleX
# X: loop id, start from 0

    List1 = [] #Left, plot1 
    List2 = [] #Right, plot2
    List3 = [] #Both, plot3
    List4 = [] #Mid to Left, plot4
    List5 = [] #Mid and Left, plot5
    List6 = [] #Mid to Right, plot6
    List7 = [] #Mid and Right, plot7
    
    # Extract valid fragments
    Lmostsite = Loop.iloc[0]['S']
    Rmostsite = Loop.iloc[1]['E']
    LoopID = Loop.iloc[0]['Name'][5:]
#     ipdb.set_trace()
    Meta_ValidFrags = ChIA_Drop[(ChIA_Drop['start']>Lmostsite) & (ChIA_Drop['end']<Rmostsite) & (ChIA_Drop['chr']==Loop.iloc[0]['chr'])]
    
#     For plot 1-3 (Left, Right, Both)
    ValidFrags = Meta_ValidFrags
    # Extract complexes name
    ComplexName = ValidFrags['Complex_name'].drop_duplicates()
    
    for ii in range(len(ComplexName)):
        complexName = ComplexName.iloc[ii]
        Complex = ValidFrags[ValidFrags['Complex_name'] == complexName]
        Complex = Complex.sort_values(by = ['start'])
        if len(Complex)<2:
            continue
            
        if Checkintersect(Complex.iloc[0]['start'],Loop.iloc[0]['S'],Complex.iloc[0]['end'],Loop.iloc[0]['E']) and Checkintersect(Complex.iloc[-1]['start'],Loop.iloc[1]['S'],Complex.iloc[-1]['end'],Loop.iloc[1]['E']):
            Tmpstr = str(Complex.iloc[0]['chr'])+' '+str(Complex.iloc[0]['start'])+' '+str(Complex.iloc[0]['end'])+' '+complexName+' '+str(len(Complex))+' '
            for ff in range(1,len(Complex)-1):
                Tmpstr += str(Complex.iloc[ff]['start'])+','+str(Complex.iloc[ff]['end'])+','
            Tmpstr += str('-1,-1')+' '+str(Complex.iloc[len(Complex)-1]['start'])+' '+str(Complex.iloc[len(Complex)-1]['end'])+ '\n'
            List3.append(Tmpstr)
        elif Checkintersect(Complex.iloc[0]['start'],Loop.iloc[0]['S'],Complex.iloc[0]['end'],Loop.iloc[0]['E']):
            Tmpstr = str(Complex.iloc[0]['chr'])+' '+str(Complex.iloc[0]['start'])+' '+str(Complex.iloc[0]['end'])+' '+complexName+' '+str(len(Complex))+' '
            for ff in range(1,len(Complex)-1):
                Tmpstr += str(Complex.iloc[ff]['start'])+','+str(Complex.iloc[ff]['end'])+','
            Tmpstr += str('-1,-1')+' '+str(Complex.iloc[len(Complex)-1]['start'])+' '+str(Complex.iloc[len(Complex)-1]['end'])+ '\n'
            List1.append(Tmpstr)
        elif Checkintersect(Complex.iloc[-1]['start'],Loop.iloc[1]['S'],Complex.iloc[-1]['end'],Loop.iloc[1]['E']):
            Tmpstr = str(Complex.iloc[0]['chr'])+' '+str(Complex.iloc[0]['start'])+' '+str(Complex.iloc[0]['end'])+' '+complexName+' '+str(len(Complex))+' '
            for ff in range(1,len(Complex)-1):
                Tmpstr += str(Complex.iloc[ff]['start'])+','+str(Complex.iloc[ff]['end'])+','
            Tmpstr += str('-1,-1')+' '+str(Complex.iloc[len(Complex)-1]['start'])+' '+str(Complex.iloc[len(Complex)-1]['end'])+ '\n'
            List2.append(Tmpstr)

#   Plot 4-5, mid --- left
    ValidFrags = Meta_ValidFrags[(Meta_ValidFrags['start']>Lmostsite) & (Meta_ValidFrags['end']<Loop.iloc[2]['E'])]
    ComplexName = ValidFrags['Complex_name'].drop_duplicates()
    for ii in range(len(ComplexName)):
        complexName = ComplexName.iloc[ii]
        Complex = ValidFrags[ValidFrags['Complex_name'] == complexName]
        Complex = Complex.sort_values(by = ['start'])
        if len(Complex)<2:
            continue
        if Checkintersect(Complex.iloc[0]['start'],Loop.iloc[0]['S'],Complex.iloc[0]['end'],Loop.iloc[0]['E']) and Checkintersect(Complex.iloc[-1]['start'],Loop.iloc[2]['S'],Complex.iloc[-1]['end'],Loop.iloc[2]['E']):
            Tmpstr = str(Complex.iloc[0]['chr'])+' '+str(Complex.iloc[0]['start'])+' '+str(Complex.iloc[0]['end'])+' '+complexName+' '+str(len(Complex))+' '
            for ff in range(1,len(Complex)-1):
                Tmpstr += str(Complex.iloc[ff]['start'])+','+str(Complex.iloc[ff]['end'])+','
            Tmpstr += str('-1,-1')+' '+str(Complex.iloc[len(Complex)-1]['start'])+' '+str(Complex.iloc[len(Complex)-1]['end'])+ '\n'
            List5.append(Tmpstr)
        elif Checkintersect(Complex.iloc[-1]['start'],Loop.iloc[2]['S'],Complex.iloc[-1]['end'],Loop.iloc[2]['E']):
            Tmpstr = str(Complex.iloc[0]['chr'])+' '+str(Complex.iloc[0]['start'])+' '+str(Complex.iloc[0]['end'])+' '+complexName+' '+str(len(Complex))+' '
            for ff in range(1,len(Complex)-1):
                Tmpstr += str(Complex.iloc[ff]['start'])+','+str(Complex.iloc[ff]['end'])+','
            Tmpstr += str('-1,-1')+' '+str(Complex.iloc[len(Complex)-1]['start'])+' '+str(Complex.iloc[len(Complex)-1]['end'])+ '\n'
            List4.append(Tmpstr)

#   Plot 6-7, mid---right  
    ValidFrags = Meta_ValidFrags[(Meta_ValidFrags['start']>Loop.iloc[2]['S']) & (Meta_ValidFrags['end']<Rmostsite)]
    ComplexName = ValidFrags['Complex_name'].drop_duplicates()
    for ii in range(len(ComplexName)):
        complexName = ComplexName.iloc[ii]
        Complex = ValidFrags[ValidFrags['Complex_name'] == complexName]
        Complex = Complex.sort_values(by = ['start'])
        if len(Complex)<2:
            continue
        if Checkintersect(Complex.iloc[0]['start'],Loop.iloc[2]['S'],Complex.iloc[0]['end'],Loop.iloc[2]['E']) and Checkintersect(Complex.iloc[-1]['start'],Loop.iloc[1]['S'],Complex.iloc[-1]['end'],Loop.iloc[1]['E']):
            Tmpstr = str(Complex.iloc[0]['chr'])+' '+str(Complex.iloc[0]['start'])+' '+str(Complex.iloc[0]['end'])+' '+complexName+' '+str(len(Complex))+' '
            for ff in range(1,len(Complex)-1):
                Tmpstr += str(Complex.iloc[ff]['start'])+','+str(Complex.iloc[ff]['end'])+','
            Tmpstr += str('-1,-1')+' '+str(Complex.iloc[len(Complex)-1]['start'])+' '+str(Complex.iloc[len(Complex)-1]['end'])+ '\n'
            List7.append(Tmpstr)
        elif Checkintersect(Complex.iloc[0]['start'],Loop.iloc[2]['S'],Complex.iloc[0]['end'],Loop.iloc[2]['E']):
            Tmpstr = str(Complex.iloc[0]['chr'])+' '+str(Complex.iloc[0]['start'])+' '+str(Complex.iloc[0]['end'])+' '+complexName+' '+str(len(Complex))+' '
            for ff in range(1,len(Complex)-1):
                Tmpstr += str(Complex.iloc[ff]['start'])+','+str(Complex.iloc[ff]['end'])+','
            Tmpstr += str('-1,-1')+' '+str(Complex.iloc[len(Complex)-1]['start'])+' '+str(Complex.iloc[len(Complex)-1]['end'])+ '\n'
            List6.append(Tmpstr)
    
#     ipdb.set_trace()
    Tempstr = ''.join(List1)            
    FinalGEMs = BedTool(Tempstr,from_string=True).saveas(savebedpath+'Region'+str(LoopID)+'_'+'Left'+'.bed')
    Tempstr = ''.join(List2)            
    FinalGEMs = BedTool(Tempstr,from_string=True).saveas(savebedpath+'Region'+str(LoopID)+'_'+'Right'+'.bed')
    Tempstr = ''.join(List3)            
    FinalGEMs = BedTool(Tempstr,from_string=True).saveas(savebedpath+'Region'+str(LoopID)+'_'+'Both'+'.bed')
    Tempstr = ''.join(List4)            
    FinalGEMs = BedTool(Tempstr,from_string=True).saveas(savebedpath+'Region'+str(LoopID)+'_'+'M2L'+'.bed')
    Tempstr = ''.join(List5)            
    FinalGEMs = BedTool(Tempstr,from_string=True).saveas(savebedpath+'Region'+str(LoopID)+'_'+'MaL'+'.bed')
    Tempstr = ''.join(List6)            
    FinalGEMs = BedTool(Tempstr,from_string=True).saveas(savebedpath+'Region'+str(LoopID)+'_'+'M2R'+'.bed')
    Tempstr = ''.join(List7)            
    FinalGEMs = BedTool(Tempstr,from_string=True).saveas(savebedpath+'Region'+str(LoopID)+'_'+'MaR'+'.bed')
    
    return 

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

def bedfn(i,Loops,ChIA_Drop,savebedpath):
    Loop = Loops.iloc[i:(i+3)] 
#     ipdb.set_trace()
    ProcessRegion(ChIA_Drop,Loop,savebedpath)
    return None

def mainfunc(path1,path2,savebedpath,MultiThread = 20):
# path1 = 'GM12878-CTCF-pooled_comp_FDR_0.2_PASS_motifext4kbboth.region.PEanno' (CTCF, RNAPII or Cohesin file)
# path2 = 'RNAPII-peaks-anchor_mt-forward_TSS-forward_sem-annot_20200711.bed' (9 regions file)
# savebedpath: path for saving extracted GEMs in .bed
# MultiThread: Number of multithread

    header_list = ['chr', 'start','end','?','Complex_name','??']
    ChIA_Drop = pd.read_csv(path1,sep = '\t',names = header_list)

    header_list = ['chr', 'S','E','Name']
    Loops = pd.read_csv(path2,sep = '\t',names=header_list)
    
    test_bedfn = functools.partial(bedfn,Loops,ChIA_Drop = ChIA_Drop,Loops=Loops,savebedpath = savebedpath)
    
    with Pool(MultiThread) as p:
        r = list(tqdm(p.imap(test_bedfn, range(0,len(Loops),3)), total = int(len(Loops)/3)))
    ProcessRegion(ChIA_Drop,Loop,savebedpath) 
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--path1',type = str)
    parser.add_argument('--path2',type = str)
    parser.add_argument('--savebedpath',type = str)
    
    parser.add_argument('--MultiThread',type = int)
    args = parser.parse_args()
    
    path1 = args.path1
    path2 = args.path2
    savebedpath = args.savebedpath
    MultiThread = args.MultiThread

    mainfunc(path1,path2,savebedpath,MultiThread)