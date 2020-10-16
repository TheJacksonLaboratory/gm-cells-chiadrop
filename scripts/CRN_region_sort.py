import pybedtools
from pybedtools import BedTool
import numpy as np
from tqdm import tqdm
import pandas as pd
import multiprocessing
from multiprocessing import Pool
import functools
# import ipdb
import argparse

# These two lines let the pybedtool know where to put the temp files.
# cleanup() will remove all the temp file create from this session in temp file folder 
# Thread = '0'

# Process the given region, retrive all valid GEMs
def ProcessRegion(ChIA_Drop,Loop,savebedpath):
    
    List1 = [] #Left, plot1
    List2 = [] #Right, plot2
    List3 = [] #Both, plot3
    List4 = [] #left of Left, plot4
    List5 = [] #right of Right, plot5
    
    # Extract valid fragments
    Lmostsite = Loop['First_start']-Loop['Second_start']+Loop['First_end']
    Rmostsite = Loop['Second_end']+Loop['Second_start']-Loop['First_end']
#     ipdb.set_trace()
    Meta_ValidFrags = ChIA_Drop[(ChIA_Drop['start']>Lmostsite) & (ChIA_Drop['end']<Rmostsite) & (ChIA_Drop['chr']==Loop['First_chr'])]
    
#     For plot 1-3 (Left, Right, Both)
    ValidFrags = Meta_ValidFrags[(Meta_ValidFrags['start']>Loop['First_start']) & (Meta_ValidFrags['end']<Loop['Second_end'])]
    # Extract complexes name
    ComplexName = ValidFrags['Complex_name'].drop_duplicates()
    
    for ii in range(len(ComplexName)):
        complexName = ComplexName.iloc[ii]
        Complex = ValidFrags[ValidFrags['Complex_name'] == complexName]
        Complex = Complex.sort_values(by = ['start'])
        if len(Complex)<2:
            continue
            
        if Checkintersect(Complex.iloc[0]['start'],Loop['First_start'],Complex.iloc[0]['end'],Loop['First_end']) and Checkintersect(Complex.iloc[-1]['start'],Loop['Second_start'],Complex.iloc[-1]['end'],Loop['Second_end']):
            Tmpstr = str(Complex.iloc[0]['chr'])+' '+str(Complex.iloc[0]['start'])+' '+str(Complex.iloc[0]['end'])+' '+complexName+' '+str(len(Complex))+' '
            for ff in range(1,len(Complex)-1):
                Tmpstr += str(Complex.iloc[ff]['start'])+','+str(Complex.iloc[ff]['end'])+','
            Tmpstr += str('-1,-1')+' '+str(Complex.iloc[len(Complex)-1]['start'])+' '+str(Complex.iloc[len(Complex)-1]['end'])+ '\n'
            List3.append(Tmpstr)
        elif Checkintersect(Complex.iloc[0]['start'],Loop['First_start'],Complex.iloc[0]['end'],Loop['First_end']):
            Tmpstr = str(Complex.iloc[0]['chr'])+' '+str(Complex.iloc[0]['start'])+' '+str(Complex.iloc[0]['end'])+' '+complexName+' '+str(len(Complex))+' '
            for ff in range(1,len(Complex)-1):
                Tmpstr += str(Complex.iloc[ff]['start'])+','+str(Complex.iloc[ff]['end'])+','
            Tmpstr += str('-1,-1')+' '+str(Complex.iloc[len(Complex)-1]['start'])+' '+str(Complex.iloc[len(Complex)-1]['end'])+ '\n'
            List1.append(Tmpstr)
        elif Checkintersect(Complex.iloc[-1]['start'],Loop['Second_start'],Complex.iloc[-1]['end'],Loop['Second_end']):
            Tmpstr = str(Complex.iloc[0]['chr'])+' '+str(Complex.iloc[0]['start'])+' '+str(Complex.iloc[0]['end'])+' '+complexName+' '+str(len(Complex))+' '
            for ff in range(1,len(Complex)-1):
                Tmpstr += str(Complex.iloc[ff]['start'])+','+str(Complex.iloc[ff]['end'])+','
            Tmpstr += str('-1,-1')+' '+str(Complex.iloc[len(Complex)-1]['start'])+' '+str(Complex.iloc[len(Complex)-1]['end'])+ '\n'
            List2.append(Tmpstr)

#   Plot 4, left element to left  
    ValidFrags = Meta_ValidFrags[(Meta_ValidFrags['start']>Lmostsite) & (Meta_ValidFrags['end']<Loop['First_end'])]
    ComplexName = ValidFrags['Complex_name'].drop_duplicates()
    for ii in range(len(ComplexName)):
        complexName = ComplexName.iloc[ii]
        Complex = ValidFrags[ValidFrags['Complex_name'] == complexName]
        Complex = Complex.sort_values(by = ['start'])
        if len(Complex)<2:
            continue
        if Checkintersect(Complex.iloc[0]['start'],Lmostsite,Complex.iloc[0]['end'],Loop['First_start']) and Checkintersect(Complex.iloc[-1]['start'],Loop['First_start'],Complex.iloc[-1]['end'],Loop['First_end']):
            Tmpstr = str(Complex.iloc[0]['chr'])+' '+str(Complex.iloc[0]['start'])+' '+str(Complex.iloc[0]['end'])+' '+complexName+' '+str(len(Complex))+' '
            for ff in range(1,len(Complex)-1):
                Tmpstr += str(Complex.iloc[ff]['start'])+','+str(Complex.iloc[ff]['end'])+','
            Tmpstr += str('-1,-1')+' '+str(Complex.iloc[len(Complex)-1]['start'])+' '+str(Complex.iloc[len(Complex)-1]['end'])+ '\n'
            List4.append(Tmpstr)

#   Plot 5, right element to right  
    ValidFrags = Meta_ValidFrags[(Meta_ValidFrags['start']>Loop['Second_start']) & (Meta_ValidFrags['end']<Rmostsite)]
    ComplexName = ValidFrags['Complex_name'].drop_duplicates()
    for ii in range(len(ComplexName)):
        complexName = ComplexName.iloc[ii]
        Complex = ValidFrags[ValidFrags['Complex_name'] == complexName]
        Complex = Complex.sort_values(by = ['start'])
        if len(Complex)<2:
            continue
        if Checkintersect(Complex.iloc[0]['start'],Loop['Second_start'],Complex.iloc[0]['end'],Loop['Second_end']) and Checkintersect(Complex.iloc[-1]['start'],Loop['Second_end'],Complex.iloc[-1]['end'],Rmostsite):
            Tmpstr = str(Complex.iloc[0]['chr'])+' '+str(Complex.iloc[0]['start'])+' '+str(Complex.iloc[0]['end'])+' '+complexName+' '+str(len(Complex))+' '
            for ff in range(1,len(Complex)-1):
                Tmpstr += str(Complex.iloc[ff]['start'])+','+str(Complex.iloc[ff]['end'])+','
            Tmpstr += str('-1,-1')+' '+str(Complex.iloc[len(Complex)-1]['start'])+' '+str(Complex.iloc[len(Complex)-1]['end'])+ '\n'
            List5.append(Tmpstr)
    
#     ipdb.set_trace()
    Tempstr = ''.join(List1)            
    FinalGEMs = BedTool(Tempstr,from_string=True).saveas(savebedpath+str(Loop['RAID_Name'])+'_'+str(Loop['First_Name'])+'_'+str(Loop['Second_Name'])+'_'+'Left'+'.bed')
    Tempstr = ''.join(List2)            
    FinalGEMs = BedTool(Tempstr,from_string=True).saveas(savebedpath+str(Loop['RAID_Name'])+'_'+str(Loop['First_Name'])+'_'+str(Loop['Second_Name'])+'_'+'Right'+'.bed')
    Tempstr = ''.join(List3)            
    FinalGEMs = BedTool(Tempstr,from_string=True).saveas(savebedpath+str(Loop['RAID_Name'])+'_'+str(Loop['First_Name'])+'_'+str(Loop['Second_Name'])+'_'+'Both'+'.bed')
    Tempstr = ''.join(List4)            
    FinalGEMs = BedTool(Tempstr,from_string=True).saveas(savebedpath+str(Loop['RAID_Name'])+'_'+str(Loop['First_Name'])+'_'+str(Loop['Second_Name'])+'_'+'LoL'+'.bed')
    Tempstr = ''.join(List5)            
    FinalGEMs = BedTool(Tempstr,from_string=True).saveas(savebedpath+str(Loop['RAID_Name'])+'_'+str(Loop['First_Name'])+'_'+str(Loop['Second_Name'])+'_'+'RoR'+'.bed')
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
    Loop = Loops.iloc[i] 
    ProcessRegion(ChIA_Drop,Loop,savebedpath)
    return None

def mainfunc(path1,path2,savebedpath,Multithread=20):
# path1: path for GEMs (i.e. 'GM12878-cohesin-pooled_comp_FDR_0.2_PASS_motifext4kbboth.region.PEanno')
# path2: path for Region (i.e. 'CRN-loops.G100kb.20200705.ext4kbboth.bedpe')
# savebedpath: path for saving extracted GEMs in .bed (i.e. './bedfiles/')
# Multithread: number of multithreading

    header_list = ['chr', 'start','end','?','Complex_name','??']
    ChIA_Drop = pd.read_csv(path1,sep = '\t',names = header_list)

    header_list = ['First_chr', 'First_start','First_end',
               'Second_chr','Second_start','Second_end',
               'First_Name','Second_Name','RAID_Name']
    Loops = pd.read_csv(path2,sep = '\t',names=header_list)


    test_bedfn = functools.partial(bedfn,Loops = Loops,ChIA_Drop = ChIA_Drop,savebedpath = savebedpath,)
    with Pool(Multithread) as p:
        r = list(tqdm(p.imap(test_bedfn, range(len(Loops))), total = len(Loops)))    
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
