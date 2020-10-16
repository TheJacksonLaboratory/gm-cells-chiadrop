import ipdb
import fnmatch
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--path',type = str)
    parser.add_argument('--savepath',type = str)
    args = parser.parse_args()
    path = args.path
    savepath = args.savepath

# Preprocess the bedfiles
# path: the file that contain all region_cr (i.e.:'example_cohesin_regions_20200424.bed')
# savepath: path to save processed region_cr files (i.e.:'./Regions_exp/')

    File = pd.read_csv(path,sep = '\t',header=None)

    All_Interval_name = File.iloc[:][3]
    All_cr = fnmatch.filter(All_Interval_name, '*_S')
    for i in range(len(All_cr)):
        All_cr[i] = All_cr[i][:-2]

    All_dict = {}
    Now_idx = 0
    for CR in All_cr:
    #     Format: chr, L_start, R_end, region_name, L_end, R_start
        CR_DF = pd.DataFrame([],columns = ['chr', 'L_start', 'R_end', 'region_name', 'L_end', 'R_start'])
        Flag = True
        M_count = 0
        while Flag:
            try:
    #             ipdb.set_trace()
    #             Now_row = File.iloc[Now_idx]
                if len(fnmatch.filter([File.iloc[Now_idx][3]], CR+'_S')) > 0:
                    Row_S = File.iloc[Now_idx]
                    Now_idx += 1
                elif len(fnmatch.filter([File.iloc[Now_idx][3]], CR+'_E')) > 0:
                    Row_E = File.iloc[Now_idx]
                    Now_idx += 1
                    DF_now = pd.DataFrame([[Row_S[0],Row_S[1],Row_E[2],CR+'_SE',Row_S[2],Row_E[1]]],columns = ['chr', 'L_start', 'R_end', 'region_name', 'L_end', 'R_start'])
                    CR_DF = CR_DF.append(DF_now,ignore_index=True)
                elif len(fnmatch.filter([File.iloc[Now_idx][3]], CR+'_M*')) > 0:
                    M_count += 1
                    Row_M = File.iloc[Now_idx]
                    Now_idx += 1
        #             S to M
                    DF_now = pd.DataFrame([[Row_S[0],Row_S[1],Row_M[2],CR+'_M:'+str(M_count)+'S',Row_S[2],Row_M[1]]],columns = ['chr', 'L_start', 'R_end', 'region_name', 'L_end', 'R_start'])
                    CR_DF = CR_DF.append(DF_now,ignore_index=True)
        #            M to E
                    DF_now = pd.DataFrame([[Row_S[0],Row_M[1],Row_E[2],CR+'_M:'+str(M_count)+'E',Row_M[2],Row_E[1]]],columns = ['chr', 'L_start', 'R_end', 'region_name', 'L_end', 'R_start'])
                    CR_DF = CR_DF.append(DF_now,ignore_index=True)
                else:
                    Flag = False
                    CR_DF.to_csv(savepath+'Region_'+CR,sep = '\t',header=False, index = False)
            except:
                CR_DF.to_csv(savepath+'Region_'+CR,sep = '\t',header=False, index = False)
                break



    
            