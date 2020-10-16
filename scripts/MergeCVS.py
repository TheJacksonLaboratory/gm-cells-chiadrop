#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--Dir_org',type = int)
    parser.add_argument('--savepath',type = str)

    args = parser.parse_args()
    
    Dir_org = args.Dir_org # path to parse all .csv files (i.e. ./directory1/)
    savepath = args.savepath # path to save results, including the file name (i.e. ./directory1/AllResults.csv)
    
    extension = 'csv'
    all_filenames = [i for i in glob.glob(Dir_org+'*.{}'.format(extension))]
#     print(all_filenames)
    combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames],ignore_index = True, axis=0)
    # combined_csv.insert(0,'LoopID',['Left_0','Left_1','Right_0','Right_1','Both_0','Both_1','None_0','None_1','Total','Left Intensity','Right Intensity','Left motif strand','Right motif strand'])
    # print(combined_csv)
    combined_csv.to_csv( savepath, index=False, encoding='utf-8-sig')

