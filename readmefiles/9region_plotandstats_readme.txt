# Load Python 
module load python/3.6.6
module load numpy/1.18.1
module load pandas/1.0.1
module load matlibplot/3.1.3
module load pybedtools/0.8.1
module load tqdm/4.41.1
module load scipy/1.3.2
module load ipdb/0.12.3
module load glob
module load time
module load functools
module load warnings

Input:
# bedfiledir: savebedpath in 9region_sort.py
# saveplotpath: path to save plots
# regionfilepath: 9 region file
# savecsvpath: path to save tables
# Fname: save plots/table name (i.e. CTCF + one of the 9 region name)
# Tname: title of plots (i.e. CTCF + one of the 9 region name)
# Colorcode: colorcode for ploting fragments

Output:
1. plots in the given saveplotpath.
2. Stats (.csv file) in the give savecsvpath.