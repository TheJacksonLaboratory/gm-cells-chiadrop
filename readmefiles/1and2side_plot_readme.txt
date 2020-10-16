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
#     bedfiledir: path for 1 side sorted complexes files. i.e. './bedfiles/'
#     bedfiledir: path for 2 side sorted complexes files. i.e. './2sides_bedfiles/'
#     saveplotpath: path to save plots. i.e. './Plots/'
#     regionfiledir: path of processed region_cr files.
#     Multithread: number of multithreading
#     MaxCr: maximum number of cr_id

Output:
Both 1 side and 2 side plots to the saveplotpath.