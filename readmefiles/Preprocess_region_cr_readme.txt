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
# path: the file that contain all region_cr (i.e.:'example_cohesin_regions_20200424.bed')
# savepath: path to save processed region_cr files (i.e.:'./Regions_exp/')

Output:
region_cr files 