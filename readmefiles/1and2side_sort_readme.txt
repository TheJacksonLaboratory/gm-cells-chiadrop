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
# path1: PEanno file. i.e.: 'GM12878-RNAPII-pooledv2_comp_FDR_0.2_ALL_promext1kbboth.region.PEanno'
# path2_dir: path of region_cr files. i.e.: './Multiregion_exp/Regions/'
# savebedpath: path to save sorted bedfiles. i.e.: './2sides_bedfiles/'
# savecsvpath: path to save tables. i.e.: './tables/'
# tmpfilepath: path to put temp files. i.e.: './tempfiles'
# Start_pos: determine the start position
# End_pos: determine the end position. For processing the entire range, set to -1
# Side: 1 side plot or 2 side plot (should be 1 or 2, integer)
# Multithread: number of multithreading
# MaxCr: maximum number of cr_id 

Output:
1. Sorted complexes files
2. Stats file (.csv). If Side==1: will have List1.csv and List2.csv. If Side == 2: will have List3.csv.