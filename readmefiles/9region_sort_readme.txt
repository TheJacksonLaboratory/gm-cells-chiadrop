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
# path1: CTCF, RNAPII or Cohesin file (i.e.: './GM12878-CTCF-pooled_comp_FDR_0.2_PASS_motifext4kbboth.region.PEanno')
# path2: 9 regions file (i.e.:'./RNAPII-peaks-anchor_mt-forward_TSS-forward_sem-annot_20200711.bed')
# savebedpath: path for saving extracted GEMs in .bed
# MultiThread: Number of multithread

Output:
bedfiles in name: savebedpath+'Region'+str(LoopID)+'_'+Type+'.bed'