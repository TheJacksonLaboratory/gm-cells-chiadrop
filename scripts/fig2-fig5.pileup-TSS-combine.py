import cooler, sys, pickle, matplotlib
matplotlib.use('Agg')
from multiprocess import Pool
import numpy as np
import matplotlib.pyplot as plt
from skimage import transform
from matplotlib.colors import Normalize

class PiecewiseNorm(Normalize):
    def __init__(self, levels, clip=False):
        # the input levels
        self._levels = np.sort(levels)
        # corresponding normalized values between 0 and 1
        self._normed = np.linspace(0, 1, len(levels))
        Normalize.__init__(self, None, None, clip)

    def __call__(self, value, clip=None):
        # linearly interpolate to get the normalized value
        return np.ma.masked_array(np.interp(value, self._levels, self._normed))
        

def parse_tad(fil, pre):

    tad = {}
    with open(fil, 'r') as source:
        for line in source:
            parse = line.rstrip().split()
            chrom = parse[0].lstrip(pre)
            start, end = int(parse[1]), int(parse[2])
            if not chrom in tad:
                tad[chrom] = []
            tad[chrom].append((start, end))
    
    for c in tad:
        tad[c].sort()
    
    return tad

def _expected_core(pool_args):
    
    hic_pool, c, maxdis, balance = pool_args

    expected = {} # average over each genomic distance
    hic = hic_pool.matrix(balance=balance, sparse=True).fetch(c)

    tmp = np.array(hic.sum(axis=0)).ravel() > 2 # filter out gap regions
    n = hic.shape[0]
    maxdis = min(n-1, maxdis)
    # Assign values for each genomic distance
    for i in range(1, maxdis+1):
        valid = tmp[:-i] * tmp[i:]
        current = hic.diagonal(i)[valid]
        if current.size > 0:
            expected[i] = [current.sum(), current.size]
    
    return expected

def calculate_expected(uri, chroms, maxdis=2000000, balance=True, nproc=1):

    # B: Block Bias, constant for each copy number pair
    hic_pool = cooler.Cooler(uri)
    res = hic_pool.binsize
    maxdis = maxdis // res
    args = []
    for c in chroms:
        args.append((hic_pool, c, maxdis, balance))
    
    # Allocate processes
    if nproc==1:
        results = list(map(_expected_core, args))
    else:
        pool = Pool(nproc)
        results = pool.map(_expected_core, args)
        pool.close()
        pool.join()

    expected = {}
    for i in range(1, maxdis+1):
        nume = 0
        denom = 0
        for extract in results:
            if i in extract:
                nume += extract[i][0]
                denom += extract[i][1]
        if nume > 0:
            expected[i] = nume / denom
    
    return expected


def distance_normalize(uri, tadfil, expected, remove_pre='chr', max_apart=2000000, min_size=200000,
    fixed_shape=(20, 20), reverse=False):

    lib = cooler.Cooler(uri)
    tads = parse_tad(tadfil, remove_pre)
    chroms = tads.keys()

    arrs = []
    for c in chroms:
        for t in tads[c]:
            if t[1] - t[0] < min_size:
                continue
            if t[1] - t[0] > max_apart:
                continue
            M = lib.matrix(balance=False, sparse=False).fetch((c, t[0], t[1]))
            idx = np.arange(len(M))
            M[idx, idx] = 1 # mask the diagonal
            for i in range(1, len(M)):
                if i in expected:
                    exp = expected[i]
                else:
                    exp = expected[max(expected)]
                M[idx[:-i], idx[i:]] = M[idx[:-i], idx[i:]] / exp
                M[idx[i:], idx[:-i]] = M[idx[i:], idx[:-i]] / exp
            
            new = transform.resize(M, output_shape=fixed_shape)
            new = new /  new.mean() * M.mean()
            if reverse:
                new = new[::-1, ::-1]
            arrs.append(new)
    
    return arrs

def pileup(arrs):

    pool = np.r_[arrs]
    avg = pool.mean(axis=0)

    return avg

uri = sys.argv[1]
tad_fil1 = sys.argv[2] # forward
tad_fil2 = sys.argv[3] # reverse
outfig = sys.argv[4]
chroms = [str(i) for i in range(1, 23)] + ['X']
max_apart = 3000000
min_size = 100000
n_process = 8
out_shape = (61, 61)
expected = calculate_expected(uri, chroms, maxdis=max_apart, balance=False, nproc=n_process)
arrs1 = distance_normalize(uri, tad_fil1, expected, remove_pre='chr', max_apart=max_apart, min_size=min_size, fixed_shape=out_shape, reverse=False)
arrs2 = distance_normalize(uri, tad_fil2, expected, remove_pre='chr', max_apart=max_apart, min_size=min_size, fixed_shape=out_shape, reverse=True)
avg = pileup(arrs1+arrs2)

print(avg.min(), avg.max())
#plt.imshow(avg, interpolation='none', cmap='seismic', vmin=0.5, vmax=3) # CTCF
plt.imshow(avg, interpolation='none', cmap='seismic', vmin=0.5, vmax=5) # RNAPII
#plt.imshow(avg, interpolation='none', cmap='seismic', vmin=1, vmax=9)
# build a non-linear colormap
#plt.imshow(avg, interpolation='none', cmap='RdYlBu_r', norm=PiecewiseNorm([0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.4, 2, 2.8, 3.5]), vmax=3.5, vmin=0.5) # 1-4, ctcf
#plt.imshow(avg, interpolation='none', cmap='RdYlBu_r', norm=PiecewiseNorm([0.63, 0.76, 0.89, 1.02, 1.15, 1.28, 1.7, 2.4, 3.3, 4.5]), vmax=4.5, vmin=0.63) # 1-4, rnapii
#plt.imshow(avg, interpolation='none', cmap='RdYlBu_r', norm=PiecewiseNorm([0.97, 1.27, 1.58, 1.89, 2.19, 2.5, 3.5, 4.9, 6.9, 9.6]), vmax=9.6, vmin=0.97)
plt.tick_params(axis='both', bottom=False, top=False, left=False, right=False,
                labelbottom=False, labeltop=False, labelleft=False, labelright=False)
plt.colorbar()
plt.savefig(outfig, bbox_inches='tight')
plt.close()