#! /bin/sh
#
# Copyright (C) 2020 jianhao <jianhao2@illinois.edu>
#
# Distributed under terms of the MIT license.
#

start=37
end=38
pos=right
p2loop_annot=/projects/ruan-lab/GM12878_CTCF_ChIA-Drop/fake_loop_annotataion_1000_loops_5_random_sample_per_loop.csv
p2complex_master=/projects/ruan-lab/GM12878_CTCF_ChIA-Drop/miasig_results/CDH0002NR_hg38_CTCF_EnrichTest_FDR_0.1/CDH0002NR_hg38_CTCF_FDR_0.1_pseudoGEM_5000_enrichTest_master.txt
p2saved_df=/projects/kimm/jianhao_chiadrop/fake_loop_samples

python extract_complexes_only_based_on_motif_v4.py --start ${start} \
    --end ${end} \
    --pos ${pos} \
    --p2loop_annot ${p2loop_annot} \
    --p2complex_master ${p2complex_master} \
    --p2saved_df ${p2saved_df}

