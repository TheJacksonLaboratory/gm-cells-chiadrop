#! /bin/sh
#
# run_generate_random_loop_annotation_file_for_NULL_hypothesis.sh
# Copyright (C) 2020 qizai <jianhao2@illinois.edu>
#
# Distributed under terms of the MIT license.
#


python generate_random_loop_annotation_file_for_NULL_hypothesis.py \
    --p2real_annot ../../../ChIA_PET/LHG0052H.e500.clusters.cis.bothanchint_G250.PETcnt_G9.motifannot_added_labels \
    --p2chrom_size ../../../ChIA_PET/hg38.chrom.sizes \
    --nloops 10 \
    --repeat 5 \
    --p2pseudo_annot pseudo_loop_annot \
