#! /bin/sh
#
# run_loop_tag_and_aggregation.sh

python loop_tag_and_aggregation.py \
    --p2loop_file ../../ChIA_PET/LHG0052H.e500.clusters.cis.bothanchint_G250.PETcnt_G9.motifannot \
    --p2bedgraph ../../data/LHG0052H.for.BROWSER.sorted.bedgraph \
    --p2save_loop_tag loop_tag \
    --nbins 1000 \
    --p2chrom_size ../../ChIA_PET/hg38.chrom.sizes \
    --p2agg_stats ChIA_PET_region_file_aggregation_stats \
    --p2binned_intensity_per_loop binned_intensity_per_loop_chia_pet
