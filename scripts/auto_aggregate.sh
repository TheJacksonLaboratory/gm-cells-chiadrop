#! /bin/sh
#
# auto_aggregate.sh
# Copyright (C) 2020 jianhao2 <jianhao2@illinois.edu>
#
# Distributed under terms of the MIT license.
#

trap "exit" INT TERM ERR
trap "kill 0" EXIT

# expr_name=ChIP-Seq.MED1.rep1.GSM2443457.SRR5139372.rmdup.q30.hg38.bedgraph
p2intensity=/Users/qizai/Dropbox/projects/Minji_project_summary/new_data_0422_2020/intensity_bedgraph/ 
p2loop_tag=/Users/qizai/Dropbox/projects/Minji_project_summary/new_data_0422_2020/original_data/loop_tags
p2chrom=/Users/qizai/Dropbox/projects/Minji_project_summary/new_data_0422_2020/original_data/hg38.chrom.sizes
p2save_dir=/Users/qizai/Dropbox/projects/Minji_project_summary/new_data_0422_2020/new_results_0526
nbins=1000
pseudo=$1

if [ $pseudo = 1 ]
then
    p2loop_annot=/Users/qizai/Dropbox/projects/Minji_project_summary/new_data_0422_2020/original_data/pseudo_loop_annotation_1000_loops_30_random_sample_per_loop_selected.csv
else
    p2loop_annot=/Users/qizai/Dropbox/projects/Minji_project_summary/new_data_0422_2020/original_data/LHG0052H.e500.clusters.cis.bothanchint_G250.PETcnt_G9.motifannot
fi

expr_name=LHG0052H.for.BROWSER.sorted.bedgraph
python aggregate_intensity.py --p2bedgraph $p2intensity \
    --expr_name $expr_name \
    --p2chrom $p2chrom \
    --p2loop_annot $p2loop_annot \
    --p2loop_tag $p2loop_tag \
    --p2save_dir $p2save_dir \
    --nbins $nbins \
    --pseudo $pseudo &
# for expr_name in $( ls -1 $p2intensity )
# do
#     echo $expr_name
#     time python aggregate_intensity.py --p2bedgraph $p2intensity \
#         --expr_name $expr_name \
#         --p2chrom $p2chrom \
#         --p2loop_annot $p2loop_annot \
#         --p2loop_tag $p2loop_tag \
#         --p2save_dir $p2save_dir \
#         --nbins $nbins \
#         --pseudo $pseudo &
# done

wait
