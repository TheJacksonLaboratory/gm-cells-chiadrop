#! /bin/sh
#
# auto_aggregate.sh
# Copyright (C) 2020 jianhao2 <jianhao2@illinois.edu>
#
# Distributed under terms of the MIT license.
#


trap "exit" INT TERM ERR
trap "kill 0" EXIT

# p2save_dir=../results/
p2save_dir=/Users/qizai/Dropbox/projects/Minji_project_summary/new_data_0422_2020/new_results_0526/
nbins=1000

expr_name=LHG0052H.for.BROWSER.sorted.bedgraph
python from_bin_matrix_to_plots.py --expr_name $expr_name \
    --p2save_dir $p2save_dir \
    --nbins $nbins \
    --norm 1 &


# for expr_name in $( cat /Users/qizai/Dropbox/projects/Minji_project_summary/new_data_0422_2020/new_results_0526/new_datasets )
# do
#     echo $expr_name
#     python from_bin_matrix_to_plots.py --expr_name $expr_name \
#         --p2save_dir $p2save_dir \
#         --nbins $nbins \
#         --norm 1 &
# done

wait
