# expr_name=SHG0180-181-182NR_hg38_cohesin_FDR_0.1_PASS
# expr_name=CDH0002NR_hg38_CTCF_FDR_0.1_PASS_pseudo_v4
expr_name=$1

pseudo=$2

data_path=/Users/qizai/Dropbox/projects/Minji_project_summary/new_data_0422_2020/

if [ $pseudo = 1 ]
then
    # p2loop_annot=../../results/pseudo_loop_annotataion_1000_loops_5_random_sample_per_loop.csv
    p2loop_annot=${data_path}/original_data/pseudo_loop_annotation_1000_loops_30_random_sample_per_loop.csv
    p2loop_tag=.
    p2intersected_complex=$data_path
    p2binning_folder=${data_path}/new_results_0514/pseudo_loop_results/bin_complex_for_left_middle_right_aggregate/
    # expr_name=${expr_name}_pseudo/
else
    p2loop_annot=${data_path}/original_data/LHG0052H.e500.clusters.cis.bothanchint_G250.PETcnt_G9.motifannot
    p2loop_tag=${data_path}/original_data/loop_tags
    p2intersected_complex=${data_path}
    p2binning_folder=${data_path}/new_results_0514/bin_complex_for_left_middle_right_aggregate/
fi

python wrappter_for_bin_found_complexes_and_aggregate_matrix.py \
    --p2loop_annot $p2loop_annot \
    --p2loop_tag $p2loop_tag \
    --p2intersected_complex_folder $p2intersected_complex \
    --p2binning_results_saved_folder $p2binning_folder \
    --expr_name $expr_name \
    --pseudo $pseudo
