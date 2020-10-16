import os
import pandas as pd
import numpy as np
import time
import argparse

def find_gem_intersection_with_region_add_tag(coord_str, region_start, region_end, 
                                              left_motif_region = None, right_motif_region = None,
                                             gem_id = None, df_region = None, fast_mode = True):
    '''
    similar to find_gem_intersection_with_region.
    first find all complexes with intersected fragments.
    Then determine if the complex has overlap with either side of the motif.
    '''
    # 1. find interesected complexes and selected only the interesected fragments.
    # split by fragments, keep the list_of_frag_coord for later concatenation.
    list_of_frag_coord = coord_str.split(';')
    
    # get rid of chr<>, split number into a tuple.
    tmp_list = [x.split(':')[1] for x in list_of_frag_coord]
    tmp_list = [x.split('-') for x in tmp_list]
    
    # order the fragment by starting site.
    tmp_list.sort(key = lambda x: int(x[0]))
    
    # convert fragment to integer. 
    # Note that end site has '(x)' at the end.
    tmp_list = [[int(x[0]), int(x[1].split('(')[0])] for x in tmp_list]
    frag_coord_filtered_array = np.array(tmp_list)
    
    # test condidtion: >= region_start, and <= region_end. any frag has start/end falls in
    # region is considered as intersected.
    intersect_indicator = (frag_coord_filtered_array >= region_start) & (frag_coord_filtered_array <= region_end)
    intersect_indicator = intersect_indicator.any(axis = 1)
    
    # check if at least two fragment intersected.
    if intersect_indicator.sum() >= 2:
        # selected_frag_idx will look like: (array([idx1, idx2]),)
        selected_frag_idx = np.where(intersect_indicator == 1)[0]
    else:
        selected_frag_idx = []
        
    
    # selected_frag_coord used for numerical computation.
    selected_frag_coord_filtered_array = [frag_coord_filtered_array[x] for x in selected_frag_idx]
    # select_list_of_frag_coord is a string.
    select_list_of_frag_coord = ';'.join([list_of_frag_coord[x] for x in selected_frag_idx])
    num_intersected_frag = len(selected_frag_idx)
    
    if select_list_of_frag_coord == '':
        return select_list_of_frag_coord, num_intersected_frag, ''
    
    # 2. When intersected. add tag to complex based on overlapping with motif.
    # there are basically 2 different approaches. 
    # i) fast_mode == True, look at left/right_motif_region, find overlap with first and last fragment.
    # ii) fast_mode == False, query from df_region by GEM_ID. and get overlapping condition of all fragments.
    left_most_frag = selected_frag_coord_filtered_array[0]
    right_most_frag = selected_frag_coord_filtered_array[-1]
    if fast_mode == True:
        if left_motif_region is None or right_motif_region is None:
            raise ValueError('Fast mode require left/right motif region in [start, end]')
        # find overlap of two intervals.
        # <---->
        #    <----->
        # or
        # <-------->
        #   <--->
        # exception:  <---->
        #                     <--->        
        # =======
        # EXTEND left/right motif region.
        expand_1kbs = np.array([-1000, 1000])
        left_motif_region_expand = left_motif_region + 4 * expand_1kbs
        right_motif_region_expand = right_motif_region + 4 * expand_1kbs

        left_side_flag = ((min(left_most_frag[1], left_motif_region_expand[1]) - 
                          max(left_most_frag[0], left_motif_region_expand[0])) > 0)       
        right_side_flag = ((min(right_most_frag[1], right_motif_region_expand[1]) - 
                          max(right_most_frag[0], right_motif_region_expand[0])) > 0)
    else:
        # query the df_region file. SUPER SLOW
        if gem_id is None or df_region is None:
            raise ValueError('non-Fast mode require gem_id and df_region file(fragment with annotation)')
        all_frag_of_gem = df_region[df_region['GEM_ID'] == gem_id]
        left_most_frag_annot = all_frag_of_gem[all_frag_of_gem['start'] == left_most_frag[0]]
        right_most_frag_annot = all_frag_of_gem[all_frag_of_gem['start'] == right_most_frag[0]]
        
        left_side_flag = (left_most_frag_annot['CTCT_annot'] == 'P').values[0]
        right_side_flag = (right_most_frag_annot['CTCT_annot'] == 'P').values[0]
        
    if left_side_flag & right_side_flag:
        complex_direction_tag = 'Both'
    elif left_side_flag:
        complex_direction_tag = 'Left'
    elif right_side_flag:
        complex_direction_tag = 'Right'
    else:
        complex_direction_tag = 'None'
    
    return select_list_of_frag_coord, num_intersected_frag, complex_direction_tag
    


def find_intersect_for_list_of_region(p2loop_annot, p2chiadrop_complex, p2saved_df, 
                                      test_idx_list, position = 'middle'):
    # Read loop annotation file for loop region and motif region.
    df_loop_annot = pd.read_csv(p2loop_annot, sep = '\t',
                               names = ['left_chr', 'left_start', 'left_end', 
                                   'right_chr', 'right_start', 'right_end',
                              'PET count', 'left_max_intensity', 'right_max_intensity', 
                              'left_max_index', 'right_max_index', 'loop_ID',
                               'left_motif_chr', 'left_motif_start', 'left_motif_end', 
                               'left_motif_strand', 'left_distance',
                              'right_motif_chr', 'right_motif_start', 'right_motif_end', 
                              'right_motif_strand', 'right_distance'])
    
    # Read chiadrop master file. Each line is a complex, with list of fragments at the end.
    column_of_interest = ['GEM_ID', 'GEM_coord', 'GEM_span', 
            'Frag_number', 'List_of_frag_coord']
    df_complex = pd.read_csv(p2chiadrop_complex, sep = '\t')[column_of_interest]

    region_count = 0
    test_region_list = ['lp{}'.format(x) for x in test_idx_list]
    
    assert position in ['middle', 'right', 'left']

    t1 = t_start = time.time()
    for test_region in test_region_list:
        saved_file_name = '{}_{}_100_region.csv'.format(test_region, position)
        #p2saved_df = os.path.join(p2saved_df, saved_file_name)

        region_loop_idx = df_loop_annot['loop_ID'] == test_region
        # region_left_anchor_start = df_loop_annot[region_loop_idx]['left_start'].values[0]
        # region_left_anchor_end = df_loop_annot[region_loop_idx]['left_end'].values[0]
        # region_right_anchor_start = df_loop_annot[region_loop_idx]['right_start'].values[0]
        # region_right_anchor_end = df_loop_annot[region_loop_idx]['right_end'].values[0]
        region_chrom_name = df_loop_annot[region_loop_idx]['left_chr'].values[0]

        left_motif_region = df_loop_annot[region_loop_idx][['left_motif_start', 'left_motif_end']].values[0]
        right_motif_region = df_loop_annot[region_loop_idx][['right_motif_start', 'right_motif_end']].values[0]

        if np.all(left_motif_region == [-1, -1]) or np.all(right_motif_region == [-1, -1]):
            # if any side of the motif are missing. continue.
            continue
        
        # Change left/right_motif_region to be a evenly spreaded range around 
        # region start and region end.
        motif_buffer = 4000
        extension = motif_buffer * np.array([-1, 1])
        left_motif_ext_start, left_motif_ext_end = left_motif_region + extension
        right_motif_ext_start, right_motif_ext_end = right_motif_region + extension
        region_span = right_motif_ext_end - left_motif_ext_start
        
        if position == 'left':
            # when motif_region is reversed, there is no way to find overlapping.
            # min(.., region[1]) - max(..., region[0]) > 0
            # here, min(..., region[1]) <= 0. max(..., region[0]) >= 1. 
            # the condition is always negative.
            right_motif_region = left_motif_region
            left_motif_region = [4001,-4000]
            right_motif_ext_end = left_motif_ext_end
            left_motif_ext_start = right_motif_ext_end - region_span
            
        elif position == 'right':
            left_motif_region = right_motif_region
            right_motif_region = [4001,-4000]
            left_motif_ext_start = right_motif_ext_start
            right_motif_ext_end = left_motif_ext_start + region_span
        
#         ipdb.set_trace()
        selected_GEM_idx = (df_complex['GEM_coord'].apply(lambda x: x.split(':')[0]) == region_chrom_name)
        df_selected_complex = df_complex[selected_GEM_idx]

        fast_filter_fn = lambda x: find_gem_intersection_with_region_add_tag(x, left_motif_ext_start, 
                                                                             right_motif_ext_end,
                                                                            left_motif_region, right_motif_region,
                                                                            fast_mode=True)
        
        frag_coord, num_frag, complex_tag = zip(*df_selected_complex['List_of_frag_coord'].apply(fast_filter_fn))
        df_selected_complex_all = df_selected_complex.assign(**{'Intersected_list_of_frag_coordinates': frag_coord,
                                                                'Intersected_fragment_number': num_frag,
                                                               'complex_direction': complex_tag})
        selected_idx = (df_selected_complex_all['Intersected_list_of_frag_coordinates'] != '')
        df_selected_complex_all = df_selected_complex_all[selected_idx]
        # if position != 'middle':
        #     overlap_idx = (df_selected_complex_all['complex_direction'] != 'None')
        #     df_selected_complex_all = df_selected_complex_all[overlap_idx]
        df_selected_complex_all.to_csv(p2saved_df+saved_file_name)

        region_count += 1
        if region_count % 100 == 0:
            t2 = time.time()
            print('   100 regions finished in {:.2f}s'.format(t2 - t1))
            t1 = time.time()
    t_end = time.time()
    print('-' * 7)
    print('in total {} regions finished in {:.2f}min'.format(region_count, (t_end - t_start)/60))


if __name__ == '__main__':
    random_idx_list = np.random.permutation(range(10000))[:10]
    region_list = ['lp{}'.format(x) for x in random_idx_list]

    parser = argparse.ArgumentParser()
    parser.add_argument('--start', type = int,
            default = 0)
    parser.add_argument('--end', type = int,
            default = 99)
    parser.add_argument('--pos', type = str,
            choices = ['left', 'right', 'middle'],
            default = 'middle')
    parser.add_argument('--p2loop_annot', type = str,
            default = '../data/LHG0052H.e500.clusters.cis.bothanchint_G250.PETcnt_G9.motifannot')
    parser.add_argument('--p2complex_master', type = str,
            default = '../data/miasig_master/CDH0002NR_hg38_CTCF_FDR_0.1_pseudoGEM_5000_enrichTest_master.txt')
    parser.add_argument('--p2saved_df', type = str)

    args = parser.parse_args()
    
    start, end = args.start, args.end
    p2loop_annot = args.p2loop_annot
    p2chiadrop_complex = args.p2complex_master
    p2saved_df = args.p2saved_df
    position = args.pos
    idx_list = range(start, end)

    try:
        if not os.path.isdir(p2saved_df):
            raise ValueError('{} is not a path'.format(p2saved_df))
    except:
        raise TypeError('p2saved_df should be a string type')

    find_intersect_for_list_of_region(p2loop_annot, p2chiadrop_complex, p2saved_df, idx_list, position)
