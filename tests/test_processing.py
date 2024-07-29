import pandas as pd
import networkx as nx


from dram2_viz.processing.process_annotations import (
    get_module_step_coverage, make_module_coverage_df, make_module_coverage_frame,
    pairwise, first_open_paren_is_all, split_into_steps, is_ko, make_module_network,
    get_module_coverage, make_etc_coverage_df, make_functional_df, get_phylum_and_most_specific,
    make_product_df, fill_product_dfs, get_annotation_ids_by_row, get_ordered_uniques
)


def test_get_ordered_uniques():
    assert get_ordered_uniques([1, 2, 3]) == [1, 2, 3]
    assert get_ordered_uniques([1, 1, 2, 3]) == [1, 2, 3]
    assert get_ordered_uniques([1, 2, 1, 3]) == [1, 2, 3]


def test_get_ids_from_row():
    in_data = pd.concat([pd.DataFrame({'ko_id': 'K00001,K00003'}, index=['id_set1']),
                         pd.DataFrame({'kegg_hit': 'Some text and then [EC:0.0.0.0]; also [EC:1.1.1.1]'}, index=['id_set2']),
                         pd.DataFrame({'peptidase_family': 'ABC1;BCD2'}, index=['id_set3']),
                         pd.DataFrame({'cazy_best_hit': 'GH4'}, index=['id_set4'])
                        ])
    out_data = get_annotation_ids_by_row(in_data)
    assert out_data['id_set1'] == {'K00001', 'K00003'}
    assert out_data['id_set2'] == {'EC:0.0.0.0', 'EC:1.1.1.1'}
    assert out_data['id_set3'] == {'ABC1', 'BCD2'}
    assert out_data['id_set4'] == {'GH4'}


def test_build_module_net(test_module_net):
    real_module_net = nx.DiGraph(num_steps=0, module_id='M12345', module_name='a module name')
    real_module_net.add_edges_from([('0,0', 'end_step_0'), ('end_step_0', '1,0'),
                                    ('1,0', 'end_step_1'), ('end_step_1', '2,0'),
                                    ('2,0', 'end_step_2')])
    assert nx.is_isomorphic(test_module_net, real_module_net)


def test_get_module_step_coverage(test_module_net):
    test_coverages1 = get_module_step_coverage(set([]), test_module_net)
    assert test_coverages1 == (3, 0, 0, [])
    test_coverages2 = get_module_step_coverage({'K00001', 'K00003'}, test_module_net)
    assert test_coverages2 == (3, 2, 2/3, ['K00001', 'K00003'])
    test_coverages2 = get_module_step_coverage({'K00001', 'K00003', 'K00002', 'K12345'}, test_module_net)
    assert test_coverages2 == (3, 3, 1, ['K00001', 'K00002', 'K00003'])


def test_make_module_coverage_df(test_annotations_df, test_module_net):
    test_module_coverage_df = make_module_coverage_df(test_annotations_df, {'M12345': test_module_net})
    module_coverage_df = pd.DataFrame([['a module name', 3, 1, 1/3, 1, 'K00001', 'gene_3']],
                                      index=['M12345'],
                                      columns=['module_name', 'steps', 'steps_present', 'step_coverage', 'ko_count',
                                               'kos_present', 'genes_present'])
    pd.testing.assert_frame_equal(test_module_coverage_df, module_coverage_df)


def test_make_module_coverage_frame(test_annotations_df, test_module_net, module_coverage_frame):
    test_module_coverage_frame = make_module_coverage_frame(test_annotations_df, {'M12345': test_module_net},
                                                            groupby_column='scaffold')
    pd.testing.assert_frame_equal(test_module_coverage_frame, module_coverage_frame)


def test_pairwise():
    assert list(pairwise([1, 2, 3])) == [(1, 2), (2, 3)]


def test_first_open_paren_is_all():
    assert first_open_paren_is_all('()')
    assert first_open_paren_is_all('(K1+K2-(K3+K5-K4))')
    assert not first_open_paren_is_all('(K1+K2-)K3+K5-K4)')


def test_split_into_steps():
    true_steps = ['K00330', 'K00331+K00332,K00331+K13378,K13380']
    assert split_into_steps('K00330+(K00331+K00332,K00331+K13378,K13380)', '+') == true_steps


def test_is_ko():
    assert is_ko('K00000')
    assert not is_ko('K1')


def test_make_module_network(module_network):
    test_network, _ = make_module_network('K00001+(K00002,K00003+K00013)+K00004')
    assert nx.is_isomorphic(test_network, module_network)


def test_get_module_coverage(module_network):
    module_network.add_edge('K00004', 'end')
    assert (3, 1, 1/3, {'K00001'}, {'K00002', 'K00004'}) == get_module_coverage(module_network, {'K00001'})
    assert (4, 1, 1/4, {'K00003'}, {'K00001', 'K00013', 'K00004'}) == get_module_coverage(module_network, {'K00003'})
    assert (3, 1, 1/3, {'K00002'}, {'K00001', 'K00004'}) == get_module_coverage(module_network, {'K00002', 'K99999'})


def test_make_etc_coverage_df(test_annotations_ids_by_row_df, etc_module_df, etc_coverage_df):
    test_etc_coverage_df = make_etc_coverage_df(etc_module_df, test_annotations_ids_by_row_df, 'scaffold')
    pd.testing.assert_frame_equal(test_etc_coverage_df, etc_coverage_df)


def test_make_functional_df(test_annotations_ids_by_row_df, function_heatmap_form, functional_df):
    test_functional_df = make_functional_df(test_annotations_ids_by_row_df, function_heatmap_form, 'scaffold')
    pd.testing.assert_frame_equal(test_functional_df, functional_df)


# TODO: actually test that the frames are correct, already done above
def test_fill_product_dfs(test_annotations_df, test_module_net, etc_module_df, function_heatmap_form, test_annotations_ids_by_row_df):
    module_nets = {'M12345': test_module_net}
    liquor_dfs = fill_product_dfs(test_annotations_df, module_nets, etc_module_df, function_heatmap_form, test_annotations_ids_by_row_df, groupby_column='scaffold')
    assert len(liquor_dfs) == 3
    assert type(liquor_dfs[0]) is pd.DataFrame
    assert type(liquor_dfs[1]) is pd.DataFrame
    assert type(liquor_dfs[2]) is pd.DataFrame


def test_make_product_df(module_coverage_frame, etc_coverage_df, functional_df):
    liquor_df = pd.DataFrame([[1/3, 1/3, False, True]],
                             index=pd.Index(['scaffold_1'], name='genome'),
                             columns=['a module name', 'Complex I: oxidoreductase', 'Category1: A function',
                                      'Category1: B function'])
    test_liquor_df = make_product_df(module_coverage_frame, etc_coverage_df, functional_df)
    pd.testing.assert_frame_equal(test_liquor_df, liquor_df)


def test_get_phylum_and_most_specific():
    assert get_phylum_and_most_specific('d__Bacteria;p__Bacteroidota;c__;o__;f__;g__;s__') == \
           'p__Bacteroidota;c__'
    assert get_phylum_and_most_specific('d__Archaea;p__;c__;o__;f__;g__;s__') == 'd__Archaea;p__'
    assert get_phylum_and_most_specific('d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae;'
                                        'g__Alistipes;s__') == 'p__Bacteroidota;g__Alistipes'
    assert get_phylum_and_most_specific('d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Enterococcaceae;'
                                        'g__Enterococcus_D;s__Enterococcus_D gallinarum') == \
        'p__Firmicutes;s__Enterococcus_D gallinarum'
