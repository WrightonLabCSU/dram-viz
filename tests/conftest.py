import pytest
import pandas as pd
import networkx as nx

from dram2_viz.definitions import DBSETS_COL

from dram2_viz.processing.process_annotations import build_module_net, get_annotation_ids_by_row


@pytest.fixture()
def annotations():
    return pd.DataFrame([['genome', 'K00001'],
                         ['genome', pd.np.NaN]],
                        index=['genome_scaffold_1_1', 'genome_scaffold_1_2'],
                        columns=['fasta', 'ko_id'])


@pytest.fixture()
def genome_summary_frame():
    return pd.DataFrame(pd.DataFrame([['K00001', 'description', 'module1', 'main', 'header1', 'subheader1'],
                                      ['K12345', 'description2', 'module2', 'main', 'header1', 'subheader1']],
                                     columns=['gene_id', 'gene_description', 'module', 'sheet', 'header', 'subheader']))


@pytest.fixture()
def summarized_genomes():
    return pd.DataFrame([['K00001', 'description', 'module1', 'main', 'header1', 'subheader1', 1],
                         ['K12345', 'description2', 'module2', 'main', 'header1', 'subheader1', 0]],
                        columns=['gene_id', 'gene_description', 'module', 'sheet', 'header', 'subheader', 'genome'])


@pytest.fixture()
def test_module_net():
    module_frame = pd.DataFrame([['0,0', 'M12345', 'a module name', 'K00001'],
                                 ['1,0', 'M12345', 'a module name', 'K00002'],
                                 ['2,0', 'M12345', 'a module name', 'K00003']],
                                columns=['path', 'module', 'module_name', 'ko'])
    test_module_net = build_module_net(module_frame)
    return test_module_net


@pytest.fixture()
def test_annotations_df():
    return pd.DataFrame([['', 'scaffold_1'],
                         ['K12345', 'scaffold_1'],
                         ['K00001', 'scaffold_1']],
                        index=['gene_1', 'gene_2', 'gene_3'], columns=['ko_id', 'scaffold'])

@pytest.fixture()
def test_annotations_ids_by_row_df(test_annotations_df):
    db_id_sets: pd.Series = get_annotation_ids_by_row(
        test_annotations_df
    )
    test_annotation_ids_by_row = test_annotations_df.copy()
    test_annotation_ids_by_row[DBSETS_COL] = db_id_sets
    return test_annotation_ids_by_row


@pytest.fixture()
def module_coverage_frame():
    return pd.DataFrame([['scaffold_1', 'M12345', 'a module name', 3, 1, 1 / 3, 1, 'K00001', 'gene_3']],
                        columns=['genome', 'module', 'module_name', 'steps', 'steps_present', 'step_coverage',
                                 'ko_count', 'kos_present', 'genes_present'])


@pytest.fixture()
def module_network():
    network = nx.DiGraph()
    network.add_edges_from([('start', 'K00001'), ('K00001', 'K00002'), ('K00001', 'K00003'), ('K00002', 'K00004'),
                            ('K00003', 'K00013'), ('K00013', 'K00004')])
    return network


@pytest.fixture()
def etc_module_df():
    return pd.DataFrame([['K00001+(K00002,K00003+K00013)+K00004', 'Complex I', 'oxidoreductase', 'M00000']],
                        columns=['definition', 'complex', 'module_name', 'module_id'])


@pytest.fixture()
def etc_coverage_df():
    return pd.DataFrame([['M00000', 'oxidoreductase', 'I', 'scaffold_1', 3, 1, 1/3, 'K00001',
                          'K00002,K00004', 'Complex I: oxidoreductase']],
                        columns=['module_id', 'module_name', 'complex', 'genome', 'path_length',
                                 'path_length_coverage', 'percent_coverage', 'genes', 'missing_genes',
                                 'complex_module_name'])


@pytest.fixture()
def function_heatmap_form():
    return pd.DataFrame([['Category1', 'SubCategory1', 'A function', 'K00001, K99999', 'A long function name', ''],
                         ['Category1', 'SubCategory1', 'A function', 'K00002', 'A second long function name', ''],
                         ['Category1', 'SubCategory2', 'B function', 'K12345', 'A longer function name', '']],
                        columns=['category', 'subcategory', 'function_name', 'function_ids', 'long_function_name',
                                 'gene_symbol'])


@pytest.fixture()
def functional_df():
    return pd.DataFrame([['Category1', 'SubCategory1', 'A function', 'K00001',
                          'A long function name; A second long function name', '', 'scaffold_1', False,
                          'Category1: A function'],
                         ['Category1', 'SubCategory2', 'B function', 'K12345', 'A longer function name', '',
                          'scaffold_1', True, 'Category1: B function']],
                        columns=['category', 'subcategory', 'function_name', 'function_ids', 'long_function_name',
                                 'gene_symbol', 'genome', 'present', 'category_function_name'])
