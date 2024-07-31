import panel as pn
from bokeh.models import Plot

from dram2_viz.apps.heatmap import Dashboard, make_product_heatmap


def build_tax_tree_selected_recurse(tax_tree, flat_tree=None):
    if flat_tree is None:
        flat_tree = []
    for node in tax_tree:
        flat_tree.append(node["id"])
        if "children" in node:
            build_tax_tree_selected_recurse(node["children"], flat_tree)
    return flat_tree


def test_make_product_heatmap(module_coverage_frame, etc_coverage_df, functional_df):
    plots = make_product_heatmap(module_coverage_frame, etc_coverage_df, functional_df)
    assert isinstance(plots, list)
    for p in plots:
        assert isinstance(p, Plot)


def test_heatmap_dashboard_instantiation(module_coverage_frame, etc_coverage_df, functional_df):
    app = Dashboard(module_coverage_frame, etc_coverage_df, functional_df)
    assert isinstance(app, Dashboard)
    assert isinstance(app.view, pn.template.BaseTemplate)


def test_heatmap_dashboard_no_filter_taxonomy_returns_all(
    module_coverage_df_from_file, etc_coverage_df_from_file, function_df_from_file, taxonomy_tree
):
    module_df_length = len(module_coverage_df_from_file)
    etc_df_length = len(etc_coverage_df_from_file)
    function_df_length = len(function_df_from_file)

    app = Dashboard(module_coverage_df_from_file, etc_coverage_df_from_file, function_df_from_file, taxonomy_tree)

    # We have to do this manually because the taxonomy filter is a panel object and requires a browser to
    # instantiate fully
    app.taxonomy_filter.value = build_tax_tree_selected_recurse(taxonomy_tree)

    m2, e2, f2 = app.filter_by_taxonomy(module_coverage_df_from_file, etc_coverage_df_from_file, function_df_from_file)

    # check that the dataframes are the same length as the input dataframes
    assert len(m2) == module_df_length
    assert len(e2) == etc_df_length
    assert len(f2) == function_df_length
