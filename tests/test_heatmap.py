import panel as pn
from bokeh.models import Plot
from dram2_viz.apps.heatmap import make_product_heatmap, Dashboard


def test_make_product_heatmap(module_coverage_frame, etc_coverage_df, functional_df):
    plots = make_product_heatmap(module_coverage_frame, etc_coverage_df, functional_df)
    assert isinstance(plots, list)
    for p in plots:
        assert isinstance(p, Plot)


def test_heatmap_dashboard_instantiation(module_coverage_frame, etc_coverage_df, functional_df):
    app = Dashboard(module_coverage_frame, etc_coverage_df, functional_df)
    assert isinstance(app, Dashboard)
    assert isinstance(app.view, pn.template.BaseTemplate)
