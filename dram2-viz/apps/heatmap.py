from __future__ import annotations

from math import pi
from typing import Optional

import pandas as pd
import panel as pn
from bokeh.core.property.vectorization import Field
from bokeh.models import Plot, Legend, LegendItem, ColorBar, LinearColorMapper
from bokeh.palettes import BuGn, Cividis256
from bokeh.plotting import figure
from bokeh.transform import linear_cmap, factor_cmap
import param
from panel_jstree import Tree


try:
    from ..definitions import DEFAULT_GROUPBY_COLUMN, BACKUP_GROUPBY_COLUMN, HEATMAP_MODULES, TAXONOMY_RANKS_REGEX, NO_TAXONOMY_RANKS
except ImportError:
    import os
    import sys

    sys.path.insert(0, os.path.abspath(".."))
    from viz.definitions import DEFAULT_GROUPBY_COLUMN, BACKUP_GROUPBY_COLUMN, HEATMAP_MODULES, TAXONOMY_RANKS_REGEX, NO_TAXONOMY_RANKS

pn.extension('tabulator', 'katex', template='bootstrap')

PALETTE_CATEGORICAL = BuGn
PALETTE_CONTINUOUS = Cividis256

HEATMAP_CELL_HEIGHT = 15
HEATMAP_CELL_WIDTH = 15


def make_heatmap_groups(df: pd.DataFrame, groupby: Optional[str] = None, title: Optional[str] = None, **kwargs):
    """
    Generate a list of heatmaps based on the given DataFrame and grouping.

    Parameters:
    - df (pd.DataFrame): The DataFrame containing the data for the heatmaps.
    - groupby (Optional[str]): The column name to group the data by. If not provided, a single heatmap will be generated for the entire DataFrame.
    - title (Optional[str]): The title for the heatmaps. If provided, it will be used as the title for all heatmaps. If not provided, the groups will be used as titles.
    - **kwargs: Additional keyword arguments to be passed to the heatmap function.

    Returns:
    - List[plt.Figure]: A list of heatmaps, each represented as a matplotlib Figure object.
    """
    if title:
        kwargs["title"] = title
    if not groupby:
        return [heatmap(df, **kwargs)]
    # if not title, use groups as titles
    return [heatmap(frame, **{"title": maybe_title, **kwargs}) for maybe_title, frame in df.groupby(groupby)]


def add_legend(p_orig: Plot | list[Plot], labels: str | list[str], side="right", index: Optional[int] = None):
    """
    Add a legend to a Bokeh plot.

    Parameters:
    - p_orig (Plot | list[Plot]): The Bokeh plot or a list of Bokeh plots to add the legend to.
    - labels (str | list[str]): The label(s) to display in the legend. Can be a single string or a list of strings.
    - side (str): The side of the plot where the legend should be placed. Default is "right".
    - index (Optional[int]): The index of the plot in the list of plots (if p_orig is a list). Default is None.

    Returns:
    - p_orig: The original Bokeh plot(s) with the legend added.
    """

    if isinstance(p_orig, list) and index is None:
        raise ValueError("If p is a list, i must be an integer")
    if isinstance(p_orig, list):
        if len(p_orig) == 0:
            return p_orig
        p = p_orig[index]
    else:
        p = p_orig
    if isinstance(labels, str):
        label = labels
        legend = Legend(items=[LegendItem(label=Field(field=label), renderers=p.renderers)])
        p.add_layout(legend, side)
    else:
        legend = Legend(items=[LegendItem(label=Field(field=label), renderers=[p.renderers[i]]) for i, label in enumerate(labels)])
        p.add_layout(legend, side)
    return p_orig


def format_chart_group(chart_group: list, title: str = ""):
    """
    Make a formatted chart group

    Parameters
    ----------
    chart_group : list
        A list of charts
    title : str
        The title of the chart group

    Returns
    -------
    pn.Column
        A panel column of charts
    """
    if not isinstance(chart_group, list) or isinstance(chart_group, tuple):
        chart_group = [chart_group]
    return pn.Column(
        pn.pane.Markdown(f"## {title}", align="center"),
        pn.Row(*chart_group)
    )


def add_colorbar(p_orig: Plot | list[Plot], index: Optional[int] = None):
    """
    Add a colorbar to a plot
    """
    if isinstance(p_orig, list) and index is None:
        raise ValueError("If p is a list, i must be an integer")
    if isinstance(p_orig, list):
        if len(p_orig) == 0:
            return p_orig
        p = p_orig[index]
    else:
        p = p_orig
    color_bar = ColorBar(color_mapper=LinearColorMapper(palette=tuple(reversed(PALETTE_CONTINUOUS)), low=0, high=1),
                         height=HEATMAP_CELL_HEIGHT * 30)
    p.add_layout(color_bar, 'right')
    return p_orig


def heatmap(df, y_col, tooltip_cols, title="", rect_kw=None, c_min=0, c_max=1, c_col: str = None, x_col: str = None, x_cols: list[str] = None, extra_y_col: str = None,  **fig_kwargs, ):
    """
    Make a heatmap from a dataframe

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe to make the heatmap from
    x_col : str
        The column to use for the x-axis
    y_col : str
        The column to use for the y-axis
    c_col : str
        The column to use for the color
    tooltip_cols : list
        A list of tuples of columns to use for the tooltips
    title : str
        The title of the plot
    rect_kw : dict
        Keyword arguments for the rect
    c_min : float
        The minimum value for the color
    c_max : float
        The maximum value for the color
    extra_y_col: str
        An extra column to use for the y-axis

    Returns
    -------
    Plot
        The heatmap plot
    """

    if x_cols:
        df = pd.melt(df, id_vars=y_col, value_vars=x_cols, var_name="x_col").drop_duplicates()
        x_col = "x_col"
        c_col = "value"
        tooltip_cols = [y_col, "value"]

    rect_kw = rect_kw or {}
    # df = df.sort_values(by=[y_col], ascending=False)

    tooltips = []
    for col in tooltip_cols:
        if isinstance(col, tuple):
            tooltips.append(col)
        else:
            tooltips.append((col.replace("_", " ").title(), f"@{col}"))

    # half_ttl_ln = len(title) / 2
    # if half_ttl_ln > len(
    #         df[x_col].unique()) and "min_border_left" not in fig_kwargs and "min_border_right" not in fig_kwargs:
    #     left_over = (half_ttl_ln - len(df[x_col].unique()))
        # fig_kwargs["min_border_left"] = int(left_over / 1.2) * HEATMAP_CELL_WIDTH
        # fig_kwargs["min_border_right"] = int(left_over / 1.2) * HEATMAP_CELL_WIDTH

    p = figure(
        frame_width=HEATMAP_CELL_WIDTH * len(df[x_col].unique()),
        frame_height=HEATMAP_CELL_WIDTH * len(df[y_col].unique()),

        x_range=sorted(list(df[x_col].unique())),
        y_range=list(df[y_col].unique()),
        tools="hover",
        toolbar_location=None,
        tooltips=tooltips,
        # title="\n".join(char for char in title),
        title=title,
        # title_location="right",
        **fig_kwargs
    )
    #
    # if extra_y_col:
    #     # p.extra_y_ranges[extra_y_col] = FactorRange(factors=list(df[extra_y_col].unique()))
    #     p.extra_y_ranges[extra_y_col] = Range1d(-0.1, len(df[extra_y_col].unique()))
    #
    #     ax2 = LinearAxis(y_range_name=extra_y_col, axis_label=extra_y_col)
    #     # ax2.axis_label_text_color = "navy"
    #     p.add_layout(ax2, 'left')


    # if x_col:
    if df[c_col].dtype == float:
        palette = tuple(reversed(PALETTE_CONTINUOUS))
        fill_color = linear_cmap(c_col, palette=palette, low=c_min, high=c_max)
    else:
        df[c_col] = df[c_col].astype(str)
        factors = sorted(df[c_col].unique())
        max_factors = max(PALETTE_CATEGORICAL.keys())
        palette = PALETTE_CATEGORICAL[max(len(factors), 3)] if len(factors) <= max_factors else PALETTE_CONTINUOUS
        fill_color = factor_cmap(c_col, palette=tuple(reversed(palette)), factors=factors)
    p.rect(x=x_col, y=y_col,
           width=0.9, height=0.9,
           source=df,
           fill_alpha=0.9,
           color=fill_color,
           **rect_kw
           )
    # else:
    #     for x_col in x_cols:
    #         palette = tuple(reversed(PALETTE_CONTINUOUS))
    #         fill_color = linear_cmap(x_col, palette=palette, low=c_min, high=c_max)
    #
    #         p.rect(x=x_col, y=y_col,
    #                width=0.9, height=0.9,
    #                source=df,
    #                fill_alpha=0.9,
    #                color=fill_color,
    #                **rect_kw
    #                )

    p.title.align = "left"

    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "14px"
    p.xaxis.major_label_orientation = pi / 2

    return p


def make_product_heatmap(
        module_df: pd.DataFrame,
        etc_df: pd.DataFrame,
        function_df: pd.DataFrame,
        y_col: str = "genome",
        taxonomy_label: pd.Series | None = None,
):
    """
    Make a product heatmap group from the module_coverage_df, etc_coverage_df, and functional_df

    Parameters
    ----------
    module_df : pd.DataFrame
        A dataframe of module coverage information
    etc_df : pd.DataFrame
        A dataframe of ETC coverage information
    function_df : pd.DataFrame
        A dataframe of functional coverage information
    y_col : str
        The column to use for the y-axis in the heatmaps, must be present in all dataframes
        default: "genome"
    taxonomy_label : pd.Series
        The current taxonomy label to use for the y-axis in the heatmaps

    Returns
    -------
    pn.Row
        A panel row of heatmaps
    """

    extra_tooltip_cols = []
    if "taxonomy" in module_df.columns:
        extra_tooltip_cols.append("taxonomy")
    if "Completeness" in module_df.columns:
        extra_tooltip_cols.append("Completeness")
    if "Contamination" in module_df.columns:
        extra_tooltip_cols.append("Contamination")

    fig1_kw = {}
    if taxonomy_label is not None:
        module_df["label"] = module_df[taxonomy_label] + " | " + module_df[y_col]
        # fig1_kw["extra_y_col"] = taxonomy_label
    else:
        module_df["label"] = module_df[y_col]
    completeness_charts = [*make_heatmap_groups(module_df, x_cols=["Contamination"], y_col="label",
                                        tooltip_cols=[y_col],
                                                **fig1_kw
                                                ),
                           *make_heatmap_groups(module_df, x_cols=["Completeness"], y_col=y_col,
                                               tooltip_cols=[y_col],
                                               y_axis_location=None)
                           ]
    if taxonomy_label is not None:
        fig1_kw["extra_y_col"] = taxonomy_label
    module_charts = make_heatmap_groups(module_df, x_col="module_name", y_col=y_col, c_col="step_coverage",
                                        tooltip_cols=["genome", "module_name", "steps", "steps_present", *extra_tooltip_cols],
                                        y_axis_location=None,
                                        **fig1_kw,
                                        title="Module")
    # etc_charts = add_colorbar(make_heatmap_groups(etc_df, x_col="module_name", y_col=y_col, c_col="percent_coverage",
    #                                  groupby="complex",
    #                                  tooltip_cols=["genome", "module_name", "path_length", "path_length_coverage",
    #                                                "genes", "missing_genes", *extra_tooltip_cols],
    #                                               y_axis_location=None,),
    #                           index=-1)
    etc_charts = make_heatmap_groups(etc_df, x_col="module_name", y_col=y_col, c_col="percent_coverage",
                                     groupby="complex",
                                     tooltip_cols=["genome", "module_name", "path_length", "path_length_coverage",
                                                   "genes", "missing_genes", *extra_tooltip_cols],
                                                  y_axis_location=None,)
    #
    function_charts = add_legend(make_heatmap_groups(function_df, x_col="function_name", y_col=y_col, c_col="present",
                                          groupby="category",
                                          tooltip_cols=["genome", "category", "subcategory",
                                                        ("Function IDs", "@function_ids"),
                                                        "function_name", "long_function_name", "gene_symbol",
                                                        *extra_tooltip_cols],
                                          y_axis_location=None,),
                                "present", side="right", index=-1)
    return [*completeness_charts, *module_charts, *etc_charts, *function_charts]

    charts = [
        format_chart_group([p for p in completeness_charts]),
        format_chart_group([p for p in module_charts]),
        format_chart_group([p for p in etc_charts], title="ETC Complexes"),
        format_chart_group([p for p in function_charts]),
    ]
    return charts

    # plot = pn.Row(*charts)
    # return plot



class Dashboard(pn.viewable.Viewer):
    """
    A class representing a dashboard for visualizing data.

    Parameters:
    -----------
    module_df : pd.DataFrame
        The module dataframe.
    etc_df : pd.DataFrame
        The etc dataframe.
    function_df : pd.DataFrame
        The function dataframe.
    tax_tree_data : optional
        The taxonomy tree data.
    """


    min_coverage = param.Number(default=0, bounds=(0, 1), label="Minimum Coverage")
    view = param.ClassSelector(class_=pn.template.FastListTemplate)
    plot_view = param.ClassSelector(class_=pn.Row)

    def __init__(self, module_df: pd.DataFrame, etc_df: pd.DataFrame, function_df: pd.DataFrame, tax_tree_data=None):
        super().__init__()
        self.module_df = module_df
        self.etc_df = etc_df
        self.function_df = function_df
        self.tax_tree_data = tax_tree_data
        self.plot_view = pn.Row()

        self.redraw_button = pn.widgets.Button(name="Redraw", button_type="primary")
        self.reset_button = pn.widgets.Button(name="Reset Filters", button_type="primary")
        self.reset_button.on_click(self.reset_filters)

        self.redraw_button.on_click(self.make_plot)

        self.tax_axis_filter = pn.widgets.Checkbox(name="Show Taxonomy", value=False)
        self.tax_axis_rank = pn.widgets.Select(name='Taxonomy Label', options=list(TAXONOMY_RANKS_REGEX), visible=False, value="genus")
        self.show_tax_box = pn.Column(
            self.tax_axis_filter,
            self.tax_axis_rank
        )
        pn.bind(self.set_taxonomy_axis_filter, self.tax_axis_filter, watch=True)

        if "taxonomy" in self.module_df.columns:
            self.taxonomy_filter = Tree(data=self.tax_tree_data, show_icons=False, cascade=True)
            sort_options = ["genome", *list(TAXONOMY_RANKS_REGEX.keys())]
        else:
            self.taxonomy_filter = None
            sort_options = ["genome"]
        self._taxonomy_filter_initiated = False
        if "Completeness" in self.module_df.columns:
            sort_options.append("Completeness")
        if "Contamination" in self.module_df.columns:
            sort_options.append("Contamination")

        self.sort_by = pn.widgets.MultiChoice(name="Sort By", options=sort_options)

        self.make_plot()

    def make_plot(self, event=None):
        """
        Make the product plot
        """
        additional_sidebar = []
        module_df = self.module_df.copy()
        etc_df = self.etc_df.copy()
        function_df = self.function_df.copy()

        if event:  # don't filter on initial load (These are filters after browser load)
            # This is a hack to make sure the taxonomy filter tree.value is set.
            # It is not set on the first load, even though we pass in the data as selected,
            # Those aren't currently back propagated to the python side
            if not self._taxonomy_filter_initiated and "taxonomy" in self.module_df.columns:
                self.reset_taxonomy()
                self._taxonomy_filter_initiated = True

            if self.min_coverage > 0:
                print(self.min_coverage)
                # module_df = self.module_df[self.module_df["step_coverage"] >= self.min_coverage]
                module_df.loc[self.module_df["step_coverage"] < self.min_coverage, "step_coverage"] = 0
                # etc_df = self.etc_df[self.etc_df["percent_coverage"] >= self.min_coverage]
                etc_df.loc[self.etc_df["percent_coverage"] < self.min_coverage, "percent_coverage"] = 0
            print("making plot")

            module_df, etc_df, function_df = self.filter_by_taxonomy(module_df, etc_df, function_df)
            module_df, etc_df, function_df = self.get_sorted_dfs(module_df, etc_df, function_df, by=self.sort_by.value)

        charts = make_product_heatmap(module_df, etc_df, function_df, taxonomy_label=None if not self.tax_axis_filter.value else self.tax_axis_rank.value)

        self.plot_view[:] = charts

        if "taxonomy" in self.module_df.columns:
            additional_sidebar.append("## Taxonomy Filter")
            additional_sidebar.append(self.taxonomy_filter)

        self.view = pn.template.FastListTemplate(
            title="DRAM2 Product Visualization",
            # main=[self.plot_view],
            main=[
                pn.Tabs(
                    ("Heatmap", self.plot_view),
                    ("Module Coverage DF", pn.widgets.Tabulator(module_df, page_size=50)),
                    ("ETC Coverage DF", pn.widgets.Tabulator(etc_df, page_size=50)),
                    ("Function DF", pn.widgets.Tabulator(self.function_df, page_size=50)),
                )
            ],
            sidebar=[
                pn.Row(self.redraw_button, self.reset_button),
                self.show_tax_box,
                self.sort_by,
                self.param.min_coverage,
                *additional_sidebar,
            ]
        )

    def reset_filters(self, event=None):
            """
            Resets the filters applied to the heatmap.

            Parameters:
            - event (optional): The event that triggered the reset. Defaults to None.
            """

            self.min_coverage = self.param.min_coverage.default

            if self.taxonomy_filter is not None:
                self.reset_taxonomy()

            self.sort_by.value = []

    def reset_taxonomy(self):
        """
        Resets the taxonomy filter to its initial state.

        This method sets the value of the `taxonomy_filter` attribute to the list of IDs of all nodes in the flat tree.
        """
        self.taxonomy_filter.value = [node["id"] for node in self.taxonomy_filter.flat_tree]

    def filter_by_taxonomy(self, module_df, etc_df, function_df):
        """
        Filter the dataframes by selected taxonomy ranks
        """
        if self.taxonomy_filter is None:
            return module_df, etc_df, function_df
        selected = self.taxonomy_filter.value
        leaves = [node for node in selected if len(node.split(";")) == NO_TAXONOMY_RANKS]
        module_df = module_df.loc[module_df["taxonomy"].isin(leaves)]
        etc_df = etc_df.loc[etc_df["taxonomy"].isin(leaves)]
        function_df = function_df.loc[function_df["taxonomy"].isin(leaves)]

        return module_df, etc_df, function_df

    def set_taxonomy_axis_filter(self, event=None):
        """
        Set the taxonomy filter
        """
        if self.tax_axis_filter.value:
            self.tax_axis_rank.visible = True
            return
        self.tax_axis_rank.visible = False

    def get_sorted_dfs(self, module_df, etc_df, function_df, by="genome"):
        """
        Sort the dataframes by taxonomy
        """
        return module_df.sort_values(by=by), etc_df.sort_values(by=by), function_df.sort_values(by=by)
