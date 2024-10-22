from __future__ import annotations

from math import pi
from pathlib import Path
from typing import Optional

import pandas as pd
import panel as pn
import param
from bokeh.core.property.vectorization import Field
from bokeh.models import ColorBar, Legend, LegendItem, LinearColorMapper, Plot
from bokeh.palettes import BuGn, Cividis256
from bokeh.plotting import figure
from bokeh.resources import INLINE
from bokeh.transform import factor_cmap, linear_cmap
from panel_jstree import Tree

from dram_viz.definitions import NO_TAXONOMY_RANKS, TAXONOMY_RANKS_REGEX

pn.extension("tabulator", "katex", template="bootstrap")

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
    return [
        heatmap(frame, **{"title": maybe_title, **kwargs}) for maybe_title, frame in df.groupby(groupby, sort=False)
    ]


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
        legend = Legend(
            items=[LegendItem(label=Field(field=label), renderers=[p.renderers[i]]) for i, label in enumerate(labels)]
        )
        p.add_layout(legend, side)
    return p_orig


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
    color_bar = ColorBar(
        color_mapper=LinearColorMapper(palette=tuple(reversed(PALETTE_CONTINUOUS)), low=0, high=1),
        height=HEATMAP_CELL_HEIGHT * 30,
    )
    p.add_layout(color_bar, "right")
    return p_orig


def heatmap(
    df,
    y_col,
    tooltip_cols,
    title="",
    rect_kw=None,
    c_min=0,
    c_max=1,
    c_col: str = None,
    x_col: str = None,
    x_cols: list[str] = None,
    extra_y_col: str = None,
    **fig_kwargs,
):
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
        **fig_kwargs,
    )

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
    p.rect(x=x_col, y=y_col, width=0.9, height=0.9, source=df, fill_alpha=0.9, color=fill_color, **rect_kw)

    p.title.text_font_size = "8pt"

    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "12px"
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
    list
        A list of heatmaps
    """

    extra_tooltip_cols = []
    if "taxonomy" in module_df.columns:
        extra_tooltip_cols.append("taxonomy")

    first_charts_kw = {}
    if taxonomy_label is not None:
        module_df["label"] = module_df[taxonomy_label] + " | " + module_df[y_col]
        # fig1_kw["extra_y_col"] = taxonomy_label
    else:
        module_df["label"] = module_df[y_col]

    completeness_charts = []
    if "Completeness" in module_df.columns:
        completeness_charts.extend(
            make_heatmap_groups(module_df, x_cols=["Contamination"], y_col="label", tooltip_cols=[y_col])
        )
        extra_tooltip_cols.append("Completeness")
        # The first chart we see (completeness, contamination, module) will have the y-axis on the left
        first_charts_kw["y_axis_location"] = None
    if "Contamination" in module_df.columns:
        completeness_charts.extend(
            make_heatmap_groups(
                module_df, x_cols=["Completeness"], y_col=y_col, tooltip_cols=[y_col], **first_charts_kw
            )
        )
        extra_tooltip_cols.append("Contamination")
        # The first chart we see (completeness, contamination, module) will have the y-axis on the left
        first_charts_kw["y_axis_location"] = None

    if taxonomy_label is not None:
        first_charts_kw["extra_y_col"] = taxonomy_label
    module_charts = make_heatmap_groups(
        module_df,
        x_col="module_name",
        y_col=y_col,
        c_col="step_coverage",
        tooltip_cols=["genome", "module_name", "steps", "steps_present", *extra_tooltip_cols],
        **first_charts_kw,
        title="Module",
    )
    # etc_charts = add_colorbar(make_heatmap_groups(etc_df, x_col="module_name", y_col=y_col, c_col="percent_coverage",
    #                                  groupby="complex",
    #                                  tooltip_cols=["genome", "module_name", "path_length", "path_length_coverage",
    #                                                "genes", "missing_genes", *extra_tooltip_cols],
    #                                               y_axis_location=None,),
    #                           index=-1)
    etc_charts = make_heatmap_groups(
        etc_df,
        x_col="module_name",
        y_col=y_col,
        c_col="percent_coverage",
        groupby="complex",
        tooltip_cols=[
            "genome",
            "module_name",
            "path_length",
            "path_length_coverage",
            "genes",
            "missing_genes",
            *extra_tooltip_cols,
        ],
        y_axis_location=None,
    )
    #
    function_charts = add_legend(
        make_heatmap_groups(
            function_df,
            x_col="function_name",
            y_col=y_col,
            c_col="present",
            groupby="category",
            tooltip_cols=[
                "genome",
                "category",
                "subcategory",
                ("Function IDs", "@function_ids"),
                "function_name",
                "long_function_name",
                "gene_symbol",
                *extra_tooltip_cols,
            ],
            y_axis_location=None,
        ),
        "present",
        side="right",
        index=-1,
    )
    return [*completeness_charts, *module_charts, *etc_charts, *function_charts]


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

    def __init__(
        self,
        module_df: pd.DataFrame,
        etc_df: pd.DataFrame,
        function_df: pd.DataFrame,
        tax_tree_data=None,
        selected_tax_tree=None,
        output_dir=None,
    ):
        super().__init__()
        self.module_df = module_df
        self.etc_df = etc_df
        self.function_df = function_df
        self.tax_tree_data = tax_tree_data
        self._output_dir = output_dir or Path.cwd()
        self.plot_view = pn.Row()
        self.download_button = pn.widgets.Button(name="Download Heatmap", button_type="primary")
        self.download_button.on_click(self.download_heatmap)

        self.redraw_button = pn.widgets.Button(name="Redraw", button_type="primary")
        self.reset_button = pn.widgets.Button(name="Reset Filters", button_type="warning")
        self.reset_button.on_click(self.reset_filters)

        self.redraw_button.on_click(self.make_plot)

        self.tax_axis_filter = pn.widgets.Checkbox(name="Show Taxonomy on Y Axis", value=False)
        self.tax_axis_rank = pn.widgets.Select(
            name="Taxonomy Label", options=list(TAXONOMY_RANKS_REGEX), visible=False, value="genus"
        )
        self.show_tax_box = pn.Column(self.tax_axis_filter, self.tax_axis_rank)
        pn.bind(self.reveal_tax_axis_rank_selector, self.tax_axis_filter, watch=True)

        if "taxonomy" in self.module_df.columns:
            self.taxonomy_filter = Tree(data=self.tax_tree_data, show_icons=False, cascade=True)

            # hack to make sure the taxonomy filter tree.value is set since it isn't set on the first load
            # TODO: remove maybe when this is put into panel
            self.taxonomy_filter.value = selected_tax_tree or []

            sort_options = ["genome", *list(TAXONOMY_RANKS_REGEX.keys())]
        else:
            self.taxonomy_filter = None
            sort_options = ["genome"]
        if "Completeness" in self.module_df.columns:
            sort_options.append("Completeness")
        if "Contamination" in self.module_df.columns:
            sort_options.append("Contamination")

        self.sort_by = pn.widgets.MultiChoice(name="Sort By", options=sort_options)

        self._init_view()
        self.download_heatmap()

    def __panel__(self, *args, **kwargs):
        return self.view

    def _init_view(self):
        additional_sidebar = []
        module_df = self.module_df.copy()
        etc_df = self.etc_df.copy()
        function_df = self.function_df.copy()

        charts = make_product_heatmap(
            module_df,
            etc_df,
            function_df,
            taxonomy_label=None if not self.tax_axis_filter.value else self.tax_axis_rank.value,
        )

        self.plot_view[:] = charts

        if "taxonomy" in self.module_df.columns:
            additional_sidebar.append(self.show_tax_box)
            additional_sidebar.append("## Taxonomy Filter")
            additional_sidebar.append(self.taxonomy_filter)

        self.view = pn.template.FastListTemplate(
            title="DRAM Product Visualization",
            # main=[self.plot_view],
            main=[
                pn.Tabs(
                    ("Heatmap", self.plot_view),
                    ("Module Coverage DF", pn.widgets.Tabulator(self.module_df, page_size=50)),
                    ("ETC Coverage DF", pn.widgets.Tabulator(self.etc_df, page_size=50)),
                    ("Function DF", pn.widgets.Tabulator(self.function_df, page_size=50)),
                )
            ],
            sidebar=[
                pn.Row(self.redraw_button, self.reset_button),
                self.download_button,
                self.sort_by,
                self.param.min_coverage,
                *additional_sidebar,
            ],
        )

    def make_plot(self, event=None):
        """
        Make the product plot
        """
        additional_sidebar = []
        module_df = self.module_df.copy()
        etc_df = self.etc_df.copy()
        function_df = self.function_df.copy()

        if self.min_coverage > 0:
            print(self.min_coverage)
            # module_df = self.module_df[self.module_df["step_coverage"] >= self.min_coverage]
            module_df.loc[self.module_df["step_coverage"] < self.min_coverage, "step_coverage"] = 0
            # etc_df = self.etc_df[self.etc_df["percent_coverage"] >= self.min_coverage]
            etc_df.loc[self.etc_df["percent_coverage"] < self.min_coverage, "percent_coverage"] = 0

        module_df, etc_df, function_df = self.filter_by_taxonomy(module_df, etc_df, function_df)
        module_df, etc_df, function_df = self.get_sorted_dfs(module_df, etc_df, function_df, by=self.sort_by.value)

        charts = make_product_heatmap(
            module_df,
            etc_df,
            function_df,
            taxonomy_label=None if not self.tax_axis_filter.value else self.tax_axis_rank.value,
        )

        self.plot_view[:] = charts

    def reset_filters(self, event=None):
        """
        Resets the filters applied to the heatmap.

        Parameters:
        - event (optional): The event that triggered the reset. Defaults to None.
        """

        self.min_coverage = self.param.min_coverage.default

        if self.taxonomy_filter is not None:
            self.taxonomy_filter.value = [node["id"] for node in self.taxonomy_filter.flat_tree]
            self.tax_axis_filter.value = False
            self.tax_axis_rank.visible = False
            self.tax_axis_rank.value = "genus"

        self.sort_by.value = []

    def filter_by_taxonomy(self, module_df, etc_df, function_df):
        """
        Filter the dataframes by selected taxonomy ranks
        """
        if self.taxonomy_filter is None:
            return module_df, etc_df, function_df
        selected = self.taxonomy_filter.value
        # leaves = [node for node in selected if len(node.split(";")) == NO_TAXONOMY_RANKS]
        # maybe we don't need this replace, but leaving in for now to be sure we match the data
        leaves = [node.replace("; ", ";") for node in selected if len(node.split(";")) == NO_TAXONOMY_RANKS]
        module_df = module_df.loc[module_df["taxonomy"].isin(leaves)]
        etc_df = etc_df.loc[etc_df["taxonomy"].isin(leaves)]
        function_df = function_df.loc[function_df["taxonomy"].isin(leaves)]

        return module_df, etc_df, function_df

    def reveal_tax_axis_rank_selector(self, event=None, tax_axis_filter_value: bool = None):
        """
        Set the taxonomy filter
        """
        if isinstance(tax_axis_filter_value, bool) and tax_axis_filter_value:
            self.tax_axis_filter.value = tax_axis_filter_value
        if self.tax_axis_filter.value:
            self.tax_axis_rank.visible = True
            return
        self.tax_axis_rank.visible = False

    def get_sorted_dfs(self, module_df, etc_df, function_df, by="genome"):
        """
        Sort the dataframes by taxonomy
        """
        return module_df.sort_values(by=by), etc_df.sort_values(by=by), function_df.sort_values(by=by)

    def download_heatmap(self, event=None, output_dir=None):
        """
        Save the heatmap to a file
        """
        output_dir = output_dir or self._output_dir
        self.plot_view.save(output_dir / "product.html", resources=INLINE)
