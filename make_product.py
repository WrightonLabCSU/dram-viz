"""
===================
DRAM Visualizations
===================

Script that generates a product visualization from the DRAM output.
"""
from __future__ import annotations
import argparse
import logging
from pathlib import Path
from typing import Optional

from bokeh.resources import INLINE
import pandas as pd
import panel as pn
import param
from panel_jstree import Tree

try:
    from .definitions import DEFAULT_GROUPBY_COLUMN, BACKUP_GROUPBY_COLUMN
    from .proccessing.process_annotations import (DBSETS_COL, MODULE_STEPS_FORM_TAG, FUNCTION_HEATMAP_FORM_TAG,
        ETC_MODULE_DF_TAG, FILES_NAMES, build_module_net, fill_product_dfs, make_product_df,
        get_phylum_and_most_specific,
        make_strings_no_repeats, get_annotation_ids_by_row, rename_genomes_to_taxa)
    from .apps.heatmap import make_product_heatmap
except ImportError:
    from definitions import DEFAULT_GROUPBY_COLUMN, BACKUP_GROUPBY_COLUMN
    from processing.process_annotations import (DBSETS_COL, MODULE_STEPS_FORM_TAG, FUNCTION_HEATMAP_FORM_TAG,
        ETC_MODULE_DF_TAG, FILES_NAMES, build_module_net, fill_product_dfs, make_product_df,
        get_phylum_and_most_specific,
        make_strings_no_repeats, get_annotation_ids_by_row, rename_genomes_to_taxa)
    from apps.heatmap import make_product_heatmap

__authors__ = ["Madeline Scyphers", "Rory Flynn"]
__copyright__ = "Copyright 2024, Wrighton Lab"

logger = logging.getLogger("dram2_log.viz")

pn.extension('tabulator', 'katex', template='bootstrap')

HEATMAP_MODULES = [
    "M00001",
    "M00004",
    "M00008",
    "M00009",
    "M00012",
    "M00165",
    "M00173",
    "M00374",
    "M00375",
    "M00376",
    "M00377",
    "M00422",
    "M00567",
]

# regex portion not used right now, but could be useful in the future
# captures in named pandas columns that remove the prefix string (e.g. "d__")
TAXONOMY_RANKS_REGEX = {
    "domain": r"(?:;?d__)(?P<domain>.*?)",
    "phylum": r"(?:;?p__)(?P<phylum>.*?)",
    "class": r"(?:;?c__)(?P<class>.*?)",
    "order": r"(?:;?o__)(?P<order>.*?)",
    "family": r"(?:;?f__)(?P<family>.*?)",
    "genus": r"(?:;?g__)(?P<genus>.*?)",
    "species": r"(?:;?s__)(?P<species>.*)"
}

NO_TAXONOMY_RANKS = len(TAXONOMY_RANKS_REGEX)


class Dashboard(pn.viewable.Viewer):
    """
    A parameterized class for the product visualization
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
        Reset the taxonomy filter
        """
        self.min_coverage = self.param.min_coverage.default

        if self.taxonomy_filter is not None:
            self.reset_taxonomy()

        self.sort_by.value = []

    def reset_taxonomy(self):
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


def build_tree(edge_df, source_col: str = "source", target_col: str = "target", state: dict = None, id_cb=None):
    def recurse_tree(source, id_cb, parent_id=None):
        target_nodes = edge_df.loc[edge_df[source_col] == source, target_col]
        # if id_cb is not None:
        js_ls = []
        for target in target_nodes:
            id_ = id_cb(source, target, parent_id)

            # children = recurse_tree(target, parent_id=id_, id_cb=id_cb)
            # state = {"opened": False, "disable": False, "selected": False} if
            js_ls.append({"text": target, "children": recurse_tree(target, parent_id=id_, id_cb=id_cb), "state": state, "id": id_})
            # js_ls.append({"text": target, "children": recurse_tree(target, id_cb=id_cb), "state": state})
        return js_ls
            # return [{"text": child, "children": recurse_tree(child, id_cb=id_cb), "state": state, "id": id_cb(source, child)} for child in children]
        # return [{"text": child, "children": recurse_tree(child, id_cb=id_cb), "state": state} for child in children]
    state = state or {}

    roots = edge_df.loc[~edge_df[source_col].isin(edge_df[target_col]), source_col].unique()
    tree_data = [{"id": root, "text": root, "children": recurse_tree(root, parent_id=root, id_cb=id_cb), "state": state} for root in roots]
    # tree_data = [{"id": root, "text": root, "children": recurse_tree(root, parent_id=root, id_cb=id_cb), "state": {"opened": False, "disable": False, "selected": False}} for root in roots]
    # tree_data = [{"text": root, "children": recurse_tree(root, parent_id=root, id_cb=id_cb), "state": state} for root in roots]
    return tree_data


def main(annotations_tsv_path,
         groupby_column=DEFAULT_GROUPBY_COLUMN,
         output_dir=None,
         module_steps_form: Optional[Path] = None,
         etc_steps_form: Optional[Path] = None,
         function_steps_form: Optional[Path] = None,
         dashboard=False
         ):
    """
    Make a product visualization

    Parameters
    ----------
    annotations_tsv_path : str
        The path to the annotations tsv file
    groupby_column : str
        The column to group by
    output_dir : str
        The output directory
    module_steps_form : str
        The module step database tsv
    etc_steps_form : str
        The etc step database tsv
    function_steps_form : str
        The function step database tsv
    show : bool
        Launch as dashboard
    """

    output_dir = output_dir or Path.cwd().resolve()
    annotations = pd.read_csv(annotations_tsv_path, sep="\t", index_col=0)

    db_id_sets: pd.Series = get_annotation_ids_by_row(
        annotations
    )
    annotation_ids_by_row = annotations.copy()
    annotation_ids_by_row[DBSETS_COL] = db_id_sets

    module_steps_form = pd.read_csv(module_steps_form or FILES_NAMES[MODULE_STEPS_FORM_TAG], sep="\t")
    etc_module_df = pd.read_csv(etc_steps_form or FILES_NAMES[ETC_MODULE_DF_TAG], sep="\t")
    function_heatmap_form = pd.read_csv(function_steps_form or FILES_NAMES[FUNCTION_HEATMAP_FORM_TAG], sep="\t")

    if groupby_column not in annotations.columns:
        if BACKUP_GROUPBY_COLUMN in annotations.columns:
            groupby_column = BACKUP_GROUPBY_COLUMN
        elif BACKUP_GROUPBY_COLUMN in annotations.columns:
            groupby_column = BACKUP_GROUPBY_COLUMN
        else:
            raise ValueError(f"Groupby column {groupby_column} not found in annotations")

    # make product
    if "bin_taxonomy" in annotations:
        # if gtdb format then get phylum and most specific
        if all(
                [
                    i[:3] == "d__" and len(i.split(";")) == 7
                    for i in annotations["bin_taxonomy"].fillna("")
                ]
        ):
            taxa_str_parser = get_phylum_and_most_specific
        # else just throw in what is there
        else:

            def taxa_str_parser(x):
                return x

        labels = make_strings_no_repeats(
            {
                row[groupby_column]: taxa_str_parser(row["bin_taxonomy"])
                for _, row in annotations.iterrows()
            }
        )
    else:
        labels = None

    # make module coverage frame
    module_nets = {
        module: build_module_net(module_df)
        for module, module_df in module_steps_form.groupby("module")

        if module in HEATMAP_MODULES
    }

    # module_coverage_df = pd.read_csv(output_dir / "module_coverage.tsv", sep="\t")
    # etc_coverage_df = pd.read_csv(output_dir / "etc_coverage.tsv", sep="\t")
    # function_df = pd.read_csv(output_dir / "function_coverage.tsv", sep="\t")

    module_coverage_df, etc_coverage_df, function_df = fill_product_dfs(
        annotations=annotations,
        module_nets=module_nets,
        etc_module_df=etc_module_df,
        function_heatmap_form=function_heatmap_form,
        annotation_ids_by_row=annotation_ids_by_row,
        groupby_column=groupby_column
    )

    tax_tree_data = None
    if "taxonomy" in annotations:
        cols = [ groupby_column, "taxonomy"]
        if "Completeness" in annotations.columns:
            cols.append("Completeness")
        if "Contamination" in annotations.columns:
            cols.append("Contamination")
        tax_df = annotations[cols].drop_duplicates()
        tax_df.rename(columns={groupby_column: "genome"}, inplace=True)

        # regex = "".join([regex for regex in TAXONOMY_RANKS_REGEX.values()])  # regex that might be useful later
        # tree = tax_df["taxonomy"].str.extractall(regex)
        tree = pd.DataFrame(tax_df["taxonomy"].str.split(";").to_list(), columns=["domain", "phylum", "class", "order", "family", "genus", "species"], index=tax_df.index)
        tax_df = tax_df.merge(tree, left_index=True, right_index=True, validate="1:1")
        tax_df["taxonomy"] = tax_df["domain"] + "; " + tax_df["phylum"] + "; " + tax_df["class"] + "; " + tax_df["order"] + "; " + tax_df["family"] + "; " + tax_df["genus"] + "; " + tax_df["species"]

        tax_edge_df = pd.concat(
            [tree[["domain", "phylum"]].rename(columns={"domain": "source", "phylum": "target"}),
             tree[["phylum", "class"]].rename(columns={"phylum": "source", "class": "target"}),
             tree[["class", "order"]].rename(columns={"class": "source", "order": "target"}),
             tree[["order", "family"]].rename(columns={"order": "source", "family": "target"}),
             tree[["family", "genus"]].rename(columns={"family": "source", "genus": "target"}),
             tree[["genus", "species"]].rename(columns={"genus": "source", "species": "target"})]
        ).drop_duplicates().reset_index(drop=True)

        tax_tree_data = build_tree(tax_edge_df, state={"opened": False, "selected": True}, id_cb=lambda source, child, parent_id: f"{parent_id}; {child}")

        module_coverage_df = tax_df.merge(module_coverage_df, on="genome", how="left")
        etc_coverage_df = tax_df.merge(etc_coverage_df, on="genome", how="left")
        function_df = tax_df.merge(function_df, on="genome", how="left")

    product_df = make_product_df(module_coverage_df, etc_coverage_df, function_df)

    if labels is not None:
        function_df = rename_genomes_to_taxa(function_df, labels)

    app = Dashboard(module_coverage_df, etc_coverage_df, function_df, tax_tree_data=tax_tree_data)

    app.plot_view.save(output_dir / "product.html", resources=INLINE)
    if dashboard:
        pn.serve(app.view, port=5006)
        # plot.show(port=5007)
    product_df.to_csv(output_dir / "product.tsv", sep="\t", index=False)
    logger.info("Completed visualization")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a product visualization from the DRAM output.")
    parser.add_argument("--annotations", help="Path to the annotations tsv file")
    parser.add_argument("--groupby_column", help="Column to group by", default=DEFAULT_GROUPBY_COLUMN)
    parser.add_argument("--output_dir", help="Output directory", default=Path.cwd().resolve())
    parser.add_argument("--module_steps_form", help="Module Step Database TSV",
                        default=FILES_NAMES[MODULE_STEPS_FORM_TAG])
    parser.add_argument("--etc_steps_form", help="ETC Step Database TSV", default=FILES_NAMES[ETC_MODULE_DF_TAG])
    parser.add_argument("--function_steps_form", help="Function Step Database TSV",
                        default=FILES_NAMES[FUNCTION_HEATMAP_FORM_TAG])
    parser.add_argument("--dashboard", help="Launch as dashboard", action='store_true')
    args = parser.parse_args()

    main(annotations_tsv_path=args.annotations, groupby_column=args.groupby_column, output_dir=args.output_dir,
         module_steps_form=args.module_steps_form, etc_steps_form=args.etc_steps_form,
         function_steps_form=args.function_steps_form, dashboard=args.dashboard)
