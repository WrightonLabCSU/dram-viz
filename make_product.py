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

import pandas as pd
import panel as pn
import param


try:
    from .definitions import DEFAULT_GROUPBY_COLUMN, BACKUP_GROUPBY_COLUMN
    from .proccessing.process_annotations import (DBSETS_COL, MODULE_STEPS_FORM_TAG, FUNCTION_HEATMAP_FORM_TAG,
        ETC_MODULE_DF_TAG, FILES_NAMES, build_module_net, fill_product_dfs, make_product_df,
        get_phylum_and_most_specific,
        make_strings_no_repeats, get_annotation_ids_by_row, rename_genomes_to_taxa)
    from .models.heatmap import make_product_heatmap
except ImportError:
    from definitions import DEFAULT_GROUPBY_COLUMN, BACKUP_GROUPBY_COLUMN
    from processing.process_annotations import (DBSETS_COL, MODULE_STEPS_FORM_TAG, FUNCTION_HEATMAP_FORM_TAG,
        ETC_MODULE_DF_TAG, FILES_NAMES, build_module_net, fill_product_dfs, make_product_df,
        get_phylum_and_most_specific,
        make_strings_no_repeats, get_annotation_ids_by_row, rename_genomes_to_taxa)
    from models.heatmap import make_product_heatmap

__authors__ = ["Madeline Scyphers", "Rory Flynn"]
__copyright__ = "Copyright 2024, Wrighton Lab"

logger = logging.getLogger("dram2_log.viz")

pn.extension('tabulator', template='bootstrap')

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


# def make_plot(
#         module_df: pd.DataFrame,
#         etc_df: pd.DataFrame,
#         function_df: pd.DataFrame,
# ):
#     """
#     Make a product heatmap group from the module_coverage_df, etc_coverage_df, and functional_df
#
#     Parameters
#     ----------
#     module_df : pd.DataFrame
#         A dataframe of module coverage information
#     etc_df : pd.DataFrame
#         A dataframe of ETC coverage information
#     function_df : pd.DataFrame
#         A dataframe of functional coverage information
#
#     Returns
#     -------
#     pn.Row
#         A panel row of heatmaps
#     """
#     # chart_type = pn.widgets.Select(name="Chart Type", options=["Heatmap", "Barchart"], value="Heatmap")
#
#     coverage_slider = pn.widgets.RangeSlider(name="Coverage Range", start=0, end=1, value=(0, 1), step=0.01)
#
#     widgets = [
#         # chart_type
#     ]
#
#     p = make_product_heatmap(module_df, etc_df, function_df)
#     # return p
#     # # plot.save(output_dir / "product.html", resources=INLINE_RESOURCES)
#     return pn.template.FastListTemplate(
#         title="DRAM2 Product Visualization",
#         main=[p],
#         sidebar=[
#             *widgets
#         ]
#     )

TAXONOMY_RANKS = {
    "domain": r"d__.*?(?=p__)",
    "phylum": r"p__.*?(?=c__)",
    "class": r"c__.*?(?=o__)",
    "order": r"o__.*?(?=f__)",
    "family": r"f__.*?(?=g__)",
    "genus": r"g__.*?(?=s__)",
    "species": r"s__.*"
}


class Dashboard(pn.viewable.Viewer):
    """
    A parameterized class for the product visualization
    """

    min_coverage = param.Number(default=0, bounds=(0, 1), label="Minimum Coverage")
    y_axis_col = param.Selector(objects=["genome", "taxonomy"], default="genome", label="Y Axis Column")
    taxonomy_ranks = param.ListSelector(
        # objects={"domain": "d__", "phylum": "p__", "class": "c__", "order": "o__", "family": "f__", "genus": "g__", "species": "s__"},
        objects=TAXONOMY_RANKS,
        default=list(TAXONOMY_RANKS.values()),
        # default=["domain", "phylum", "class", "order", "family", "genus", "species"],
        label="Taxonomy Group By Ranks",
    )

    view = param.ClassSelector(class_=pn.template.FastListTemplate)

    plot_view = param.ClassSelector(class_=pn.Row)

    def __init__(self, module_df: pd.DataFrame, etc_df: pd.DataFrame, function_df: pd.DataFrame):
        super().__init__()
        self.module_df = module_df
        self.etc_df = etc_df
        self.function_df = function_df
        self.plot_view = pn.Row()

        self.redraw_button = pn.widgets.Button(name="Redraw", button_type="primary")
        self.reset_button = pn.widgets.Button(name="Reset Filters", button_type="primary")
        self.reset_button.on_click(self.reset_filters)

        self.redraw_button.on_click(self.make_plot)

        regex = "".join([f"(?P<{rank}>{regex})" for rank, regex in TAXONOMY_RANKS.items()])
        # self.module_df["taxonomy"].str.extractall(regex)

        self.module_df = pd.concat([self.module_df, self.module_df["taxonomy"].str.extractall(regex).reset_index(drop=True)], axis=1)
        self.etc_df = pd.concat([self.etc_df, self.etc_df["taxonomy"].str.extractall(regex).reset_index(drop=True)], axis=1)
        self.function_df = pd.concat([self.function_df, self.function_df["taxonomy"].str.extractall(regex).reset_index(drop=True)], axis=1)

        if "taxonomy" in self.module_df.columns:

            self.taxonomy_filter = pn.Tabs(
                *[pn.widgets.MultiChoice(name=rank,
                                         value=self.module_df[rank].unique().tolist(),
                                         options=self.module_df[rank].unique().tolist(),
                                         )
                  for rank in TAXONOMY_RANKS]
            )
        else:
            self.taxonomy_filter = None

        self.make_plot()


    # @param.depends("min_coverage", "y_axis_col", "taxonomy_ranks", watch=True)
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
        print("making plot")

        taxonomy_labels = None
        if self.taxonomy_ranks != list(TAXONOMY_RANKS.values()):
            if "taxonomy" not in self.module_df.columns:
                additional_sidebar.append(pn.pane.Alert("Taxonomy column not found in supplied data", alert_type='warning'))
            else:
                regex = "|".join([regex for regex in self.param.taxonomy_ranks.objects if regex not in self.taxonomy_ranks])
                module_df["taxonomy"] = module_df["taxonomy"].str.replace(regex, "", regex=True)
                etc_df["taxonomy"] = etc_df["taxonomy"].str.replace(regex, "", regex=True)
                function_df["taxonomy"] = function_df["taxonomy"].str.replace(regex, "", regex=True)

        module_df, etc_df, function_df = self.filter_by_taxonomy(module_df, etc_df, function_df)

        charts = make_product_heatmap(module_df, etc_df, function_df, y_col=self.y_axis_col, taxonomy_labels=taxonomy_labels)

        self.plot_view[:] = charts

        # return p
        # # plot.save(output_dir / "product.html", resources=INLINE_RESOURCES)
        if "taxonomy" in self.module_df.columns:
            additional_sidebar.append(self.param.taxonomy_ranks)
            additional_sidebar.append(self.taxonomy_filter)

        self.view = pn.template.FastListTemplate(
            title="DRAM2 Product Visualization",
            # main=[self.plot_view],
            main=[
                pn.Tabs(
                    ("Heatmap", self.plot_view),
                    ("Module Coverage DF", pn.widgets.Tabulator(module_df, page_size=50)),
                    ("ETC Coverage DF", pn.widgets.Tabulator(etc_df, page_size=50)),
                    ("Function DF", pn.widgets.Tabulator(self.function_df, page_size=50))
                )
            ],
            sidebar=[
                pn.Row(self.redraw_button, self.reset_button),
                self.param.min_coverage,
                self.param.y_axis_col,
                *additional_sidebar,
            ]
        )

    def reset_filters(self, event=None):
        """
        Reset the taxonomy filter
        """
        self.min_coverage = self.param.min_coverage.default
        self.y_axis_col = self.param.y_axis_col.default
        self.taxonomy_ranks = self.param.taxonomy_ranks.default

        if self.taxonomy_filter is not None:
            for multiselect in self.taxonomy_filter:
                multiselect.value = multiselect.options

    def filter_by_taxonomy(self, module_df, etc_df, function_df):
        """
        Filter the dataframes by selected taxonomy ranks
        """
        if self.taxonomy_filter is None:
            return module_df, etc_df, function_df
        for multiselect in self.taxonomy_filter:
            module_df = module_df.loc[module_df[multiselect.name].isin(multiselect.value)]
            etc_df = etc_df.loc[etc_df[multiselect.name].isin(multiselect.value)]
            function_df = function_df.loc[function_df[multiselect.name].isin(multiselect.value)]
            # print(multiselect.name)

        for multiselect in self.taxonomy_filter:
            # check if all values in multiselect are still valid after filtering
            if not all([i in module_df[multiselect.name].unique() for i in multiselect.value]):
                multiselect.value = [i for i in multiselect.value if i in module_df[multiselect.name].unique()]

        return module_df, etc_df, function_df


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

    if "taxonomy" in annotations:
        tax_df = annotations[[groupby_column, "taxonomy"]].drop_duplicates()
        tax_df.rename(columns={groupby_column: "genome"}, inplace=True)

        module_coverage_df =tax_df.merge(module_coverage_df, on="genome", how="left")
        etc_coverage_df = tax_df.merge(etc_coverage_df, on="genome", how="left")
        function_df = tax_df.merge(function_df, on="genome", how="left")

    # module_coverage_df.to_csv(output_dir / "module_coverage.tsv", sep="\t", index=False)
    # etc_coverage_df.to_csv(output_dir / "etc_coverage.tsv", sep="\t", index=False)
    # function_df.to_csv(output_dir / "function_coverage.tsv", sep="\t", index=False)
    product_df = make_product_df(module_coverage_df, etc_coverage_df, function_df)

    if labels is not None:
        function_df = rename_genomes_to_taxa(function_df, labels)

    app = Dashboard(module_coverage_df, etc_coverage_df, function_df)

    # plot = make_plot(
    #     module_coverage_df,
    #     etc_coverage_df,
    #     function_df,
    # )

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
