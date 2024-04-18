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
    from .definitions import DEFAULT_GROUPBY_COLUMN
    from .proccessing.process_annotations import (DBSETS_COL, MODULE_STEPS_FORM_TAG, FUNCTION_HEATMAP_FORM_TAG,
        ETC_MODULE_DF_TAG, FILES_NAMES, build_module_net, fill_product_dfs, make_product_df,
        get_phylum_and_most_specific,
        make_strings_no_repeats, get_annotation_ids_by_row, rename_genomes_to_taxa)
    from .models.heatmap import make_product_heatmap
except ImportError:
    from definitions import DEFAULT_GROUPBY_COLUMN
    from processing.process_annotations import (DBSETS_COL, MODULE_STEPS_FORM_TAG, FUNCTION_HEATMAP_FORM_TAG,
        ETC_MODULE_DF_TAG, FILES_NAMES, build_module_net, fill_product_dfs, make_product_df,
        get_phylum_and_most_specific,
        make_strings_no_repeats, get_annotation_ids_by_row, rename_genomes_to_taxa)
    from models.heatmap import make_product_heatmap

__authors__ = ["Madeline Scyphers", "Rory Flynn"]
__copyright__ = "Copyright 2024, Wrighton Lab"

logger = logging.getLogger("dram2_log.viz")

pn.extension(template='bootstrap')

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


class Dashboard(pn.viewable.Viewer):
    """
    A parameterized class for the product visualization
    """

    min_coverage = param.Number(default=0, bounds=(0, 1), label="Minimum Coverage")

    view = param.ClassSelector(class_=pn.template.FastListTemplate)

    def __init__(self, module_df: pd.DataFrame, etc_df: pd.DataFrame, function_df: pd.DataFrame):
        super().__init__()
        self.module_df = module_df
        self.etc_df = etc_df
        self.function_df = function_df

        # self.view = pn.template.FastListTemplate(
        #             title="DRAM2 Product Visualization",
        #             main=[p],
        #             sidebar=[
        #                 self.param.min_coverage
        #             ]
        #         )
        self.make_plot()

    @pn.depends("min_coverage", watch=True)
    def make_plot(self):
        """
        Make the product plot
        """
        if self.min_coverage > 0:
            print(self.min_coverage)
            module_df = self.module_df[self.module_df["step_coverage"] >= self.min_coverage]
            etc_df = self.etc_df[self.etc_df["percent_coverage"] >= self.min_coverage]
        else:
            module_df = self.module_df
            etc_df = self.etc_df
        print("making plot")
        p = make_product_heatmap(module_df, etc_df, self.function_df)



        # return p
        # # plot.save(output_dir / "product.html", resources=INLINE_RESOURCES)
        self.view = pn.template.FastListTemplate(
            title="DRAM2 Product Visualization",
            main=[p],
            sidebar=[
                self.param.min_coverage
            ]
        )


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

    module_coverage_df, etc_coverage_df, function_df = fill_product_dfs(
        annotations=annotations,
        module_nets=module_nets,
        etc_module_df=etc_module_df,
        function_heatmap_form=function_heatmap_form,
        annotation_ids_by_row=annotation_ids_by_row,
        groupby_column=groupby_column
    )
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
