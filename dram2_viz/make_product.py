"""
===================
DRAM Visualizations
===================

Script that generates a product visualization from the DRAM output.
"""
from __future__ import annotations
import logging
from pathlib import Path
from typing import Optional

import click
from bokeh.resources import INLINE
import pandas as pd
import panel as pn

from dram2_viz.definitions import DEFAULT_GROUPBY_COLUMN, BACKUP_GROUPBY_COLUMN, HEATMAP_MODULES, DBSETS_COL, \
    MODULE_STEPS_FORM_TAG, FUNCTION_HEATMAP_FORM_TAG, ETC_MODULE_DF_TAG, FILES_NAMES
from dram2_viz.processing.process_annotations import (build_module_net, fill_product_dfs, make_product_df,
                                                      get_phylum_and_most_specific,
                                                      make_strings_no_repeats, get_annotation_ids_by_row,
                                                      rename_genomes_to_taxa, build_tree, build_tax_edge_df, build_taxonomy_df)
from dram2_viz.apps.heatmap import Dashboard


logger = logging.getLogger("dram2_log.viz")


@click.command()
@click.option(
    "--annotations",
    "-a",
    type=Path,
    help="Path to the annotations tsv file"
)
@click.option(
    "--groupby-column",
    "-g",
    type=str,
    default=DEFAULT_GROUPBY_COLUMN,
    help="Column to group by"
)
@click.option(
    "--output-dir",
    "-o",
    type=Path,
    help="Path to the output directory",
    default=Path.cwd().resolve()
)
@click.option(
    "--module-steps-form",
    type=Path,
    help="Path to Module Step Database TSV",
    default=FILES_NAMES[MODULE_STEPS_FORM_TAG]
)
@click.option(
    "--etc-steps-form",
    type=Path,
    help="Path to ETC Step Database TSV",
    default=FILES_NAMES[ETC_MODULE_DF_TAG]
)
@click.option(
    "--function-steps-form",
    type=Path,
    help="Path to Function Step Database TSV",
    default=FILES_NAMES[FUNCTION_HEATMAP_FORM_TAG]
)
@click.option(
    "--dashboard",
    "-d",
    is_flag=True,
)
def main(annotations,
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
    annotations = pd.read_csv(annotations, sep="\t", index_col=0)

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
        annotations_df=annotations,
        module_nets=module_nets,
        etc_module_df=etc_module_df,
        function_heatmap_form=function_heatmap_form,
        annotation_ids_by_row=annotation_ids_by_row,
        groupby_column=groupby_column
    )

    tax_tree_data = None
    if "taxonomy" in annotations:
        tax_df = build_taxonomy_df(annotations, groupby_column)

        tax_edge_df = build_tax_edge_df(tax_df)

        tax_tree_data = build_tree(tax_edge_df, state={"opened": False, "selected": True}, id_cb=lambda source, child, parent_id: f"{parent_id}; {child}")

        module_coverage_df = tax_df.merge(module_coverage_df, on="genome", how="left")
        etc_coverage_df = tax_df.merge(etc_coverage_df, on="genome", how="left")
        function_df = tax_df.merge(function_df, on="genome", how="left")

    product_df = make_product_df(module_coverage_df, etc_coverage_df, function_df)

    if labels is not None:
        function_df = rename_genomes_to_taxa(function_df, labels)

    app = Dashboard(module_coverage_df, etc_coverage_df, function_df, tax_tree_data=tax_tree_data)

    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)
    app.plot_view.save(output_dir / "product.html", resources=INLINE)
    product_df.to_csv(output_dir / "product.tsv", sep="\t", index=False)
    if dashboard:
        pn.serve(app.view, port=5006)
        # plot.show(port=5007)
    logger.info("Completed visualization")


if __name__ == "__main__":
    main()
