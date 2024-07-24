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

try:
    from .definitions import DEFAULT_GROUPBY_COLUMN, BACKUP_GROUPBY_COLUMN, HEATMAP_MODULES
    from .proccessing.process_annotations import (DBSETS_COL, MODULE_STEPS_FORM_TAG, FUNCTION_HEATMAP_FORM_TAG,
        ETC_MODULE_DF_TAG, FILES_NAMES, build_module_net, fill_product_dfs, make_product_df,
        get_phylum_and_most_specific,
        make_strings_no_repeats, get_annotation_ids_by_row, rename_genomes_to_taxa)
    from .apps.heatmap import make_product_heatmap, Dashboard
except ImportError:
    from definitions import DEFAULT_GROUPBY_COLUMN, BACKUP_GROUPBY_COLUMN, HEATMAP_MODULES
    from processing.process_annotations import (DBSETS_COL, MODULE_STEPS_FORM_TAG, FUNCTION_HEATMAP_FORM_TAG,
        ETC_MODULE_DF_TAG, FILES_NAMES, build_module_net, fill_product_dfs, make_product_df,
        get_phylum_and_most_specific,
        make_strings_no_repeats, get_annotation_ids_by_row, rename_genomes_to_taxa)
    from apps.heatmap import make_product_heatmap, Dashboard

__authors__ = ["Madeline Scyphers", "Rory Flynn"]
__copyright__ = "Copyright 2024, Wrighton Lab"

logger = logging.getLogger("dram2_log.viz")


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
