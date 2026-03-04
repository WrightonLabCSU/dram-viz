"""
===================
DRAM Visualizations
===================

Script that generates a product visualization from the DRAM output.
"""

from __future__ import annotations

import logging
from pathlib import Path

import click
import panel as pn
import polars as pl

from dram_viz.apps.heatmap import Dashboard
from dram_viz.definitions import (
    BACKUP_ALIAS_COLUMNS,
    BACKUP_FASTA_COLUMNS,
    BACKUP_GROUPBY_COLUMNS,
    BACKUP_LABEL_COLUMNS,
    DEFAULT_ALIAS_COLUMN,
    DEFAULT_FASTA_COLUMN,
    DEFAULT_GROUPBY_COLUMN,
    DEFAULT_LABEL_COLUMN,
)
from dram_viz.processing.process_annotations import (
    build_tax_edge_df_pl,
    build_tax_tree_selected_recurse,
    build_taxonomy_df_pl,
    build_tree_pl,
)
from dram_viz.rule_parser.src.rules import (
    ID_EXPR_DICT,
    CompiledRules,
    build_present_map,
    evaluate_cycles,
)

logger = logging.getLogger("dram.viz")

pn.config.reuse_sessions = True
pn.config.global_loading_spinner = True


def get_column_name(
    name: str,
    backup_names: list[str],
    name_of_data: str,
    df: pl.LazyFrame | pl.DataFrame = None,
    cols: list[str] = None,
    raise_error: bool = True,
) -> str:
    assert df is not None or cols is not None, (
        "Must provide either a dataframe or a list of columns to search for the column name"
    )
    if cols is None:
        if isinstance(df, pl.LazyFrame):
            cols = df.collect_schema().names()
        else:
            cols = df.columns
    if name in cols:
        return name
    for backup_name in backup_names:
        if backup_name in cols:
            return backup_name
    if raise_error:
        raise ValueError(
            f"None of the following column names were found in the {name_of_data}: {[name] + backup_names}"
        )
    else:
        return name


@click.command()
@click.option(
    "--rules_tsv",
    "-r",
    type=click.Path(exists=True),
    default=Path(__file__).parent / "data/rules.tsv",
    help="This is an optional path to a rules file with strict formatting. It will over write the original rules file that is stored with the script.",
)
@click.option("--annotations", "-a", type=Path, help="Path to the annotations tsv file")
@click.option(
    "--fasta_column",
    "-f",
    type=str,
    default=DEFAULT_FASTA_COLUMN,
    help="Column from annotations file to use as fasta names",
)
@click.option(
    "--output_dir",
    "-o",
    type=Path,
    help="Path to the output directory",
    default=Path.cwd().resolve(),
)
@click.option(
    "--mapping",
    "-m",
    type=Path,
    help="Path to mapping file of gene abundances in samples. Must have a column named 'Geneid' that matches the query_id column in the annotations file.",
)
@click.option(
    "--label_column",
    "-l",
    type=str,
    default=DEFAULT_LABEL_COLUMN,
    help="Column in rules file to use as labels in visualization",
)
@click.option(
    "--alias_column",
    "-ac",
    type=str,
    default=DEFAULT_ALIAS_COLUMN,
    help="Column from rules file to use as alias in rules. This allows you to alias a rule on one line and then use that alias in a different line.",
)
@click.option(
    "--group_colunm",
    "-g",
    type=str,
    help="Column in rules file to group by in visualization (creates separate heatmap for each unique value in this column)",
    default=DEFAULT_GROUPBY_COLUMN,
)
@click.option(
    "--dashboard",
    "-d",
    is_flag=True,
    show_default=True,
    default=False,
)
@click.option(
    "--save-dataframes",
    "-sd",
    is_flag=True,
    show_default=True,
    default=False,
)
def main(
    rules_tsv,
    annotations,
    fasta_column=DEFAULT_FASTA_COLUMN,
    output_dir=None,
    mapping=None,
    label_column=DEFAULT_LABEL_COLUMN,
    alias_column=DEFAULT_ALIAS_COLUMN,
    group_colunm=DEFAULT_GROUPBY_COLUMN,
    dashboard=False,
    save_dataframes=False,
):
    """
    Make a product heatmap visualization from the DRAM output.
    """
    import time

    s = time.time()

    rules_lf = pl.scan_csv(
        rules_tsv, separator="\t", infer_schema_length=None
    ).fill_null("")

    rules_cols = rules_lf.collect_schema().names()
    group_colunm = get_column_name(
        group_colunm, BACKUP_GROUPBY_COLUMNS, name_of_data="rules", cols=rules_cols
    )
    label_column = get_column_name(
        label_column, BACKUP_LABEL_COLUMNS, name_of_data="rules", cols=rules_cols
    )
    alias_column = get_column_name(
        alias_column,
        BACKUP_ALIAS_COLUMNS,
        name_of_data="rules",
        cols=rules_cols,
        raise_error=False,
    )
    if "long_name" not in rules_cols:
        rules_lf = rules_lf.with_columns(long_name=pl.col(label_column))

    rules_lf = rules_lf.with_columns(
        pl.col("long_name").fill_null(pl.col(label_column)).alias("long_name")
    )

    anno = pl.read_csv(
        annotations,
        separator="\t",
        infer_schema_length=10_000,
        # columns=list(ID_EXPR_DICT.keys()) + [groupby_column]
    )
    fasta_column = get_column_name(
        fasta_column, BACKUP_FASTA_COLUMNS, name_of_data="annotations", df=anno
    )
    anno = anno.rename({fasta_column: "genome"})

    kw = dict(
        rules=rules_lf,
        label_col=label_column,
        alias_col=alias_column,
        rules_col="rule",
        allow_visualize_functions=True,
    )
    # kw = dict(rules_path=rules_path, label_col="module", parent_col=alias_column, rules_col="rule")
    compiled = CompiledRules.from_rules(**kw)
    logger.info(f"Compiled rules in {time.time() - s} seconds")
    if mapping:
        mapping_df = pl.read_csv(mapping, separator="\t", ignore_errors=True).fill_null(
            0
        )
        sample_names = mapping_df.columns[1:]
        mapping_df = mapping_df.with_columns(
            summed_sample_abundance=pl.sum_horizontal(sample_names)
        )
        anno = anno.join(
            mapping_df.select(["Geneid", "summed_sample_abundance"]).rename(
                {"Geneid": "query_id"}
            ),
            on=["query_id"],
            how="left",
        )
        del mapping_df

        samples, present_map, anno_df = build_present_map(
            anno,
            sample_col="genome",
            besthit_cols=list(ID_EXPR_DICT.keys()),
            needed_features=compiled.needed_features,
            additional_cols=["summed_sample_abundance"],
        )
    else:
        samples, present_map = build_present_map(
            anno,
            sample_col="genome",
            besthit_cols=list(ID_EXPR_DICT.keys()),
            needed_features=compiled.needed_features,
        )

    logger.info(f"Built present map in {time.time() - s} seconds")

    dfs = evaluate_cycles(
        compiled=compiled,
        samples=samples,
        present_map=present_map,
        annotations=anno,
        sample_col="genome",
        additional_cols=["long_name"],
        group_col=group_colunm,
    )
    logger.info("Evaluated all rules in:")
    logger.info(time.time() - s)

    if mapping:
        mapped_dfs = {}
        for group, frame in compiled.df.group_by(group_colunm, maintain_order=True):
            group = group[0]
            mapped_dfs[group] = []
            for rn in frame.select(label_column).unique().to_series():
                mapped_dfs[group].append(
                    anno_df.filter(pl.col("hit").is_in(compiled.features_by_rules[rn]))
                    .group_by("genome")
                    .agg(pl.col("summed_sample_abundance").sum())
                    .with_columns(name=pl.lit(rn))
                )
            df = (
                dfs[group]
                .join(
                    pl.concat(mapped_dfs[group]),
                    on=["genome", label_column],
                    how="left",
                )
                .with_columns(pl.col("summed_sample_abundance").fill_null(0))
            )
            dfs[group] = df
            # dfs[group] = dfs[group].join(pl.concat(mapped_dfs[group]), on=["genome", label_column], how="left").with_columns(pl.col("summed_sample_abundance").fill_null(0))

    extra_cols = [
        col for col in ["Completeness", "Contamination"] if col in anno.columns
    ]
    if extra_cols:
        df = anno.unpivot(
            index="genome", on=extra_cols, variable_name=label_column
        ).unique()
        # df = df.rename({groupby_column: "genome"})
        # We do this reorder the dfs dictionary to ensure that the metadata dataframe is the first one
        d = {"Meta": df}
        d.update(dfs)
        dfs = d

    tax_tree_data = None
    selected_tax_tree = None
    if "taxonomy" in anno:
        tax_df = build_taxonomy_df_pl(anno, "genome")

        tax_edge_df, tax_df = build_tax_edge_df_pl(tax_df)

        tax_tree_data = build_tree_pl(
            tax_edge_df,
            state={"opened": False, "selected": True},
            id_cb=lambda source, child, parent_id: f"{parent_id};{child}",
        )
        selected_tax_tree = build_tax_tree_selected_recurse(tax_tree_data)

        for key, df in dfs.items():
            dfs[key] = tax_df.join(df, on="genome", how="left")

    del anno

    if save_dataframes:
        for key, df in dfs.items():
            df.write_csv(output_dir / f"{key}_df.tsv", separator="\t")
            logger.info(f"Saved {key} dataframe to {output_dir / f'{key}_df.tsv'}")

    kw = dict(
        dfs=dfs,
        taxanomy_tree_data=tax_tree_data,
        selected_tax_tree=selected_tax_tree,
        mapping=bool(mapping),
    )
    logger.info(
        f"Finished all processing in: {time.time() - s} seconds, starting visualization"
    )
    if dashboard:
        pn.serve(
            lambda: Dashboard(**kw),
            port=5006,
        )
    else:
        Dashboard(**kw)


if __name__ == "__main__":
    main()
