from __future__ import annotations

import logging
import re
from collections import Counter, defaultdict
from itertools import chain, tee
from typing import Optional

import networkx as nx
import numpy as np
import pandas as pd

from dram_viz.definitions import (
    DBSETS_COL,
    DEFAULT_GROUPBY_COLUMN,
    ETC_COVERAGE_COLUMNS,
    ID_FUNCTION_DICT,
    KO_REGEX,
    TAXONOMY_LEVELS,
)

logger = logging.getLogger("dram.viz")


def build_module_net(module_df):
    """Starts with a data from including a single module"""
    # build net from a set of module paths
    num_steps = max([int(i.split(",")[0]) for i in set(module_df["path"])])
    module_net = nx.DiGraph(
        num_steps=num_steps,
        module_id=list(module_df["module"])[0],
        module_name=list(module_df["module_name"])[0],
    )
    # go through all path/step combinations
    for module_path, frame in module_df.groupby("path"):
        split_path = [int(i) for i in module_path.split(",")]
        step = split_path[0]
        module_net.add_node(module_path, kos=set(frame["ko"]))
        # add incoming edge
        if step != 0:
            module_net.add_edge("end_step_%s" % (step - 1), module_path)
        # add outgoing edge
        module_net.add_edge(module_path, "end_step_%s" % step)
    return module_net


def get_module_step_coverage(kos, module_net):
    # prune network based on what kos were observed
    pruned_module_net = module_net.copy()
    module_kos_present = set()
    for node, data in module_net.nodes.items():
        if "kos" in data:
            ko_overlap = data["kos"] & kos
            if len(ko_overlap) == 0:
                pruned_module_net.remove_node(node)
            else:
                module_kos_present = module_kos_present | ko_overlap
    # count number of missing steps
    missing_steps = 0
    for node, data in pruned_module_net.nodes.items():
        if ("end_step" in node) and (pruned_module_net.in_degree(node) == 0):
            missing_steps += 1
    # get statistics
    num_steps = pruned_module_net.graph["num_steps"] + 1
    num_steps_present = num_steps - missing_steps
    coverage = num_steps_present / num_steps
    return num_steps, num_steps_present, coverage, sorted(module_kos_present)


def make_module_coverage_df(annotation_df, module_nets):
    kos_to_genes = defaultdict(list)
    ko_id: Optional[str] = None
    ko_id_names: list[str] = ["kegg_id", "kofam_id", "ko_id"]
    for id in ko_id_names:
        if id in annotation_df:
            ko_id = id
            break
    if ko_id is None:
        raise ValueError(
            f"""
            No KEGG or KOfam id column could be found.
            These names were tried: {', '.join(ko_id_names)}.
            """
        )
    for gene_id, ko_list in annotation_df[ko_id].items():
        if type(ko_list) is str:
            for ko in ko_list.split(","):
                kos_to_genes[ko].append(gene_id)
    coverage_dict = {}
    for i, (module, net) in enumerate(module_nets.items()):
        (
            module_steps,
            module_steps_present,
            module_coverage,
            module_kos,
        ) = get_module_step_coverage(set(kos_to_genes.keys()), net)
        module_genes = sorted([gene for ko in module_kos for gene in kos_to_genes[ko]])
        coverage_dict[module] = [
            net.graph["module_name"],
            module_steps,
            module_steps_present,
            module_coverage,
            len(module_kos),
            ",".join(module_kos),
            ",".join(module_genes),
        ]
    coverage_df = pd.DataFrame.from_dict(
        coverage_dict,
        orient="index",
        columns=[
            "module_name",
            "steps",
            "steps_present",
            "step_coverage",
            "ko_count",
            "kos_present",
            "genes_present",
        ],
    )
    return coverage_df


def pairwise(iterable):
    """s -> (s0, s1), (s1, s2), (s2, s3), ..."""
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def first_open_paren_is_all(str_):
    """Go through string and return true"""
    curr_level = 1
    for char in str_[1:-1]:
        if char == ")":
            curr_level -= 1
        elif char == "(":
            curr_level += 1
        if curr_level == 0:
            return False
    return True


def split_into_steps(definition, split_char=" "):
    """Very fancy split on string of chars"""
    curr_level = 0
    step_starts = [-1]
    for i, char in enumerate(definition):
        if char == "(":
            curr_level += 1
        if char == ")":
            curr_level -= 1
        if (curr_level == 0) and (char in split_char):
            step_starts.append(i)
    step_starts.append(len(definition))
    steps = list()
    for a, b in pairwise(step_starts):
        step = definition[a + 1 : b]
        if step.startswith("(") and step.endswith(")"):
            if first_open_paren_is_all(step):
                step = step[1:-1]
        steps.append(step)
    return steps


def is_ko(ko):
    return re.match(KO_REGEX, ko) is not None


def make_module_network(definition, network: nx.DiGraph = None, parent_nodes=("start",)):
    # TODO: Figure out how to add 'end' step to last step at end
    if network is None:
        network = nx.DiGraph()
    last_steps = []
    for step in split_into_steps(definition, ","):
        prev_steps = parent_nodes
        for substep in split_into_steps(step, "+"):
            if is_ko(substep):
                for prev_step in prev_steps:
                    network.add_edge(prev_step, substep)
                prev_steps = [substep]
            else:
                network, prev_steps = make_module_network(substep, network, prev_steps)
        last_steps += prev_steps
    return network, last_steps


def get_module_coverage(module_net: nx.DiGraph, genes_present: set):
    max_coverage = -1
    max_coverage_genes = list()
    max_coverage_missing_genes = list()
    max_path_len = 0
    for net_path in nx.all_simple_paths(module_net, source="start", target="end"):
        net_path = set(net_path[1:-1])
        overlap = net_path & genes_present
        coverage = len(overlap) / len(net_path)
        if coverage > max_coverage:
            max_coverage = coverage
            max_coverage_genes = overlap
            max_coverage_missing_genes = net_path - genes_present
            max_path_len = len(net_path)
    return (
        max_path_len,
        len(max_coverage_genes),
        max_coverage,
        max_coverage_genes,
        max_coverage_missing_genes,
    )


def make_module_coverage_frame(annotations_df, module_nets, groupby_column=DEFAULT_GROUPBY_COLUMN):
    # go through each scaffold to check for modules
    module_coverage_dict = dict()
    for group, frame in annotations_df.groupby(groupby_column, sort=False):
        module_coverage_dict[group] = make_module_coverage_df(frame, module_nets)
    module_coverage = pd.concat(module_coverage_dict)
    module_coverage.index = module_coverage.index.set_names(["genome", "module"])
    return module_coverage.reset_index()


def make_etc_coverage_df(
    etc_module_df,
    annotation_ids_by_row: pd.DataFrame,
    groupby_column=DEFAULT_GROUPBY_COLUMN,
):
    etc_coverage_df_rows = list()
    for _, module_row in etc_module_df.iterrows():
        definition = module_row["definition"]
        # remove optional subunits
        definition = re.sub(r"-K\d\d\d\d\d", "", definition)
        module_net, _ = make_module_network(definition)
        # add end node
        no_out = [node for node in module_net.nodes() if module_net.out_degree(node) == 0]
        for node in no_out:
            module_net.add_edge(node, "end")
        # go through each genome and check pathway coverage
        for group, frame in annotation_ids_by_row.groupby(groupby_column):
            # get annotation genes
            grouped_ids = set(get_all_annotation_ids(frame).keys())
            (
                path_len,
                path_coverage_count,
                path_coverage_percent,
                genes,
                missing_genes,
            ) = get_module_coverage(module_net, grouped_ids)
            complex_module_name = "Complex %s: %s" % (
                module_row["complex"].replace("Complex ", ""),
                module_row["module_name"],
            )
            etc_coverage_df_rows.append(
                [
                    module_row["module_id"],
                    module_row["module_name"],
                    module_row["complex"].replace("Complex ", ""),
                    group,
                    path_len,
                    path_coverage_count,
                    path_coverage_percent,
                    ",".join(sorted(genes)),
                    ",".join(sorted(missing_genes)),
                    complex_module_name,
                ]
            )
    return pd.DataFrame(etc_coverage_df_rows, columns=ETC_COVERAGE_COLUMNS)


def make_functional_df(
    annotation_ids_by_row,
    function_heatmap_form,
    groupby_column=DEFAULT_GROUPBY_COLUMN,
):
    # clean up function heatmap form
    function_heatmap_form = function_heatmap_form.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    function_heatmap_form = function_heatmap_form.fillna("")
    # build dict of ids per genome
    genome_to_id_dict = dict()
    for genome, frame in annotation_ids_by_row.groupby(groupby_column, sort=False):
        id_list = get_all_annotation_ids(frame).keys()
        genome_to_id_dict[genome] = set(id_list)
    # build long from data frame
    rows = list()
    for function, frame in function_heatmap_form.groupby("function_name", sort=False):
        for bin_name, id_set in genome_to_id_dict.items():
            presents_in_bin = list()
            functions_present = set()
            for _, row in frame.iterrows():
                function_id_set = set([i.strip() for i in row.function_ids.strip().split(",")])
                present_in_bin = id_set & function_id_set
                functions_present = functions_present | present_in_bin
                presents_in_bin.append(len(present_in_bin) > 0)
            function_in_bin = np.all(presents_in_bin)
            row = frame.iloc[0]
            rows.append(
                [
                    row.category,
                    row.subcategory,
                    row.function_name,
                    ", ".join(functions_present),
                    "; ".join(get_ordered_uniques(frame.long_function_name)),
                    "; ".join(get_ordered_uniques(frame.gene_symbol)),
                    bin_name,
                    function_in_bin,
                    "%s: %s" % (row.category, row.function_name),
                ]
            )
    return pd.DataFrame(
        rows,
        columns=list(function_heatmap_form.columns) + ["genome", "present", "category_function_name"],
    )


# TODO: refactor this to handle splitting large numbers of genomes into multiple heatmaps here
def fill_product_dfs(
    annotations_df,
    module_nets,
    etc_module_df,
    function_heatmap_form,
    annotation_ids_by_row: pd.DataFrame,
    groupby_column=DEFAULT_GROUPBY_COLUMN,
):
    module_coverage_frame = make_module_coverage_frame(annotations_df, module_nets, groupby_column)

    # make ETC frame
    etc_coverage_df = make_etc_coverage_df(etc_module_df, annotation_ids_by_row, groupby_column)

    # make functional frame
    function_df = make_functional_df(
        annotation_ids_by_row,
        function_heatmap_form,
        groupby_column,
    )

    return module_coverage_frame, etc_coverage_df, function_df


def make_product_df(module_coverage_frame, etc_coverage_df, function_df):
    liquor_df = pd.concat(
        [
            module_coverage_frame.pivot(index="genome", columns="module_name", values="step_coverage"),
            etc_coverage_df.pivot(
                index="genome",
                columns="complex_module_name",
                values="percent_coverage",
            ),
            function_df.pivot(
                index="genome",
                columns="category_function_name",
                values="present",
            ),
        ],
        axis=1,
        sort=False,
    )
    return liquor_df


def get_phylum_and_most_specific(taxa_str):
    taxa_ranks = [i[3:] for i in taxa_str.split(";")]
    phylum = taxa_ranks[1]
    most_specific_rank = TAXONOMY_LEVELS[sum([len(i) > 0 for i in taxa_ranks]) - 1]
    most_specific_taxa = taxa_ranks[sum([len(i) > 0 for i in taxa_ranks]) - 1]
    if most_specific_rank == "d":
        return "d__%s;p__" % most_specific_taxa
    if most_specific_rank == "p":
        return "p__%s;c__" % most_specific_taxa
    else:
        return "p__%s;%s__%s" % (phylum, most_specific_rank, most_specific_taxa)


def make_strings_no_repeats(genome_taxa_dict: dict):
    labels = dict()
    seen = Counter()
    for genome, taxa_string in genome_taxa_dict.items():
        final_taxa_string = "%s_%s" % (taxa_string, str(seen[taxa_string]))
        seen[taxa_string] += 1
        labels[genome] = final_taxa_string
    return labels


def get_annotation_ids_by_row(data):
    functions = {i: j for i, j in ID_FUNCTION_DICT.items() if i in data.columns}
    missing = [i for i in ID_FUNCTION_DICT if i not in data.columns]
    logger.info(
        "Note: the fallowing id fields "
        f"were not in the annotations file and are not being used: {missing},"
        f" but these are {list(functions.keys())}"
    )
    out = data.apply(
        lambda x: {i for k, v in functions.items() if not pd.isna(x[k]) for i in v(str(x[k])) if not pd.isna(i)},
        axis=1,
    )
    return out


def get_all_annotation_ids(data):
    annotation_series = data[DBSETS_COL]
    annotation_series.apply(list)
    out = Counter(chain(*annotation_series.values))
    return out


def rename_genomes_to_taxa(function_df, labels):
    function_df = function_df.copy()
    new_genome_column = [labels[i] for i in function_df["genome"]]
    function_df["genome"] = pd.Series(new_genome_column, index=function_df.index)
    return function_df


class DramUsageError(Exception):  # TODO maybe remove or make more specific
    "Raised when dram is not used corectly, usally it means you are missing a step"

    pass


def get_ordered_uniques(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x) or pd.isna(x))]


def build_taxonomy_df(annotations_df: pd.DataFrame, groupby_column=DEFAULT_GROUPBY_COLUMN):
    cols = [groupby_column, "taxonomy"]
    if "Completeness" in annotations_df.columns:
        cols.append("Completeness")
    if "Contamination" in annotations_df.columns:
        cols.append("Contamination")
    tax_df = annotations_df[cols].drop_duplicates()
    tax_df.rename(columns={groupby_column: "genome"}, inplace=True)
    return tax_df


def build_tax_edge_df(
    tax_df,
):
    # regex = "".join([regex for regex in TAXONOMY_RANKS_REGEX.values()])  # regex that might be useful later
    # tree = tax_df["taxonomy"].str.extractall(regex)
    tree = pd.DataFrame(
        tax_df["taxonomy"].str.split(";").to_list(),
        columns=["domain", "phylum", "class", "order", "family", "genus", "species"],
        index=tax_df.index,
    )
    tax_df = tax_df.merge(tree, left_index=True, right_index=True, validate="1:1")

    tax_edge_df = (
        pd.concat(
            [
                tree[["domain", "phylum"]].rename(columns={"domain": "source", "phylum": "target"}),
                tree[["phylum", "class"]].rename(columns={"phylum": "source", "class": "target"}),
                tree[["class", "order"]].rename(columns={"class": "source", "order": "target"}),
                tree[["order", "family"]].rename(columns={"order": "source", "family": "target"}),
                tree[["family", "genus"]].rename(columns={"family": "source", "genus": "target"}),
                tree[["genus", "species"]].rename(columns={"genus": "source", "species": "target"}),
            ]
        )
        .drop_duplicates()
        .reset_index(drop=True)
    )

    return tax_edge_df, tax_df


def build_tree(edge_df, source_col: str = "source", target_col: str = "target", state: dict = None, id_cb=None):
    """
    Builds a tree structure from an edge DataFrame.

    Parameters:
    - edge_df (DataFrame): The edge DataFrame containing the source and target nodes.
    - source_col (str): The name of the column in edge_df that represents the source nodes. Default is "source".
    - target_col (str): The name of the column in edge_df that represents the target nodes. Default is "target".
    - state (dict): A dictionary representing the state of the tree nodes. Default is None.
    - id_cb (callable): A callback function that generates unique IDs for the tree nodes. Default is None.

    Returns:
    - tree_data (list): A list of dictionaries representing the tree structure.

    """

    def recurse_tree(source, id_cb, parent_id=None):
        """
        Recursively builds the tree structure.

        Parameters:
        - source: The current source node.
        - id_cb (callable): A callback function that generates unique IDs for the tree nodes.
        - parent_id: The ID of the parent node. Default is None.

        Returns:
        - js_ls (list): A list of dictionaries representing the tree structure.

        """
        target_nodes = edge_df.loc[edge_df[source_col] == source, target_col]
        js_ls = []
        for target in target_nodes:
            id_ = id_cb(source, target, parent_id)
            js_ls.append(
                {
                    "text": target,
                    "children": recurse_tree(target, parent_id=id_, id_cb=id_cb),
                    "state": state,
                    "id": id_,
                }
            )
        return js_ls

    state = state or {}
    roots = edge_df.loc[~edge_df[source_col].isin(edge_df[target_col]), source_col].unique()
    tree_data = [
        {"id": root, "text": root, "children": recurse_tree(root, parent_id=root, id_cb=id_cb), "state": state}
        for root in roots
    ]
    return tree_data


def build_tax_tree_selected_recurse(tax_tree, flat_tree=None):
    if flat_tree is None:
        flat_tree = []
    for node in tax_tree:
        flat_tree.append(node["id"])
        if "children" in node:
            build_tax_tree_selected_recurse(node["children"], flat_tree)
    return flat_tree
