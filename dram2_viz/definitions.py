from __future__ import annotations

import re
from pathlib import Path

DEFAULT_GROUPBY_COLUMN = "sample"
BACKUP_GROUPBY_COLUMN = "fasta"

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
DBSETS_COL = "db_id_sets"
MODULE_STEPS_FORM_TAG = "module_step_form"
FUNCTION_HEATMAP_FORM_TAG = "function_heatmap_form"
ETC_MODULE_DF_TAG = "etc_module_database"
FILES_NAMES: dict[str, Path] = {
    MODULE_STEPS_FORM_TAG: Path(__file__).parent.resolve() / "data/module_step_form.tsv",
    FUNCTION_HEATMAP_FORM_TAG: Path(__file__).parent.resolve() / "data/function_heatmap_form.tsv",
    ETC_MODULE_DF_TAG: Path(__file__).parent.resolve() / "data/etc_module_database.tsv",
}
ID_FUNCTION_DICT = {
    'kegg_genes_id': lambda x: [x],
    'ko_id': lambda x: [j for j in x.split(',')],
    'kegg_id': lambda x: [j for j in x.split(',')],
    'kegg_hit': lambda x: [i[1:-1] for i in
                           re.findall(r'\[EC:\d*.\d*.\d*.\d*\]', x)],
    'peptidase_family': lambda x: [j for j in x.split(';')],
    'cazy_best_hit': lambda x: [x.split('_')[0]],
    'pfam_hits': lambda x: [j[1:-1].split('.')[0]
                            for j in re.findall(r'\[PF\d\d\d\d\d.\d*\]', x)],
    'camper_id': lambda x: [x],
    'fegenie_id': lambda x: [x],
    'sulfur_id': lambda x: [x],
    'methyl_id': lambda x: [i.split(' ')[0].strip() for i in x.split(',')]
}
KO_REGEX = r"^K\d\d\d\d\d$"
ETC_COVERAGE_COLUMNS = [
    "module_id",
    "module_name",
    "complex",
    "genome",
    "path_length",
    "path_length_coverage",
    "percent_coverage",
    "genes",
    "missing_genes",
    "complex_module_name",
]
TAXONOMY_LEVELS = ["d", "p", "c", "o", "f", "g", "s"]
