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
