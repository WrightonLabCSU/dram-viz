from src.rules import (
    evaluate_rules_on_anno,
    CompiledRules,
    Name,
    # Number,
    And,
    Or,
    Steps,
    Call,
)
import polars as pl
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"

BASIC_NEEDS = {
    "K00399",
    "K00400",
    "K00401",
    "K00402",
    "K03421",
    "K03422",
    "K02591",
    "K02586",
    "K02588",
}

COMPOUND_NEEDS = {
    "K00400",
    "K00399",
    "K03422",
    "K02591",
    "K13380",
    "K02588",
    "K00401",
    "K00402",
    "K00333",
    "K03421",
    "K02586",
    "K14126",
}

BASIC_PARSE = {
    "A": And(
        a=Or(
            a=Or(
                a=Or(
                    a=Or(
                        a=Or(
                            a=Name(name="K00399", db=None),
                            b=Name(name="K00400", db=None),
                        ),
                        b=Name(name="K00401", db=None),
                    ),
                    b=Name(name="K00402", db=None),
                ),
                b=Name(name="K03421", db=None),
            ),
            b=Name(name="K03422", db=None),
        ),
        b=Call(
            name="not",
            args=(
                Or(
                    a=Or(
                        a=Name(name="K02591", db=None),
                        b=Name(name="K02586", db=None),
                    ),
                    b=Name(name="K02588", db=None),
                ),
            ),
        ),
    )
}

COMPOUND_PARSE = {
    "A": And(
        a=Or(
            a=Or(
                a=Or(
                    a=Or(
                        a=Or(
                            a=Name(name="K00399", db=None),
                            b=Name(name="K00400", db=None),
                        ),
                        b=Name(name="K00401", db=None),
                    ),
                    b=Name(name="K00402", db=None),
                ),
                b=Name(name="K03421", db=None),
            ),
            b=Name(name="K03422", db=None),
        ),
        b=Call(
            name="not",
            args=(
                Or(
                    a=Or(
                        a=Name(name="K02591", db=None), b=Name(name="K02586", db=None)
                    ),
                    b=Name(name="K02588", db=None),
                ),
            ),
        ),
    ),
    "B": Or(
        a=Call(
            name="percent",
            args=(
                Name(name="50", db=None),
                Steps(
                    parts=(
                        Name(name="K00399", db=None),
                        And(
                            a=And(
                                a=And(
                                    a=Name(name="K00401", db=None),
                                    b=Name(name="K00402", db=None),
                                ),
                                b=Or(
                                    a=Name(name="K00333", db=None),
                                    b=Name(name="K03421", db=None),
                                ),
                            ),
                            b=Or(
                                a=Name(name="K03422", db=None),
                                b=Name(name="K13380", db=None),
                            ),
                        ),
                        Name(name="K00400", db=None),
                        Name(name="K00401", db=None),
                        Name(name="K00402", db=None),
                        Name(name="K03421", db=None),
                        Name(name="K03422", db=None),
                    )
                ),
            ),
        ),
        b=Or(
            a=Or(a=Name(name="K02591", db=None), b=Name(name="K02586", db=None)),
            b=Name(name="K02588", db=None),
        ),
    ),
    "C": And(
        a=Call(
            name="not",
            args=(
                Call(
                    name="column_contains",
                    args=(
                        Name(name="kegg_description", db=None),
                        Name(name="NitrateReducer", db=None),
                    ),
                ),
            ),
        ),
        b=Call(
            name="column_count_values",
            args=(
                Name(name="heme_regulatory_motif_count", db=None),
                Name(name="ge", db=None),
                Name(name="4", db=None),
                Name(name="ge", db=None),
                Name(name="3", db=None),
            ),
        ),
    ),
    "D": Call(
        name="percent",
        args=(
            Name(name="50", db=None),
            Steps(
                parts=(
                    Name(name="K00399", db=None),
                    Name(name="K14126", db=None),
                    Name(name="K02588", db=None),
                )
            ),
        ),
    ),
}


def test_basic_compiled_rules():
    compiled = CompiledRules.from_rules(rules_path=DATA_DIR / "basic_rules.tsv")
    assert compiled.needed_features == BASIC_NEEDS

    assert compiled.rules == BASIC_PARSE


def test_basic_rules_pass_annotation():

    df = evaluate_rules_on_anno(
        rules_path=DATA_DIR / "basic_rules.tsv",
        annotations_path=DATA_DIR / "basic_anno.tsv",
        sample_col="input_fasta",
    )

    true_df = pl.DataFrame(
        {
            "input_fasta": ["bin1", "bin2", "bin3", "bin4", "bin5", "bin6"],
            "A": [True, False, False, False, False, False],
        }
    )
    assert df.equals(true_df)


def test_compoun_compiled_rules():
    compiled = CompiledRules.from_rules(rules_path=DATA_DIR / "compound_rules.tsv")
    assert compiled.needed_features == COMPOUND_NEEDS

    assert compiled.rules == COMPOUND_PARSE


def test_compound_rules_pass_annotation():

    df = evaluate_rules_on_anno(
        rules_path=DATA_DIR / "compound_rules.tsv",
        annotations_path=DATA_DIR / "basic_anno.tsv",
        sample_col="input_fasta",
    )

    true_df = pl.DataFrame(
        {
            "input_fasta": ["bin1", "bin2", "bin3", "bin4", "bin5", "bin6"],
            "A": [True, False, False, False, False, False],
            "B": [False, False, False, False, False, True],
            "C": [False, False, False, False, True, False],
            "D": [True, False, False, False, False, False],
        }
    )
    assert df.equals(true_df)
