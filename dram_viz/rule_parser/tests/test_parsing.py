from src.rules import (
    evaluate_rules_on_anno,
    CompiledRules,
    Name,
    PipeChain,
    And,
    Or,
    Steps,
    Call,
    String,
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
        parts=(
            Or(
                parts=(
                    Name(value="K00399", db=None),
                    Name(value="K00400", db=None),
                    Name(value="K00401", db=None),
                    Name(value="K00402", db=None),
                    Name(value="K03421", db=None),
                    Name(value="K03422", db=None),
                )
            ),
            Call(
                value="not",
                args=(
                    Or(
                        parts=(
                            Name(value="K02591", db=None),
                            Name(value="K02586", db=None),
                            Name(value="K02588", db=None),
                        )
                    ),
                ),
            ),
        )
    )
}

COMPOUND_PARSE = {
    "A": And(
        parts=(
            Or(
                parts=(
                    Name(value="K00399", db=None),
                    Name(value="K00400", db=None),
                    Name(value="K00401", db=None),
                    Name(value="K00402", db=None),
                    Name(value="K03421", db=None),
                    Name(value="K03422", db=None),
                )
            ),
            Call(
                value="not",
                args=(
                    Or(
                        parts=(
                            Name(value="K02591", db=None),
                            Name(value="K02586", db=None),
                            Name(value="K02588", db=None),
                        )
                    ),
                ),
            ),
        )
    ),
    "B": Or(
        parts=(
            Call(
                value="percent",
                args=(
                    Name(value="50", db=None),
                    Steps(
                        parts=(
                            Name(value="K00399", db=None),
                            And(
                                parts=(
                                    Name(value="K00401", db=None),
                                    Name(value="K00402", db=None),
                                    Or(
                                        parts=(
                                            Name(value="K00333", db=None),
                                            Name(value="K03421", db=None),
                                        )
                                    ),
                                    Or(
                                        parts=(
                                            Name(value="K03422", db=None),
                                            Name(value="K13380", db=None),
                                        )
                                    ),
                                )
                            ),
                            Name(value="K00400", db=None),
                            Name(value="K00401", db=None),
                            Name(value="K00402", db=None),
                            Name(value="K03421", db=None),
                            Name(value="K03422", db=None),
                        )
                    ),
                ),
            ),
            Or(
                parts=(
                    Name(value="K02591", db=None),
                    Name(value="K02586", db=None),
                    Name(value="K02588", db=None),
                )
            ),
        )
    ),
    "C": PipeChain(
        calls=(
            Call(
                value="not",
                args=(
                    Call(
                        value="filter_contains",
                        args=(
                            Name(value="kegg_description", db=None),
                            String(value="nitrate reductase"),
                        ),
                    ),
                ),
            ),
            Call(
                value="column_count_values",
                args=(
                    Name(value="heme_regulatory_motif_count", db=None),
                    Name(value="ge", db=None),
                    Name(value="2", db=None),
                    Name(value="ge", db=None),
                    Name(value="2", db=None),
                ),
            ),
        )
    ),
    "D": Call(
        value="percent",
        args=(
            Name(value="50", db=None),
            Steps(
                parts=(
                    Name(value="K00399", db=None),
                    Name(value="K14126", db=None),
                    Name(value="K02588", db=None),
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
            "C": [True, False, False, False, False, False],
            "D": [True, False, False, False, False, False],
        }
    )
    assert df.equals(true_df)
