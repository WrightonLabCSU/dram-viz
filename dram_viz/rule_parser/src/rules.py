from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple, List, Optional, Set
import operator
import os
from pathlib import Path

import numpy as np
import polars as pl
from lark import Lark, Transformer

OP_TO_EXPR = {
    "gt": operator.gt,
    "ge": operator.ge,
    "lt": operator.lt,
    "le": operator.le,
    "eq": operator.eq,
    "ne": operator.ne,
}
ALLOWED_CMPOPS = set(OP_TO_EXPR.keys())


class RuleError(Exception):
    pass


ID_EXPR_DICT = {
    "camper_id": pl.col("camper_id").cast(pl.Utf8).cast(pl.List(pl.Utf8)),
    "fegenie_id": pl.col("fegenie_id").cast(pl.Utf8).cast(pl.List(pl.Utf8)),
    "sulfur_id": pl.col("sulfur_id").cast(pl.Utf8).cast(pl.List(pl.Utf8)),
    "kegg_genes_id": pl.col("kegg_genes_id").cast(pl.Utf8).cast(pl.List(pl.Utf8)),
    "kofam_genes_id": pl.col("kofam_genes_id").cast(pl.Utf8).cast(pl.List(pl.Utf8)),
    "ko_id": pl.col("ko_id").str.split(","),
    "kegg_id": pl.col("kegg_id").str.split("/"),
    "kofam_id": pl.col("kofam_id").str.split("/"),
    "peptidase_family": pl.col("peptidase_family").str.split(";"),
    "merops_family": pl.col("merops_family").str.split(";"),
    "camper_EC": pl.col("camper_EC").str.extract_all(r"\[EC:(\d+\.\d+\.\d+\.\d+)\]"),
    "kegg_EC": pl.col("kegg_EC").str.extract_all(r"\[EC:(\d+\.\d+\.\d+\.\d+)\]"),
    "kofam_EC": pl.col("kofam_EC").str.extract_all(r"\[EC:(\d+\.\d+\.\d+\.\d+)\]"),
    "pfam_hits": pl.col("pfam_hits").str.extract_all(r"\[(PF\d{5})\.\d+\]"),
    "pfam_id": pl.col("pfam_id").str.extract_all(r"\[(PF\d{5})\.\d+\]"),
    "cazy_hits": (
        pl.col("cazy_hits").cast(pl.Utf8)
        # capture just the EC number part inside "(EC ...)"
        .str.extract_all(r"\(EC\s+([\d\.\-]+)\)")
        # prefix each captured number with "EC:"
        .list.eval(pl.concat_str([pl.lit("EC:"), pl.element()]))
    ),
    "cazy_subfam_ec": (
        pl.col("cazy_subfam_ec")
        .cast(pl.Utf8)
        .str.extract_all(r"([\d\.\-]+)")
        .list.eval(pl.concat_str([pl.lit("EC:"), pl.element()]))
    ),
    "dbcan_EC": (
        pl.col("dbcan_EC")
        .cast(pl.Utf8)
        .str.extract_all(r"([\d\.\-]+)")
        .list.eval(pl.concat_str([pl.lit("EC:"), pl.element()]))
    ),
    "cazy_best_hit": pl.col("cazy_best_hit").cast(pl.Utf8).str.split("_"),
    "dbcan_id": pl.col("dbcan_id").cast(pl.Utf8).str.split("_"),
    "methyl_id": (
        pl.col("methyl_id")
        .cast(pl.Utf8)
        .str.split(",")
        .list.eval(pl.element().str.strip_chars().str.split(" ").list.first())
    ),
}

CALL_FUNCTIONS = {
    "not",
    "percent",
    "at_least",
    "column_count_values",
    "column_sum_values",
    "column_contains",
}

# AST node definitions and Grammar definition


@dataclass(frozen=True)
class Expr: ...


@dataclass(frozen=True)
class Name(Expr):
    name: str
    db: str | None = None


@dataclass(frozen=True)
class Number(Expr):
    value: float


@dataclass(frozen=True)
class And(Expr):
    a: Expr
    b: Expr


@dataclass(frozen=True)
class Or(Expr):
    a: Expr
    b: Expr


@dataclass(frozen=True)
class Steps(Expr):
    parts: Tuple[Expr, ...]  # list-of-parts (coverage/count units)


@dataclass(frozen=True)
class Call(Expr):
    name: str
    args: Tuple[Expr, ...]

    def __post_init__(self):  # validation
        if self.name not in CALL_FUNCTIONS:
            raise RuleError(f"Unknown function: {self.name}")
        n_args = len(self.args)
        match self.name:
            case "not":
                if n_args != 1:
                    raise RuleError("not(...) expects 1 arg")
            case "percent" | "at_least":
                if n_args != 2:
                    raise RuleError(f"{self.name}(n, Group) expects 2 args")
            case "column_count_values":
                if n_args != 5:
                    raise RuleError(
                        "column_count_values(column, val_op, val_threshold, count_op, count_thrreshold)"
                        " expects 5 args"
                    )
                for op in (self.args[1], self.args[3]):
                    op = _as_str(op)
                    if op not in ALLOWED_CMPOPS:
                        raise RuleError(
                            f"column_count_values ops must be one of {sorted(ALLOWED_CMPOPS)}; got {op}"
                        )
            case "column_sum_values":
                if n_args != 3:
                    raise RuleError(
                        "column_sum_values(column, op, threshold) expects 3 args"
                    )
                op = _as_str(self.args[1])
                if op not in ALLOWED_CMPOPS:
                    raise RuleError(
                        f"column_sum_values op must be one of {sorted(ALLOWED_CMPOPS)}; got {op}"
                    )
            case "column_contains":
                if n_args != 2:
                    raise RuleError("column_contains(column, value) expects 2 args")


def _as_str(e: Expr) -> str:
    try:
        return e.name
    except AttributeError:
        raise RuleError(f"Could not convert node to string, got {e}")


def _as_int(e: Expr) -> int:
    if isinstance(e, Number) and float(e.value).is_integer():
        return int(e.value)
    if isinstance(e, Name) and e.name.isdigit():
        return int(e.name)
    raise RuleError(f"Expected integer literal, got {e}")


def _as_float(e: Expr) -> float:
    if isinstance(e, Number):
        return float(e.value)
    try:
        return float(e.name)
    except (ValueError, AttributeError):
        raise RuleError(
            f"Could not convert node to gloat. Expected numeric literal, got {e}"
        )


class ASTTransformer(Transformer):
    def simple_name(self, items):
        return Name(name=str(items[0]), db=None)

    def qualified_name(self, items):
        db, name = items
        return Name(name=str(name), db=str(db))

    def number(self, items):
        return Number(float(str(items[0])))

    def and_(self, items):
        return And(items[0], items[1])

    def or_(self, items):
        return Or(items[0], items[1])

    def group(self, items):
        # square brackets are just grouping, not a node
        return items[0]

    def step_(self, items):
        return Steps(tuple(items))

    def call(self, items):
        name = str(items[0])
        args = tuple(items[1:])
        return Call(name, args)


@dataclass(frozen=True)
class CompiledRules:
    rules: Dict[str, Expr]  # top-level rules with macros expanded
    needed_features: Set[str]

    @classmethod
    def from_rules(cls, *args, **kwargs) -> CompiledRules:
        definitions, rules = load_rules(*args, **kwargs)

        # Expand macros inside definitions too (so GroupName -> Steps(...) can nest)
        defs_expanded = {
            k: expand_macros(v, definitions) for k, v in definitions.items()
        }

        # Expand rules using expanded defs
        # we need to hit again in case defs is empty (no parent col)
        # and we still need to add needed features from rules
        needed_features: Set[str] = set()
        rules_expanded = {
            k: expand_macros(v, defs_expanded, needed_features=needed_features)
            for k, v in rules.items()
        }

        return cls(rules=rules_expanded, needed_features=needed_features)


def load_rules(
    rules_path: str = None,
    rules: pl.LazyFrame = None,
    label_col: str = "name",
    parent_col: str = "parent",
    rules_col: str = "child",
) -> Tuple[Dict[str, Expr], Dict[str, Expr]]:
    """
    Assumes TSV has columns at least: name, parent, child
    Convention used here:
      - if `name` is non-empty: this is an OUTPUT RULE whose expression is in `child`
      - if `name` is empty and `parent` is non-empty: this is a DEFINITION macro: parent := child
    If your file uses a slightly different convention, adjust this function only.
    """
    assert (rules_path is not None) != (
        rules is not None
    ), "Either rules_path or rules DataFrame must be provided, but not both."
    if rules_path:
        lf = pl.scan_csv(
            rules_path, separator="\t", infer_schema_length=None
        ).fill_null("")
    else:
        lf = rules.fill_null("")

    # Normalize columns
    cols = set(lf.collect_schema().names())
    required = {label_col, rules_col}
    if not required.issubset(cols):
        raise ValueError(
            f"rules TSV missing required columns {required}. Found: {lf.columns}"
        )

    has_parent_col = parent_col and (parent_col in cols)

    with open(Path(__file__).parent.absolute() / "rules.lark") as f:
        parser = Lark(f, parser="lalr", transformer=ASTTransformer())

    lf = lf.with_columns(
        [
            pl.col(rules_col)
            .str.strip_chars()
            .map_elements(parser.parse, return_dtype=pl.Object)
        ]
    )

    lf = lf.with_columns(pl.col(pl.String).replace("", None))

    rules = {
        a: b
        for a, b in lf.filter(~pl.col(label_col).is_null())
        .select([pl.col(label_col), pl.col(rules_col)])
        .collect()
        .iter_rows()
    }

    if has_parent_col:
        definitions = {
            a: b
            for a, b in lf.filter(~pl.col(parent_col).is_null())
            .select([pl.col(parent_col), pl.col(rules_col)])
            .collect()
            .iter_rows()
        }
    else:
        definitions = {}

    return definitions, rules


def expand_macros(
    expr: Expr, definitions: Dict[str, Expr], needed_features: Set[str] = None
) -> Expr:
    """Expand recursively macros in expr using definitions"""
    memo: Dict[Expr, Expr] = {}
    stack: list[str] = []
    skip_needed = False
    if needed_features is None:
        skip_needed = True

    def recurse(e: Expr, add_name_to_needed: bool = True) -> Expr:
        if e in memo:
            return memo[e]

        if isinstance(e, Name):
            name = e.name
            if name in definitions:
                if name in stack:
                    raise RuleError(f"Cycle detected: {' -> '.join(stack + [name])}")
                stack.append(name)
                out = recurse(definitions[name])
                stack.pop()
            else:
                out = e
                # This will slightly
                if add_name_to_needed and not skip_needed:
                    needed_features.add(name)
        elif isinstance(e, Number):
            out = e
        elif isinstance(e, And):
            out = And(recurse(e.a), recurse(e.b))
        elif isinstance(e, Or):
            out = Or(recurse(e.a), recurse(e.b))
        elif isinstance(e, Steps):
            out = Steps(tuple(recurse(p) for p in e.parts))
        elif isinstance(e, Call):
            fn, args = e.name, e.args
            call_rec = []
            # Determine which args should have their Names added to needed_features
            # since some call args are just numbers or ops
            for i, arg in enumerate(args):
                if fn in {"percent", "at_least"}:
                    add_name = 1 == i
                elif fn == "not":
                    add_name = 0 == i
                elif fn in {
                    "column_count_values",
                    "column_sum_values",
                    "column_contains",
                }:
                    add_name = False
                # unknown function; be conservative. Worse case is we count extra needed features
                # such as numbers or ops, which is not ideal, but quite unlikely to be
                # a gene id. So almost always harmless.
                else:
                    add_name = True
                call_rec.append(recurse(arg, add_name_to_needed=add_name))

            out = Call(e.name, tuple(call_rec))
        else:
            raise TypeError(e)

        memo[e] = out
        return out

    return recurse(expr)


def build_present_map(
    lf: pl.DataFrame,
    sample_col: str,
    besthit_cols: List[str],
    needed_features: Set[str],
) -> Tuple[List[str], Dict[str, np.ndarray]]:
    """Build present_map of needed gene_ids from annotations DataFrame"""
    besthit_cols = [col for col in besthit_cols if col in lf.columns]
    for col in besthit_cols:
        lf = lf.with_columns(ID_EXPR_DICT[col].alias(col)).explode(col)
    lf = lf.select([sample_col] + besthit_cols)

    # unpivot to long (sample, hit)
    hit_col = "hit"
    lf = (
        lf.unpivot(
            index=sample_col,
            on=besthit_cols,
            variable_name="db",
            value_name=hit_col,
        )
        .drop("db")
        .drop_nulls()
    )

    samples = lf.select(sample_col).unique().sort(sample_col).to_series().to_list()

    df: pl.DataFrame = lf.select([sample_col, hit_col]).filter(
        pl.col(hit_col).is_in(list(needed_features))
    )

    sample_index = {s: i for i, s in enumerate(samples)}
    n = len(samples)
    present_map: Dict[str, np.ndarray] = {
        f: np.zeros(n, dtype=bool) for f in needed_features
    }

    # Group by hit_id and set booleans
    for hit_id, sub in df.group_by(hit_col, maintain_order=False):
        hit = hit_id[0]
        arr = present_map.get(hit)
        if arr is None:
            continue
        sub_s = sub.select(sample_col).unique().to_series().to_list()
        for s in sub_s:
            arr[sample_index[s]] = True

    return samples, present_map


class Evaluator:
    def __init__(
        self,
        samples: List[str],
        present_map: Dict[str, np.ndarray],
        sample_col: str,
        annotations: Optional[pl.DataFrame] = None,
    ):
        self.samples = samples
        self.present_map = present_map
        self.annotations = annotations
        self.sample_col = sample_col
        self._memo: Dict[Expr, np.ndarray] = {}
        self._memo_list: Dict[Expr, np.ndarray] = {}
        self._order_df = pl.DataFrame(
            {
                self.sample_col: self.samples,
                "_order": range(len(self.samples)),
            }
        )

        self._all_false = np.zeros(len(samples), dtype=bool)

    def eval_bool(self, expr: Expr) -> np.ndarray:
        if expr in self._memo:
            return self._memo[expr]

        if isinstance(expr, Name):
            out = self.present_map.get(expr.name, self._all_false)
            self._memo[expr] = out
            return out

        if isinstance(expr, And):
            out = np.logical_and(self.eval_bool(expr.a), self.eval_bool(expr.b))
            self._memo[expr] = out
            return out

        if isinstance(expr, Or):
            out = np.logical_or(self.eval_bool(expr.a), self.eval_bool(expr.b))
            self._memo[expr] = out
            return out

        if isinstance(expr, Call):
            out = self.eval_call(expr)
            self._memo[expr] = out
            return out

        if isinstance(expr, Number):
            # numeric alone isn't boolean; treat as error
            raise RuleError(f"Literal used where boolean expected: {expr}")

        if isinstance(expr, Steps):
            raise RuleError(
                f"Step expressions (expression seperated by commas)"
                " logic cannot be evaluted on their own. They must be used"
                " in a supporting function such as `percent`. Offending"
                " rule: {expr}"
            )

        raise RuleError(
            "Something failed to parse properly."
            f" Tried to evaluate truth value of {expr}"
        )

    def eval_cycle(self, expr: Steps) -> np.ndarray:
        if expr in self._memo_list:
            return self._memo_list[expr]

        parts = [self.eval_bool(p) for p in expr.parts]
        mat = (
            np.stack(parts, axis=1)
            if parts
            else np.zeros((len(self.samples), 0), dtype=bool)
        )
        self._memo_list[expr] = mat
        return mat

    def eval_call(self, call: Call) -> np.ndarray:
        fn = call.name
        args = call.args

        match fn:
            case "not":
                return self.not_(self.eval_bool(args[0]))
            case "percent":
                return self.percent(_as_int(args[0]), self.eval_cycle(args[1]))
            case "at_least":
                return self.at_least(_as_int(args[0]), self.eval_cycle(args[1]))
            case "column_count_values":
                return self.column_count_values(
                    col=_as_str(args[0]),
                    val_op=_as_str(args[1]),
                    val_thr=_as_float(args[2]),
                    count_op=_as_str(args[3]),
                    count_thr=_as_float(args[4]),
                )
            case "column_sum_values":
                return self.column_count_values(
                    col=_as_str(args[0]), op=_as_str(args[1]), thr=_as_float(args[2])
                )
            case "column_contains":
                return self.column_contains(col=_as_str(args[0]), val=_as_str(args[1]))
            case _:
                raise RuleError(f"Unable to parse function in rules. Function: {fn}")

    def _sort_df_to_ordered_df(self, df):
        return (
            df.join(self._order_df, on=self.sample_col, how="left")
            .sort("_order")
            .drop("_order")
        )

    # Call functions
    @staticmethod
    def not_(x: np.ndarray) -> np.ndarray:
        return np.logical_not(x)

    @staticmethod
    def percent(n: int, x: np.ndarray) -> np.ndarray:
        thr = n / 100.0
        if x.shape[1] == 0:
            return np.zeros(x.shape[0], dtype=bool)
        cov = x.mean(axis=1)
        return cov >= thr

    @staticmethod
    def at_least(k: int, x: np.ndarray) -> np.ndarray:
        return x.sum(axis=1) >= k

    def column_count_values(
        self, col: str, val_op: str, val_thr: float, count_op: str, count_thr: float
    ) -> np.ndarray:
        if col not in self.annotations.columns:
            raise RuleError(f"Missing column '{col}' for column_count_values()")
        if val_op not in OP_TO_EXPR:
            raise ValueError(
                f"Unsupported value operation={val_op!r} to compare values for column {col}. Use one of {sorted(OP_TO_EXPR)}"
            )
        val_cmp_fn = OP_TO_EXPR[val_op]
        if count_op not in OP_TO_EXPR:
            raise ValueError(
                f"Unsupported count operation={count_op!r} to compare value counts for column {col}. Use one of {sorted(OP_TO_EXPR)}"
            )
        count_cmp_fn = OP_TO_EXPR[count_op]

        df = (
            self.annotations.group_by(self.sample_col)
            # first do the val cmp fn on the column values to get a bool of which rows
            # pass the value threshold (e.g., per row: col_val >= val_thr)
            # Then sum those booleans to get a count of how many rows per sample pass
            .agg(
                count_cmp_fn(
                    val_cmp_fn(pl.col(col).cast(float), val_thr).sum(), count_thr
                )
            )
        )
        df = self._sort_df_to_ordered_df(df)
        return df.select(pl.col(col)).to_series().to_numpy()

    def column_sum_values(self, col: str, op: str, thr: float) -> np.ndarray:
        if col not in self.annotations.columns:
            raise RuleError(f"Missing column '{col}' for column_sum_values()")
        try:
            cmp_fn = OP_TO_EXPR[op]
        except KeyError:
            raise ValueError(f"Unsupported op={op!r}. Use one of {sorted(OP_TO_EXPR)}")

        df = self.annotations.group_by(self.sample_col).agg(
            cmp_fn(pl.col(col).cast(float).sum(), thr)
        )
        df = self._sort_df_to_ordered_df(df)
        return df.select(pl.col(col)).to_series().to_numpy()

    def column_contains(self, col: str, val: str) -> np.ndarray:
        if col not in self.annotations.columns:
            raise RuleError(f"Missing column '{col}' for column_contains()")

        df = self.annotations.group_by(self.sample_col).agg(
            pl.col(col).str.contains(val).any()
        )

        df = self._sort_df_to_ordered_df(df)
        return df.select(pl.col(col)).to_series().to_numpy()


def evaluate_rules(
    compiled: CompiledRules,
    samples: List[str],
    present_map: Dict[str, np.ndarray],
    annotations: Optional[pl.DataFrame] = None,
    sample_col: Optional[str] = None,
    transpose: bool = False,
) -> pl.DataFrame:
    ev = Evaluator(
        samples=samples,
        present_map=present_map,
        sample_col=sample_col,
        annotations=annotations,
    )

    rule_names = sorted(compiled.rules.keys())
    data = {}
    for rn in rule_names:
        data[rn] = ev.eval_bool(compiled.rules[rn])

    df = pl.DataFrame(data)
    df = df.with_columns(pl.Series(sample_col, samples)).select(
        [sample_col] + rule_names
    )

    if transpose:
        # rules x samples: make sample columns
        df_t = df.drop("sample").transpose(include_header=True)
        # df_t: first column is "column" (former rule name), others are sample-indexed columns
        df_t = df_t.rename({"column": "rule"})
        # add sample names as headers
        # polars transpose produces string column names "column_0", "column_1"...; rename to sample ids
        rename_map = {}
        for i, s in enumerate(samples):
            rename_map[f"column_{i}"] = s
        df_t = df_t.rename(rename_map)
        return df_t

    return df


def evaluate_rules_on_anno(
    annotations_path: os.PathLike = None,
    annotations: pl.DataFrame = None,
    sample_col: str = "input_fasta",
    *args,
    **kwargs,
):
    compiled = CompiledRules.from_rules(*args, **kwargs)
    print(f"Need {len(compiled.needed_features)} features")
    assert (annotations_path is not None) != (
        annotations is not None
    ), "Exactly one of annotations_path or annotations DataFrame must be provided."
    if annotations_path:
        try:
            annotations = pl.read_csv(
                annotations_path, separator="\t", infer_schema_length=10_000
            )
        except Exception as e:
            annotations = pl.read_csv(
                annotations_path, separator="\t", infer_schema_length=None
            )

    samples, present_map = build_present_map(
        annotations,
        sample_col=sample_col,
        besthit_cols=list(ID_EXPR_DICT.keys()),
        needed_features=compiled.needed_features,
    )

    df = evaluate_rules(
        compiled,
        samples,
        present_map,
        annotations=annotations,
        sample_col=sample_col,
        transpose=False,
    )
    return df
