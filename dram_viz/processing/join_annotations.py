import pandas as pd


def join_dataframes(df1, df2, on_column, how="inner"):
    """
    Joins two pandas DataFrames on a specified column.

    Args:
        df1 (pd.DataFrame): The first DataFrame.
        df2 (pd.DataFrame): The second DataFrame.
        on_column (str): The column name to join on.
        how (str): Type of join - 'inner', 'outer', 'left', 'right'. Default is 'inner'.

    Returns:
        pd.DataFrame: The resulting DataFrame after the join.
    """
    if on_column not in df1.columns or on_column not in df2.columns:
        raise ValueError(f"The column '{on_column}' must be present in both DataFrames.")

    return pd.merge(df1, df2, on=on_column, how=how)
