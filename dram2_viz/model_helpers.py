import panel as pn


def format_chart_group(chart_group, title=""):
    if not isinstance(chart_group, list) or isinstance(chart_group, tuple):
        chart_group = [chart_group]
    return pn.Column(
        pn.pane.Markdown(f"## {title}", align="center"),
        pn.Row(*chart_group)
    )


def chart_group(chart, df, groupby_column, group_title, **chart_kwargs):
    """
    Group the data by the groupby_column and return a list of charts
    """
    charts = [
        chart(df_group, title=f"{group_title}: {group}", **chart_kwargs)
        for group, df_group in df.groupby(groupby_column)
    ]

    return format_chart_group(charts, title=group_title)
