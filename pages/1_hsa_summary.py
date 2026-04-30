"""
Page 1 — HSA Summary

High-level summary of HSA (HLA-Specific Antibody) distributions with demographic
variables in the ~4K patient cohort.

Left plot:  Recipient / HSA count histogram by a user-selected X variable, colored by subgroup.
            Supports both raw counts and within-X-group percentage stacking.
Right plot: MFI (Mean Fluorescence Intensity) distribution histogram, colored by subgroup.

Both plots support dynamic filtering by SEX, RACE, HLA_LOCUS, HLA_MATCH, MFI_GRP,
HSA_CNTs, and AGEGRP via a multi-select dropdown.
"""

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
import plotly.express as px
import os
import dash
from dash import html, dcc, callback, Output, Input
import dash_bootstrap_components as dbc
import warnings

# Suppress FutureWarnings from pandas / plotly
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.copy_on_write = True

dash.register_page(__name__)

# ── Data Loading ──────────────────────────────────────────────────────────────
default_path = "C:/Users/tzhang/DSA"
os.chdir(default_path)

df = pd.read_csv('dsa_data1.csv', sep=',')

# Standardise column names
df.rename(columns={
    'SEX_CDE': 'SEX',
    'broad_race_mapped': 'RACE',
    'RECIPIENT_AGEDEC': 'AGEGRP',
}, inplace=True)

# Clean MFI_GRP: replace NaN / 'NAN' with 0, cast to int
df.loc[(df['MFI_GRP'] == 'NAN') | (df['MFI_GRP'].isna()), 'MFI_GRP'] = 0
df['MFI_GRP'] = df['MFI_GRP'].astype(int)

# Clean AGEGRP: cast to int, cap at 90
df['AGEGRP'] = df['AGEGRP'].astype(int)
df.loc[df['AGEGRP'] > 90, 'AGEGRP'] = 90

# Keep only the four main race groups
df = df[df['RACE'].isin(['AFA', 'API', 'CAU', 'HIS'])]

# ── Widget Option Lists ──────────────────────────────────────────────────────
# Variables available for X-axis / subgroup radio buttons
group_category = ['SEX', 'RACE', 'HLA_LOCUS', 'HLA_MATCH', 'HSA_CNTs', 'MFI_GRP', 'AGEGRP']


def _build_filter_opts(series, col_name, cap=None):
    """Build 'COL=value' filter-option strings for a column.

    Parameters
    ----------
    series : pd.Series   – the column data
    col_name : str        – prefix used in the option label
    cap : int | None      – if set, values >= cap are collapsed to cap
    """
    def _fmt(x):
        if x != x:                       # NaN check
            return "None"
        if cap is not None and x >= cap:
            return f'{col_name}={cap}'
        return f'{col_name}={x}'
    return series.apply(_fmt).unique().tolist()


# Build combined filter dropdown options across all grouping variables
# HSA_CNTs values ≥ 10 are capped at 10 for grouping purposes
filter_category = (
    _build_filter_opts(df['SEX'], 'SEX')
    + _build_filter_opts(df['RACE'], 'RACE')
    + _build_filter_opts(df['HLA_LOCUS'], 'HLA_LOCUS')
    + _build_filter_opts(df['HLA_MATCH'], 'HLA_MATCH')
    + _build_filter_opts(df['MFI_GRP'], 'MFI_GRP')
    + _build_filter_opts(df['HSA_CNTs'], 'HSA_CNTs', cap=10)
    + _build_filter_opts(df['AGEGRP'], 'AGEGRP')
)


# ── Shared Helper: apply user-selected 'COL=value' filters ──────────────────
_INVALID = {None, "None", "NA", "na", "nan", ""}


def _apply_filters(temp_df, filter_chosen):
    """Apply 'COL=value' style filters to *temp_df*.

    For numeric columns the comparison uses ``>`` (greater-than);
    for string columns it uses ``==`` (equality).
    Handles both list and single-value inputs from the dropdown.
    """
    if filter_chosen is None:
        return temp_df

    items = filter_chosen if isinstance(filter_chosen, list) else [filter_chosen]
    for itm in items:
        if itm in _INVALID:
            continue
        col, val = itm.split('=', 1)
        if is_numeric_dtype(temp_df[col]):
            temp_df = temp_df.loc[temp_df[col] > int(val)]
        else:
            temp_df = temp_df.loc[temp_df[col] == val]
    return temp_df


# ── Shared Helper: set ordered Categorical dtype ─────────────────────────────
def _set_categorical(series, var_name):
    """Convert *series* to an ordered Categorical with sensible ordering.

    HLA_MATCH uses a fixed clinical order; SEX / RACE / HLA_LOCUS use
    alphabetical; everything else is sorted numerically.
    """
    s = series.astype(str)
    if var_name in ('SEX', 'RACE', 'HLA_LOCUS'):
        return pd.Categorical(s)
    if var_name == 'HLA_MATCH':
        return pd.Categorical(s, categories=['8/8', '7/8', '6/8'], ordered=True)
    cats = series.sort_values().unique().astype(str)
    return pd.Categorical(s, categories=cats, ordered=True)


# ── Page Layout ──────────────────────────────────────────────────────────────
WIDGET_STYLE = {
    'width': '20%', 'display': 'inline-block',
    'font-size': '12px', 'height': '20vh', 'padding': '0 10',
}
PLOT_STYLE = {
    'width': '49%', 'display': 'inline-block',
    'font-size': '12px', 'height': '50vh', 'padding': '0 10',
}

layout = dbc.Container([
    html.Hr(),

    # ── Control widgets row ──
    dbc.Row([
        # X-axis variable selector (left plot only)
        html.Div([
            dbc.Label("Variables for X axis of the plot1"),
            dbc.RadioItems(options=group_category, value=group_category[2], id='control_x_plot1'),
        ], style=WIDGET_STYLE),

        # Subgroup / colour variable selector (both plots)
        html.Div([
            dbc.Label("Variables for subgroup of both plots"),
            dbc.RadioItems(options=group_category, value=group_category[0], id='control_group'),
        ], style=WIDGET_STYLE),

        # Y-axis metric selector (raw count vs percentage)
        html.Div([
            dbc.Label("Variables for Y axis of both plots"),
            dbc.RadioItems(
                options=[
                    'Total_reciepient_count',
                    'Total_reciepient_count_percentage',
                    'Total_HSA_count',
                    'Total_HSA_count_percentage',
                ],
                value='Total_reciepient_count',
                id='control_count',
            ),
        ], style=WIDGET_STYLE),

        # Multi-select filter dropdown
        html.Div([
            dcc.Dropdown(id="control_filter", options=filter_category, multi=True),
        ], style={
            'width': '20%', 'color': 'blue', 'font-size': '12px',
            'height': '20vh', 'margin': '20px', 'padding': '0 10px',
        }),
    ]),

    # ── Plot row ──
    dbc.Row([
        # Left: HSA / recipient count histogram
        html.Div([dcc.Graph(id='HSA_count')], style=PLOT_STYLE),
        # Right: MFI distribution histogram
        html.Div([dcc.Graph(id='mfi_count')], style=PLOT_STYLE),
    ]),
])


# ── Callback: HSA / Recipient Count Histogram (left plot) ───────────────────
@callback(
    Output('HSA_count', 'figure'),
    Input('control_x_plot1', 'value'),
    Input('control_group', 'value'),
    Input('control_count', 'value'),
    Input('control_filter', 'value'),
)
def update_HSA_graph(x_plot1_chosen, subgroup_chosen, count_chosen, filter_chosen):
    """Build the left histogram.

    Depending on *count_chosen* the plot shows either:
    * raw counts  → grouped histogram (Total_reciepient_count / Total_HSA_count)
    * percentages → stacked bar normalised within each X-group
    """
    temp_df = df.copy()
    # Cap extreme values for cleaner visualisation
    temp_df.loc[temp_df['HSA_CNTs'] > 10, 'HSA_CNTs'] = 10
    temp_df.loc[temp_df['MFI_VALUE_NUM'] > 20000, 'MFI_VALUE_NUM'] = 20000
    temp_df['subgroup'] = temp_df[subgroup_chosen]

    # Apply user-selected filters
    temp_df = _apply_filters(temp_df, filter_chosen)

    # Deduplicate: one row per recipient for recipient-level counts,
    # otherwise one row per (recipient, subgroup) combination
    if count_chosen == 'Total_reciepient_count':
        temp_df_unique = temp_df.drop_duplicates(subset=['RECIPIENT_GUID'])
    else:
        temp_df_unique = temp_df.drop_duplicates(subset=['RECIPIENT_GUID', subgroup_chosen])

    # Drop rows with NaN in the axis / subgroup columns
    temp_df_unique = temp_df_unique[
        temp_df_unique[x_plot1_chosen].notna() & temp_df_unique['subgroup'].notna()
    ]

    # Apply ordered categorical types for consistent axis ordering
    temp_df_unique['subgroup'] = _set_categorical(temp_df_unique['subgroup'], subgroup_chosen)
    temp_df_unique[x_plot1_chosen] = _set_categorical(temp_df_unique[x_plot1_chosen], x_plot1_chosen)
    temp_df_unique = temp_df_unique.sort_values([x_plot1_chosen, 'subgroup'])

    # ── Build figure ──
    if count_chosen in ('Total_reciepient_count', 'Total_HSA_count'):
        # Raw-count histogram
        fig = px.histogram(
            temp_df_unique, x=x_plot1_chosen, color='subgroup',
            color_discrete_sequence=px.colors.qualitative.Alphabet,
            histfunc='sum', opacity=0.6,
        ).update_layout(font=dict(size=12), yaxis_title='Recipient counts')
        fig.update_yaxes(range=[0, 10000])
    else:
        # Percentage stacked bar (normalised within each X-group)
        def _normalizer(x):
            return x / x.sum()

        temp_df_unique['count'] = (
            temp_df_unique['HSA_CNTs']
            if count_chosen == 'Total_HSA_count_percentage'
            else 1
        )
        temp_df_unique['count'] = temp_df_unique.groupby(x_plot1_chosen)['count'].transform(_normalizer)
        agg = temp_df_unique.groupby([x_plot1_chosen, 'subgroup'])['count'].sum().reset_index()

        fig = px.bar(
            agg, x=x_plot1_chosen, y='count', color='subgroup',
            labels={'subgroup': str(subgroup_chosen)},
            color_discrete_sequence=px.colors.qualitative.Alphabet,
            barmode='stack',
        ).update_layout(
            yaxis=dict(tickformat='.0%', title='Recipient counts'),
            font=dict(size=12),
        )

    # Common styling
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', categoryorder='category ascending')
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black')
    fig.update_layout(legend_title_text=subgroup_chosen, height=450, plot_bgcolor='white')
    return fig


# ── Callback: MFI Distribution Histogram (right plot) ───────────────────────
@callback(
    Output('mfi_count', 'figure'),
    Input('control_group', 'value'),
    Input('control_count', 'value'),
    Input('control_filter', 'value'),
)
def update_mfi_graph(subgroup_chosen, count_chosen, filter_chosen):
    """Build the right histogram showing MFI value distributions.

    X axis: MFI_VALUE_NUM (capped at 20 000), 100 bins.
    Colour: user-selected subgroup variable.
    """
    temp_df = df.copy()
    temp_df.loc[temp_df['HSA_CNTs'] > 10, 'HSA_CNTs'] = 10
    temp_df.loc[temp_df['MFI_VALUE_NUM'] > 20000, 'MFI_VALUE_NUM'] = 20000
    temp_df['subgroup'] = temp_df[subgroup_chosen]

    # Apply user-selected filters
    temp_df = _apply_filters(temp_df, filter_chosen)

    # Deduplicate
    if count_chosen == 'Total_reciepient_count':
        temp_df_unique = temp_df.drop_duplicates(subset=['RECIPIENT_GUID'])
    else:
        temp_df_unique = temp_df.drop_duplicates(subset=['RECIPIENT_GUID', subgroup_chosen])

    temp_df_unique = temp_df_unique.sort_values('subgroup')

    fig = px.histogram(
        temp_df_unique, x='MFI_VALUE_NUM', color='subgroup',
        color_discrete_sequence=px.colors.qualitative.Alphabet,
        histfunc='avg', opacity=0.6, nbins=100,
    ).update_layout(font=dict(size=12), yaxis_title='MFI counts')

    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', range=[0, 20000])
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', range=[0, 200])
    fig.update_layout(height=450, plot_bgcolor='white')
    return fig
