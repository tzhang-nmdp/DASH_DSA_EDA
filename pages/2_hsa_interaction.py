"""
Page 2 — HSA Interaction Heatmaps

Detailed 4-way interaction analysis of demographic variables on HSA distribution.

Users select 4 categorical variables to map onto a grid of heatmap subplots:
  - Upper columns & Outer rows → define the subplot grid
  - Inner rows & Lower columns → define the axes within each heatmap cell
  - Optional normalization (row-wise, column-wise, or none)
  - MFI threshold range slider to re-bin MFI groups dynamically
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
import dash
from dash import html, dcc, callback, Output, Input
import dash_bootstrap_components as dbc
import warnings

pd.options.mode.copy_on_write = True
warnings.simplefilter(action='ignore', category=FutureWarning)

dash.register_page(__name__)

# ── Data Loading ──────────────────────────────────────────────────────────────
default_path = "C:/Users/tzhang/DSA"
os.chdir(default_path)

# Pre-aggregated DSA count table with demographic cross-tabulations
df = pd.read_csv('dsa_data2.csv', sep=',', header=0)

# ── Ordered category definitions for consistent axis ordering ────────────────
CATEGORY_ORDERS = {
    'HLA_MATCH': ['HLA_MATCH=8/8', 'HLA_MATCH=7/8', 'HLA_MATCH=6/8'],
    'HSA_CNTs': [f'HSA_CNTs={i}' for i in range(1, 11)],
    'RACE': ['RACE=AFA', 'RACE=API', 'RACE=CAU', 'RACE=HIS', 'RACE=OTH'],
    'AGR_GRP': [f'AGR_GRP={d}' for d in range(0, 100, 10)],
    'HLA_LOCUS': [
        'HLA_LOCUS=A', 'HLA_LOCUS=B', 'HLA_LOCUS=C',
        'HLA_LOCUS=DRB1', 'HLA_LOCUS=DRB3', 'HLA_LOCUS=DRB4', 'HLA_LOCUS=DRB5',
        'HLA_LOCUS=DPA1', 'HLA_LOCUS=DPB1', 'HLA_LOCUS=DQA1', 'HLA_LOCUS=DQB1',
    ],
}

# Dropdown options: all columns except the last 3 (count + metadata columns)
dropdown_options = [{"label": x, "value": x} for x in df.columns[:-3]]

# ── Page Layout ──────────────────────────────────────────────────────────────
WIDGET_STYLE = {
    'width': '25%', 'display': 'inline-block',
    'font-size': '12px', 'height': '10vh', 'padding': '0 5',
}

layout = dbc.Container([
    dbc.Row([
        # Upper column variable selector
        html.Div([
            dbc.Label("Interaction variables at upper columns"),
            dcc.Dropdown(id="upper_col", options=dropdown_options, value="HSA_CNTs", clearable=False, style={"width": "100%"}),
        ], style=WIDGET_STYLE),

        # Outer row variable selector
        html.Div([
            dbc.Label("Interaction variables at outer rows"),
            dcc.Dropdown(id="outter_row", options=dropdown_options, value="SEX", clearable=False, style={"width": "100%"}),
        ], style=WIDGET_STYLE),

        # Inner row variable selector
        html.Div([
            dbc.Label("Interaction variables at inner rows"),
            dcc.Dropdown(id="inner_row", options=dropdown_options, value="RACE", clearable=False, style={"width": "100%"}),
        ], style=WIDGET_STYLE),

        # Lower column variable selector
        html.Div([
            dbc.Label("Interaction variables at lower columns"),
            dcc.Dropdown(id="lower_col", options=dropdown_options, value="HLA_LOCUS", clearable=False, style={"width": "100%"}),
        ], style=WIDGET_STYLE),

        # Normalization toggle (row-wise, column-wise, or none)
        html.Div([
            dbc.Label("Normalization axis"),
            dbc.RadioItems(
                id='normalization_factors',
                options=[{"label": x, "value": x} for x in ['X-axis', 'Y-axis', 'No']],
                value='No', inline=True,
            ),
        ], style={'width': '20%', 'display': 'inline-block', 'font-size': '12px', 'height': '10vh', 'padding': '0 5'}),

        # MFI threshold range slider for dynamic MFI group binning
        html.Div([
            dbc.Label("MFI cutoffs"),
            dcc.RangeSlider(
                0, 20000, value=[500, 20000],
                marks={v: str(v) for v in [0, 100, 500, 1000, 2000, 5000, 10000, 20000]},
                id='mfi_cutoff',
            ),
        ], style={'width': '80%', 'display': 'inline-block', 'font-size': '12px', 'height': '10vh', 'padding': '0 5'}),
    ]),

    # Heatmap output
    dbc.Row([
        html.Div([
            dcc.Graph(id="dsa-heatmaps-graph"),
        ], style={'width': '100%', 'display': 'inline-block', 'font-size': '12px', 'height': '50vh', 'padding': '0 56'}),
    ]),
])


# ── Helpers ──────────────────────────────────────────────────────────────────
def _is_selected(var_name, *selected_vars):
    """Return True if *var_name* is among the 4 selected interaction variables."""
    return var_name in selected_vars


def _apply_category_order(data, var_name, selected_vars):
    """Apply an ordered Categorical type to *var_name* if it is selected and
    has a predefined ordering in CATEGORY_ORDERS."""
    if not _is_selected(var_name, *selected_vars):
        return
    if var_name in CATEGORY_ORDERS:
        data[var_name] = pd.Categorical(
            data[var_name], categories=CATEGORY_ORDERS[var_name], ordered=True,
        )


# ── Callback: Generate the 4-way interaction heatmap grid ───────────────────
@callback(
    Output("dsa-heatmaps-graph", "figure"),
    Input("upper_col", "value"),
    Input("outter_row", "value"),
    Input("inner_row", "value"),
    Input("lower_col", "value"),
    Input("normalization_factors", "value"),
    Input("mfi_cutoff", "value"),
)
def dsa_filter_heatmap(upper_col_var, outter_row_var, inner_row_var,
                       lower_col_var, normalization_decision, mfi_threshold):
    """Build a grid of heatmap subplots showing DSA count interactions.

    The subplot grid is defined by upper_col_var (columns) × outter_row_var (rows).
    Within each subplot cell, inner_row_var defines Y and lower_col_var defines X.
    """
    selected = (upper_col_var, outter_row_var, inner_row_var, lower_col_var)

    # Select only the 4 interaction columns + the count column (col index 8)
    data = df[[upper_col_var, inner_row_var, outter_row_var, lower_col_var, df.columns[8]]].copy()

    # ── Dynamic MFI group re-binning based on slider thresholds ──
    if _is_selected('MFI_GRP', *selected):
        data = data[data['MFI_GRP'] != 'NAN']
        mfi_l, mfi_h = int(mfi_threshold[0]), int(mfi_threshold[1])
        data['MFI_GRP'] = data['MFI_GRP'].apply(
            lambda x: (
                f'MFI=<{mfi_l}' if 0 < int(x) <= mfi_l
                else (f'MFI={mfi_l}~{mfi_h}' if mfi_l < int(x) <= mfi_h
                      else (f'MFI=>{mfi_h}' if int(x) > mfi_h else 0))
            )
        )
        data['MFI_GRP'] = pd.Categorical(
            data['MFI_GRP'],
            categories=[f'MFI=<{mfi_l}', f'MFI={mfi_l}~{mfi_h}', f'MFI=>{mfi_h}'],
            ordered=True,
        )

    # ── Apply ordered categories for consistent axis ordering ──
    for var in CATEGORY_ORDERS:
        _apply_category_order(data, var, selected)

    # ── Aggregate counts by all 4 variables ──
    data = (
        data
        .groupby([upper_col_var, inner_row_var, outter_row_var, lower_col_var])[df.columns[8]]
        .sum()
        .reset_index()
        .sort_values([upper_col_var, inner_row_var, outter_row_var, lower_col_var])
    )

    # ── Determine subplot grid dimensions ──
    upper_col_idx = len(data.iloc[:, 0].unique())   # number of subplot columns
    inner_row_idx = len(data.iloc[:, 1].unique())   # Y-axis levels inside each cell
    outter_row_idx = len(data.iloc[:, 2].unique())  # number of subplot rows
    lower_col_idx = len(data.iloc[:, 3].unique())   # X-axis levels inside each cell

    # ── Build subplot grid ──
    fig = make_subplots(
        rows=outter_row_idx, cols=upper_col_idx,
        shared_xaxes=True, shared_yaxes=True,
        horizontal_spacing=0.002, vertical_spacing=0.006,
    )

    n = 0  # counter — only the first heatmap shows the colour-scale legend
    for i in range(upper_col_idx):
        for j in range(outter_row_idx):
            # Extract the subset for this (upper_col, outter_row) combination
            mask = np.where(
                (data.iloc[:, 0] == data.iloc[:, 0].unique()[i]) &
                (data.iloc[:, 2] == data.iloc[:, 2].unique()[j])
            )
            # Reshape count values into a 2-D matrix (inner_row × lower_col)
            z = data.iloc[np.r_[mask], 4].to_numpy().reshape((inner_row_idx, lower_col_idx))

            # Apply normalization if requested (safe against division by zero)
            if normalization_decision == 'X-axis':
                row_sums = z.sum(axis=1, keepdims=True)
                z = z / np.where(row_sums == 0, 1, row_sums)
            elif normalization_decision == 'Y-axis':
                col_sums = z.sum(axis=0, keepdims=True)
                z = z / np.where(col_sums == 0, 1, col_sums)

            # Extract axis labels from the reshaped data
            x = data.iloc[np.r_[mask], 3].to_numpy().reshape((inner_row_idx, lower_col_idx))[0]
            y = data.iloc[np.r_[mask], 1].to_numpy().reshape((inner_row_idx, lower_col_idx)).T[0]

            # Add heatmap trace (only the first subplot shows the colour scale)
            heatmap = go.Heatmap(z=z, x=x, y=y, colorscale='Viridis', showscale=(n == 0))
            fig.add_trace(heatmap, row=j + 1, col=i + 1)
            n += 1

    # ── Label the outer grid axes ──
    for i in range(upper_col_idx):
        fig.update_xaxes(title_text=data.iloc[:, 0].unique()[i], row=1, col=i + 1, side='top')
    for j in range(outter_row_idx):
        fig.update_yaxes(title_text=data.iloc[:, 2].unique()[j], row=j + 1, col=1)

    fig.update_traces(textfont_size=10)
    fig.update_layout(
        title_text="Exploratory Analysis of Clinical Variables Interactions in H(D)SA Data",
        width=1200, height=600,
    )
    return fig
