"""
Page 3 — DSA Match

Exploratory analyses of HLA matching characteristics between transplants
with and without Donor-Specific Antibodies (DSA).

Left plot:  Variable comparison (violin / bar) between DSA-positive and DSA-negative
            groups, with statistical test annotations (t-test or proportions z-test).
Right plot: Per-sample DSA allele scatter plot faceted by HLA locus, with marker size
            proportional to allele frequency and shape indicating MFI data availability.
"""

import scipy.stats as stats
from statsmodels.stats.proportion import proportions_ztest
import numpy as np
import pandas as pd
import plotly.express as px
import os
import math
import json
import dash
from dash import html, dcc, callback, Output, Input
import dash_bootstrap_components as dbc
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.copy_on_write = True

dash.register_page(__name__)

# ── Data Loading ──────────────────────────────────────────────────────────────
default_path = "C:/Users/tzhang/DSA"
os.chdir(default_path)

# Antigen + allele level antibody data with demographics, frequencies, and CIWD status
hsa_antigen_allele_freq_ciwd_df = pd.read_csv(
    'dsa_data3.csv', sep=',', header=0,
)

# HCT Essentials transplant-level clinical data
hla_save_df = pd.read_csv(
    'dsa_hla_data.csv', sep=',',
)
hla_save_df.fillna('nan', inplace=True)
hla_save_df = hla_save_df.drop_duplicates(['CRID', 'TXNUM'])

# Derive HLA match level from the 4 core loci (A, B, C, DRB1)
hla_save_df['HLA_MATCH_LEVEL_8'] = hla_save_df[
    ['A_RES_MG', 'B_RES_MG', 'C_RES_MG', 'DRB1_RES_MG']
].apply(
    lambda x: [m for allele in x for m in str(allele).split(':') if pd.notna(allele)],
    axis=1,
)
hla_save_df['HLA_MATCH_LEVEL'] = hla_save_df['HLA_MATCH_LEVEL_8'].apply(
    lambda x: x.count('Hmm') + x.count('Imm') + 2 if x != [] else 'NAN'
)

# DSA dictionary: {RECIPIENT_GUID -> "allele1&allele2&..."} for patients with confirmed DSA
with open('hla_dsa_dict.json', "r") as jsonfile:
    dsa_dict = json.load(jsonfile)

# Recipients with DSA targeting the 8-locus match panel (A, B, C, DRB1)
HLA_MATCH_8_LOCI = ('A', 'B', 'C', 'DRB1')
HLA_save_dsa_hm8_RID = [
    key for key, value in dsa_dict.items()
    if any(allele.startswith(HLA_MATCH_8_LOCI) for allele in value.split('&'))
]

# All recipients with any DSA
HLA_save_dsa_RID = list(dsa_dict.keys())

# ── Rename columns for consistency ──────────────────────────────────────────
hsa_antigen_allele_df = hsa_antigen_allele_freq_ciwd_df.rename(columns={
    'SEX_CDE': 'SEX', 'broad_race_mapped': 'RACE', 'HLA_LOCUS_CDE': 'HLA',
})
hsa_antigen_allele_df = hsa_antigen_allele_df.sort_values(['SEX', 'RACE', 'HLA'])

# ── Build filter dropdown options ───────────────────────────────────────────
filter_category = (
    hsa_antigen_allele_df['SEX'].apply(lambda x: f'SEX={x}' if x == x else "None").unique().tolist()
    + hsa_antigen_allele_df['RACE'].apply(lambda x: f'RACE={x}' if x == x else "None").unique().tolist()
    + hsa_antigen_allele_df['HLA'].apply(lambda x: f'HLA={x}' if x == x else "None").unique().tolist()
)

# ── Constants ───────────────────────────────────────────────────────────────
HLA_LOCUS_ORDER = ['A', 'B', 'C', 'DRB1', 'DRB3', 'DRB4', 'DRB5', 'DPA1', 'DPB1', 'DQA1', 'DQB1']
MFI_SHAPES = {0: 'circle', 1: 'arrow', 2: 'star'}

# Common plot styling
PLOT_KWARGS = dict(
    plot_bgcolor='white',
    font=dict(family="Arial, monospace", size=10, color="black"),
    margin=dict(l=20, r=20, t=40, b=20),
)
AXIS_LINE = dict(showline=True, linecolor='black')


# ── Shared Helper: apply 'COL=value' filters ────────────────────────────────
def _apply_filters(data, variable_filters):
    """Parse and apply 'COL=value' style filters from the dropdown."""
    if variable_filters is None:
        return data
    items = [str(x) for x in variable_filters] if isinstance(variable_filters, list) else [str(variable_filters)]
    for vf in items:
        col, val = vf.split('=', 1)
        data = data[data[col] == val]
    return data


def _normalizer(x):
    """Normalise a group so values sum to 1 (for percentage bar charts)."""
    s = x.sum()
    return x / s if s != 0 else x


def _proportions_test_annotation(df_tmp_2, Variable_sel):
    """Compute a proportions z-test on the most frequent category and return
    the p-value annotation text."""
    main_cat = df_tmp_2.sort_values('count', ascending=False)[Variable_sel].iloc[0]
    count = [
        int(df_tmp_2.loc[(df_tmp_2[Variable_sel] == main_cat) & (df_tmp_2['DSA'] == 'Y'), 'count'].sum()),
        int(df_tmp_2.loc[(df_tmp_2[Variable_sel] == main_cat) & (df_tmp_2['DSA'] == 'N'), 'count'].sum()),
    ]
    nobs = [
        int(df_tmp_2.loc[df_tmp_2['DSA'] == 'Y', 'count'].sum()),
        int(df_tmp_2.loc[df_tmp_2['DSA'] == 'N', 'count'].sum()),
    ]
    _, p_value = proportions_ztest(count, nobs)
    return f"Proportions z-test p-value: {p_value:.2e}"


def _add_pvalue_annotation(fig, text):
    """Add a p-value annotation at the top-centre of the figure."""
    fig.add_annotation(
        x=0.5, y=1.01, xref="paper", yref="paper",
        text=text, showarrow=False,
        font=dict(family="Arial, monospace", size=12, color="black"),
    )


# ── Page Layout ──────────────────────────────────────────────────────────────
CTRL_STYLE = {
    'width': '25%', 'display': 'inline-block',
    'font-size': '12px', 'height': '10vh', 'padding': '0 10',
}

layout = dbc.Container([
    html.Hr(),

    # ── Control widgets ──
    dbc.Row([
        # Filter dropdown (SEX, RACE, HLA)
        html.Div([
            dbc.Label("Variable for filtering"),
            dcc.Dropdown(
                options=[{'label': c, 'value': c} for c in filter_category if c == c and c != 'None'],
                value=None, multi=True, id='Variable_filters',
            ),
        ], style=CTRL_STYLE),

        # Variable selector for the comparison plot (left)
        html.Div([
            dbc.Label("Variable for comparisons"),
            dcc.Dropdown(
                options=[{'label': c, 'value': c} for c in hla_save_df.columns],
                value='LOW_DONOR_MATCH_CATEGORY', multi=False, id='Variable_selected',
            ),
        ], style=CTRL_STYLE),

        # Variable type toggle
        html.Div([
            dbc.Label("Variable type for comparisons"),
            dcc.Dropdown(
                options=['continuous', 'categorical_base_value', 'categorical_percentage'],
                value='categorical_percentage', multi=False, id='Variable_type',
            ),
        ], style=CTRL_STYLE),

        # DSA grouping mode: colour by DSA or put DSA on X-axis
        html.Div([
            dbc.Label("DSA and variable setting"),
            dbc.RadioItems(options=['DSA_group', 'DSA_x_axis'], value='DSA_group', id='DSA_setting'),
        ], style=CTRL_STYLE),
    ]),

    # MFI cutoff slider
    dbc.Row([
        html.Div([
            dbc.Label("MFI cutoffs"),
            dcc.RangeSlider(
                0, 5000, value=[0],
                marks={v: str(v) for v in [0, 100, 500, 1000, 2000, 5000]},
                id='mfi_cutoff',
            ),
        ], style={'width': '90%', 'display': 'inline-block', 'font-size': '12px', 'height': '10vh', 'padding': '0 5'}),
    ]),

    # ── Plot row ──
    dbc.Row([
        # Left: variable comparison plot
        html.Div([dcc.Graph(id='var_comp_plot')],
                 style={'width': '40%', 'display': 'inline-block', 'font-size': '12px', 'height': '50vh', 'padding': '0 10'}),
        # Right: per-sample DSA scatter plot
        html.Div([dcc.Graph(id='dsa_sample_plot')],
                 style={'width': '60%', 'display': 'inline-block', 'font-size': '12px', 'height': '55vh', 'padding': '0 10'}),
    ]),
])


# ── Callback 1: Variable Comparison Plot (left) ─────────────────────────────
@callback(
    Output('var_comp_plot', 'figure'),
    Input('Variable_selected', 'value'),
    Input('Variable_filters', 'value'),
    Input('Variable_type', 'value'),
    Input('DSA_setting', 'value'),
    Input('mfi_cutoff', 'value'),
)
def var_comp_plot(Variable_sel, Variable_filters, Variable_ty, DSA_set, mfi_ctf):
    """Compare a user-selected clinical variable between DSA and non-DSA groups.

    Supports:
    - continuous  → violin + box with Welch's t-test
    - categorical_base_value → grouped bar with proportions z-test
    - categorical_percentage → stacked percentage bar with proportions z-test

    DSA_set controls whether DSA status is the colour grouping or the X-axis.
    mfi_ctf filters DSA recipients by minimum MFI threshold.
    """
    hsa_tmp = hsa_antigen_allele_df.copy()

    # Determine DSA recipient list based on MFI cutoff
    if int(mfi_ctf[0]) == 0:
        dsa_rid_mfi = HLA_save_dsa_RID
    else:
        # Only count recipients whose DSA alleles exceed the MFI threshold
        dsa_rid_mfi = hsa_antigen_allele_df.loc[
            (hsa_antigen_allele_df['MFI_VALUE_NUM'].notna()) &
            (hsa_antigen_allele_df['MFI_VALUE_NUM'] > int(mfi_ctf[0])) &
            (hsa_antigen_allele_df['DSA_HAR'] == 'Y'),
            'RECIPIENT_GUID',
        ].unique()

    # Apply demographic filters and restrict HLA-SAVE data accordingly
    hsa_tmp = _apply_filters(hsa_tmp, Variable_filters)
    if Variable_filters is not None:
        save_df = hla_save_df[hla_save_df['RECIPIENT_GUID'].isin(hsa_tmp['RECIPIENT_GUID'])].copy()
    else:
        save_df = hla_save_df.copy()

    # Label each recipient as DSA-positive or -negative
    save_df['DSA'] = save_df['RECIPIENT_GUID'].apply(lambda x: 'Y' if x in dsa_rid_mfi else 'N')
    save_df = save_df.drop_duplicates(subset=['RECIPIENT_GUID']).copy()

    # Clean: drop NaN / 'nan' in the comparison variable
    df2 = save_df.dropna(subset=[Variable_sel, 'DSA']).copy()
    df2 = df2[df2[Variable_sel] != 'nan'].copy()

    # ── Continuous variable: violin + Welch's t-test ──
    if Variable_ty == 'continuous':
        df2['DSA'] = pd.Categorical(df2['DSA'], categories=['Y', 'N'], ordered=True)
        dsa_y = df2[df2['DSA'].astype(str) == 'Y'][Variable_sel].astype(float)
        dsa_n = df2[df2['DSA'].astype(str) == 'N'][Variable_sel].astype(float)
        _, p_val = stats.ttest_ind(dsa_y, dsa_n, equal_var=False, nan_policy='omit')

        fig = px.violin(df2, y=Variable_sel, color='DSA', box=True, title=f"{Variable_sel} comparison")
        _add_pvalue_annotation(fig, f"T-test p-value: {p_val:.2e}")
        fig.update_layout(
            legend_title_text='nonDSA vs DSA',
            xaxis_title=Variable_sel, yaxis_title='Count', **PLOT_KWARGS,
        )

    # ── Categorical variable ──
    else:
        df2 = df2.value_counts([Variable_sel, 'DSA']).reset_index()
        df2['DSA'] = pd.Categorical(df2['DSA'], categories=['Y', 'N'], ordered=True)
        df2[Variable_sel] = pd.Categorical(
            df2[Variable_sel].astype(str),
            categories=df2[Variable_sel].sort_values().unique().astype(str),
            ordered=True,
        )
        p_text = _proportions_test_annotation(df2, Variable_sel)

        # Determine X-axis and colour based on DSA_set
        if DSA_set == 'DSA_group':
            x_col, color_col, legend = Variable_sel, 'DSA', 'nonDSA vs DSA'
        else:
            x_col, color_col, legend = 'DSA', Variable_sel, Variable_sel

        # Apply percentage normalization if requested
        if 'percentage' in Variable_ty:
            group_key = Variable_sel if DSA_set == 'DSA_group' else 'DSA'
            df2['count'] = df2.groupby(group_key)['count'].transform(_normalizer)
            y_title = 'Percentage'
        else:
            y_title = 'Count'

        fig = px.bar(df2, x=x_col, y='count', color=color_col, title=f"{Variable_sel} comparison")
        _add_pvalue_annotation(fig, p_text)
        fig.update_layout(legend_title_text=legend, yaxis_title=y_title, **PLOT_KWARGS)

    fig.update_xaxes(**AXIS_LINE)
    fig.update_yaxes(**AXIS_LINE)
    return fig


# ── Callback 2: Per-Sample DSA Allele Scatter (right) ───────────────────────
@callback(
    Output('dsa_sample_plot', 'figure'),
    Input('Variable_filters', 'value'),
    Input('mfi_cutoff', 'value'),
)
def dsa_sample_plot(Variable_filters, mfi_ctf):
    """Scatter plot showing individual DSA alleles for all DSA-positive recipients.

    - Faceted by HLA locus
    - Colour: DSA_HAR (Y = confirmed DSA allele, N = non-DSA allele)
    - Size: proportional to HLA allele frequency (log-scaled)
    - Shape: arrow if MFI data exists, circle otherwise
    """
    hsa_tmp = _apply_filters(hsa_antigen_allele_df.copy(), Variable_filters)

    # Select DSA-positive recipients with valid antigen annotations
    cols = ['RECIPIENT_GUID', 'SEX', 'RACE', 'HLA', 'ANTIGEN_HAR', 'ALLELE_HAR',
            'MFI_VALUE_NUM', 'HLA_FREQ', 'CIWD3', 'CIWD2', 'DSA', 'DSA_HAR']
    dsa_df = hsa_tmp.loc[
        (hsa_tmp['RECIPIENT_GUID'].isin(HLA_save_dsa_RID)) &
        (hsa_tmp['ANTIGEN_HAR'].notna()) &
        (hsa_tmp['ANTIGEN_HAR'] != 'nan'),
        cols,
    ]

    # Apply MFI threshold filter
    if int(mfi_ctf[0]) != 0:
        dsa_df = dsa_df.loc[
            (dsa_df['MFI_VALUE_NUM'].notna()) &
            (dsa_df['MFI_VALUE_NUM'] > int(mfi_ctf[0]))
        ]

    dsa_df = dsa_df.reset_index(drop=True)

    # Map CIWD frequency categories to short labels (C / I / WD / nan)
    # NOTE: the original code had `.value_counts()` appended to the CIWD2 mapping
    # which replaced the column with a frequency Series — this is now fixed.
    dsa_df['CIWD2'] = dsa_df['CIWD2'].astype(str).apply(
        lambda x: 'C' if x == "['Common']" else ('WD' if x == "['Well-Documented']" else 'nan')
    )
    dsa_df['CIWD3'] = dsa_df['CIWD3'].astype(str).apply(
        lambda x: 'C' if x == "['C']" else ('I' if x == "['I']" else ('WD' if x == "['WD']" else 'nan'))
    )
    # Prefer CIWD3 over CIWD2; fall back to 'nan'
    dsa_df['CIWD_FREQ'] = dsa_df[['CIWD3', 'CIWD2']].apply(
        lambda x: x.iloc[0] if x.iloc[0] != 'nan' else (x.iloc[1] if x.iloc[1] != 'nan' else 'nan'),
        axis=1,
    )

    # Flag whether MFI data exists (for marker shape encoding)
    dsa_df['MFI_check'] = dsa_df['MFI_VALUE_NUM'].apply(lambda x: 1 if x == x and x != 'nan' else 0)

    # Prepare plot data
    plot_df = dsa_df.fillna('nan').copy()
    plot_df['RECIPIENT_GUID'] = plot_df['RECIPIENT_GUID'].astype(str)
    plot_df['HLA'] = pd.Categorical(plot_df['HLA'], categories=HLA_LOCUS_ORDER, ordered=True)
    plot_df = plot_df.sort_values(['HLA', 'DSA_HAR'], ascending=True).reset_index(drop=True)

    # Size: log-scaled allele frequency (larger = more common allele)
    hover_cols = ['SEX', 'RACE', 'ANTIGEN_HAR', 'HLA_FREQ', 'CIWD_FREQ', 'MFI_VALUE_NUM']
    fig = px.scatter(
        plot_df, x="RECIPIENT_GUID", y="ANTIGEN_HAR", facet_col="HLA",
        color=plot_df['DSA_HAR'],
        size=plot_df['HLA_FREQ'].astype(float).apply(
            lambda x: 5 * (7 + round(math.log10(x))) if x != 0 else 1
        ),
        hover_data=hover_cols,
    )

    fig.update_layout(
        xaxis_type='category', yaxis_type='category',
        legend_title_text='DSA_loci', title='DSA distribution',
        xaxis_title='index', yaxis_title='DSA allele',
        legend=dict(itemsizing='trace'),
        **PLOT_KWARGS,
    )
    fig.update_xaxes(showticklabels=False, title=None, **AXIS_LINE)
    fig.update_yaxes(**AXIS_LINE, gridcolor='lightgrey')

    # Use arrow markers for alleles with MFI data, circles otherwise
    fig.update_traces(
        marker=dict(symbol=[MFI_SHAPES[x] for x in plot_df['MFI_check']]),
        selector=dict(mode="markers"),
    )
    # Rotate facet labels for readability
    for annotation in fig['layout']['annotations']:
        annotation['textangle'] = 315

    return fig
