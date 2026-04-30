"""
Page 4 — DSA Outcome / Survival Analysis

Transplant outcome analyses regarding DSA status.  Fits a Cox Proportional Hazards
model with user-selected continuous and categorical covariates, then displays:

  - Left-top:    Kaplan-Meier survival curves (DSA=0 vs DSA=1) with 95 % CI bands
  - Left-bottom: Model summary statistics rendered as an image
  - Right:       Interactive CoxPH coefficient table (sortable, filterable, CSV export)
"""

import numpy as np
import pandas as pd
import os
import re
import json

import plotly.graph_objects as go
import plotly.express as px

import dash
from dash import html, dcc, callback, Output, Input
import dash_bootstrap_components as dbc
import dash_table as dt

from lifelines import KaplanMeierFitter, CoxPHFitter
from PIL import Image, ImageDraw, ImageFont

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.copy_on_write = True

dash.register_page(__name__)

# ── Data Loading ──────────────────────────────────────────────────────────────
default_path = "C:/Users/tzhang/DSA"
os.chdir(default_path)

# HCT Essentials: transplant-level clinical + outcome data
hct_df = pd.read_csv('dsa_hct_data.csv', sep=',')
death_dict = pd.read_csv('dsa_death_dict.csv', sep=',')
disease_dict = pd.read_csv('dsa_disease_dict.csv', sep=',')

# DSA dictionary: {RECIPIENT_GUID -> "allele1&allele2&..."}
with open('hla_dsa_dict.json', "r") as jsonfile:
    dsa_dict = json.load(jsonfile)
HLA_save_dsa_RID = list(dsa_dict.keys())

# Recipients with DSA targeting the 8-locus match panel (A, B, C, DRB1)
HLA_MATCH_8_LOCI = ('A', 'B', 'C', 'DRB1')
HLA_save_dsa_hm8_RID = [
    key for key, value in dsa_dict.items()
    if any(allele.startswith(HLA_MATCH_8_LOCI) for allele in value.split('&'))
]

# ── Data Cleaning ────────────────────────────────────────────────────────────
# Identify all death-related columns
death_idx = [col for col in hct_df.columns if 'DEATH' in col]

# Deduplicate to one row per (CRID, TXNUM) and exclude fully matched 8/8
hct_df = hct_df.drop_duplicates(['CRID', 'TXNUM'])
hct_df = hct_df[hct_df['LOW_DONOR_MATCH_CATEGORY'] != '8_of_8']

# Fill missing primary death cause with 'Over_all_survival' (i.e. alive / censored)
hct_df['PRI_DEATH'] = hct_df['PRI_DEATH'].fillna('Over_all_survival')

# B-leader match status: fill NaN as False, cast to int (0/1)
hct_df['B_LEADER_MATCHED_STATUS'] = hct_df['B_LEADER_MATCHED_STATUS'].fillna(False).astype(int)

# Derive HLA match level from the 4 core loci
hct_df['HLA_MATCH_LEVEL_8'] = hct_df[
    ['A_RES_MG', 'B_RES_MG', 'C_RES_MG', 'DRB1_RES_MG']
].apply(
    lambda x: [m for allele in x for m in str(allele).split(':') if pd.notna(allele)],
    axis=1,
)
hct_df['HLA_MATCH_LEVEL'] = hct_df['HLA_MATCH_LEVEL_8'].apply(
    lambda x: x.count('Hmm') + x.count('Imm') + 2 if x != [] else 'NAN'
)

# Convert survival time from "X days" string to months; cap missing at 120 months
hct_df['SURV_TIME'] = hct_df['SURV_TIME'].apply(
    lambda x: 120 if x != x else int(int(re.sub(' days', '', x)) / 30.44)
)

# ── Outcome dropdown options ────────────────────────────────────────────────
filter_category = {}
filter_category['PRI_DEATH'] = hct_df['PRI_DEATH'].unique().tolist()
# Add composite death categories
filter_category['PRI_DEATH'].extend(['DEATH_INFECT', 'DEATH_ORGAN_FAILURE', 'DEATH_MALIGNANCY'])

# ── Page Layout ──────────────────────────────────────────────────────────────
WIDGET_STYLE = {
    'width': '25%', 'display': 'inline-block',
    'font-size': '12px', 'height': '10vh', 'padding': '0 5',
}

layout = dbc.Container([
    html.Hr(),

    # ── Row 1: Covariate & outcome selectors ──
    dbc.Row([
        # Continuous covariates multi-select
        html.Div([
            dbc.Label("Clinical covariates - continuous"),
            dcc.Dropdown(
                id="clinical_var_con",
                options=[{'label': c, 'value': c} for c in hct_df.columns],
                value=['AGE', 'D_AGE', 'SEX', 'YEARTX', 'NUMNMDP8'],
                clearable=False, multi=True, style={"width": "100%"},
            ),
        ], style=WIDGET_STYLE),

        # Categorical covariates multi-select
        html.Div([
            dbc.Label("Clinical covariates - categorical"),
            dcc.Dropdown(
                id="clinical_var_cat",
                options=[{'label': c, 'value': c} for c in hct_df.columns],
                value=['DISEASE', 'RACE', 'GRAFTYPE', 'CONDCLAS_ALL', 'LOW_DONOR_MATCH_CATEGORY', 'SOURCE'],
                clearable=False, multi=True, style={"width": "100%"},
            ),
        ], style=WIDGET_STYLE),

        # Outcome variable selector
        html.Div([
            dbc.Label("Outcome variable"),
            dcc.Dropdown(
                id="pri_outcome",
                options=[{"label": x, "value": x} for x in filter_category['PRI_DEATH']],
                value="Graft_rejection_or_failure",
                clearable=False, multi=False, style={"width": "100%"},
            ),
        ], style=WIDGET_STYLE),

        # HLA loci scope toggle
        html.Div([
            dbc.Label("HLA loci setting"),
            dbc.RadioItems(options=['All', 'A-B-C-DRB1'], value='All', id='hla_loci_set'),
        ], style=WIDGET_STYLE),
    ]),

    # ── Row 2: Survival time cutoff slider ──
    dbc.Row([
        html.Div([
            dbc.Label("Survival time cutoff (months)"),
            dcc.RangeSlider(
                0, 60, value=[60],
                marks={v: str(v) for v in [0, 1, 2, 3, 6, 12, 24, 36, 48, 60]},
                id='survival_time_cutoff',
            ),
        ], style={'width': '75%', 'display': 'inline-block', 'font-size': '12px', 'height': '10vh', 'padding': '0 5'}),
    ]),

    # ── Row 3: Output plots & table ──
    dbc.Row([
        # Left column: KM curve + model summary image
        dbc.Col([
            html.Div([
                dcc.Graph(id='coxph_survival_cruve_plot'),
            ], style={'width': '100%', 'display': 'inline-block', 'font-size': '12px', 'height': '40vh', 'padding': '0 0'}),
            html.Div([
                dcc.Graph(id='coxph_survival_summary', style={'height': '20vh'}),
            ], style={'width': '100%', 'display': 'inline-block', 'font-size': '12px', 'height': '20vh', 'padding': '0 0'}),
        ]),
        # Right column: interactive CoxPH coefficient table
        dbc.Col([
            html.Div(
                id='coxph_survival_summary_table',
                style={'width': '100%', 'display': 'inline-block', 'font-size': '15px', 'height': '100vh', 'padding': '0 0'},
            ),
        ]),
    ]),
])


# ── Helper: format CoxPH summary for display ────────────────────────────────
def _to_sci(x, prec=2):
    """Convert a value to scientific notation string; pass through NaN / empty."""
    if pd.isnull(x) or x == '':
        return x
    try:
        return f"{float(x):.{prec}e}"
    except Exception:
        return x


# ── Main Callback ────────────────────────────────────────────────────────────
@callback(
    Output('coxph_survival_cruve_plot', 'figure'),
    Output('coxph_survival_summary', 'figure'),
    Output('coxph_survival_summary_table', 'children'),
    Input('clinical_var_con', 'value'),
    Input('clinical_var_cat', 'value'),
    Input('pri_outcome', 'value'),
    Input('hla_loci_set', 'value'),
    Input('survival_time_cutoff', 'value'),
)
def coxph_survival_cruve_plot(clinical_var_con, clinical_var_cat,
                              pri_outcome, hla_loci_set, survival_time_cutoff):
    """Fit a CoxPH model and produce three outputs:

    1. Kaplan-Meier survival curves (DSA=0 green vs DSA=1 red) with 95 % CI
    2. Model summary statistics rendered as an image (concordance, AIC, LR test)
    3. Interactive DataTable of CoxPH coefficients
    """
    tmp = hct_df.copy()
    cutoff = int(survival_time_cutoff[0])

    # ── Optionally restrict DSA to 8-locus panel ──
    if hla_loci_set == 'A-B-C-DRB1':
        tmp['DSA_8'] = 0
        tmp.loc[tmp['RECIPIENT_GUID'].isin(HLA_save_dsa_hm8_RID), 'DSA_8'] = 1
        tmp['DSA'] = tmp['DSA_8']

    # ── Build covariate lists ──
    coxph_vars = []
    cat_drop = []

    # Continuous covariates
    con_list = (
        [str(x) for x in clinical_var_con] if isinstance(clinical_var_con, list)
        else ([str(clinical_var_con)] if clinical_var_con else [])
    )

    # Categorical covariates → one-hot encode, drop sparse / reference levels
    if clinical_var_cat is not None:
        cat_input = [str(x) for x in clinical_var_cat] if isinstance(clinical_var_cat, list) else [str(clinical_var_cat)]
        for cat_var in cat_input:
            tmp = pd.get_dummies(tmp, columns=[cat_var], prefix=f'Cat_{cat_var}')
            cat_cols = [c for c in tmp.columns if f'Cat_{cat_var}' in c]
            tmp[cat_cols] = tmp[cat_cols].astype(int)

            # Drop levels with < 10 DSA-positive observations (or keep only the largest)
            dsa_sums = tmp.loc[tmp['DSA'] == 1, cat_cols].sum()
            if len(cat_cols) > 2:
                keep = dsa_sums[dsa_sums >= 10].index
                if len(keep) > 0:
                    cat_drop.extend(dsa_sums[dsa_sums < 10].index)
                else:
                    cat_drop.extend(dsa_sums.sort_values().index[1:])
            else:
                # Binary: drop the less frequent level as reference
                cat_drop.append(dsa_sums.sort_values().index[0])

    cat_list = [c for c in tmp.columns if 'Cat_' in c and c not in cat_drop]

    # Assemble final covariate list: DSA + outcome + time + continuous + categorical
    coxph_vars = ['DSA', pri_outcome, 'SURV_TIME'] + con_list + cat_list

    # ── Define the event indicator ──
    if pri_outcome == 'Over_all_survival':
        # Any death cause → event
        tmp[pri_outcome] = (tmp[death_idx[:7]].fillna(0) != 0).any(axis=1).astype(bool)
    elif pri_outcome in ('DEATH_INFECT', 'DEATH_ORGAN_FAILURE', 'DEATH_MALIGNANCY'):
        tmp[pri_outcome] = tmp[pri_outcome].astype(int).astype(bool)
    else:
        # Specific cause of death
        tmp[pri_outcome] = (tmp[death_idx[:7]].fillna(0) == pri_outcome).any(axis=1).astype(bool)

    # Subset to model columns and apply time cutoff
    tmp = tmp[coxph_vars].copy()
    tmp.loc[tmp['SURV_TIME'] > cutoff, pri_outcome] = False          # censor beyond cutoff
    tmp['SURV_TIME'] = tmp['SURV_TIME'].clip(upper=cutoff)
    tmp = tmp.dropna(subset=coxph_vars).reset_index(drop=True)

    # ── Fit CoxPH model ──
    cph = CoxPHFitter(penalizer=0.01)
    cph.fit(tmp, "SURV_TIME", pri_outcome)

    # ── Fit Kaplan-Meier curves for DSA=0 and DSA=1 ──
    mask_no_dsa = (tmp["DSA"] == 0)
    T, E = tmp["SURV_TIME"], tmp[pri_outcome]

    kmf0 = KaplanMeierFitter()
    kmf0.fit(T[mask_no_dsa], event_observed=E[mask_no_dsa], label="0")
    kmf1 = KaplanMeierFitter()
    kmf1.fit(T[~mask_no_dsa], event_observed=E[~mask_no_dsa], label="1")

    # ── Build KM survival curve figure ──
    fig = go.Figure()

    # DSA=0 (green) — survival line + 95 % CI band
    fig.add_trace(go.Scatter(
        x=kmf0.survival_function_.index, y=kmf0.survival_function_['0'],
        line=dict(shape='hv', width=3, color='rgb(36, 180, 36)'),
        legendgroup='DSA=0', name='DSA=0', showlegend=True,
    ))
    fig.add_trace(go.Scatter(
        x=kmf0.confidence_interval_.index,
        y=kmf0.confidence_interval_['0_upper_0.95'],
        line=dict(shape='hv', width=0), showlegend=False,
    ))
    fig.add_trace(go.Scatter(
        x=kmf0.confidence_interval_.index,
        y=kmf0.confidence_interval_['0_lower_0.95'],
        line=dict(shape='hv', width=0),
        fill='tonexty', fillcolor='rgba(36, 180, 36, 0.4)', showlegend=False,
    ))

    # DSA=1 (red) — survival line + 95 % CI band
    fig.add_trace(go.Scatter(
        x=kmf1.survival_function_.index, y=kmf1.survival_function_['1'],
        line=dict(shape='hv', width=3, color='rgb(180, 36, 36)'),
        legendgroup='DSA=1', name='DSA=1', showlegend=True,
    ))
    fig.add_trace(go.Scatter(
        x=kmf1.confidence_interval_.index,
        y=kmf1.confidence_interval_['1_upper_0.95'],
        line=dict(shape='hv', width=0), showlegend=False,
    ))
    fig.add_trace(go.Scatter(
        x=kmf1.confidence_interval_.index,
        y=kmf1.confidence_interval_['1_lower_0.95'],
        line=dict(shape='hv', width=0),
        fill='tonexty', fillcolor='rgba(180, 36, 36, 0.4)', showlegend=False,
    ))

    fig.update_xaxes(showline=True, linewidth=2, linecolor='black')
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black')
    fig.update_layout(
        xaxis_title="Months since Transplantation",
        yaxis_title="Survival probability",
        margin=dict(r=0, t=10, l=0),
        font_size=14, xaxis_title_font_size=18, yaxis_title_font_size=18,
        plot_bgcolor='white',
        legend=dict(x=1.01, y=0.9, title="DSA Group", font=dict(size=12)),
    )

    # ── Build model summary statistics image ──
    sr = cph.log_likelihood_ratio_test()
    headers = {
        "model": "'lifelines.CoxPHFitter'",
        "duration col": "'SURV_TIME'",
        "event col": f"'{pri_outcome}'",
        "baseline estimation": "breslow",
        "number of observations": f"{tmp.shape[0]:g}",
        "number of events observed": f"{tmp[tmp[pri_outcome] != 0].shape[0]:g}",
        "Concordance": f"{cph.concordance_index_:.2f}",
        "Partial AIC": f"{cph.AIC_partial_:.2f}",
        "log-likelihood ratio test": f"{sr.test_statistic:.2f} on {sr.degrees_freedom} df",
    }
    headers_df = pd.DataFrame.from_dict(headers, orient='index')
    headers_df.columns = [''] * len(headers_df.columns)

    # Render summary text as a PIL image, then display via px.imshow
    summary_text = headers_df.to_string()
    font = ImageFont.truetype("arial.ttf", 15)

    # Measure text size, create appropriately sized image
    tmp_img = Image.new("RGB", (1, 1), "white")
    tmp_draw = ImageDraw.Draw(tmp_img)
    bbox = tmp_draw.textbbox((10, 10), summary_text, font=font)
    img_w, img_h = bbox[2] - bbox[0] + 20, bbox[3] - bbox[1] + 40

    image = Image.new("RGB", (img_w, img_h), "white")
    draw = ImageDraw.Draw(image)
    draw.text((10, 10), summary_text, font=font, fill="black")

    # Crop to content area
    image_cropped = image.crop((10, 20, (img_w + img_w) // 2, (img_h + img_h) // 2))
    summary_fig = px.imshow(image_cropped, width=600, height=300)
    summary_fig.update_xaxes(visible=False, showticklabels=False)
    summary_fig.update_yaxes(visible=False, showticklabels=False)

    # ── Build CoxPH coefficient table ──
    summary = cph.summary.copy()
    # Format coef, exp(coef), se(coef) to 2 decimal places
    for i in np.r_[0:3, 8:9]:
        summary.iloc[:, i] = summary.iloc[:, i].apply(lambda x: f"{x:.2f}")
    # Format p-values in scientific notation
    summary.iloc[:, 9] = summary.iloc[:, 9].apply(lambda x: _to_sci(x, prec=2))

    # Select key columns: coef, exp(coef), se(coef), z, p
    coef_df = summary.iloc[:, np.r_[0:3, 8:10]].copy()
    coef_df['variable'] = coef_df.index
    col_order = ['variable'] + coef_df.columns[:-1].tolist()
    coef_df = coef_df.reindex(columns=col_order)

    table = html.Div(dt.DataTable(
        data=coef_df.to_dict('records'),
        columns=[{"name": c, "id": c} for c in coef_df.columns],
        page_size=20,
        filter_action="native",
        sort_action="native",
        export_format="csv",
        export_headers='all',
    ))

    return fig, summary_fig, table
