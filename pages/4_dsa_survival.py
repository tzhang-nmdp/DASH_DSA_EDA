import datetime
import numpy as np
import pandas as pd
from pandas.api.types import is_string_dtype
from pandas.api.types import is_numeric_dtype
import plotly
import plotly.graph_objects as go
import plotly.express as px
from ipywidgets import widgets
import nbformat
import os, sys, math
import dash
from dash import Dash, html, dash_table, dcc, callback, Output, Input
import plotly.express as px
import dash_mantine_components as dmc
from plotly.subplots import make_subplots
import dash_bootstrap_components as dbc
import dash_html_components as html
import warnings
# Code from: https://github.com/plotly/dash-labs/tree/main/docs/demos/multi_page_example1
dash.register_page(__name__)

# reading the dsa data file
df = pd.read_csv('c:/Users/tzhang/Desktop/Project/Registry model/dsa/temp_df_dsa_cnt_table_hla_loci.csv',sep='\t')
filter_category = df['SEX'].apply(lambda x: 'SEX' + "=" + str(x) if x == x else "None").unique().tolist() + df['RACE'].apply(lambda x: 'RACE' + "=" + str(x) if x == x else "None").unique().tolist() + df['HLA_LOCUS'].apply(lambda x: 'HLA_LOCUS' + "=" + str(x) if x == x else "None").unique().tolist() + df['HLA_MATCH10'].apply(lambda x: 'HLA_MATCH10' + "=" + str(x) if x == x else "None").unique().tolist()

layout =dbc.Container([ 
        dbc.Row([               
        html.Div(children=[
        dbc.Label("Interation variables at upper columns"),
        dcc.Dropdown(
            id="test",
            options=[{"label": x, "value": x} for x in df.columns],
            value="HLA_MATCH10",
            clearable=False,
            style={"width": "100%"}
        )], style={'width': '25%', 'display': 'inline-block', 'font-size': '12px', 'height':'10vh', 'padding': '0 5'})
    ])
    ])        
        