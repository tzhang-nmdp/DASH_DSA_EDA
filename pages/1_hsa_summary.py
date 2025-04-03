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
pd.options.mode.copy_on_write = True
import warnings
# Suppress all FutureWarnings
warnings.simplefilter(action='ignore', category=FutureWarning)
dash.register_page(__name__)

# change the working directory
default_path="c:/Users/tzhang/Desktop/Project/Registry model/dsa"
os.chdir(default_path)
cwd = os.getcwd()
print(cwd)

# reading the dsa data file
df = pd.read_csv('Patient_Antibody_Report_original_data_allele.sim.csv',sep=',')

# Set the subgrouping / filtering variable list 
group_category = [ 'SEX', 'RACE', 'HLA_LOCUS','HLA_MATCH10','HLA_MATCH8', 'HSA_CNT', 'MFI_GRP']
group_category

filter_category = df['SEX'].apply(lambda x: 'SEX' + "=" + str(x) if x == x else "None").unique().tolist() + df['RACE'].apply(lambda x: 'RACE' + "=" + str(x) if x == x else "None").unique().tolist() + df['HLA_LOCUS'].apply(lambda x: 'HLA_LOCUS' + "=" + str(x) if x == x else "None").unique().tolist() + df['HLA_MATCH10'].apply(lambda x: 'HLA_MATCH10' + "=" + str(x) if x == x else "None").unique().tolist() + df['HLA_MATCH8'].apply(lambda x: 'HLA_MATCH8' + "=" + str(x) if x == x else "None").unique().tolist() + df['MFI_GRP'].apply(lambda x: 'MFI_GRP' + "=" + str(x) if x == x else "None").unique().tolist()  + df['HSA_CNT'].apply(lambda x: 'HSA_CNT' + "=" + str(x) if x == x and x < 10 else ( 'HSA_CNT' + '=10' if x == x and x >=10 else"None")).unique().tolist()
# filter_category = [x for x in filter_category if x == x]
# filter_category = filter_category + ['None']
filter_category

layout=dbc.Container([
# title    
            # html.Div(style={'backgroundColor': colors['background'],'color':'blue', 'margin':'0px', 'padding':'20px'},children="DSA exploratory analyses"),
            html.Hr(),

# x_axis widget (plot 1 only)      
            dbc.Row(
                [html.Div(children=[ 
            dbc.Label("Variable selection for X axis of the plot1"),                                        
            dbc.RadioItems(options=group_category,value=group_category[2],id='control_x_plot1'),   
    ], style={'width': '20%', 'display': 'inline-block',  'font-size': '12px',  'height':'20vh', 'padding': '0 10'}),  
                        
# subgroup widget (both plots)               
            html.Div(children=[
            dbc.Label("Variable selection for subgroup of both plots"),                  
                dbc.RadioItems(options=group_category,value=group_category[0],id='control_group'), 
    ], style={'width': '20%', 'display': 'inline-block', 'font-size': '12px', 'height':'20vh', 'padding': '0 10'}),  
                           
# analytic set widget (both plots)               
            html.Div(children=[
            dbc.Label("DSA count analytic set for both plots"),                  
                dbc.RadioItems(options=['Total_HSA_count', 'Total_reciepient_count'],value='Total_reciepient_count',id='control_count'), 
    ], style={'width': '20%', 'display': 'inline-block', 'font-size': '12px', 'height':'20vh', 'padding': '0 5'}),
 
                                    
# filter widget                                  
            html.Div(children=[dcc.Dropdown(
            id="control_filter",
            #value=filter_category[0],
            options=filter_category,
            multi=True,
        )],style={'width': '20%', 'color':'blue', 'fontColor': 'blue', 'font-size': '12px', 'height':'20vh',
          'margin':'20px', 'padding':'0 10px'})  
            ]), 
                        
# plots with DSA and MFI distributions            
            html.Div(children=[
        dcc.Graph(
            id='HSA_count'
        ),  
    ], style={'width': '49%', 'display': 'inline-block', 'font-size': '12px', 'height':'10vh', 'padding': '0 10'}),                      
            html.Div(children=[        
        dcc.Graph(
            id='mfi_count'
        ),        
    ], style={'width': '49%', 'display': 'inline-block', 'font-size': '12px', 'height':'10vh','padding': '0 10'}),
            ])

# call the input and output for DSA distributions 
@callback(
    Output(component_id='HSA_count', component_property='figure'),
    Input(component_id='control_x_plot1', component_property='value'),
    Input(component_id='control_group', component_property='value'),   
    Input(component_id='control_count', component_property='value'),     
    Input(component_id='control_filter', component_property='value')    
)
# define plot for DSA distributions 
def update_HSA_graph( x_plot1_chosen, subgroup_chosen,count_chosen,  filter_chosen):   
    temp_df=df.copy() 
    temp_df.loc[temp_df['HSA_CNT'] > 10,'HSA_CNT'] = 10       
    temp_df.loc[temp_df['MFI_VALUE_NUM'] > 20000,'MFI_VALUE_NUM'] = 20000      
    temp_df['subgroup'] = temp_df[subgroup_chosen]    
     
    if isinstance(filter_chosen, list): 
        for filter_itm in filter_chosen: 
            if filter_itm not in [None, "None", "NA", "na","nan", ""] and is_numeric_dtype(temp_df[filter_itm.split('=')[0]]):
                temp_df = temp_df.iloc[np.r_[temp_df[filter_itm.split('=')[0]] ==  int(filter_itm.split('=')[1])],:]   
            elif filter_itm not in [None, "None", "NA", "na","nan", ""] and is_string_dtype(temp_df[filter_itm.split('=')[0]]):
                temp_df = temp_df.iloc[np.r_[temp_df[filter_itm.split('=')[0]] ==  filter_itm.split('=')[1]],:]                    
    else:
        if filter_chosen not in [None, "None", "NA", "na","nan", ""] and is_numeric_dtype(temp_df[filter_chosen.split('=')[0]]):
            temp_df = temp_df.iloc[np.r_[temp_df[filter_chosen.split('=')[0]] ==  int(filter_chosen.split('=')[1])],:]    
        elif filter_chosen not in [None, "None", "NA", "na","nan", ""] and is_string_dtype(temp_df[filter_chosen.split('=')[0]]):
            temp_df = temp_df.iloc[np.r_[temp_df[filter_chosen.split('=')[0]] ==  filter_chosen.split('=')[1]],:]                       
      
    if count_chosen == 'Total_reciepient_count':
        temp_df_unique = temp_df.sort_values('MFI_GRP').drop_duplicates(subset=['RID'], keep='last')  
    elif x_plot1_chosen=='HLA_LOCUS' or subgroup_chosen=='HLA_LOCUS':
        temp_df_unique = temp_df.drop_duplicates(subset=['RID', 'ALLELE_HAR'])  
    else:
        temp_df_unique = temp_df.drop_duplicates(subset=['RID', subgroup_chosen])      
    temp_df_unique=temp_df_unique.sort_values([x_plot1_chosen,'subgroup'])   
    fig=px.histogram(temp_df_unique,x = x_plot1_chosen, color = 'subgroup', histfunc='sum', opacity = 0.6 ).update_layout(font=dict(size=12), yaxis_title = 'DSA counts')
    return fig

# call the input and output for MFI distributions 
@callback(
    Output(component_id = 'mfi_count', component_property = 'figure'),
    Input(component_id = 'control_group', component_property = 'value'),
    Input(component_id='control_count', component_property='value'),          
    Input(component_id = 'control_filter', component_property = 'value')    
)
# define plot for MFI distributions 
def update_mfi_graph(subgroup_chosen, count_chosen, filter_chosen): 
    temp_df=df.copy()  
    temp_df.loc[temp_df['HSA_CNT'] > 10,'HSA_CNT'] = 10       
    temp_df.loc[temp_df['MFI_VALUE_NUM'] > 20000,'MFI_VALUE_NUM'] = 20000      
    temp_df['subgroup'] = temp_df[subgroup_chosen]    
        
    print(subgroup_chosen)    
    print(filter_chosen)
    if isinstance(filter_chosen, list): 
        for filter_itm in filter_chosen: 
            if filter_itm not in [None, "None", "NA", "na","nan", ""] and is_numeric_dtype(temp_df[filter_itm.split('=')[0]]):
                temp_df = temp_df.iloc[np.r_[temp_df[filter_itm.split('=')[0]] ==  int(filter_itm.split('=')[1])],:]   
            elif filter_itm not in [None, "None", "NA", "na","nan", ""] and is_string_dtype(temp_df[filter_itm.split('=')[0]]):
                temp_df = temp_df.iloc[np.r_[temp_df[filter_itm.split('=')[0]] ==  filter_itm.split('=')[1]],:]                    
    else:
        if filter_chosen not in [None, "None", "NA", "na","nan", ""] and is_numeric_dtype(temp_df[filter_chosen.split('=')[0]]):
            temp_df = temp_df.iloc[np.r_[temp_df[filter_chosen.split('=')[0]] ==  int(filter_chosen.split('=')[1])],:]    
        elif filter_chosen not in [None, "None", "NA", "na","nan", ""] and is_string_dtype(temp_df[filter_chosen.split('=')[0]]):
            temp_df = temp_df.iloc[np.r_[temp_df[filter_chosen.split('=')[0]] ==  filter_chosen.split('=')[1]],:]            
         
    if count_chosen == 'Total_reciepient_count':
        temp_df_unique = temp_df.sort_values('MFI_GRP').drop_duplicates(subset=['RID'], keep='last')
    else:
        temp_df_unique = temp_df.drop_duplicates(subset=['RID', 'ALLELE_HAR'])    
    temp_df_unique=temp_df_unique.sort_values('subgroup')      
    fig=px.histogram(temp_df_unique,x = 'MFI_VALUE_NUM', color = 'subgroup',histfunc = 'avg' , opacity = 0.6, nbins=100 ).update_layout(font=dict(size=12), yaxis_title = 'MFI counts')
    fig.update_xaxes(range=[0,20000])
    return fig
