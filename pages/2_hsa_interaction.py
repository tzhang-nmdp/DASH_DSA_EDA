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
# Code from: https://github.com/plotly/dash-labs/tree/main/docs/demos/multi_page_example1
dash.register_page(__name__)

# reading the dsa data file
df = pd.read_csv('c:/Users/tzhang/Desktop/Project/Registry model/dsa/temp_df_dsa_cnt_table_hla_loci_combined.csv',sep=',')
print(df.columns)

layout =dbc.Container([ 
        dbc.Row([               
        html.Div(children=[
        dbc.Label("Interation variables at upper columns"),
        dcc.Dropdown(
            id="upper_col",
            options=[{"label": x, "value": x} for x in df.columns[:-2]],
            value="HLA_MATCH10",
            clearable=False,
            style={"width": "100%"}
        )], style={'width': '25%', 'display': 'inline-block', 'font-size': '12px', 'height':'10vh', 'padding': '0 5'}), 
        
        html.Div(children=[
        dbc.Label("Interation variables at outter rows"),         
        dcc.Dropdown(
            id="outter_row",
            options=[{"label": x, "value": x} for x in df.columns[:-2]],
            value="SEX",
            clearable=False,
            style={"width": "100%"}
        )], style={'width': '25%', 'display': 'inline-block', 'font-size': '12px', 'height':'10vh', 'padding': '0 5'}),  
        
        html.Div(children=[
        dbc.Label("Interation variables at inner rows"),        
        dcc.Dropdown(
            id="inner_row",
            options=[{"label": x, "value": x} for x in df.columns[:-2]],
            value="RACE",
            clearable=False,
            style={"width": "100%"}
        )], style={'width': '25%', 'display': 'inline-block', 'font-size': '12px', 'height':'10vh', 'padding': '0 5'}),  
        
        html.Div(children=[
        dbc.Label("Interation variables at lower columns"),        
        dcc.Dropdown(
            id="lower_col",
            options=[{"label": x, "value": x} for x in df.columns[:-2]],
            value="HLA_LOCUS",
            clearable=False,
            style={"width": "100%"}
        )], style={'width': '25%', 'display': 'inline-block', 'font-size': '12px', 'height':'10vh', 'padding': '0 5'}), 
        
        # html.Div(children=[
        # dbc.Label("Filtering variables(subgroups)"),    
        # dcc.Dropdown(
        #     id="filter_factors",
        #     value=filter_category[0],
        #     options=filter_category,
        #     multi=True,
        # )], style={'width': '40%', 'display': 'inline-block', 'font-size': '12px', 'height':'10vh', 'padding': '0 5'}), 
        
        html.Div(children=[
        dbc.Label("Normalization axis"),         
        dbc.RadioItems(
        id='normalization_factors',   
        options=[{"label": x, "value": x} for x in ['X-axis', 'Y-axis','No']],
        value='No',
        inline=True
        )], style={'width': '20%', 'display': 'inline-block', 'font-size': '12px', 'height':'10vh', 'padding': '0 5'}),  
        
        # html.Div(children=[
        # dbc.Label("Base-count variables"),         
        # dbc.RadioItems(
        # id='base_variable',   
        # options=[{"label": x, "value": x} for x in ['Patient_no', 'HSA_CNT','MFI']],
        # value='100'  
        # )], style={'width': '20%', 'display': 'inline-block', 'font-size': '12px', 'height':'10vh', 'padding': '0 5'}), 
 
        html.Div(children=[     
        dcc.RangeSlider(0, 10000, value=[500, 10000], marks={0: "0", 100: "100", 500: "500", 1000: "1000", 2000: "2000", 5000: "5000", 10000: "10000"}, id='mfi_cutoff')
        ], style={'width': '80%', 'display': 'inline-block', 'font-size': '12px', 'height':'10vh', 'padding': '0 5'}),    

    #     html.Div(children=[
    #     dbc.Label("MFI cutoff(H-L)"),         
    #     dbc.RadioItems(
    #     id='mfi_cutoff',   
    #     options=[{"label": x, "value": x} for x in ['10000-5000' , '10000-1000' , '10000-500' , '5000-1000' ,'5000-500', '1000-500']],
    #     value='10000-5000',
    #     inline=True 
    #     )], style={'width': '80%', 'display': 'inline-block', 'font-size': '12px', 'height':'10vh', 'padding': '0 5'}),    

        html.Div(children=[                        
        dcc.Graph(id="dsa-heatmaps-graph")], 
        style={'width': '100%', 'display': 'inline-block', 'font-size': '12px', 'height':'20vh', 'padding': '0 56'})        
        ])   
    ])

# call the input and output for DSA distributions 
@callback(Output("dsa-heatmaps-graph", "figure"), 
          Input("upper_col", "value"),         
          Input("outter_row", "value"),
          Input("inner_row", "value"),
          Input("lower_col", "value"),  
          Input("normalization_factors", "value"),         
          Input("mfi_cutoff", "value")                                                           
          )
# def update_output(value):
#     return 'You have selected "{}"'.format(value[5])

# # define plot for DSA distributions 
def dsa_filter_heatmap(upper_col_var, outter_row_var, inner_row_var, lower_col_var, normalization_decision, mfi_threshold):
    # Create subplots
    if upper_col_var=='HLA_LOCUS' or outter_row_var=='HLA_LOCUS' or inner_row_var=='HLA_LOCUS' or lower_col_var=='HLA_LOCUS':
        temp_df_dsa_cnt_table_cr=df[df['HLA_LOCI_combined']==1].copy()
    else:
        temp_df_dsa_cnt_table_cr=df[df['HLA_LOCI_combined']==0].copy()       
    
    # reorder the data according to the upper/outter and lower/inner variable selections
    temp_df_dsa_cnt_table_cr=temp_df_dsa_cnt_table_cr[[upper_col_var, inner_row_var, outter_row_var, lower_col_var,temp_df_dsa_cnt_table_cr.columns[7]]]
    # set MFI_GRP by mfi_threshold
    if upper_col_var=='MFI_GRP' or outter_row_var=='MFI_GRP' or inner_row_var=='MFI_GRP' or lower_col_var=='MFI_GRP':   
        mfi_threshold_l=int(mfi_threshold[0])
        mfi_threshold_h=int(mfi_threshold[1])    
        temp_df_dsa_cnt_table_cr['MFI_GRP_new'] = temp_df_dsa_cnt_table_cr['MFI_GRP'].apply(lambda x: 'MFI=<'+str(mfi_threshold_l) if x<=mfi_threshold_l and x!=0 else ( 'MFI='+str(mfi_threshold_l) + '~'+ str(mfi_threshold_h) if x>mfi_threshold_l and x<=mfi_threshold_h else ( 'MFI=>'+str(mfi_threshold_h) if x>mfi_threshold_h else 0)))
        temp_df_dsa_cnt_table_cr['MFI_GRP'] = temp_df_dsa_cnt_table_cr['MFI_GRP_new']
        temp_df_dsa_cnt_table_cr=temp_df_dsa_cnt_table_cr.drop(['MFI_GRP_new'],axis=1)
        temp_df_dsa_cnt_table_cr=temp_df_dsa_cnt_table_cr.loc[temp_df_dsa_cnt_table_cr['MFI_GRP']!=0,:] 
        temp_df_dsa_cnt_table_cr['MFI_GRP'] = pd.Categorical(temp_df_dsa_cnt_table_cr['MFI_GRP'], categories=['MFI=<'+str(mfi_threshold_l),'MFI='+str(mfi_threshold_l) + '~'+ str(mfi_threshold_h),'MFI=>'+str(mfi_threshold_h)], ordered=True)      
    if upper_col_var=='HLA_MATCH8' or outter_row_var=='HLA_MATCH8' or inner_row_var=='HLA_MATCH8' or lower_col_var=='HLA_MATCH8':
        temp_df_dsa_cnt_table_cr['HLA_MATCH8'] = pd.Categorical(temp_df_dsa_cnt_table_cr['HLA_MATCH8'], categories=['HLA_MATCH8=8/8','HLA_MATCH8=7/8','HLA_MATCH8=6/8'], ordered=True)
    if upper_col_var=='HLA_MATCH10' or outter_row_var=='HLA_MATCH10' or inner_row_var=='HLA_MATCH10' or lower_col_var=='HLA_MATCH10':    
        temp_df_dsa_cnt_table_cr['HLA_MATCH10'] = pd.Categorical(temp_df_dsa_cnt_table_cr['HLA_MATCH10'], categories=['HLA_MATCH10=10/10','HLA_MATCH10=9/10','HLA_MATCH10=8/10'], ordered=True)
    if upper_col_var=='HSA_CNTs' or outter_row_var=='HSA_CNTs' or inner_row_var=='HSA_CNTs' or lower_col_var=='HSA_CNTs':
        temp_df_dsa_cnt_table_cr['HSA_CNTs'] = pd.Categorical(temp_df_dsa_cnt_table_cr['HSA_CNTs'], categories=['HSA_CNTs=1','HSA_CNTs=2','HSA_CNTs=3','HSA_CNTs=4','HSA_CNTs=5','HSA_CNTs=6','HSA_CNTs=7','HSA_CNTs=8','HSA_CNTs=9','HSA_CNTs=10'], ordered=True)
    if upper_col_var=='RACE' or outter_row_var=='RACE' or inner_row_var=='RACE' or lower_col_var=='RACE':
        temp_df_dsa_cnt_table_cr['RACE'] = temp_df_dsa_cnt_table_cr['RACE'].apply(lambda x: 'RACE=AFA' if x=='RACE=Black' else ('RACE=API' if x=='RACE=Asian' else ('RACE=HIS' if x=='RACE=Hispanic' else ('RACE=CAU' if x=='RACE=White'else 'RACE=OTH'))))
        temp_df_dsa_cnt_table_cr['RACE'] = pd.Categorical(temp_df_dsa_cnt_table_cr['RACE'], categories=['RACE=AFA','RACE=API','RACE=CAU','RACE=HIS','RACE=OTH'], ordered=True)
    if upper_col_var=='HLA_LOCUS' or outter_row_var=='HLA_LOCUS' or inner_row_var=='HLA_LOCUS' or lower_col_var=='HLA_LOCUS':
        temp_df_dsa_cnt_table_cr['HLA_LOCUS'] = pd.Categorical(temp_df_dsa_cnt_table_cr['HLA_LOCUS'], categories=['HLA_LOCUS=A','HLA_LOCUS=B','HLA_LOCUS=C','HLA_LOCUS=DRB1','HLA_LOCUS=DRB3','HLA_LOCUS=DRB4','HLA_LOCUS=DRB5','HLA_LOCUS=DPA1','HLA_LOCUS=DPB1','HLA_LOCUS=DQA1','HLA_LOCUS=DQB1'], ordered=True)
    temp_df_dsa_cnt_table_cr=pd.DataFrame(temp_df_dsa_cnt_table_cr.groupby([upper_col_var, inner_row_var, outter_row_var, lower_col_var]).sum()).reset_index()
    temp_df_dsa_cnt_table_cr=temp_df_dsa_cnt_table_cr.sort_values([upper_col_var, inner_row_var, outter_row_var, lower_col_var])

    # define the axis dimensions of upper/outter and lower/inner dataframe(subplots) by the levels of the corresponding variables    
    upper_col_idx = len(temp_df_dsa_cnt_table_cr.iloc[:,0].unique())
    inner_row_idx = len(temp_df_dsa_cnt_table_cr.iloc[:,1].unique())
    outter_row_idx = len(temp_df_dsa_cnt_table_cr.iloc[:,2].unique())
    lower_col_idx = len(temp_df_dsa_cnt_table_cr.iloc[:,3].unique())  
    # dsa_cnt_min = temp_df_dsa_cnt_table_cr.iloc[:,4].min()
    # dsa_cnt_max = temp_df_dsa_cnt_table_cr.iloc[:,4].max() 
    # dsa_cnt_colorscale= [
    #         [dsa_cnt_min, 'rgb(68, 1, 84)'],        #0
    #         [round((dsa_cnt_max+dsa_cnt_min)/2), 'rgb(63, 174, 85)'],  #100
    #         [dsa_cnt_max, 'rgb(255, 253, 141)']           #100000
    #     ]
    # dsa_cnt_colorscale = [
    #     [0, 'rgb(68, 1, 84)'],   # Red at the lowest value
    #     [0.5, 'rgb(63, 174, 85)'],   # Yellow in the middle
    #     [1, 'rgb(255, 253, 141)']    # Green at the highest value
    # ]

    # initialize the subplots with the axis dimensions              
    fig = make_subplots(rows=outter_row_idx, cols=upper_col_idx, shared_xaxes=True, shared_yaxes=True,horizontal_spacing=0.002,vertical_spacing=0.006)
    n=0
    # add heatmaps to subplots by looping the upper/outter columns and the upper/outter rows
    for i in range(upper_col_idx):
        for j in range(outter_row_idx):
            # extract subset of dsa data by the indices of the upper/outter columns and the upper/outter rows
            m=np.where((temp_df_dsa_cnt_table_cr.iloc[:,0]==temp_df_dsa_cnt_table_cr.iloc[:,0].unique()[i]) & (temp_df_dsa_cnt_table_cr.iloc[:,2]==temp_df_dsa_cnt_table_cr.iloc[:,2].unique()[j]))
            # reshape the subset dataframe by the axis dimensions of lower/inner columns and the lower/inner rows
            z=temp_df_dsa_cnt_table_cr.iloc[np.r_[m],4].to_numpy().reshape((inner_row_idx,lower_col_idx))
            # z[outter_row_idx-1,0]=dsa_cnt_min
            # z[outter_row_idx-1,lower_col_idx-1]=dsa_cnt_max
            # set normalization decision
            if normalization_decision == 'X-axis':
                row_sums = z.sum(axis=1, keepdims=True)
                row_sums = row_sums if row_sums.all() != 0 else 1
                z = z / row_sums
            elif normalization_decision == 'Y-axis':
                col_sums = z.sum(axis=0, keepdims=True)
                col_sums = col_sums if col_sums.all() != 0 else 1
                z = z / col_sums                
            elif normalization_decision == 'No':
                pass            
            # assign the lower/inner labels of x and y axis by lower/inner columns and the lower/inner rows
            x=temp_df_dsa_cnt_table_cr.iloc[np.r_[m],3].to_numpy().reshape((inner_row_idx,lower_col_idx))[0]
            y=temp_df_dsa_cnt_table_cr.iloc[np.r_[m],1].to_numpy().reshape((inner_row_idx,lower_col_idx)).transpose()[0]
            # add the heatmap subplots using the subset dataframe for each loop
            if n==0:
                heatmap = go.Heatmap(z=z, x=x, y=y,   colorscale = 'Viridis', showscale = True)
            else:
                heatmap = go.Heatmap(z=z, x=x, y=y,   colorscale = 'Viridis', showscale = False)
            fig.add_trace(heatmap, row=j+1, col=i+1) 
            n=n+1  # counter for heatmaps       

    # assign axis titles by the upper/outter columns and the upper/outter rows
    for i in range(upper_col_idx):
        fig.update_xaxes(title_text=temp_df_dsa_cnt_table_cr.iloc[:,0].unique()[i], row=1, col=i+1, side='top')
    for j in range(outter_row_idx):
        fig.update_yaxes(title_text=temp_df_dsa_cnt_table_cr.iloc[:,2].unique()[j], row=j+1,  col=1)
        
    # update the figure size and font size
    fig.update_traces(  textfont_size=10)
    fig.update_layout(  title_text="Exploratory analyse of Clinical Variables Interactions in H(D)SA Data",
                        width=1250,  
                        height=600 )
    return fig
