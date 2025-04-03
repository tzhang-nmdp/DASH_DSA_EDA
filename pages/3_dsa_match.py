import datetime
import numpy as np
from numpy._typing._array_like import NDArray
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
import json
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
hsa_allele_df = pd.read_csv('snowflake/allele_df_cwid.csv',sep=',')
hsa_antigen_df =  pd.read_csv('snowflake/dsa/Patient_Antibody_Report_original_data_antigen.txt',sep='\t')
hla_save_df = pd.read_csv('snowflake/hla_save/Tao 2024 1106_HLA_SAVE.csv',sep=',')
hla_save_df.fillna('nan',inplace=True)
# dsa_dict=pd.read_json('snowflake/hla_save_d.json', orient='index')
with open('snowflake/hla_save_d_2f.json', "r") as jsonfile:
    dsa_dict = json.load(jsonfile)
HLA_save_dsa_RID=[int(key) for key in dsa_dict.keys()]

# Set the subgrouping / filtering variable list 
# group_category = [ 'SEX', 'RACE', 'HLA_LOCUS','HLA_MATCH10','HLA_MATCH8', 'DSA_CNT', 'MFI_GRP']
# group_category
hsa_allele_df=hsa_allele_df.rename(columns={'SEX_CDE':'SEX','ORG_RACE_GROUP_SHORT_NME':'RACE','HLA_LOCUS_CDE':'HLA'})
hsa_antigen_df=hsa_antigen_df.rename(columns={'SEX_CDE':'SEX','ORG_RACE_GROUP_SHORT_NME':'RACE','HLA_LOCUS_CDE':'HLA'})
hsa_allele_df=hsa_allele_df.sort_values(['SEX','RACE','HLA'])
filter_category = hsa_allele_df['SEX'].apply(lambda x: 'SEX' + "=" + str(x) if x == x else "None").unique().tolist() + hsa_allele_df['RACE'].apply(lambda x: 'RACE' + "=" + str(x) if x == x else "None").unique().tolist() + hsa_allele_df['HLA'].apply(lambda x: 'HLA' + "=" + str(x) if x == x else "None").unique().tolist() 
filter_category

layout=dbc.Container([
# title    
            # html.Div(style={'backgroundColor': colors['background'],'color':'blue', 'margin':'0px', 'padding':'20px'},children="DSA exploratory analyses"),
            html.Hr(),

# x_axis widget (for both main plots)      
            dbc.Row([
                html.Div(children=[ 
            dbc.Label("Selected variable for filtering"),                                          
            dcc.Dropdown(options=[{'label': col, 'value': col} for col in filter_category if col==col and col!='None'],value=None, multi=True, id='Variable_filters')  
    ], style={'width': '20%', 'display': 'inline-block',  'font-size': '12px',  'height':'10vh', 'padding': '0 10'}),  

# x_axis widget (1st side plots)      
                html.Div(children=[ 
            dbc.Label("Selected variable for comparisons"),                                        
            dcc.Dropdown(options=[{'label': col, 'value': col} for col in hla_save_df.columns],value='R Age', multi=False , id='Variable_selected')  
    ], style={'width': '20%', 'display': 'inline-block',  'font-size': '12px',  'height':'10vh', 'padding': '0 10'}),  

# x_axis widget (1st side plots)      
                html.Div(children=[ 
            dbc.Label("Selected variable type for comparisons"),                                        
            dcc.Dropdown(options=['continuous', 'categorical'],value='continuous', multi=False ,id='Variable_type')
    ], style={'width': '20%', 'display': 'inline-block',  'font-size': '12px',  'height':'10vh', 'padding': '0 10'}),   

# x_axis widget (2nd side plots)      
                html.Div(children=[ 
            dbc.Label("Selected sample for dsa distribution"),                                        
            dcc.Dropdown(options=[{'label': col, 'value': col} for col in HLA_save_dsa_RID],value=HLA_save_dsa_RID, multi=True , id='Sample_selected')
    ], style={'width': '40%', 'display': 'inline-block',  'font-size': '12px',  'height':'10vh', 'padding': '0 10'}),        

            ]), 

# plots with DSA and MFI distributions   
        dbc.Row([
                html.Div(children=[
            dcc.Graph(
                id='donor_cnt_plot'
            ),  
        ], style={'width': '60%', 'display': 'inline-block', 'font-size': '12px', 'height':'50vh', 'padding': '0 10'}),   

            html.Div(children=[
        dcc.Graph(
            id='var_comp_plot'
        ),  
        ], style={'width': '40%', 'display': 'inline-block', 'font-size': '12px', 'height':'50vh', 'padding': '0 10'})
            ]),  
        dbc.Row([
            html.Div(children=[        
        dcc.Graph(
            id='dsa_cnt_plot'
        ),        
        ], style={'width': '60%', 'display': 'inline-block', 'font-size': '12px', 'height':'40vh','padding': '0 10'}),

            html.Div(children=[        
        dcc.Graph(
            id='dsa_sample_plot'
        ),        
        ], style={'width': '40%', 'display': 'inline-block', 'font-size': '12px', 'height':'40vh','padding': '0 10'})
         ])
         ])

# call the input and output for DSA distributions 
@callback(
    Output(component_id='donor_cnt_plot', component_property='figure'),
    Input(component_id='Variable_filters', component_property='value') 
)
# define plot for DSA distributions 
# scatter plot for MFI distributions with hla antigen types at sample level
def donor_cnt_plot(Variable_filters): 
    hsa_allele_df_tmp=hsa_allele_df.copy()
    if Variable_filters is not None:
        if isinstance(Variable_filters,list):
            Variable_filter_list=[str(x) for x in Variable_filters]
        else:
            Variable_filter_list=[str(Variable_filters)]
        for Variable_filter in Variable_filter_list:        
            filter_col=Variable_filter.split('=')[0]
            filter_category=Variable_filter.split('=')[1]
            hsa_allele_df_tmp =  hsa_allele_df_tmp[hsa_allele_df_tmp[filter_col]==filter_category]
        Donor_cnt_RID_df=hsa_allele_df_tmp[['RID','DONOR_10_OF_10_GT_75_AV','DONOR_9_OF_10_GT_75_AV','DONOR_8_OF_8_GT_75_AV','DONOR_7_OF_8_GT_75_AV','DONOR_8_OF_8_GT_75_35U_AV','DONOR_7_OF_8_GT_75_35U_AV']].copy()
    else:
        Donor_cnt_RID_df =  hsa_allele_df[['RID','DONOR_10_OF_10_GT_75_AV','DONOR_9_OF_10_GT_75_AV','DONOR_8_OF_8_GT_75_AV','DONOR_7_OF_8_GT_75_AV','DONOR_8_OF_8_GT_75_35U_AV','DONOR_7_OF_8_GT_75_35U_AV']]
    Donor_cnt_RID_df['DSA'] =  Donor_cnt_RID_df['RID'].apply(lambda x: 'Y' if x in HLA_save_dsa_RID else 'N')
    Donor_cnt_RID_df = Donor_cnt_RID_df.drop_duplicates(subset=['RID']).copy()
    Donor_cnt_RID_df=Donor_cnt_RID_df.rename(columns={'DONOR_10_OF_10_GT_75_AV':'HLA10/10','DONOR_9_OF_10_GT_75_AV':'HLA9/10','DONOR_8_OF_8_GT_75_AV':'HLA8/8','DONOR_7_OF_8_GT_75_AV':'HLA7/8','DONOR_8_OF_8_GT_75_35U_AV':'HLA8/8_35U','DONOR_7_OF_8_GT_75_35U_AV':'HLA7/8_35U'})
    Donor_cnt_RID_df.fillna(0,inplace=True)
    colors = {'Y': 'red','N': 'blue'}
    sides = {'Y': 'positive','N': 'negative'}
    donor_types=['HLA10/10','HLA9/10','HLA8/8','HLA7/8','HLA8/8_35U','HLA7/8_35U']
    pointposs={'Y': 0.5,'N': -0.5}
    # Create the scatter plot
    fig = go.Figure()
    i=0
    for donor_type in donor_types:
        for name, group in Donor_cnt_RID_df.groupby('DSA'):
            df_tmp1=Donor_cnt_RID_df[Donor_cnt_RID_df['DSA']==name]
            if i<=1:
                fig.add_trace(go.Violin(x=[donor_type]*df_tmp1.shape[0], y=df_tmp1[donor_type].apply(lambda x: math.log(x,10) if x!=0 else 0), side=sides[name],line_color= colors[name], name='DSA='+name))
            else:
                fig.add_trace(go.Violin(x=[donor_type]*df_tmp1.shape[0], y=df_tmp1[donor_type].apply(lambda x: math.log(x,10) if x!=0 else 0), side=sides[name], line_color= colors[name], showlegend=False))
            # , points='all', pointpos=pointposs[name], jitter=0.1, scalemode='width')) 
            i+=1
    fig.update_traces(meanline_visible=True)   
    fig.update_xaxes(
        # mirror=True,
        # ticks='outside',
        showline=True,
        linecolor='black',
        # gridcolor='lightgrey'
    )
    fig.update_yaxes(
        # mirror=True,
        # ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey'
    ) 
    fig.update_layout(violingap=0, violinmode='overlay', 
                    xaxis_type='category',legend_title_text='nonDSA vs DSA',title='Donor CNT distribution',xaxis_title='Donor type',yaxis_title='Donor counts (log)', plot_bgcolor='white',
                    font=dict(family="Arial , monospace",size=10, color="black"),
                    margin=dict(l=20, r=20, t=40, b=20))
    return fig

# call the input and output for DSA distributions 
@callback(
    Output(component_id='var_comp_plot', component_property='figure'),
    Input(component_id='Variable_selected', component_property='value'),
    Input(component_id='Variable_filters', component_property='value'),    
    Input(component_id='Variable_type', component_property='value')
)
# define plot for MFI distributions 
def var_comp_plot(Variable_sel, Variable_filters, Variable_ty): 
    hsa_allele_df_tmp=hsa_allele_df.copy()
    if Variable_filters is not None:
        if isinstance(Variable_filters,list):
            Variable_filter_list=[str(x) for x in Variable_filters]
        else:
            Variable_filter_list=[str(Variable_filters)]
        for Variable_filter in Variable_filter_list:        
            filter_col=Variable_filter.split('=')[0]
            filter_category=Variable_filter.split('=')[1]
            hsa_allele_df_tmp =  hsa_allele_df_tmp[hsa_allele_df_tmp[filter_col]==filter_category]
        HLA_save_dsa_RID_df=hla_save_df[hla_save_df['Nmdp Rid'].isin(hsa_allele_df_tmp['RID'])].copy()
    else:    
        HLA_save_dsa_RID_df =  hla_save_df.copy()
    HLA_save_dsa_RID_df['DSA'] =  HLA_save_dsa_RID_df['Nmdp Rid'].apply(lambda x: 'Y' if x in HLA_save_dsa_RID else 'N')  
    HLA_save_dsa_RID_df=HLA_save_dsa_RID_df.drop_duplicates(subset=['Nmdp Rid']).copy()

    # Create a histogram for the selected column
    if Variable_ty =='continuous':  # Check if it's numeric
        df_tmp_2=HLA_save_dsa_RID_df
        df_tmp_2['DSA'] = pd.Categorical(df_tmp_2['DSA'], categories=['Y','N'], ordered=True)
        fig = px.violin(df_tmp_2, y=Variable_sel,  color='DSA', box=True, title=f"{Variable_sel} comparison")
    else:
        df_tmp_2=HLA_save_dsa_RID_df.value_counts([Variable_sel,'DSA']).reset_index()
        df_tmp_2['DSA'] = pd.Categorical(df_tmp_2['DSA'], categories=['Y','N'], ordered=True)
        fig=px.bar(df_tmp_2, x=df_tmp_2[Variable_sel], 
                            y=df_tmp_2['count'], color='DSA', title=f"{Variable_sel} comparison")
        
    fig.update_layout(plot_bgcolor='white',legend_title_text='nonDSA vs DSA',
                    font=dict(family="Arial , monospace",size=10, color="black"),
                    margin=dict(l=20, r=20, t=40, b=20))   
    fig.update_xaxes(
        # mirror=True,
        # ticks='outside',
        showline=True,
        linecolor='black',
        # gridcolor='lightgrey'
    )
    fig.update_yaxes(
        # mirror=True,
        # ticks='outside',
        showline=True,
        linecolor='black',
        # gridcolor='lightgrey'
    )  
    return fig


# call the input and output for DSA distributions 
@callback(
    Output(component_id='dsa_cnt_plot', component_property='figure'),
    Input(component_id='Variable_filters', component_property='value')
)

# scatter plot for MFI distributions with hla antigen types at sample level
def dsa_cnt_plot(Variable_filters):
    hsa_allele_df_tmp=hsa_allele_df.copy()
    hsa_antigen_df_tmp=hsa_antigen_df.copy()
    if Variable_filters is not None:
        if isinstance(Variable_filters,list):
            Variable_filter_list=[str(x) for x in Variable_filters]
        else:
            Variable_filter_list=[str(Variable_filters)]
        for Variable_filter in Variable_filter_list:        
            filter_col=Variable_filter.split('=')[0]
            filter_category=Variable_filter.split('=')[1]
            hsa_allele_df_tmp =  hsa_allele_df_tmp.loc[hsa_allele_df_tmp[filter_col]==filter_category]
            hsa_antigen_df_tmp =  hsa_antigen_df_tmp.loc[hsa_antigen_df_tmp[filter_col]==filter_category]            
        hdsa_RID_df=pd.concat([hsa_antigen_df_tmp[['RID','HLA','DONOR_10_OF_10_GT_75_AV','DONOR_9_OF_10_GT_75_AV','DONOR_8_OF_8_GT_75_AV','DONOR_7_OF_8_GT_75_AV','DONOR_8_OF_8_GT_75_35U_AV','DONOR_7_OF_8_GT_75_35U_AV']],hsa_allele_df_tmp[['RID','HLA','DONOR_10_OF_10_GT_75_AV','DONOR_9_OF_10_GT_75_AV','DONOR_8_OF_8_GT_75_AV','DONOR_7_OF_8_GT_75_AV','DONOR_8_OF_8_GT_75_35U_AV','DONOR_7_OF_8_GT_75_35U_AV']]],axis=0)
    else:
        hdsa_RID_df =  pd.concat([hsa_antigen_df[['RID','HLA','DONOR_10_OF_10_GT_75_AV','DONOR_9_OF_10_GT_75_AV','DONOR_8_OF_8_GT_75_AV','DONOR_7_OF_8_GT_75_AV','DONOR_8_OF_8_GT_75_35U_AV','DONOR_7_OF_8_GT_75_35U_AV']],hsa_allele_df[['RID','HLA','DONOR_10_OF_10_GT_75_AV','DONOR_9_OF_10_GT_75_AV','DONOR_8_OF_8_GT_75_AV','DONOR_7_OF_8_GT_75_AV','DONOR_8_OF_8_GT_75_35U_AV','DONOR_7_OF_8_GT_75_35U_AV']]],axis=0) 
    hdsa_RID_df['DSA'] =  hdsa_RID_df['RID'].apply(lambda x: 'Y' if x in HLA_save_dsa_RID else 'N')
    hdsa_RID_df_cnt: pd.Series[int]=hdsa_RID_df.value_counts(['DSA', 'HLA', 'RID']).sort_index()
    hdsa_RID_df_cnt=hdsa_RID_df_cnt.reset_index()
    colors = {'A': 'green','B': 'blue','C' : 'red','DRB1' : 'yellow','DRB3' : 'yellow','DRB4': 'yellow','DRB5': 'yellow','DPA1': 'brown','DPB1': 'brown','DQA1': 'brown','DQB1': 'brown'}
    colors = {'Y': 'red','N': 'blue'}
    sides = {'Y': 'positive','N': 'negative'}
    pointposs={'Y': 0.5,'N': -0.5}
    # Create the scatter plot
    fig = go.Figure()
    for name, group in hdsa_RID_df_cnt.groupby('DSA'):
        df_tmp3=hdsa_RID_df_cnt[hdsa_RID_df_cnt['DSA']==name]
        fig.add_trace(go.Violin(x=df_tmp3['HLA'], y=df_tmp3['count'].apply(lambda x: math.log(x,10)), side=sides[name], line_color= colors[name], name='DSA='+name))# , points='all', pointpos=pointposs[name], jitter=0.1, scalemode='width')) 
    fig.update_traces(meanline_visible=True)   
    fig.update_xaxes(
        # mirror=True,
        # ticks='outside',
        showline=True,
        linecolor='black',
        # gridcolor='lightgrey'
    )
    fig.update_yaxes(
        # mirror=True,
        # ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey'
    ) 
    fig.update_layout(violingap=0, violinmode='overlay', 
                    xaxis_type='category',legend_title_text='nonDSA vs DSA',title='H(D)SA CNT distribution',xaxis_title='HLA Loci',yaxis_title='H(D)SA counts (log)', plot_bgcolor='white',
                    font=dict(family="Arial , monospace",size=10, color="black"),
                    margin=dict(l=20, r=20, t=40, b=20))
    return fig


# call the input and output for DSA distributions 
@callback(
    Output(component_id='dsa_sample_plot', component_property='figure'), 
    Input(component_id='Variable_filters', component_property='value'),   
    Input(component_id='Sample_selected', component_property='value')  
)

def dsa_sample_plot(Variable_filters, Sample_sel):
    hsa_allele_df_tmp=hsa_allele_df.copy()
    if Variable_filters is not None:
        if isinstance(Variable_filters,list):
            Variable_filter_list=[str(x) for x in Variable_filters]
        else:
            Variable_filter_list=[str(Variable_filters)]
        for Variable_filter in Variable_filter_list:        
            filter_col=Variable_filter.split('=')[0]
            filter_category=Variable_filter.split('=')[1]
            hsa_allele_df_tmp =  hsa_allele_df_tmp[hsa_allele_df_tmp[filter_col]==filter_category]
        dsa_sample_df =  hsa_allele_df_tmp.loc[(hsa_allele_df_tmp['RID'].isin(HLA_save_dsa_RID)) & (hsa_allele_df_tmp['ANTIGEN_HAR_allele'].notna()) & (hsa_allele_df_tmp['ANTIGEN_HAR_allele']!='nan'),['RID','SEX','RACE','HLA','ANTIGEN_HAR_allele','ALLELE_HAR','MFI_VALUE_NUM','ciwd3','ciwd2']]
    else:    
        dsa_sample_df =  hsa_allele_df_tmp.loc[(hsa_allele_df_tmp['RID'].isin(HLA_save_dsa_RID)) & (hsa_allele_df_tmp['ANTIGEN_HAR_allele'].notna()) & (hsa_allele_df_tmp['ANTIGEN_HAR_allele']!='nan'),['RID','SEX','RACE','HLA','ANTIGEN_HAR_allele','ALLELE_HAR','MFI_VALUE_NUM','ciwd3','ciwd2']]
    dsa_sample_df=dsa_sample_df.copy()
    dsa_sample_df['DSA']=dsa_sample_df[['RID','ANTIGEN_HAR_allele']].apply((lambda x: 'Y' if x.iloc[1] in dsa_dict[str(x.iloc[0])].split('&') else 'N'),axis=1)
    dsa_sample_df['ciwd2']= dsa_sample_df['ciwd2'].astype(str).apply(lambda x: 'C' if x=="['Common']" else 'WD' if x=="['Well-Documented']" else 'nan').value_counts()
    dsa_sample_df['ciwd3']= dsa_sample_df['ciwd3'].astype(str).apply(lambda x: 'C' if x=="['C']" else 'I' if x=="['I']" else 'WD' if x=="['WD']" else 'nan')
    dsa_sample_df['ALLELE_FREQ']=dsa_sample_df[['ciwd3','ciwd2']].apply((lambda x: x[0] if x[0]!='nan' else x[1] if x[1]!='nan' else 'nan'),axis=1)
    dsa_sample_df['MFI_check']=dsa_sample_df['MFI_VALUE_NUM'].apply(lambda x: 1 if x!='nan' and x==x else 0)
    print(dsa_sample_df)    
    # dsa_sample_df=HLA_save_dsa_RID_df.value_counts(['RID', 'HLA', 'ANTIGEN_HAR_allele','ALLELE_HAR']).sort_index()
    # dsa_sample_df=dsa_sample_df.reset_index()
    colors = {'A': 'green','B': 'blue','C' : 'red','DRB1' : 'yellow','DRB3' : 'yellow','DRB4': 'yellow','DRB5': 'yellow','DPA1': 'brown','DPB1': 'brown','DQA1': 'brown','DQB1': 'brown'}
    sizes={'C':50,'I':20,'WD':10,'nan':5}
    shapes={0:'circle',1:'arrow'}
    dsa={'Y':'red','N':'green'}
    hover_data_list=['RID','SEX','RACE','DSA','ALLELE_HAR','ALLELE_FREQ','MFI_VALUE_NUM']
    # Create a bubble for the selected sample
    # fig = go.Figure()
    if isinstance(Sample_sel,list):
        Sample_sel_list=[int(x) for x in Sample_sel]
    else:
        Sample_sel_list=[int(Sample_sel)]
    dsa_sample_df_tmp=dsa_sample_df[dsa_sample_df['RID'].isin(Sample_sel_list)].copy()
    dsa_sample_df_tmp=dsa_sample_df_tmp.fillna('nan')
    # dsa_sample_df_tmp=dsa_sample_df_tmp.rename(columns={'HLA':'HLA'})
    dsa_sample_df_tmp['RID'] = dsa_sample_df_tmp['RID'].astype(str)
    dsa_sample_df_tmp['HLA'] = pd.Categorical(dsa_sample_df_tmp['HLA'], categories=['A','B','C','DRB1','DRB3','DRB4','DRB5','DPA1','DPB1','DQA1','DQB1'], ordered=True)
    dsa_sample_df_tmp2=dsa_sample_df_tmp.sort_values(['HLA','DSA'],ascending=True).copy()
    dsa_sample_df_tmp2=dsa_sample_df_tmp2.reset_index(drop=True)
    fig = px.scatter(dsa_sample_df_tmp2,x="RID", y="ALLELE_HAR",facet_col="HLA", color=dsa_sample_df_tmp2['DSA'], size=dsa_sample_df_tmp2['ALLELE_FREQ'].apply(lambda x: sizes[x]),hover_data=hover_data_list)
    fig.update_layout(xaxis_type='category', yaxis_type='category', legend_title_text='DSA_loci',title='DSA distribution',xaxis_title='index',yaxis_title='DSA allele',font=dict(family="Arial , monospace",size=10, color="black"), plot_bgcolor='white')
    fig.update_xaxes(
        # mirror=True,
        # ticks='outside',
        showticklabels=False,
        title=None,
        showline=True,
        linecolor='black',
        # gridcolor='lightgrey'
    )
    fig.update_yaxes(
        # mirror=True,
        # ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='lightgrey'
    )   
    fig.update_traces(
    marker=dict(
    symbol=[shapes[x] for x in dsa_sample_df_tmp2['MFI_check']]
    ),
    selector=dict(mode="markers"),
    )
    for annotation in fig['layout']['annotations']:
        annotation['textangle'] = 315  
    return fig
