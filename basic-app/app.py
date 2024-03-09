from shiny.express import input, render, ui
import plotly.express as px
import plotly.graph_objects as go
from shinywidgets import render_widget
import numpy as np
import functions
import pandas as pd

ui.page_opts(title="SARS-CoV-2 Chronic Infection Calculator")

with ui.sidebar():
    ui.input_select("var", "Select Bin Size", 
                    choices= [int(50), int(500), int(1000)]),

chronic = functions.parse_mutation_files('/Users/egill/Projects/chronic_infection_python/basic-app/data/chronicnucl.tsv')
deer = functions.parse_mutation_files('/Users/egill/Projects/chronic_infection_python/basic-app/data/deernucl.tsv')
global_ = functions.parse_mutation_files('/Users/egill/Projects/chronic_infection_python/basic-app/data/globalnucl.tsv')


@render_widget
def hist1():
    x0 = chronic
    x1 = global_
    x2 = deer

    fig = go.Figure()
    fig.add_trace(go.Histogram(
    x=x0,
    histnorm='percent',
    name='chronic', # name used in legend and hover labels
    xbins=dict(
        start=1,
        end=30001,
        size = input.var()),
    # marker_color='#EB89B5',
    opacity=0.75
    ))
    fig.add_trace(go.Histogram(
    x=x1,
    histnorm='percent',
    name='global',
    xbins=dict(
        start=1,
        end=30001,
        size = input.var()),
    #marker_color='#330C73',
    opacity=0.75
    ))

    fig.add_trace(go.Histogram(
    x=x2,
    histnorm='percent',
    name='deer',
    xbins=dict(
        start=1,
        end=30001,
        size = input.var()),
    #marker_color='#330C73',
    opacity=0.75
    ))
    fig.update_layout(
    title_text='Distribution of Mutations\nAcross Genome', # title of plot
    xaxis_title_text='Genome Position (nt)', # xaxis label
    yaxis_title_text='Normalized Mutation Count', # yaxis label
    bargap=0.2, # gap between bars of adjacent location coordinates
    bargroupgap=0.1 # gap between bars of the same location coordinates
    )
    return fig
    
