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
                    choices= [int(50), int(500), int(1000), 'gene', 'genes_split'])

chronic, total_chronic = functions.parse_mutation_files('/Users/egill/Projects/chronic_infection_python/basic-app/data/chronicnucl.tsv')
deer, total_deer = functions.parse_mutation_files('/Users/egill/Projects/chronic_infection_python/basic-app/data/deernucl.tsv')
global_, total_global = functions.parse_mutation_files('/Users/egill/Projects/chronic_infection_python/basic-app/data/globalnucl.tsv')


@render_widget
def hist1():
    x0 = chronic
    x1 = global_
    x2 = deer

    counts0, bins0 = functions.make_bins(x0,input.var())
    counts1, bins1 = functions.make_bins(x1,input.var())
    counts2, bins2 = functions.make_bins(x2,input.var())
    fig = go.Figure()
    fig.add_trace(go.Bar(
    x=bins0,
    y=[x/total_chronic for x in counts0], # normalize bin counts by total number of mutations
    # histnorm='percent',
    name='chronic', # name used in legend and hover labels,
    # marker_color='#EB89B5',
    opacity=0.75
    ))
    fig.add_trace(go.Bar(
    x=bins0,
    y=[x/total_global for x in counts1], # normalize bin counts by total number of mutations
    # histnorm='percent',
    name='global', # name used in legend and hover labels,
    # marker_color='#EB89B5',
    opacity=0.75
    ))

    fig.add_trace(go.Bar(
    x=bins0,
    y=[x/total_deer for x in counts2], # normalize bin counts by total number of mutations
    # histnorm='percent',
    name='deer', # name used in legend and hover labels,
    # marker_color='#EB89B5',
    opacity=0.75
    ))
    fig.update_layout(
    title_text='Distribution of Mutations\nAcross Genome', # title of plot
    xaxis_title_text='Genome Position', # xaxis label
    yaxis_title_text='Proportion of Mutations', # yaxis label
    bargap=0.2, # gap between bars of adjacent location coordinates
    bargroupgap=0.1 # gap between bars of the same location coordinates
    )
    return fig
    
