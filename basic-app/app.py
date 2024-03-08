from shiny.express import input, render, ui
import plotly.express as px
from shinywidgets import render_plotly

ui.page_opts(title="SARS-CoV-2 Chronic Infection Calculator")

with ui.sidebar():
    ui.input_select("binsize", "Select Bin Size", 
                    choices= [50, 500, 'gene', 'gene - split spike']),


