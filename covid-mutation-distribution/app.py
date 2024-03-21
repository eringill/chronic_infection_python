from functools import partial
from shiny.express import input, render, ui
from shiny import reactive
from shiny.ui import page_navbar
import plotly.graph_objects as go
from shinywidgets import render_widget
import functions
import re

ui.page_opts(
    title="SARS-CoV-2 Chronic Infection Calculator",
    page_fn=partial(page_navbar, id="page"),
    fillable=True
)


with ui.nav_panel("Home"):

    with ui.layout_columns(col_widths=(4, 8)):
        with ui.card():
            ui.input_select("var", "Select Bin Size", 
                    choices= ['genes_split', 'gene', int(500), int(1000)])
            (ui.input_text_area("var2", "Please enter a comma-separated list of nucleotide positions where mutations occur here", ""),)
            ui.input_select("var3", "Select Palette",
                    choices= ["viridis", "inferno", "plasma"])

        with ui.card():
            chronic, total_chronic = functions.parse_mutation_files('covid-mutation-distribution/chronicnucl.tsv')
            deer, total_deer = functions.parse_mutation_files('covid-mutation-distribution/deernucl.tsv')
            global_, total_global = functions.parse_mutation_files('covid-mutation-distribution/globalnucl.tsv')

            @reactive.calc
            def plot_user_input():
                mutated_nucleotide_list = input.var2().split(',')
                try: 
                    int_nuc_list = [re.sub('\D', '', i) for i in mutated_nucleotide_list]
                    mut_nuc_list = [int(i) for i in int_nuc_list]
                except: return [[0,0,0], [1,1,1], 1]
                counts, bins0 = functions.make_bins(mut_nuc_list, input.var())
                total_counts = sum(counts)
                return counts, bins0, total_counts

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
                name='chronic', # name used in legend and hover labels,
                marker_color=functions.select_palette(input.var3())[0],
                opacity=0.75
                ))
                fig.add_trace(go.Bar(
                x=bins0,
                y=[x/total_global for x in counts1], # normalize bin counts by total number of mutations
                name='global', # name used in legend and hover labels,
                marker_color=functions.select_palette(input.var3())[1],
                opacity=0.75
                ))

                fig.add_trace(go.Bar(
                x=bins0,
                y=[x/total_deer for x in counts2], # normalize bin counts by total number of mutations
                name='deer', # name used in legend and hover labels,
                marker_color=functions.select_palette(input.var3())[2],
                opacity=0.75
                ))
                fig.add_trace(go.Bar(
                x=bins0,
                y=[x/plot_user_input()[2] for x in plot_user_input()[0]], # normalize bin counts by total number of mutations
                name='user_input', # name used in legend and hover labels,
                marker_color=functions.select_palette(input.var3())[3],
                opacity=0.75
                ))
                fig.update_layout(
                title_text='Distribution of Mutations\nAcross Genome', # title of plot
                xaxis_title_text='Genome Position', # xaxis label
                yaxis_title_text='Proportion of Mutations', # yaxis label
                bargap=0.2, # gap between bars of adjacent location coordinates
                bargroupgap=0.1, # gap between bars of the same location coordinates
                plot_bgcolor='white'
                )
                fig.update_yaxes(
                mirror=True,
                ticks='outside',
                showline=True,
                linecolor='black',
                gridcolor='lightgrey'
                )
                fig.update_xaxes(
                mirror=True,
                ticks='outside',
                showline=True,
                linecolor='black',
                gridcolor='white'
                )

                return fig
            
            @reactive.calc
            def calc_likelihoods():
                likelihood_list, most_likely = functions.most_likely(input.var(), global_, chronic, deer, input.var2())
                return likelihood_list, most_likely

            @render.text
            def txt():
                return f'The log likelihoods of your sequence fitting the mutation distributions above are as follows:'

            @render.text
            def txt1():
                try:
                    return f'{calc_likelihoods()[0][0][1]}: {calc_likelihoods()[0][0][0]:.2f}'
                except:
                    return f'Please enter a comma-separated list of integer nucleotide positions in the box on the left to see your results.'
            

            @render.text
            def txt2():
                try:
                    return f'{calc_likelihoods()[0][1][1]}: {calc_likelihoods()[0][1][0]:.2f}'
                except:
                    pass

            @render.text
            def txt3():
                try:
                    return f'{calc_likelihoods()[0][2][1]}: {calc_likelihoods()[0][2][0]:.2f}'
                except:
                    pass

            @render.text
            def txt4():
                try:
                    return f'Your sequence best fits the distribution of {calc_likelihoods()[1][1]} mutations. ({functions.times_more_likely(calc_likelihoods()[0])} times more likely.)'
                except:
                    pass

with ui.nav_panel("Application Notes"):
    ui.markdown(
        '''
        # Overview
        This application was developed by the Computational Analysis, Modelling and Evolutionary Outcomes 
        ([CAMEO](https://covarrnet.ca/computational-analysis-modelling-and-evolutionary-outcomes-cameo/)) 
        pillar of Canada's Coronavirus Variants Rapid Response Network ([CoVaRR-Net](https://covarrnet.ca/)). 
        Data analysis, code and maintenance of the application are conducted by Erin E. Gill, Fiona S.L. Brinkman,
        and Sarah Otto. 
        
        # Background
        This application accepts a list of genomic positions where mutations occur vs. the Wuhan strain of SARS-CoV-2
        (reference sequence NC_045512.2 in GenBank).
        
        '''
    )