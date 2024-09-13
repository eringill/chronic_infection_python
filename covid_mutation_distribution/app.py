# imports
from functools import partial # for adding a navbar
from shiny.express import input, render, ui # interactivity
from shiny import reactive, render # reactivity (i.e. calculations)
from shiny.ui import page_navbar # for adding a navbar
from shiny.types import FileInfo
import plotly.graph_objects as go # graph
from shinywidgets import render_widget # rendering graph
import functions # functions from functions.py
import nextcladefunctions
import re # regex
from pathlib import Path
import faicons
from htmltools import HTML
import os
import time
import subprocess
import shutil
import stat

ui.page_opts(
    title="SMDP: SARS-CoV-2 Mutation Distribution Profiler",
    page_fn=partial(page_navbar, id="page"),
    fillable=True
)
# add dark mode switch to navbar
ui.nav_spacer()
with ui.nav_control():
        ui.input_dark_mode() # << 

css_file = Path(__file__).parent / "css" / "styles.css"

# Name of application tab
with ui.nav_panel("Home"):
    with ui.card():
        with ui.accordion(id="acc", open="Quick Start Information"):  
            with ui.accordion_panel("Quick Start Information"):  
                ui.include_css(css_file)
                ui.markdown(
                '''
                
                <!-- Google tag (gtag.js) -->
                <script async src="https://www.googletagmanager.com/gtag/js?id=G-BRMKPZKYHQ"></script>
                <script>
                window.dataLayer = window.dataLayer || [];
                function gtag(){dataLayer.push(arguments);}
                gtag('js', new Date());

                gtag('config', 'G-BRMKPZKYHQ');
                </script>

                <div>
                <img src=https://drive.google.com/thumbnail?id=10krtHUH8xbWrfkcVLfZTu6JJAFKtDK6b>
                </div>
                <p class="opening_paragraph">
                Given a user-provided set of SARS-CoV-2 nucleotide mutations or genome consensus sequence, this application compares the probability of generating this set from the following four distributions:</p>
                <ul class="unordered_list">
                    <li>Mutations observed during the first nine months of the pandemic (pre-VoC) (<b>global pre-VoC distribution</b>)
                    <li>Mutations observed during the Omicron era (<b>global Omicron distribution</b>)
                    <li>Mutations observed in chronic infections (<b>chronic distribution</b>)
                    <li>Mutations observed in zoonotic spillovers from humans to white-tailed deer (<b>deer distribution</b>)
                </ul>
                <p class="opening_paragraph">In addition, the application will inform the user if the mutation pattern is:</p>
                <ul class="unordered_list">
                    <li>Consistent with <b>molnupiravir use</b> (via examination of the transition:transversion ratio)
                    <li>A <b>mutator lineage</b> (contains a mutation in nsp14/exonuclease that is known to increase the mutation rate of the lineage)
                </ul>
                <p class="opening_paragraph">See <b>Application Notes</b> tab for more information.</p>
                
                '''    
                )

    # layout of columns on first tab
    with ui.layout_columns(col_widths=(4, 8)):
        # first column (or "card")
        with ui.card():
            # variables defined by user input
            # bin size
            with ui.tooltip(id="btn_tooltip", placement="right"):
                ui.input_select("var", "Select Bin Size", 
                    choices= ['genes_split', 'gene', int(500), int(1000)])
                'This is the number and type of segments that the genome will be divided into when plotting mutations and calculating likelihoods.'
            # nucleotide positions where mutations occur - example is shown by default
            ui.input_radio_buttons(
                'var2',
                'Please select a lineage whose mutation distribution you would like to visualize, then click "Submit"',
                {'C897A, G3431T, A7842G, C8293T, G8393A, G11042T, C12789T, T13339C, T15756A, A18492G, ins21608, C21711T, G21941T, T22032C, C22208T, A22034G, C22295A, C22353A, A22556G, G22770A, G22895C, T22896A, G22898A, A22910G, C22916T, del23009, G23012A, C23013A, T23018C, T23019C, C23271T, C23423T, A23604G, C24378T, C24990T, C25207T, A26529C, A26610G, C26681T, C26833T, C28958A':'BA.2.86 (Omicron lineage with chronic-like mutation profile)',
                 'G3692T, G4181T, G5617T, C5822T, C6402T, C6638T, C6990T, C7124T, C7926T, C8733T, G9053T, C9611T, C10029T, G11083T, A11201G, C13665T, T14014G, C16466T, C19011T, C19572T, C20589T, C21618G, C21627T, G22028, T22917G, A23403G, C23604G, G24410A, G24815A, C25427T, C25469T, G25793A, T26767C, C27509T, T27638C, C27752T, T28072, T28092, A28247, A28461G, G28881T, G28916T, G29402T':HTML("Lineage from white-tailed deer (Sample 4205 from (<a href=' https://www.cell.com/iscience/fulltext/S2589-0042(23)02396-9 '>Kotwa et al., 2023</a>)) "),
                 'G4460A, G11071A, G3004A, T724C, C11300T, G22186A, G20493A, C2638T, G9128A, C24133T, C12445T, T25150C, G14743A, G18025A, A22633G, C12789T, G28325A, A6626G, T9007C, A15775G, A1844G, C5621T, G12761A, G22899A, C6606T': HTML('Molnupiravir-induced mutation signature (Patient D from (<a href="https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(23)00393-2/fulltext">Fountain-Jones et al. 2024</a>))'),
                 '1':'I want to enter my own list of lineage-defining mutations',
                 '2':HTML('I want to enter sequence data in FASTA format (lineage-defining mutations calculated via <a href=" https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli/index.html ">NextClade CLI</a>)'),
                 }
            )
            with ui.panel_conditional("input.var2 === '1'"):
                with ui.tooltip(id="cond_tooltip", placement="right"):
                    ui.input_text_area("var4", "Please enter a comma-separated list of the lineage-defining mutations (using genomic nucleotide positions, example shown)", 
                                        "C897A, G3431T, A7842G, C8293T, G8393A, G11042T, C12789T, T13339C, T15756A, A18492G, ins21608, C21711T, G21941T, T22032C, C22208T, A22034G, C22295A, C22353A, A22556G, G22770A, G22895C, T22896A, G22898A, A22910G, C22916T, del23009, G23012A, C23013A, T23018C, T23019C, C23271T, C23423T, A23604G, C24378T, C24990T, C25207T, A26529C, A26610G, C26681T, C26833T, C28958A",autoresize=True,)
                    'Power analyses suggest that a minimum of 10 lineage-defining mutations are needed for accurate results.'
            with ui.panel_conditional("input.var2 === '2'"):
                #with ui.tooltip(id="cond_tooltip2", placement="right"):
                ui.input_file("file1", "Please select a file that contains a SARS-CoV-2 genome consensus sequence (FASTA header required, U must be converted to T)", accept=['.fasta', '.FASTA', '.fa'], multiple = False,)
                #'You must include a SINGLE FASTA header and all U in the sequence should be converted to T.'
            
            @reactive.calc
            def parsed_file():
                err_msg = "Error"
                if not input.file1():
                    return
                if not os.path.exists("./data/results/"):
                    os.makedirs("./data/results/")
                try:
                    with open(input.file1()[0]["datapath"], "r") as f:
                        content = f.read()
                        err_msg += "read"
                    timestr = time.strftime("%m%d-%H%M%S")
                    file_path = Path(__file__).parent / f"./data/results/{timestr}.fasta" #need unique name
                    os.system("chown -R shiny:shiny ./data/results/")
                    err_msg += "own"
                    # change permissions for the directory to 777
                    subprocess.run(['chmod', '0777', "./data/results/"])
                    err_msg += "mod"
                    # make node for file in directory
                    fh1 = os.open (f"./data/results/{timestr}.fasta", os.O_CREAT, 0o777)
                    os.close (fh1)
                    err_msg += "touch"
                    with open(file_path, "w") as f2:
                        err_msg += "open"
                        f2.write(content)
                        err_msg += "write"
                    return content, file_path
                except:
                    return err_msg, err_msg
            
            @reactive.effect
            @reactive.event(input.file1)
            def prompt_submit():
                contents, file_path = parsed_file()
 
            # colour palette
            with ui.tooltip(id="btn_tooltip2", placement="right"):
                ui.input_select("var3", "Select Color Palette",
                    choices= ["plasma", "viridis", "inferno", "seaborn"])
                'You can change the colors of the plot here.'
                
            with ui.tooltip(id="btn_tip_submit", placement="right"):
                ui.input_action_button("submit", "Submit", class_="btn-success")
                'Click here to analyze your list of mutations.'

        # second column (or "card")
        with ui.card():
            private_muts = reactive.value(None)
            
            @reactive.effect
            @reactive.event(input.file1)
            def _():
                results, file_path = parsed_file()
                if results.startswith("Error"):
                    private_muts.set("Error")
                else:
                    private_muts.set(nextcladefunctions.execute_nextclade(file_path))
            @reactive.calc
            def number_of_mutations():
                # return the number of mutations that the user has entered
                if private_muts.get():
                    if private_muts.get() == "Error":
                        return 
                    return len(functions.parse_user_input(private_muts.get()))
                elif input.var2() != '1':
                    return len(functions.parse_user_input(input.var2()))
                elif input.var2() == '1':
                    return len(functions.parse_user_input(input.var4()))

            @render.text
            @reactive.event(input.submit, ignore_none=False)
            def print_mutations():
                if private_muts.get():
                    if private_muts.get() == "Error":
                        private_muts.set(None)
                        return 'Please double check your input to ensure that it includes only numeric nucleotide positions between 1 and 30000 (no commas inside digits) and either zero, one or two of the nucleotides A, C, T, G or U. Optionally, each list item may start OR end with "ins", "del" or "indel". Please enter ONLY the first nucleotide at which an insertion, deletion or indel occurs (e.g. del28248). Do not use "_" characters. If you uploaded a file, make sure that the file contains at least 100 nucleotides and that it is the correct file format.'
                    transitions, transversions = functions.transition_or_transversion(private_muts.get())
                elif input.var2() != '1':
                    transitions, transversions = functions.transition_or_transversion(input.var2())                                   
                elif input.var2() == '1':
                    transitions, transversions = functions.transition_or_transversion(input.var4())
                if transversions == False:
                    return 'Please double check your input to ensure that it includes only numeric nucleotide positions between 1 and 30000 (no commas inside digits) and either zero, one or two of the nucleotides A, C, T, G or U. Optionally, each list item may start OR end with "ins", "del" or "indel". Please enter ONLY the first nucleotide at which an insertion, deletion or indel occurs (e.g. del28248). Do not use "_" characters. If you uploaded a file, make sure that the file contains at least 100 nucleotides and that it is the correct file format.'
                if number_of_mutations() == 1:
                    return f'You have entered {number_of_mutations()} mutation.'
                else:
                    return f'You have entered {number_of_mutations()} mutations.'
            with ui.layout_column_wrap(width=1/2):
                with ui.card():
                    @reactive.calc
                    def get_transition_transversion_ratio():
                                # gui accepts input as a string, so it first needs to be split into a list 
                                # splits occur wherever there is a comma
                                # then pass to function defined in functions.py
                                # to get number of transitions, transversions
                        if private_muts.get():
                            transitions, transversions = functions.transition_or_transversion(private_muts.get())
                        elif input.var2() != '1':
                            transitions, transversions = functions.transition_or_transversion(input.var2())                                   
                        elif input.var2() == '1':
                            transitions, transversions = functions.transition_or_transversion(input.var4())
                        return transitions, transversions
                    
                    with ui.tooltip(id="btn_tooltip3", placement="right"):        
                        @render_widget
                        @reactive.event(input.submit, ignore_none=False)
                                # function to plot transition/transversion ratio heatmap
                        def heatmap():
                            transitions, transversions = get_transition_transversion_ratio()
                            fig2 = go.Figure()
                            config = {'displayModeBar': False}
                            if transversions == False:
                                fig2.update_yaxes(showticklabels=False)
                                fig2.update_xaxes(showticklabels=False)
                                fig2.update_layout(height=150, width=300)
                                fig2.update_layout(title_text='Transition to Transversion Ratio') # title of plot
                                return fig2
                            fig2.add_trace(go.Heatmap(
                                z=[[float(transitions)/float(transversions)]],
                                text=[[f'{float(transitions)/transversions:.2f}']],
                                texttemplate='%{text}',
                                colorscale='RdBu',
                                textfont={'size':20},
                                zmax=16, zmin=0,
                                hovertemplate='Transition - Transversion Ratio: %{z}',
                                colorbar=dict(
                                    title="Ratio",
                                    titleside="top",
                                    tickmode="array",
                                    tickvals=[2, 8, 14],
                                    labelalias={2: "Typical", 14: "Molnupiravir-induced"},
                                    ticks="outside"
                                )))
                            
                            fig2.update_yaxes(showticklabels=False)
                            fig2.update_xaxes(showticklabels=False)
                            fig2.update_layout(height=150, width=300)
                            fig2.update_layout(title_text='Transition to Transversion Ratio') # title of plot
                            
                            # return figure
                            return fig2
                        'The transition:transversion ratio of SARS-CoV-2 is typically ~2:1, while molnupiravir induces a ratio of between 9:1 and 14:1 (Gruber et al. 2024).'
                    ui.markdown(
                            '''
                            <p class="opening_paragraph">(see <b><a href="https://movbranchapp.streamlit.app/">movbranch</a></b> for additional molnupiravir analyses)</p>
                            '''
                            )
                with ui.card():
                    with ui.tooltip(id="btn_tooltip4", placement="right"):
                        with ui.value_box(
                                    showcase=faicons.icon_svg("dna", width="80px"),
                                    theme="bg-gradient-blue-purple"
                                ):
                                    "Changes at known mutator sites:"
                                    
                                    @render.ui
                                    @reactive.event(input.submit, ignore_none=False)
                                    def mut_lineage():
                                        if private_muts.get():
                                            transitions, transversions = functions.transition_or_transversion(private_muts.get())
                                        elif input.var2() != '1':
                                            transitions, transversions = functions.transition_or_transversion(input.var2())                                   
                                        elif input.var2() == '1':
                                            transitions, transversions = functions.transition_or_transversion(input.var4())
                                        if transversions == False:
                                            return ''
                                        if private_muts.get():
                                            if (functions.mut_lineage_parsing(private_muts.get())[0] == '') and (functions.mut_lineage_parsing(private_muts.get())[1] == ''):
                                                return 'NO'
                                            else:
                                                return f'Confirmed: {functions.mut_lineage_parsing(private_muts.get())[0]}'
                                        if input.var2() != '1':
                                            if (functions.mut_lineage_parsing(input.var2())[0] == '') and (functions.mut_lineage_parsing(input.var2())[1] == ''):
                                                return 'NO'
                                            else:
                                                return f'Confirmed: {functions.mut_lineage_parsing(input.var2())[0]}'
                                        elif input.var2() == '1':
                                            if (functions.mut_lineage_parsing(input.var4())[0] == '') and (functions.mut_lineage_parsing(input.var4())[1] == ''):
                                                return 'NO'
                                            else:
                                                return f'Confirmed: {functions.mut_lineage_parsing(input.var4())[0]}'

                                    @render.ui
                                    @reactive.event(input.submit, ignore_none=False)
                                    def potential_mut_lineage():
                                        if private_muts.get():
                                            if (functions.mut_lineage_parsing(private_muts.get())[0] == '') and (functions.mut_lineage_parsing(private_muts.get())[1] == ''):
                                                return ''
                                            else:
                                                return f'Potential: {functions.mut_lineage_parsing(private_muts.get())[1]}'
                                        elif input.var2() != '1':
                                            if (functions.mut_lineage_parsing(input.var2())[0] == '') and (functions.mut_lineage_parsing(input.var2())[1] == ''):
                                                return ''
                                            else:
                                                return f'Potential: {functions.mut_lineage_parsing(input.var2())[1]}'
                                        elif input.var2() == '1':
                                            if (functions.mut_lineage_parsing(input.var4())[0] == '') and (functions.mut_lineage_parsing(input.var4())[1] == ''):
                                                return ''
                                            else:
                                                return f'Potential: {functions.mut_lineage_parsing(input.var4())[1]}'
                        'See Application Notes table for a list of Confirmed and Potential mutator sites.'

            with ui.card():
                # shiny won't use file paths in quotes, you have to use pathlib
                # define mutation distributions and total number of mutations for chronic sequences
                chronic_data = Path(__file__).parent / "./data/chronicnucl.tsv"
                chronic, total_chronic = functions.parse_mutation_files(chronic_data)
                # for deer sequences
                deer_data = Path(__file__).parent / "./data/deernucl.tsv"
                deer, total_deer = functions.parse_mutation_files(deer_data)
                # for global sequences
                global_data = Path(__file__).parent / "./data/globalnucl.tsv"
                global_, total_global = functions.parse_mutation_files(global_data)
                
                # for late global sequences
                globallate_data = Path(__file__).parent / "./data/globallatenucl.tsv"
                global_late, total_lateglobal = functions.parse_mutation_files(globallate_data)

                # once nucleotide positions where mutations occur are entered into the text box, these
                # calculations occur reactively
                @reactive.calc
                # function to parse user mutation data input for plotting
                def plot_user_input():
                    # gui accepts input as a string, so it first needs to be split into a list 
                    # splits occur wherever there is a comma
                    if private_muts.get():
                        mutated_nucleotide_list = functions.parse_user_input(private_muts.get())
                    elif input.var2() != '1':
                        mutated_nucleotide_list = functions.parse_user_input(input.var2())
                    elif input.var2() == '1':
                        mutated_nucleotide_list = functions.parse_user_input(input.var4())
                    # try to remove non-digit characters, then convert each string in list into
                    # a digit
                    try: 
                        int_nuc_list = [re.sub('\D', '', i) for i in mutated_nucleotide_list]
                        digit_nuc_list = [int(i) for i in int_nuc_list]
                        mut_nuc_list = [i for i in digit_nuc_list if (i < 30001 and i > 0)]
                    # if this fails, return some dummy data so the plot is still rendered
                    except ValueError: 
                        int_nuc_list = [re.sub('\D', '', i) for i in mutated_nucleotide_list]
                        digit_nuc_list = []
                        for i in int_nuc_list:
                            try:
                                digit_nuc_list.append(int(i))
                                mut_nuc_list = [i for i in digit_nuc_list if (i < 30001 and i > 0)]
                            except ValueError:
                                return [[0,0,0,0], [1,1,1,1], 1]
                        if len(mut_nuc_list) == 0:
                            return [[0,0,0,0], [1,1,1,1], 1]               
                    except: return [[0,0,0,0], [1,1,1,1], 1]
                    if private_muts.get():
                        transitions, transversions = functions.transition_or_transversion(private_muts.get())
                    elif input.var2() != '1':
                        transitions, transversions = functions.transition_or_transversion(input.var2())                                   
                    elif input.var2() == '1':
                        transitions, transversions = functions.transition_or_transversion(input.var4())
                    if transversions == False:
                        return [[0,0,0,0], [1,1,1,1], 1]
                    # otherwise, parse the list of mutation positions into bins based on the size
                    # specified by the user
                    counts, bins0 = functions.make_bins(mut_nuc_list, input.var())
                    # determine the total number of mutations that are in the list supplied by the
                    # user
                    total_counts = sum(counts)
                    if total_counts == 0:
                        total_counts += 1
                    # return everything
                    
                    return counts, bins0, total_counts

                # plot out mutation distributions
                @render_widget
                @reactive.event(input.submit, ignore_none=False)
                # function to plot graph
                def hist1():
                    # assign x variables
                    x0 = global_
                    x1 = global_late
                    x2 = chronic
                    x3 = deer
                    opacity = 1.0
                    
                    # calculate mutations per user-specified bin size in histogram for each distribution
                    # chronic
                    counts0, bins0 = functions.make_bins(x0,input.var())
                    # global
                    counts1, bins1 = functions.make_bins(x1,input.var())
                    # global late
                    counts2, bins2 = functions.make_bins(x2,input.var())
                    # deer
                    counts3, bins3 = functions.make_bins(x3,input.var(), deer = True)
                    # instatiate figure
                    fig = go.Figure()
                    if private_muts.get():
                        transitions, transversions = functions.transition_or_transversion(private_muts.get())
                    elif input.var2() != '1' or '2':
                        transitions, transversions = functions.transition_or_transversion(input.var2())                                   
                    elif input.var2() == '1':
                        transitions, transversions = functions.transition_or_transversion(input.var4())
                    if transversions == False:
                        return fig
                    # add plot of nucleotide positions specified by user
                    fig.add_trace(go.Bar(
                    x=bins0,
                    y=[x/plot_user_input()[2] for x in plot_user_input()[0]], # normalize bin counts by total number of mutations
                    name='user input', # name used in legend and hover labels,
                    marker_color=functions.select_palette(input.var3())[0], # user specifies colour palette
                    opacity=opacity
                    ))
                    # add plot of early global
                    fig.add_trace(go.Bar(
                    x=bins0,
                    y=[x/total_global for x in counts0], # normalize bin counts by total number of mutations
                    name='global pre-VoC', # name used in legend and hover labels,
                    marker_color=functions.select_palette(input.var3())[1], # user specifies colour palette
                    opacity=opacity
                    ))
                    # add plot of late global distribution
                    fig.add_trace(go.Bar(
                    x=bins0,
                    y=[x/total_lateglobal for x in counts1], # normalize bin counts by total number of mutations
                    name='global Omicron', # name used in legend and hover labels,
                    marker_color=functions.select_palette(input.var3())[2], # user specifies colour palette
                    opacity=opacity
                    ))
                    # add plot of chronic distribution
                    fig.add_trace(go.Bar(
                    x=bins0,
                    y=[x/total_chronic for x in counts2],
                    name='chronic',
                    marker_color=functions.select_palette(input.var3())[3],
                    opacity=opacity    
                    ))
                    # add plot of deer distribution
                    fig.add_trace(go.Bar(
                    x=bins0,
                    y=[x/total_deer for x in counts3], # normalize bin counts by total number of mutations
                    name='deer', # name used in legend and hover labels,
                    marker_color=functions.select_palette(input.var3())[4], # user specifies colour palette
                    opacity=opacity
                    ))
                    fig.update_layout(
                    title_text='Distribution of Mutations\nAcross Genome', # title of plot
                    xaxis_title_text='Genome Position', # xaxis label
                    yaxis_title_text='Proportion of Mutations', # yaxis label
                    bargap=0.2, # gap between bars of adjacent location coordinates
                    bargroupgap=0.1, # gap between bars of the same location coordinates
                    plot_bgcolor='white' # specify white background
                    )
                    fig.update_yaxes( # make y axes and ticks look pretty
                    mirror=True,
                    ticks='outside',
                    showline=True,
                    linecolor='black',
                    gridcolor='lightgrey'
                    )
                    fig.update_xaxes( # make x axes and ticks look pretty
                    mirror=True,
                    ticks='outside',
                    showline=True,
                    linecolor='black',
                    gridcolor='white'
                    )
                    # return figure
                    return fig
            
            with ui.card():
                @render.text
                @reactive.event(input.submit, ignore_none=False)
                def txt():
                    return f'The log likelihoods of your sequence fitting the mutation distributions above are as follows: (higher is better)'
                with ui.layout_column_wrap(width=1/2):
                # once nucleotide positions where mutations occur are entered into the text box, these
                # calculations occur reactively
                    @reactive.calc
                    # function to calculate log likelihoods of user's mutation distribution fitting each 
                    # of the specified mutation distributions
                    def calc_likelihoods():
                        # input user's bin size selection, global mutations, chronic mutations, deer mutations, user's mutations
                        if private_muts.get():
                            likelihood_list, most_likely = functions.most_likely(input.var(), global_, global_late, chronic, deer, private_muts.get())
                        elif input.var2() != '1':
                            likelihood_list, most_likely = functions.most_likely(input.var(), global_, global_late, chronic, deer, input.var2())
                        elif input.var2() == '1':
                            likelihood_list, most_likely = functions.most_likely(input.var(), global_, global_late, chronic, deer, input.var4())
                        # return a list of tuples: [(global_likelihood, 'global'), (global_late_likelihood, 'global_late'),(chronic_likelihood, 'chronic'), (deer_likelihood, 'deer')]
                        # and the name of the distribution that the user's list of mutations fits best (e.g. 'chronic')
                        return likelihood_list, most_likely
                    with ui.card():
                        # print text out for the user

                        with ui.value_box(
                            
                            showcase=faicons.icon_svg("globe", width='50px'),
                            theme=ui.value_box_theme(name = 'pre_voc', fg='white', bg='#e16462'),
                            showcase_layout="left center", 
                            max_height='90px'
                        ):
                            "global pre-VoC"
                            
                            @render.ui
                            @reactive.event(input.submit, ignore_none=False)
                            def txt1():
                            # if reactive calculations have been performed (i.e. likelihoods have been calculated),
                            # display likelihoods, otherwise prompt user to enter a list of mutated nucleotide positions
                                if private_muts.get():
                                    transitions, transversions = functions.transition_or_transversion(private_muts.get())
                                elif input.var2() != '1':
                                    transitions, transversions = functions.transition_or_transversion(input.var2())                                   
                                elif input.var2() == '1':
                                    transitions, transversions = functions.transition_or_transversion(input.var4())
                                if transversions == False:
                                    return ''
                                try:
                                    return f'{calc_likelihoods()[0][0][0]:.2f}'
                                except:
                                    pass
                        # print text out for the user

                        with ui.value_box(
                            
                            showcase=faicons.icon_svg("earth-americas", width="50px"),
                            theme=ui.value_box_theme(name = 'omicron', fg='white', bg='#b12a90'),
                            sshowcase_layout="left center", 
                            max_height='90px'
                        ):
                            "global Omicron"
                            @render.ui
                            @reactive.event(input.submit, ignore_none=False)
                            def txt2():
                                if private_muts.get():
                                    transitions, transversions = functions.transition_or_transversion(private_muts.get())
                                elif input.var2() != '1':
                                    transitions, transversions = functions.transition_or_transversion(input.var2())                                   
                                elif input.var2() == '1':
                                    transitions, transversions = functions.transition_or_transversion(input.var4())
                                if transversions == False:
                                    return ''
                                # if reactive calculations have been performed (i.e. likelihoods have been calculated),
                                # display likelihoods, otherwise don't do anything 
                                try:
                                    return f'{calc_likelihoods()[0][1][0]:.2f}'
                                except:
                                    pass
                    with ui.card():
                        # print text out for the user
                        with ui.value_box(
                            
                            showcase=faicons.icon_svg("head-side-virus", width="50px"),
                            theme=ui.value_box_theme(name = 'chronic', fg='white', bg='#6a00a8'),
                            sshowcase_layout="left center", 
                            max_height='90px'
                        ):
                            "chronic"
                            @render.ui
                            @reactive.event(input.submit, ignore_none=False)
                            def txt3():
                                if private_muts.get():
                                    transitions, transversions = functions.transition_or_transversion(private_muts.get())
                                elif input.var2() != '1':
                                    transitions, transversions = functions.transition_or_transversion(input.var2())                                   
                                elif input.var2() == '1':
                                    transitions, transversions = functions.transition_or_transversion(input.var4())
                                if transversions == False:
                                    return ''
                                # if reactive calculations have been performed (i.e. likelihoods have been calculated),
                                # display likelihoods, otherwise don't do anything 
                                try:
                                    return f'{calc_likelihoods()[0][2][0]:.2f}'
                                except:
                                    pass

                        # print text out for the user
                        with ui.value_box(
                            
                            showcase=faicons.icon_svg("virus-covid", width="50px"),
                            theme=ui.value_box_theme(name = 'deer', fg='white', bg='#0d0887'),
                            showcase_layout="left center", 
                            max_height='90px'
                        ):
                            "deer"            
                            @render.ui
                            @reactive.event(input.submit, ignore_none=False)
                            def txt4():
                                # if reactive calculations have been performed (i.e. likelihoods have been calculated),
                                # display likelihoods, otherwise don't do anything
                                if private_muts.get():
                                    transitions, transversions = functions.transition_or_transversion(private_muts.get())
                                elif input.var2() != '1':
                                    transitions, transversions = functions.transition_or_transversion(input.var2())                                   
                                elif input.var2() == '1':
                                    transitions, transversions = functions.transition_or_transversion(input.var4())
                                if transversions == False:
                                    return ''
                                try:
                                    return f'{calc_likelihoods()[0][3][0]:.2f}'
                                except:
                                    pass
                        
            # print text out for the user
            with ui.value_box(
                showcase=faicons.icon_svg("check", width="50px"),
                theme="blue",
            ):
                "Your sequence best fits the following distribution:"
                @render.ui
                @reactive.event(input.submit, ignore_none=False)
                def txt5():
                    if private_muts.get():
                        transitions, transversions = functions.transition_or_transversion(private_muts.get())
                    elif input.var2() != '1':
                        transitions, transversions = functions.transition_or_transversion(input.var2())                                   
                    elif input.var2() == '1':
                        transitions, transversions = functions.transition_or_transversion(input.var4())
                    if transversions == False:
                        return ''
                    # if reactive calculations have been performed (i.e. likelihoods have been calculated),
                    # display likelihoods, otherwise don't do anything
                    if calc_likelihoods()[0][0][0] == float(0):
                        return ''
                    try:
                        return f'{calc_likelihoods()[1][1].replace("_", " ")}'
                        
                    except:
                        pass
                @render.ui
                @reactive.event(input.submit, ignore_none=False)
                def txt6():
                    if private_muts.get():
                        transitions, transversions = functions.transition_or_transversion(private_muts.get())
                    elif input.var2() != '1':
                        transitions, transversions = functions.transition_or_transversion(input.var2())                                   
                    elif input.var2() == '1':
                        transitions, transversions = functions.transition_or_transversion(input.var4())
                    if transversions == False:
                        return ''
                    try:
                        if int(functions.times_more_likely(calc_likelihoods()[0])[0]) > 99999:     
                            more_likely = functions.sci_notation(functions.times_more_likely(calc_likelihoods()[0])[0], sig_fig=1)
                        else:
                            more_likely = f'{functions.times_more_likely(calc_likelihoods()[0])[0]:.2f}'
                        dist = functions.times_more_likely(calc_likelihoods()[0])[1]
                        if private_muts.get():
                            private_muts.set(None)
                        return f'({more_likely} times more likely than the {dist.replace("_", " ")} distribution.)'
                    except:
                        private_muts.set(None)
                        return f'Please enter a list of nucleotide positions or upload a FASTA file to calculate likelihoods.'
                    
                        
            
    ui.markdown(
        '''
        <p class="footer">If you use this tool, please cite the following: <a href="https://doi.org/10.48550/arXiv.2407.11201"><b>Gill, E.E. et al.</b> SMDP: SARS-CoV-2 Mutation Distribution Profiler for rapid estimation of mutational histories of unusual lineages. <i>arXiv</i> 2024; 2407.11201v2</a></p>
        '''
    )                

                
# name of notes tab 
with ui.nav_panel("Application Notes"):
    # markdown of text to appear on second tab page
    ui.markdown(
'''
### Background
SARS-CoV-2 evolution exhibits a strong clock-like signature with mutational changes accumulating over time, but this pattern is punctuated by “saltational changes”, where lineages appear with a higher number of mutations than expected from their divergence time from other lineages ([Neher (2022)](https://academic.oup.com/ve/article/8/2/veac113/6887176)). Such unusual lineages are thought to reflect long passage times within immunocompromised individuals, sharing many of the same signatures seen in chronic infections ([Harari et al. (2022)](https://www.nature.com/articles/s41591-022-01882-4)). 

When unusual lineages arise, however, it is challenging to know the evolutionary history leading to the observed genomic changes.  Other processes, including passage through animals, ([Bashor et al. 2021](https://www.pnas.org/doi/full/10.1073/pnas.2105253118), [Naderi et al. (2023)](https://elifesciences.org/articles/83685)) mutator lineages with error-prone polymerases ([Takada et al. (2023)](https://doi.org/10.1016/j.isci.2023.106210)), and exposure to mutagens such as molnupiravir ([Gruber et al. (2024)](https://onlinelibrary.wiley.com/doi/10.1002/jmv.29642)), can also leave unusual genomic signatures. 

Given a user-provided set of nucleotide mutations defining an unusual lineage of SARS-CoV-2 or a SARS-CoV-2 genome consensus sequence (from which lineage-defining mutations are derived via the [NextClade CLI](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli/index.html)), this application compares the probability of generating this set from the following four distributions:
- The list of mutations observed during the first nine months of the pandemic, prior to the spread of VoC [Harari et al. (2022)](https://www.nature.com/articles/s41591-022-01882-4). (**global pre-VoC distribution**)
- The list of mutations observed in Omicron-era sequences by Harari et al., included submission dates only up to 25 May 2022. (**global Omicron distribution**)
- The list of mutations compiled from 27 chronic infections of immunocompromised individuals [Harari et al. (2022)](https://www.nature.com/articles/s41591-022-01882-4). (**chronic distribution**)
- The list of mutations inferred from 109 separate zoonotic spillovers from humans to white-tailed deer [Feng et al. (2023)](https://www.nature.com/articles/s41467-023-39782-x). (**deer distribution**)

In the first paper, the authors demonstrate that specific lineage-defining mutation patterns occur in SARS-CoV-2 genomes that are sequenced from chronic infections vs. mutations that occurred in SARS-CoV-2 genomes sequenced around the globe at the start of the pandemic (before the rise of Variants of Concern (VOCs)). They also analyzed lineage-defining mutation patterns in VOCs, and concluded that “mutations in chronic infections are predictive of lineage-defining mutations of VOCs”.

Feng et al. sequenced hundreds of SARS-CoV-2 samples obtained from white-tailed deer in the United States. They observed Alpha, Gamma, Delta and Omicron VOCs and determined that the deer infections arose from a minimum of 109 separate transmission events from humans. In addition, the deer were then able to transmit the virus to each other. Deer infections resulted in three documented human zoonoses. The SARS-CoV-2 virus displayed specific adaptation patterns in deer, which differ from adaptations seen in humans. 

In addition, the app informs the user whether the data contain signals consistent with:
- **Past molnupiravir use:** The transition-to-transversion ratio of mutations is calculated in the focal lineage and compared to a background ratio of ~2:1 for SARS-CoV-2 and to case-control cohort studies indicate a ratio of ~14:1 under molnupiravir treatment ([Gruber et al. (2024)](https://onlinelibrary.wiley.com/doi/10.1002/jmv.29642)). A high ratio may thus suggest past exposure to molnupiravir or a similar factor inducing transitions. A sample molnupiravir-induced mutation distribution is taken from [Fountain-Jones et al. (2024)](https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(23)00393-2/fulltext#:~:text=We%20found%20that%20as%20early,patients%20not%20treated%20with%20molnupiravir)
- **Mutator lineages:** Mutator alleles may contribute to the unusual features of a lineage by increasing the rate and type of mutation. Known mutators have been observed in nsp14 within the ExoN proofreading domain of SARS-CoV-2.  P203L in nsp14 was shown to have an elevated substitution rate in phylogenetic analyses, which was confirmed to double the mutation rate when passaged through hamsters ([Takada et al. (2023)](https://doi.org/10.1016/j.isci.2023.106210)). Sites F60S and C39F in nsp14 were associated with a 22-fold and 6-fold higher substitution rate in phylogenetic analyses ([Mack et al. (2023)](https://link.springer.com/article/10.1186/s12967-020-02344-6)). We considered mutations at sites 39, 60, and 203 in nsp14 to be known mutators and mutations in sites 90, 92, 191, 268, and 273, which fall within the ExoN proofreading domain of nsp14, to be potential mutators.

**Table 1: Mutator Sites.** Known and Potential mutator sites (denoted by “Confirmed” and “Potential” in the “Site Type” column, respectively) are listed in the table below. Known sites have been confirmed experimentally, and the specific amino acid / nucleotide changes leading to mutator phenotypes are shown. Potential sites lie within the ExoN proofreading domain of nsp14 (as shown in Mack et al. 2023). The wild type amino acids, their positions within the mature nsp14 protein, encoding nucleotides and genomic locations are shown for these sites, but changes that would lead to mutator phenotypes have not been confirmed.
<table>
    <tr>
        <th>Gene</th>
        <th>Amino Acid Change</th>
        <th>Nucleotide Change</th>
        <th>Site Type</th>
        <th>Reference</th>
    </tr>
    <tr>
        <td>nsp14 / exonuclease</td>
        <td>C39F</td>
        <td>G18,155T</td>
        <td>Confirmed</td>
        <td>(Mack et al. 2023)</td>
    </tr>
    <tr>
        <td>nsp14 / exonuclease</td>
        <td>F60S</td>
        <td>T18,218C</td>
        <td>Confirmed</td>
        <td>(Takada et al. 2023)</td>
    </tr>
    <tr>
        <td>nsp14 / exonuclease</td>
        <td>P203L</td>
        <td>C18,647T</td>
        <td>Confirmed</td>
        <td>(Mack et al. 2023)</td>
    </tr>
    <tr>
        <td>nsp14 / exonuclease</td>
        <td>D90</td>
        <td>18,307-18,309 (GAT)</td>
        <td>Potential</td>
        <td>(Mack et al. 2023)</td>
    </tr>
    <tr>
        <td>nsp14 / exonuclease</td>
        <td>E92</td>
        <td>18,313-18,315 (GAG)</td>
        <td>Potential</td>
        <td>(Mack et al. 2023)</td>
    </tr>
    <tr>
        <td>nsp14 / exonuclease</td>
        <td>E191</td>
        <td>18,610-18,612 (GAG)</td>
        <td>Potential</td>
        <td>(Mack et al. 2023)</td>
    </tr>
    <tr>
        <td>nsp14 / exonuclease</td>
        <td>H268</td>
        <td>18,841-18,843 (CAT)</td>
        <td>Potential</td>
        <td>(Mack et al. 2023)</td>
    </tr>
    <tr>
        <td>nsp14 / exonuclease</td>
        <td>D273</td>
        <td>18,856-18,858 (GAT)</td>
        <td>Potential</td>
        <td>(Mack et al. 2023)</td>
    </tr>
</table>
<br>

### Application Use
This application accepts a list of comma separated nucleotide positions in a SARS-CoV-2 genome where lineage-defining mutations occur. **Lineage-defining mutations are the subset of mutations in a lineage that have occurred since divergence from the larger SARS-CoV-2 tree.** A list of lineage-defining mutations (the “mutation set”) for [pangolin-designated SARS-CoV-2 lineages](https://en.wikipedia.org/wiki/Phylogenetic_Assignment_of_Named_Global_Outbreak_Lineages) can be found [here](https://github.com/cov-lineages/pango-designation?tab=readme-ov-file). The tool will also accept a FASTA file containing a **SINGLE** SARS-CoV-2 genome consensus sequence. In this case, the [NextClade CLI](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli/index.html) is used to determine lineage-defining mutations (called private mutations in NextClade).

The application determines the likelihood of observing the mutation set as a random draw from each distribution (chronic infection, deer-specific mutations, global (pre-VOC) and global (Omicron era)). The log likelihood of observing the mutation set from each distribution is displayed (in natural log units).

Because the mutational data sets are sparse, the method bins sites across the genome when calculating likelihoods. The user can define the bin of interest: genes, genes splitting the spike protein into regions of interest, genome split into 500 nucleotide windows, or genome split into 1000 nucleotide windows. For a given bin choice, the log-likelihood of drawing the user-defined mutation set from each distribution is calculated from the multinomial distribution as:

<p class="equation"><i>Σ(log(((distribution bin counts + 1) / Σ(distribution bin counts + 1))<sup>user bin counts</sup>))</i></p>

The addition of one to each bin ensures that there are no bins lacking data.

### Notes on Input and Useful Tools
- Your list can be formatted **with** or **without** nucleotide abbreviations. e.g. `C897A, G3431T, A7842G, C8293T,...`  OR `897, 3431, 7842, 8293,...`
- These coordinates MUST be **genomic** coordinates, **not gene** coordinates like `S:G107Y`
- Indels should be reported by including the first position only e.g. `ins21608` or `del28248` **NOT** `ins21608TCATGCCGCTGT` or `del28248_28250`
- If you would like to convert gene coordinates to nucleotide coordinates, try using Theo Sanderson’s [tool](https://codon2nucleotide.theo.io/)
- FASTA files must contain a single sequence with a canonical header (e.g. `>genome_sequence`), have one of the following suffixes: `.FASTA`, `.fasta` or `.fa` and **ALL** U nucleotides must be converted to T before upload

### Additional Information
More details are available in the [arXiv preprint](https://doi.org/10.48550/arXiv.2407.11201).

### Example Mutation Distributions
The example mutation distributions available for analysis on the main page are as follows:
- **BA.2.86:** C897A, G3431T, A7842G, C8293T, G8393A, G11042T, C12789T, T13339C, T15756A, A18492G, ins21608, C21711T, G21941T, T22032C, C22208T, A22034G, C22295A, C22353A, A22556G, G22770A, G22895C, T22896A, G22898A, A22910G, C22916T, del23009, G23012A, C23013A, T23018C, T23019C, C23271T, C23423T, A23604G, C24378T, C24990T, C25207T, A26529C, A26610G, C26681T, C26833T, C28958A
- **White-tailed deer sample:** G3692T, G4181T, G5617T, C5822T, C6402T, C6638T, C6990T, C7124T, C7926T, C8733T, G9053T, C9611T, C10029T, G11083T, A11201G, C13665T, T14014G, C16466T, C19011T, C19572T, C20589T, C21618G, C21627T, G22028, T22917G, A23403G, C23604G, G24410A, G24815A, C25427T, C25469T, G25793A, T26767C, C27509T, T27638C, C27752T, T28072, T28092, A28247, A28461G, G28881T, G28916T, G29402T
- **Molnupiravir-induced signature:** G4460A, G11071A, G3004A, T724C, C11300T, G22186A, G20493A, C2638T, G9128A, C24133T, C12445T, T25150C, G14743A, G18025A, A22633G, C12789T, G28325A, A6626G, T9007C, A15775G, A1844G, C5621T, G12761A, G22899A, C6606T

<p class="footer">If you use this tool, please cite the following: <a href="https://doi.org/10.48550/arXiv.2407.11201"><b>Gill, E.E. et al.</b> SMDP: SARS-CoV-2 Mutation Distribution Profiler for rapid estimation of mutational histories of unusual lineages. <i>arXiv</i> 2024; 2407.11201v2</a></p>
        
'''
    )
          
# name of notes tab 
with ui.nav_panel("Contact"):
    # markdown of text to appear on second tab page
    ui.markdown(
'''
### Acknowledgements
This application was developed by the **Computational Analysis, Modelling and Evolutionary Outcomes** ([CAMEO](https://covarrnet.ca/computational-analysis-modelling-and-evolutionary-outcomes-cameo/)) pillar of Canada's **Coronavirus Variants Rapid Response Network** ([CoVaRR-Net](https://covarrnet.ca/)). Data analysis, code and maintenance of the application are conducted by Erin E. Gill, Sheri Harari, Aijing Feng, Fiona S.L. Brinkman, and Sarah Otto. 

Funding was gratefully provided by CoVaRR-Net, which is supported by Genome Canada, Innovation, Science and Economic Development Canada (ISED) and CIHR (grant #ARR-175622). This project was supported by funding from a CoVaRR-Net Rapid Response Research Grant.

### Feedback, Issues and Feature Requests
We're pleased to accept any feedback you have. You can submit an issue on the issues page of the [GitHub repository](https://github.com/eringill/chronic_infection_python). 

You can also email questions, comments or suggestions to Erin Gill at erin.gill81(at)gmail.com.

<p class="footer">If you use this tool, please cite the following: <a href="https://doi.org/10.48550/arXiv.2407.11201"><b>Gill, E.E. et al.</b> SMDP: SARS-CoV-2 Mutation Distribution Profiler for rapid estimation of mutational histories of unusual lineages. <i>arXiv</i> 2024; 2407.11201v2</a></p>

'''
    )           


    