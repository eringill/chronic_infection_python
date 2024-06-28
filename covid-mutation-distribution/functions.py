# functions to be used in app.py

# imports
import pandas as pd # processing dataframes
import numpy as np # numbers are important!
import re # regex
import math # math is important!
from pathlib import Path

# for test purposes only
# example_mutation_list = [897, 3431, 7842, 8293, 8393, 11042, 12789, 13339, 15756, 18492, 21608, 21711, 21941, 22032, 22208, 22034, 22295, 22353, 22556, 22770, 22895, 22896, 22898, 22910, 22916, 23009, 23012, 23013, 23018, 23019, 23271, 23423, 23604, 24378, 24990, 25207, 26529, 26610, 26681, 26833, 28958]

# function to parse nucleotide mutation files
def parse_mutation_files(filename):
    '''
    input: name of file to be opened
    
    outputs: mut_list-list of mutations contained in file. If multiple mutations are observed in the same
    location, the location appears multiple times in the list. total_mutations-total number of mutations in file
    '''
    # open the file
    df = pd.read_csv(filename, sep='\t')
    # rename the columns
    df.columns = ['position', 'counts']
    # instantiate the list of mutations
    mut_list = []
    # add each nucleotide position where there is a mutation to the list. 
    # for each genome position x where there is a mutation, the position is added n times, 
    # where n is the integer in the 'counts' column
    for x, y in zip(df.counts.tolist(), df.position.tolist()):
        mut_list.extend([y] * x)
    # calculate the total number of mutations in the list
    total_mutations = sum(df.counts.tolist())
    return mut_list, total_mutations

# function to parse gene files
# gene bins from Wuhan reference sequence NC_045512.2
def parse_gene_files(filename):
    '''
    input: user-selected 'gene' or 'genes_split' as bin size
    
    outputs: genelist-list of nucleotide gene start positions, 
    names-list of gene names
    '''
    # if the user selects "gene" as bin size
    if filename == 'gene':
        gene_data = Path(__file__).parent / "./data/genes.csv"
        df = pd.read_csv(gene_data)
        # make a list of nucleotide gene start coordinates
        genelist = df['start'].tolist()
        # make a list of gene names
        names = df['gene'].tolist()
        names.pop()
    # if the user selects "genes_split" as bin size
    # this option splits the spike protein up into three sections:
    # NTD, RBD and postRBD
    elif filename == 'genes_split':
        split_gene_data = Path(__file__).parent / "./data/genes_split.csv"
        df = pd.read_csv(split_gene_data)
        # make a list of nucleotide gene start coordinates
        genelist = df['start'].tolist()
        # make a list of gene names
        names = df['gene'].tolist()
        names.pop()
    return genelist, names

# function to make bins based on either genes or a specific number of nucleotides,
# depending on what the user selects. Mutation positions are then put into bins.
def make_bins(x, binsize):
    x = list(set(x))
    '''
    inputs: x-list of nucleotide positions where mutations occur, binsize-user-defined bin size
    for plotting and likelihood calculations
    
    outputs: counts-a list of the number of mutations that fall into each bin, bins0- the names
    of the bins
    '''
    # first see if the user has selected an integer bin size
    try:
        int(binsize)
        # if this is the case, make a list of the number of mutations that
        # fall into each bin (based on genome size of 30,000) and then make a list
        # of the centers of each of the bins for plotting 
        # np.histogram gives you the bin edges https://numpy.org/doc/stable/reference/generated/numpy.histogram.html
        counts, bins0 = np.histogram(x, bins=range(1,30001,int(binsize)))
        bins0 = 0.5 * (bins0[:-1] + bins0[1:])
    except ValueError:
        y = [i for i in x if i > 265 and i < 30001]
        if binsize == 'gene':
            # first get the gene start positions and gene names from file using
            # parse_gene_files() function
            genebins, names = parse_gene_files('gene')
        else:
            # first get the gene start positions and gene names from file using
            # parse_gene_files() function
            genebins, names = parse_gene_files('genes_split')
        # then make a list of the number of mutations that fall into each bin (gene)
        counts, bins = np.histogram(y, bins=genebins)
        # the last name is n/a, so remove it from the list of gene names
        bins0 = names
    return counts, bins0

# function to calculate likelihood of user's mutation list belonging to specified distributions
def get_likelihood(existing_bin_counts, test_bin_counts):
    '''
    inputs: existing_bin_counts-number of mutations from existing distribution (global, chronic, deer) 
    that fall into each user-specified bin, 
    test_bin_counts-number of mutations from user's list that fall into each user-specified bin
    
    output: log likelihood of user's list of mutations fitting the existing distribution of mutations
    '''
    # convert both inputs into numpy arrays for easier math operations
    existing_bin_counts = np.array(existing_bin_counts)
    test_bin_counts = np.array(test_bin_counts)
    # perform the likelihood calculation
    return np.sum(np.log(((existing_bin_counts + 1)/np.sum(existing_bin_counts + 1)) ** test_bin_counts))

# function to determine most likely distribution for the user's list of mutations
def most_likely(binsize, global_, global_late, chronic, deer, mutated_nucleotide_list):
    '''
    inputs: binsize-user-selected binsize, global_-list of mutated nucleotide positions in global distribution,
    chronic-list of mutated nucleotide positions in chronic distribution, deer-list of mutated nucleotide positions
    in deer distribution, mutated_nucleotide_list-user-specified list of mutated nucleotide positions
    
    outputs: zipped-a list of tuples containing the likelihood that the user's mutation distribution fits each of the
    existing distributions in the following format [(global_likelihood, 'global'), (global_late_likelihood, 'global_late'),
    (chronic_likelihood, 'chronic'), (deer_likelihood, 'deer')],
    best_fit: the name of the distribution that the user's list of mutations fits best (e.g. 'chronic')
    '''
    # first try to see if user input of mutated nucleotides can be processed
    try:
        # gui accepts input as a string, so it first needs to be split into a list 
        # splits occur wherever there is a comma 
        mutated_nucleotide_list = mutated_nucleotide_list.rstrip(',').rstrip().split(',') 
        # try to remove non-digit characters, then convert each string in list into
        # an integer
        int_nuc_list = [re.sub('\D', '', i) for i in mutated_nucleotide_list]
        digit_nuc_list = [int(i) for i in int_nuc_list]
        mut_nuc_list = [i for i in digit_nuc_list if i > 0 and i < 30001]
    except ValueError:
        # if this fails, return None and exit function - there will be a message printed on the
        # screen prompting the user to enter appropriate input
        int_nuc_list = [re.sub('\D', '', i) for i in mutated_nucleotide_list]
        mut_nuc_list = []
        for i in int_nuc_list:
            try:
                if int(i) < 30001 and int(i) > 0:
                    mut_nuc_list.append(int(i))
                else: 
                    pass 
            except:
                pass
    except:    
        names = ['', '', '', '']
        dummy_likelihoods = ['','','','']
        # zip the two lists together
        dummy_zipped = list(zip(dummy_likelihoods, names))
        dummy_fit = ['','']
        return dummy_zipped, dummy_fit
    
    # if the user's input is processed successfully, split the mutated nucleotide positions
    # into bins
    mut_counts, mut_bins = make_bins(mut_nuc_list, binsize)
    # get bins for global, chronic and deer
    global_counts, global_bins = make_bins(global_,binsize)
    global_late_counts, global_late_bins = make_bins(global_late, binsize)
    chronic_counts, chronic_bins = make_bins(chronic,binsize)
    deer_counts, deer_bins = make_bins(deer,binsize)
    
    # calculate all likelihoods using the number of mutations per bin in the user's input and
    # in existing distributions
    global_likelihood = get_likelihood(global_counts, mut_counts)
    global_late_likelihood = get_likelihood(global_late_counts, mut_counts)
    chronic_likelihood = get_likelihood(chronic_counts, mut_counts)
    deer_likelihood = get_likelihood(deer_counts, mut_counts)
    
    # make a list of all likelihoods
    likelihood_list = [global_likelihood, global_late_likelihood, chronic_likelihood, deer_likelihood]
    # create a matching list of names for the list above
    names = ['global_preVoC', 'global_Omicron', 'chronic', 'deer']
    # zip the two lists together
    zipped = list(zip(likelihood_list, names))
    # find the name of the distribution that best fits the user's input
    best_fit = max(zipped)
    return zipped, best_fit

# function to figure out how many times more likely the best fit distribution is than the default (global)
def times_more_likely(zipped_likelihood_list):
    '''
    input: zipped_likelihood_list-a list of tuples containing the likelihood that the user's mutation distribution 
    fits each of the existing distributions in the following format [(global_likelihood, 'global'), (global_late_likelihood, 'global_late'),
    (chronic_likelihood, 'chronic'), (deer_likelihood, 'deer')]
    
    output: the number of times that the user's mutation distribution is better explained by the best
    fit distribution than the next best fit distribution
    '''
    # first sort the input list so that most likely is first
    sorted_list = sorted(zipped_likelihood_list, key=lambda x: x[0])
    # first unzip the input list, keep only the likelihoods (global, chronic, deer)
    unzipped_names = [j for (i, j) in sorted_list]
    unzipped_nums = [i for (i, j) in sorted_list]
    # convert items in the unzipped number list to floats
    unzipped_nums_float = [float(i) for i in unzipped_nums]
    if unzipped_nums_float[0] == float(0):
        return '', 'Please enter a list of nucleotide positions in order to calculate likelihoods.'
    # return the number of times that the user's mutation distribution is better explained by the best
    # fit distribution than the global distribution
    else:
        return math.exp(unzipped_nums_float[3] - unzipped_nums_float[2]), unzipped_names[2]

# function to select colour palettes for the plot
def select_palette(palette_name):
    '''
    input: palette_name-user-specified palette name
    
    output: colour_list-list of hex codes for colours to use when plotting graph
    '''
    if palette_name == "viridis":
        colour_list = ['#7ad151', '#22a884', '#2a788e', '#414487', '#440154']
    elif palette_name == "inferno":
        colour_list = ['#fca50a', '#dd513a', '#932667', '#420a68', '#000004']
    elif palette_name == "plasma":
        colour_list = ['#fca636', '#e16462', '#b12a90', '#6a00a8', '#0d0887']
    elif palette_name == "seaborn":
        colour_list = ['#0173B2', '#029E73', '#D55E00', '#CC78BC', '#ECE133']
    return colour_list

# function to parse user's input into a list and make sure each entry is unique
def parse_user_input(input):
    '''
    input:
    input - string of comma-separated nucleotide positions where mutations occur (optionally
    flanked by nucleotides and or "ins", "del" or "indel")
    output - list of unique nucleotide positions where mutations occur with trailing commas and
    whitespace removed
    '''
    # first strip trailing commas and whitespace from end of string input
    # then split by commas
    input_split = input.rstrip(',').rstrip().split(',')
    # remove any additional tabs, newlines, returns or whitespace from list
    input_cleaned = [i.strip(' \t\n\r') for i in input_split]
    # remove any list entries that are empty
    input_nonempty = [i for i in input_cleaned if i != '']
    # remove any duplicate entries
    return list(set(input_nonempty))

# function to check if user input matches specified format
def check_for_standard_nucleotides(nuc_list):
    '''
    input:
    nuc_list - list of unique nucleotide positions where mutations occur optionally flanked by 
    uppercase nucleotides and /or "ins", "del" or "indel" with trailing commas and whitespace removed
    output:
        if each element in the list matches the regular expression, True else False
        regex matches the following:
            optionally "ins", "del" or "indel" optionally followed by "A", "C", "T" or "G", 
            followed by a digit between 1 and 30000, optionally followed by "A", "C", "T" or "G", 
            optionally followed by "ins", "del" or "indel" (but ONLY if "ins", "del" or "indel" 
            is NOT at the start of the capture group)
    '''
    # instantiate the regular expression
    pattern = re.compile(r"^(?:INS|DEL|INDEL)?[ACTG]?\d{1,5}[ACTG]?(?:(?<!^)INS|DEL|INDEL)?$")
    # for each item in the user's input list, check to see if it matches the regular expression
    for i in nuc_list:
        if re.match(pattern, i):
            pass
        else:
            # if any items don't match, return False
            return False
    # if all items match, return True
    return True

# function to count the number of transitions and transversions in a list of mutations with
# associated nucleotides
def transition_or_transversion(nuc_pos_list):
    '''
    input:
    nuc_pos_list - string of comma-separated nucleotide positions where mutations occur (optionally
    flanked by nucleotides and or "ins", "del" or "indel")
    output: 
    transitions - count of transitions in user's list of mutations
    transversions - count of transversions in user's list of mutations
    '''
    # first generate list of unique nucleotide positions where mutations occur with trailing commas and
    # whitespace removed
    nuc_pos_list_blank_removed = parse_user_input(nuc_pos_list)
    # replace uracil with thymidine
    nuc_list_standard = [s.replace('u', 'T').replace('U', 'T') for s in nuc_pos_list_blank_removed]
    # convert all alphabetic characters to uppercase
    nuc_list_upper = [i.upper() for i in nuc_list_standard]
    # make sure user's input matches expected pattern
    # if not, instead of returning counts of transitions and transversions, return
    # False, False
    if check_for_standard_nucleotides(nuc_list_upper) == False:
        return False, False
    # remove digits from list of mutations so that we can count transitions and transversions
    nuc_pos_list_no_digits = [re.sub('\d+', '', i) for i in nuc_list_upper]
    # instantiate transition and transversion counts
    transitions = 0
    transversions = 0
    # remove 'INS', 'DEL' and 'INDEL' from list, only keep list items with 2 nucleotides
    nuc_list_parsed = [i for i in nuc_pos_list_no_digits if (len(i) == 2)]
    # iterate through each item in list of mutated nucleotides
    # classify each as a transition or transversion
    for i in nuc_list_parsed:
        if i[0] == 'A':
            if i[1] == 'G':
                transitions += 1
            elif i[1] in ['C', 'T']:
                transversions += 1
        elif i[0] == 'C':
            if i[1] == 'T':
                transitions += 1
            elif i[1] in ['A', 'G']:
                transversions += 1
        elif i[0] == 'G':
            if i[1] == 'A':
                transitions += 1
            elif i[1] in ['C', 'T']:
                transversions += 1
        elif i[0] == 'T':
            if i[1] == 'C':
                transitions += 1
            elif i[1] in ['A', 'G']:
                transversions += 1
    # if there are no transversions, add one to the count to avoid a division by 0 error
    if transversions == 0:
        transversions += 1
    # return the counts of transitions and transversions
    return transitions, transversions

                
# function to count mutations that confer a mutator phenotype
def mut_lineage_parsing(nuc_pos_list):
    '''
    input:
    nuc_pos_list - string of comma-separated nucleotide positions where mutations occur (optionally
    flanked by nucleotides and or "ins", "del" or "indel")
    output:
    mutator_text - list of user-entered mutations conferring mutator phenotype in string format
    potential_mutator_text - list of user-entered mutations potentially conferring mutator phenotype
    in string format
    '''
    try:
        # parse user input to list of unique nucleotide positions where mutations occur with trailing commas and
        # whitespace removed
        nuc_pos_list_parsed = parse_user_input(nuc_pos_list)
        # remove any non-digit characters from list of nucleotide positions
        int_nuc_list = [re.sub('\D', '', i) for i in nuc_pos_list_parsed]
        # convert each list entry to an integer
        mut_nuc_list = [int(i) for i in int_nuc_list]
        # instantiate lists of mutator and potential mutator mutations
        mutator_list = [18155, 18218, 18647]
        potential_mutator_list = [18307, 18308, 18309, 18313, 18314, 18315, 18610, 18611, 18612, 18841, 18842, 18843, 18856, 18857, 18858]
        # check to see if the user has entered a nucleotide position that confers a mutator
        # phenotype
        mutator_exists = set(mut_nuc_list) & set(mutator_list)
        # check to see if the user has entered a nucleotide position that potentially 
        # confers a mutator phenotype
        potential_mutator_exists = set(mut_nuc_list) & set(potential_mutator_list)
        # join entries in list of user's mutator mutations into a single string
        mutator_text = ', '.join(str(e) for e in (list(mutator_exists)))
        # join entries in list of user's potential mutator mutations into a single string
        potential_mutator_text = ', '.join(str(e) for e in (list(potential_mutator_exists)))
    except:
        # if the above block doesn't work, provide dummy text
        mutator_text = ''
        potential_mutator_text = ''
    return mutator_text, potential_mutator_text
    
# function to parse numeric data into scientific notation
def sci_notation(number, sig_fig=2):
    '''
    input:
    number - the number to convert into scientific notation
    sig_fig - the number of significant figures to include in the scientific notation
    output:
    number in scientific notation
    '''
    # format the number using built-in .format function
    # looks like 8.2e+7
    ret_string = "{0:.{1:d}e}".format(number, sig_fig)
    # split the formatted number into two parts (the number itself and the exponent)
    a, b = ret_string.split("e")
    # remove leading "+" and strip leading zeros
    b = int(b)
    return a + " * 10^" + str(b)
        