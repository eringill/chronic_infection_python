# functions to be used in app.py

# imports
import pandas as pd # processing dataframes
import numpy as np # numbers are important!
import re # regex
import math # math is important!

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
        df = pd.read_csv('covid-mutation-distribution/genes.csv')
        # make a list of nucleotide gene start coordinates
        genelist = df['start'].tolist()
        # make a list of gene names
        names = df['gene'].tolist()
    # if the user selects "genes_split" as bin size
    # this option splits the spike protein up into three sections:
    # NTD, RBD and postRBD
    elif filename == 'genes_split':
        df = pd.read_csv('covid-mutation-distribution/genes_split.csv')
        # make a list of nucleotide gene start coordinates
        genelist = df['start'].tolist()
        # make a list of gene names
        names = df['gene'].tolist()
    return genelist, names

# function to make bins based on either genes or a specific number of nucleotides,
# depending on what the user selects. Mutation positions are then put into bins.
def make_bins(x, binsize):
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
    # if the user has specified 'gene' or 'genes_split' instead
    except ValueError:
        if binsize == 'gene':
            # first get the gene start positions and gene names from file using
            # parse_gene_files() function
            genebins, names = parse_gene_files('gene')
            # then make a list of the number of mutations that fall into each bin (gene)
            counts, bins = np.histogram(x, bins=genebins)
            # the last name is n/a, so remove it from the list of gene names
            names.pop()
            bins0 = names
        else:
            # first get the gene start positions and gene names from file using
            # parse_gene_files() function
            genebins, names = parse_gene_files('genes_split')
            # then make a list of the number of mutations that fall into each bin (gene)
            counts, bins = np.histogram(x, bins=genebins)
            # the last name is n/a, so remove it from the list of gene names
            names.pop()
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
def most_likely(binsize, global_, chronic, deer, mutated_nucleotide_list):
    '''
    inputs: binsize-user-selected binsize, global_-list of mutated nucleotide positions in global distribution,
    chronic-list of mutated nucleotide positions in chronic distribution, deer-list of mutated nucleotide positions
    in deer distribution, mutated_nucleotide_list-user-specified list of mutated nucleotide positions
    
    outputs: zipped-a list of tuples containing the likelihood that the user's mutation distribution fits each of the
    existing distributions in the following format [('global', global_likelihood), ('chronic', chronic_likelihood), ('deer', deer_likelihood)],
    best_fit: the name of the distribution that the user's list of mutations fits best (e.g. 'chronic')
    '''
    # first try to see if user input of mutated nucleotides can be processed
    try:
        # gui accepts input as a string, so it first needs to be split into a list 
        # splits occur wherever there is a comma 
        mutated_nucleotide_list = mutated_nucleotide_list.split(',') 
        # try to remove non-digit characters, then convert each string in list into
        # an integer
        int_nuc_list = [re.sub('\D', '', i) for i in mutated_nucleotide_list]
        mut_nuc_list = [int(i) for i in int_nuc_list]
    except:
        # if this fails, return None and exit function - there will be a message printed on the
        # screen prompting the user to enter appropriate input
        return None
    
    # if the user's input is processed successfully, split the mutated nucleotide positions
    # into bins
    mut_counts, mut_bins = make_bins(mut_nuc_list, binsize)
    # get bins for global, chronic and deer
    global_counts, global_bins = make_bins(global_,binsize)
    chronic_counts, chronic_bins = make_bins(chronic,binsize)
    deer_counts, deer_bins = make_bins(deer,binsize)
    
    # calculate all likelihoods using the number of mutations per bin in the user's input and
    # in existing distributions
    global_likelihood = get_likelihood(global_counts, mut_counts)
    chronic_likelihood = get_likelihood(chronic_counts, mut_counts)
    deer_likelihood = get_likelihood(deer_counts, mut_counts)
    
    # make a list of all likelihoods
    likelihood_list = [global_likelihood, chronic_likelihood, deer_likelihood]
    # create a matching list of names for the list above
    names = ['global', 'chronic', 'deer']
    # zip the two lists together
    zipped = list(zip(likelihood_list, names))
    # find the name of the distribution that best fits the user's input
    best_fit = max(zipped)
    return zipped, best_fit

# function to figure out how many times more likely the best fit distribution is than the default (global)
def times_more_likely(zipped_likelihood_list):
    '''
    input: zipped_likelihood_list-a list of tuples containing the likelihood that the user's mutation distribution 
    fits each of the existing distributions in the following format [('global', global_likelihood), 
    ('chronic', chronic_likelihood), ('deer', deer_likelihood)]
    
    output: the number of times that the user's mutation distribution is better explained by the best
    fit distribution than the global distribution
    '''
    # first unzip the input list, keep only the likelihoods (global, chronic, deer)
    unzipped = [i for (i, j) in zipped_likelihood_list]
    # convert items in the unzipped list to floats
    unzipped_nums = [float(i) for i in unzipped]
    # return the number of times that the user's mutation distribution is better explained by the best
    # fit distribution than the global distribution
    return math.exp(max(unzipped_nums) - unzipped_nums[0])

# function to select colour palettes for the plot
def select_palette(palette_name):
    '''
    input: palette_name-user-specified palette name
    
    output: colour_list-list of hex codes for colours to use when plotting graph
    '''
    if palette_name == "viridis":
        colour_list = ['#5ec962', '#21918c', '#3b528b', '#440154']
    elif palette_name == "inferno":
        colour_list = ['#f98e09', '#bc3754', '#57106e', '#000004']
    elif palette_name == "plasma":
        colour_list = ['#f89540', '#cc4778', '#7e03a8', '#0d0887']
    return colour_list