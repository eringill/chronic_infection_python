import pandas as pd

# function to parse nucleotide mutation files
def parse_mutation_files(filename):
    df = pd.read_csv(filename, sep='\t')
    df.columns = ['position', 'counts']
    mut_list = []
    for x, y in zip(df.counts.tolist(), df.position.tolist()):
        mut_list.extend([y] * x)
    return mut_list

# function to make mutation dataframe
def mutation_df():
    chronic = parse_mutation_files('/Users/egill/Projects/chronic_infection_python/basic-app/data/chronicnucl.tsv')
    deer = parse_mutation_files('/Users/egill/Projects/chronic_infection_python/basic-app/data/deernucl.tsv')
    global_mut = parse_mutation_files('/Users/egill/Projects/chronic_infection_python/basic-app/data/globalnucl.tsv')
    
    mut_dict = {'chronic': chronic,
                'deer': deer,
                'global': global_mut}
    
    df = pd.DataFrame.from_dict(mut_dict, orient='index').T.reset_index(drop=True)
    return df


# function to parse gene files
def parse_gene_files(filename):
    df = pd.read_csv(filename)
    return df

# function to make bins
def make_bins(binsize):
    if binsize == 50:
        bins = [b for b in range(1, 30001, 50)]
        bins = 0.5 * (bins[:-1] + bins[1:])
    elif binsize == 500:
        bins = [b for b in range(1, 30001, 500)]
    elif binsize == 'gene':
        df = parse_gene_files("/Users/egill/Projects/chronic_infection_python/basic-app/data/genes_split.csv")
        bins = df['start'].tolist()
    else:
        df = parse_gene_files("/Users/egill/Projects/chronic_infection_python/basic-app/data/genes_split.csv")
        bins = df['start'].tolist()
    return bins

# function to put mutations in bins
def get_counts(df):
    pass