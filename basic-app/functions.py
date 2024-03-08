import pandas as pd

# function to parse nucleotide mutation files
def parse_mutation_files(filename):
    df = pd.read_csv(filename, sep='\t')
    return df

# function to parse gene files
def parse_gene_files(filename):
    df = pd.read_csv(filename)
    return df

# function to make bins
def make_bins(binsize):
    if binsize == 50:
        bins = [b for b in range(1, 30001, 50)]
    elif binsize == 500:
        bins = [b for b in range(1, 30001, 500)]
    elif binsize == 'gene':
        df = parse_gene_files("data/genes.csv")
        bins = df['start'].tolist()
    else:
        df = parse_gene_files("data/genes_split.csv")
        bins = df['start'].tolist()
    return bins

# function to put mutations in bins
