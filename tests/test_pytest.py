import pytest
from pathlib import Path
from covid_mutation_distribution import functions

mut_list = [897, 3431, 7842, 8293, 8393, 11042, 12789, 13339, 15756, 18492, 21608, 21711, 21941, 22032, 22208, 22034, 22295, 22353, 22556, 22770, 22895, 22896, 22898, 22910, 22916, 23009, 23012, 23013, 23018, 23019, 23271, 23423, 23604, 24378, 24990, 25207, 26529, 26610, 26681, 26833, 28958]
test_dist = Path(__file__).parent/ "test_dist.tsv"
test_genes = Path(__file__).parent/ "test_genes.csv"

def test_parse_mutation_files_list():
    mut_list, total_mutation = functions.parse_mutation_files(test_dist)
    assert len(mut_list) == 282
    
def test_parse_mutation_files_count():
    mut_list, total_mutation = functions.parse_mutation_files(test_dist)
    assert total_mutation == 282
    
def test_parse_gene_files():
    genelist, names = functions.parse_gene_files('gene')
    assert names[0] == 'nsp1'
    
def test_parse_gene_files_list():
    genelist, names = functions.parse_gene_files('gene')
    assert len(genelist) == len(names) + 1

def test_make_bins_1000():
    counts, bins0 = functions.make_bins(functions.parse_mutation_files(test_dist)[0], 1000)
    assert len(counts) == 30
    
def test_make_bins_equal_bins_counts():
    counts, bins0 = functions.make_bins(functions.parse_mutation_files(test_dist)[0], 1000)
    assert len(counts) == len(bins0)
    
def test_make_bins_gene():
    counts, bins0 = functions.make_bins(functions.parse_mutation_files(test_dist)[0], 'gene')
    assert len(counts) == len(bins0)
    
def test_make_bins_gene_names():
    counts, bins0 = functions.make_bins(functions.parse_mutation_files(test_dist)[0], 'gene')
    assert bins0[-1] == 'ORF10'
    
def test_make_bins_genes_split():
    counts, bins0 = functions.make_bins(functions.parse_mutation_files(test_dist)[0], 'genes_split')
    assert len(counts) == len(bins0)
    
def test_make_bins_genes_split_names():
    counts, bins0 = functions.make_bins(functions.parse_mutation_files(test_dist)[0], 'genes_split')
    assert bins0[-1] == 'ORF10'