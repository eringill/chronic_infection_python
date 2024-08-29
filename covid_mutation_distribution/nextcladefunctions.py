# functions to be used in app.py

# imports
import pandas as pd # processing dataframes
import numpy as np # numbers are important!
import re # regex
import math # math is important!
from pathlib import Path
from subprocess import call

def generate_alignment_script():
    ref_seqs = ['wuhan', 'BA2', 'BA286', 'XBB']
    path = Path(__file__).parent / "./data/results/nextcladerun.sh"
    results = Path(__file__).parent / "data/results/"
    test_file = Path(__file__).parent / "./data/reference_seqs/testFASTA/test.fasta"
    path.touch()
    with path.open('w') as file:
        for i in ref_seqs:
            file.write(f'nextclade run {test_file} --output-tsv {results}/{i}results.tsv --input-dataset /Users/egill/Projects/chronic_infection_python/covid_mutation_distribution/data/reference_seqs/{i}/\n')

def generate_alignments():        
    generate_alignment_script()
    path = Path(__file__).parent / "./data/results/nextcladerun.sh"
    with open(path, 'rb') as file:
        script = file.read()
    call(script, shell=True)
    
def get_best_reference():
    generate_alignments()
    ref_seqs = ['wuhan', 'BA2', 'BA286', 'XBB']
    score_list = []
    for i in ref_seqs:
        file = Path(__file__).parent / "data/results/"
        results = pd.read_csv(f'{file}/{i}results.tsv', sep= '\t')
        score = results['alignmentScore'].tolist()
        score_list.append(int(score[0]))
    m = max(score_list)
    best_index = score_list.index(m)
    return ref_seqs[best_index]

print(get_best_reference())