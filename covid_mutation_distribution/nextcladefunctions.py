# functions to be used in app.py

# imports
import pandas as pd # processing dataframes
import numpy as np # numbers are important!
import re # regex
import math # math is important!
from pathlib import Path
from subprocess import call
import os
import shutil

def generate_alignment_script(path):
    sep = '.'
    path = f"{path}"
    path_stripped = path.split(sep, 1)[0]
    ref_seqs = ['wuhan', 'BA2', 'BA286', 'XBB']
    script_path = Path(__file__).parent / f"{path_stripped}nextcladerun.sh"
    test_file = Path(__file__).parent / "./data/reference_seqs/"
    script_path.touch()
    with script_path.open('w') as file:
        for i in ref_seqs:
            file.write(f'nextclade run {path} --output-tsv {path_stripped}{i}results.tsv --input-dataset {test_file}/{i}/\n')
    return script_path

def generate_alignments(script_path):   
    with open(script_path, 'rb') as file:
        script = file.read()
    call(script, shell=True)

def get_best_reference(script_path):
    sep = 'nextcladerun.sh'
    path_stripped = script_path.split(sep, 1)[0]
    generate_alignments(script_path)
    ref_seqs = ['wuhan', 'BA2', 'BA286', 'XBB']
    score_list = []
    for i in ref_seqs:
        path = Path(__file__).parent / path_stripped
        results = pd.read_csv(f'{path}/{i}results.tsv', sep= '\t')
        score = results['alignmentScore'].tolist()
        try:
            score_list.append(int(score[0]))
        except:
            return "Error"
    m = max(score_list)
    best_index = score_list.index(m)
    return ref_seqs[best_index]

def get_private_mutations(script_path):
    sep = 'nextcladerun.sh'
    script_path = f"{script_path}"
    path_stripped = script_path.split(sep, 1)[0]
    best_reference = get_best_reference(script_path)
    if best_reference == "Error":
        return "Error"
    tsv = path_stripped + best_reference + "results.tsv"
    file = Path(__file__).parent / tsv
    df = pd.read_csv(file, sep = '\t')
    return ((f'{df["privateNucMutations.reversionSubstitutions"][0]},{df["privateNucMutations.labeledSubstitutions"][0]},{df["privateNucMutations.unlabeledSubstitutions"][0]}').split(','))

def parse_private_mutations(script_path):
    mutation_list = get_private_mutations(script_path)
    # sep = 'nextcladerun.sh'
    # path_stripped = script_path.split(sep, 1)[0]
    if mutation_list == "Error":
        # path = Path(__file__).parent / path_stripped
        # os.system('rm -rf %s/*' % path)
        return "Error"
    pattern = r'(\|.*)|(^-.*)'
    fixed_list = []
    for item in mutation_list:
        item = re.sub(pattern, '', item)
        if len(item) > 0:
            fixed_list.append(item)        
    # path = Path(__file__).parent / path_stripped
    # os.system('rm -rf %s/*' % path)
    return (','.join(fixed_list))
        
def execute_nextclade(path):
    return parse_private_mutations(generate_alignment_script(path))
#print(execute_nextclade('/Users/egill/Desktop/testFASTA/test.fasta'))
