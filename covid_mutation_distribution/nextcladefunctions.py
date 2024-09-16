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
import subprocess

def generate_alignment_script(input_path):
    path_string = str(input_path)
    path_root = path_string.split('.')[0]
    file_to_copy = Path(__file__).parent / "./nextclade"
    path = Path(__file__).parent / "./data/results/"
    shutil.copy2(file_to_copy, path)
    os.environ['PATH'] += os.pathsep + str(path)
    ref_seqs = ['wuhan', 'BA2', 'BA286', 'XBB']
    script_path = Path(__file__).parent / "./data/results/nextcladerun.sh"
    # results = Path(__file__).parent / "./data/results/"
    test_file = Path(__file__).parent / "./data/reference_seqs/"
    user_file = Path(__file__).parent / input_path
    script_path.touch()
    with script_path.open('w') as file:
        for i in ref_seqs:
            file.write(f'nextclade run {user_file} --output-tsv {path_root}{i}results.tsv --input-dataset {test_file}/{i}/\n')

def generate_alignments():
    subprocess.run(['chmod', '0777', "./data/results/nextclade"])        
    path = Path(__file__).parent / "./data/results/nextcladerun.sh"
    with open(path, 'rb') as file:
        script = file.read()
    call(script, shell=True)

def get_best_reference(input_path):
    generate_alignments()
    path_string = str(input_path)
    path_root = path_string.split('.')[0]
    ref_seqs = ['wuhan', 'BA2', 'BA286', 'XBB']
    score_list = []
    for i in ref_seqs:
        results = pd.read_csv(f'{path_root}{i}results.tsv', sep= '\t')
        score = results['alignmentScore'].tolist()
        try:
            score_list.append(int(score[0]))
        except:
            return "Error"
    m = max(score_list)
    best_index = score_list.index(m)
    return ref_seqs[best_index]

def get_private_mutations(input_path):
    best_reference = get_best_reference(input_path)
    path_string = str(input_path)
    path_root = path_string.split('.')[0]
    if best_reference == "Error":
        return "Error"
    tsv = best_reference + "results.tsv"
    df = pd.read_csv(f'{path_root}{tsv}', sep = '\t')
    return ((f'{df["privateNucMutations.reversionSubstitutions"][0]},{df["privateNucMutations.labeledSubstitutions"][0]},{df["privateNucMutations.unlabeledSubstitutions"][0]}').split(','))

def parse_private_mutations(input_path):
    mutation_list = get_private_mutations(input_path)
    if mutation_list == "Error":
        path = Path(__file__).parent / "data/results"
        os.system('rm -rf %s/*' % path)
        return "Error"
    pattern = r'(\|.*)|(^-.*)'
    fixed_list = []
    for item in mutation_list:
        item = re.sub(pattern, '', item)
        if len(item) > 0:
            fixed_list.append(item)        
    path = Path(__file__).parent / "data/results"
    os.system('rm -rf %s/*' % path)
    return (','.join(fixed_list))
        
def execute_nextclade(input_path):
    generate_alignment_script(input_path)
    return parse_private_mutations(input_path)
#print(execute_nextclade(Path(__file__).parent / '/Users/egill/Desktop/testFASTA/test.fasta'))
