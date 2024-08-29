# functions to be used in app.py

# imports
import pandas as pd # processing dataframes
import numpy as np # numbers are important!
import re # regex
import math # math is important!
from pathlib import Path

def generate_alignments():
    ref_seqs = ['wuhan', 'BA2', 'BA286', 'XBB']
    for i in ref_seqs:
        