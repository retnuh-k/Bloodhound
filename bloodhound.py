import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import os
import sys
from sklearn import preprocessing
from hmmlearn import hmm 
import pickle
from acahmms import *

def main():
    if not os.path.exists("aca_struct_hmms.pkl"):
        acaHmms = buildAcaHmms()
    else:
        with open("aca_struct_hmms.pkll", "rb") as file: acaHmms = pickle.load(file)
    
if __name__ == '__main__':
    main()  # next section explains the use of sys.exit
    sys.exit()