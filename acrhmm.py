"""
#Learn an HMM for Acr Peptides


*   Alphabet: amino acids
*   States: chem_states
"""

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


acrProtDF = pd.read_csv("acrproteins.csv",header=0)
acrseqdf = acrProtDF["Seq"]

def toNums(lis):
  n =[]
  for x in lis:
   n.append(ord(x))
  return n

#make a matrix of seqs
acrs = []
acrs2 = []
for i in range(acrseqdf.shape[0]):
  #print([list(acrseqdf[i])])
  acrs.append([toNums(list(acrseqdf[i]))])
  acrs2.append(toNums(list(acrseqdf[i])))

lengths = [len(x) for x in acrs]

plt.hist(lengths, bins = 20)
plt.show()


def buildAcrHmm():
  acr_general_hmm = hmm.MultinomialHMM(n_components=6,init_params="ste")
  acr_general_hmm.emissionprob_ = np.random.rand(6,20)
  acr_general_hmm.transmat_ = np.random.rand(6,6)
  start_p = [1,0,0,0,0,0]

  acr_general_hmm.startprob_ = start_p

  # concatenate acr proteins to one X for hmm.fit
  pep_input = np.array(list(np.sum(acrseqdf.to_numpy().flatten())))
  le_pep = preprocessing.LabelEncoder()
  le_pep.fit_transform(pep_input)
  acr_general_hmm.fit(le_pep.fit_transform(pep_input).reshape(-1,1))


peps = pd.read_csv("/content/drive/MyDrive/01 University/04 Senior/Computational Genomics/Project/peps.csv")

candidate_seqs = peps.seq_aa
candidate_seqs

#test_input = np.array(list(np.sum(candidate_seqs.to_numpy().flatten())))

#print(test_input)
scores = []
for i in range(candidate_seqs.shape[0]):
  test_input = list(candidate_seqs[i])
  le_pep_test = preprocessing.LabelEncoder()
  le_pep_test.fit_transform(test_input)
  scores.append(acr_general_hmm.score(le_pep_test.fit_transform(test_input).reshape(-1,1)))

peps.insert(loc = 2, column='scores', value=scores)

peps

mask = (peps['seq_aa'].str.len() < 270)
pep_cands = peps.loc[mask]

pep_cands.hist('scores')

mask2 = (pep_cands['scores'] > -180)

pep_cands = pep_cands[mask2]

pep_cands

def get_cand_coords(candsFile):
  coords = []
  for txt in candsFile.seq_name:
        #txt = candsFile.loc[i, 'seq_name']
        a = [int(s) for s in re.findall(r'\b\d+\b', txt)]
        coords.append(list(a))
  return coords

locations = np.asarray(get_cand_coords(pep_cands))