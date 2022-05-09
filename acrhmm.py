"""
#Learn an HMM for Acr Peptides


*   Alphabet: amino acids
*   States: chem_states
"""

import numpy as np
import pandas as pd
import re
from sklearn import preprocessing
from hmmlearn import hmm 
from acahmms import *

le_pep = preprocessing.LabelEncoder()
le_pep.fit(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y','X', 'B', 'Z','J', ' '])

def toNums(lis):
  n =[]
  for x in lis:
   n.append(ord(x))
  return n

def build_acr_hmm():
  acrProtDF = pd.read_csv("acrproteins.csv",header=0)
  acrseqdf = acrProtDF["Seq"]
  acr_general_hmm = hmm.MultinomialHMM(n_components=6,init_params="ste")
  acr_general_hmm.emissionprob_ = np.random.rand(6,20)
  acr_general_hmm.transmat_ = np.random.rand(6,6)
  start_p = [1,0,0,0,0,0]

  acr_general_hmm.startprob_ = start_p

  # concatenate acr proteins to one X for hmm.fit
  pep_input = np.array(list(np.sum(acrseqdf.to_numpy().flatten())))
  acr_general_hmm.fit(le_pep.transform(pep_input).reshape(-1,1))
  return acr_general_hmm

def get_cand_coords(candsFile):
  coords = []
  for txt in candsFile.seq_name:
        # Find location in genemark name
        a = [int(s) for s in re.findall(r'\b\d+\b', txt)]
        coords.append(list(a))
  return coords

def get_acr_cands(pep_df, acrHmm):
  """Returns a dataframe with protein seqs above a score determined threshold"""
  seqs = pep_df.seq_aa
  scores = []
  for i in range(seqs.shape[0]):
    s = list(seqs[i])
    le_pep.transform(s)
    scores.append(acrHmm.score(le_pep.transform(s).reshape(-1,1)))

  pep_df.insert(loc = 2, column='scores', value=scores)

  mask = (pep_df['seq_aa'].str.len() < 270) # Discard proteins longer than expected
  pep_cands = pep_df.loc[mask]
  mask2 = (pep_cands['scores'] > -180) # Score threshhold
  pep_cands = pep_cands[mask2]
  return pep_cands