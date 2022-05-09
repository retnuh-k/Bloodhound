from hmmlearn import hmm
from sklearn import preprocessing
import numpy as np
import pickle
import os

"""ACA secondary structure sequences"""

# NOTE: aca1 and aca2 have known secondary structures, the rest were predicted with jpred
aca1struct = np.array(list("----------HHH----HHHHHHHHHHH---HHHHHHHH---HHHHHHHH------------HHHHHHHHHHHH-----"))
aca2struct = np.array(list("--HHHHHHHHHH----HHHHHHH-----HHHHHHHHH------HHHHHHHHHHHHHHHHHHHHHHHHHHHH--EEEEE----HHHHHHHH----HHHHHHHHHHHHHHHH---EEEE--------"))
aca3struct = np.array(list("------------HHHHHHHHHH----HHHHHHHH---HHHHHHHHH------------HHHHHHHHHH--"))
aca4struct = np.array(list("--HHHHHHHHHHHH------HHHHHHHH-----HHHHHHHH---HHHHHHHHHHHHHHHHHHHHH--"))
aca5struct = np.array(list("--HHHHHHHH----HHHHHHH----HHHHHHHHH--EEEEE-----EEEE----------"))
aca6struct = np.array(list("--HHHHHHHH----HHHHHHHH---HHHHHHHH---------HHHHHHHHHHHHHHHH-------"))
aca7struct = np.array(list("-------------HHHHHHHHHH----HHHHHHHH--HHHHHHH------HHHHHHHHHHHHHH----"))
aca9struct = np.array(list("--EEEEE---EEEEE----HHHHHHHHHH----HHHHHHHH--HHHHHHHHHH--"))
aca10struct = np.array(list("-------HHHHHHHHHH----HHHHHHHH---HHHHHHHH-------HHHHHHHHHHHH------"))
aca11struct = np.array(list("---HHHHHHHHH----HHHHHH-----HHHHHHHH----E-----HHHHHHHHHHHHHH----"))
aca12struct = np.array(list("----HHHHHHHHH-----HHHHHHH----HHHHHHHH---HH-----HHHHHHHHHHHH----EHHHHHHHHHHHHHHHHHHHHH--EEEEE-"))
aca13struct = np.array(list("----HHHHHHHHHHHHH--HHHHHHHHHHH-HHHHHHHH--HHHHHHHHHHHHHHHH--"))
acaStructs = np.array([aca1struct,aca2struct,aca3struct,aca4struct,aca5struct,aca6struct,aca7struct,aca9struct,aca10struct,aca11struct,aca12struct,aca13struct])

# Encodes domain chars to integers, '-'=0, 'E'=1, 'H'=2
domains = ['-','E','H']
le_struct = preprocessing.LabelEncoder()
le_struct.fit(domains)

def readJnetInput(JNET_folder):
  all_jpred = []
  for filename in os.listdir(JNET_folder):
    if filename.endswith("jnet"):
      f = open(JNET_folder + '/' + filename, 'r')
      line1 = f.readline()
      line2 = f.readline()
      temp = list(line2[9:line2.index('\n')])
      count = temp.count(',')
      for c in range(count):
        temp.remove(',')
      all_jpred.append(temp)
  return all_jpred

#helper method from https://stackoverflow.com/questions/1883980/find-the-nth-occurrence-of-substring-in-a-string
def find_nth(haystack, needle, n):
  start = haystack.find(needle)
  while start >= 0 and n > 1:
    start = haystack.find(needle, start+len(needle))
    n -= 1
  return start

def positionsOfInterest(num_greater, folder_dir):
  gene_pos = []
  curr_pos = 0
  for filename in os.listdir(folder_dir):
    if filename.endswith("jnet"):
      if num_greater[curr_pos] == 12:
        first = find_nth(filename, '|', 4)
        second = find_nth(filename, '|', 5)
        third = find_nth(filename, '_', 3)
        start = filename[first+1:second]
        end = filename[second+1:third]
        gene_pos.append((start,end))
      curr_pos += 1
  return gene_pos

def normalizeScore(lengths, scores):
  return np.log(scores)/lengths

def get_domain_tmat(seq):
  probs = np.zeros(shape=(len(domains),len(domains)))
  for i in range(len(seq)-1):
    probs[seq[i],seq[i+1]] += 1
  counts = [np.count_nonzero(seq[:-1]==n) for n in range(len(domains))]
  probs = (probs.T/counts).T
  probs[np.isnan(probs)] = 0
  for j in range(probs.shape[0]):
    if sum(probs[j,:]) == 0: probs[j,j] = 1 
  return probs

def build_aca_hmms():
    acaHmms = []
    for aca in acaStructs:
        a = le_struct.transform(aca)
        acahmm = hmm.MultinomialHMM(n_components=len(domains),init_params="ste")
        acahmm.emissionprob_ = np.eye(3)
        acahmm.transmat_ = get_domain_tmat(a)
        start_p = np.zeros(len(domains))
        start_p[a[0]] = 1
        acahmm.startprob_ = start_p
        acaHmms.append(acahmm)
    with open("aca_struct_hmms.pkl", "wb") as file: pickle.dump(acaHmms, file)
    return acaHmms

def aca_cand_locs(acaHmms, JNET_folder):
  likelihoods = np.zeros((acaStructs.shape[0], len(acaHmms)))
  lengths = np.zeros((acaStructs.shape[0], len(acaHmms)))
  a_ind = 0
  h_ind = 0
  for a in acaStructs:
    for h in acaHmms:
      likelihoods[a_ind, h_ind] = (np.exp(h.score(le_struct.transform(a).reshape(-1,1))))
      lengths[a_ind, h_ind] = len(a)
      h_ind += 1
    h_ind =  0
    a_ind += 1

  all_jpred = readJnetInput(JNET_folder)

  likelihoods_sample = np.zeros((len(all_jpred), len(acaHmms)))
  lengths_sample = np.zeros((len(all_jpred), len(acaHmms)))
  j_ind = 0
  h_ind = 0
  for j in all_jpred:
    for h in acaHmms:
      likelihoods_sample[j_ind, h_ind] = (np.exp(h.score(le_struct.transform(j).reshape(-1,1))))
      lengths_sample[j_ind, h_ind] = len(j)
      h_ind += 1
    h_ind = 0
    j_ind += 1

  normalized_aca = normalizeScore(lengths, likelihoods)
  normalized_sample = normalizeScore(lengths_sample, likelihoods_sample)
  compare_aca = np.diag(normalized_aca)
  comparison = np.zeros(lengths_sample.shape)

  num_greater = np.zeros(normalized_sample.shape[0])
  for x in range(normalized_sample.shape[0]):
    #for l in ghao:
    comparison[x,:] = normalized_sample[x,:] > compare_aca
    num_greater[x] = comparison[x,:].sum()
  
  return positionsOfInterest(num_greater, JNET_folder)