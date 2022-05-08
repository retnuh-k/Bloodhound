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
from acrhmm import *

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

def main():

  if len(sys.argv) != 2:
    print("Command line format: python bloodhound.py <jnet_output_dir>")
    exit()
  
  JNET_folder = sys.argv[1]

  if not os.path.exists("aca_struct_hmms.pkl"):
    acaHmms = buildAcaHmms()
  else:
    with open("aca_struct_hmms.pkl", "rb") as file: acaHmms = pickle.load(file)
  
  if not os.path.exists("acr_general_hmm.pkl"):
    acrHmm = buildAcrHmm()
  else:
    with open("acr_general_hmm.pkl", "rb") as file: acrHmm = pickle.load(file)

  #FIND MOST LIKELY ACA REGION
  likelihoods = np.zeros((acaStructs.shape[0], len(acaHmms)))
  lengths = np.zeros((acaStructs.shape[0], len(acaHmms)))
  a_ind = 0
  h_ind = 0
  for a in acaStructs:
    for h in acaHmms:
      likelihoods[a_ind, h_ind] = (np.exp(h.score(le.transform(a).reshape(-1,1))))
      lengths[a_ind, h_ind] = len(a)
      h_ind += 1
  h_ind = 0
  a_ind += 1

  all_jpred = readJnetInput(JNET_folder)

  likelihoods_sample = np.zeros((len(all_jpred), len(acaHmms)))
  lengths_sample = np.zeros((len(all_jpred), len(acaHmms)))
  j_ind = 0
  h_ind = 0
  for j in all_jpred:
    for h in acaHmms:
      likelihoods_sample[j_ind, h_ind] = (np.exp(h.score(le.transform(j).reshape(-1,1))))
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

  aca_pos = positionsOfInterest(num_greater, JNET_folder)
  




if __name__ == '__main__':
  main()  # next section explains the use of sys.exit
  sys.exit()