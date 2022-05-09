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
from matplotlib.collections import LineCollection

def main():

  if len(sys.argv) != 3:
    print("Command line format: python bloodhound.py <genemark_csv> <jnet_output_dir>")
    exit()
  
  genemark_csv = sys.argv[1]
  JNET_folder = sys.argv[2]
  peps = pd.read_csv(genemark_csv)

  # BUILD HMMS
  acaHmms = build_aca_hmms()
  acrHmm = build_acr_hmm()

  # GET CANDIDATE ACA REGIONS
  aca_cands = aca_cand_locs(acaHmms, JNET_folder)

  # GET CANDIDATE ACR REGIONS
  acr_cands = get_acr_cands(peps, acrHmm)
  acr_cand_locs = get_cand_coords(acr_cands)

  # PRINT CANDIDATE REGIONS
  print(f"Acr candidate regions:\n{acr_cand_locs}\n")
  print(f"Aca candidate regions:\n{aca_cands}")

  # PLOT REGIONS

  coords = []
  # Get range of input for plot
  for txt in peps.seq_name:
        # Find location in genemark name
        a = [int(s) for s in re.findall(r'\b\d+\b', txt)]
        coords.append(list(a))
  coords = sorted(coords)

  fig, ax = plt.subplots()
  ax.plot([int(coords[0][0]), int(coords[-1][-1])], [1, 1], color='dimgrey', linewidth = 10)
  ax.plot([int(coords[0][0]), int(coords[-1][-1])], [2, 2], color='dimgrey', linewidth = 10)
  ax.set_xticklabels([])
  #Y = 1: ACA
  #Y = 2: ACR

  for x in range(len(acr_cand_locs)):
    ax.plot([acr_cand_locs[x][0]-10000, acr_cand_locs[x][1]+10000], [2, 2], color='red', linewidth = 10)

  for x in range(len(aca_cands)):
    ax.plot([int(aca_cands[x][0]) - 10000, int(aca_cands[x][1]) +10000], [1, 1], color='blue', linewidth = 10)

  plt.ylim([0, 3])
  plt.show()


if __name__ == '__main__':
  main()
  sys.exit()