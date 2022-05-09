import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import sys
from hmmlearn import hmm 
from acahmms import *
from acrhmm import *

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

  # GET TABLE
  aca_loci_df = pd.DataFrame(columns=["loci_start","loci_end"])
  aca_loci_df["loci_start"] = [loci[0] for loci in aca_cands]
  aca_loci_df["loci_end"] = [loci[1] for loci in aca_cands]
  aca_loci_df.sort_values("loci_start",inplace=True)

  acr_loci_df = pd.DataFrame(columns=["loci_start","loci_end"])
  acr_loci_df["loci_start"] = [loci[0] for loci in acr_cand_locs]
  acr_loci_df["loci_end"] = [loci[1] for loci in acr_cand_locs]
  acr_loci_df.sort_values("loci_start",inplace=True)
  
  print(f"Aca candidate loci:\n{aca_loci_df}\n")
  print(f"Acr candidate loci:\n{acr_loci_df}")

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
    ax.plot([int(acr_cand_locs[x][0])-75, int(acr_cand_locs[x][1])+75], [2, 2], color='red', linewidth = 10)

  for x in range(len(aca_cands)):
    ax.plot([int(aca_cands[x][0])-75, int(aca_cands[x][1])+75], [1, 1], color='blue', linewidth = 10)

  plt.ylim([0, 3])
  plt.show()


if __name__ == '__main__':
  main()
  sys.exit()