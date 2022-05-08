from hmmlearn import hmm
from sklearn import preprocessing
import numpy as np
import pickle

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
le = preprocessing.LabelEncoder()
le.fit(domains)

def getDomainTmat(seq):
  print("test")
  probs = np.zeros(shape=(len(domains),len(domains)))
  for i in range(len(seq)-1):
    probs[seq[i],seq[i+1]] += 1
  counts = [np.count_nonzero(seq[:-1]==n) for n in range(len(domains))]
  probs = (probs.T/counts).T
  probs[np.isnan(probs)] = 0
  for j in range(probs.shape[0]):
    if sum(probs[j,:]) == 0: probs[j,j] = 1 
  return probs

def buildAcaHmms():
    acaHmms = []
    for aca in acaStructs:
        a = le.transform(aca)
        acahmm = hmm.MultinomialHMM(n_components=len(domains),init_params="ste")
        acahmm.emissionprob_ = np.eye(3)
        acahmm.transmat_ = getDomainTmat(a)
        start_p = np.zeros(len(domains))
        start_p[a[0]] = 1
        acahmm.startprob_ = start_p
        acaHmms.append(acahmm)
    with open("aca_struct_hmms.pkl", "wb") as file: pickle.dump(acaHmms, file)
    return acaHmms