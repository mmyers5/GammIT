import model_plotato_2 as mplo2
import numpy as np
import matplotlib.pyplot as plt

def read_emcee(infile):
  with open(infile) as f:
    fileList = f.readlines()

  z_masses = np.zeros((len(fileList),2))
  for i,fileName in enumerate(fileList):
    result = mplo2.Result(fileName.strip('\n'))
    theta_max = result.theta_max()
    z_masses[i,1] = theta_max[0]
    z_masses[i,0] = result.obs['zred']

  return z_masses
