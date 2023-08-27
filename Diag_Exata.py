import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
from time import time
from scipy.linalg import eigvalsh

from Classe_rolha import ED
from StatMec import Mecstat as MS

t0 = time()

N = 6
            
Nt = N
St = N%2

t = 1
U = 5*t
mu = U/2

HS1 = ED(N)

labels1 = HS1.Gera_ind()
labelsS = HS1.Sym(labels1, St, Nt)

ns = len(labelsS)

H1 = HS1.H_mu(ns, mu, labelsS)
H2 = HS1.H_hop(ns, t, labelsS)
H3 = HS1.H_int(ns, U, labelsS)

H = H1 + H2 + H3

"""
config, List_labels, List_H = HS1.Config(labels1, mu, t, U)

eigsL, eigsvL = HS1.Find_Blocks(labels1, mu, t, U)

configGS, GS, GSE = HS1.Find_GS(eigsL, eigsvL, config)

indexes = HS1.Find_ind(configGS, config)

iGS = indexes[1]
iGSC = indexes[0]
iGSCd = indexes[2]

from math import ceil

site = ceil(N/2) - 1

ClabelsM1 = HS1.C_Cdag(List_labels[iGSC], site) # Ne = 2 -> Ne = 3
ClabelsP1 = HS1.C_Cdag(List_labels[iGS], site) # Ne = 3 -> Ne = 4

E_NS, Z_NS = sp.linalg.eigh(List_H[iGS])
E_NSp1, Z_NSp1 = sp.linalg.eigh(List_H[iGSCd])
E_NSm1, Z_NSm1 = sp.linalg.eigh(List_H[iGSC])

Crep, Cdagrep = HS1.Op_rep(List_labels[iGSCd], ClabelsP1)
CrepP, CdagrepP = HS1.Op_rep(List_labels[iGS], ClabelsM1)
"""

#t, GF = HS1.Time_ev(500, 0.1, Crep, Cdagrep, CrepP, CdagrepP, List_H[iGSCd], List_H[iGS], List_H[iGSC], GS)

#rho, om = HS1.FFTGF(t, GF)

en = sp.linalg.eigvalsh(H)

tf = time()

print()
print(f'Tempo de processamento: {round((tf-t0), 1)}s')
print()

plt.plot(en, 'o')
plt.show()

#en = sp.linalg.eigvalsh(List_H[iGS])


#fig = plt.figure(figsize = (12, 8))

#HS1.plottGF(t, GF)

#plt.plot(om, np.abs(rho))
#for e in en:
  #plt.axvline(e, c = 'green')
#plt.xlim(-15, 15)

#plt.show()
