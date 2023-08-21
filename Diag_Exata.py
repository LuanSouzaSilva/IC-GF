import matplotlib.pyplot as plt
import scipy as sp
import numpy as np

from Classe_rolha import ED
from StatMec import Mecstat as MS

N = 4
            
HS1 = ED(N)
labels1 = HS1.Gera_ind()

Nt = N
St = N%2

states = HS1.Sym(labels1, St, Nt)
n_states = len(states)
H1 = HS1.H_hop(n_states, 1, states)
H2 = HS1.H_mu(n_states, 2.5, states)
H3 = HS1.H_int(n_states, 5, states)

H =  H1 + H2 + H3 + (5*N/4)*np.identity(len(H1))

#plt.imshow(H)
eigs = sp.linalg.eigvalsh(H)
plt.plot(eigs, 'o')
plt.show()