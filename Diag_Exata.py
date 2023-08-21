import matplotlib.pyplot as plt
import numpy as np

from Classe_rolha import ED
from StatMec import Mecstat as MS

N = 2
            
HS1 = ED(3)
labels1 = HS1.Gera_ind()

Nt = N
St = N%2

states = HS1.Sym(labels1, St, Nt)
n_states = len(states)
H = HS1.H_hop(n_states, 1, states) + HS1.H_mu(n_states, 2.5, states) + HS1.H_int(n_states, 5, states)

plt.imshow(H)
plt.show()