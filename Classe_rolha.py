import matplotlib.pyplot as plt
from tqdm.notebook import tqdm
import itertools
import numpy as np
import scipy as sp
import warnings
warnings.filterwarnings("ignore")

pi = np.pi

class ED():
    def __init__(self, N):
        self.N = N

    def Gera_ind(self):
        labels = []
        for label in itertools.product(range(1, 5), repeat=self.N): #Este loop faz um loop recursivo, que serve para automatizar o processo de escolha do número de sítios
            labels.append(list(label))

        labels = np.array(labels)

        return labels
    
    def Sym(self, labels, S, N_e):
        labelsS = np.zeros((len(labels), len(labels[0])))
        labelsN = np.zeros((len(labels), len(labels[0])))
        for i in range(len(labels)):
            for j in range(self.N):
                match labels[i][j]:
                    case 1:
                        labelsS[i][j] = 0
                        labelsN[i][j] = 0
                    case 2:
                        labelsS[i][j] = 1
                        labelsN[i][j] = 1
                    case 3:
                        labelsS[i][j] = -1
                        labelsN[i][j] = 1
                    case 4:
                        labelsS[i][j] = 0
                        labelsN[i][j] = 2


        NS = []
        for i in range(len(labels)):
            NS.append([sum(labelsN[i]), sum(labelsS[i])])

        NS = np.array(NS)

        states = []
        labels = list(labels)

        for i in range(len(NS)):
            if NS[i][0] == N_e and NS[i][1] == S:
                states.append(list(labels[i]))

        return states
    
    def C(self, sit, spin, state):
        new_state = np.zeros(self.N)

        for i in range(self.N):
            if i == sit:
                match state[sit], spin:
                    case 2, 2:
                        new_state[sit] = 1
                    case 4, 2:
                        new_state[sit] = 3
                    case 4, 3:
                        new_state[sit] = 2
                    case 3, 3:
                        new_state[sit] = 1
                    case 4, 4:
                        new_state[sit] = 1
                    case 3, 2:
                        new_state[sit] = 0
                    case 2, 3:
                        new_state[sit] = 0
                    case 2, 4:
                        new_state[sit] = 0
                    case 3, 4:
                        new_state[sit] = 0
                    case 1:
                        new_state[sit] = 0
            else:
                new_state[i] = state[i]

        return new_state

    def C_dag(self, sit, spin, state):
        new_state1 = np.zeros(self.N)

        for i in range(self.N):
            if i == sit:
                match state[sit], spin:
                    case 1, 2:
                        new_state1[sit] = 2
                    case 3, 2:
                        new_state1[sit] = 4
                    case 1, 3:
                        new_state1[sit] = 3
                    case 2, 3:
                        new_state1[sit] = 4
                    case 1, 4:
                        new_state1[sit] = 4
                    case 4:
                        new_state1[sit] = 0
                    case 2, 2:
                        new_state1[sit] = 0
            else:
                new_state1[i] = state[i]


        return new_state1
    
    def H_mu(self, n_states, mu, labels):

        H1 = np.zeros((n_states, n_states))

        for k in range(self.N):
            for i in range(len(labels)):
                C0 = self.C(k, 2, labels[i])
                Cdag = self.C_dag(k, 2, C0)

                for j in range(len(labels)):
                    count = 0
                    for l in range(self.N):
                        if labels[j][l] == Cdag[l]:
                            count += 1
                    if count == self.N:
                        H1[i][j] += -mu
                C0 = self.C(k, 3, labels[i])
                Cdag = self.C_dag(k, 3, C0)

                for j in range(len(labels)):
                    count = 0
                    for l in range(self.N):
                        if labels[j][l] == Cdag[l]:
                            count += 1
                    if count == self.N:
                        H1[i][j] += -mu
        return H1

    def H_hop(self, n_states, t, labels):
        H2 = np.zeros((n_states, n_states))

        for k in range(self.N-1):
            for i in range(len(labels)):
                C0 = self.C(k, 2, labels[i])
                Cdag = self.C_dag(k+1, 2, C0)

                for j in range(len(labels)):
                    count = 0
                    for l in range(self.N):
                        if labels[j][l] == Cdag[l]:
                            count += 1
                    if count == self.N:
                        H2[i][j] += -t

                C0 = self.C(k, 3, labels[i])
                Cdag = self.C_dag(k+1, 3, C0)

                for j in range(len(labels)):
                    count = 0
                    for l in range(self.N):
                        if labels[j][l] == Cdag[l]:
                            count += 1
                    if count == self.N:
                        H2[i][j] += -t


        for k in range(1, self.N):
            for i in range(len(labels)):
                C0_ = self.C(k, 2, labels[i])
                Cdag_ = self.C_dag(k-1, 2, C0_)

                for j in range(len(labels)):
                    count = 0
                    for l in range(self.N):
                        if labels[j][l] == Cdag_[l]:
                            count += 1
                    if count == self.N:
                        H2[i][j] += -t

                C0_ = self.C(k, 3, labels[i])
                Cdag_ = self.C_dag(k-1, 3, C0_)

                for j in range(len(labels)):
                    count = 0
                    for l in range(self.N):
                        if labels[j][l] == Cdag_[l]:
                            count += 1
                    if count == self.N:
                        H2[i][j] += -t
        return H2
    
    def H_int(self, n_states, U, labels):

        H3 = np.zeros((n_states, n_states))

        for k in range(self.N):
            for i in range(len(labels)):
                C0 = self.C(k, 2, labels[i])
                Cdag = self.C_dag(k, 2, C0)

                C0_ = self.C(k, 3, Cdag)
                Cdag_ = self.C_dag(k, 3, C0_)

                for j in range(len(labels)):
                    count = 0
                    for l in range(self.N):
                        if labels[j][l] == Cdag_[l]:
                            count += 1
                    if count == self.N:
                        H3[i][j] += U

        return H3
    def Config(self, lab, mu, t, U):
        config = []
        for i in range(len(lab)):
            countNe = 0
            counts = 0
            for j in range(self.N):
                match lab[i][j]:
                    case 1:
                        countNe += 0
                        counts += 0
                    case 2:
                        countNe += 1
                        counts += 1
                    case 3:
                        countNe += 1
                        counts += -1
                    case 4:
                        countNe += 2
                        counts += 0
                config.append([countNe, counts])

        config.sort()
        config = list(config for config,_ in itertools.groupby(config))

        List_labels = []

        for conf in config:
            List_labels.append(self.Sym(lab, self.N, conf[1], conf[0]))

        List_H = []

        for lab1 in List_labels:
            H1 = self.H_hop(len(lab1), t, lab1)
            H2 = self.H_int(len(lab1), U, lab1)
            H3 = self.H_mu(len(lab1), mu, lab1)

            HU = H1 + H2 + H3
            List_H.append(HU)

        return np.array(config), np.array(List_labels), np.array(List_H)
    
    def Find_Blocks(self, lab, mu, t, U):

        config, List_labels, List_H = self.Config(lab, mu, t, U)

        eigsL, eigsvL = [], []

        for H in List_H:
            eigs, eigsv = sp.linalg.eigh(H)
            eigsv = np.transpose(eigsv)

            eigsL.append(eigs)
            eigsvL.append(eigsv)

        return np.array(eigsL), np.array(eigsvL)
    
    def Find_GS(eigvals, eigvecs, conf):
        count = 0
        GSE = min(eigvals[0])
        GS = eigvecs[count][0]
        eigp = eigvals[0]
        configGS = conf[count]
        for eigs in eigvals:
            if min(eigs) < GSE:
                GSE = min(eigs)
                GS = eigvecs[count][0]
                eigp = eigs
                configGS = conf[count]
                count += 1

        return configGS, GS, GSE
    
    def Find_ind(confGS, configs):
        confgs = np.abs(confGS)

        confC = confgs - np.array([1, 1])
        confCd = confgs + np.array([1, 1])

        indexes = []

        count = 0
        for conf in configs:
            if (conf==confC).all() == True:
                indexes.append(count)
            elif (conf==confgs).all() == True:
                indexes.append(count)
            elif (conf==confCd).all() == True:
                indexes.append(count)
            count += 1

        indexes = np.array(indexes)

        return indexes
    
    def C_Cdag(labs, site): #Ne -> Ne+1
        Clabels1 = np.array(labs)

        for i in range(len(labs)):
            if labs[i][site] == 2:
                Clabels1[i][site] = 0
            elif labs[i][site] == 4:
                Clabels1[i][site] = 0
            elif labs[i][site] == 3:
                Clabels1[i][site] = 4
            elif labs[i][site] == 1:
                Clabels1[i][site] = 2
            else:
                pass
        return Clabels1
    
    def Op_rep(self, labs, Clabs):
        Cdagrep = []
        for i in range(len(labs)):
            aux = []
            for j in range(len(Clabs)):
                count = 0
                for k in range(self.N):
                    if labs[i][k] == Clabs[j][k]:
                        count += 1
                if count == self.N:
                    aux.append((-1)**(i+1))
                else:
                    aux.append(0)

                Cdagrep.append(aux)

        Cdagrep = np.array(Cdagrep)
        Crep = np.transpose(Cdagrep)

        return Crep, Cdagrep
    
class Mecstat():
    def __init__(self, eigE):
        self.eigE = eigE

    def Part_fun(self, beta):
        Z = sum(np.exp(-beta*self.eigE))
        return Z
    
    def E_mean(self, beta):
        Em = sum(self.eigE*np.exp(-beta*self.eigE))/self.Part_fun(beta)
        E2m = sum((self.eigE**2)*np.exp(-beta*self.eigE))/self.Part_fun(beta)
        Evar = E2m - Em**2
        return Em, E2m, Evar
    
    def Sh_entr(self, beta):
        P = np.exp(-beta*self.eigE)/self.Part_fun(beta)
        Shannon_entropy = -sum(P*np.log(P))
        return Shannon_entropy

    

from scipy import linalg
            
HS1 = ED(3)
labels1 = HS1.Gera_ind()

states = HS1.Sym(labels1, 0, 4)
n_states = len(states)
H = HS1.H_hop(n_states, 0, states) + HS1.H_mu(n_states, 2.5, states) + HS1.H_int(n_states, 5, states)

eig = linalg.eigvalsh(H)

T = np.arange(0.05, 5, 0.1)
beta = 1/T

Stat = Mecstat(eig)

Em, Estd = [], []
S = []

for bet in beta:
    em, _, evar = Stat.E_mean(bet)
    Em.append(em)
    Estd.append(np.sqrt(evar))
    S.append(Stat.Sh_entr(bet)/4)

Em = np.array(Em)
Estd = np.array(Estd)

plt.plot(np.abs(Estd/Em))
plt.show()

#plt.imshow(H)
#plt.colorbar()
#plt.show()
print('Terminei')