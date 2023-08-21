import numpy as np

class Mecstat():
    def __init__(self, eigE):
        self.eigE = eigE

    def Part_fun(self, beta): #Calcula a função de particao através de um dado beta = 1/T e dos niveis de energia
        Z = sum(np.exp(-beta*self.eigE))
        return Z
    
    def E_mean(self, beta):# Calcula a media da energia, da energia quadratica e a variancia da energia
        Em = sum(self.eigE*np.exp(-beta*self.eigE))/self.Part_fun(beta)
        E2m = sum((self.eigE**2)*np.exp(-beta*self.eigE))/self.Part_fun(beta)
        Evar = E2m - Em**2
        return Em, E2m, Evar
    
    def Sh_entr(self, beta):#Calcula a entropia shannon
        P = np.exp(-beta*self.eigE)/self.Part_fun(beta)
        Shannon_entropy = -sum(P*np.log(P))
        return Shannon_entropy