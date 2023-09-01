import Pkg; #Pkg.add("LaTeXStrings")#Pkg.add("ITensors") #Instala os pacotes
#Pkg.add("ProgressBars")


#Importa os pacotes
using LaTeXStrings
using ITensors
using Plots
using ProgressBars

import ITensors: op

print(" Comecei ")#Comecei

function op(::OpName"expτSS", ::SiteType"S=1/2", s1::Index, s2::Index; tau)
  h =
    1 / 2 * op("S+", s1) * op("S-", s2) +
    1 / 2 * op("S-", s1) * op("S+", s2) +
    op("Sz", s1) * op("Sz", s2)
  return exp(tau * h)
end

function main1(; N=10, cutoff=1E-10, dt=0.1, ttotal=50.0)
    #Numero de passos
    Nsteps = Int(ttotal / dt)
    #Escolhe o sitio como sendo o do meio
    c = div(N, 2)
    t = 0.0
  
    # Make an array of 'site' indices
    s = siteinds("S=1/2", N; conserve_qns=true)
  
    #Faz os gates em ordem, e depois na ordem reversa
    gates = ops([("expτSS", (n, n + 1), (tau=-dt * im / 2,)) for n in 1:(N - 1)], s)
    append!(gates, reverse(gates))
  
    #Inicializa um MPS com spins alternados
    psi0 = productMPS(s, n -> isodd(n) ? "Up" : "Dn")

    psi = copy(psi0) #Faz uma copia de psi0
    s = siteind(psi,c)
    newpsi= op(s,"Sz")*psi[c] #Aplica Sz_5 em psi
    noprime!(newpsi)
    psi[c]= newpsi

    #Inicializa variaveis
    corr = zeros(ComplexF64, 0)
    time = zeros(Float64, 0)
  
    #Evolucao temporal
    for step in tqdm(1:Nsteps)
      psi = apply(gates, psi; cutoff)
    
      ket = copy(psi)
      s = siteind(ket,c)
      newket= op(s,"Sz")*ket[c]
      noprime!(newket)
      ket[c]= newket

      append!(corr, inner(psi, ket))
      append!(time, t)
      
      t += dt
    end
  
    return time, corr
  end

time, corr = main1()


p = plot(time, real(corr), label = L"\langle S_5(t)S_5(0)\rangle", xlabel = "t", ylabel = L"\langle S_5(t)S_5(0)\rangle", lw = 2) #Plota o grafico
savefig(p, "Spin_corr.pdf") 
display(p) #Abre o grafico

print(" Terminei ") #Terminei