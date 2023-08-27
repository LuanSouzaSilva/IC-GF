import Pkg; #Pkg.add("ProgressBars")
using ITensors
using Plots
using ProgressBars

import ITensors: op

print(" Comecei ")

let
#function main(; N=5, cutoff1=1E-8, δt=0.1, ttotal=1.0)
    # Compute the number of steps to do

    #Nsteps = Int(ttotal / δt)
    N = 5
    cutoff1 = 1E-8
    ttotal=1.0
  
    Npart = N
    t = 1.0
    U = 2.0
    mu = U/2
    tau = 0.1

    sites = siteinds("Electron", N; conserve_qns=true)

    os = OpSum()
    for b in 1:(N - 1)
    os += -t, "Cdagup", b, "Cup", b + 1
    os += -t, "Cdagup", b + 1, "Cup", b
    os += -t, "Cdagdn", b, "Cdn", b + 1
    os += -t, "Cdagdn", b + 1, "Cdn", b

    end

    for i in 1:N
    os += U, "Nupdn", i
    end
    for i in 1:N
    os += -mu, "Nup", i
    os += -mu, "Ndn", i
    end

    H = MPO(os, sites)

    nsweeps = 6
    maxdim = [50, 100, 200, 400, 800, 800]
    cutoff = [1E-12]
 
    state = ["Emp" for n in 1:N]
    p = Npart
    for i in N:-1:1
    if p > i
        #println("Doubly occupying site $i")
        state[i] = "UpDn"
        p -= 2
    elseif p > 0
        #println("Singly occupying site $i")
        state[i] = (isodd(i) ? "Up" : "Dn")
        p -= 1
    end
    end
    
    # Initialize wavefunction to be bond 
    # dimension 10 random MPS with number
    # of particles the same as `state`
    psi0 = randomMPS(sites, state, N)

    # Check total number of particles:
    #@show flux(psi0)

    # Start DMRG calculation:
    energy, psi1 = dmrg(H, psi0; nsweeps, maxdim, cutoff)

    println("\nGround State Energy = $energy")

    c = div(N, 2) # center site

    # Make gates (1,2),(2,3),(3,4),...
    gates = ITensor[]
    for j in 1:(N - 1)
        s1 = sites[j]
        s2 = sites[j + 1]
        hj =
        -t*op("Cdagup", s1) * op("Cup", s2) 
        -t*op("Cup", s1) * op("Cdagup", s2) 
        -t*op("Cdagdn", s1) * op("Cdn", s2) 
        -t*op("Cdn", s1) * op("Cdagdn", s2)
        -mu*(op("Nup", s1) + op("Ndn", s1))
        +U*op("Nupdn", s1)

        Gj = exp(-im * (tau / 2) * hj)
        push!(gates, Gj)
    end
    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    append!(gates, reverse(gates))
  
    # Do the time evolution by applying the gates
    # for Nsteps steps
    psi = copy(psi1)
    s = siteind(psi, c)
    newpsi= op(s,"Cdagup")*psi[c]
    noprime!(newpsi)
    psi[c]= newpsi
    println("First value is ", inner(psi, psi1))

    print(typeof(psi1))

    corr = zeros(ComplexF64, 0)
    time = zeros(Float64, 0)
  
    # Now do the same evolution with an MPO
    for t in 0.0:tau:ttotal
      psi = apply(gates, psi; cutoff1)#, apply_dag=true)
      normalize!(psi)
    
      ket = copy(psi)
      s = siteind(ket,c)
      newket= op(s,"Cup")*ket[c]
      noprime!(newket)
      ket[c]= newket

      append!(corr, inner(psi, ket))
      append!(time, t)
      
    end
  
    return time, corr
end

#time, corr = main()

p = plot(time, real(corr))
display(p)


print(" Terminei ")