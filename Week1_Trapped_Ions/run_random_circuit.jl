#!/usr/bin/env julia

using PastaQ

function PastaQ.gate(::GateName"R"; theta::Real, phi::Real)
    [
        cos(theta/2)    (-im * exp(-im * phi) * sin(theta/2))
        (-im * exp(im * phi) * sin(theta/2))     cos(theta/2)
    ]
end

function PastaQ.gate(::GateName"M"; Theta::Real)
    [
        cos(Theta)    0    0    (-im * sin(Theta))
        0    cos(Theta)    (-im * sin(Theta))    0
        0    (-im * sin(Theta))    cos(Theta)    0
        (-im * sin(Theta))    0    0    cos(Theta)
    ]
end

function run(N, depth)
    # Random circuit.
    gates = Vector{Tuple}[]

    for i in 1:depth
        one_qubit_layer = Tuple[]
        two_qubit_layer = Tuple[]

        for j in 1:N
            gate = ("R", j, (theta=2pi*rand(), phi=2pi*rand()))
            push!(one_qubit_layer, gate)
        end

        # Alternate start qubit for pairs.
        idx_first = i % 2 + 1

        for j in idx_first:2:(N-1)
            gate = ("M", (j, j+1), (Theta=2pi*rand(),))
            push!(two_qubit_layer, gate)
        end

        push!(gates, one_qubit_layer)
        push!(gates, two_qubit_layer)
    end

    psi = runcircuit(N, gates)
end

#Input Parameters
N = 4
depth = 2

psi = run(N, depth) #Store the MPS state

using ITensors
sites = siteinds("S=1/2",N)
states = [isodd(n) ? "Up" : "Dn" for n=1:N]

A = zeros(1, 2^N) # for storing theoretical probability of each qubit bit string
for x in 0:2^N-1
    y = digits(x, base=2, pad=N)    #Convert Decimal Integer to Binary Number

for n in 1:N
    if y[n] == 0
    states[n] = "Dn" # Quantum State "Dn" when corresponding bit is zero
    else
    states[n] = "Up" # Quantum State "Up" when corresponding bit is one
    end
end
    xpsi = productMPS(sites,states) # Quantum state corresponding to Classical Bit String
    a = inner(xpsi, psi) # Probability Amplitude of xpsi in psi
    b = a * a'          # Probability of occurance of xpsi upon measurement
    A[1,x+1] = b
    
end
println(A)
