using Convex, Plots, Random, LinearAlgebra

N_Ant = 5;
δ = [0;5];
no_δ = length(δ)
1im*π*(0:N_Ant-1)'.* sin.(δ./180*π)
Aθ1(DOA) = exp.(1im*π*(0:N_Ant-1)'.*sin.(DOA[:]./180*π) )
A = Aθ1(δ)

N_snap = 50
SNR = 10

S = 1/√2 * (randn(no_δ,N_snap) + 1im*randn(no_δ,N_snap))
corrFactor = 0
Rss = I(no_δ)
Rss[1,2],Rss[2,1] = corrFactor,corrFactor'

σ = 10^(-SNR/20) / √2 * (randn(N_Ant,N_snap)+1im*randn(N_Ant,N_snap))
A'*sqrt(Rss)
Y = A'*sqrt(Rss) * S + σ

NGrid = 1800
Grid = range(-90,stop=90,length=NGrid)
AGrid = Aθ1(Grid)