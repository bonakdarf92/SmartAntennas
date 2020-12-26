using Plots
using Convex
using Gurobi
using LinearAlgebra
#using CuArrays

println("Test")
n = 100
σ = 0.1

A = randn(n,n)
x0 = randn(n)
x0[rand(collect(1:10),Int(n*σ))] .= 0 
x = Variable(n)

y = A*x0 + 0.1*randn(n)
pr = minimize(0.5*sumsquares(A*x - y) + 2*norm_1(x))

solve!(pr,Gurobi.Optimizer())

plot(x.value)
plot!(x0)

kelvin = 6841
blue_red = 1.118
lr_b = [0.3449440921233201, 0.31299511102256444, 0.12060422558367614,
0.18144144645243712, 0.11946633025088182, 0.04842304913535485,
0.04467635530343231, 0.16012872207334228, 0.1784220440072845,
0.18005930568683895, 0.19027196151588324, 0.191061711722997,
0.05121485252925914, 0.10626070626400436, 0.5332429027564216,
0.2891554197967443, 0.009749579987027667, 0.010365695207488415]

hr_b = [0.15342246294255762, 0.16269296865913963, 0.15156155312398212,
0.1810052495578579, 0.14052534906052036, 0.11880743370300048,
0.11523861551790471, 0.20553346505702125, 0.21868965225932188,
0.21877059206369076, 0.2059354681419188, 0.19738626151226807,
0.1020488037731007, 0.10626070626400436, 0.17132562307058136,
0.12867729367776415, 0.025679649372667382, 0.03158527671286214]

function linearReg(measurements)
    t = size(measurements)
    l = 1:t[1]
    β = (l' * l) \ (l' * measurements[:,1])
    return β
end


function regression(measurements)

    
end
