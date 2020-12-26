using Plots
using Convex
using DSP,LibSndFile, SampledSignals, PortAudio
using MAT, FileIO
println("Lab 1 - 4.4")

plot(0:0.1:2π, sin.(0:0.1:2π))

A(r,h) = 2π*(r^2 + r*h)
V(r,h) = π * (r^2*h)
h1(r) = 330 / (π*r^2)
eq(r) = A(r,h1(r))

s = 0.5:0.1:10;
radii = eq.(s);
plot(s, radii)

println("Minimum ", minimum(radii), " at ", s[argmin(radii)] , " cm")

println("Lab 1 - 4.5")

glob_warm = matread("globwarm.mat")
year, T_a = glob_warm["year"], glob_warm["Ta"]
plot(year,T_a)

function moving_avg(x, m)
    N = length(x)
    out = zeros(N,1)
    for k in 1:length(x)
        if k < m + 1
            out[k,1] = sum(x[1:k+m]) / (k+m)
        elseif k > length(x) - m
            out[k,1] = sum(x[(k-m):end] ) / (k+m)
        else
            out[k,1] = sum(x[k-m : k+m]) / (2*m+1)
        end
    end
    return out 
end

mTa2 = moving_avg(T_a, 2)
mTa3 = moving_avg(T_a, 3);
mTa5 = moving_avg(T_a, 5);
mTa7 = moving_avg(T_a, 7);
mTa9 = moving_avg(T_a, 9);
plot(year,T_a, label="Orig")
plot!(year, mTa2, label="m=2")
plot!(year, mTa3, label="m=3")
plot!(year, mTa5, label="m=5")
plot!(year, mTa7, label="m=7")
plot!(year, mTa9, label="m=9")

## You see that the curve gets smoother by increasing mat
println("Lab 1 - 4.6")

n = 0:100
F = 1 # Hz
T = 0.05 # Sample Time
s = sin.(2π*F*n*T)
plot(n*T,s,label="Sin")
plot(n*T,s,label="Sin Stem", line=:stem)

using SampledSignals
using FileIO, LibSndFile
stream = PortAudioStream()
buf = read(stream.source,Int(10*stream.samplerate))
save(joinpath(homedir(), "Desktop", "myvoice.wav"),buf)
plot(buf)
sp = spectrogram(buf.data[:,1],fs=Int(buf.samplerate))
plot(sp.freq,DSP.pow2db.(sp.power))
using ControlSystems
h = 0.1
Aa = [1 h; 0 1]
B = [0; 1]
C = [1 0]
sys = ss(Aa,B,C,0, h)
using LinearAlgebra
Q, R = I, I 
L = lqr(sys,Q,R)
u(x,t) = -L * x .+ 1.5(t[t .>= 2.5])
t = 0:h:5
x0 = [1,0]
y,t,x,uout = lsim(sys,u,t,x0=x0)

