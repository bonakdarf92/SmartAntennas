using Plots, Convex, LinearAlgebra
#using Gurobi
using Random
using Peaks
using Polynomials
Na = 10;
T = 500;
Pa = 3;
θ = [-23, 10, 42];
Power = [30, 30, 30]

A = exp.((im * π) * (0:Na-1) .* sin.(π .* θ'./180))
#S = 1/√2 .* (randn(Pa, T) .+ im *randn(Pa, T));
#ρ = 0;
#Rss = diagm(Power[:])

#N = 1/√2 .* 0.01*(randn(Na, T) .+ im*randn(Na, T));

# Y = A*sqrt(Rss)*S + N;
#Y = A*S + N;
#R = 1/size(Y)[2] * Y * conj(Y)';
#UU,l,_ = svd(R);

#doa = range(-90,stop=90,length=1800)
#GG = UU[:,Pa+1 : Na]' *  exp.( (im * π) * (0:Na-1) .* sin.(π .* doa' ./ 180) )
#f = real(1 ./ sum(conj(GG) .* GG, dims=1));
#f ./= maximum(real(f))


function generate_signal(A, antennas, sources, snapshots, noise, rng)
    S = 1/√2 .* (randn(rng[Threads.threadid()],sources, snapshots,) .+ im *randn(rng[Threads.threadid()] ,sources, snapshots));
    N = 1/√2 .* noise * (randn(rng[Threads.threadid()] ,antennas, snapshots) .+ im*randn(rng[Threads.threadid()] ,antennas, snapshots));
    Y = A*S + N;
    return S, N, Y
end


function music(Y, antennas, sources, doa_grid)
    R = 1/size(Y)[2] * Y * Y';
    U,Λ,_ = svd(R);
    G = U[:,sources+1 : antennas]' *  exp.((im * π) * (0:antennas-1) .* sin.(π .* doa_grid' ./ 180) )
    fθ = real(1 ./ sum(conj(G) .* G, dims=1));
    fθ ./= maximum(real(fθ))
    # fθ .*= maximum(abs.(R))
    ind, vals = peakprom(Maxima(),20*log10.(fθ[:]);minprom=10);
    θs = sort(doa_grid[ind])
    if length(θs) < sources
        #println("Anomalie mit nur ", θs, " gefunden" )
        difference = sources - length(θs);
        for dDiff = 1:difference
            pushfirst!(θs,90)
        end
        #println("Neuer θ ", θs)
    end
    return θs
end


function rootMusic(Y, antennas, sources)
    R̂ = 1/size(Y)[2] * Y * Y';
    U,Λ,_ = svd(R̂);
    Ĝ = U[:,sources+1 : antennas] * U[:,sources+1 : antennas]';

    sqRT = zeros(Complex,2*antennas-1)
    for k in -(antennas-1):(antennas-1)
        sqRT[k+antennas] = sum(diag(Ĝ, k));
    end
    rP = Polynomials.Polynomial(sqRT);
    r_sol = roots(rP);
    roots_ = r_sol[sortperm(r_sol,by=abs, rev=true)];
    roots_2 = roots_[Int(length(r_sol)/2)+1 : end];
    fθ = asin.(angle.(roots_2[1:sources])./π) * 180/π;
    return sort(fθ)
end

function crb(A, T, θ, nu)
    D = (pi / 180) .* (im * π) * (0:Na-1) .* cos.(π .* θ'./180) .* A;#exp.((im * π) * (0:Na-1) .* sin.(π .* θ'./180))
    P = I(size(A)[2])
    R = A * P * A' + nu * I(size(A)[1]);

    M = P * A' / R * A * P
    ΠA_ = I(10) - A * pinv(A);
    crb = real((D' * ΠA_ * D) .* transpose(M));
    crb = inv(crb) * nu / 2 / T 
    return sum(diag(crb)), crb
end

noise_levels = 20
rng = MersenneTwister.(1:Threads.nthreads());
power_noise = 10 .^ range(-2.5, stop=1, length=noise_levels)
rmse = zeros(3,noise_levels);
counter = 1;
MonteCarlo = 1000;
for s in power_noise
    println("Noise level: ", s)
    current_rmse_rootmusic = 0;
    current_rmse_music = 0;
    Threads.@threads for t = 1:MonteCarlo
        S,N,Y = generate_signal(A, Na, Pa, T, s, rng);
        θ_hat_music = music(Y, Na, Pa, doa);
        θ_hat_rootmusic = rootMusic(Y, Na, Pa);
        current_rmse_music += sum((θ - θ_hat_music).^2);
        current_rmse_rootmusic += sum((θ - θ_hat_rootmusic).^2);
    end
    rmse[1,counter] = sqrt(current_rmse_music / MonteCarlo / 3)
    rmse[2,counter] = sqrt(current_rmse_rootmusic / MonteCarlo / 3)
    rmse[3,counter], _ = crb(A, T, θ, s);
    counter += 1
end


function plot_snr(power_noise, rmse)
    plot(10*log10.(1 ./ power_noise), 10*log10.(rmse[1,:]))
    plot!(10*log10.(1 ./ power_noise), 10*log10.(rmse[2,:]))
    plot!(10*log10.(1 ./ power_noise), 10*log10.(rmse[3,:]))
end

function maxLik(Y, grid, antenna, sources)
    Θ_all = grid
    result = zeros(length(grid))
    R̂ = 1/size(Y)[2] * Y * Y';
    counter = 1;
    for k in Θ_all
    a = exp.((im * π) * (0:antenna-1)' .* sin.(π .* k ./ 180) )'
    ΠA_ = I(size(Y)[1]) - a' * inv((conj(a')' * a')) * (conj(a')')
    val[counter] = tr(ΠA_* R̂)
    println(val[k])
    counter += 1;
    end
    return val
end


function PR_DML(Y,antennas, sources, grid)
    U, Λ = eig(1/size(Y)[2] * Y * Y');
    #Λ = diag(real(Λ));
    #[lambda_hat, order] = sort(lambda_hat, 'descend');
    #U = U(:, order);
    absZ = abs.(sqrt(Λ) .* (U'*repeat(grid, antennas)));
    bf_spec = sum(absZ.^2, 1);
end

function omp(Y,sources, doa_grid, MaxIter=10)
    Ω = zeros(sources);
    antennas = size(Y)[1]
    A = zeros(Complex,antennas,1);
    Q = Y;
    a = exp.((im * π) * (0:antennas-1)' .* sin.(π .* doa_grid ./ 180) )'
    for k in 1:sources
        p = findmax(diag(abs.(conj((conj(Q') * a)') * (conj(Q')*a))))[2]
        Ω[k] = doa_grid[p];
        if k == 1
            A[:] = exp.((im * π) * (0:antennas-1) .* sin.(π .* doa_grid[p] ./ 180) );
        else
            A = [A exp.((im * π) * (0:antennas-1) .* sin.(π .* doa_grid[p] ./ 180) )];
        end
        S = pinv(A)*Y
        Q = Y - A*S;#(I(antennas) - A * pinv(A)) * Y
    end
    return -Ω
end

λ(q,X) = eigvals(X)[q]

function PR_OMP(Y, sources, doa_grid, Maxiter=10)
    R = Y;
    Ω = zeros(sources);
    antennas = size(Y)[1];
    A = zeros(Complex, antennas, 1);
    counter = 1;
    for k in 1:antennas
        if k == 1
            A[:] = exp.((im * π) * (0:antennas-1) .* sin.(π .* doa_grid[p] ./ 180) );
        else
            A = [A exp.((im * π) * (0:antennas-1) .* sin.(π .* doa_grid[p] ./ 180) )];
        end
        ΠA_ =  (I(size(Y)[1]) - A*pinv(A)) * (X*conj(X'));
        p = Complex(0);
        for q = sources:antennas
            p += λ(k, ΠA_);
        end
    end
    return Ω

end



# test, l = music(Y, Na, P, doa, Rss)
# plot(doa,20*log10.(test'))
# plot!(θ',seriestype=:vline)

# Aθ = maxLik(Y,3)
