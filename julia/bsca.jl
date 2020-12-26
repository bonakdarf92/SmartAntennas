module bsca
using SparseArrays, LinearAlgebra


N = 2000;
K = 4000;
I = 4000;
ρ_real = 5;
ρ = 10; 
Sample = 20;

function desLemma(A,y,μ,K,x0,algset,discout)
    objval = zeros(1,1001)
    cputime = zeros(1,1001)
    x = x0
    A_norm = sum(A.^2,1)
    L = 3*sum(A_norm.^2) + A_norm * abs(y)

    L      = discount * L;
    lambda = 1. / L;
    
    ATx  = A' * x;
    loss = ATx.^2 - y;

    objval[1]  = 0.25*sum(loss.^2) + μ * sum(abs(x));

    return objval
end

function descentLemma(A, y, mu, K, N, x0, algsetup, discount)
    
    # algType       = algsetup.algorithm;
    # MaxIter_inner = algsetup.MaxIter_inner;
    # stepsize      = algsetup.stepsize;
    # MaxIter_outer = algsetup.MaxIter_outer;
    
    MaxIter_outer = 500
    objval  = zeros(1, MaxIter_outer + 1);
    cputime = zeros(1, MaxIter_outer + 1);
    
    x   = x0;
    
    
    Anorm = sum(A.^2, 1);
    L     = 3 * sum(Anorm.^2) + Anorm * abs(y);
    
    # % for n = 1: 1: N
    # %     An = A(:, n);
    # %     normAn = An' * An;
    # %     L = L + 3 * normAn^2 + normAn * abs(y(n));
    # % end
    
    L      = discount * L;
    lambda = 1 / L;
    
    ATx  = A' * x;
    loss = ATx.^2 - y;
    
    # cputime(1) = toc;
    
    objval[1]  = 0.25*sum(loss.^2) + mu * sum(abs(x));
    
    # println(algsetup.legend);
    println(["iteration ",0, ", value " , objval(1)]);
    
    for t = 1: MaxIter_outer
        # tic;
        
        gradient_g = A * (ATx .* loss);
        gradient_h = (x' * x) * x + x;
        p          = lambda * gradient_g - gradient_h;
        vx         = max(p - lambda * mu * ones(K,1), zeros(K,1)) - max(-p - lambda * mu * ones(K,1), zeros(K, 1));
    
        t_star     = root_polynomial3(vx'*vx);
    
        x          = -t_star * vx;
    
        ATx  = A' * x;
        loss = ATx.^2 - y;
        
        # cputime[t+1] = cputime[t] + toc;
    
        objval[t+1] = 0.25 * sum(loss.^2) + mu * sum(abs(x));
    
    # %     disp(['iteration ' num2str(t+1) ', value ' num2str(objval(t+1)) ', and judge ' num2str(objval(t+1)<=objval(t))]);
    # %     if objval(t+1) > objval(t)
    # %         disp(['root value ' num2str(t_star) ', accuracy ' num2str(value_polynomial3(vx'*vx, t_star))]);
    # %     end
    end
    
    # sol.objval  = objval;
    # sol.x       = x;
    # sol.cputime = cputime;
    
    end
    
    
    function value_polynomial3(c4, t)
        return t^3 * c4 + t - 1;
    end


    function root_polynomial3(c4)
        # % 1/4 * c4 * t^4 + 1/2 * t^2 - t
        Sigma1 = 1/2/c4;
        Sigma2 = 1/3/c4;
        Sigma3 = Sigma1^2 + Sigma2^3;
        Sigma3_sqrt = sqrt(Sigma3);
        gamma = nthroot(3,Sigma1 + Sigma3_sqrt) + nthroot(3,Sigma1 - Sigma3_sqrt);  
        return gamma 
    
    # %     % performing the exact line search the stepsize by Matlab
    # %     f       = @(x) c4/4 * x^4 + 1/2 * x^2 - x;
    # %     options = optimoptions('fmincon','Display','off','Algorithm','sqp');
    # %     gamma   = fmincon(f,0,[],[],[],[],0,[],[],options);
    
        # gamma = max(gamma, 0); % for numerical stability
    end


    function nthroot(n::Integer, r::Real)
        r < 0 || n == 0 && throw(DomainError())
        n < 0 && return 1 / nthroot(-n, r)
        r > 0 || return 0
        x = r / n
        prevdx = r
        while true
            y = x ^ (n - 1)
            dx = (r - y * x) / (n * y)
            abs(dx) ≥ abs(prevdx) && return x
            x += dx
            prevdx = dx
        end
    end


    function FUN_QI_parameter( M, N, density )
        
        x_ori = sprandn(N,1,density);
        
        sigma2 = 0e-4;
        
        A = randn(M, N);
        # % A = orth(A')';
        for n = 1:M
            A[n,:] = A[n,:] / norm(A[n,:]);
        end
        
        # % % ------to make the diagonal of (A'*A) equal to 1
        # % d_AtA = sum(A.^2, 1)'; % diagonal of AtA
        # % d_AtA_InvSR = d_AtA.^-0.5; % inverse square root
        # % % A = A * diag(d_AtA_InvSR);
        # % for n = 1: 1: N,
        # %     A(:,n) = A(:,n) * d_AtA_InvSR(n);
        # % end;
        # % % ------to make the diagonal of (A'*A) equal to 1
        
        b = abs.(abs.((A * x_ori).^2 + sqrt(sigma2) * randn(M, 1)));
        μ = 0.05 * max.( abs.( A' * b ) ); 
        A = A';

        return A,b,μ,x_ori
    end


        
        
        # % maximum number of iterations
        MaxIter_pSCA = 50; sol_pSCA = FUN_initialize(Sample, MaxIter_pSCA);
        MaxIter_bSCA = 50; sol_bSCA = FUN_initialize(Sample, MaxIter_bSCA);
        MaxIter_admm = 50; sol_admm = FUN_initialize(Sample, MaxIter_admm);
        
        sol_pSCA1 = FUN_initialize(Sample, MaxIter_pSCA);
        sol_bSCA1 = FUN_initialize(Sample, MaxIter_bSCA);
        sol_admm1 = FUN_initialize(Sample, MaxIter_admm);
        
        for s = 1: Sample
            println("sample " ,s);
            
            # % generate the data
            global D; D = randn(N, I);
            for n = 1: N
                D(n, :) = D(n, :) / norm(D(n, :));
            end
        # %     D = zeros(N,I);    
        # %     for i = 1: 1: I
        # %         D(randi(N), i) = 1;
        # %     end
            # S0 = sprandn(I, K, 0.05); % density   
            
            P0 = sqrt(100/I) * randn(N, ρ_real);
            Q0 = sqrt(100/K) * randn(ρ_real, K);
            X0 = P0 * Q0; 
            sigma = 0.01;
            V = sigma * randn(N, K); 
            
            Y = X0 + D * S0 + V; #% observation
        
            # % own parameters
            c_λ = 2.5 * 10^-1;
            λ = c_λ * norm(Y); #%spectral norm
            c_μ = 2 * 10^-3; 
            μ = c_μ / 10 * norm(D'*(Y), inf);
        
            # % initial point (common for all algorithms)
            initial_P = randn(N, rho);
            initial_Q = randn(rho, K);
            initial_S = zeros(I,K);
        
            # % initial value (same initial point for all algorithms)
            val0 = FUN_objval(Y, initial_P, initial_Q, D, initial_S, lambda, mu);
            
            # %% BSCA algorithm
            sol_bSCA = FUN_bSCA(initial_P, initial_Q, initial_S, MaxIter_bSCA, s, sol_bSCA);
        
            # %% fully parallel SCA algorithm
            sol_pSCA = FUN_pSCA(initial_P, initial_Q, initial_S, MaxIter_pSCA, s, sol_pSCA);
            
            # %% ADMM algorithm
            sol_admm = FUN_admm(initial_P, initial_Q, initial_S, MaxIter_admm, s, sol_admm);
        
            # % new initial point (common for all algorithms)
            initial_P = sqrt(100/I)*randn(N,ρ);
            initial_Q = sqrt(100/K)*randn(ρ,K);
            initial_S = zeros(I,K);
        
            # % initial value (same initial point for all algorithms)
            val0 = FUN_objval(Y, initial_P, initial_Q, D, initial_S, lambda, mu);
            
            # %% BSCA algorithm
            sol_bSCA1 = FUN_bSCA(initial_P, initial_Q, initial_S, MaxIter_bSCA, s, sol_bSCA1);
        
            # %% fully parallel SCA algorithm
            sol_pSCA1 = FUN_pSCA(initial_P, initial_Q, initial_S, MaxIter_pSCA, s, sol_pSCA1);
            
            # %% ADMM algorithm
            sol_admm1 = FUN_admm(initial_P, initial_Q, initial_S, MaxIter_admm, s, sol_admm1);
        end
        
        figure;
        semilogy(mean(sol_pSCA.time, 1), mean(sol_pSCA.val, 1), 'd-',  'LineWidth', 1.5);
        # hold on; box on; grid on;
        semilogy(mean(sol_bSCA.time, 1), mean(sol_bSCA.val, 1), 'o-',  'LineWidth', 1.5);
        semilogy(mean(sol_admm.time, 1), mean(sol_admm.val, 1), 'kx-', 'LineWidth', 1.5);
        legend('The parallel SCA algorithm (state-of-the-art)', 'The BSCA algorithm (proposed)', 'ADMM (state-of-the-art)');
        xlabel('CPU time (seconds)'); 
        ylabel('objective function value');
        
        function sol = FUN_initialize(Sample, MaxIter)
        sol.val    = zeros(Sample, MaxIter + 1);
        sol.time = zeros(Sample, MaxIter + 1);
        end
        
        function val = FUN_objval(Y, P, Q, D, S, lambda, mu)
        val = 0.5 * norm(Y - P * Q - D * S, 'fro') ^ 2 + 0.5 * lambda * (norm(P, 'fro') ^ 2 + norm(Q, 'fro') ^ 2) + mu * norm(S(:), 1);
        end
        
        function sol_bSCA = FUN_bSCA(initial_P, initial_Q, initial_S, MaxIter_bSCA, s, sol_bSCA)
            % block SCA algorithm
            global val0;
            global Y;
            global D;
            global lambda;
            global mu;
            global rho;
            global I;
            global K;
        
            P = initial_P; 
            Q = initial_Q; 
            S = initial_S;    
            sol_bSCA.val(s,1)  = val0; 
            sol_bSCA.time(s,1) = 0;
            d_DtD = diag(diag(D' * D));
            disp(['Block SCA algorithm, iteration ' num2str(0) ', time ' num2str(0) ', value ' num2str(sol_bSCA.val(s,1))]);
        
            for t = 1: 1: MaxIter_bSCA
                tic;
                
                Y_DS = Y - D * S;
                
                P = Y_DS * Q' * (Q * Q' + lambda * eye(rho)) ^ -1;
        
                Q = (P' * P + lambda * eye(rho)) ^ -1 * P' * Y_DS;
                
                G = d_DtD * S - D' * (P * Q - Y_DS); clear Y_DS
                S_new = d_DtD ^ -1 * (max(G - mu * ones(I, K), zeros(I, K)) - max(-G - mu * ones(I, K), zeros(I, K))); clear G
                cS = S_new - S;
                
                # %-------------------- to calculate the stepsize by exact line search----------------
                B = D*cS;
                C = P*Q+D*S-Y;
                
                c = sum(sum(B.^2,1));
                d = sum(sum(B.*C,1))+mu*(norm(S_new(:),1)-norm(S(:),1));
        
                # clear B C
                # % calculating the stepsize by closed-form expression
                gamma = max(0, min(-d / c, 1));
                clear a b c d
        
                # % variable update
                S = S + gamma * cS; %clear cS S_new
                
                sol_bSCA.time(s,t+1) = toc + sol_bSCA.time(s,t);                        
                sol_bSCA.val(s,t+1) = FUN_objval(Y, P, Q, D, S, lambda, mu);
        
                println("Block SCA algorithm, iteration ", t ", time " , sol_bSCA.time(s,t+1)) ', value ' num2str(sol_bSCA.val(s,t+1))]);
            end
        
            X_bSCA = P * Q;
            S_bSCA = S;
        
            clear P Q S gamma
        end
        
        function sol_pSCA = FUN_pSCA(initial_P, initial_Q, initial_S, MaxIter_pSCA, s, sol_pSCA)
            % parallel SCA algorithm
            global val0;
            global Y;
            global D;
            global lambda;
            global mu;
            global rho;
            global I;
            global K;
        
            P = initial_P; 
            Q = initial_Q; 
            S = initial_S;    
            sol_pSCA.val(s,1)  = val0; 
            sol_pSCA.time(s,1) = 0;
            d_DtD = diag(diag( D' * D));
            disp(['Parallel SCA algorithm, iteration ' num2str(0) ', time ' num2str(0) ', value ' num2str(sol_pSCA.val(s,1))]);
        
            for t = 1: 1: MaxIter_pSCA
                tic;
                
                Y_DS = Y - D * S;
                
                P_new = Y_DS * Q' * (Q * Q' + lambda * eye(rho)) ^ -1;
                cP = P_new - P;
                
                Q_new = (P' * P + lambda * eye(rho)) ^ -1 * P' * Y_DS;
                cQ = Q_new - Q;
                
                G = d_DtD * S - D' * (P * Q - Y_DS); clear Y_DS
                S_new = d_DtD ^ -1 * (max(G - mu * ones(I, K), zeros(I, K)) - max(-G - mu*ones(I, K), zeros(I, K))); clear G
                cS = S_new - S;
                
                %-------------------- to calculate the stepsize by exact line search----------------
                A = cP * cQ;
                B = P * cQ + cP * Q + D * cS;
                C = P * Q + D * S - Y;
                
                a = 2 * sum(sum(A.^2, 1));
                b = 3 * sum(sum(A.*B, 1));
                c = sum(sum(B.^2, 1)) + 2 * sum(sum(A.*C, 1)) + lambda * sum(sum(cP.^2, 1)) + lambda * sum(sum(cQ.^2, 1));
                d = sum(sum(B.*C, 1)) + lambda * sum(sum(cP.*P, 1)) + lambda * sum(sum(cQ.*Q, 1)) + mu * (norm(S_new(:), 1) - norm(S(:), 1));
        
                clear A B C
                % calculating the stepsize by closed-form expression
                Sigma1 = (-(b / 3 / a) ^ 3 + b * c / 6 / a^2 - d / 2 / a);
                Sigma2 = c / 3 / a - (b / 3 / a) ^ 2;
                Sigma3 = Sigma1 ^ 2 + Sigma2 ^ 3;
                Sigma3_sqrt = sqrt(Sigma3);
                if Sigma3 >= 0
                    gamma = nthroot(Sigma1 + Sigma3_sqrt, 3)...
                        + nthroot(Sigma1 - Sigma3_sqrt, 3)...
                        - b / 3 / a;
                else
                    C1 = 1; C1(4) = - (Sigma1 + Sigma3_sqrt);
                    C2 = 1; C2(4) = - (Sigma1 - Sigma3_sqrt);
                    R = real(roots(C1) + roots(C2)) - b/3/a * ones(3,1);
                    gamma = min(R(R>0));
                    clear C1 C2 R;
                end
                clear Sigma1 Sigma2 Sigma3 Sigma3_sqrt
                clear a b c d
                gamma = max(0, min(gamma, 1)); 
        
        %         % performing the exact line search the stepsize by Matlab
        %         f=@(x) a/4*x^4+b/3*x^3+c/2*x^2+d*x;
        %         options = optimoptions('fmincon','Display','off','Algorithm','sqp');
        %         gamma = fmincon(f,0,[],[],[],[],0,1,[],options);
                
                % variable update
                P = P + gamma * cP; clear cP P_new
                Q = Q + gamma * cQ; clear cQ Q_new
                S = S + gamma * cS; clear cS S_new
                
                sol_pSCA.time(s,t+1) = toc + sol_pSCA.time(s,t);                        
                sol_pSCA.val(s,t+1) = FUN_objval(Y, P, Q, D, S, lambda, mu);
        
                disp(['Parallel SCA algorithm, iteration ' num2str(t) ', time ' num2str(sol_pSCA.time(s,t+1)) ', value ' num2str(sol_pSCA.val(s,t+1))...
                    ', stepsize ' num2str(gamma)]);
            end
        
            X_pSCA = P * Q;
            S_pSCA = S;
        
            clear P Q S gamma
        end
        
        function sol_admm = FUN_admm(initial_P, initial_Q, initial_S, MaxIter_admm, s, sol_admm)
            % ADMM algorithm
            global val0;
            global Y;
            global D;
            global lambda;
            global mu;
            global rho;
            global I;
            global K;
        
            P  = initial_P; 
            E  = initial_S; 
            Pi = zeros(I,K); 
            c  = 10^4;
            sol_admm.val(s,1)  = val0;
            sol_admm.time(s,1) = 0;
        
            for t = 1:1:MaxIter_admm
                tic;
        
                Q  = (P' * P + lambda * eye(rho)) ^ -1 * P' * (Y - D * E);
                F  = max(E + Pi / c - mu / c * ones(I, K), zeros(I, K)) - max(-E - Pi / c - mu / c * ones(I, K), zeros(I, K));
                P  = (Y - D * E) * Q' * (Q * Q' + lambda * eye(rho)) ^ -1;
                E  = (D' * D + c * eye(I)) ^ -1 * (c * F - D' * (P * Q - Y) - Pi);
                Pi = Pi + c * (E - F);
        
                sol_admm.val(s, t + 1)  = FUN_objval(Y, P, Q, D, E, lambda, mu);
                sol_admm.time(s, t + 1) = toc + sol_admm.time(s, t);
        
                disp(['ADMM algorithm, iteration ' num2str(t) ', time ' num2str(sol_admm.time(s,t+1)) ', value ' num2str(sol_admm.val(s,t+1))]);
            end
        
            X_admm = P * Q;
            S_admm = E;
        
            clear P Q E F Pi
        end

end