function [W,Theta,Rsum]= active_RIS_precoding(M,K,N,Ps_max,Pr_max,sigma2,sigmar2,eta_k,Theta,W,h_k,f_k,G)
iteration=30;


for Q=1:iteration

w_k = w_k_generate(K,M,W);
K = double(K);
M = double(M);
N = int32(N);
disp(['Class of K: ', class(K), ', Value: ', num2str(K)]);
disp(['Class of M: ', class(M), ', Value: ', num2str(M)]);
disp(['Class of N: ', class(N), ', Value: ', num2str(N)]);

[~,gamma_k] = SINR_calculate(K,M,N,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);

H_k = H_k_generate(K,M,N,h_k,f_k,G,Theta);

Rho_k=gamma_k;

eps_k=eps_update(K,M,N,Rho_k,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);

[V,A]=v_A_k_generate(K,M,N,Rho_k,eta_k,eps_k,h_k,f_k,G,Theta);

W = w_k2W(K,M,w_k);
W=cvx_solve_W(M,K,G,Theta,V,A,W,Ps_max,Pr_max,sigmar2);
w_k = w_k_generate(K,M,W);

[Rsum(2*Q-1),gamma_k] = SINR_calculate(K,M,N,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);

eps_k=eps_update(K,M,N,Rho_k,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);

[nu,Lam] = nu_Lam_generate(K,M,N,Rho_k,eta_k,eps_k,h_k,f_k,G,w_k,sigmar2);
theta=Theta*ones(N,1);
theta= cvx_solve_theta(N,K,M,theta,nu,Lam,w_k,G,Pr_max,sigmar2);
Theta=diag(theta);

[Rsum(2*Q),gamma_k] = SINR_calculate(K,M,N,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);

if Q>1
    if (Rsum(2*Q)-Rsum(2*Q-1))/Rsum(2*Q-1)<0.01
        break;
    end
end

end

%save('Phi1.mat','theta');

%plot(Rsum)
Rsum=max(Rsum);
end

