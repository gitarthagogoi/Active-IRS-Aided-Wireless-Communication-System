function [W, Theta, Rsum]= random_RIS_precoding(M,K,N,Ps_max,sigma2,eta_k,Theta,W,h_k,f_k,G)
iteration=30;


for Q=1:iteration

w_k = w_k_generate(K,M,W);
[~,gamma_k] = SINR_calculate(K,M,N,eta_k,h_k,f_k,G,Theta,w_k,sigma2,0);

H_k = H_k_generate(K,M,N,h_k,f_k,G,Theta);

Rho_k=gamma_k;

eps_k=eps_update(K,M,N,Rho_k,eta_k,h_k,f_k,G,Theta,w_k,sigma2,0);

[V,A]=v_A_k_generate(K,M,N,Rho_k,eta_k,eps_k,h_k,f_k,G,Theta);

W = w_k2W(K,M,w_k);
W = cvx_solve_W_for_passiveRIS(M,K,G,Theta,V,A,W,Ps_max);
w_k = w_k_generate(K,M,W);

[Rsum(Q),gamma_k] = SINR_calculate(K,M,N,eta_k,h_k,f_k,G,Theta,w_k,sigma2,0);

%eps_k=eps_update(K,M,N,Rho_k,eta_k,h_k,f_k,G,Theta,w_k,sigma2,0);

%[nu,~] = nu_Lam_generate(K,M,N,Rho_k,eta_k,eps_k,h_k,f_k,G,w_k,0);
%theta=Theta*ones(N,1);

%theta= Passive_RIS_cvx_solve_theta(N,K,M,theta,nu,w_k,G);

%Theta=diag(theta);

%[Rsum(2*Q),gamma_k] = SINR_calculate(K,M,N,eta_k,h_k,f_k,G,Theta,w_k,sigma2,0);

if Q>1
    if (Rsum(Q)-Rsum(Q-1))/Rsum(Q-1)<0.001
        break;
    end
end

end

%plot(Rsum)
Rsum=max(Rsum);
end

