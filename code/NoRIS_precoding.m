function [W, Rsum]= NoRIS_precoding(M,K,N,Ps_max,sigma2,eta_k,W,h_k,f_k,G)
iteration=30;

Theta=zeros(N,N);
sigmar2=0;

for Q=1:iteration

w_k = w_k_generate(K,M,W);
[~,gamma_k] = SINR_calculate(K,M,N,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);

Rho_k=gamma_k;

eps_k=eps_update(K,M,N,Rho_k,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);

[V,A]=v_A_k_generate(K,M,N,Rho_k,eta_k,eps_k,h_k,f_k,G,Theta);

W = w_k2W(K,M,w_k);

W=cvx_solve_W_for_noRIS(M,K,G,Theta,V,A,W,Ps_max);

w_k = w_k_generate(K,M,W);

%eps_k=eps_update(K,M,N,Rho_k,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);

[Rsum(Q),gamma_k] = SINR_calculate(K,M,N,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2);

if Q>1

    if (Rsum(Q)-Rsum(Q-1))/Rsum(Q-1)<0.005
        break;
    end
end

end

%plot(Rsum);


Rsum=max(Rsum);
end

