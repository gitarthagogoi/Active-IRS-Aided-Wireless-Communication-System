function [W,Rsum,Rsum_no]= active_RIS_precoding_SI(M,K,N,Pr_max,sigma2,sigmar2,Theta,W,h_k,f_k,G,SI_factor)

w_k = w_k_generate(K,M,W);

U=zeros(N,N);
for k=1:K
    w_k_temp=reshape(w_k(k,:),M,1);
    U=U+diag(G*w_k_temp)*(diag(G*w_k_temp))';
end
U=U+sigmar2*eye(N);

theta=Theta*ones(N,1);

H = SI_factor*1/sqrt(2)*(randn(N,N)+1j*randn(N,N));

H_k = zeros(K,N,N);
for k=1:K   
    H_k(k,:,:)=diag(reshape(f_k(K,:),N,1)')*H*(diag(reshape(f_k(K,:),N,1)')^-1);
end

Phi_opt = Interference_cancellation_multiuser(theta,H_k,0.001,10,150);


theta_no = theta;
[Rsum_no,~] = SINR_calculate_SI(K,M,N,H_k,h_k,f_k,G,theta_no,w_k,sigma2,sigmar2,Pr_max);

theta_opt = Phi_opt;
[Rsum,~] = SINR_calculate_SI(K,M,N,H_k,h_k,f_k,G,theta_opt,w_k,sigma2,sigmar2,Pr_max);

end

