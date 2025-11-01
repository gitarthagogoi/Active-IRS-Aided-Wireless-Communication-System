function [nu,Lam] = nu_Lam_generate(K,M,N,Rho_k,eta_k,eps_k,h_k,f_k,G,w_k,sigmar2)

nu_k=zeros(K,N);
Lam_k=zeros(K,N,N);

for k=1:K
    temp_f_k=reshape(f_k(k,:),N,1);
    temp_h_k=reshape(h_k(k,:),M,1);
    temp_w_k=reshape(w_k(k,:),M,1);
    
    temp1=sqrt(eta_k(k)*(1+Rho_k(k)))*eps_k(k)'*diag(temp_f_k')*G*temp_w_k;
    
    temp2=zeros(M,M);
    for j=1:K
        temp3=reshape(w_k(j,:),M,1);
        temp3=temp3*temp3';
        temp2=temp2+temp3;
    end
    
    nu_k(k,:)=temp1-abs(eps_k(k))^2*diag(temp_f_k')*G*temp2*temp_h_k;
    
    Lam_k(k,:,:)=abs(eps_k(k))^2*diag(temp_f_k')*diag(temp_f_k)*sigmar2+abs(eps_k(k))^2*diag(temp_f_k')*G*temp2*G'*diag(temp_f_k);
end

nu=zeros(N,1);
Lam=zeros(N,N);

for k=1:K
    nu=nu+reshape(nu_k(k,:),N,1);
    Lam=Lam+reshape(Lam_k(k,:,:),N,N);
end

end

