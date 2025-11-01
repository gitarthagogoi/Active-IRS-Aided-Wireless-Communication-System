function Rho_k = Rho_k_update(M,K,H_k,eps_k,eta_k,w_k)
    Rho_k=zeros(K,1);
    for k=1:K
        temp1=reshape(H_k(k,:),M,1);
        temp2=reshape(w_k(k,:),M,1);
        temp=real(eps_k(k)'*temp1'*temp2)/sqrt(eta_k(k));
        Rho_k(k)=(temp^2+temp*sqrt(temp^2+4))/2;
    end
    
end

