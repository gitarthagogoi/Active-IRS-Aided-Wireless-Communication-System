function eps_k=eps_update(K,M,N,Rho_k,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2)
eps_k=zeros(K,1);
for k=1:K
    h_k_temp=reshape(h_k(k,:),M,1);
    f_k_temp=reshape(f_k(k,:),N,1);
    temp1=h_k_temp+G'*Theta*f_k_temp;
    
    temp2=reshape(w_k(k,:),M,1);
    temp3=sqrt(eta_k(k)*(1+Rho_k(k)))*temp1'*temp2;
    
    temp4=0;
    for j=1:K
        w_j=reshape(w_k(j,:),M,1);
        temp4=temp4+norm(temp1'*w_j)^2;
    end
    temp4=temp4+norm(f_k_temp'*Theta')^2*sigmar2+sigma2;    
    eps_k(k)=temp3/temp4;
end


end

