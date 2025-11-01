function [Rsum,gamma_k] = SINR_calculate(K,M,N,eta_k,h_k,f_k,G,Theta,w_k,sigma2,sigmar2)

gamma_k=zeros(K,1);

for k=1:K
    h_k_temp=reshape(h_k(k,:),M,1);
    f_k_temp=reshape(f_k(k,:),N,1);
    disp(['K: ', num2str(K), ', M: ', num2str(M), ', N: ', num2str(N)]);
    disp(['Class of N: ', class(N), ', Value: ', num2str(N)]);
    temp1=h_k_temp+G'*Theta*f_k_temp;
    
    temp2=reshape(w_k(k,:),M,1);
    temp3=temp1'*temp2;
    
    temp4=0;
    for j=1:K
        w_j=reshape(w_k(j,:),M,1);
        temp4=temp4+norm(temp1'*w_j)^2;
    end
    temp4=temp4-norm(temp3)^2+norm(f_k_temp'*Theta')^2*sigmar2+sigma2;    
    gamma_k(k)=norm(temp3)^2/temp4;
end
Rsum=0;
for k=1:K
    Rsum=Rsum+eta_k(k)*log2(1+gamma_k(k));
end
end

