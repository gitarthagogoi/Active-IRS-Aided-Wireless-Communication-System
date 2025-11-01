function W = w_k2W(K,M,w_k)
W=zeros(K*M,1);

for k=1:K
    W(M*(k-1)+1:1:M*k)=reshape(w_k(k,:),M,1);
end

end
