function w_k = w_k_generate(K,M,W)
w_k=zeros(K,M);

for k=1:K
    w_k(k,:)=reshape(W(M*(k-1)+1:1:M*k),M,1);
end

end

