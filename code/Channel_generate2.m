function [ h_k,f_k,G] = Channel_generate2(K,N,M,large_fading_AI,large_fading_DI,Dis_BStoRIS,Dis_BStoUser,Dis_RIStoUser,f_c)
h_k=zeros(K,M);
f_k=zeros(K,N);
G=zeros(N,M);

for k=1:K
	h_k(k,:)=channel_H2(M,1,Dis_BStoUser(k),large_fading_AI,f_c);             
end


for k=1:K
	f_k(k,:)=channel_F(N,1,Dis_RIStoUser(k),large_fading_DI,f_c);
end

G(:,:)=channel_G(N,M,Dis_BStoRIS,large_fading_DI,f_c);

end

