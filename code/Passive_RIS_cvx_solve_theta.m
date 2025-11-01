function theta= Passive_RIS_cvx_solve_theta(N,K,M,theta,nu,w_k,G)

U=zeros(N,N);
for k=1:K
    w_k_temp=reshape(w_k(k,:),M,1);
    U=U+diag(G*w_k_temp)*(diag(G*w_k_temp))';
end

U=0.5*(U+U');

U = U;
nu = nu;

% cvx_begin quiet
%     variable theta(N,1) complex;
%     minimize ((theta')*U*theta-2*real((theta')*nu))
%     subject to
%     for n = 1:N
%         abs(theta(n))<=1;
%     end
% cvx_end

theta=MMAlgorithm(U, -nu, theta, 1000, 0.000001);

end