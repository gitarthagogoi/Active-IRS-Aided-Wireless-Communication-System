function theta= cvx_solve_theta(N,K,M,theta,nu,Lam,w_k,G,Pr_max,sigmar2)

U=zeros(N,N);
for k=1:K
    w_k_temp=reshape(w_k(k,:),M,1);
    U=U+diag(G*w_k_temp)*(diag(G*w_k_temp))';
end
U=U+sigmar2*eye(N);
Lam=0.5*(Lam+Lam');

U=0.5*(U+U');

lambda=0;
% cvx_begin quiet
%     variable theta(N,1) complex;
%     maximize(-(theta')*Lam*theta+2*real((theta')*nu))
%     subject to
%     theta'*U*theta<=Pr_max;
% cvx_end    

theta=inv(Lam+lambda*U)*nu;

lambda_left=0;
lambda_right=10;

theta_left=inv(Lam+lambda_left*U)*nu;
theta_right=inv(Lam+lambda_right*U)*nu;
P_left=theta_left'*U*theta_left;
P_right=theta_right'*U*theta_right;

if P_left<Pr_max
    theta=theta_left;
    return;
end

while P_left>Pr_max && P_right>Pr_max
	lambda_left=lambda_right;
	lambda_right=lambda_right*10;
	theta_left=inv(Lam+lambda_left*U)*nu;
	theta_right=inv(Lam+lambda_right*U)*nu;
	P_left=theta_left'*U*theta_left;
	P_right=theta_right'*U*theta_right;
end  

while abs(lambda_left-lambda_right)>0.0001

lambda_middle=(lambda_left+lambda_right)/2;
theta_middle=inv(Lam+lambda_middle*U)*nu;
P_middle=theta_middle'*U*theta_middle;

if P_middle<Pr_max
    lambda_right=lambda_middle;
    theta_right=inv(Lam+lambda_right*U)*nu;
    P_right=theta_right'*U*theta_right;
else if P_middle>Pr_max
    lambda_left=lambda_middle;
	theta_left=inv(Lam+lambda_left*U)*nu;
    P_left=theta_left'*U*theta_left;
    else 
        theta=theta_middle;
        return;
end
    
end

theta=theta_right;

end