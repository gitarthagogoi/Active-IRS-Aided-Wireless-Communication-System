function W=cvx_solve_W_for_passiveRIS(M,K,G,Theta,V,A,W,Ps_max)

A=0.5*(A+A');
%A=A+10^(-50);

A = A*10^5;
V = V*10^5;

% cvx_begin quiet
%     cvx_precision low
%     variable W(M*K,1) complex;
%     minimize((W')*A*W-2*real((V')*W))
%     subject to
%     W'*W<=Ps_max;
% cvx_end

W=(A+0*eye(size(A,1)))^(-1)*V;
if W'*W<=Ps_max
    return;
end

lambda_left = 0;
lambda_right = 10;

W_left=(A+lambda_left*eye(size(A,1)))^(-1)*V;
W_right=(A+lambda_right*eye(size(A,1)))^(-1)*V;
while W_right'*W_right>Ps_max
    lambda_right=2*lambda_right;
    W_right=(A+lambda_right*eye(size(A,1)))^(-1)*V;
end

lambda_center = (lambda_left+lambda_right)/2;
W_center=(A+lambda_center*eye(size(A,1)))^(-1)*V;
    
while (lambda_right-lambda_left)/lambda_right>0.0000001
    
    if W_center'*W_center < Ps_max
        lambda_right = lambda_center;
    else
        lambda_left = lambda_center;
    end
    
    lambda_center = (lambda_left+lambda_right)/2;
    W_center=(A+lambda_center*eye(size(A,1)))^(-1)*V;
end

W = W_center;

end