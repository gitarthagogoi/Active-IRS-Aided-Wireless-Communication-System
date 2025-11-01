function W=cvx_solve_W(M,K,G,Theta,V,A,W,Ps_max,Pr_max,sigmar2)

temp=(G'*Theta)*(Theta'*G);

theta=Theta*ones(size(Theta,1),1);

D=kron(eye(K),temp);
D=0.5*(D+D');

A=0.5*(A+A')+10^(-50)*eye(size(A,1));
%A=A+10^(-50);

cvx_begin quiet
    cvx_precision low
    variable W(M*K,1) complex;
    minimize((W')*A*W-2*real((V')*W))
    subject to
    W'*W<=Ps_max;
    W'*D*W<=Pr_max-theta'*theta*sigmar2;
cvx_end

% [~ ,W] = QCQP2_solver(A,V,eye(size(A,1)),D,Ps_max,Pr_max-theta'*theta*sigmar2);



end