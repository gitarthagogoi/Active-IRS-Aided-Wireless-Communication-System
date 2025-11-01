% Reference:
% Pan, Cunhua, et al. 
% "Multicell MIMO communications relying on intelligent reflecting surfaces."
% IEEE Transactions on Wireless Communications (2020).
%
% Author:
% Kunzan Liu, July 2020
%
% Optimization Problem:
% min_phi phi'*U*phi + 2*real(phi'*v)
% s.t. abs(phi_m) = 1, phi = [phi_1, ..., phi_M]^T
function phi = MMAlgorithm(U, v, phi, iter, accuracy)
% U(M*M) - input matrix
% v(M*1) - input vector
% phi(M*1) - initialization
% accuracy - threshold
% iter - maximum iteration
M = length(v);
lambda = max(eigs(U));
t = 1;
f_phi_new = phi'*U*phi+2*real(phi'*v);
while t <= iter
    f_phi = f_phi_new;
    q = (lambda*eye(M)-U)*phi-v;
    phi = exp(1i*angle(q));
    f_phi_new = phi'*U*phi+2*real(phi'*v);
    if abs(f_phi-f_phi_new)/abs(f_phi_new) <= accuracy
        break;
    end
    t = t+1;
end
end