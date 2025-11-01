function [Val ,x] = QCQP2_solver(H,f,Q,D,Pmax,Dmax)

H_max = max([max(max(abs(H))) max(max(abs(f)))]);
Q_max = max(max(abs(Q)));
D_max = max(max(abs(D)));

scale = 1;

H = H/H_max*scale;
f = f/H_max*scale;
Q = Q/Q_max*scale;
D=D/D_max*scale;

Pmax = Pmax/Q_max*scale;
Dmax = Dmax/D_max*scale;

N = size(H,2);      %维度

lambda = 1000000;        %拉格朗日系数
lambda2 = 1000000;        %拉格朗日系数

step_length = 0.001;     %迭代步长

eta = 100000;

Iteration = N*1000;     %迭代次数

x = H^(-1)*f;          %随机初始化

%x = randn(N,1)+1j*randn(N,1);

% Min_eva = Inf;

%f_best = Inf;

for p = 1: Iteration
    
    x_gradient = H.'*conj(x)-conj(f) + (lambda+eta*(x'*Q*x - Pmax)*((x'*Q*x)>Pmax))*(Q.'*conj(x)*((x'*Q*x)>Pmax)) ...
                       + (lambda2+eta*(x'*D*x - Dmax)*((x'*D*x)>Dmax))*(D.'*conj(x)*((x'*D*x)>Dmax))   ;
    
    x = x - step_length * conj(x_gradient)/max(abs(x_gradient));

    lambda_gradient = (x'*Q*x - Pmax); 
    lambda_gradient = lambda_gradient*((x'*Q*x)>Pmax);

    lambda2_gradient = (x'*D*x - Dmax);
    lambda2_gradient = lambda2_gradient*((x'*D*x)>Dmax);

    lambda = lambda + step_length * lambda_gradient;
    lambda2 = lambda2 + step_length * lambda2_gradient;


    Eva(p) = x'*H*x - 2*real(f'*x) + lambda*(x'*Q*x - Pmax)*((x'*Q*x)>Pmax)+lambda2*(x'*D*x - Dmax)*((x'*D*x)>Dmax) ...
              + eta/2*norm(x'*Q*x - Pmax)^2 + eta/2*norm(x'*D*x - Dmax)^2 ;
            
    Eva_without(p) = x'*H*x - 2*real(f'*x);


    step_length = 100/(100+p);

%    if (Eva_without(p)<Min_eva) && (x'*Q*x <= Pmax) && (x'*D*x <= Dmax)
%        Min_x = x;
%    end
end

figure;
hold on;
box on;
grid on;
plot(real(Eva));
plot(real(Eva_without));

% x = Min_x;

H = H*H_max/scale;
f = f*H_max/scale;
Q = Q*Q_max/scale;
Pmax = Pmax*Q_max/scale;

D=D*D_max/scale;
Dmax = Dmax*D_max/scale;

Val = real(x'*H*x - 2*real(f'*x));

end