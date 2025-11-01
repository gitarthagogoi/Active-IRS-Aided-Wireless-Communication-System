function [H] = channel_H2(N,M,dis,large_fading,f_c)
% N number of receiver
% M number of transmitter
f_c=f_c*10^9;
Lambda=3*10^8/f_c;

H = zeros(N,M);
for aa=1:N
    for bb=1:M
       H(aa,bb) =( randn()+1j*randn())/sqrt(2);
    end
end

% H = sqrt(1/10)*exp(1j*2*pi*rand())+sqrt(9/10)*H;
H = sqrt(1/2)*exp(1j*2*pi*rand(N,M))+sqrt(1/2)*H;
%H = sqrt(10/11)*exp(1j*2*pi*rand())+sqrt(1/11)*H;

a = Lambda/4/pi/dis;

a = 10^(-4.12)/(dis^2.8);
a= sqrt(a);

H = a*H;
end
