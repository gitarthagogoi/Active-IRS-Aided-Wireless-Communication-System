function [G] = channel_G(N,M,dis,large_fading,f_c)
% N number of receiver
% M number of transmitter
f_c=f_c*10^9;
Lambda=3*10^8/f_c;

G = zeros(N,M);
for aa=1:N
    for bb=1:M
       G(aa,bb) =( randn()+1j*randn())/sqrt(2);
    end
end

% G = sqrt(10/11)*exp(1j*2*pi*rand())+sqrt(1/11)*G;
G = sqrt(1/2)*exp(1j*2*pi*rand(N,M))+sqrt(1/2)*G;

G = exp(1j*2*pi*rand(N,M));

A=Lambda^2/4;

a=sqrt(A/4/pi)/dis;

a = 10^(-3.73)/(dis^large_fading);
a= 2*sqrt(a);   %3dB element gain

G = a*G;
end

