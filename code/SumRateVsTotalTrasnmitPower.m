clear workspace;
clear all;
clc;
tic;

user_X=300;
ExNumber=60;

M=4;  % Base station antennas

Ps_max_dB=[-10:5:40];
Ps_max=10.^(Ps_max_dB/10);
Pr_max=Ps_max;

BW=25.12*1000000;    % Bandwidth (25MHz)
sigma2=10^(-10);
sigmar2=sigma2;

f_c=5;  % Carrier frequency (GHz)
K=4;    % Number of users
N=512;  % Number of RIS elements
eta_k=ones(K,1); % Weight factors

large_fading_AI=2.2;
large_fading_DI=2.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rsum=zeros(length(Ps_max_dB),ExNumber);
Rsum_noRIS=zeros(length(Ps_max_dB),ExNumber);
Rsum_random=zeros(length(Ps_max_dB),ExNumber);
Rsum_passive=zeros(length(Ps_max_dB),ExNumber);
Rsum_practical=zeros(length(Ps_max_dB),ExNumber);

 % Number of quantization bits for practical phase model

for a=1:length(Ps_max_dB)
    fprintf('Iteration %d\n',a);
    parfor b=1:ExNumber
        [Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser]=Position_generate_2(K,user_X);
        [h_k,f_k,G] = Channel_generate2(K,N,M,large_fading_AI,large_fading_DI,Dis_BStoRIS,Dis_BStoUser,Dis_RIStoUser,f_c);
        
        Theta=diag(exp(1j*2*pi*rand(N,1))); % Ideal Phase
        W=exp(1j*2*pi*rand(K*M,1))*sqrt(Ps_max(a)/K/M);
        % Apply practical phase model
        Q=4;
        quantized_phases = round(angle(diag(Theta))/(2*pi)*2^Q)*(2*pi/2^Q);   % Quantize the phase angles
        Theta_practical = diag(exp(1j*quantized_phases));  % Construct the practical Theta with quantized phases
        
        disp(size(Theta)); % Should be [N, N]
        disp(size(quantized_phases));
        disp(size(Theta_practical)); % Should be [N, N]
        [W, Rsum_noRIS(a,b)]= NoRIS_precoding(M,K,N,Ps_max(a),sigma2,eta_k,W,h_k,f_k,G);
        [W, ~, Rsum_random(a,b)]= random_RIS_precoding(M,K,N,Ps_max(a),sigma2,eta_k,Theta,W,h_k,f_k,G);
        [W, Theta, Rsum_passive(a,b)]= passive_RIS_precoding(M,K,N,Ps_max(a),sigma2,eta_k,Theta,W,h_k,f_k,G);
        
        
        [W,Theta,Rsum(a,b)]= active_RIS_precoding(M,K,N,Ps_max(a)*0.99,Pr_max(a)*0.01,sigma2,sigmar2,eta_k,Theta,W,h_k,f_k,G);
        [W,Theta_practical,Rsum_practical(a,b)]= active_RIS_precoding(M,K,N,Ps_max(a)*0.99,Pr_max(a)*0.01,sigma2,sigmar2,eta_k,Theta_practical,W,h_k,f_k,G);
    end
end

Rsum_mean=mean(Rsum,2);
Rsum_noRIS_mean=mean(Rsum_noRIS,2);
Rsum_passive_mean=mean(Rsum_passive,2);
Rsum_random_mean=mean(Rsum_random,2);
Rsum_practical_mean=mean(Rsum_practical,2);

figure;
hold on;
box on;
grid on;
plot(Ps_max_dB,Rsum_mean,'-r^','LineWidth',1.5);
plot(Ps_max_dB,Rsum_passive_mean,'-bo','LineWidth',1.5);
plot(Ps_max_dB,Rsum_random_mean,'-m^','LineWidth',1.5);
plot(Ps_max_dB,Rsum_noRIS_mean,'--k','LineWidth',1.5);
plot(Ps_max_dB,Rsum_practical_mean,'-gs','LineWidth',1.5);

xlabel('Total transmit power $P^{\rm max}$ (dBm)','Interpreter','latex');
ylabel('Sum-rate (bps/Hz)','Interpreter','latex');
title('SumRate Vs Total Transmit Power(Strong LOS)');
set(gca,'FontName','Times','FontSize',12);
legend('Active RIS (Ideal Phase)','Passive RIS','Random phase shift','Without RIS','Active RIS (Practical Phase)','Interpreter','latex','FontSize',12);

save('main_2_practical.mat','Rsum_mean','Rsum_passive_mean','Rsum_random_mean','Rsum_noRIS_mean','Rsum_practical_mean');

toc