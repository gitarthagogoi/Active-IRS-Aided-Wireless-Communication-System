clear workspace;
clear all;
clc;
tic;
user_X=[100:50:600];
ExNumber=10; %No of Exp


M=4;  %»Number of Base Station Antennas
Ps_max=10;        %Mx power for base station
Pr_max=10;        %Mx power for active RIS

sigma2 = 10^(-10); %Noise power

sigmar2=sigma2;

%sigma2=10^(-10);    %ÓÃ»§ÒýÈëÈÈÔëÉù
%sigmar2=10^(-10);   %active RISÒýÈëÈÈÔëÉù

f_c=5;      %Carrier Freq
L = 2; % Number of quantization levels (example value)
K=4;            %No of users
N=512;          %No of RIS Elements
eta_k=ones(K,1);%weight factor

large_fading_AI=2.2;
large_fading_DI=2.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rsum=zeros(length(user_X),ExNumber);
Rsum_noRIS=zeros(length(user_X),ExNumber);
Rsum_passive=zeros(length(user_X),ExNumber);
Rsum_random=zeros(length(user_X),ExNumber);

for a=1:length(user_X)
    fprintf('Processing Distance: %d m\n',user_X(a));
    parfor b=1:ExNumber
        %Generate positions and channels
        [Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser]=Position_generate(K,user_X(a));    %»ùÕ¾RISÓÃ»§Î»ÖÃÉèÖÃ
        [ h_k,f_k,G] = Channel_generate2(K,N,M,large_fading_AI,large_fading_DI,Dis_BStoRIS,Dis_BStoUser,Dis_RIStoUser,f_c);
        
        %Ideal Phase
        theta_ideal=2*pi*rand(N,1);%random phases
        Theta=diag(exp(1j*theta_ideal));%ideal phase shift matrix
        
        
        %Practical Phase
        Q=4; %number of discrete phase lev0els
        %quantized_phases =round(rand(N,1)*(Q-1))*((2*pi)/(Q-1)); %quantized Phases
        %Theta_practical =diag(exp(1j*quantized_phases));%quantized theta
        
        %Q=4; %number of discrete phase levels
        %m=randi([0, Q-1],N,1);
        quantized_phases =round(theta_ideal/(2*pi)*2^Q)*(2*pi/2^Q); %quantized Phases
        Theta_practical =diag(exp(1j*quantized_phases));%quantized theta
        disp(size(Theta));           % Should be [N, N]
        disp(size(Theta_practical)); % Should be [N, N]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate random phase shifts from uniform distribution between [0, 2*pi)
        %theta_ideal = 2*pi*rand(N, 1);

        % Quantize the phases to L discrete levels
        %quantized_phases = round(theta_ideal * (L - 1) / (2*pi)) * (2*pi / (L - 1));

        % Create the phase shift matrix Theta
        %Theta = diag(exp(1j * quantized_phases));  % Diagonal matrix with quantized phases
        %Q = 4; % Number of discrete phase levels
        %quantized_phases = round(rand(N, 1) * (Q - 1)) * (2 * pi / Q); % Quantize phases
        %Theta = diag(exp(1j * quantized_phases)); % Create practical phase shift matrix
        %Theta=diag(exp(1j*2*pi*rand(N,1)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %without IRS
        W=exp(1j*2*pi*rand(K*M,1))*sqrt(Ps_max/(K*M));
        [W, Rsum_noRIS(a,b)]= NoRIS_precoding(M,K,N,Ps_max,sigma2,eta_k,W,h_k,f_k,G);
        
        [W, ~, Rsum_random(a,b)]= random_RIS_precoding(M,K,N,Ps_max,sigma2,eta_k,Theta,W,h_k,f_k,G);
%       Rsum_noRIS(a,b)=0;
        %passive IRS
        [W, Theta, Rsum_passive(a,b)]= passive_RIS_precoding(M,K,N,Ps_max,sigma2,eta_k,Theta,W,h_k,f_k,G);
        Theta=100*Theta;
%       [W,Theta,Rsum(a,b)]= active_RIS_precoding(M,K,N,Ps_max/2,Pr_max/2,sigma2,sigmar2,eta_k,Theta,W,h_k,f_k,G);
        [W,Theta,Rsum(a,b)]= active_RIS_precoding(M,K,N,Ps_max*0.99,Pr_max*0.01,sigma2,sigmar2,eta_k,Theta,W,h_k,f_k,G);
        
        [W,Theta_practical,Rsum_practical(a,b)]= active_RIS_precoding(M,K,N,Ps_max*0.99,Pr_max*0.01,sigma2,sigmar2,eta_k,Theta_practical,W,h_k,f_k,G);
        
        
    end
end

Rsum_mean=mean(Rsum,2);
Rsum_practical_mean=mean(Rsum_practical,2);
Rsum_noRIS_mean=mean(Rsum_noRIS,2);
Rsum_passive_mean=mean(Rsum_passive,2);
Rsum_random_mean=mean(Rsum_random,2);

figure;
hold on;
box on;
grid on;
plot(user_X,Rsum_mean,'-r^','LineWidth',1.5);
plot(user_X,Rsum_practical_mean,'-g^','LineWidth',1.5);
%plot(user_X,Rsum_passive_mean,'-bo','LineWidth',1.5);
plot(user_X,Rsum_random_mean,'-m^','LineWidth',1.5);
plot(user_X,Rsum_noRIS_mean,'--k','LineWidth',1.5);
ylabel('Sum-rate (bps/Hz)','Interpreter','latex');
xlabel('Distance $L$ (m)','Interpreter','latex');
title('SumRate Vs Distance(Weak LOS)');

set(gca,'FontName','Times','FontSize',12);
%ylim([0, max([Rsum_mean;Rsum_practical_mean;Rsum_noRIS_mean;Rsum_passive_mean;Rsum_random_mean])*1.1]);
%ylim([0 80]);
%legend('Active RIS~~\,~($P_{\rm BS}^{\rm max}=5$ W, $P_{\rm A}^{\rm max}=5$ W)','Passive RIS~\,~($P_{\rm BS}^{\rm max}=10$ W)','Random phase shift ($P_{\rm BS}^{\rm max}=10$ W)','Without RIS ($P_{\rm BS}^{\rm max}=10$ W)','Interpreter','latex','FontSize',12);
legend('Active IRS (Ideal)','Active IRS (Practical)','Random phase shift','Without RIS','Interpreter','latex','FontSize',12);

%legend('Active RIS~~\,~($P_{\rm BS}=5$ W, $P_{\rm A}=5$ W)','Passive RIS~\,~($P_{\rm BS}=10$ W)','Without RIS ($P_{\rm BS}=10$ W)','Interpreter','latex','FontSize',12);
save('main_1_2.mat','Rsum_mean','Rsum_passive_mean','Rsum_random_mean','Rsum_noRIS_mean');

toc;
%title('{\bf Typical communication scenario}','Interpreter','latex');
%title('{\bf Unusual communication scenario}','Interpreter','latex');
