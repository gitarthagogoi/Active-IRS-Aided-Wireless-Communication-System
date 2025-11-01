clear workspace;
clear all;
clc;
tic;
user_X = [100:50:600];
ExNumber = 2; % No of Exp

M = 4;  % Number of Base Station Antennas
Ps_max = 10;        % Max power for base station
Pr_max = 10;        % Max power for active RIS

sigma2 = 10^(-10); % Noise power
sigmar2 = sigma2;

f_c = 5;      % Carrier Frequency
K = 4;            % No of users
N = 512;          % No of RIS Elements
eta_k = ones(K, 1); % Weight factor

large_fading_AI = 2.2;
large_fading_DI = 2.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rsum = zeros(length(user_X), ExNumber);
Rsum_noRIS = zeros(length(user_X), ExNumber);
Rsum_passive = zeros(length(user_X), ExNumber);
Rsum_random = zeros(length(user_X), ExNumber);

for a = 1:length(user_X)
    fprintf('Processing Distance: %d m\n', user_X(a));
    parfor b = 1:ExNumber
        % Generate positions and channels
        [Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser] = Position_generate(K, user_X(a));
        [h_k, f_k, G] = Channel_generate2(K, N, M, large_fading_AI, large_fading_DI, Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser, f_c);
        
        % Ideal Phase
        theta_ideal = 2 * pi * rand(N, 1); % Random phases
        Theta = diag(exp(1j * theta_ideal)); % Ideal phase shift matrix

        % Practical Phase
        Q = 4; % Number of discrete phase levels
        quantized_phases = round(theta_ideal / (2 * pi) * 2^Q) * (2 * pi / 2^Q); % Quantized phases
        Theta_practical = diag(exp(1j * quantized_phases)); % Quantized theta

        Q1 = 3; % Number of discrete phase levels
        quantized_phases3 = round(theta_ideal / (2 * pi) * 2^Q) * (2 * pi / 2^Q); % Quantized phases
        Theta_practical3 = diag(exp(1j * quantized_phases3)); % Quantized theta

        Q2 = 2; % Number of discrete phase levels
        quantized_phases2 = round(theta_ideal / (2 * pi) * 2^Q) * (2 * pi / 2^Q); % Quantized phases
        Theta_practical2 = diag(exp(1j * quantized_phases2)); % Quantized theta

        % Passive IRS
        W = exp(1j * 2 * pi * rand(K * M, 1)) * sqrt(Ps_max / (K * M));
        [W, Theta, Rsum_passive(a, b)] = passive_RIS_precoding(M, K, N, Ps_max, sigma2, eta_k, Theta, W, h_k, f_k, G);
        
        % Active IRS
        [W, Theta, Rsum(a, b)] = active_RIS_precoding(M, K, N, Ps_max * 0.99, Pr_max * 0.01, sigma2, sigmar2, eta_k, Theta, W, h_k, f_k, G);
        [W, Theta_practical, Rsum_practical(a, b)] = active_RIS_precoding(M, K, N, Ps_max * 0.99, Pr_max * 0.01, sigma2, sigmar2, eta_k, Theta_practical, W, h_k, f_k, G);
        [W, Theta_practical3, Rsum_practical3(a, b)] = active_RIS_precoding(M, K, N, Ps_max * 0.99, Pr_max * 0.01, sigma2, sigmar2, eta_k, Theta_practical3, W, h_k, f_k, G);
        [W, Theta_practical2, Rsum_practical2(a, b)] = active_RIS_precoding(M, K, N, Ps_max * 0.99, Pr_max * 0.01, sigma2, sigmar2, eta_k, Theta_practical2, W, h_k, f_k, G);
    end
end

Rsum_mean = mean(Rsum, 2);
Rsum_practical_mean = mean(Rsum_practical, 2);
Rsum_practical_mean3 = mean(Rsum_practical3, 2);
Rsum_practical_mean2 = mean(Rsum_practical2, 2);
Rsum_passive_mean = mean(Rsum_passive, 2);

% Plot results
figure;
plot(user_X, Rsum_mean, '-r^', 'LineWidth', 1.5, 'DisplayName', 'Active IRS (Ideal)');
hold on;
plot(user_X, Rsum_practical_mean, '-g^', 'LineWidth', 1.5, 'DisplayName', 'Active IRS (Practical, 4)');
plot(user_X, Rsum_practical_mean3, '-m^', 'LineWidth', 1.5, 'DisplayName', 'Active IRS (Practical, 3)');
plot(user_X, Rsum_practical_mean2, '--k', 'LineWidth', 1.5, 'DisplayName', 'Active IRS (Practical, 2)');
%plot(user_X, Rsum_passive_mean, '-bo', 'LineWidth', 1.5, 'DisplayName', 'Passive IRS');
hold off;

% Add labels, legend, and grid
xlabel('Distance $L$ (m)', 'Interpreter', 'latex');
ylabel('Sum-rate (bps/Hz)', 'Interpreter', 'latex');
legend('Location', 'best');
grid on;
title('Sum-rate vs. Distance (Weak LOS)');

% Magnify y-axis by setting custom limits
y_min = min([Rsum_mean; Rsum_practical_mean; Rsum_practical_mean3; Rsum_practical_mean2; Rsum_passive_mean]);
y_max = max([Rsum_mean; Rsum_practical_mean; Rsum_practical_mean3; Rsum_practical_mean2; Rsum_passive_mean]);

% Apply a margin for better visualization
margin = (y_max - y_min) * 0.4; % 20% margin
ylim([y_min - margin, y_max + margin]); % Adjust y-axis limits

toc;