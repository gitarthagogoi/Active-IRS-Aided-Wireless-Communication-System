clear workspace;
clear all;
clc;
tic;

% Parameters
user_X = 300;
ExNumber = 60;

M = 4; % Number of base station antennas
N = [100:100:900]; % Number of RIS elements

Ps_max = 10; % Maximum base station power
Pr_max = Ps_max; % Maximum RIS power

sigma2 = 10^(-10); % User noise variance
sigmar2 = sigma2; % RIS noise variance

f_c = 5; % Carrier frequency

K = 4; % Number of users
eta_k = ones(K, 1); % User weights

large_fading_AI = 2.2; % Large fading for AI
large_fading_DI = 2.2; % Large fading for DI

% Initialize result matrices
Rsum = zeros(length(N), ExNumber);
Rsum_noRIS = zeros(length(N), ExNumber);
Rsum_passive = zeros(length(N), ExNumber);
Rsum_random = zeros(length(N), ExNumber);
Rsum_practical = zeros(length(N), ExNumber);

% Main loop
for a = 1:length(N)
    fprintf('Iteration %d: N = %d\n', a, N(a));
    
    parfor b = 1:ExNumber
        % Generate positions and channels
        [Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser] = Position_generate_2(K, user_X);
        [h_k, f_k, G] = Channel_generate(K, N(a), M, large_fading_AI, large_fading_DI, Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser, f_c);

        % Initialize RIS phase shifts and precoding vectors
        Theta = diag(exp(1j * 2 * pi * rand(N(a), 1)));
        W = exp(1j * 2 * pi * rand(K * M, 1)) * sqrt(Ps_max / K / M);

        % Quantize RIS phase shifts (practical implementation)
        Q = 4;
        quantized_phases = round(angle(diag(Theta)) / (2 * pi) * 2^Q) * (2 * pi / 2^Q);
        Theta_practical = diag(exp(1j * quantized_phases));

        % Debug: Print sizes of key variables
        fprintf('Debug: Sizes of key variables for N = %d\n', N(a));
        fprintf('Size of h_k: %s\n', mat2str(size(h_k)));
        fprintf('Size of f_k: %s\n', mat2str(size(f_k)));
        fprintf('Size of G: %s\n', mat2str(size(G)));
        fprintf('Size of Theta: %s\n', mat2str(size(Theta)));
        fprintf('Size of W: %s\n', mat2str(size(W)));

        % No RIS case
        [W, Rsum_noRIS(a, b)] = NoRIS_precoding(M, K, N(a), Ps_max, sigma2, eta_k, W, h_k, f_k, G);

        % Random RIS phase shifts
        [W, ~, Rsum_random(a, b)] = random_RIS_precoding(M, K, N(a), Ps_max, sigma2, eta_k, Theta, W, h_k, f_k, G);

        % Passive RIS
        [W, Theta, Rsum_passive(a, b)] = passive_RIS_precoding(M, K, N(a), Ps_max, sigma2, eta_k, Theta, W, h_k, f_k, G);

        % Active RIS (ideal phase shifts)
        [W, Theta, Rsum(a, b)] = active_RIS_precoding(M, K, N(a), Ps_max * 0.99, Pr_max * 0.01, sigma2, sigmar2, eta_k, Theta, W, h_k, f_k, G);

        % Active RIS (practical phase shifts)
        [W, Theta_practical, Rsum_practical(a, b)] = active_RIS_precoding(M, K, N(a), Ps_max * 0.99, Pr_max * 0.01, sigma2, sigmar2, eta_k, Theta_practical, W, h_k, f_k, G);
    end
end

% Compute mean values
Rsum_mean = mean(Rsum, 2);
Rsum_noRIS_mean = mean(Rsum_noRIS, 2);
Rsum_passive_mean = mean(Rsum_passive, 2);
Rsum_random_mean = mean(Rsum_random, 2);
Rsum_practical_mean = mean(Rsum_practical, 2);

% Plot results
figure;
hold on;
box on;
grid on;
plot(N, Rsum_mean, '-r^', 'LineWidth', 1.5);
plot(N, Rsum_passive_mean, '-bo', 'LineWidth', 1.5);
plot(N, Rsum_random_mean, '-m>', 'LineWidth', 1.5);
plot(N, Rsum_noRIS_mean, '--k', 'LineWidth', 1.5);
plot(N, Rsum_practical_mean, '-gs', 'LineWidth', 1.5); % Practical phase model
xlabel('Number of RIS elements $N$', 'Interpreter', 'latex');
ylabel('Sum-rate (bps/Hz)', 'Interpreter', 'latex');
title('Sum-Rate Vs Number of Elements (Strong LOS)');
set(gca, 'FontName', 'Times', 'FontSize', 12);
legend('Active RIS (Ideal Phase)', 'Passive RIS', 'Random phase shift', 'Without RIS', 'Active RIS (Practical Phase)', 'Interpreter', 'latex', 'FontSize', 12);

% Save results
save('main_3.mat', 'Rsum_mean', 'Rsum_passive_mean', 'Rsum_random_mean', 'Rsum_noRIS_mean', 'Rsum_practical_mean');

toc;