clear all; close all; clc;
tic;

user_X = 300;
ExNumber = 2;

M = 4;  % Base station antennas
Ps_max_dB = -10:5:40;
Ps_max = 10.^(Ps_max_dB / 10);
Pr_max = Ps_max;

BW = 25.12e6; % 25.12 MHz
sigma2 = 1e-10;
sigmar2 = sigma2;

f_c = 5; % Carrier frequency in GHz
K = 4;
N = 64;
eta_k = ones(K, 1);

P_element = 17.566e-3; % W
P_ris_hw = N * P_element;

EE = zeros(length(Ps_max_dB), ExNumber);
EE_quantized = zeros(length(Ps_max_dB), ExNumber, 4);
EE_noRIS = zeros(length(Ps_max_dB), ExNumber);

for a = 1:length(Ps_max_dB)
    fprintf('Iteration %d\n', a);
    parfor b = 1:ExNumber
        [Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser] = Position_generate_2(K, user_X);
        [h_k, f_k, G] = Channel_generate2(K, N, M, 2.2, 2.2, Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser, f_c);

        Theta = diag(exp(1j * 2 * pi * rand(N, 1)));
        W = exp(1j * 2 * pi * rand(K * M, 1)) * sqrt(Ps_max(a) / K / M);

        % Active RIS EE Maximization
        [W, Theta, EE(a, b)] = active_RIS_precoding_EE(M, K, N, Ps_max(a)*0.99, Pr_max(a)*0.01, ...
            sigma2, sigmar2, eta_k, Theta, W, h_k, f_k, G, BW, P_ris_hw);

        % Quantized RIS Phases
        for Q = 1:4
            quantized_phases = round(angle(diag(Theta)) / (2 * pi) * 2^Q) * (2 * pi / 2^Q);
            Theta_q = diag(exp(1j * quantized_phases));
            [W, Theta_q, EE_quantized(a, b, Q)] = active_RIS_precoding_EE(M, K, N, Ps_max(a)*0.99, Pr_max(a)*0.01, ...
                sigma2, sigmar2, eta_k, Theta_q, W, h_k, f_k, G, BW, P_ris_hw);
        end

        % No RIS EE (direct transmission only)
        [W_noRIS, Rsum_noRIS] = NoRIS_precoding(M, K, N, Ps_max(a), sigma2, eta_k, W, h_k, f_k, G);
        EE_noRIS(a, b) = Rsum_noRIS / Ps_max(a);  % No hardware power from RIS
    end
end

% Averages
EE_mean = mean(EE, 2);
EE_quantized_mean = squeeze(mean(EE_quantized, 2));
EE_noRIS_mean = mean(EE_noRIS, 2);

% Plot only ideal active RIS and quantized RIS
figure;
hold on; box on; grid on;

plot(Ps_max_dB, EE_mean, '-r^', 'LineWidth', 1.5, 'DisplayName', 'Ideal Active RIS');

colors = ['b', 'g', 'm', 'c'];
markers = ['o', 's', 'd', 'v'];

for Q = 1:4
    plot(Ps_max_dB, EE_quantized_mean(:, Q), ...
        [colors(Q) '-' markers(Q)], ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('%d-bit Quantized RIS', Q));
end

xlabel('Total transmit power $P^{\rm max}$ (dBm)', 'Interpreter', 'latex');
ylabel('Energy Efficiency (bps/Hz/W)', 'Interpreter', 'latex');
title('EE vs Transmit Power (Active RIS)', 'Interpreter', 'latex');

legend('Ideal Active IRS','1-bit Quantized IRS','2-bit Quantized IRS','3-bit Quantized IRS','4-bit Quantized IRS','Interpreter', 'latex', 'FontSize', 12, 'Location', 'best');
set(gca, 'FontName', 'Times', 'FontSize', 12);

save('EE_maximization_results.mat', 'EE_mean', 'EE_quantized_mean', 'EE_noRIS_mean');
toc;
