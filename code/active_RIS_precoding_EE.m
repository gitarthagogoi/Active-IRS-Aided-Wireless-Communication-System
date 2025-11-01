% ========== Modified active_RIS_precoding_EE for Energy Efficiency Maximization ==========
function [W, Theta, EE] = active_RIS_precoding_EE(M, K, N, Ps_max, Pr_max, sigma2, sigmar2, eta_k, Theta, W, h_k, f_k, G, BW, P_ris_hw)
    iteration = 30;
    EE_record = zeros(1, 2 * iteration);

    for Q = 1:iteration
        w_k = w_k_generate(K, M, W);
        K = double(K); M = double(M); N = int32(N);

        [~, gamma_k] = SINR_calculate(K, M, N, eta_k, h_k, f_k, G, Theta, w_k, sigma2, sigmar2);
        H_k = H_k_generate(K, M, N, h_k, f_k, G, Theta);

        Rho_k = gamma_k;
        eps_k = eps_update(K, M, N, Rho_k, eta_k, h_k, f_k, G, Theta, w_k, sigma2, sigmar2);
        [V, A] = v_A_k_generate(K, M, N, Rho_k, eta_k, eps_k, h_k, f_k, G, Theta);

        W = w_k2W(K, M, w_k);
        W = cvx_solve_W(M, K, G, Theta, V, A, W, Ps_max, Pr_max, sigmar2);
        w_k = w_k_generate(K, M, W);

        [Rsum_1, gamma_k] = SINR_calculate(K, M, N, eta_k, h_k, f_k, G, Theta, w_k, sigma2, sigmar2);
        Rsum_1 = BW * log2(1 + sum(gamma_k));

        eps_k = eps_update(K, M, N, Rho_k, eta_k, h_k, f_k, G, Theta, w_k, sigma2, sigmar2);
        [nu, Lam] = nu_Lam_generate(K, M, N, Rho_k, eta_k, eps_k, h_k, f_k, G, w_k, sigmar2);

        theta = Theta * ones(N, 1);
        theta = cvx_solve_theta(N, K, M, theta, nu, Lam, w_k, G, Pr_max, sigmar2);
        Theta = diag(theta);

        [Rsum_2, gamma_k] = SINR_calculate(K, M, N, eta_k, h_k, f_k, G, Theta, w_k, sigma2, sigmar2);
        Rsum_2 = BW * log2(1 + sum(gamma_k));

        total_power = Ps_max + P_ris_hw;
        EE_record(2 * Q - 1) = Rsum_1 / total_power;
        EE_record(2 * Q)     = Rsum_2 / total_power;

        if Q > 1 && abs(EE_record(2 * Q) - EE_record(2 * Q - 1)) / EE_record(2 * Q - 1) < 0.01
            break;
        end
    end

    EE = max(EE_record);
end
