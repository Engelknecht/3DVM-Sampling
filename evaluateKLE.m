function [E_nodes, Y_nodes, a_nodes] = evaluateKLE( ...
        Phi, lambda_sorted, ...        % KLE-Eigenbasis
        mu_E, mu_Y, mu_a, ...          % Mittelwerte
        scale_E, scale_Y, scale_a, ... % rel. Standard­abweichungen
        rho_EY, rho_Ea, rho_Ya)        % Korrelationen

% --- 1. Zufallsvektor ξ ~ N(0,Σ_kle) ziehen ------------------------------
M = numel(lambda_sorted);
Sigma_kle = [ eye(M),        rho_EY*eye(M), rho_Ea*eye(M); ...
              rho_EY*eye(M), eye(M),        rho_Ya*eye(M); ...
              rho_Ea*eye(M), rho_Ya*eye(M), eye(M)       ];

xi_all = chol(Sigma_kle,'lower') * randn(3*M,1);
xi_E   = xi_all(1:M);           % Teilvektoren
xi_Y   = xi_all(M+1:2*M);
xi_a   = xi_all(2*M+1:end);

% --- 2. Varianz-skalieren ----------------------------------------------
sigma_E = scale_E * mu_E;
sigma_Y = scale_Y * mu_Y;
sigma_a = scale_a * mu_a;

fac     = 1 / sqrt(sum(lambda_sorted));   % = 1/√λ_sum
fac_E   = sigma_E * fac;
fac_Y   = sigma_Y * fac;
fac_a   = sigma_a * fac;

% --- 3. Knotenfelder aufbauen ------------------------------------------
sqrtlam = sqrt(lambda_sorted);
E_nodes = mu_E + fac_E * (Phi * (sqrtlam .* xi_E));
Y_nodes = mu_Y + fac_Y * (Phi * (sqrtlam .* xi_Y));
a_nodes = mu_a + fac_a * (Phi * (sqrtlam .* xi_a));
end
