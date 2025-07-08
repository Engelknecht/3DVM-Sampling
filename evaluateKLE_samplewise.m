function [E_nodes,Y_nodes,a_nodes] = evaluateKLE_samplewise(M, Phi,lambda_sorted,mu_E,mu_Y,mu_a,scale_E,scale_Y,scale_a, n_n,xi_E, xi_Y, xi_a)
    % Stochastische Realisierung für eine Stichprobe xi

% KLE-Felder
E_nodes = zeros(n_n,1);
Y_nodes = zeros(n_n,1);
a_nodes = zeros(n_n,1);
for i = 1:M
  E_nodes = E_nodes + sqrt(lambda_sorted(i)) * Phi(:,i) * xi_E(i);
  Y_nodes = Y_nodes + sqrt(lambda_sorted(i)) * Phi(:,i) * xi_Y(i);
  a_nodes = a_nodes + sqrt(lambda_sorted(i)) * Phi(:,i) * xi_a(i);
end

% Streuung kontrollieren:
E_nodes = E_nodes / std(E_nodes);  % Standardisieren
Y_nodes = Y_nodes / std(Y_nodes);
a_nodes = a_nodes / std(a_nodes);

% Skalieren auf gewünschte Streuung (z. B. 5 %) um Mittelwert

E_nodes = mu_E + scale_E * mu_E * E_nodes;
Y_nodes = mu_Y + scale_Y * mu_Y * Y_nodes;
a_nodes = mu_a + scale_a * mu_a * a_nodes;
end

