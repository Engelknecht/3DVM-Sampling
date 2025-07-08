function [Phi, lambda_sorted] = solveKLE(COORD, CovMat, M, rho_EY, rho_Ea, rho_Ya, n_n,x,y)

[Phi, Lambda] = eig(CovMat);
[lambda_sorted, idx] = sort(diag(Lambda), 'descend');
Phi = Phi(:, idx(1:M));
lambda_sorted = lambda_sorted(1:M);
% Beziehung 
Sigma_kle = [eye(M), rho_EY*eye(M), rho_Ea*eye(M);
             rho_EY*eye(M), eye(M), rho_Ya*eye(M);
             rho_Ea*eye(M), rho_Ya*eye(M), eye(M)];

L_kle = chol(Sigma_kle, 'lower');
Z = randn(3*M, 1);
xi_all = L_kle * Z;
xi_E = xi_all(1:M);
xi_Y = xi_all(M+1:2*M);
xi_a = xi_all(2*M+1:end);

% KLE-Felder
E_nodes = zeros(n_n,1);
Y_nodes = zeros(n_n,1);
a_nodes = zeros(n_n,1);

for i = 1:M
  E_nodes = E_nodes + sqrt(lambda_sorted(i)) * Phi(:,i) * xi_E(i);
  Y_nodes = Y_nodes + sqrt(lambda_sorted(i)) * Phi(:,i) * xi_Y(i);
  a_nodes = a_nodes + sqrt(lambda_sorted(i)) * Phi(:,i) * xi_a(i);
end
% Eigenwerte x-Achse
y_vals = unique(round(COORD(2,:), 2));
counts = arrayfun(@(y) sum(abs(COORD(2,:) - y) < 1e-3), y_vals);
[~, max_idx] = max(counts);
line_y = y_vals(max_idx);
line_nodes = abs(COORD(2,:) - line_y) < 1e-3;
x_line = COORD(1, line_nodes);
Phi_line = Phi(line_nodes, 1:M);

% Sortieren nach x für Überischt 
[x_sorted, idx] = sort(x_line);
Phi_sorted = Phi_line(idx, :);
styles = cell(1, M);  % Pre-allokieren
base_styles = {'-','--','-.',':'};
markers = {'o','s','x','^','v','+','*','d','.'};

for i = 1:M
    ls = base_styles{mod(i-1, length(base_styles)) + 1};
    mk = markers{mod(i-1, length(markers)) + 1};
    styles{i} = [ls mk];  % z. B. '--o'
end
figure; hold on;
for i = 1:M
    plot(x_sorted, Phi_sorted(:,i), styles{i}, 'LineWidth', 2, ...
         'DisplayName', sprintf('\\lambda_{%d} = %.3f', i, lambda_sorted(i)));
end
xlabel('x'); ylabel('\phi_i(x)');
legend('show'); title(sprintf('KLE-Eigenfunktionen auf y = %.2f', line_y));
grid on;
%
% Alle Eigenwerte
% Abbildung Eigenwerte KLE Gesamt
figure; hold on;
colors = lines(5); % 5 verschiedene Farbenxx
styles = {'-o','-s','-','--','-.'}; % verschiedene Linienstile

for i = 1:M
    plot(x, Phi(:,i), styles{i}, 'LineWidth', 2, 'DisplayName', ...
        sprintf('\\lambda_{%d} = %.3f', i, lambda_sorted(i)));
end

xlabel('x'); ylabel('\phi_i(x)');
legend('show', 'Location', 'best');
title('Erste 5 Eigenfunktionen der KLE mit zugehörigen Eigenwerten');
grid on;
%
end

