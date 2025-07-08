function [xi_all] = sampleXi(Phi,inputArg2)
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
end

