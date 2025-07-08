function [Xi_samples, w_samples] = generateCollocationPoints(Xi_dim, p_order)
%% Generates Gauss-Hermite collocation points and weights (tensor grid)
% INPUT:
%   M         - number of random variables (dimensions)
%   p_order   - polynomial order of PCE
% OUTPUT:
%   Xi_samples - matrix of collocation points (P x M)
%   w_samples  - associated quadrature weights (P x 1)

% Univariate Gauss-Hermite rule
[xi_1D, w_1D] = hermquad(p_order+1);  % e.g., n = p+1 for exactness

% Build tensor product (full tensor grid)
grid = cell(1, Xi_dim);
wgt = cell(1, Xi_dim);
[grid{:}] = ndgrid(xi_1D);
[wgt{:}] = ndgrid(w_1D);

% Reshape to (P × M)
Xi_samples = [];
for i = 1:Xi_dim
    Xi_samples = [Xi_samples, grid{i}(:)];
end

% Tensor-product weights
w_samples = wgt{1}(:);
for i = 2:Xi_dim
    w_samples = w_samples .* wgt{i}(:);
end
end
