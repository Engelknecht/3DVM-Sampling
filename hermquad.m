function [x, w] = hermquad(n)
%% Compute 1D Gauss-Hermite quadrature nodes and weights
%   x - nodes
%   w - weights
    i = (1:n-1)';
    a = zeros(n,1);
    b = sqrt(i/2);
    T = diag(a) + diag(b,1) + diag(b,-1);
    [V,D] = eig(T);
    x = diag(D);
    w = sqrt(pi) * V(1,:).^2';
end