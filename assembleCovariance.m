function [CovMat, cov_func] = assembleCovariance(points, lc, n_n)
cov_func = @(p1, p2) exp(-norm(p1 - p2)^2 / lc^2); % Gauß
CovMat = zeros(n_n);
for i = 1:n_n
  for j = i:n_n
    val = cov_func(points(i,:), points(j,:));
    CovMat(i,j) = val;
    CovMat(j,i) = val;
  end
end
end

