function plotCollocationPoints(Xi_samples)
%PLOTCOLLOCATIONPOINTS  Visualize collocation samples for each parameter.
%
%   Xi_samples : (P x d) matrix of collocation points
%
%   Creates one figure per parameter showing the sample distribution.

[P, d] = size(Xi_samples);
for j = 1:d
    figure;
    scatter(1:P, Xi_samples(:,j), 25, 'filled');
    xlabel('Sample index');
    ylabel(sprintf('\\xi_{%d}', j));
    title(sprintf('Collocation points for parameter %d', j));
    grid on;
end
end

