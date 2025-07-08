function coeffs = computePCEcoeffs(QoI_samples, Xi_samples, p_order)
% COMPUTEPCECOEFFS - Berechnet die PCE-Koeffizienten via Least-Squares
%
% QoI_samples : [N x K] Stichproben der Zielgrößen (QoIs)
% Xi_samples  : [N x d] Stichproben der Zufallsgrößen
% p_order     : Ordnung der PCE
% coeffs      : [P x K] Matrix der PCE-Koeffizienten

[N, d] = size(Xi_samples);
K = size(QoI_samples, 2);

multi_idx = multi_index(d, p_order);
P = size(multi_idx, 1);

Psi_eval = zeros(N, P);
for n = 1:N
    Psi_eval(n, :) = evaluate_orthopoly(multi_idx, Xi_samples(n, :));
end

coeffs = (Psi_eval' * Psi_eval) \ (Psi_eval' * QoI_samples);

end
