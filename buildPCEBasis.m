function [Xi_samples, w_samples, alpha, Psi_eval, PsiSqNorm] = buildPCEBasis(Xi_dim, p_order)
%% Kombinierte Hermite-PCE-Basis + Collocation-Punkte für nicht-intrusive PCE
% INPUT:
%   M        - # KLE-Terme = Dimension des Zufallsraums
%   p_order  - PCE-Ordnung
% OUTPUT:
%   Xi_samples  - Collocation-Punkte (P x M)
%   w_samples   - Gewichte (P x 1)
%   alpha       - Multiindex (P x M)
%   Psi_eval    - Auswertung aller Psi_j(Xi_i) (P x P)
%   PsiSqNorm   - Normen von Psi_j

    % Erzeuge PCE-Basis
    [alpha, ~, Psi_p, PsiSqNorm, P] = Hermite_PC(Xi_dim, p_order);

    % Collocation-Punkte (Gauss-Hermite Quadratur)
    [Xi_samples, w_samples] = generateCollocationPoints(Xi_dim, p_order);

    % Evaluiere alle Psi_j an allen Xi_i
    Psi_eval = zeros(size(Xi_samples,1), P);
    for j = 1:P
        for i = 1:size(Xi_samples,1)
            Psi_eval(i,j) = evaluateMultivariateHermite(Psi_p{j}, Xi_samples(i,:));
        end
    end
end
