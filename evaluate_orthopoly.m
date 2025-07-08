function psi_val = evaluate_orthopoly(multi_idx, xi)
% EVALUATE_ORTHOPOLY - Bewertet multivariate orthonormale Hermite-Polynome
%   multi_idx : [P x d] Matrix der Multiindizes
%   xi        : [1 x d] Punkt in Zufallsraum
%   psi_val   : [1 x P] Auswertung aller Polynome am Punkt xi

P = size(multi_idx, 1);  % Anzahl Polynome
d = length(xi);          % Dimension
psi_val = ones(1, P);    % Vektor initialisieren

for k = 1:P
    for j = 1:d
        deg = multi_idx(k, j);
        psi_val(k) = psi_val(k) * hermite_he_prob(deg, xi(j));
    end
end

end
