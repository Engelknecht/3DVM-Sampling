function alpha = multi_index(M, p_order)
% MULTI_INDEX - Erzeugt einen Multiindex für M-dimensionale Polynome bis zur Ordnung p_order

alpha = cell(p_order + 1, 1);
alpha{1} = zeros(1, M);  % Ordnung 0

switch M
    case 1
        for q = 1:p_order
            alpha{q + 1} = q;
        end
    otherwise
        for q = 1:p_order
            s  = nchoosek(1:M + q - 1, M - 1);
            s1 = zeros(size(s, 1), 1);
            s2 = (M + q) + s1;
            alpha{q + 1} = flipud(diff([s1 s s2], 1, 2)) - 1;
        end
end

alpha = cell2mat(alpha);  % als Matrix zurückgeben

end
