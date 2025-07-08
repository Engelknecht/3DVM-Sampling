function H = hermite_he_prob(n, x)
% HERMITE_HE_PROB berechnet probabilistische Hermite-Polynome H_n(x)
% Rekursion: H_n(x) = x H_{n-1}(x) - (n-1) H_{n-2}(x)

if n == 0
    H = 1;
elseif n == 1
    H = x;
else
    Hm2 = 1;
    Hm1 = x;
    for k = 2:n
        H = x .* Hm1 - (k - 1) * Hm2;
        Hm2 = Hm1;
        Hm1 = H;
    end
end

end
