function val = evaluateMultivariateHermite(poly_coeffs, xi)
% poly_coeffs – Koeffizientenvektor (1xN) einer multivariaten Hermite
% xi          – Punkt im Zufallsraum (1xM)
    val = 1;
    dim = length(xi);
    for d = 1:dim
        val = val * polyvalHermite(poly_coeffs, xi(d));
    end
end
