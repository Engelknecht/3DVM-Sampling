function Postprocessing(results, COORD, SURF, elem_type, ...
                         size_xy, size_z, size_hole, ...
                         Y_nodes, a_nodes, E_nodes, ...
                         M, Phi, lambda_sorted, x, y)
% Postprocessing – Visualisiert FEM-Ergebnisse aus dem results-Struct

% --- Netz darstellen
draw_mesh(COORD, SURF, elem_type);

% --- Lastpfad (Hysterese)
figure;
plot(results.alpha, results.zeta, 'x-'); hold on;
plot(results.alpha([11 21 31 41]), results.zeta([11 21 31 41]), 'ro');
hold off;
axis tight; enlarge_axis(0.1, 0.1);
xlabel('Arbeit der äußeren Kräfte'); ylabel('Belastungsskala');
legend('alle Zeitschritte', 'visualisierte Schritte', 'Location', 'northwest');
title('Lastpfad');

% --- Verschiebungen visualisieren
U_total = sqrt(sum(results.U.^2, 1));
draw_quantity(COORD, SURF, 10 * results.U, U_total, elem_type, size_xy, size_z, size_hole);
title('Gesamtverschiebung');

% --- Abhängigkeit: Rechenzeit vs. plastische Punkte
figure;
x = results.assembly(1:results.assembly_step, 1);
y = results.assembly(1:results.assembly_step, 2);
x_ext = results.n_int;
y_ext_meas = results.assembly_elast_time;
x_r = x(x > 0); y_r = y(x > 0);

X = [ones(length(x_r), 1), x_r];
b = X \ y_r;  % lineare Regression
x_e = [x; x_ext];
y_e = [ones(length(x_e), 1), x_e] * b;
y_ext = [1, x_ext] * b;

plot(x, y, 'bx', x_e, y_e, 'b-', x_ext, y_ext, 'bo', x_ext, y_ext_meas, 'ro');
legend('plastisch – Messwerte', 'plastisch – Regression', ...
       'plastisch – Extrapolation', 'elastisch', 'Location', 'northwest');
xlabel('Anzahl plastischer Integrationspunkte');
ylabel('Rechenzeit (s)');
axis tight; enlarge_axis(0.1, 0.1);
title('Rechenzeit in Abhängigkeit von plastischer Verteilung');

% --- KLE-Feldvisualisierungen (E, Y, a)
figure; scatter3(COORD(1,:), COORD(2,:), Y_nodes, 15, Y_nodes, 'filled');
title('Fließbedingung über Knoten'); colorbar;

figure; scatter3(COORD(1,:), COORD(2,:), a_nodes, 15, a_nodes, 'filled');
title('Kinematische Verfestigung über Knoten'); colorbar;

figure; scatter3(COORD(1,:), COORD(2,:), E_nodes, 15, E_nodes, 'filled');
title('E-Modul über Knoten'); colorbar;

% --- KLE-Eigenfunktionen auf ausgewählter Linie
y_vals = unique(round(COORD(2,:), 2));
counts = arrayfun(@(yy) sum(abs(COORD(2,:) - yy) < 1e-3), y_vals);
[~, max_idx] = max(counts);
line_y = y_vals(max_idx);
line_nodes = abs(COORD(2,:) - line_y) < 1e-3;
x_line = COORD(1, line_nodes);
Phi_line = Phi(line_nodes, 1:M);

[x_sorted, idx] = sort(x_line);
Phi_sorted = Phi_line(idx, :);
styles = {'-', '--', '-.', ':'};

figure; hold on;
for i = 1:M
    plot(x_sorted, Phi_sorted(:, i), styles{mod(i-1, length(styles)) + 1}, ...
         'LineWidth', 2, ...
         'DisplayName', sprintf('\\phi_{%d} (\\lambda=%.2e)', i, lambda_sorted(i)));
end
xlabel('x'); ylabel('\phi_i(x)');
legend('show'); title(sprintf('KLE-Eigenfunktionen auf y = %.2f', line_y));
grid on;

end
