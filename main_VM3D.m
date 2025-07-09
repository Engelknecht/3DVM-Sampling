%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  main_VM3D.m – Top‑Level Script (modular Version)           %
%  ---------------------------------------------------------  %
%  Calls the building‑block functions created from your       %
%  original monolithic script.  Each step is now a clear      %
%  function call so dass du einzelne Teile leichter testen,   %
%  auslagern oder ersetzen kannst.                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% 1 | Parameter‑Struktur laden
[elem_type, level, mu_E, mu_a, mu_Y, poisson,               ...
 traction_force, size_xy, size_z, size_hole,                ...
 M, lc, scale_E, scale_Y, scale_a,                          ...
 rho_EY, rho_Ea, rho_Ya,                                    ...
 p_order, Xi_dim, xi_cdf] = setupParameter();

%% 2 | 3‑D‑Mesh erzeugen (bleibt dein alter Generator)
switch elem_type
    case 'P1'; [COORD,ELEM,SURF,NEUMANN,Q] = mesh_P1(level,size_xy,size_z,size_hole);
    case 'P2'; [COORD,ELEM,SURF,NEUMANN,Q] = mesh_P2(level,size_xy,size_z,size_hole);
    case 'Q1'; [COORD,ELEM,SURF,NEUMANN,Q] = mesh_Q1(level,size_xy,size_z,size_hole);
    case 'Q2'; [COORD,ELEM,SURF,NEUMANN,Q] = mesh_Q2(level,size_xy,size_z,size_hole);
    otherwise ;  error('Unknown element type »%s«',elem_type);
end
n_n   = size(COORD,2);                % # nodes
%points = COORD(1:2,:)';               % (n_n × 2) – nur x,y‑Koordinaten
x = COORD(1,:)';     % x Koordinate 
y = COORD(2,:)';     % y Koordinate
points = [x y];

%% 3 | Kovarianzmatrix aufbauen
[CovMat, ~] = assembleCovariance(points, lc, n_n);

%% 4 | Karhunen–Loève‐Entwicklung
[Phi, lambda_sorted] = solveKLE(COORD, CovMat, M, rho_EY, rho_Ea, rho_Ya, n_n,x,y);

%% 5 | Zufällige Feld‑Realisierung gemäß multivariater KLE
Sigma_kle = [eye(M), rho_EY*eye(M), rho_Ea*eye(M);
             rho_EY*eye(M), eye(M), rho_Ya*eye(M);
             rho_Ea*eye(M), rho_Ya*eye(M), eye(M)];
L_kle = chol(Sigma_kle, 'lower');

%% 5 | Collocation Punkte für PCE
[Xi_samples, w_samples] = generateCollocationPoints(Xi_dim, p_order);   % Gauss-Hermite Punkte
P = size(Xi_samples,1);
% Visualisiere die Collocation-Punkte für jede Zufallsvariable separat
plotCollocationPoints(Xi_samples);

QoI_samples = zeros(P, 5);  % Speicher für 5 QoIs
results = struct();         % letztes erfolgreiches Resultat

for p = 1:P
    fprintf('\n--- Berechne Stichprobe %d von %d ---\n', p, P);
    xi = Xi_samples(p,:)';

    % 5.1 | Feld-Realisierung
    xi_all = L_kle * xi;
    xi_E = xi_all(1:M);
    xi_Y = xi_all(M+1:2*M);
    xi_a = xi_all(2*M+1:end);
    [E_nodes, Y_nodes, a_nodes] = evaluateKLE_samplewise(M, Phi, lambda_sorted, ...
                             mu_E, mu_Y, mu_a, scale_E, scale_Y, scale_a, ...
                             n_n, xi_E, xi_Y, xi_a);

    % 5.2 | Knoten → Elemente
 [E_elem, Y_elem, a_elem] = nodes2elem(ELEM, E_nodes, Y_nodes, a_nodes);

    if any(E_elem < 1e3) || any(Y_elem <= 0) || any(a_elem < 0)
        warning('Sample %d: Unphysical material properties -> skipped.', p);
        continue;
    end


    % 5.3 | FEM-Simulation
    sample_results = runLoadStep(COORD, SURF, ELEM, NEUMANN, Q, ...
                       elem_type, poisson, traction_force, ...
                       E_elem, Y_elem, a_elem, ...
                       size_xy, size_z, size_hole, p);

    if isfield(sample_results, 'failed') && sample_results.failed
        continue; % oder QoI_samples(p,:) = NaN;
    end

    results = sample_results;  % nur bei Erfolg überschreiben

    % 5.4 | QoI-Auswertung
    QoI_samples(p,1) = max(vecnorm(results.U));           % maximale Verschiebung
    QoI_samples(p,2) = sum(results.n_plast);              % Anzahl plastischer Punkte
    QoI_samples(p,3) = max(results.stress_vm);            % max. v. Mises-Spannung
    QoI_samples(p,4) = mean(results.strain_eqv);          % mittlere Äquivalenzdehnung
    QoI_samples(p,5) = norm(results.U);                   % L2-Norm der Verschiebung
  

     end
[Xi_samples, w_samples, alpha, Psi_eval, PsiSqNorm] = buildPCEBasis(Xi_dim, p_order);

%% 6 | PCE-Approximation auswerten (z. B. über Regression oder quadr. Projekt)
coeffs = computePCEcoeffs(QoI_samples, Xi_samples, p_order);

% Tabelle der QoI-Stichproben ausgeben
printQoITable(QoI_samples);
% Histogramme jeder QoI-Spalte anzeigen
plotQoIHistograms(QoI_samples);

%% 6 | PCE-Approximation auswerten (z. B. über Regression oder quadr. Projekt)
coeffs = computePCEcoeffs(QoI_samples, Xi_samples, p_order);



%% 7 | Visualisierung
plotQoIStatistics(coeffs, p_order, QoI_samples);

%% 8 | Post‑Processing (Plots, Statistiken …)
Postprocessing(results, COORD, SURF, elem_type, ...
                         size_xy, size_z, size_hole, ...
                         Y_nodes, a_nodes, E_nodes, ...
                         M, Phi, lambda_sorted, x, y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Ende main_VM3D.m                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
