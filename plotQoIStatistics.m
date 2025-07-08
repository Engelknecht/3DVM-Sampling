function plotQoIStatistics(coeffs, p_order, QoI_samples)

% Erwartungswert und Varianz aus PCE-Koeffizienten
mean_PCE = coeffs(1,:);                          
var_PCE  = sum(coeffs(2:end,:).^2, 1);           
std_PCE  = sqrt(var_PCE);                        

% MC-Stichproben
[N, K] = size(QoI_samples);
mean_MC = mean(QoI_samples,1);                  

% Labels für QoIs
labels = {'max. Verschiebung', 'Anz. plastische Punkte', ...
          'max. v. Mises-Spannung', 'mittlere Dehnung', ...
          'L2-Norm der Verschiebung'};

% Plot für jede QoI separat
for k = 1:K
    figure; hold on;
    
    % Balken: Erwartungswert & Std
    bar(1, mean_PCE(k), 0.4, 'FaceColor', [0.2 0.6 1], 'DisplayName', 'PCE-Erwartungswert');
    bar(2, std_PCE(k), 0.4, 'FaceColor', [1.0 0.4 0.3], 'DisplayName', 'PCE-Std');

    % Streudiagramm: Stichproben
    scatter(ones(N,1)*3, QoI_samples(:,k), 20, 'k', 'filled', ...
            'MarkerFaceAlpha', 0.3, 'DisplayName', 'Stichproben');

    % Punkt: MC-Mittelwert
    plot(4, mean_MC(k), 'x', 'MarkerSize', 12, 'MarkerFaceColor', 'y', ...
         'MarkerEdgeColor', 'k', 'DisplayName', 'MC-Mittelwert');
     
    % Achsen & Layout
    xlim([0.5 4.5]);
    xticks(1:4);
    xticklabels({'PCE-Mean', 'PCE-Std', 'Samples', 'MC-Mean'});
    ylabel('Wert');
    title(sprintf('%s – Statistik aus PCE & Stichproben (Ordnung %d)', labels{k}, p_order));
    grid on; box on;
    legend('Location','northwest');
end

end