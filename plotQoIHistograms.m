function plotQoIHistograms(QoI_samples)
%PLOTQOIHISTOGRAMS  Plot histogram of each QoI sample column separately.
%
%   QoI_samples : (N x K) matrix of QoI values

labels = {'max. Verschiebung', 'Anz. plastische Punkte', ...
          'max. v. Mises-Spannung', 'mittlere Dehnung', ...
          'L2-Norm der Verschiebung'};

[~, K] = size(QoI_samples);
for k = 1:K
    figure;
    histogram(QoI_samples(:,k));
    xlabel('Wert');
    ylabel('H\xE4ufigkeit');
    title(sprintf('Verteilung von %s', labels{k}));
    grid on;
end
end

