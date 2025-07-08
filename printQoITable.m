function printQoITable(QoI_samples)
%PRINTQOITABLE  Nicely print QoI sample results in tabular form.
%
%  QoI_samples : [N x 5] matrix of QoI values for N samples

labels = {'Sample','u_{max}','n_{plast}','sigma_{vm,max}','strain_{mean}','||u||_2'};
header = sprintf('%6s | %12s %12s %12s %12s %12s\n', labels{:});
sep    = repmat('-',1, length(header));
fprintf('\n%s%s', header, sep);

for i = 1:size(QoI_samples,1)
    fprintf('%6d | %12.4g %12.0f %12.4g %12.4g %12.4g\n', i, QoI_samples(i,1), ...
            QoI_samples(i,2), QoI_samples(i,3), QoI_samples(i,4), QoI_samples(i,5));
end
end