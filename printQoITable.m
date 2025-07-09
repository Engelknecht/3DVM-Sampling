function T = printQoITable(QoI_samples)
%PRINTQOITABLE  Display QoI sample results as MATLAB table.
%
%  QoI_samples : [N x 5] matrix of QoI values for N samples
%
%  T : returned table of results

labels = {'Sample','u_max','n_plast','sigma_vm_max','strain_mean','norm_u'};
T = array2table([(1:size(QoI_samples,1))' , QoI_samples], ...
                'VariableNames', labels);
disp(T);
end
