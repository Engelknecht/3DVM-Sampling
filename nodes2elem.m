function [E_elem, Y_elem, a_elem] = nodes2elem(ELEM, E_nodes, Y_nodes, a_nodes)                                       
% -----------------------------------------------------------------------

% --- basic checks -------------------------------------------------------
if nargin ~= 4
    error('nodes2elem needs exactly four input arguments.');
end

nElem = size(ELEM,2);

E_elem = zeros(1,nElem);
Y_elem = zeros(1,nElem);
a_elem = zeros(1,nElem);

% --- loop over elements -------------------------------------------------
for k = 1:nElem
    idx = ELEM(:,k);           % global node indices of element k
    E_elem(k) = mean(E_nodes(idx));
    Y_elem(k) = mean(Y_nodes(idx));
    a_elem(k) = mean(a_nodes(idx));
end
end