function [IX, IX_per_xi, x2, x2_per_ti] = get_data_in_ti(x,xi,opt)
% x - 1Xm vector of any length
% xi - nX2 matrix of start/end values
% opt - optional input, tradeoff between run-time(1) and memory(2) efficiency 
% TODO: change order of output (and change all places I used this
% function!!!!!!!!) also change 

% set default option
if ~exist('opt', 'var')
    opt=1;
end
% if user asked for IX per interval, let's just run with for loop...
if nargout>1
    opt=2;
end

switch opt
case 1
    IX = find(any(x>xi(:,1)&x<xi(:,2),1));
case 2
    for ii = 1:size(xi,1)
        IX{ii} = find(x>xi(ii,1) & x<xi(ii,2));
    end
    IX_per_xi = IX;
    IX=cat(2,IX{:});
end

if nargout>=3
    x2 = x(IX);
end

if nargout>=4
    x2_per_ti = cellfun(@(x)(x(x)), IX_per_xi, 'UniformOutput',false);
end

end