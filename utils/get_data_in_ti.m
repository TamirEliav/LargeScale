function [IX, IX_per_ti, t2, t2_per_ti] = get_data_in_ti(t,ti,opt)
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
    IX = find(any(t>ti(:,1)&t<ti(:,2),1));
case 2
    IX = {}; % init
    for ii = 1:size(ti,1)
        IX{ii} = find(t>ti(ii,1) & t<ti(ii,2));
    end
    IX_per_ti = IX;
    IX=cat(2,IX{:});
end

if nargout>=3
    t2 = t(IX);
end

if nargout>=4
    t2_per_ti = cellfun(@(x)(t(x)), IX_per_ti, 'UniformOutput',false);
end

end