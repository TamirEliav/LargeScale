function IX = get_data_in_ti(t,ti,opt)
% t - 1Xm vector of any length
% ti - nX2 matrix of start/end timestamps
% opt - optional input, tradeoff between run-time(1) and memory(2) efficiency 

if ~exist('opt', 'var')
    opt=1;
end

switch opt
    case 1
        IX = find(any(t>ti(:,1)&t<ti(:,2),1));
    case 2
        for ii = 1:size(ti,1)
            IX{ii} = find(t>ti(ii,1) & t<ti(ii,2));
        end
        IX=cat(2,IX{:});
end


end