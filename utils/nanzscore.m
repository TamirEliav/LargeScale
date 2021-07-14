function [z,mu,sigma] = nanzscore(x,flag,dim)
% same as zscore, but can handle (ignore) nan

% [] is a special case for std and mean, just handle it out here.
if isequal(x,[]), z = x; return; end

if nargin < 2
    flag = 0;
end
if nargin < 3
    % Figure out which dimension to work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

% Validate flag
if ~(isequal(flag,0) || isequal(flag,1) || isempty(flag))
    error(message('stats:trimmean:BadFlagReduction'));
end

% Compute X's mean and sd, and standardize it
mu = nanmean(x,dim);
sigma = nanstd(x,flag,dim);
sigma0 = sigma;
sigma0(sigma0==0) = 1;
z = (x-mu) ./ sigma0;

