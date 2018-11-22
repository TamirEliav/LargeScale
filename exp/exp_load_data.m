function exp_data = exp_load_data(varargin)

%% parse input
if nargin == 0
    error('cell ID not specified')
    return;
end
load_all = nargin==1;
exp_ID = varargin{1};
exp_data=struct();

%% details
if any(contains(varargin, 'details')) | load_all
    details = load(fullfile('L:\Analysis\Results\exp\details',[exp_ID '_exp_details']));
    exp_data.details = details.details;
end

%% position
% if any(contains(varargin, 'pos')) | load_all
%     position = load(fullfile('L:\Analysis\Results\exp\position',[exp_ID '_exp_position']));
%     exp_data.pos = position.position;
% end


end