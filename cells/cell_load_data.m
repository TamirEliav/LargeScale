function cell_data = cell_load_data(varargin)

%% parse input
if nargin == 0
    error('cell ID not specified')
    return;
end
load_all = nargin==1;
cell_ID = varargin{1};
cell_data=struct();

%% details
if any(contains(varargin, 'details')) | load_all
    details = load(fullfile('L:\Analysis\Results\Cells\details',[cell_ID '_cell_details']));
    cell_data.details = details.details;
end

%% spikes
if any(contains(varargin, 'spikes')) | load_all
    spikes = load(fullfile('L:\Analysis\Results\Cells\spikes',[cell_ID '_cell_spikes']));
    cell_data.spikes = spikes.spikes;
end


end