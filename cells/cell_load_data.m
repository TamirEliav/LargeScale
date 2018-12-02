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
    details = load(fullfile('L:\Analysis\Results\cells\details',[cell_ID '_cell_details']));
    cell_data.details = details.details;
end

%% spikes
if any(contains(varargin, 'spikes')) | load_all
    spikes = load(fullfile('L:\Analysis\Results\cells\spikes',[cell_ID '_cell_spikes']));
    cell_data.spikes = spikes.spikes;
end

%% RecStability
if any(contains(varargin, 'RecStability')) | load_all
    RecStability = load(fullfile('L:\Analysis\Results\cells\RecStability',[cell_ID '_cell_RecStability']));
    cell_data.RecStability = RecStability.RecStability;
end

%% FE (flight epochs)
if any(contains(varargin, 'FE')) | load_all
    FE = load(fullfile('L:\Analysis\Results\cells\FE',[cell_ID '_cell_FE']));
    cell_data.FE = FE.FE;
end

%% Firing rate map 
if any(contains(varargin, 'FR_map')) | load_all
    FR_map = load(fullfile('L:\Analysis\Results\cells\FR_map',[cell_ID '_cell_FR_map']));
    cell_data.FR_map = FR_map.FR_map;
end

%% Ipos
if any(contains(varargin, 'Ipos')) | load_all
    Ipos = load(fullfile('L:\Analysis\Results\cells\Ipos',[cell_ID '_cell_Ipos']));
    cell_data.Ipos = Ipos.Ipos;
end


end







