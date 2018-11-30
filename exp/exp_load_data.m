function exp = exp_load_data(varargin)

%% parse input
if nargin == 0
    error('cell ID not specified')
    return;
end
load_all = nargin==1;
exp_ID = varargin{1};
exp = struct();

%% path
if any(contains(varargin, 'path')) | load_all
    exp.path = DS_get_exp_path(exp_ID);
end

%% details
if any(contains(varargin, 'details')) | load_all
    details = load(fullfile('L:\Analysis\Results\exp\details',[exp_ID '_exp_details']));
    exp.details = details.details;
end

%% position
if any(contains(varargin, 'position')) | load_all
    position = load(fullfile('L:\Analysis\Results\exp\position',[exp_ID '_exp_position']));
    exp.pos = position.pos;
end

%% flight
if any(contains(varargin, 'flight')) | load_all
    flight = load(fullfile('L:\Analysis\Results\exp\flight',[exp_ID '_exp_flight']));
    exp.flight = flight.flight;
end


end