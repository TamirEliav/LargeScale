function exp = exp_load_data(varargin)

%% parse input
if nargin == 0
    error('exp_ID not specified')
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

%% pos_y_std
if any(contains(varargin, 'pos_y_std')) | load_all
    pos_y_std = load(fullfile('L:\Analysis\Results\exp\pos_y_std',[exp_ID '_exp_pos_y_std']));
    exp.pos_y_std = pos_y_std.pos_y_std;
end

%% LM (landmarks)
if any(contains(varargin, 'LM')) | load_all
    %% TODO: temp, make it pretty....
    exp_path = exp_load_data(exp_ID,'path');
    load(exp_path.path.calib_tunnel_file);
    clear LM
    LM_file = 'L:\Analysis\Code\calib\Landmarks.xlsx';
    T = readtable(LM_file);
    pos_proj = POS_calc_linearized([T.pos_X T.pos_Y], calib_tunnel);
    T.pos_proj = pos_proj(:,1);
    T.pos_proj_y = pos_proj(:,2);
    LM = table2struct(T);
    exp.LM = LM;
end


end





