function exp = exp_load_data(exp_ID,varargin)

%% options to load
field_options = {
    'path',         0
    'details',      1
    'pos',          1
    'flight',       1
    'flight_6m',    1
    'LM',           0
    'csc_raw_stats',1
    'wingbeat',     1
    'ripples',      1
    'MUA',          1
    'PE',           1
    'rest',         1
    'MUA_FR_map',   1
    };

%% parse input
load_all = nargin==1; % we will use this as a permitting value, i.e. load all AVAILABLE fields
if load_all
    fields2load = field_options(:,1);
else
	fields2load = varargin;
end

%%
exp = struct();
for ii_field = 1:length(fields2load)
    field_name = fields2load{ii_field};
    switch field_name
        case field_options([field_options{:,2}]==1, 1)
            file2load = fullfile('L:\Analysis\Results\exp', field_name, [exp_ID '_exp_' field_name '.mat']);
            if exist(file2load,'file')
                field_data = load(file2load);
                exp.(field_name) = field_data.(field_name);
            elseif ~load_all
                % raise error only if the user asked for this field
                % specifically, otherwise just load all the available
                % fields
%                 error(sprintf('field %s data does not exist for exp %s',field_name, exp_ID))
            end
            
        % other fields (specific handle...)
        case 'path'
            exp.path = DS_get_exp_path(exp_ID);
        case 'LM'
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
        otherwise
            error(sprintf('field name %s not recognized',field_name));
    end
end

end





