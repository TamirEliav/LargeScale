function Nlx_detect_spikes_CSC2(dir_IN,dir_OUT,params)

data = struct();
data = get_files(data, dir_IN);
params = validate_params(params);
data = load_ref_ch(data,params);

for TT = params.TT_to_use
    data = TT_detect_spikes(data,params,TT);
end

data = check_coincidence_detection(data,params);
data = export_data(data,params);

end


%% local functions
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%
function data = arrange_files(data, dir_IN, dir_OUT)
    files_raw = dir( fullfile(dir_IN,['*_TT*_ch*']) );
    TT_files_raw = {};
    TT_ch_exist = [];
    for ii_file = 1:length(files_raw)
        file_name = files_raw(ii_file).name;
        TT_str = regexp(file_name, '_TT([\d+])','tokens','once');
        TT_str = TT_str{1};
        TT_num = str2num(TT_str);
        ch_str = regexp(file_name, '_ch([\d+])','tokens','once');
        ch_str = ch_str{1};
        ch_num = str2num(ch_str);
        TT_files_raw{TT_num,ch_num} = file_name;
    end
    TT_ch_exist = cellfun(@(x)(~isempty(x)), TT_files_raw);
    num_TTs = size(TT_files_raw,1);
    
    %% export to data struct
    data.dir_IN = dir_IN;
    data.dir_OUT = dir_OUT;
    data.files = TT_files_raw;
end
%%
function params = validate_params(params)
    if length(params.thr_uV) == 1
        params.thr_uV = repmat(params.thr_uV, size(TT_ch_exist));
    end
end
%%
function data = load_ref_ch(data,params)
    data.ref.signal = [];
    if ~isempty(params.ref_ch)
        filename = data.files{params.ref_ch(1),params.ref_ch(2)};
        [data.ref.signal,~] = Nlx_csc_read(fullfile(data.dir_IN,filename), params.t_start_end);
    end
end
%%
function data = TT_detect_spikes(data,params,TT)
    TT_data = struct();
    TT_data.TT = TT;
    TT_data = TT_load_data(data,TT_data,params);
    TT_data = TT_detect_events(data,TT_data,params);
    TT_data = TT_merge_events(data,TT_data,params);
    TT_data = TT_extract_wfrms(data,TT_data,params);
    TT_data = TT_corr_lib(data,TT_data,params);
    TT_data = TT_extract_valid(data,TT_data,params);
    % save results to data
    data.TT_data(TT) = TT_data;
end
%%
function TT_data = TT_load_data(data,TT_data,params);
    csc = {};
    timestamps = {};
    TT = TT_data.TT;
    act_ch = find(params.active_TT_channels(TT, :));
    non_act_ch = find(~params.active_TT_channels(TT, :));
    for ii_ch = act_ch
        filename = data.files{TT,ii_ch};
        [csc{ii_ch}, timestamps{ii_ch}] = Nlx_csc_read(fullfile(dir_IN,filename),params.t_start_end);
        if ~isempty(params.ref_ch) && params.ref_ch(1) ~= TT
            csc{ii_ch} = csc{ii_ch} - ref_ch_csc;
        end
    end
    % sanity check - make sure all csc files are in the same length
    if range([cellfun(@length, timestamps(act_ch))]) == 0
        timestamps = timestamps{act_ch(1)};
    else
        error('CSC file of channels from the same TT are different in length!!!')
    end
    % place zeros for disconnected channels
    for ii_ch = non_act_ch
        csc{ii_ch} = zeros(size(timestamps));
    end
    
    % save to TT_data
    TT_data.csc = csc;
    TT_data.timestamps = timestamps;
    TT_data.TT = TT;
end
%%
function TT_data = TT_detect_events(data,TT_data,params)
    
end
%%
function data = check_coincidence_detection(data,params)

end


%%
