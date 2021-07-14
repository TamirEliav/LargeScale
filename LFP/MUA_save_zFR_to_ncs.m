function MUA_save_zFR_to_ncs(exp_ID,forcecalc)

%% defaults
if nargin==1; forcecalc = 0; end

%% get exp info
exp = exp_load_data(exp_ID,'details','path','MUA');

%% check if file already exist
file_name = fullfile(exp.path.LFP_bands,'ripple', [exp_ID '_MUA_zFR' '.ncs']);
if exist(file_name,'file')
    if forcecalc
        % delete existing file
        warning('MUA FR file already exist and you chose to override it, deleting old file!')
        delete(file_name);
    else
        warning('MUA FR file already exist, use forcecalc to override it!');
        return;
    end
end
    
%% write
nlx_csc_write(file_name, exp.MUA.zFR, exp.MUA.t, exp.MUA.fs);

end