function  Nlx_filter_CSC2(file_IN, file_OUT, t_start_end, filter_params)

% file_IN
% file_OUT
% t_start_end
% filter_params.passband        2X1 vector of filter bandpass freq
% filter_params.fwin            in seconds (default is 2 minutes)
% filter_params.resample_fs     new sample rate
% filter_params.type            'fir1' is the default, no other types are supported now.
% filter_params.order           default is 300


%% sanity checks
if exist(file_OUT, 'file')
    error('File already exist, this code should not overwrite/append existing .ncs files');
end

%% check and set params
if ~isfield(filter_params,'order')
    filter_params.order = 300;
end
if ~isfield(filter_params,'fwin')
    filter_params.fwin = 2 * 60;
end
if ~isfield(filter_params,'type')
    filter_params.type = 'fir1';
end

%% read file
[signal, ts, fs] = Nlx_csc_read(file_IN, t_start_end);

%% filter!
Wn   = filter_params.passband / (fs/2);           
filt = fir1(filter_params.order, Wn,'bandpass' );       
signal2  = filtfilt(filt,1,signal); 
ts2 = ts;

%% re-sample the filtered data if needed
if isfield(filter_params,'resample_fs')
    fs3 = filter_params.resample_fs;
    
    % check for nyquist first!
    if ( fs3/2 < filter_params.passband(2) )
        error('resampling will create aliasing, use nyquist role!!!')
    end
    T_fs3 = 1e6/fs3;
    ts3 = ts2(1) : T_fs3 : ts2(end);
    signal3 = interp1(ts2,signal2,ts3);
end

%% write filtered data to file
if isfield(filter_params,'resample_fs')
    nlx_csc_write(file_OUT, signal3, ts3, fs3);
else
    nlx_csc_write(file_OUT, signal2, ts2, fs);
end



end





