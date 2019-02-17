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
[signal, ts, fs, params] = Nlx_csc_read(file_IN, t_start_end);

%% filter!
switch filter_params.type 
    case 'fir1'
        Wn   = filter_params.passband / (fs/2);
        b = fir1(filter_params.order, Wn,'bandpass');
        a = 1;
        signal = filtfilt(b,a,signal);
    case 'highpassfir1'
        Wn   = filter_params.passband / (fs/2);
        b = fir1(filter_params.order, Wn,'high');
        a = 1;
        signal = filtfilt(b,a,signal);
    case 'butter'
        Wn   = filter_params.passband / (fs/2);
        [b,a] = butter(filter_params.order,Wn,'bandpass');
        signal = filtfilt(b,a,signal);
    case 'firpm'
        Wn   = filter_params.passband / (fs/2);
        freq = [0 filter_params.passband fs/2];
        amp  = [0 1 1 0   ];
        Wn   = freq / (fs/2);
        b = firpm(filter_params.order,Wn,amp);
        a = 1;
        signal = filtfilt(b,a,signal);
    case 'custom'
        custom_filt = designfilt(filter_params.designfilt_input{:},'SampleRate', fs);
        signal = filtfilt(custom_filt,signal);
    otherwise
            error('no such filter type supported!')
end

%% re-sample the filtered data if needed
if isfield(filter_params,'resample_fs')
    fs_resample = filter_params.resample_fs;
    
    % check for nyquist first!
    if ( fs_resample/2 < filter_params.passband(2) )
        error('resampling will create aliasing, use nyquist rule!!!')
    end
    dt_resample = 1e6/fs_resample;
    ts_resample = ts(1) : dt_resample : ts(end);
    signal = interp1(ts,signal,ts_resample);
end

%% write filtered data to file
% TODO: add filter params to NLX header (as a comment?!)
if isfield(filter_params,'resample_fs')
    nlx_csc_write(file_OUT, signal, ts_resample, fs_resample, params.header)
else
    nlx_csc_write(file_OUT, signal, ts, fs, params.header)
end



end





