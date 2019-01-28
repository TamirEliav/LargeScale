function Nlx_detect_spikes_CSC(dir_IN,dir_OUT,params)

%% 
% Tamir,
% 01/2018
% Adopted from Didi, Arseny and Michael 
% 
% Here we detect the spikes for each tetrode save them.
% Output is 2 NTT files containing:
% 1. detected spikes
% 2. spikes that were detected but thrown away after library comparison

%% open log file
log_name_str = ['spikes_detection_' datestr(clock, 'yyyy-mm-dd HH-MM-SS') '.txt'];
log_name_out = fullfile(dir_OUT, log_name_str );
if ~exist(dir_OUT,'dir'), mkdir(dir_OUT); end 
diary off; diary(log_name_out); diary on

%%
% % % forcerecalc = 0;
% % % 
% % % %%
% % % if exist(dir_OUT,'dir') && ~forcerecalc
% % %     disp(['Detect-spikes was already done in: ',dir_OUT]);
% % %     return; 
% % % end

%% default params
% % % % dir_IN = 'L:\Analysis\pre_proc\0148\20170607\spikes_raw';
% % % % dir_OUT = 'L:\Analysis\pre_proc\0148\20170607\spikes_detection';
% % % % params.thr_uV = 40;
% % % % params.ref_TT = 1; % []
% % % % params.ref_ch = 1; % []
% % % % params.TT_to_use = [1 2 3 4];
% params.t_start_end = [];
% params.merge_thr_crs_width = 4;
% params.use_neg_thr = 0;
% params.lib_spike_shapes = 'library_of_acceptable_spike_shapes.mat';
% params.lib_corr_thr = 0.8;
params.nSamples = 32;
params.AlignSample = 8;

%% Get files
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

%% arrange params
if length(params.thr_uV) == 1
    params.thr_uV = repmat(params.thr_uV, size(TT_ch_exist))
end

if length(params.use_neg_thr) == 1
    params.use_neg_thr = repmat(params.use_neg_thr, size(TT_ch_exist))
end

%% init time measure and stats
time_measure = struct();
stats = struct();

%% load ref channel (if exist)
tic 
if ~isempty(params.ref_ch)
    filename = TT_files_raw{params.ref_ch(1),params.ref_ch(2)};
    [ref_ch_csc,~] = Nlx_csc_read(fullfile(dir_IN,filename),params.t_start_end);
    % remove this TT from the tetrodes list to use
%     params.TT_to_use( params.TT_to_use == params.ref_ch(1) ) = [];
end
time_measure.load_ref_ch = toc;

%% detect spikes
for TT = params.TT_to_use

    disp('')
    disp(['running detection on TT' num2str(TT)])
    
    %% sanity checks
    if sum(TT_ch_exist(TT,:)) == 0
        disp(['There are no data for TT ' TT ', skipping...'])
        continue;
    end
    
%     if params.ref_ch(1) == TT
%         disp(['TT' TT ' is using for reference, skipping...'])
%         continue;
%     end
    
    %% load raw data (CSCs)
    tic
    csc = {};
    timestamps  = {};
    act_ch = find(params.active_TT_channels(TT, :));
    non_act_ch = find(~params.active_TT_channels(TT, :));
    for ii_ch = act_ch
        filename = TT_files_raw{TT,ii_ch};
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
    time_measure.load_TT(TT) = toc;
        
    %% detect threshold crossing 
    fprintf('\t Detect spike (thr crossing) \n')
    tic
    thr_cross_vec = {};
    thr_cross_IX = {};
    SPK_start = {}; 
    SPK_End = {}; 
    
    Last_Spike_IX = 0; % Initialize
    SPK_timestamp = [];% Initialize

    for ii_ch = act_ch
        thr = params.thr_uV(TT,ii_ch);
        use_neg_thr = params.use_neg_thr(TT,ii_ch);
        % report stats (choosen thr vs. median/std)
        stats.thr_div_abs_median(ii_ch,TT) = thr / median(abs(csc{ii_ch}));
        stats.thr_div_std(ii_ch,TT) = thr / std(csc{ii_ch});
        
        if use_neg_thr
            thr_cross_IX{ii_ch} = find(csc{ii_ch}>thr | csc{ii_ch}<-thr);
        else
            thr_cross_IX{ii_ch} = find(csc{ii_ch}>thr);
        end
        
        % place '1' at threshold crossing:
        thr_cross_vec{ii_ch} = zeros(1,length(csc{ii_ch}));
        thr_cross_vec{ii_ch}(thr_cross_IX{ii_ch})=1; % place '1' at threshold crossing

        % Find the start and end segment of each spiking event on each channel of the tetrode:
        diff_thres_cross_vec{ii_ch} = diff(thr_cross_vec{ii_ch});
        SPK_start{ii_ch} = find(diff_thres_cross_vec{ii_ch} ==1)+1;
        SPK_End{ii_ch} = find(diff_thres_cross_vec{ii_ch} ==(-1))+1;
        
        if isempty(SPK_End{1})
            continue
        end
        
        %Check for the case that csc start/end point already crossed the thr
        if SPK_End{ii_ch}(1) < SPK_start{ii_ch}(1)
            SPK_End{ii_ch}(1) = [];
        end
        if SPK_start{ii_ch}(end) > SPK_End{ii_ch}(end)
            SPK_start{ii_ch}(end) = [];
        end

        % correct the spikes index according to maximum
        temp_SPK_events = []; % Initialize
        SPK_max{ii_ch} = [];
        SPK_IX{ii_ch} = []; % with regard to the maximum value of events that crossed thr
        for jj = 1:length(SPK_End{ii_ch})
            temp_IX_vec = SPK_start{ii_ch}(jj):SPK_End{ii_ch}(jj);
            temp_SPK_events = csc{ii_ch}(temp_IX_vec);
            % using abs for the negative thr option (if no negative thr 
            % used it will not affect the positive, unless there is a sharp 
            % transition from positive to negative or vice versa without 
            % passing below abs thr)
            [temp_SPK_max,IX] = max(abs(temp_SPK_events)); 
%             [temp_SPK_max,IX] = max(temp_SPK_events);
            max_IX = temp_IX_vec(IX); % Correct for the real timestamp relative to the entire recording session
            SPK_IX{ii_ch}(1,jj) = max_IX;
        end
    end
    clear diff_thres_cross_vec
    clear thres_cross_vec
    clear thres_cross_IX
    time_measure.thr_crs(TT) = toc;
	stats.detect_nEvents(:,TT) = cellfun(@length,SPK_IX);
    
    %% Merge spikes from different channels
    % Merge spikes from all channels of the same TT
    % If two spikes are detected on differnt channel but they are close
    % enough in time (<params.merge_thr_crs_width), we merge them to a single
    % spike
    fprintf('\t Merge spikes from different channels \n')
    tic
    
    % set '1' at events per ch and sum them up
    events_all_ch = zeros(4,length(csc{1}));
    for ii_ch = 1:4
        events_all_ch( ii_ch, SPK_IX{ii_ch} ) = 1;
%         events_all_ch( ii_ch, SPK_IX{ii_ch} ) = abs(csc{ii_ch}(SPK_IX{ii_ch}));
    end
    events_ch_combined = sum(events_all_ch);
    clear events_all_ch;
    % convolve with a rectangular kernal with width of the merging event 
    rect_kernel = ones(1,params.merge_thr_crs_width);
    events_ch_combined_filtered = filtfilt(rect_kernel, 1, events_ch_combined);
    % find the peaks in the convolved signal, representing the different
    % merged events (across channels)
    [~,SPK_IXs_All_Ch_combined_sorted] = findpeaks(events_ch_combined_filtered);
    
    clear events_ch_combined
    time_measure.merge_ch(TT) = toc;
    stats.merge_nEvents(TT) = length(SPK_IXs_All_Ch_combined_sorted);
    
    %% Extract spike waveforms
    % Extract the waveform and timestamps of each spike (as in Neuralynx, the peak
    % will be the 8th sample out of over all 32 samples of
    %the spike shape vec.
    
    fprintf('\t Extract waveforms \n')
    tic
    
    if (SPK_IXs_All_Ch_combined_sorted(1) < params.AlignSample)
        SPK_IXs_All_Ch_combined_sorted(1) = [];
    end
    if (SPK_IXs_All_Ch_combined_sorted(end) + params.nSamples-params.AlignSample > length(timestamps))
        SPK_IXs_All_Ch_combined_sorted(end) = [];
    end
    
    SPK_waveforms = zeros(params.nSamples,4,length(SPK_IXs_All_Ch_combined_sorted));% Initialize
    % IX relative to spikes peak
    trigger_IX = repmat( [(1-params.AlignSample):(params.nSamples-params.AlignSample)]',1,length(SPK_IXs_All_Ch_combined_sorted));
    % IX relative to csc signal
    trigger_IX = trigger_IX + repmat(SPK_IXs_All_Ch_combined_sorted,[params.nSamples 1]);
    for ch = 1:4
        SPK_waveforms(:,ch,:) = csc{ch}(trigger_IX);
    end
    % get spikes ts
    SPK_timestamp = timestamps(SPK_IXs_All_Ch_combined_sorted);
    
    time_measure.extract_wvfrm(TT) = toc;
    stats.extract_nWvfrm(:,TT) = size(SPK_waveforms);
    stats.extract_nTimestamps(TT) = length(SPK_timestamp);
    
    %% library of acceptable spike shapes
    % Clean artifacts using the library of acceptable spike shapes.
    % Now we will throw away noisy theshold crossing using the library of
    % acceptable spike shapes as reference

    fprintf('\t library of acceptable spike shapes \n')
    tic
    load(params.lib_spike_shapes);
    
    % For each event take the channel with the largest peak (8th point)
    [~,max_ch_IX] = max(squeeze(abs(SPK_waveforms(params.AlignSample,:,:))),[],1);
    spikes_waveforms = zeros(size(SPK_waveforms,1),size(SPK_waveforms,3));
    for ch=1:4
        IX = find(max_ch_IX == ch);
        spikes_waveforms(:,IX) = squeeze(SPK_waveforms(:,ch,IX));
    end
    
    % calc corr
    xxx_lags_shifts = [1:30; 2:31; 3:32];
    ccc = [];
    rrr = [];
    for ii_shift = 1:size(xxx_lags_shifts,1)
        xxx_lags = xxx_lags_shifts(ii_shift, :);
        ccc = corr(spikes_waveforms(xxx_lags,:), library_of_acceptable_spike_shapes(:,2:end-1)');
        rrr(ii_shift,:) = max(ccc,[],2);
    end
    rrr = max(rrr,[],1);
    vector_of_accepted_spikes = ( rrr >=  params.lib_corr_thr );
    
    figure('Units','normalized','Position',[0 0 1 1]);
    hold on
    h= histogram(rrr);
    h.NumBins = h.NumBins * 5;
    plot(repelem(params.lib_corr_thr,2), get(gca,'ylim'), 'r', 'LineWidth',2);
    xlabel('max lib corr')
    ylabel('Counts')
    title(sprintf('TT %d',TT));
    lib_corr_figname = fullfile(dir_OUT, sprintf('lib_corr_TT_%d',TT));
    saveas(gcf,lib_corr_figname,'tif')
    close gcf
    
    time_measure.lib_corr(TT) = toc;
    stats.lib_nAccepted(TT) = sum(vector_of_accepted_spikes);
    stats.lib_nNotAccepted(TT) = sum(~vector_of_accepted_spikes);
    
    %% get spikes data for accepted spikes
    % TODO: enable saving the non-accepted waveforms (or at least report
    % how many)
    % Find the IXs of the accepted spike waveforms and store the accepted and
    % not-accepted waveforms speratly.
    tic
%     Spike_waveforms_accepted = SPK_waveforms(:,vector_of_accepted_spikes,:);
%     Spike_waveforms_NO_accepted = SPK_waveforms(:,~vector_of_accepted_spikes,:);
    
    Timestamps_accepted_spikes_TT{TT} = SPK_timestamp(vector_of_accepted_spikes);
    spikes_TT{TT} = SPK_waveforms(:,:,vector_of_accepted_spikes);
    time_measure.divide_lib_corr_results(TT) = toc;
    stats.extract_lib_accepted_nWvfrm(:,TT) = size(spikes_TT{TT});
    stats.extract_lib_accepted_nTimestamps(TT) = length(Timestamps_accepted_spikes_TT{TT});

    %% Clear un-needed variables
    clear IX_NO_accepted IX_NO_accepted_sleep IX_accepted_sleep
    clear SPK_All_Ch_combined
    clear Spike_waveforms_NO_accepted
    clear Spike_waveforms_NO_accepted_sleep Spike_waveforms_accepted_sleep Timestamps_accepted_spikes_sleep
    clear VT VT_Parameters ccc spikes_sleep
    clear vector_of_accepted_spikes vector_of_accepted_spikes_sleep
    clear csc
%     clear timestamps 

end %end looping over Tetrodes


%% Coincidence detection - Tamir's version
do_coincidence_detection = 1;
tic
if do_coincidence_detection && (length(params.TT_to_use) > 1)

    disp('-------------------------------------------')
    disp('Coincidence-Detection (eliminate artifacts)')
    
    %% find the invalid ts (with CD)
    % first, mark 1 for event for each TT
    TTs_events = zeros(length(Timestamps_accepted_spikes_TT), length(timestamps));
    for ii_TT = 1:length(Timestamps_accepted_spikes_TT)
        spikes_ts_IX = ismember(timestamps,Timestamps_accepted_spikes_TT{ii_TT});
        TTs_events(ii_TT,spikes_ts_IX) = 1;
    end
    mask = conv(sum(TTs_events), ones(1,params.CD_detect_win_len), 'same'); % sum across TT and convolve with CD detection window
    mask = mask >= params.CD_n_TT_thr; % apply min number of TT for detection
    mask = conv(mask, ones(1,params.CD_invalid_win_len), 'same'); % convolve dections with the invalidation window (we are after detection!)
    TTs_events(:,mask>0) = 0; % invalidate events
    
    % Last, update accepted spikes (ts+shape)
    for ii_TT = 1:length(params.TT_to_use)
        TT = params.TT_to_use(ii_TT);
        % find the invalid spikes in each TT
        old_spikes_ts = Timestamps_accepted_spikes_TT{TT};
        new_spikes_ts = timestamps( TTs_events(TT,:)==1 );
        invalid_spikes = ~ismember(old_spikes_ts, new_spikes_ts);
        Timestamps_accepted_spikes_TT{TT}(invalid_spikes) = [];
        spikes_TT{TT}(:,:,invalid_spikes) = [];
        stats.CD_nEventRemove(TT) = length(invalid_spikes);
    end
    
end
time_measure.CD = toc;

%%
stats.rec_len.nSamples = length(timestamps);
stats.rec_len.totalTimeHours = range(timestamps)*1e-6/60/60;

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the data in Neuralynx NTT files for three cases:
% (1) All the data.
tic
disp('-------------------------')
disp('Write results (NTT files)')
for TT = params.TT_to_use

    % make sure there are data from this TT
    if sum(TT_ch_exist(TT,:)) == 0
        disp(['There are no data for TT ' TT ', skipping...'])
        continue;
    end
    
    CSC_filename = TT_files_raw{TT,find( TT_ch_exist(TT,:)==1, 1 ,'first' )};
    CSC_filename = regexprep(CSC_filename, '_ch([\d+])', '');
    NTT_filename = regexprep(CSC_filename, '.ncs', '.NTT');
    
    filename_out = fullfile(dir_OUT,NTT_filename);
   
    header_file = 'Nlx_header_NTT.txt';
    header = textread(header_file, '%s', 'delimiter', '\n', 'whitespace', '');

    Timestamps = Timestamps_accepted_spikes_TT{TT};
    CellNumbers = zeros(size(Timestamps));
    Samples = spikes_TT{TT};
    fs = 1e6/median(diff(timestamps));
    ADMaxValue = 32767;
    InputRange = max(Samples(:));
    ADC = InputRange / ADMaxValue / 1e6;
    Samples = Samples ./ ADC ./ 1e6;
    ADC_str = sprintf('%.24f',ADC);
    InputRange_str = sprintf('%g',InputRange);
    ADC_str_IX = contains(header, 'ADBitVolts');
    InputRange_str_IX = contains(header, 'InputRange');
    header{ADC_str_IX} = sprintf('-ADBitVolts %s %s %s %s', ADC_str, ADC_str, ADC_str, ADC_str);
    header{InputRange_str_IX} = sprintf('-InputRange %s %s %s %s', InputRange_str, InputRange_str, InputRange_str, InputRange_str);
    header{contains(header, '-SamplingFrequency')} = sprintf('-SamplingFrequency %g',fs);
    
    Mat2NlxSpike(filename_out, 0, 1, [], [1 0 1 0 1 1], ...
        Timestamps, CellNumbers, Samples, header);
    
end % end looping over tetrodes
time_measure.write_NTTs = toc;

%% report timing
time_measure_fields = fieldnames(time_measure);
total_runtime = 0;
for ii_field = 1:length(time_measure_fields)
    total_runtime = total_runtime + sum(time_measure.(time_measure_fields{ii_field}));
end
time_measure.total_runtime = total_runtime;
fn_structdisp(time_measure)

%% report used parmas
fn_structdisp(stats)

%% report results stats (num spikes detected / CD / lib / ...)
fn_structdisp(stats)

%% save some meta-data to .mat file (params,stats)
params_fileout = fullfile(dir_OUT,'params');
save(params_fileout, 'params');
stats_fileout = fullfile(dir_OUT,'stats');
save(stats_fileout, 'stats');
time_measure_fileout = fullfile(dir_OUT,'time_measure');
save(time_measure_fileout, 'time_measure');

%% FIN
disp('spikes detection FINISH!!!')

%% close log file
diary off





