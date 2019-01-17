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
% % % % params.TT_to_use = [4];
% params.t_start_end = [];
% params.merge_thr_crs_width = 4;
% params.use_neg_thr = 0;
% params.lib_spike_shapes = 'library_of_acceptable_spike_shapes.mat';
% params.lib_corr_thr = 0.8;

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

%% arrange data/params
if length(params.thr_uV) == 1
    params.thr_uV = repmat(params.thr_uV, size(TT_ch_exist))
end

%% load ref channel (if exist)
if ~isempty(params.ref_ch)
    filename = TT_files_raw{params.ref_ch(1),params.ref_ch(2)};
    [ref_ch_csc,~] = Nlx_csc_read(fullfile(dir_IN,filename),params.t_start_end);
    % remove this TT from the tetrodes list to use
%     params.TT_to_use( params.TT_to_use == params.ref_ch(1) ) = [];
end


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
    csc = {};
    timestamps  = {};
    % TODO: for disconnected channels place zeros as the signal
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
    
        
    %% detect threshold crossing 
    disp('------------------------------------')
    disp('Detect spike (thr crossing)         ')
    tic
    thr_cross_vec = {};
    thr_cross_IX = {};
    SPK_start = {}; 
    SPK_End = {}; 
    
    Last_Spike_IX = 0; % Initialize
    SPK_timestamp = [];% Initialize

    % TODO: run only over valid channels
    for ii_ch = act_ch
        thr = params.thr_uV(TT,ii_ch);
        thr_cross_vec{ii_ch} = zeros(1,length(csc{ii_ch})); % TODO: initiate first!
        if params.use_neg_thr
            thr_cross_IX{ii_ch} = find(csc{ii_ch}>thr | csc{ii_ch}<-thr);
        else
            thr_cross_IX{ii_ch} = find(csc{ii_ch}>thr);
        end;
        
        % place '1' at threshold crossing:
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
            [temp_SPK_max,IX] = max(temp_SPK_events);
            max_IX = temp_IX_vec(IX);
            SPK_IX{ii_ch}(1,jj) = max_IX; % Correct for the real timestamp relative to the entire recording session
        end
        
    end
    clear diff_thres_cross_vec
    clear thres_cross_vec
    clear thres_cross_IX
    toc
  
    
    %% Merge spikes from different channels
    % Merge spikes from all channels of the same TT
    % If two spikes are detected on differnt channel but they are close
    % enough in time (<params.merge_thr_crs_width), we merge them to a single
    % spike
    disp('------------------------------------')
    disp('Merge spikes from different channels')
    tic
    
% % % % % % % % % % % % % % % % % % % % % % % %     
    events_all_ch = zeros(4,length(csc{1}));
    for ii_ch = 1:4
        events_all_ch( ii_ch, SPK_IX{ii_ch} ) = 1;
    end
    events_ch_combined = sum(events_all_ch);
    clear events_all_ch;
    toc
% % % % % % % % % % % % % % % % % % % % % % % %     
    rect_kernel = ones(1,params.merge_thr_crs_width);
    events_ch_combined_filtered = filtfilt(rect_kernel, 1, events_ch_combined);
    toc
% % % % % % % % % % % % % % % % % % % % % % % %     
    [~,SPK_IXs_All_Ch_combined_sorted] = findpeaks(events_ch_combined_filtered);
%     [~,SPK_IXs_All_Ch_combined_sorted] = findpeaks(events_ch_combined, 'MinPeakDistance', param.spikes.x_sep_spike_thres);
    clear events_ch_combined
    toc
% % % % % % % % % % % % % % % % % % % % % % % %     

    toc
    
    %% Extract spike waveforms
    % Extract the waveform and timestamps of each spike (as in Neuralynx, the peak
    % will be the 8th sample out of over all 32 samples of
    %the spike shape vec.
    
    disp('------------------------------------')
    disp('Extract waveforms')
    tic
    
    SPK_waveforms = zeros(4,length(SPK_IXs_All_Ch_combined_sorted),32);% Initialize
    current_file_spike_counter = 0;
    for ii_spike = 1:length(SPK_IXs_All_Ch_combined_sorted)
        current_spike_max_IX = SPK_IXs_All_Ch_combined_sorted(ii_spike);
        if ((current_spike_max_IX+24<=length(timestamps))&& (current_spike_max_IX-7>0)) % i.e., we are NOT cutting the spike in the middle, TODO: no need to check this every loop, we can check only for the first and last spikes....
            current_file_spike_counter = current_file_spike_counter + 1;
            SPK_timestamp(1,Last_Spike_IX+current_file_spike_counter) = timestamps(current_spike_max_IX);
            for curr_chan_idx=1:4
                SPK_waveforms(curr_chan_idx,Last_Spike_IX+current_file_spike_counter,:) = csc{curr_chan_idx}(current_spike_max_IX-7:1:current_spike_max_IX+24);
            end
        end
    end
    Last_Spike_IX = Last_Spike_IX + current_file_spike_counter;
    
    toc
    
    %% library of acceptable spike shapes
    % Clean artifacts using the library of acceptable spike shapes.
    % Now we will throw away noisy theshold crossing using the library of
    % acceptable spike shapes as reference

    disp('------------------------------------')
    disp('library of acceptable spike shapes')
    tic
    
    % First we will need to normalize peak amplitude of each spike to a value
    % of '1' such that we can comapre it to the library of acceptable spike
    % shapes:
    
    load(params.lib_spike_shapes);
    vector_of_accepted_spikes = zeros( 1, length(SPK_waveforms) ) + NaN ; % Initialize
    vector_of_max_r_values = zeros( 1, 1 ) + NaN ;
    
    % For each event take the channel with the largest peak (8th point)
    spikes_waveforms = zeros(size(SPK_waveforms,2),32);
    for ii_spike = 1:size(SPK_waveforms,2)
        [~,ii_ch_max] = max(SPK_waveforms(:,ii_spike,8));
        spikes_waveforms(ii_spike,:) = SPK_waveforms(ii_ch_max,ii_spike,:);
    end
    
    % calc corr
    xxx_lags_shifts = [1:30; 2:31; 3:32];
    ccc = [];
    rrr = [];
    for ii_shift = 1:size(xxx_lags_shifts,1)
        xxx_lags = xxx_lags_shifts(ii_shift, :);
        ccc = corr(spikes_waveforms(:,xxx_lags)', library_of_acceptable_spike_shapes(:,2:end-1)');
        rrr(ii_shift,:) = max(ccc,[],2);
    end
    rrr = max(rrr,[],1);
    
    vector_of_accepted_spikes = ( rrr >=  params.lib_corr_thr ); % TODO: plot hist of 'r' values + thr

    % Find the IXs of the accepted spike waveforms and store the accepted and
    % not-accepted waveforms speratly.
    IX_accepted = find(vector_of_accepted_spikes ==1);
    IX_NO_accepted = find(vector_of_accepted_spikes ==0);
    
    
    toc
    
    
    
    %%
    % Extract the accepted (and not accepted) spikes and define new varialbes:
    Spike_waveforms_accepted = zeros(4,length(IX_accepted),32);
    Spike_waveforms_NO_accepted = zeros(4,length(IX_NO_accepted),32);
    
    counter_accepted = 0;
    counter_NO_accepted = 0;
    for ii = 1:length(vector_of_accepted_spikes)
        if vector_of_accepted_spikes(ii) == 1;
            counter_accepted = counter_accepted + 1;
            
            for curr_chan_idx=1:4
                Spike_waveforms_accepted(curr_chan_idx,counter_accepted,:) = SPK_waveforms(curr_chan_idx,ii,:);
            end
        else
            counter_NO_accepted = counter_NO_accepted + 1;
            for curr_chan_idx=1:4
                Spike_waveforms_NO_accepted(curr_chan_idx,counter_NO_accepted,:) = SPK_waveforms(curr_chan_idx,ii,:);
            end
        end
    end
    
    
    %%
    Timestamps_accepted_spikes = SPK_timestamp(IX_accepted);
    %rotate the timestamps to be 1 x num_records
    %timestamps = rot90(Timestamps_accepted_spikes);
    
    % Clear un-needed variables
    clear IX_NO_accepted IX_NO_accepted_sleep IX_accepted_sleep
    clear SPK_All_Ch_combined
    clear Spike_waveforms_NO_accepted
    clear Spike_waveforms_NO_accepted_sleep Spike_waveforms_accepted_sleep Timestamps_accepted_spikes_sleep
    clear VT VT_Parameters ccc spikes_sleep
    clear vector_of_accepted_spikes vector_of_accepted_spikes_sleep
    clear csc
%     clear timestamps 

    [numCh, numRec ~] = size(Spike_waveforms_accepted);
    spikes = zeros(32,4,numRec);
    for rec=1:numRec
        %rec/numRec
        for channel = 1:numCh
            current_channel_waveform = squeeze(Spike_waveforms_accepted(channel,rec,:));
            %for point=1:32
            %             spikes(point, channel, rec) = current_channel_waveform(point)*amplitude_factor;
            %end
%             spikes(:, channel, rec) = current_channel_waveform*param.spikes.amplitude_factor;
            spikes(:, channel, rec) = current_channel_waveform;
        end
    end
    
    Timestamps_accepted_spikes_TT{TT} = Timestamps_accepted_spikes;
    spikes_TT{TT} = spikes;

end %end looping over Tetrodes


%% Coincidence detection - Tamir's version
do_coincidence_detection = 1;
if do_coincidence_detection && (length(params.TT_to_use) > 1)

    disp('-------------------------------------------')
    disp('Coincidence-Detection (eliminate artifacts)')
    tic
    
    %% find the invalid ts (with CD)
    TTs_events = zeros(length(Timestamps_accepted_spikes_TT), length(timestamps));
    for ii_TT = 1:length(Timestamps_accepted_spikes_TT)
        spikes_ts_IX = ismember(timestamps,Timestamps_accepted_spikes_TT{ii_TT});
        TTs_events(ii_TT,spikes_ts_IX) = 1;
    end
    mask = conv(sum(TTs_events), ones(1,params.CD_detect_win_len), 'same');
    mask = mask >= params.CD_n_TT_thr;
    mask = conv(mask, ones(1,params.CD_invalid_win_len), 'same');
    TTs_events(:,mask>0) = 0;
    
    % Last, update accepted spikes (ts+shape)
    for ii_TT = 1:length(params.TT_to_use)
        TT = params.TT_to_use(ii_TT);
        % find the invalid spikes in each TT
        old_spikes_ts = Timestamps_accepted_spikes_TT{TT};
        new_spikes_ts = timestamps( TTs_events(TT,:)==1 );
        invalid_spikes = ~ismember(old_spikes_ts, new_spikes_ts);
        Timestamps_accepted_spikes_TT{TT}(invalid_spikes) = [];
        spikes_TT{TT}(:,:,invalid_spikes) = [];
    end
    
    toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the data in Neuralynx NTT files for three cases:
% (1) All the data.

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


%% close log file
diary off





