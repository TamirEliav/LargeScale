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

% % % %%
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
params.merge_thr_crs_width = 4;
params.use_neg_thr = 0;
params.lib_spike_shapes = 'library_of_acceptable_spike_shapes.mat';
params.lib_corr_thr = 0.8;
param.coincidence_window = 1e3;

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
    params.TT_to_use( params.TT_to_use == params.ref_ch(1) ) = [];
end


%% detect spikes
for TT = params.TT_to_use

    disp(['running detection on TT' TT])
    
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
        if ~isempty(params.ref_ch)
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


% % %     % take as a start point all the spikes from ch1
% % %     SPK_All_Ch_combined = SPK_IX{1};
% % % 
% % %     % merge into this list spikes from all other channels
% % %     for curr_chan_idx=2:4
% % %         % First loop over the spikes detected on the first channel and compare their seperation from those detected on the second channel
% % %         % Save only those which are seperated by a minimal number of bins to avoid counting the same spike twice
% % %         tic
% % %         for ii_spike = 1:length(SPK_All_Ch_combined)
% % %             current_spike_IX = SPK_All_Ch_combined(ii_spike);
% % %             shared_IXs = find(abs(SPK_IX{curr_chan_idx} - current_spike_IX)<= param.spikes.x_sep_spike_thres);
% % %             if ~isempty(shared_IXs) % i.e., the same spike is detected twice
% % %                 SPK_IX{curr_chan_idx}(shared_IXs) = [];
% % %             end
% % %         end
% % %         toc
% % %         % Now that we do not have the same spikes on two channels we can merge the spike IXs of both channels,
% % %         % as those represent different spikes
% % %         SPK_All_Ch_combined = [SPK_All_Ch_combined,SPK_IX{curr_chan_idx}];
% % %     end
% % %     
% % %     % sort them to preserve their temporal order:
% % %     SPK_IXs_All_Ch_combined_sorted = sort(SPK_All_Ch_combined);
    
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
        ccc = corr(spikes_waveforms(:,xxx_lags)', library_of_acceptable_spike_shapes(:,xxx_lags)');
        rrr(ii_shift,:) = max(ccc,[],2);
    end
    rrr = max(rrr,[],1);
    
    vector_of_accepted_spikes = ( rrr >=  params.lib_corr_thr ); % TODO: plot hist of 'r' values + thr
    
% % % %     for ii_spike = 1:size(SPK_waveforms,2)
% % % %         %(ii_spike/size(SPK_waveforms,2))*100
% % % %         spike_shape_4channels = zeros(32,4);
% % % %         for jj_channel = 1:size(SPK_waveforms,1)
% % % %             spike_shape_4channels(:,jj_channel) = SPK_waveforms{jj_channel,ii_spike}';
% % % %         end
% % % %         % Choose the channel # for which the spike has the largest height:
% % % %         [ stam  idx_channel_max_height ] = max( max( spike_shape_4channels ) );
% % % %         spike_shape = spike_shape_4channels( :, idx_channel_max_height )' ;
% % % %         
% % % %         if ( std( spike_shape(2:end-1) ) == 0 ), % If this is a completely FLAT "spike", I cannot compute CORRCOEF, so I will set r = 0
% % % %             vector_of_max_r_values( ii_spike ) = 0 ;  % Set r = 0 in this case
% % % %         else % If this spike DOES have some shape (this is the case basically for ALL the recorded waveforms)
% % % %             % Compute the correlation coefficients with all the acceptable spike shapes -- lag 0 :
% % % %             xxx_lags = 2 : 31 ; % Use CENTRAL 30 points
% % % %             ccc = corrcoef([ spike_shape(2:end-1) ; library_of_acceptable_spike_shapes(:,xxx_lags)]' ); % Correlation coefficients
% % % %             rrr_vec_lag_0 = ccc( 1, 2:end ); % All the correlation coefficients with the acceptable spike shapes
% % % %             % Compute the correlation coefficients with all the acceptable spike shapes -- lag (+1) :
% % % %             xxx_lags = 1 : 30 ; % Use FIRST 30 points (RIGHT shift of the "Acceptable Spike Shapes matrix")
% % % %             ccc = corrcoef([ spike_shape(2:end-1) ; library_of_acceptable_spike_shapes(:,xxx_lags)]' ); % Correlation coefficients
% % % %             rrr_vec_lag_plus1 = ccc( 1, 2:end ); % All the correlation coefficients with the acceptable spike shapes
% % % %             % Compute the correlation coefficients with all the acceptable spike shapes -- lag (-1) :
% % % %             xxx_lags = 3 : 32 ; % Use LAST 30 points (LEFT shift of the "Acceptable Spike Shapes matrix")
% % % %             ccc = corrcoef([ spike_shape(2:end-1) ; library_of_acceptable_spike_shapes(:,xxx_lags)]' ); % Correlation coefficients
% % % %             rrr_vec_lag_minus1 = ccc( 1, 2:end ); % All the correlation coefficients with the acceptable spike shapes
% % % %             % Save the MAXIMAL r value -- use the maximal correlation among the 3 lags (-1,0,+1):
% % % %             vector_of_max_r_values( ii_spike ) = max( [ rrr_vec_lag_0  rrr_vec_lag_plus1  rrr_vec_lag_minus1 ] );
% % % %         end
% % % %         
% % % %         % Determine if this spike should be Accepted (set value to '1') or Rejected (set value to '0'):
% % % %         vector_of_accepted_spikes( ii_spike ) = ...
% % % %             vector_of_max_r_values( ii_spike )  >=  r_threshold ;
% % % %         % Accept the spike shape ('1') if its correlation with ANY of the acceptable shapes
% % % %         % is >= r_threshold ; else, reject the spike ('0').
% % % %         % FOR DBG: plot waveform in 'r' if rejected and 'b' if accepted
% % % %         %     if (vector_of_sccepted_spikes( ii_spike )== 1)
% % % %         %     plot(spike_shape)
% % % %         %     else
% % % %         %     plot(spike_shape,'r')
% % % %         %     end
% % % %         %     title(['spike num = ',num2str(ii_spike), ' ,r val = ' num2str(vector_of_max_r_values( ii_spike ))])
% % % %         %     waitforbuttonpress
% % % %     end % End "Loop over spikes extracted from the Ntt file"
    

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
    clear timestamps 

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
    
    % first, build the indeces of overlapping timestamps between all TT pairs
    CD1_IX = {};
    for ii_TT1 = 1:length(params.TT_to_use)
        for ii_TT2 = 1:length(params.TT_to_use)
            if ii_TT1 == ii_TT2 
                continue;
            end
            ts_TT1 = Timestamps_accepted_spikes_TT{ params.TT_to_use(ii_TT1) };
            ts_TT2 = Timestamps_accepted_spikes_TT{ params.TT_to_use(ii_TT2) };
            IX = find_vec_in_vec( ts_TT1, ts_TT2);
            IX_TT1 = find( abs(ts_TT1 - ts_TT2(IX)) <= param.coincidence_window );
            IX_TT2 = IX(IX_TT1);
            CD1_IX{ii_TT1,ii_TT2,1} = IX_TT1;
            CD1_IX{ii_TT1,ii_TT2,2} = IX_TT2;
        end
    end
    
    % second, check for 1st order coincidence
    CD2_IX = {};
    for ii_TT1 = 1:size(CD1_IX,1)
        
        other_TTs_IX = 1:size(CD1_IX,1);
        other_TTs_IX(other_TTs_IX==ii_TT1) = []; % ix of other tts
        sdf = CD1_IX(ii_TT1,other_TTs_IX,1);
        IX_TT1 = mintersect( sdf{:} ); % indeces of current tt that coincide with ALL other tts event
        for ii_TT2 = other_TTs_IX
            temp_IX = ismember( CD1_IX{ii_TT1,ii_TT2,1}, IX_TT1 );
            CD2_IX{ii_TT1, ii_TT2, 1} = CD1_IX{ii_TT1,ii_TT2,1}(temp_IX);
            CD2_IX{ii_TT1, ii_TT2, 2} = CD1_IX{ii_TT1,ii_TT2,2}(temp_IX);
        end
    end
    
    % third, check for 2nd order coincidence 
    CD3_IX = {};
    for ii_TT1 = 1:size(CD1_IX,1)
        
        other_TTs_IX = 1:size(CD1_IX,1);
        other_TTs_IX(other_TTs_IX==ii_TT1) = []; % ix of other tts
        sdf = CD1_IX(other_TTs_IX,ii_TT1,2); % note the swapping
        CD3_IX{ii_TT1} = mintersect( sdf{:} );
    end
    
    % Last, update accepted spikes (ts+shape)
    for ii_TT = 1:length(params.TT_to_use)
        TT = params.TT_to_use(ii_TT);
        Timestamps_accepted_spikes_TT{TT}(CD3_IX{ii_TT}) = [];
        spikes_TT{TT}(:,:,CD3_IX{ii_TT}) = [];
    end
    
    toc
end
% %         time_diff = abs(bsxfun(@minus, ts_A, ts_B'));


%%  Read Timestamp data and use Coincidence-Detection across Tetrodes to eliminate artifacts (just as we do in the wired case): --------
% -------------------------------------------- OLD version ----------------
do_coincidence_detection = 0;
if do_coincidence_detection
    
    disp('-------------------------------------------')
    disp('Coincidence-Detection (eliminate artifacts)')
    tic
    
    idx_coincidence_vec{length(param.tetrodes.use_tetrodes)} = [];  % Initialize this variable (for later)
    
    % Find coincidence-detection events = coincidence-detection on a millisecond-scale:
    for ii_spikes = 1:length(Timestamps_accepted_spikes_TT{1}), % Loop over the spikes of the FIRST tetrode
        t_spike = Timestamps_accepted_spikes_TT{1}(ii_spikes);
        idx{1} = ii_spikes; % Index of this spike
        
        for ii_file = 2:length(param.tetrodes.use_tetrodes) % Loop over the other tetrodes, finding the coincidnce-detection
            idx{ii_file} = find( abs( Timestamps_accepted_spikes_TT{ii_file} - t_spike ) <= param.spikes.coincidence_window ); % THE COINCIDENCE DETECTION
        end
        
        % Check that the coincidence-detection occurred on ALL tetrodes:
        test_variable = 1;
        for ii_file =1: length(param.tetrodes.use_tetrodes) % Loop over all tetrode files
            if ( isempty( idx{ii_file} ) ), % If there is NO coincidence-detection
                test_variable = 0 ;
            end
        end
        
        if length(param.tetrodes.use_tetrodes)>=3 %If there are at least 3 tetrodes
            % Check that the temporal separation between the OTHER two tetrodes also meets the coincidence-detection criterion
            % (this is needed since the spikes on the OTHER two tetrodes may be up to twice-the-time apart, in principle!!! ):
            if ( test_variable == 1 ), % If we passed the first test
                if ( abs( Timestamps_accepted_spikes_TT{2}(idx{2}(1)) - Timestamps_accepted_spikes_TT{3}(idx{3}(1)) ) > param.spikes.coincidence_window  ),
                    test_variable = 0 ; % Reset test_variable if the time between spikes on the OTHER two tetrodes is too large
                end
            end
        end;
        
        % Save the indexes of the coincidence-detection "spikes" to be REMOVED:
        if ( test_variable == 1 ), % If all indexes exist = there IS a coincidence detection on ALL tetrodes
            for ii_file =1:length( param.tetrodes.use_tetrodes) % Loop over all tetrode files
                idx_coincidence_vec{ii_file} = [ idx_coincidence_vec{ii_file}  idx{ii_file} ];
            end
        end
%         if mod(ii_spikes,1000) == 0
%             disp(['Processing coindicedence detection across TT- ' num2str((ii_spikes/length(Timestamps_accepted_spikes_TT{tetrodes(1)}))*100)])
%         end
    end %end looping over spikes for coincidence detection
    
    toc
    
else
    idx_coincidence_vec{1}=[];idx_coincidence_vec{2}=[];idx_coincidence_vec{3}=[];idx_coincidence_vec{4}=[];
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the data in Neuralynx NTT files for three cases:
% (1) All the data.

if ~exist(dir_OUT,'dir'), mkdir(dir_OUT); end 
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
   
    % load generic header to add to NTT files
    load('ntt_generic_header.mat');
    
    Timestamps_accepted_spikes=Timestamps_accepted_spikes_TT{TT};
%     Timestamps_accepted_spikes(idx_coincidence_vec{TT})=[];
    spikes=spikes_TT{TT};
%     spikes(:,:,idx_coincidence_vec{TT})=[];

    %only export timestamps and data points
    FieldSelection = [1 0 0 0 1 0];
    Mat2NlxSpike(filename_out, 0, 1, [], FieldSelection, Timestamps_accepted_spikes, spikes);
    
    
end % end looping over tetrodes

