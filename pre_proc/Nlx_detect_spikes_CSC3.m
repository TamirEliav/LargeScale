function Nlx_detect_spikes_CSC3(dir_IN,dir_OUT,params)

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
% params.is_save_artifacts = 0;

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
    fprintf('\t Loading raw data \n')
    tic
    act_ch = find(params.active_TT_channels(TT, :));
    non_act_ch = find(~params.active_TT_channels(TT, :));
    csc = [];
    timestamps = [];
    for ch = act_ch
        filename = TT_files_raw{TT,ch};
        [csc(ch,:), timestamps, fs] = Nlx_csc_read(fullfile(dir_IN,filename),params.t_start_end);
        if ~isempty(params.ref_ch) && params.ref_ch(1) ~= TT
            csc(ch,:) = csc(ch,:) - ref_ch_csc;
        end
    end
    % place zeros for disconnected channels
    for ch = non_act_ch
        csc(ch,:) = zeros(size(timestamps));
    end
    time_measure.load_TT(TT) = toc;

    %% detect threshold crossing events
    fprintf('\t Detect spike (thr crossing) \n')
    tic
    TT_ch_events_IX = {};
    for ch = act_ch
        thr = params.thr_uV(TT,ch);
        use_neg_thr = params.use_neg_thr(TT,ch);
        % report stats (choosen thr vs. median/std)
        stats.TT_ch_std(TT,ch) = std(csc(ch,:));
        stats.TT_ch_abs_median(TT,ch) = median(abs(csc(ch,:)));
        stats.thr_div_abs_median(ch,TT) = thr / stats.TT_ch_abs_median(TT,ch);
        stats.thr_div_std(ch,TT) = thr / stats.TT_ch_std(TT,ch);
        
        events_pos_pks_IX = [];
        events_neg_pks_IX = [];
        % detect positive events
        tic
        [~,events_pos_pks_IX] = findpeaks(csc(ch,:),'MinPeakHeight',thr);
        toc
        % detect negative events
        tic
        if params.use_neg_thr
            [~,events_neg_pks_IX] = findpeaks(-csc(ch,:),'MinPeakHeight',thr);
        end
        toc
        tic
        events_pks_IX = sort([events_pos_pks_IX events_neg_pks_IX], 'ascend');
        toc
        tic
        % remove spikes from beginning or end of data
        events_pks_IX(events_pks_IX < params.AlignSample) = [];
        events_pks_IX(events_pks_IX + params.nSamples-params.AlignSample > length(timestamps)) = [];
        toc
        TT_ch_events_IX{ch} = events_pks_IX;
    end
    time_measure.thr_crs(TT) = toc;
	stats.detect_nEvents(:,TT) = cellfun(@length,TT_ch_events_IX);
    
    %% merge event from mall channels
    fprintf('\t Merge spikes from different channels \n')
    tic
    M = zeros(4,length(timestamps));
    for ch = act_ch
        M(ch,TT_ch_events_IX{ch}) = csc(ch,TT_ch_events_IX{ch});
    end
    M = max(abs(M));
    [~,TT_events_IX] = findpeaks(M,'MinPeakHeight',min(params.thr_uV(TT,:)));
    clear M
    time_measure.merge_ch(TT) = toc;
    stats.merge_nEvents(TT) = length(TT_events_IX);
    
    %% extract waveforms
    fprintf('\t Extract waveforms \n')
    tic
    wvfrms = zeros(params.nSamples,4,length(TT_events_IX)); % Initialize
    % IX relative to spikes peak
    trigger_IX = repmat( [(1-params.AlignSample):(params.nSamples-params.AlignSample)]',1,length(TT_events_IX));
    % IX relative to csc signal
    trigger_IX = trigger_IX + repmat(TT_events_IX,[params.nSamples 1]);
    for ch = 1:4
        temp = csc(ch,:); % TODO: how we can add another dimension and avoid this stupid loop and temp...?!
        wvfrms(:,ch,:) = temp(trigger_IX);
    end
    
    time_measure.extract_wvfrm(TT) = toc;
    stats.extract_nWvfrm(:,TT) = size(wvfrms);
    
    %% library of acceptable spike shapes
    % Clean artifacts using the library of acceptable spike shapes.
    % Now we will throw away noisy theshold crossing using the library of
    % acceptable spike shapes as reference

    fprintf('\t library of acceptable spike shapes \n')
    tic
    load(params.lib_spike_shapes);
    
    % For each event take the channel with the largest peak (8th point)
    [~,max_ch_IX] = max(squeeze(abs(wvfrms(params.AlignSample,:,:))),[],1);
    largest_waveforms = zeros(size(wvfrms,1),size(wvfrms,3));
    for ch=1:4
        IX = find(max_ch_IX == ch);
        largest_waveforms(:,IX) = squeeze(wvfrms(:,ch,IX));
    end
    
    % calc corr
    xxx_lags_shifts = [1:30; 2:31; 3:32];
    ccc = [];
    rrr = [];
    for ii_shift = 1:size(xxx_lags_shifts,1)
        xxx_lags = xxx_lags_shifts(ii_shift, :);
        ccc = corr(largest_waveforms(xxx_lags,:), library_of_acceptable_spike_shapes(:,2:end-1)');
        rrr(ii_shift,:) = max(ccc,[],2);
    end
    rrr = max(rrr,[],1);
    TT_events_lib_thrded = ( rrr >=  params.lib_corr_thr );
    
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
    stats.lib_nAccepted(TT) = sum(TT_events_lib_thrded);
    stats.lib_nNotAccepted(TT) = sum(~TT_events_lib_thrded);
    
    % create accepted spikes list
    TT_spikes_IX = TT_events_IX(TT_events_lib_thrded);
    TT_invalid_lib_IX = TT_events_IX(~TT_events_lib_thrded);
    
    %% detect possible artifact (by high voltage)
    for ch = act_ch
        artifact_thr = 6 .* stats.TT_ch_std(TT,ch);
        TT_high_amp_art_IX{ch} =  find(abs(csc(ch,:)) > artifact_thr);
    end
    
    %% TODO: consider adding minimal window seperation
    
    %% arrange data from current TT
    all_TT.TT_events_IX{TT} = TT_events_IX;
    all_TT.wvfrms{TT} = wvfrms;
    all_TT.TT_events_lib_thrded{TT} = TT_events_lib_thrded;
    all_TT.TT_spikes_IX{TT} = TT_spikes_IX;
    all_TT.TT_invalid_lib_IX{TT} = TT_invalid_lib_IX;
    all_TT.TT_high_amp_art_IX{TT} = TT_high_amp_art_IX;
    
end

%% Coincidence detection - Tamir's version
do_coincidence_detection = 1;
tic
if do_coincidence_detection && (length(params.TT_to_use) > 1)

    disp('-------------------------------------------')
    disp('Coincidence-Detection (eliminate artifacts)')
    
    %% find the invalid ts (with CD)
    % first, mark 1 for event for each TT
    high_amp_art_IX = cat(2,all_TT.TT_high_amp_art_IX{:});
    mask = zeros(size(timestamps));
    for ii = 1:length(high_amp_art_IX)
        mask_temp = zeros(size(timestamps));
        mask_temp(high_amp_art_IX{ii}) = 1;
        mask = mask+mask_temp;
    end
% %     mask = zeros(length(params.TT_to_use), length(timestamps));
% %     % TODO: replace with summing up the ones in the for loop (to avoid
% %     % large memoryu allocation for 16 TT case)
% %     for ii_TT = 1:length(params.TT_to_use)
% %         TT = params.TT_to_use(ii_TT);
% % %         mask(ii_TT, all_TT.TT_spikes_IX{TT}) = 1;
% %         mask(ii_TT, all_TT.TT_events_IX{TT}) = 1; % NOTE: we are using all the detected events!! (not only the expected spikes)
% %     end
% %     mask = sum(mask);
    mask = conv(mask, ones(1,params.CD_detect_win_len), 'same'); % sum across TT and convolve with CD detection window
%     mask = mask >= params.CD_n_TT_thr; % apply min number of TT for detection
    mask = mask >= params.CD_n_ch_thr; % apply min number of ch for detection
    mask = conv(mask, ones(1,params.CD_invalid_win_len), 'same'); % convolve dections with the invalidation window (we are after detection!)
    mask = mask>0; % now in mask we have 1 for CD invalidation
    
    % Last, update accepted spikes (ts+shape)
    for ii_TT = 1:length(params.TT_to_use)
        TT = params.TT_to_use(ii_TT);
        % find the valid/invalid spikes in each TT
        old_TT_spikes_IX = all_TT.TT_spikes_IX{TT};
        spikes_CD_valid_IX = old_TT_spikes_IX(find(~mask(old_TT_spikes_IX)));
        spikes_CD_invalid_IX = old_TT_spikes_IX(find(mask(old_TT_spikes_IX)));
        
        % update the global variable
        all_TT.TT_spikes_IX{TT} = spikes_CD_valid_IX;
        all_TT.spikes_CD_invalid_IX{TT} = spikes_CD_invalid_IX;
        
        stats.CD_nEventRemove(TT) = length(spikes_CD_valid_IX);
    end
    
end
time_measure.CD = toc;

%%
stats.rec_len.nSamples = length(timestamps);
stats.rec_len.totalTimeHours = range(timestamps)*1e-6/60/60;

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save the data in Neuralynx NTT files
% (1) valid spikes
% (2) invalid events by library
% (3) invalid events by CD
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
    
    % 1. only valid spikes
    filename_out = fullfile(dir_OUT,NTT_filename);
    [~,events_IX] = ismember(all_TT.TT_spikes_IX{TT}, all_TT.TT_events_IX{TT});
    Timestamps = timestamps(all_TT.TT_spikes_IX{TT});
    Samples = all_TT.wvfrms{TT}(:,:,events_IX);
    write_NTT_file(filename_out, Timestamps, Samples, [], fs)
    
    if params.is_save_artifacts
        % 2. all detected events, marked by units: 
        %   0 valid spikes
        %   1 lib removed
        %   2 CD removed
        filename_out = fullfile(dir_OUT,strrep(NTT_filename,'.NTT','_with_artifacts.NTT'));
        Timestamps = timestamps(all_TT.TT_events_IX{TT});
        Samples = all_TT.wvfrms{TT};
        [~,valid_spikes_events_IX] = ismember(all_TT.TT_spikes_IX{TT},          all_TT.TT_events_IX{TT});
        [~,lib_removed_events_IX]  = ismember(all_TT.TT_invalid_lib_IX{TT},     all_TT.TT_events_IX{TT});
        [~,CD_removed_events_IX]   = ismember(all_TT.spikes_CD_invalid_IX{TT},  all_TT.TT_events_IX{TT});
        CellNumbers = zeros(size(Timestamps));
        CellNumbers(valid_spikes_events_IX) = 0;
        CellNumbers(lib_removed_events_IX) = 1;
        CellNumbers(CD_removed_events_IX) = 2;
        write_NTT_file(filename_out, Timestamps, Samples, CellNumbers, fs);
    end
    
    
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

end



%%
function write_NTT_file(filename_out, Timestamps, Samples, CellNumbers, fs)
    
header_file = 'Nlx_header_NTT.txt';
header = textread(header_file, '%s', 'delimiter', '\n', 'whitespace', '');

if isempty(CellNumbers)
    CellNumbers = zeros(size(Timestamps));
end
%     csc header
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

end


%%








