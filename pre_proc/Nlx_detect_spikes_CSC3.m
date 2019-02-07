function Nlx_detect_spikes_CSC3(dir_IN, dir_OUT, params, forcecalc)

%% 
% Tamir,
% 01/2018
% Adopted from Didi, Arseny and Michael 
% 
% Here we detect the spikes for each tetrode and save them.
% Output is 2 NTT files containing:
% 1. detected spikes
% 2. detected spikes + events that were detected but removed along the process

%% input validation
if nargin<4; forcecalc = 0; end
if ~exist(dir_IN,'dir')
    error('input folder does not exist');
end
if exist(dir_OUT,'dir')
    if forcecalc
        % delete existing output dir
        warning('spikes detection output dir already existing and you chose to override it, deleting old spikes detection dir!');
        rmdir(dir_OUT,'s')
    else
        error('spikes detection output folder already exist, use forcecalc to override it!');
    end
end
% at this point we should not have the output dir, so let's create it!
mkdir(dir_OUT)

%% open log file
log_name_str = ['spikes_detection_' datestr(clock, 'yyyy-mm-dd HH-MM-SS') '.txt'];
log_name_out = fullfile(dir_OUT, log_name_str );
if ~exist(dir_OUT,'dir'), mkdir(dir_OUT); end 
diary off; diary(log_name_out); diary on
start_runtime_tic = tic;

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
% "inflate" thr for all channels
if size(params.thr,1) == 1
    % same for all TTs!
    params.thr = repmat(params.thr, size(TT_ch_exist));
elseif size(params.thr,2) == 1
    % per TT (same for all channels in that TT...)
    params.thr = repmat(params.thr, [1 size(TT_ch_exist,2)]);
elseif any(size(params.thr)~=size(TT_ch_exist))
    error('wrong input dimensions!');
end
params.thr_uV     = nan(size(TT_ch_exist)); % will be filled later 
params.thr_median = nan(size(TT_ch_exist)); % will be filled later 

if size(params.use_neg_thr,1) == 1
    % same for all TTs!
    params.use_neg_thr = repmat(params.use_neg_thr, size(TT_ch_exist));
elseif size(params.use_neg_thr,2) == 1
    % per TT (same for all channels in that TT...)
    params.use_neg_thr = repmat(params.use_neg_thr, [1 size(TT_ch_exist,2)]);
elseif any(size(params.use_neg_thr)~=size(TT_ch_exist))
    error('wrong input dimensions!');
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
    
    %% load raw data (CSCs)
    fprintf('\t Loading raw data \n')
    tic
    act_ch = find(params.active_TT_channels(TT, :));
    non_act_ch = find(~params.active_TT_channels(TT, :));
    csc = []; % TODO: consider pre-allocating this large array! (need to dummy read one channel to know the size)
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
        % calc channel stats
        stats.TT_ch_std(TT,ch) = std(csc(ch,:));
        stats.TT_ch_abs_median(TT,ch) = median(abs(csc(ch,:)));
        
        switch params.thr_type 
            case 'uVolt'
                params.thr_uV(TT,ch) = params.thr(TT,ch);
                params.thr_median(TT,ch) = params.thr(TT,ch) / stats.TT_ch_abs_median(TT,ch);
            case 'median'
                params.thr_median(TT,ch) = params.thr(TT,ch);
                params.thr_uV(TT,ch) = params.thr_median(TT,ch) .* stats.TT_ch_abs_median(TT,ch);
            otherwise
                error('unsupported thr type!');
        end
        thr = params.thr_uV(TT,ch);
        
        events_pos_pks_IX = [];
        events_neg_pks_IX = [];
        % detect positive events
        [~,events_pos_pks_IX] = findpeaks(csc(ch,:),'MinPeakHeight',thr);
        % detect negative events
        if params.use_neg_thr
            [~,events_neg_pks_IX] = findpeaks(-csc(ch,:),'MinPeakHeight',thr);
        end
        events_pks_IX = sort([events_pos_pks_IX events_neg_pks_IX], 'ascend');
        % remove spikes from beginning or end of data
        events_pks_IX(events_pks_IX < params.AlignSample) = [];
        events_pks_IX(events_pks_IX + params.nSamples-params.AlignSample > length(timestamps)) = [];
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
    
    %% remove duplicates (within window)
    M = zeros(4,length(timestamps));
    for ch = act_ch
        M(ch,TT_spikes_IX) = csc(ch,TT_spikes_IX);
    end
    M = max(abs(M));
    
    [~,new_spikes_IX] = findpeaks(M, 'MinPeakDistance',params.min_sep_events, 'MinPeakHeight',min(params.thr_uV(TT,:)));
    TT_min_sep_win_removed = setdiff(TT_spikes_IX,new_spikes_IX);
    TT_spikes_IX = new_spikes_IX;
    clear M;
    
    %% detect possible artifact (by high voltage) - later for coincidence detection
    for ch = act_ch
%         artifact_thr = 6 .* stats.TT_ch_std(TT,ch);
        artifact_thr = params.thr_uV(TT,ch);
        TT_high_amp_art_IX{ch} =  find(abs(csc(ch,:)) > artifact_thr);
    end
    
    %% TODO: consider adding minimal window seperation
    
    %% arrange data from current TT
    all_TT.TT_events_IX{TT} = TT_events_IX;
    all_TT.wvfrms{TT} = wvfrms;
    all_TT.TT_events_lib_thrded{TT} = TT_events_lib_thrded;
    all_TT.TT_invalid_lib_IX{TT} = TT_invalid_lib_IX;
    all_TT.TT_min_sep_win{TT} = TT_min_sep_win_removed;
    all_TT.TT_high_amp_art_IX{TT} = TT_high_amp_art_IX;
    all_TT.TT_spikes_IX{TT} = TT_spikes_IX;
    
end

%% Coincidence detection - Tamir's version
do_coincidence_detection = 1;
tic
if do_coincidence_detection && (length(params.TT_to_use) > 1)

    disp('-------------------------------------------')
    disp('Coincidence-Detection (eliminate artifacts)')
    
    %% find the invalid ts (with CD)
    % first, mark 1 for event for each TT
    
    % option 1 - using all detected event pooled for each TT from all channels -
    % uVolt thr different from spikes
    mask = zeros(size(timestamps));
    for ii_TT = 1:length(all_TT.TT_high_amp_art_IX)
        TT_high_amp_art_IX = all_TT.TT_high_amp_art_IX{ii_TT};
        high_amp_art_IX = cat(2, TT_high_amp_art_IX{:});
        mask_temp = zeros(size(timestamps));
        mask_temp(high_amp_art_IX) = 1;
        mask_temp = conv(mask_temp, ones(1,params.CD_detect_win_len), 'same'); % convolve each source with some window
        mask_temp(mask_temp>1) = 1;
        mask = mask+mask_temp;
    end
    nSourceEvent_thr = params.CD_n_TT_thr;
    
% %     % option 2 - using all detected event from all channels - uVolt thr
% %     % seperate from spikes
% %     high_amp_art_IX = cat(2,all_TT.TT_high_amp_art_IX{:});
% %     mask = zeros(size(timestamps));
% %     for ii = 1:length(high_amp_art_IX)
% %         mask_temp = zeros(size(timestamps));
% %         mask_temp(high_amp_art_IX{ii}) = 1;
% %         mask = mask+mask_temp;
% %     end
% %     nSourceEvent_thr = params.CD_n_ch_thr;
    
% %     % option 3 - using detected events per TT (all events/only lib
% %     % detected) - uVolt thr as for spikes  (original code)
% %     mask = zeros(length(params.TT_to_use), length(timestamps));
% %     % TODO: replace with summing up the ones in the for loop (to avoid
% %     % large memoryu allocation for 16 TT case)
% %     for ii_TT = 1:length(params.TT_to_use)
% %         TT = params.TT_to_use(ii_TT);
% % %         mask(ii_TT, all_TT.TT_spikes_IX{TT}) = 1;
% %         mask(ii_TT, all_TT.TT_events_IX{TT}) = 1; % NOTE: we are using all the detected events!! (not only the expected spikes)
% %     end
% %     mask = sum(mask);
% %     nSourceEvent_thr = params.CD_n_TT_thr;

    mask = mask >= nSourceEvent_thr ; % apply min number of sources for detection (sources = TT or ch)
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
    
    % 2. all detected events, marked by units: 
    %   0 valid spikes
    %   1 lib removed
    %   2 duplicates removed (min_sep_win)
    %   3 CD removed
    Timestamps = timestamps(all_TT.TT_events_IX{TT});
    Samples = all_TT.wvfrms{TT};
    [~,valid_spikes_events_IX] = ismember(all_TT.TT_spikes_IX{TT},          all_TT.TT_events_IX{TT});
    [~,lib_removed_events_IX]  = ismember(all_TT.TT_invalid_lib_IX{TT},     all_TT.TT_events_IX{TT});
    [~,min_sep_win_events_IX]  = ismember(all_TT.TT_min_sep_win{TT},        all_TT.TT_events_IX{TT});
    [~,CD_removed_events_IX]   = ismember(all_TT.spikes_CD_invalid_IX{TT},  all_TT.TT_events_IX{TT});
    CellNumbers = zeros(size(Timestamps));
    CellNumbers(valid_spikes_events_IX) = 0;
    CellNumbers(lib_removed_events_IX) = 1;
    CellNumbers(min_sep_win_events_IX) = 2;
    CellNumbers(CD_removed_events_IX) = 3;
    
    %% Generate figures for the detected spikes+artifacts 
    % (even if user didn't want to save those in NTT file)
    unit_colors = {0.7.*[1 1 1],'r','g','b'};
    unit_labels = {'detected','lib','min sep','Coincidence'};
    unit_numbers = [0 1 2 3];
    features_pairs = [1 2;1 3;1 4;2 3;2 4;3 4];
    max_points_plot_wvfrm = 2000;
    max_points_plot_cluster = 50000;
    max_volt_plot = 500;
    %% waveforms
    figure('Units','normalized','Position',[0 0 1 1]);
    pnl = panel();
    pnl.pack(4,4);
    pnl.margin = [30 30 20 20];
    pnl.de.margin = 10;
    h=pnl.title(NTT_filename);
    h.Interpreter = 'none';
    h.Position = [0.5 1.03];
    h.FontSize = 14;
    for ii_unit = 1:length(unit_numbers)
        unit = unit_numbers(ii_unit);
        h=pnl(ii_unit).ylabel(unit_labels{ii_unit});
        h.Position = [-0.03 0.5];
        h.FontSize = 13;
        for ch=1:4
            pnl(ii_unit,ch).select();
            IX = find(CellNumbers==unit);
            rng(0);
            n = length(IX);
            k = min(max_points_plot_wvfrm,n);
            subset_IX = randsample(n,k);
            plot(squeeze(Samples(:,ch,IX(subset_IX))), 'Color',unit_colors{ii_unit});
            xlimits = get(gca,'xlim');
            ylimits = get(gca,'ylim');
            xlimits = min(abs(xlimits),max_volt_plot).*sign(xlimits);
            ylimits = min(abs(ylimits),max_volt_plot).*sign(ylimits);
            set(gca,'xlim',xlimits);
            set(gca,'ylim',ylimits);
            text(1,1, sprintf('n=%d',length(IX)), 'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top');
            if ii_unit == length(unit_numbers)
                xlabel(sprintf('ch%d',ch))
            end
        end
    end
    pnl(length(unit_numbers),1).select();
    xlabel('Samples')
    ylabel('{\mu}Volt')
    linkaxes( pnl.de.axis , 'x');
    pnl(1,4).select();
    text(1.1,1.1,sprintf('max waveforms=%d',max_points_plot_wvfrm),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top');
    filename_out = fullfile(dir_OUT,strrep(NTT_filename,'.NTT','_artifacts_waveforms'));
    saveas(gcf, filename_out, 'tif')
    close(gcf)
    %% clusters
    figure('Units','normalized','Position',[0 0 1 1]);
    pnl = panel();
    pnl.pack(4,6);
    pnl.margin = [30 30 20 20];
    pnl.de.margin = 15;
    h=pnl.title(NTT_filename);
    h.Interpreter = 'none';
    h.Position = [0.5 1.03];
    h.FontSize = 14;
    for ii_unit = 1:length(unit_numbers)
        unit = unit_numbers(ii_unit);
        h=pnl(ii_unit).ylabel(unit_labels{ii_unit});
        h.Position = [-0.03 0.5];
        h.FontSize = 13;
        for ii_pair = 1:size(features_pairs,1)
            pnl(ii_unit,ii_pair).select(); hold on;
            IX = find(CellNumbers==unit);
            rng(0);
            n = length(IX);
            k = min(max_points_plot_cluster,n);
            subset_IX = randsample(n,k);
            x = squeeze(Samples(params.AlignSample, features_pairs(ii_pair,1), IX(subset_IX)));
            y = squeeze(Samples(params.AlignSample, features_pairs(ii_pair,2), IX(subset_IX)));
            plot(x,y, '.', 'Color',unit_colors{ii_unit});
            xlimits = get(gca,'xlim');
            ylimits = get(gca,'ylim');
            xlimits = min(abs(xlimits),max_volt_plot).*sign(xlimits);
            ylimits = min(abs(ylimits),max_volt_plot).*sign(ylimits);
            set(gca,'xlim',xlimits);
            set(gca,'ylim',ylimits);
            text(1,1, sprintf('n=%d',length(IX)), 'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top');
            xlabel(sprintf('ch%d (%s)', features_pairs(ii_pair,1), '{\mu}V'))
            ylabel(sprintf('ch%d (%s)', features_pairs(ii_pair,2), '{\mu}V'))
        end
    end
    pnl(1,size(features_pairs,1)).select();
    text(1.1,1.1,sprintf('max points=%d',max_points_plot_cluster),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top');
    filename_out = fullfile(dir_OUT,strrep(NTT_filename,'.NTT','_artifacts_clusters'));
    saveas(gcf, filename_out, 'tif')
    close(gcf)
    
    %% save artifacts in NTT file
    if params.is_save_artifacts
        filename_out = fullfile(dir_OUT,strrep(NTT_filename,'.NTT','_with_artifacts.NTT'));
        write_NTT_file(filename_out, Timestamps, Samples, CellNumbers, fs);
    end
    
    
end % end looping over tetrodes
time_measure.write_NTTs = toc;

%% total run-time
time_measure.total_runtime = toc(start_runtime_tic);

%% report parmas/stats/run-time
fn_structdisp(stats)
fn_structdisp(time_measure)
fn_structdisp(params)

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
InputRange = max(abs(Samples(:)));
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








