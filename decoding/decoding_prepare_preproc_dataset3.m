%%
clear
clc

%% get exp info
% main_data_dir = "D:\sequences\DATA";
data_proc_dir = "D:\sequences\seq_uri_eden\proc";
% exp_ID = 'b9861_d180524';
% exp_ID = 'b0034_d180313';
exp_ID = 'b9861_d180526';
% exp_ID = 'b9861_d180527';
exp = exp_load_data(exp_ID);
details = exp.details;

%% get sleep ts
sleep_ti = exp_get_sessions_ti(exp_ID, 'Sleep1','Sleep2');
sleep_ti(any(isnan(sleep_ti),2),:) = []; % remove nan in case of missing sessions

%% prepare NTT file names
NTT_files = {};
tt_to_use = ones(1,exp.details.numTT);
tt_to_use = tt_to_use & exp.details.TT_to_use;
tt_to_use = tt_to_use & contains(exp.details.TT_loc,{'CA1','CA3'});

% data_dir = fullfile(main_data_dir,num2str(details.batNum,'%.4d'),datestr(exp.details.date,'yyyymmdd'),'spikes_detection');
filenames = "spikes_"+exp_ID+"_TT"+find(tt_to_use)'+".ntt";
NTT_files = fullfile(exp.path.decoding_spikes_detection,filenames);

%% load spikes data
clear data;
for ii_file = 1:length(NTT_files)
    %% read NTT file
    NTT_file = NTT_files{ii_file};
    [Timestamps, CellNumbers, Samples, Header] = ...
         Nlx2MatSpike(NTT_file, [1 0 1 0 1], 1, 1, [] );
    %% parse header
    ADBitVolts = sscanf(Header{contains(Header,'ADBitVolts')},'-ADBitVolts %f %f %f %f');
    %% convert bits to uVolts
    Samples = Samples .* ADBitVolts' .* 1e6;
    %% calc features
    features = squeeze(range(Samples,1));
    %% create spikes structure
    spikes.ts = Timestamps;
    spikes.features = features;    
    spikes.NTT_file = NTT_file;
    %%
    NTTs(ii_file) = spikes;
end

%% arrange data
fs=500;
dt=1e6/fs;
tstart = min(exp.details.session_ts(:));
tend = max(exp.details.session_ts(:));
edges = tstart : dt: tend;
t = edges2centers(edges);

% behavioral
position = interp1(exp.pos.proc_1D.ts, exp.pos.proc_1D.pos, t);
direction = interp1(exp.pos.proc_1D.ts, sign(exp.pos.proc_1D.vel_csaps), t);

FE = exp.flight.FE;
FE = FE([FE.distance] > 100);
% FE = FE_remove_low_speed(FE,0.85);
FE_ti = [[FE.start_ts];[FE.end_ts]]';
FE_dist = [FE.distance];

is_FE = any(t>FE_ti(:,1)&t<FE_ti(:,2),1);
is_sleep = any(t>sleep_ti(:,1)&t<sleep_ti(:,2),1);

rest_ti = exp.rest.ti;

% neural
multiunits = nan(length(t),4,length(NTTs));
multiunits_spikes = false(length(t),length(NTTs));
for ii_tt = 1:length(NTTs)
    % we might have two spikes in the same bin. in such case we consider 
    % the features of the last spike
    [N,~,BIN] = histcounts(NTTs(ii_tt).ts,edges);    
    multiunits_spikes(N>0,ii_tt) = true;
    features = NTTs(ii_tt).features;
    valid = BIN~=0;
    BIN(~valid)=[];
    features(:,~valid)=[];
    multiunits(BIN,:,ii_tt) = features'; 
end

%% arrange data into one simple struct
data = struct();
data.exp_ID = exp.details.exp_ID;
data.fs = fs;
data.dt = dt;
data.t = t;
data.multiunits_spikes = multiunits_spikes;
data.multiunits = multiunits;
data.position = position;
data.direction = direction;
data.FE_ti = FE_ti;
data.sleep_ti = sleep_ti;
data.PE_ts = [exp.PE.thr.peak_ts];
data.PE_ti = [exp.PE.thr.start_ts; exp.PE.thr.end_ts]';
data.rest_ti = rest_ti;

%% save data to mat file
mkdir(data_proc_dir);
file_out = fullfile(data_proc_dir, details.exp_ID);
save(file_out, '-struct', 'data');


%% check what is removed when we remove low-speed
FE = exp.flight.FE;
FE = FE([FE.distance] > 100);
figure
hold on
plot([FE.pos],[FE.vel],'.k');
FE = FE_remove_low_speed(FE,0.85);
plot([FE.pos],[FE.vel],'.r');
xlabel('position (m)')
ylabel('Velocity (m/s)')
legend('original','after removing low speed','Location','Best')
title('remove low speed tails from flights')
file_out = fullfile(data_proc_dir, [details.exp_ID '_FE_remove_low_speed']);
saveas(gcf, file_out, 'jpeg');

%%
function FE = FE_remove_low_speed(FE,thr_prc)
    for ii_FE = 1:length(FE)
        fe = FE(ii_FE);
        thr = abs(thr_prc * median([fe.vel]));
        start_relative_IX = find(abs([fe.vel])>thr,1, 'first');
        end_relative_IX = find(abs([fe.vel])>thr,1, 'last');
        fe.ts = fe.ts(start_relative_IX:end_relative_IX);
        fe.pos = fe.pos(start_relative_IX:end_relative_IX);
        fe.vel = fe.vel(start_relative_IX:end_relative_IX);
        fe.start_ts = fe.ts(1);
        fe.end_ts = fe.ts(end);
        fe.duration = range(fe.ts);
        fe.distance = range(fe.pos);
        FE(ii_FE) = fe;
    end
end




%% ===============================================================
%% check how many spikes we loose if we increase the time bins
% FE = exp.flight.FE;
% FE = FE([FE.distance] > 100);
% FE_ti = [[FE.start_ts];[FE.end_ts]]';
% 
% dt_options = 1:10;
% fs_options = 1000./dt_options;
% lost_spikes_n = zeros(3,length(fs_options));
% lost_spikes_prc = zeros(3,length(fs_options));
% ii_tt = 2;
% mean_fr = 1e6*length(NTTs(2).ts) / range(t);
% for ii_opt = 1:length(fs_options)
%     fs = fs_options(ii_opt);
%     dt=1e6/fs;
%     tstart = min(exp.details.session_ts(:));
%     tend = max(exp.details.session_ts(:));
%     edges = tstart : dt: tend;
%     t = edges2centers(edges);
% 
%     is_FE = any(t>FE_ti(:,1)&t<FE_ti(:,2),1);
%     is_sleep = any(t>sleep_ti(:,1)&t<sleep_ti(:,2),1);
%     
%     [N,~,BIN] = histcounts(NTTs(ii_tt).ts,edges);
%     
%     sdf1 = (           N>1) .* (N-1);
%     sdf2 = (is_FE    & N>1) .* (N-1);
%     sdf3 = (is_sleep & N>1) .* (N-1);
%     lost_spikes_n(1, ii_opt) = sum(sdf1);
%     lost_spikes_n(2, ii_opt) = sum(sdf2);
%     lost_spikes_n(3, ii_opt) = sum(sdf3);
%     n1 = sum(N);
%     n2 = sum(N .* is_FE);
%     n3 = sum(N .* is_sleep);
%     lost_spikes_prc(1, ii_opt) = sum(sdf1) / n1;
%     lost_spikes_prc(2, ii_opt) = sum(sdf2) / n2;
%     lost_spikes_prc(3, ii_opt) = sum(sdf3) / n3;
% end
% 
% figure
% subplot(121)
% plot(dt_options,lost_spikes_n','.-')
% xlabel('time bin (ms)')
% ylabel('No. of lost spikes (counts)')
% legend('Entire recording','flight','sleep')
% subplot(122)
% plot(dt_options,100.*lost_spikes_prc','.-')
% xlabel('time bin (ms)')
% ylabel('Percentage of lost spikes (%)')
% sgtitle({ 'Quantifying how many spikes we loose if we increase time bins';
%         "analyzing exp " + exp_ID + " using a tetrode with high mean firing rate: " + mean_fr + "Hz"},...
%         'Interpreter','None');
% legend('Entire recording','flight','sleep')

%%



