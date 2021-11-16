function decoding_prepare_exp_data(exp_ID)

%% get exp info
dir_out = "F:\sequences\proc";
% exp = exp_load_data(exp_ID,'details','path','pos','flight','PE','rest');
exp = exp_load_data(exp_ID,'details','path','pos','flight','rest');
details = exp.details;

%% get sleep ts
sleep_ti = exp_get_sessions_ti(exp_ID, 'Sleep1','Sleep2');
sleep_ti(any(isnan(sleep_ti),2),:) = []; % remove nan in case of missing sessions

%% prepare NTT file names
NTT_files = {};
tt_to_use = ones(1,exp.details.numTT);
tt_to_use = tt_to_use & ismember(1:exp.details.numTT,exp.details.TT_to_use);
tt_to_use = tt_to_use & contains(exp.details.TT_loc,{'CA1','CA3'});

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

% behavioral data
position = interp1(exp.pos.proc_1D.ts, exp.pos.proc_1D.pos, t);
direction = interp1(exp.pos.proc_1D.ts, sign(exp.pos.proc_1D.vel_csaps), t);

% remove low speed from flight data
thr_prc = 0.85;
keep_takeoff_landing = 1;
balls_loc = exp.rest.balls_loc;
max_dist_from_ball = 1;
FE = exp.flight.FE;
FE = FE([FE.distance] > 100);
FE_before = FE;
FE = FE_remove_low_speed(FE,thr_prc,keep_takeoff_landing,balls_loc,max_dist_from_ball);
FE_after = FE;
FE_ti = [[FE.start_ts];[FE.end_ts]]';
FE_dist = [FE.distance];

is_FE = any(t>FE_ti(:,1)&t<FE_ti(:,2),1);
is_sleep = any(t>sleep_ti(:,1)&t<sleep_ti(:,2),1);

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
data.rest_ti = exp.rest.ti;
% data.PE_ts = [exp.PE.thr.peak_ts];
% data.PE_ti = [exp.PE.thr.start_ts;exp.PE.thr.end_ts]';

%% save data to mat file
mkdir(dir_out);
file_out = fullfile(dir_out, exp_ID);
save(file_out, '-struct', 'data');

%% plot the spikes features used to train (during flight)
ch_pairs = [
    1 2;
    1 3;
    1 4;
    2 3;
    2 4;
    3 4];
for TT=1:size(multiunits,3)
    %%
    TF_FE = multiunits_spikes(:,TT) & is_FE';
    TF_sleep = multiunits_spikes(:,TT) & is_sleep';
    % TODO: add also data points for rest on the balls
    figure('WindowState','maximized');
    tiledlayout(2,3,'TileSpacing','compact');
    for ii_pnl = 1:size(ch_pairs,1)
        ch_pair = ch_pairs(ii_pnl,:);
        nexttile
        hold on
        plot(multiunits(TF_FE,ch_pair(1),TT), ...
             multiunits(TF_FE,ch_pair(2),TT), '.r');
         plot(multiunits(TF_sleep,ch_pair(1),TT), ...
             multiunits(TF_sleep,ch_pair(2),TT), '.b');
         xlabel("Ch "+ch_pair(1))
         ylabel("Ch "+ch_pair(2))
         legend('flight','sleep')
    end
    % TODO: TT is the INDEX of the tetorde! correct it!
    sgtitle({exp_ID;"TT "+TT},'interpreter','none')
    file_out = fullfile(dir_out, [exp_ID '_features_TT_' num2str(TT)]);
    saveas(gcf, file_out, 'jpeg');
end


%% check what is removed when we remove low-speed
figure
hold on
plot([FE_before.pos],[FE_before.vel],'.k');
plot([FE_after.pos],[FE_after.vel],'.r');
xlabel('position (m)')
ylabel('Speed (m/s)')
legend('original','after removing low speed','Location','Best')
title({exp_ID;'Remove low speed tails from flights'},'Interpreter','none')
file_out = fullfile(dir_out, [details.exp_ID '_FE_remove_low_speed']);
saveas(gcf, file_out, 'jpeg');

end


%%
function FE = FE_remove_low_speed(FE,thr_prc,keep_takeoff_landing,balls_loc,max_dist_from_ball)
    keep_takeoff_landing = 1;
    for ii_FE = 1:length(FE)
        fe = FE(ii_FE);
        thr = abs(thr_prc * median([fe.vel])); % median of the single flight vel
        start_relative_IX = find(abs([fe.vel])>thr,1, 'first');
        end_relative_IX = find(abs([fe.vel])>thr,1, 'last');
        if keep_takeoff_landing
            if min(abs(fe.pos(1) - balls_loc)) < max_dist_from_ball
                start_relative_IX = 1;
            end
            if min(abs(fe.pos(end) - balls_loc)) < max_dist_from_ball
                end_relative_IX = length(fe.vel);
            end
        end
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

%%



