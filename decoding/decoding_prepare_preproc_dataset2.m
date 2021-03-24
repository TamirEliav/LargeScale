%%
clear
clc

%%
main_data_dir = "D:\sequences\seq_uri_eden\DATA";
data_proc_dir = "D:\sequences\seq_uri_eden\proc";
exp_ID = 'b9861_d180524';
exp = exp_load_data(exp_ID);
details = exp.details;

%% prepare file names
NTT_files = {};
tt_to_use = ones(1,exp.details.numTT);
tt_to_use = tt_to_use & exp.details.TT_to_use;
tt_to_use = tt_to_use & contains(exp.details.TT_loc,{'CA1','CA3'});
tt_to_use(4)=false;

data_dir = fullfile(main_data_dir,num2str(details.batNum),'spikes_detection');
filenames = "spikes_"+exp_ID+"_TT"+find(tt_to_use)'+".ntt";
NTT_files = fullfile(data_dir,filenames);

% NTT_files = {
%     'D:\sequences\seq_uri_eden\DATA\9861\spikes_detection\spikes_b9861_d180524_TT1.NTT',...
%     'D:\sequences\seq_uri_eden\DATA\9861\spikes_detection\spikes_b9861_d180524_TT2.NTT',...
%     'D:\sequences\seq_uri_eden\DATA\9861\spikes_detection\spikes_b9861_d180524_TT3.NTT',...
%     'D:\sequences\seq_uri_eden\DATA\9861\spikes_detection\spikes_b9861_d180524_TT4.NTT',...
%     };


%% prepare behavioral data
FE = exp.flight.FE;
% FE = FE([FE.direction]==1);
FE = FE([FE.distance] > 100);
ti = [[FE.start_ts];[FE.end_ts]]';

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
    data(ii_file) = spikes;
end

%% arange data
fs=1000;
dt=1e6/fs;
t=[];
for ii=1:size(ti,1)
    t = [t ti(ii,1): dt : ti(ii,2)];
end
multiunits = nan(length(t),4,length(data));
multiunits_spikes = false(length(t),length(data));
for tt = 1:length(data)
    %%
    IX = find_nearest_point(t,data(tt).ts);
    tdiff = abs(t-data(tt).ts(IX));
    TF = tdiff < dt/2;
    multiunits_spikes(TF,tt) = true;
    multiunits(TF,:,tt) = data(tt).features(:,IX(TF))';
end
position = interp1(exp.pos.proc_1D.ts, exp.pos.proc_1D.pos, t);
direction = interp1(exp.pos.proc_1D.ts, sign(exp.pos.proc_1D.vel_csaps), t);
position = [position;direction];

%%
t2 = 1:length(t);
t2 = t2./fs;

%% save data to mat file
mkdir(data_proc_dir);
file_out = fullfile(data_proc_dir, details.exp_ID);
save(file_out, 't','t2','fs','position','multiunits_spikes','multiunits');






%%






