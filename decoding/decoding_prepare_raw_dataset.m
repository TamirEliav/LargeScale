function decoding_prepare_raw_dataset(exp_ID)

%%
% clear
% clc
% exp_ID = 'b9861_d180616';
% exp_ID = 'b9861_d180619';
% exp_ID = 'b9861_d180705';
% exp_ID = 'b9861_d180708';


%%
exp=exp_load_data(exp_ID,'details','path','pos','flight');

%% find active channels
active_channels = find(exp.details.activeChannels')-1;

%% read first ncs file to get neural data ts
file_in = fullfile(exp.path.nlx, sprintf('CSC%d.ncs',active_channels(1)) );
[signal, ts, fs, params] = Nlx_csc_read(file_in, []);

%%
FE=[exp.flight.FE];
FE([FE.distance]<100)=[];
directions=[1 -1];
for ii_dir=1:2
    disp(ii_dir)
    % get data
    FE_dir = FE([FE.direction]==directions(ii_dir));
    FE_ts=[[FE_dir.start_ts];[FE_dir.end_ts]]';
    ts_IX = get_data_in_ti(ts,FE_ts);
    
    % read raw data
    % TODO: pre-allocate memory! (here in matlab or better directly in the
    % hd5 file!!)
    raw_data=[];
    for ii_ch = 1:length(active_channels)
        disp(ii_ch)
        file_in = fullfile(exp.path.nlx, sprintf('CSC%d.ncs',active_channels(ii_ch)) );
        signal = Nlx_csc_read(file_in, []);
        raw_data(:,ii_ch) = signal(ts_IX);
    end
    
    % write data to hd5 file
    out_dir = 'L:\Analysis\Results\decoding\raw_datasets';
    out_file = fullfile(out_dir,exp_ID+"_raw_data_dir_"+ii_dir+".h5");
    h5create(out_file ,'/raw',size(raw_data));
    h5write(out_file , '/raw', raw_data);
    h5create(out_file ,'/raw_ts',size(ts_IX));
    h5write(out_file , '/raw_ts', ts(ts_IX));
    h5create(out_file ,'/raw_fs',size(fs));
    h5write(out_file , '/raw_fs', fs);
    h5create(out_file ,'/pos',size(exp.pos.proc_1D.pos_csaps));
    h5write(out_file , '/pos', exp.pos.proc_1D.pos_csaps);
    h5create(out_file ,'/vel',size(exp.pos.proc_1D.vel_csaps));
    h5write(out_file , '/vel', exp.pos.proc_1D.vel_csaps);
    h5create(out_file ,'/pos_ts',size(exp.pos.proc_1D.ts));
    h5write(out_file , '/pos_ts', exp.pos.proc_1D.ts);
    h5disp(out_file )
end



%%











%%





%%


