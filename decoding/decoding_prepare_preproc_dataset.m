function decoding_prepare_preproc_dataset(exp_ID)

%%
clear
clc
% exp_ID = 'b9861_d180616';
% exp_ID = 'b9861_d180619';
exp_ID = 'b9861_d180705';
% exp_ID = 'b9861_d180708';


%% read data
exp=exp_load_data(exp_ID,'details','path','pos');

%%
chunk_size_sec = 5;
avg_win = 1000;

%%
for ii_dir=1:2
    disp("dir "+ii_dir);
    %% read raw data
    main_dir = 'L:\Analysis\Results\decoding\raw_datasets';
    raw_file = fullfile(main_dir,exp_ID+"_raw_data_dir_"+ii_dir+".h5");
    fs = h5read(raw_file,'/raw_fs');
    ts = h5read(raw_file,'/raw_ts');
    info = h5info(raw_file);
    raw_info = info.Datasets( strcmp({info.Datasets.Name},'raw') );
    num_ch = raw_info.Dataspace.Size(2);
    num_win=floor(length(ts)/avg_win);
    
    %% create filter bank
    chunk_size_samples = round(chunk_size_sec*fs);
    nChunks = floor( length(ts)/chunk_size_samples);
    fb = cwtfilterbank( 'SignalLength', chunk_size_samples, ...
                        'Wavelet','amor',...
                        'SamplingFrequency',fs,...
                        'VoicesPerOctave',4,...
                        'FrequencyLimits',[100 10000]);
    [~,freqs] = cwt(zeros(1,chunk_size_samples),'FilterBank',fb);
    
    %% calc sub-sampled ts
    ts2 = reshape(ts(1:num_win*avg_win),avg_win,[]);
    ts2 = mean(ts2,1);
    
    %% create h5 out file
    proc_file = fullfile(main_dir,exp_ID+"_pre_proc_dir_"+ii_dir+".h5");
    % TODO: handle existing file!!!!!
    h5create(proc_file ,'/inputs/fourier_frequencies',[length(freqs)]);
    h5write(proc_file ,'/inputs/fourier_frequencies',freqs);
    h5create(proc_file ,'/inputs/ts',[length(ts)]);
    h5write(proc_file ,'/inputs/ts',ts);
    h5create(proc_file ,'/inputs/wavelets', [num_ch length(freqs) num_win]);
    h5create(proc_file ,'/outputs/position',[                  1  num_win]);
    h5create(proc_file ,'/outputs/speed',   [                  1  num_win]);
    
    
    %% write input (wavelets)
    for chnl = 1:num_ch
        %%
%         disp("channel "+chnl)
        chnl_data=zeros(length(freqs),length(ts));
        for chnk = 1:nChunks
            %%
%             disp("chunk "+chnk)
            fprintf('dir %d, channel %d, chunk %d\n',ii_dir,chnl,chnk);
            start = (chnk-1)*chunk_size_samples+1;
            count = chunk_size_samples;
            chnk_data=h5read(raw_file,'/raw',[start chnl],[count 1]);
            chnk_data = gpuArray(single(chnk_data));
            chnl_data(:,start+[0:count-1]) = gather(abs(cwt(chnk_data ,'FilterBank',fb)));
            
        end
        chnl_data = reshape(chnl_data(:,1:num_win*avg_win),length(freqs),num_win,[]);
        chnl_data = mean(chnl_data,3);
        chnl_data = shiftdim(chnl_data,-1);
        start = [chnl 1 1];
        count = size(chnl_data);
        h5write(proc_file ,'/inputs/wavelets', chnl_data, start, count);
    end
    
    %% write output variables to hd5 file
    % first, sub-sample output variables
    pos = h5read(raw_file,'/pos');
    speed = abs(h5read(raw_file,'/vel'));
    pos_ts = h5read(raw_file,'/pos_ts');
    pos2 = interp1(pos_ts,pos, ts2);
    speed2 = interp1(pos_ts,speed, ts2);
    % now, write to file
    h5write(proc_file ,'/outputs/position', pos2);
    h5write(proc_file ,'/outputs/speed', speed2);
end
                    



%%





%%
% x=chnl_data_avg;
% x=reshape(x,size(x,1)/2,2,[]);
% x=squeeze(mean(x,2));
% x=normalize(x,2);
% noise = mean(abs(diff(x')));
% [~,order]=sort(noise);
% plot(x(order,:)')
% plotbrowser

%%






%%