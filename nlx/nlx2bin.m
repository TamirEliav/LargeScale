function nlx2bin(exp_ID, dir_out)
arguments
    exp_ID
    dir_out = 'D:\DATA\raw_bin';
end

%% get exp info
exp = exp_load_data(exp_ID,'details','path');
active_channels = exp.details.activeChannels;
active_channels = active_channels';
active_channels = active_channels(:);

%% init ouput array
%  dummy read of first valid channel to get the data length
ii_ch = find(active_channels, 1, 'first');
TT = ceil(ii_ch/4);
ch_num = mod(ii_ch-1,4)+1;
dummy_file = fullfile(exp.path.nlx,['CSC_TT' num2str(TT) '_' num2str(ch_num) '.ncs']);
dummy_data = Nlx_csc_read(dummy_file,[]);
% init out array
data = zeros(length(active_channels), length(dummy_data), 'int16');
file_OUT = fullfile(dir_out, sprintf('%s_raw_neur.bin',exp_ID));

%% read data (from .ncs)
for ii_ch = 1:length(active_channels)
    TT = ceil(ii_ch/4);
    ch_num = mod(ii_ch-1,4)+1;
    file_IN = fullfile(exp.path.nlx,['CSC_TT' num2str(TT) '_' num2str(ch_num) '.ncs']);
    data(ii_ch,:) = int16(Nlx_csc_read(file_IN,[]));
end

%% write data (to .bin)
fid = fopen(file_OUT, 'w'); 
fwrite(fid, data, 'int16');
fclose(fid);

end