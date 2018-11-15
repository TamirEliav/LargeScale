function [signal, ts, fs] = Nlx_csc_read(file_in, limits_ts)



if isempty(limits_ts)
    [ts_blocks, fs, samples_blocks, header] = Nlx2MatCSC(file_in, [1 0 1 0 1], 1, 1);
else
    [ts_blocks, fs, samples_blocks, header] = Nlx2MatCSC(file_in, [1 0 1 0 1], 1, 4, limits_ts);
end

% re-arrange data samples to a vector of uVolt samples
signal = reshape(samples_blocks, 1, prod(size(samples_blocks)));
% fs = fs(1);
field_name = 'SamplingFrequency';
temp=textscan(header{find(not(cellfun('isempty', strfind(header, field_name)))) , :},'%s%f');
fieled_value = temp{2};
fs = fieled_value;

if 0
% TODO: check if the following params exist in the header, if not, use defaults
% read some header parameters values
field_name = 'ADBitVolts';
temp=textscan(header{find(not(cellfun('isempty', strfind(header, field_name)))) , :},'%s%f');
fieled_value = temp{2};
ADBitVolts = fieled_value;
field_name = 'InputInverted';
temp=textscan(header{find(not(cellfun('isempty', strfind(header, field_name)))) , :},'%s %s');
fieled_value = temp{2}{:};
is_data_inverted = strcmp('True',fieled_value );

signal = signal .* ADBitVolts; % change to volt
signal = signal .* 1e6; % change to uvolt
if is_data_inverted
    signal = -signal; % invert back
end
end

sample_time_usec = 1e6/fs;
block_size = size(samples_blocks,1);
block_relative_time = linspace(0, sample_time_usec*(block_size-1), block_size );
sdf1 = repmat(ts_blocks,block_size,1);
sdf2 = repmat(block_relative_time', 1, length(ts_blocks));
ts = sdf1 + sdf2;
ts = reshape(ts, 1, prod(size(ts)));























