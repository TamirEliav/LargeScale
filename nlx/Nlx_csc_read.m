function [signal, ts, fs, params] = Nlx_csc_read(file_in, limits_ts)



if isempty(limits_ts)
    [ts_blocks, fs, samples_blocks, header] = Nlx2MatCSC(file_in, [1 0 1 0 1], 1, 1);
else
    [ts_blocks, fs, samples_blocks, header] = Nlx2MatCSC(file_in, [1 0 1 0 1], 1, 4, limits_ts);
end

% re-arrange data samples to a vector of uVolt samples
signal = reshape(samples_blocks, 1, prod(size(samples_blocks)));

%% parse header
SamplingFrequency   = textscan(header{contains(header,'SamplingFrequency')},    '%s %f');
ADMaxValue          = textscan(header{contains(header,'ADMaxValue')},           '%s %d');
InputRange          = textscan(header{contains(header,'InputRange')},           '%s %f');
ADBitVolts          = textscan(header{contains(header,'ADBitVolts')},           '%s %f');
InputInverted       = textscan(header{contains(header,'InputInverted')},        '%s %s');
SamplingFrequency =  SamplingFrequency{2};
ADMaxValue        =   ADMaxValue{2};
InputRange        =   InputRange{2};
ADBitVolts        =   ADBitVolts{2};
InputInverted     =   InputInverted{2};

%% arrange output data
fs = SamplingFrequency;
signal = signal .* ADBitVolts; % change to volt
signal = signal .* 1e6; % change to uvolt

sample_time_usec = 1e6/fs;
block_size = size(samples_blocks,1);
block_relative_time = linspace(0, sample_time_usec*(block_size-1), block_size );
sdf1 = repmat(ts_blocks,block_size,1);
sdf2 = repmat(block_relative_time', 1, length(ts_blocks));
ts = sdf1 + sdf2;
ts = reshape(ts, 1, prod(size(ts)));

%% arrange output params
params.SamplingFrequency = SamplingFrequency;
params.ADMaxValue = ADMaxValue;
params.InputRange = InputRange;
params.ADBitVolts = ADBitVolts;
params.InputInverted = InputInverted;
params.header = header;



%%
