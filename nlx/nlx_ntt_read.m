function NTT = nlx_ntt_read(filename,ti)
arguments
    filename
    ti = [];
end

if isempty(ti)
    [Timestamps, CellNumbers, Samples, Header] = ...
         Nlx2MatSpike(filename, [1 0 1 0 1], 1, 1, [] );
else
    [Timestamps, CellNumbers, Samples, Header] = ...
         Nlx2MatSpike(filename, [1 0 1 0 1], 1, 4, ti );
end

%% parse header
ADBitVolts = sscanf(Header{contains(Header,'ADBitVolts')},'-ADBitVolts %f %f %f %f');
fs = sscanf(Header{contains(Header,'SamplingFrequency')},'-SamplingFrequency %f %f %f %f');

%% convert bits to uVolts
Samples = Samples .* ADBitVolts' .* 1e6;

%% create spikes structure
NTT.ts = Timestamps;
NTT.waveforms = Samples;
NTT.cluster_id = CellNumbers;
NTT.filename = filename;
NTT.header = Header;

end