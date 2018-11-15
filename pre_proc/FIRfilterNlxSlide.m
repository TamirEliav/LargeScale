function [Samples2,Timestamps2,SampleFrequency]= FIRfilterNlxSlide(file_IN,file_OUT,t_start_end,passband,slidwin)

%% read file
[signal, ts, fs] = Nlx_csc_read(file_IN, t_start_end);

[ts_blocks,~,~,~,Samples,NlxHeader] = Nlx2MatCSC(file_IN, FieldSelection,1,ExtractMode,ExtractionModeArray) ;

fs_row_IX = find(~cellfun(@isempty, strfind(NlxHeader, '-SamplingFrequency')));
fs_str = NlxHeader(fs_row_IX);
sdf = regexp(fs_str , ' ', 'split');
sdf = sdf{1}(2);
SampleFrequency = str2num(sdf{1});
% SampleFrequency=str2num(NlxHeader{8}(20:27));
SamplePeriod_usec = 1e6/SampleFrequency;

block_size = size(Samples,1);
block_relative_time = linspace(0, SamplePeriod_usec*(block_size-1), block_size );
sdf1 = repmat(ts_blocks,block_size,1);
sdf2 = repmat(block_relative_time', 1, length(ts_blocks));
Timestamps = sdf1 + sdf2;

Timestamps = reshape(Timestamps, 1, prod(size(Timestamps)));
Samples = reshape(Samples, 1, prod(size(Samples)));

% Timestamps  = NlxfilledTimstamps( Timestamps ); 
% Timestamps = Timestamps(:); 
% Samples = Samples(:); 
% Samples  = Samples- mean(Samples); 

fwin = slidwin*60*10^6;                                              % we will run over the data in 2-min windows but save                                                               % only the central 1-min to avoid edge effects in the filtering process (see below)
fwin = round(fwin/SamplePeriod_usec);
rwin = round(fwin/2);                                          
fwin = 2*rwin; 

v= 0:rwin:length(Samples)-fwin; 
v =repmat(v(:),1,fwin); 
t = 0:1:fwin-1; 
t = repmat(t,size(v,1),1); 
v = v+t; 
v = v+1;
 
Timestamps = Timestamps(v); 
v = Samples(v); 

filtOrder   = 300 ; 
Wn                = passband / (SampleFrequency/2);           
filt              = fir1(filtOrder, Wn,'bandpass' );       

parfor i=1:size(v,1)
    v(i,:)  = filtfilt(filt,1,v(i,:)')'; 
end

Samples2 = v(:,rwin/2:end-(rwin/2)-1); 
Timestamps2 = Timestamps(:,rwin/2:end-(rwin/2)-1); 
Samples2 = permute(Samples2,[2 1]); 
Samples2= Samples2(:); 
Timestamps2 = permute(Timestamps2,[2 1]); 
Timestamps2 = Timestamps2(:); 





