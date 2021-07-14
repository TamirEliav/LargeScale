function [LFP, ts, fs, params, ch_valid] = LFP_load(exp_ID,band,TT_or_ch_to_use)
% band is optional (as one of the following):
% delta
% theta
% low_gamma
% high_gamma
% ripple
% Third argument is optional and can be list of specific TT/ch to load
% (vector/matrix). If matrix in the form TT X ch of booleans
% out: 
%   LFP (time X TT X ch)

%% get exp info
exp = exp_load_data(exp_ID,'details','path');
prm = PARAMS_GetAll();
active_channels = exp.details.activeChannels;

%% get folder
if ~exist('band','var')
    dir_IN = exp.path.LFP;
else
    dir_IN = fullfile(exp.path.LFP_bands,band);
end

%% option for specific TT/ch
if exist('TT_or_ch_to_use','var')
    if isvector(TT_or_ch_to_use)
        TT = TT_or_ch_to_use;
        mask = false(size(active_channels));
        mask(TT,:) = true;
    elseif ismatrix(TT_or_ch_to_use)
        mask = TT_or_ch_to_use;
    else 
        error('Wrong input')
    end
    ch_valid = active_channels & mask;
else
    ch_valid = active_channels;
end

%% load using parfor - turns out it is slower... probably due to reading from the hard-drive bottelneck 
% % % tic
% % % nTT = size(ch_valid,1);
% % % nCh = size(ch_valid,2);
% % % files = fullfile(dir_IN, "LFP_"+exp_ID+"_TT"+[1:nTT]'+"_ch"+[1:nCh]+".ncs");
% % % % LFP=[]; % TODO: preallocate!
% % % % LFP = zeros(16,13388288);
% % % % LFP = zeros(16,13388288);
% % % % LFP = zeros(13388288,16);
% % % ts=[];
% % % fs=[];
% % % params=[];
% % % parfor ii_file = 1:numel(files)
% % %     if ~ch_valid(ii_file)
% % %         LFP(:,ii_file) = nan;
% % %         continue;
% % %     end
% % %     LFP_file_IN = str2mat(files(ii_file));
% % % %     fprintf('loading file %d\n',ii_file)
% % % %     [LFP(ii_file,:), ts, fs, params] = Nlx_csc_read(LFP_file_IN , []);
% % % %     [LFP(:,ii_file)] = Nlx_csc_read(LFP_file_IN , []);
% % %     Nlx_csc_read(LFP_file_IN , []);
% % % %     fprintf('file %d was loaded\n',ii_file)
% % % end
% % % toc

%% get data length
[TT,ch] = find(ch_valid,1,'first');
LFP_file_IN = fullfile(dir_IN, ['LFP_' exp_ID '_TT' num2str(TT) '_ch' num2str(ch) '.ncs']);
[~, ts, ~, ~] = Nlx_csc_read(LFP_file_IN , []); 

%%
nPoints = length(ts);
nTT = size(ch_valid,1);
nCh = size(ch_valid,2);
LFP = zeros(nPoints,nTT,nCh);
for TT = 1:nTT
    for ch = 1:nCh
        if ~ch_valid(TT,ch)
            LFP(:,TT,ch) = nan;
            continue;
        end
        LFP_file_IN = fullfile(dir_IN, ['LFP_' exp_ID '_TT' num2str(TT) '_ch' num2str(ch) '.ncs']);
        % place data in columns as this is how matlab store it continuously 
        % in the memory, so it's slightly faster
        [LFP(:,TT,ch), ts, fs, params] = Nlx_csc_read(LFP_file_IN , []); 
    end
end

end



