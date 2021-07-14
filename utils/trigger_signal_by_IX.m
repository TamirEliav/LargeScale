function [trig_signal] = trigger_signal_by_IX(signal, trigger_IX, window_num_samples, fill_value)
% This function triggers a signal by given set of timestamps.
% 
% IN:
% signal - 1xN samples
% trigger_IX - 1xK triggering points
% window - if scalar, window is symmetric, if 2 elements vector window can
% be asymmetric
% 
% OUT:
% triggerd_signal - KxN matrix 
% 
% COMMENTS:
% - For triggering points were singal is not available, nans are used.
% - Assuming continuous signal

%%
trig_signal = [];
switch numel(window_num_samples)
    case 1
        window_num_samples = [-window_num_samples window_num_samples];
    case 2
%         TODO:sanity check for positive/negative values
    otherwise
            error('wrong number of elements in triggering window!')
end
if nargin<4
    fill_value = nan;
end
    
%%
trig_win_IX_relative = window_num_samples(1):window_num_samples(2);
sdf1 = repmat(trig_win_IX_relative, length(trigger_IX), 1);
sdf2 = repmat(trigger_IX, length(trig_win_IX_relative), 1)';
trig_IX = sdf1 + sdf2;
no_data_trig_IX = find( trig_IX(:,1) < 1  | trig_IX(:,end) > length(signal));
trig_IX(no_data_trig_IX,:) = 1; % put dummy index (later fill with nans or user input)
trig_signal = signal(trig_IX);
trig_signal(no_data_trig_IX,:) = fill_value;


%%


end