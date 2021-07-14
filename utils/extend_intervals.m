function xthr_high_low = extend_intervals(xthr_high,xthr_low)

%% testing
% xthr_high = [0 0 1 0 0 0 1 0 0 0];
% xthr_low  = [0 1 1 1 0 0 1 1 0 1];

%% make sure xthr_low is true wherever xthr_high is true
% xthr_low(xthr_high) = 1;

%%
xthr_low_lbl = bwlabel(xthr_low);
xthr_low_high_lbl_valid = xthr_low_lbl .* xthr_high;
valid_lbl = unique(xthr_low_high_lbl_valid(xthr_low_high_lbl_valid>0));
xthr_high_low = xthr_low;
xthr_high_low(~ismember(xthr_low_lbl,valid_lbl))=0;

end