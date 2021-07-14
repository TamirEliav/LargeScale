function PE_detect(exp_ID)

%% get exp info
exp = exp_load_data(exp_ID,'details','path','ripples','MUA');
prm = PARAMS_GetAll();

%%
rpl = exp.ripples.all;
mua = exp.MUA.events;
t = exp.MUA.t; % we will use MUA timestamps

%%
PE = merge_ripples_and_MUA(rpl,mua,t);
PE2 = merge_ripples_and_MUA(rpl(g1==1),mua(g2==1),t);

%% save population events ts to nlx event files
Filename = fullfile(exp.path.LFP_bands,'ripple', [exp_ID '_PE' '.nev']);
Mat2NlxEV( Filename, 0, 1, [], [1 0 0 0 0 0], [PE.peak_ts]);

%% save ripples detection results
file_name = fullfile('L:\Analysis\Results\exp\PE',[exp_ID '_exp_PE']);
save(file_name,'PE');

end

function PE = merge_ripples_and_MUA(rpl,mua,t)

%%
mua_ti = [mua.start_ts; mua.end_ts]';
rpl_ti = [rpl.start_ts; rpl.end_ts]';
is_mua = any(t>=mua_ti(:,1)&t<=mua_ti(:,2),1);
is_rpl = any(t>=rpl_ti(:,1)&t<=rpl_ti(:,2),1);

%%
mua_lbl = bwlabel(is_mua);
% max(mua_lbl)
if max(mua_lbl) ~= length(mua)
    error('Something went wrong...');
end
valid_events = unique(mua_lbl .* is_rpl);
valid_events(valid_events==0)=[];
% length(valid_events) / length(mua.events)
PE = mua(valid_events);
PE = rmfield(PE,'start_IX');
PE = rmfield(PE,'end_IX');
PE = rmfield(PE,'peak_IX');
% TODO: add fields for ripples quantifications (z/ripple-gamma-ratio/...)

end