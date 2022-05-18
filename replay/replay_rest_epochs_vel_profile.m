%%
clear
clc

%% params
win_s = 1;

%% choose bats / sessions
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
clear exp_list
groupsummary(T,'bat_num')
if exist('bats_to_include','var')
    T = groupfilter(T,"bat_num",@(x)ismember(x,bats_to_include),'bat_num');
end
bats = unique(T.bat_num);

%% arrange data
exps = cellfun(@(exp_ID)(exp_load_data(exp_ID,'details','rest')), T.exp_ID);
rests = [exps.rest];
epochs = [rests.events];
fs = rests(1).fs;
win_samples = round(win_s*fs);

%% trigger vel
sdf = nan(length(epochs), 2*win_samples+1);
sdf2 = false(1,length(epochs));
ii = 0;
epoch_type = 'rest';
params_opt = 11;
for ii_exp = 1:length(rests)
    exp = exps(ii_exp);
    exp_ID = exp.details.exp_ID;
    rest = exp.rest;
    nEpochs = length(rest.events);
    
    events = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
    seqs = [events.seq_model];
    [seqs, TF] = decoding_apply_seq_inclusion_criteria(seqs);
    events(~TF)=[];
    gEarlyLate = categorical([events.time_from_epoch_start]<1,[false,true],["Late","Early"]);
    rest_epochs_num_with_early_replay_events = unique([events(gEarlyLate=="Early").epoch_num]);
    sdf2(ii+[1:nEpochs]) = ismember(1:nEpochs,rest_epochs_num_with_early_replay_events);

%     trigger_IX = start_IXs + [-win_samples:win_samples];
%     trigger_IX(trigger_IX<1) = nan;
%     sdf2 = rest.vel_smooth(trigger_IX);
    sdf(ii+[1:nEpochs], :) = trigger_signal_by_IX(abs(rest.vel_smooth), [rest.events.start_IX], win_samples);
    ii = ii + nEpochs;
end
% gHasEarlyEvents = categorical(sdf2,[true false],["with_early","without_early"]);

%%
figure
hold on
% plot( nanmean(sdf(gHasEarlyEvents=='with_early',:)) );
% plot( nanmean(sdf(gHasEarlyEvents=='without_early',:)) );
plot( nanmean(sdf( sdf2,:)) );
plot( nanmean(sdf(~sdf2,:)) );
plot( nanmean(sdf) ,'k');
rescale_plot_data('x',[1/fs win_samples+1]);
xlabel('Time (s)');
ylabel('Speed (m/s)');
legend("with early replay","without early replay","all");

%%
figure
imagesc(sdf)
% imagesc(sdf(sdf2,:))
% imagesc(sdf(~sdf2,:))
rescale_plot_data('x',[1/fs win_samples+1])
xlabel('Time (s)')
ylabel('#epoch')
colorbar
hax=gca;
hax.CLim = [0 0.5];


%%
