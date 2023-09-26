%%
clear
clc

%% choose bats / sessions
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
clear exp_list
groupsummary(T,'bat_num')
if exist('bats_to_include','var')
    T = groupfilter(T,"bat_num",@(x)ismember(x,bats_to_include),'bat_num');
end
bats = unique(T.bat_num);

%% load and prepare data
nbins = 50;
coverage_all = zeros(height(T),2,2,nbins); % [exp] X [sleep/rest] X [2 directions] X [space]
ccc_all = zeros(height(T),2); % [exp] X [2 directions]
n_seqs_all = zeros(height(T),2,2); % [exp] X [sleep/rest] X [2 directions]
for ii_exp = 1:height(T)
    % load exp data
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'details');

    epoch_type = 'sleep';
    params_opt = 11;
    [coverage_all(ii_exp,1,:,:), n_seqs_all(ii_exp,1,:)] = calc_coverage(exp_ID, epoch_type, params_opt,nbins);
    epoch_type = 'rest';
    params_opt = 11;
    [coverage_all(ii_exp,2,:,:),  n_seqs_all(ii_exp,2,:)]= calc_coverage(exp_ID, epoch_type, params_opt,nbins);
    ccc1 = corr(squeeze(coverage_all(ii_exp,:,1,:))');
    ccc2 = corr(squeeze(coverage_all(ii_exp,:,2,:))');
    ccc_all(ii_exp,1) = ccc1(2);
    ccc_all(ii_exp,2) = ccc2(2);
end

%% save processed data
dir_OUT = 'E:\Tamir\work\PROJECTS\LargeScale\paper_replay\data_prepared_for_figures';
filename = 'replay_coverage';
file_OUT = fullfile(dir_OUT,filename);
save(file_OUT,'T','coverage_all','ccc_all','n_seqs_all');

%% plot coverage
xbins = linspace(0,1,nbins);
fig=figure;
fig.Units = 'centimeters';
fig.Position([3 4]) = [15 6];
axes('Units','normalized','Position',[0.1 0.2 0.85 0.65])
hold on
for ii_exp = 1:height(T)
    exp_ID = T.exp_ID{ii_exp};
    c = squeeze(coverage_all(ii_exp,:,:,:));
    cla
    lw = 2;
    plot(xbins, squeeze(c(1,1,:)),'-b','LineWidth',lw)
    plot(xbins, squeeze(c(2,1,:)),'--b','LineWidth',lw)
    plot(xbins, squeeze(c(1,2,:)),'-r','LineWidth',lw)
    plot(xbins, squeeze(c(2,2,:)),'--r','LineWidth',lw)
    text(0.05,1.05,sprintf('#%d %s',ii_exp,exp_ID),'Interpreter','none','HorizontalAlignment','left','Units','normalized');
    xlabel('Position (norm.)')
    ylabel('Replay coverage (counts)')
    hl=legend({'sleep dir 1','rest dir 1','sleep dir 2','rest dir 2'},'NumColumns',2);
    hl.Units = 'normalized';
    hl.Position([1 2]) = [0.6 0.85];
    dir_OUT = 'F:\sequences\decoded_figs\replay_coverage';
    filename = sprintf('ii_%.3d_replay_coverage_%s',ii_exp,exp_ID);
    file_OUT = fullfile(dir_OUT,filename);
    saveas(fig,file_OUT,'jpg')
end


%%
function [coverage, n_seqs] = calc_coverage(exp_ID, epoch_type, params_opt,nbins)

coverage = zeros(2,nbins);
n_seqs = zeros(2,1);

[events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
if isempty(events)
    return;
end
seqs = [events.seq_model];
[seqs, TF] = decoding_apply_seq_inclusion_criteria(seqs);
events(~TF)=[];
if isempty(seqs)
    return;
end
seqs_edges = [seqs.start_pos_norm; seqs.end_pos_norm];
seqs_edges = [min(seqs_edges ); max(seqs_edges )]';

xbins = linspace(0,1,nbins);
directions = [1 -1];
for ii_dir = 1:2
    direction = directions(ii_dir);
    coverage(ii_dir,:) = sum( xbins>seqs_edges([seqs.state_direction]==direction,1) & ...
                              xbins<seqs_edges([seqs.state_direction]==direction,2),1);
    n_seqs(ii_dir) = sum([seqs.state_direction]==direction);
end



end