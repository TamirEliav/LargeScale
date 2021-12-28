function decoding_compare_replay_speeds(exp_ID, epoch_type, params_opts, event_type)
arguments
    %% 
    exp_ID = 'b9861_d180526'
    epoch_type {mustBeMember(epoch_type,{'sleep','rest','flight'})} = 'sleep'
    params_opts = [8:14]
    event_type {mustBeMember(event_type,{'PE','posterior','ripples','MUA'})} = 'posterior'
end

%% IN/OUT folders
dir_OUT = 'F:\sequences\events_quantification\compare_replay_speeds';
mkdir(dir_OUT);

%% load data
exp = exp_load_data(exp_ID, 'details','path');
events_per_opt={};
params_per_opt={};
for ii_opt = 1:length(params_opts)
    params_opt = params_opts(ii_opt);
    [events_per_opt{ii_opt} params_per_opt{ii_opt}] =...
        decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
end
params_per_opt = [params_per_opt{:}];
dec_params_per_opt = [params_per_opt.decode];
replay_speeds = [dec_params_per_opt.replay_speed];

%% arrange sequences quantification struct array
events_all = [events_per_opt{:}];
nEvents = cumsum(cellfun(@length, events_per_opt));
grp = interp1(nEvents,1:length(nEvents),1:nEvents(end),'next','extrap');
seqs = [events_all.seq_model];
[seqs.grp] = disperse(grp);
[seqs.mean_prob] = events_all.mean_prob;
[seqs.confidence_sparsity] = disperse([seqs.confidence_sparsity] .* [seqs.mean_prob]);
[seqs.replay_speed] = disperse(replay_speeds(grp));

%%
close all
figs={};
filters={};
[figs{end+1},filters{end+1}]= plot_features_vs_replay_speed(seqs);
[figs{end+1},filters{end+1}]= plot_features_vs_replay_speed(seqs,'duration','>',0.05,'distance','>',5);
[figs{end+1},filters{end+1}]= plot_features_vs_replay_speed(seqs,'duration','>',0.05,'distance','>',5,'confidence_KS','>',0.5);
[figs{end+1},filters{end+1}]= plot_features_vs_replay_speed(seqs,'duration','>',0.05,'distance','>',5,'confidence_KS','>',0.6);
[figs{end+1},filters{end+1}]= plot_features_vs_replay_speed(seqs,'duration','>',0.05,'distance','>',5,'confidence_KS','>',0.7);
[figs{end+1},filters{end+1}]= plot_features_vs_replay_speed(seqs,'duration','>',0.05,'distance','>',5,'confidence_KS','>',0.8);
[figs{end+1},filters{end+1}]= plot_features_vs_replay_speed(seqs,'duration','>',0.05,'distance','>',5,'confidence_sparsity','>',0.5);
[figs{end+1},filters{end+1}]= plot_features_vs_replay_speed(seqs,'duration','>',0.05,'distance','>',5,'confidence_sparsity','>',0.6);
[figs{end+1},filters{end+1}]= plot_features_vs_replay_speed(seqs,'duration','>',0.05,'distance','>',5,'confidence_sparsity','>',0.7);
[figs{end+1},filters{end+1}]= plot_features_vs_replay_speed(seqs,'duration','>',0.05,'distance','>',5,'confidence_sparsity','>',0.8);
[figs{end+1},filters{end+1}]= plot_features_vs_replay_speed(seqs,'duration','>',0.05,'distance','>',5,'score','>',0.2);
[figs{end+1},filters{end+1}]= plot_features_vs_replay_speed(seqs,'duration','>',0.05,'distance','>',5,'score','>',0.3);
[figs{end+1},filters{end+1}]= plot_features_vs_replay_speed(seqs,'duration','>',0.05,'distance','>',5,'score','>',0.4);
[figs{end+1},filters{end+1}]= plot_features_vs_replay_speed(seqs,'duration','>',0.05,'distance','>',5,'score','>',0.5);
[figs{end+1},filters{end+1}]= plot_features_vs_replay_speed(seqs,'duration','>',0.05,'distance','>',5,'score','>',0.6);
[figs{end+1},filters{end+1}]= plot_features_vs_replay_speed(seqs,'duration','>',0.05,'distance','>',5,'score','>',0.7);

compare_fig_title = {sprintf('%s - %s - %s',exp_ID,epoch_type,event_type);
                    'comparing sequence features across replay speeds'};
for ii_fig = 1:length(figs)
    fig = figs{ii_fig};
    filter = filters{ii_fig};
    figure(fig);
    h=sgtitle(compare_fig_title); h.Interpreter = 'none'; %h.Position(2) = -0.025;
    filters_str = arrayfun(@(x)(sprintf('%s_%g',x.feature,x.val)), filter, 'UniformOutput', false)';
    filters_str = strjoin(filters_str,'_');
    filename = fullfile(dir_OUT, sprintf('%s_compare_replay_speeds_%s_%s_%s',exp_ID,epoch_type,event_type,filters_str) );
    saveas(fig, [filename '.jpg'], 'jpg');
end
% close all

%% "specificity" plots
% % % % forward_lbl_map = containers.Map({true,false},{'Forward','Reverse'});
% % % % direction_lbl_map = containers.Map({-1,1},{'Inbound','Outbound'});
% % % % features = {'score';'compression';'distance';'confidence_KS';};
% % % % features_label = {'Replay score';'Compression';'Distance (m)';'Confidence (KS)';};
% % % % for ii_feature = 1:length(features)
% % % %     fig = figure;
% % % %     fig.WindowState = 'maximized';
% % % %     ht = tiledlayout(3,length(valid_grps),'TileSpacing','compact');
% % % %     feature = features{ii_feature};
% % % %     for ii_opt = valid_grps
% % % %         nexttile
% % % %         opt_events = events_per_opt{ii_opt};
% % % %         opt_seqs = [opt_events.seq];
% % % %         vals = [opt_seqs.(feature)];
% % % %         grps = [opt_seqs.forward];
% % % %         violinplot(vals,grps);
% % % %         [~,ID] = findgroups(grps);
% % % %         xticklabels(arrayfun(@(x)(forward_lbl_map(x)), ID, 'UniformOutput', false));
% % % %         ylabel(features_label(ii_feature));
% % % %         title("x"+replay_speeds(ii_opt));
% % % %     end
% % % %     for ii_opt = valid_grps
% % % %         nexttile
% % % %         opt_events = events_per_opt{ii_opt};
% % % %         opt_seqs = [opt_events.seq];
% % % %         vals = [opt_seqs.(feature)];
% % % %         grps = [opt_seqs.direction];
% % % %         violinplot(vals,grps);
% % % %         [~,ID] = findgroups(grps);
% % % %         xticklabels(arrayfun(@(x)(direction_lbl_map(x)), ID, 'UniformOutput', false));
% % % %         xlabel('Direction');
% % % %         ylabel(features_label(ii_feature));
% % % %         title("x"+replay_speeds(ii_opt));
% % % %     end
% % % %     for ii_opt = valid_grps
% % % %         nexttile
% % % %         opt_events = events_per_opt{ii_opt};
% % % %         opt_seqs = [opt_events.seq];
% % % %         vals = [opt_seqs.(feature)];
% % % %         grps = [opt_seqs.state_direction];
% % % %         violinplot(vals,grps);
% % % %         [~,ID] = findgroups(grps);
% % % %         xticklabels(arrayfun(@(x)(direction_lbl_map(x)), ID, 'UniformOutput', false));
% % % %         xlabel('State direction');
% % % %         ylabel(features_label(ii_feature));
% % % %         title("x"+replay_speeds(ii_opt));
% % % %     end
% % % %     h=suptitle({exp_ID; sprintf('comparing sequence %s across replay speeds',feature)});
% % % %     h.Interpreter = 'none';
% % % %     h.Position(2) = -0.025;
% % % %     dirname = fullfile(dir_OUT, 'compare_directions');
% % % %     mkdir(dirname);
% % % %     filename = fullfile(dirname,...
% % % %         [sprintf('%s_compare_by_directions_%s_%s_dec_prms',exp_ID,epoch_type,event_type) ...
% % % %          sprintf('_%d',params_opts) ...
% % % %          '_' feature]);
% % % %     saveas(gcf, filename, 'jpg');
% % % % end

end




%%
function [fig, filter] = plot_features_vs_replay_speed(seqs,flt_features,flt_funcs,flt_vals)
arguments 
    seqs
end
arguments (Repeating)
    flt_features
    flt_funcs
    flt_vals
end
    %% filter by features (optional)
    filter = repelem(struct('feature',[],'func',[],'val',[]),0);
    valid = true(size(seqs));
    for ii_feature = 1:length(flt_features)
        feature = flt_features{ii_feature};
        func = flt_funcs{ii_feature};
        val = flt_vals{ii_feature};
        filter(ii_feature).feature = feature;
        filter(ii_feature).func = func;
        filter(ii_feature).val = val;
        if ~isa(func,'function_handle')
            func = str2func(func);
        end
        valid = valid & func([seqs.(feature)],val);
    end
    seqs(~valid)=[];
    
    %%
    features2plot = {'duration','speed','distance','start_pos','end_pos','middle_pos','direction','state_direction','forward','compression','confidence_KS','confidence_HPD','confidence_sparsity','score'};
    features_label = {'Duration (s)','Speed (m/s)','Distance (m)','Start pos (m)','End pos (m)','Middle pos (m)','Direction','State direction','Forward','Compression','Confidence (KS)','Confidence (HPD)','Confidence (sparsity)','score'};
    fig = figure;
    fig.WindowState = 'maximized';
    ht = tiledlayout('flow','TileSpacing','compact');
    ht.XLabel.String = 'Replay speed';
    for ii_feature = 1:length(features2plot)
        nexttile
        feature = features2plot{ii_feature};
        grps = [seqs.grp];
        vals = [seqs.(feature)];
        h=violinplot(double(vals),grps); % TODO: do not use violin plot for boolean/discrete variables
        [h.ShowMean] = disperse(repelem(1,length(h)));
        set([h.MeanPlot],'LineWidth',2);
        set([h.MeanPlot],'Color','k');
        xticklabels("x"+unique([seqs.replay_speed]));
        ylabel(features_label(ii_feature));
    end
    
    %%
    nexttile
    hold on
    [G,ID] = findgroups([seqs.replay_speed]);
    median_compression = splitapply(@median, [seqs.compression], G);
    mean_compression = splitapply(@mean, [seqs.compression], G);
    plot(ID,median_compression,'o-');
    plot(ID,mean_compression,'o-');
    h=refline(1,0);h.Color = 'k';
    legend('median','mean','Location','northwest');
    ylabel('Compression')

    %% add filters details text
    filter_str = arrayfun(@(x)([x.feature x.func num2str(x.val)]), filter, 'UniformOutput', false)';
    if isempty(filter_str)
        filter_str = 'no filter';
    end
    annotation('textbox', [0.8 0.95 0.2 0.05], 'String',filter_str,...
        'LineStyle','None','Units','normalized','HorizontalAlignment','left','Interpreter','none');
    
end




%%

