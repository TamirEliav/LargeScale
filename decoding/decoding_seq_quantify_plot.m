function decoding_seq_quantify_plot(exp_ID, epoch_type, params_opt, event_type)
arguments
    %% 
    exp_ID = 'b9861_d180526'
    epoch_type {mustBeMember(epoch_type,{'sleep','rest','flight'})} = 'sleep'
    params_opt = 11;
    event_type {mustBeMember(event_type,{'PE','posterior','ripples','MUA'})} = 'posterior'
end

%% IN/OUT folders
dir_OUT = 'F:\sequences\events_quantification\scatters';
mkdir(dir_OUT);

%% load data
[events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
if length(events)==0
    fprintf('No event for %s %s dec param %d %s\n', exp_ID, epoch_type, params_opt, event_type);
    return;
end
seqs = [events.seq];

%% plot
close all
figs={}; filters={};
[figs{end+1},filters{end+1}]= plot_features_scatters(seqs);
[figs{end+1},filters{end+1}]= plot_features_scatters(seqs,'duration','>',0.05,'distance','>',5);
[figs{end+1},filters{end+1}]= plot_features_scatters(seqs,'duration','>',0.05,'distance','>',5,'confidence_KS','>',0.5);
[figs{end+1},filters{end+1}]= plot_features_scatters(seqs,'duration','>',0.05,'distance','>',5,'confidence_KS','>',0.6);
[figs{end+1},filters{end+1}]= plot_features_scatters(seqs,'duration','>',0.05,'distance','>',5,'confidence_KS','>',0.7);
[figs{end+1},filters{end+1}]= plot_features_scatters(seqs,'duration','>',0.05,'distance','>',5,'confidence_KS','>',0.8);
[figs{end+1},filters{end+1}]= plot_features_scatters(seqs,'duration','>',0.05,'distance','>',5,'confidence_sparsity','>',0.5);
[figs{end+1},filters{end+1}]= plot_features_scatters(seqs,'duration','>',0.05,'distance','>',5,'confidence_sparsity','>',0.6);
[figs{end+1},filters{end+1}]= plot_features_scatters(seqs,'duration','>',0.05,'distance','>',5,'confidence_sparsity','>',0.7);
[figs{end+1},filters{end+1}]= plot_features_scatters(seqs,'duration','>',0.05,'distance','>',5,'confidence_sparsity','>',0.8);
[figs{end+1},filters{end+1}]= plot_features_scatters(seqs,'duration','>',0.05,'distance','>',5,'score','>',0.2);
[figs{end+1},filters{end+1}]= plot_features_scatters(seqs,'duration','>',0.05,'distance','>',5,'score','>',0.3);
[figs{end+1},filters{end+1}]= plot_features_scatters(seqs,'duration','>',0.05,'distance','>',5,'score','>',0.4);
[figs{end+1},filters{end+1}]= plot_features_scatters(seqs,'duration','>',0.05,'distance','>',5,'score','>',0.5);
[figs{end+1},filters{end+1}]= plot_features_scatters(seqs,'duration','>',0.05,'distance','>',5,'score','>',0.6);
[figs{end+1},filters{end+1}]= plot_features_scatters(seqs,'duration','>',0.05,'distance','>',5,'score','>',0.7);

title_str = sprintf('%s - %s - x%d - %s',exp_ID, epoch_type, params.decode.replay_speed, event_type);
for ii_fig = 1:length(figs)
    fig = figs{ii_fig};
    filter = filters{ii_fig};
    figure(fig);
    h=sgtitle(title_str); h.Interpreter = 'none'; %h.Position(2) = -0.025;
    filters_str = arrayfun(@(x)(sprintf('%s_%g',x.feature,x.val)), filter, 'UniformOutput', false)';
    filters_str = strjoin(filters_str,'_');
    filename = fullfile(dir_OUT, sprintf('%s_events_%s_dec_prm_%d_%s_%s',exp_ID,epoch_type,params_opt,event_type,filters_str));
    saveas(fig, [filename '.jpg'], 'jpg');
end
% close all

end


%%
function [fig, filter] = plot_features_scatters(seqs,flt_features,flt_funcs,flt_vals)
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

    %% plot scatter matrix figure
    T = struct2table(seqs);
    T = T(:,["distance","duration","start_pos","end_pos","score","compression","confidence_KS","confidence_sparsity"]);
    fig = figure;
    fig.WindowState = 'maximized';
    if ~isempty(T)
        [h,ax,bigax]=gplotmatrix(table2array(T),[],[],[],[],[],[],[],T.Properties.VariableNames);
        set([ax.XLabel ax.YLabel],'Interpreter','none');
    else
        text(.5,.5,'No data points');
    end
    
    %% add filters details text
    filter_str = arrayfun(@(x)([x.feature x.func num2str(x.val)]), filter, 'UniformOutput', false)';
    if isempty(filter_str)
        filter_str = 'no filter';
    end
    annotation('textbox', [0.8 0.95 0.2 0.05], 'String',filter_str,...
        'LineStyle','None','Units','normalized','HorizontalAlignment','left','Interpreter','none');
end


%%
