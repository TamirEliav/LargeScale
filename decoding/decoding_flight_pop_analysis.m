function decoding_flight_pop_analysis(exp_ID_list,params_opt)
arguments
    %%
    exp_ID_list = {'b0184_d191128';'b0184_d191129';'b0184_d191130';'b0184_d191201';'b0184_d191202';'b0184_d191203';'b0184_d191204';}
    params_opt = 4
end

%% folders
dir_OUT = 'F:\sequences\decoded_figs\flight\population'
mkdir(dir_OUT)

%% load data
epoch_type = 'flight';
dir_IN = 'F:\sequences\decoded_figs\flight\conf_mat';
res_all = {};
for ii_exp = 1:length(exp_ID_list)
    exp_ID = exp_ID_list{ii_exp};
    exp = exp_load_data(exp_ID,'details');
    filename = fullfile(dir_IN, sprintf('%s_%s_decoding_opt_%d',exp_ID,epoch_type,params_opt) );
    clear res
    load(filename);
    res.exp_ID = exp_ID;
    res.sparsity_max = max(res.sparsity);
    res.details = exp.details;
    res_all{ii_exp} = res;
end
res_all = [res_all{:}];
% res_all = [res_all.res];

%% arrange data
num_TT_used = cellfun(@(x)(sum(contains(x,'CA'))), {details.TT_loc});
[res_all.num_TT_used] = disperse(num_TT_used);

%%
fig=figure;
hax=gca;
clr = hax.ColorOrder;
close(fig)

%% plot scatters
fig=figure;
fig.WindowState = 'maximized';
X = [];
X(:,end+1) = [res_all.pos_err_mean];
X(:,end+1) = [res_all.pos_err_median];
X(:,end+1) = [res_all.sparsity_mean];
X(:,end+1) = [res_all.sparsity_max];
X(:,end+1) = [res_all.num_TT_used];
labels = {'pos err (mean) [m]';'pos err (median) [m]';'sparsity (mean)';'sparsity (max)';'numTT'};
details = [res_all.details];
G = [details.batNum];
[g,gN] = grp2idx(G);
h=gplotmatrix(X,[],G,clr,[],[],[],[],labels);
clear hlinks
for ii_grp = 1:size(h,3)
    % get relevant plot handles
    hl_grp = findobj(h(:,:,ii_grp),'type','line');
    % link brush
    hlinks(ii_grp) = linkprop(hl_grp,'BrushData');
    % add exp_ID details
    for ii_h = 1:length(hl_grp)
        hl_grp(ii_h).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('exp ID',exp_ID_list(g==ii_grp));
    end
end
setappdata(fig,'LinkBrushData',hlinks);
brush on
sgtitle('Flight decoding populaiton analysis')
filename = fullfile(dir_OUT, 'flight_decoding_population_gmatplot');
saveas(fig, filename , 'fig');
saveas(fig, filename , 'jpg');

%% plot
fig=figure;
fig.WindowState = 'maximized';
hold on
x = [res_all.pos_err_median]';
y = [res_all.sparsity_max]';
% h = plot(x, y, 'o');
% h = scatter(x,y,[],g,'filled');
h = splitapply(@(x,y,exp_ID)(plot(x,y,'.','MarkerSize',15,'UserData',exp_ID)), x,y,exp_ID_list,g);
for ii_grp = 1:length(h)
    h(ii_grp).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('exp ID','UserData');
end
legend("bat "+gN)
xlabel('Median position error (m)')
ylabel('Max predicted sparsity')
sgtitle('Flight decoding populaiton analysis')
filename = fullfile(dir_OUT, 'flight_decoding_population_sparsity_vs_error');
saveas(fig, filename , 'fig');
saveas(fig, filename , 'jpg');

end