function decoding_flight_pop_analysis(exp_list,params_opt)
arguments
    %%
    exp_list
    params_opt = 4
end

%% folders & params
prm = PARAMS_GetAll();
dir_OUT = 'F:\sequences\decoded_figs\flight\population';
mkdir(dir_OUT);

%% arrange exp list
exp_t = DS_get_exp_summary();
exp_t(~contains(exp_t.exp_ID, exp_list),:) = [];
exp_t = add_manual_inclusion(exp_t);
    
%% load data
epoch_type = 'flight';
dir_IN = 'F:\sequences\decoded_figs\flight\conf_mat';
res_all = {};
for ii_exp = 1:height(exp_t)
    exp_ID = exp_t.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'details');
    filename = fullfile(dir_IN, sprintf('%s_%s_decoding_opt_%d',exp_ID,epoch_type,params_opt) );
    clear res res_raw
    load(filename,'res');
    res.exp_ID = exp_ID;
    res.sparsity_max = max(res.sparsity);
    res.details = exp.details;
    res_all{ii_exp} = res;
end
res_all = [res_all{:}];

%% arrange data
details = [res_all.details];
num_TT_used = cellfun(@(x)(sum(contains(x,'CA'))), {details.TT_loc});
[res_all.num_TT_used] = disperse(num_TT_used);
[res_all.pos_err_median_prc] = disperse([res_all.pos_err_median_prc].*100);

%% arrange results in a matrix
X = [];
X(:,end+1) = [res_all.pos_err_mean];
X(:,end+1) = [res_all.pos_err_median];
X(:,end+1) = [res_all.pos_err_median_prc];
X(:,end+1) = [res_all.mean_err_prob];
X(:,end+1) = [res_all.err_prob_by_predicted_max];
X(:,end+1) = [res_all.num_TT_used];
labels = {
    'pos err (mean) [m]';
    'pos err (median) [m]';
    'pos err (median) [%]';
    'mean error prob';
    'max predicted error prob';
    'numTT'};
G = [details.batNum];
[g,gN] = grp2idx(G);

%% arrange results in a table
T = table(...
    exp_t.bat_num,...
    [res_all.pos_err_mean]',...
    [res_all.pos_err_median]',...
    [res_all.pos_err_mean_prc]',...
    [res_all.pos_err_median_prc]',...
    [res_all.err_prob_by_predicted_max]',...
    [res_all.mean_err_prob]',...
    [res_all.num_TT_used]',...
    exp_t.included_manual,...
    'VariableNames',...
    {'bat_num','pos_err_mean_m','pos_err_median_m','pos_err_mean_prc','pos_err_median_prc','err_prob_by_predicted_max','mean_err_prob','num_TT_used','included_manual'},...
    'RowNames', exp_t.exp_ID);
T.Properties.DimensionNames{1} = 'exp_ID';

%% add auto inclusion criteria 
TF = true(size(T,1),1);
TF = TF & T.mean_err_prob < prm.decode.inc_criteria.err_prob;
TF = TF & T.pos_err_median_prc  < prm.decode.inc_criteria.median_err_prc;
% TF = TF & T.err_prob_by_predicted_max < prm.decode.inc_criteria.max_predict_err_prob;
T.included_auto = TF;
filename = 'F:\sequences\inclusion\exp_inc_list.xlsx';
writetable(T,filename,'WriteVariableNames',true,'WriteRowNames',true);

%% colors (TODO: get better colors from a pallete)
switch 2
    case 1
        fig=figure;
        hax=gca;
        clr = hax.ColorOrder;
        close(fig)
    case 2
        clr = brewermap(12,'set1');
end

%% plot scatters
fig=figure;
fig.WindowState = 'maximized';
hax=gca;
hax.ColorOrder = clr;
h=gplotmatrix(X,[],G,clr,[],[],[],[],labels);
clear hlinks
for ii_grp = 1:size(h,3)
    % get relevant plot handles
    hl_grp = findobj(h(:,:,ii_grp),'type','line');
    % link brush
    hlinks(ii_grp) = linkprop(hl_grp,'BrushData');
    % add exp_ID details
    for ii_h = 1:length(hl_grp)
        hl_grp(ii_h).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('exp ID',exp_list(g==ii_grp));
    end
end
setappdata(fig,'LinkBrushData',hlinks);
brush on
sgtitle('Flight decoding populaiton analysis')
filename = fullfile(dir_OUT, 'flight_decoding_pop_scatters_all');
saveas(fig, filename , 'fig');
saveas(fig, filename , 'jpg');

%% max predicted error probability vs. median error
fig=figure;
fig.WindowState = 'maximized';
hax=gca;
hax.ColorOrder = clr;
hold on
x = [res_all.pos_err_median_prc]';
y = [res_all.err_prob_by_predicted_max]';
h = splitapply(@(x,y,exp_ID)(plot(x,y,'.','MarkerSize',15,'UserData',exp_ID)), x,y,exp_list,g);
for ii_grp = 1:length(h)
    h(ii_grp).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('exp ID','UserData');
end
xline(1)
hax=gca;
hax.XScale = 'log';
hax.YScale = 'log';
legend("bat "+gN);
xlabel('Median position error (normalized to environment size) [%]')
ylabel('Max predicted error probability')
sgtitle('Flight decoding populaiton analysis')
filename = fullfile(dir_OUT, 'flight_decoding_pop_scatter_max_predict_err_prob_vs_median_error');
saveas(fig, filename , 'fig');
saveas(fig, filename , 'jpg');

%% max predicted error probability vs. mean error rate
fig=figure;
fig.WindowState = 'maximized';
hax=gca;
hax.ColorOrder = clr;
hold on
x = T.mean_err_prob;
y = T.err_prob_by_predicted_max;
TF = T.included_manual;
bats_numbers = unique(T.bat_num);
for ii_bat = 1:length(bats_numbers)
    bat_num = bats_numbers(ii_bat);
    IX1 = T.bat_num==bat_num;
    h1 = plot(x(IX1),y(IX1),'o', 'MarkerSize',8,'UserData',T.exp_ID(IX1),'LineWidth',1.5,'DisplayName',"bat "+bat_num);
    h1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('exp ID','UserData');
end
for ii_bat = 1:length(bats_numbers)
    bat_num = bats_numbers(ii_bat);
    IX2 = T.bat_num==bat_num & TF==1;
    if ~any(IX2)
        continue;
    end
    h2 = plot(x(IX2),y(IX2),'*k', 'MarkerSize',5,'UserData',T.exp_ID(IX2),'DisplayName','yes');
    h2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('exp ID','UserData');
end
for ii_bat = 1:length(bats_numbers)
    bat_num = bats_numbers(ii_bat);
    IX3 = T.bat_num==bat_num & TF==0.5;
    if ~any(IX3)
        continue;
    end
    h2 = plot(x(IX3),y(IX3),'+k', 'MarkerSize',5,'UserData',T.exp_ID(IX3),'DisplayName','maybe');
    h2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('exp ID','UserData');
end

% h1 = splitapply(@(x,y,exp_ID)(plot(x,y,'o', 'MarkerSize',5,'UserData',exp_ID,'LineWidth',1.8)), x,    y,    exp_list,    g);
% h2 = splitapply(@(x,y,exp_ID)(plot(x,y,'xk','MarkerSize',4,'UserData',exp_ID)),                 x(TF),y(TF),exp_list(TF),g(TF));
% for ii_grp = 1:length(h)
%     h1(ii_grp).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('exp ID','UserData');
%     h2(ii_grp).DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('exp ID','UserData');
% end

xline(prm.decode.inc_criteria.err_prob)
yline(prm.decode.inc_criteria.max_predict_err_prob)
hax=gca;
hax.XScale = 'log';
hax.YScale = 'log';
% legend("bat "+bats_numbers,'Location','northwest');
legend('Location','northwest')
xlabel('Error probability')
ylabel('Max predicted error probability')
sgtitle('Flight decoding populaiton analysis')
filename = fullfile(dir_OUT, 'flight_decoding_pop_scatter_max_predict_err_prob_vs_error_prob');
saveas(fig, filename , 'fig');
saveas(fig, filename , 'jpg');


%% median error vs. mean error rate
fig=figure;
fig.WindowState = 'maximized';
hax=gca;
hax.ColorOrder = clr;
hold on
x = T.mean_err_prob;
y = [res_all.pos_err_median_prc]';
TF = T.included_manual;
bats_numbers = unique(T.bat_num);
for ii_bat = 1:length(bats_numbers)
    bat_num = bats_numbers(ii_bat);
    IX1 = T.bat_num==bat_num;
    h1 = plot(x(IX1),y(IX1),'o', 'MarkerSize',8,'UserData',T.exp_ID(IX1),'LineWidth',1.5,'DisplayName',"bat "+bat_num);
    h1.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('exp ID','UserData');
end
for ii_bat = 1:length(bats_numbers)
    bat_num = bats_numbers(ii_bat);
    IX2 = T.bat_num==bat_num & TF==1;
    if ~any(IX2)
        continue;
    end
    h2 = plot(x(IX2),y(IX2),'*k', 'MarkerSize',5,'UserData',T.exp_ID(IX2),'DisplayName','yes');
    h2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('exp ID','UserData');
end
for ii_bat = 1:length(bats_numbers)
    bat_num = bats_numbers(ii_bat);
    IX3 = T.bat_num==bat_num & TF==0.5;
    if ~any(IX3)
        continue;
    end
    h2 = plot(x(IX3),y(IX3),'+k', 'MarkerSize',5,'UserData',T.exp_ID(IX3),'DisplayName','maybe');
    h2.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('exp ID','UserData');
end

xline(prm.decode.inc_criteria.err_prob);
yline(prm.decode.inc_criteria.median_err_prc);
hax=gca;
hax.XScale = 'log';
hax.YScale = 'log';
% legend("bat "+bats_numbers,'Location','northwest');
legend('Location','northwest')
xlabel('Error probability')
ylabel('Median position error (normalized to environment size) [%]')
sgtitle('Flight decoding populaiton analysis')
filename = fullfile(dir_OUT, 'flight_decoding_pop_scatter_median_err_vs_error_prob');
saveas(fig, filename , 'fig');
saveas(fig, filename , 'jpg');


end


%%
function exp_t = add_manual_inclusion(exp_t)
%%
inc_list_filename = "L:\Analysis\Code\inclusion_lists\decoding_exp_list.xlsx";
T = readtable(inc_list_filename, 'ReadRowNames',true,'ReadVariableNames',true);
exp_t = innerjoin(exp_t,T,'keys','exp_ID');

end




%%
