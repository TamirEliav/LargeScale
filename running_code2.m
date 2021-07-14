%%
exps_IDs = exp_t.exp_ID;

%%
exps = cellfun(@(exp_ID)( exp_load_data(exp_ID,'details','PE') ), exps_IDs);
PEs = {exps.PE};
PE_all = [exps.PE];
nPE = cellfun(@length, PEs);
mean_FR  = cellfun(@(PE)(mean([PE.peak_FR])), PEs);
mean_zFR  = cellfun(@(PE)(mean([PE.peak_zFR])), PEs);
order_feature = nPE .* mean_zFR;
[~,ranks]  = ismember(order_feature,flip(unique(order_feature)));
% [~,sort_IX] = sort(order_feature,'descend');
[exps.nPE] = disperse(nPE);
[exps.mean_zFR] = disperse(mean_zFR);
[exps.order_feature] = disperse(order_feature);
[exps.rank] = disperse(ranks);

%%
figure
tiledlayout('flow')
nexttile
h=scatter(nPE,mean_FR,20,order_feature,'filled');
h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('ranks',ranks);
h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('feature value',order_feature);
% lsline
xlabel('No. of events')
ylabel('Mean firing rate (Hz)')
colorbar

nexttile
h=scatter(nPE,mean_zFR,20,order_feature,'filled');
h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('ranks',ranks);
h.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('feature value',order_feature);
% lsline
xlabel('No. of events')
ylabel('Mean firing rate (z)')
colorbar

%% sort exps
% exps=exps(sort_IX);

T1 = table( arrayfun(@(exp)(exp.details.exp_ID), exps,'UniformOutput',0),...
            arrayfun(@(exp)(exp.details.batNum), exps,'UniformOutput',1),...
            datetime(arrayfun(@(exp)(exp.details.date), exps,'UniformOutput',1),'Format','yyMMdd'),...
            arrayfun(@(exp)(exp.nPE), exps,'UniformOutput',1),...
            arrayfun(@(exp)(exp.mean_zFR), exps,'UniformOutput',1),...
            arrayfun(@(exp)(exp.order_feature), exps,'UniformOutput',1),...
            arrayfun(@(exp)(exp.rank), exps,'UniformOutput',1),...
            'VariableNames',{'exp_ID','batNum','date','nPE','mean_zFR','order_feature','rank'});
T1.Properties.RowNames = string(datestr(T1.date,'yymmdd'));

%%
cells_table_filename = "L:\Analysis\Code\inclusion_lists\cells_summary.xlsx";
cells_t = readtable(cells_table_filename);
T2 = groupcounts(cells_t,'date');
T2.Properties.VariableNames(2) = {'nCells'};
T2.date = datetime(T2.date,'Format','yyMMdd');
T2.Properties.RowNames = string(datestr(T2.date,'yymmdd'));

%%
T = outerjoin(T2,T1,'Keys','Row');
T = sortrows(T,'nCells','descend');
clc
T

%%








%%
