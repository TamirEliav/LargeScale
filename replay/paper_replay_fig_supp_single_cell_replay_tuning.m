%% Replay - Fig supp XXX - single cell replay tuning
clear 
clc
close all

%% data options 
cat_IX = 3; % sleep+rest pooled

%% plotting options

%% graphics params

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Extended_Data_Fig_XXX_single cell_replay_tuning';
fig_caption_str = 'single_cell_replay_tuning';
log_name_str = [fig_name_str '_log_file' '.txt'];
log_name_str = strrep(log_name_str , ':', '-');
log_name_str = strrep(log_name_str , ' ', '_');
log_name_out = fullfile(res_dir, log_name_str);

%% open log file
diary off
diary(log_name_out)
diary on
disp('Log file');
disp(['created: ', datestr(clock)]);
disp('======================================================');
disp([fig_name_str ':' fig_caption_str]);   
disp('======================================================');
disp('');

%% create figure
% figure_size_cm = [21.0 29.7]; % ~A4
figure_size_cm = [21.6 27.9]; % ~US letter
fig = figure;
% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',7);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'DefaultAxesUnits','centimeters');
set(gcf,'PaperType','usletter')
% set(gcf,'PaperType','<custom>');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 figure_size_cm]);
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]); % position on screen...
set(gcf, 'Renderer', 'painters');
% set(gcf, 'color', 'none');
set(groot, 'defaultAxesColor','None')
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');

% create panels
panels{1}(1) = axes('position', [3 20 4 4]);

%% load data - get relevant cells
prm=PARAMS_GetAll();
cells_t = DS_get_cells_summary();
cells_t(~ismember(cells_t.bat, [184,194,2382,148,34,9861,2289] ),:) = [];
[exp_list,T] = decoding_get_inclusion_list();
cells_exp_ID = cellfun(@(c)DS_get_exp_ID_from_cell_ID(c),cells_t.cell_ID,'UniformOutput',false);
TF = ismember(cells_exp_ID,exp_list);
cells_t(~TF,:)=[];
cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',1);
details = [cells.details];
details(~contains({details.brain_area}, {'CA1','CA3'})) = [];
details(~ismember([details.ClusterQuality], [2])) = [];
cells_t = cells_t({details.cell_ID},:);
% cells_t('b2382_d190623_TT13_SS01',:) = []; % not stable - remove this cell

%% load data - get cells results
cells = cellfun(@(c)cell_load_data(c,'details','FR_map','replay_FR_map','signif','inclusion'),cells_t.cell_ID);

%% calc shuffle corr
% corr_type = 'Pearson';
corr_type = 'Spearman';
details = [cells.details];
cells_bats_num = [details.bat];
cells_inclusion = cat(1,cells.inclusion);
cells_signif = cat(1,cells.signif);
FR_maps = cat(1,cells.FR_map);
FR_maps = arrayfun(@(x)x.all,FR_maps);
PSTH_flight = [];
for ii_cell = 1:length(cells)
    for ii_dir=1:2
        PSTH_flight(ii_cell,ii_dir,:) = FR_maps(ii_cell,ii_dir).PSTH;
    end
end
replay_FR_maps = cat(1,cells.replay_FR_map);
replay_FR_maps = cat(3,replay_FR_maps.replay_PSTH_all);
PSTH_replay = [];
for ii_cell = 1:length(cells)
    for ii_dir=1:2
        PSTH_replay(ii_cell,ii_dir,:) = replay_FR_maps(cat_IX,ii_dir,ii_cell).PSTH;
    end
end
ccc_shuffle = {};
ccc_data = {};
for ii_dir = 1:2
    cells_inc_TF = true(size(cells))';
    cells_inc_TF = cells_inc_TF & [replay_FR_maps(cat_IX,ii_dir,:).valid_fields_coverage_prc]>0.5;
    cells_inc_TF = cells_inc_TF & [cells_inclusion(:,ii_dir).TF];
    cells_inc_TF = cells_inc_TF & [cells_inclusion(:,ii_dir).pyr];
    cells_inc_TF = cells_inc_TF & [cells_signif(:,ii_dir).SI_thr_shuffle];
    cells_inc_TF = cells_inc_TF & [cells_signif(:,ii_dir).SI_thr_signif];
%     cells_inc_TF = cells_inc_TF & [cells_signif(:,ii_dir).TF]; % this excludes maps without signif fields
    
    X = squeeze(PSTH_flight(cells_inc_TF,ii_dir,:))';
    Y = squeeze(PSTH_replay(cells_inc_TF,ii_dir,:))';
    ccc = corr(X,Y,'rows','pairwise','type',corr_type);
    mask_data = eye(size(ccc),'logical');
    mask_shuffle = tril(true(size(ccc)),-1);
    same_bat_TF = cells_bats_num(cells_inc_TF) == cells_bats_num(cells_inc_TF)';
    mask_shuffle = mask_shuffle & same_bat_TF;
    ccc_data{ii_dir} = ccc(mask_data);
    ccc_shuffle{ii_dir} = ccc(mask_shuffle);
end

%% plot
axes(panels{1}(1));
cla reset
hold on
histogram(cat(1,ccc_data{:}),'FaceColor',0.5*[1 1 1],'normalization','pdf','NumBins',21,'EdgeColor','None');
histogram(cat(1,ccc_shuffle{:}),'EdgeColor','k','DisplayStyle','Stairs','normalization','pdf','LineWidth',1.5);
title(corr_type)

%% add panel letters
font_size = 11;
axes(panels{1}(1))
text(-0.08,1.15, 'a', 'Units','normalized','FontWeight','bold','FontSize',font_size);

%%
fig_name = sprintf('%s',fig_name_str);
file_out = fullfile(res_dir, fig_name);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
disp('figure saved!')


%%
features = {
    'nSpikes',
    'nReplays',
    'replayTotalDuration',
    'valid_env_coverage_m',
    'valid_env_coverage_prc',
    'valid_fields_coverage_m',
    'valid_fields_coverage_prc',
    'maps_corr_rho_pearson',
    'maps_corr_rho_spearman',};
figure
M = [];
for ii_fn = 1:length(features)
    fn = features{ii_fn};
    M(:,ii_fn) = [replay_FR_maps(3,:,:).(fn)];
end
h=gplotmatrix(M,[],[],[],[],[],[],[],features)


%%
