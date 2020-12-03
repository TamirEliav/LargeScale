%% Large Scale - Fig. 4 - comparing datasets (different track sizes)

%%
clear 
clc

%% params
nboot = 10000;
alpha = 0.05; % TODO: decide if one or two sided test
% CV_err_bar = 'CI';
% CV_err_bar = 'SEM';
CV_err_bar = 'SD';
use_large_bins_results = 0;
norm_speed_ratio = 0;
norm_dist_from_center = 0;
use_only_center = 0;
center_part_len = [nan 6 300];

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'Fig_4';
fig_caption_str = 'Comparing datasets';
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

%% report figure params
fprintf('nboot = %d\n',nboot);
fprintf('CV_err_bar = %s\n',CV_err_bar);
fprintf('normalize speed ratio = %d\n',norm_speed_ratio);
fprintf('Use only the center part of the small environment = %d\n',use_only_center);
if use_only_center
    fprintf('\tBat small no change\n');
    fprintf('\tBat small center part length = %dm\n',center_part_len(2));
    fprintf('\tRat small center part length = %dm\n',center_part_len(3)/100);
end


%% create figure
% =========================================================================
% figure_size_cm = [21.0 29.7]; % ~A4
figure_size_cm = [21.6 27.9]; % ~US letter
figure ;
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
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');
pause(0.2); % workaround to solve matlab automatically changing the axes positions...

% create panels
panel_AB_size_raster = [3 1];
panel_AB_size_FR_map = [3 1];
panel_AB_pos = [2 16.6];
panel_AB = [];
for ii_dataset =1:2
    for ii_example=1:5
        offset_x = (ii_example-1)*3.6;
        offset_y = (ii_dataset-1)*5.1;
        offset = panel_AB_pos + [offset_x offset_y];
        panel_AB(ii_dataset,ii_example,1) = axes('position', [offset+[0 2.25] panel_AB_size_FR_map]);
        panel_AB(ii_dataset,ii_example,2) = axes('position', [offset+[0 1   ] panel_AB_size_raster]);
        panel_AB(ii_dataset,ii_example,3) = axes('position', [offset+[0 0   ] panel_AB_size_raster]);
    end
end
panel_AB = flipdim(panel_AB,1);
% now panel_AB is 2X5X3 = dataset X example X raster/FR_map

panel_CDE_size = [2.3 2.3];
panel_CDE_pos = [4 6.7];
panel_CDE = [];
for ii_dataset=1:3
    for ii_feature=1:3
        offset_x = (ii_feature-1)*4;
        offset_y = (ii_dataset-1)*3;
        offset = panel_CDE_pos + [offset_x offset_y];
        panel_CDE(ii_dataset,ii_feature) = axes('position', [offset panel_CDE_size]);
    end
end
panel_CDE = flipdim(panel_CDE,1);
% panel_B = permute(panel_B,[2 1]);
% now panel_B is 3X3 = dataset X feature

panel_FGHI_size = [2.5 2.8];
panel_FGHI_pos = [17 2.3];
panel_FGHI = [];
for r=1:4
    for c=1:1
        offset_x = (c-1)*3.3;
        offset_y = (r-1)*3.3;
        offset = panel_FGHI_pos + [offset_x offset_y];
        panel_FGHI(r,c) = axes('position', [offset panel_FGHI_size]);
    end
end
panel_FGHI = flipdim(panel_FGHI,1);

if norm_dist_from_center
%     norm_dist_from_center_panels(1) = axes('position', [3  2.5 3 3]);
    norm_dist_from_center_panels(2) = axes('position', [3 2 3 3]);
    norm_dist_from_center_panels(3) = axes('position', [9 2 3 3]);
end

%% load data ==============================================================
clear datasets
data_dir = 'L:\processed_data_structs';
if use_large_bins_results
    % supp fig
    datasets(1)=load(fullfile(data_dir,'cells_bat_200m.mat'));
    datasets(2)=load(fullfile(data_dir,'cells_bat_6m_binning_option3.mat'));
    datasets(3)=load(fullfile(data_dir,'cells_rat_3m_signif+FE_bin_5cm_ker_10cm.mat'));
else
    % main fig 
    datasets(1)=load(fullfile(data_dir,'cells_bat_200m.mat'));
    datasets(2)=load(fullfile(data_dir,'cells_bat_6m.mat'));
    datasets(3)=load(fullfile(data_dir,'cells_rat_3m_signif+FE.mat'));
end
datasets(1).name = {'Bat large env.'};
datasets(2).name = {'Bat small env.'};
datasets(3).name = {'Rat small env.'};

% fix struct names for 6m dataset
[datasets(2).cells.stats]  = disperse([datasets(2).cells.stats_6m]);
% [datasets(2).cells.signif] = disperse([datasets(2).cells.signif_6m]);
% [datasets(2).cells.FR_map] = disperse([datasets(2).cells.FR_map_6m]);
% [datasets(2).cells.stats] = disperse(arrayfun(@(cell)(cell.stats_6m([1 2])), datasets(2).cells, 'UniformOutput', false));
[datasets(2).cells.signif] = disperse(arrayfun(@(cell)(cell.signif_6m([1 2])), datasets(2).cells, 'UniformOutput', false));
[datasets(2).cells.FR_map] = disperse(arrayfun(@(cell)(cell.FR_map_6m([1 2])), datasets(2).cells, 'UniformOutput', false));
[datasets(2).cells.fields] = disperse(arrayfun(@(cell)(cell.fields_6m([1 2])), datasets(2).cells, 'UniformOutput', false));
[datasets(2).cells.FE] = disperse(arrayfun(@(cell)(cell.FE_6m([1 2])), datasets(2).cells, 'UniformOutput', false));
datasets(2).cells = rmfield(datasets(2).cells, 'stats_6m');
datasets(2).cells = rmfield(datasets(2).cells, 'signif_6m');
datasets(2).cells = rmfield(datasets(2).cells, 'FR_map_6m');
datasets(2).cells = rmfield(datasets(2).cells, 'fields_6m');
datasets(2).cells = rmfield(datasets(2).cells, 'FE_6m');
datasets(2).cells = rmfield(datasets(2).cells, 'cluster_quality');

%%
prm = PARAMS_GetAll();

%% cell examples options
% panel_A_opt = [3 1 8 11 2];
% panel_A_opt = [3 1 8 11 13];
panel_A_opt = [1 11 8 2 3];

% panel_B_opt = [3 4 8 1 2];
% panel_B_opt = [17:21];
% panel_B_opt = [22:26];
% panel_B_opt = [27:31];
panel_B_opt = [8 3 4 30 27];

bat_examples_options = {
    'b2311_d191218_TT3_SS08'   % 1
    'b2329_d190918_TT11_SS01'
    'b2329_d190919_TT11_SS04'
    'b2382_d190620_TT9_SS03'
    'b2382_d190623_TT13_SS05'  % 5
    'b2311_d191209_TT3_SS04'
    'b2311_d191218_TT3_SS02'
    'b2311_d191218_TT3_SS05'
    'b2329_d190923_TT2_SS02'
    'b2329_d190916_TT7_SS01'   % 10
    'b2311_d191218_TT3_SS07'
    'b2311_d191218_TT3_SS04'
    'b2329_d190918_TT7_SS01'   % 13
};
rat_examples_options = [
    386     % 1
    725
    738
    742
    1118    % 5
    1167
    1196
    1346
    1555
    1689    % 10
    1846
    2070
    2233
    2336
    2378    % 15
    2569
    2391
    1345
    1456
    1651    % 20
    1929
    1451
    1366
    913
    2657    % 25
    915
    724
    2628
    1872
    1285    % 30
    842
    ];

% report cell examples ID
fprintf('========= Panels A+B examples =========\n')
fprintf('Bat examples:\n')
fprintf('\t%s\n',bat_examples_options{panel_A_opt})
fprintf('Rat examples:\n')
fprintf('\t%d\n',rat_examples_options(panel_B_opt))

%% FR map + rasters - examples ============================================
cells2 = datasets(2).cells;
cells3 = datasets(3).cells;
details2=[cells2.details];
details3=[cells3.details];
IX2 = cellfun(@(cell_ID)find(contains({details2.cell_ID},cell_ID)), bat_examples_options(panel_A_opt))';
IX3 = arrayfun(@(id)find(ismember([details3.id],id)), rat_examples_options(panel_B_opt))';
cell_examples = [cells2(IX2);cells3(IX3)];
if 1
for ii_dataset = 1:size(cell_examples,1)
    for ii_cell = 1:size(cell_examples,2)
        %%
        cell = cell_examples(ii_dataset,ii_cell);
        c = prm.graphics.colors.flight_directions;
        
        switch ii_dataset
            case 1
                xlimits = [0 6];
%                 area_limits = [6 12];
                area_size = 0.30;
%                 rescale_gain_offset = [1 6];
                rescale_gain_offset = [1 mean(cell.stats.all.valid_speed_pos)-3];
                cell_str = cell.details.cell_ID;
                m_or_cm = 'm';
            case 2
                xlimits = [0 3];
%                 area_limits = [-150 150];
                area_size = 30;
%                 rescale_gain_offset = [1/100 -150];
                rescale_gain_offset = [1/100 mean(cell.stats.all.valid_speed_pos)-150];
                cell_str = "cell "+cell.details.id;
                m_or_cm = 'cm';
        end

        % map+fields
        axes(panel_AB(ii_dataset, ii_cell, 1));
        cla
        hold on
        maps=[cell.FR_map.all];
        x = maps(1).bin_centers;
        y = cat(1,maps.PSTH);
        m = round(max(y(:)));
        ylimits = [0 m+1];
        % low-speed area
        area_lowerval = ylimits(1) - 0.23*range(ylimits);
        area_upperval = ylimits(1) - 1e-2*range(ylimits);
%         area([area_limits(1)  cell.stats.all.valid_speed_pos(1)               ], repelem(area_upperval,2), area_lowerval, 'FaceColor',0.7*[1 1 1],'EdgeColor','none','ShowBaseLine','off','Clipping','off');
%         area([                cell.stats.all.valid_speed_pos(2) area_limits(2)], repelem(area_upperval,2), area_lowerval, 'FaceColor',0.7*[1 1 1],'EdgeColor','none','ShowBaseLine','off','Clipping','off');
        area_X = cell.stats.all.valid_speed_pos(1) - area_size;
        area_Y = repelem(area_upperval,2);
        area(cell.stats.all.valid_speed_pos(1) - [area_size 0], repelem(area_upperval,2), area_lowerval, 'FaceColor',0.7*[1 1 1],'EdgeColor','none','ShowBaseLine','off','Clipping','off');
        area(cell.stats.all.valid_speed_pos(2) + [0 area_size], repelem(area_upperval,2), area_lowerval, 'FaceColor',0.7*[1 1 1],'EdgeColor','none','ShowBaseLine','off','Clipping','off');
        h=plot(x,y,'LineWidth',1.2);
        [h.Color] = disperse(c);
        % fields
        dir_offsets = [-0.1 -0.17]+0.015;
        for ii_dir=1:2
            fields = cell.fields{ii_dir};
            if isfield(fields,'in_low_speed_area')
                fields([fields.in_low_speed_area])=[];
            end
            for ii_field = 1:length(fields)
                plot(fields(ii_field).edges_prc, ...
                     repelem(dir_offsets(ii_dir)*range(ylimits),2),...
                    'Linewidth', 2, 'Color', c{ii_dir},'Clipping','off');
                fprintf('dataset %d cell %d direction %d field %d size = %.2f%s\n',...
                    ii_dataset, ii_cell, ii_dir, ii_field, fields(ii_field).width_prc, m_or_cm);
            end
        end
        rescale_plot_data('x',rescale_gain_offset);
        box off
        h=gca;
        h.TickDir = 'out';
        h.XTick = [];
        h.YTick = [0 m];
        h.YLim = ylimits;
        h.XLim = xlimits;
        h.YRuler.TickLabelGapOffset = -1;
%         title(cell_str,'FontSize',7,'FontWeight','normal','Interpreter','none');

% % %         % cell details
% % %         cell_num_str_pos_x   = [0.50 0.50 0.50 0.50 0.50 0.45 0.50 0.50 0.50];
% % %         cell_num_str_pos_y   = [1.05 1.05 1.05 0.85 0.90 0.90 0.90 0.90 0.90];
% % %         cell_stats_str_pos_x = [0.80 0.95 0.80 0.80 0.80 0.80 0.12 0.50 0.85];
% % %         cell_stats_str_pos_y = [1.20 1.05 1.15 0.90 1.10 1.10 1.10 0.90 0.90]+0.05;
% % %         text(cell_num_str_pos_x(ii_cell), cell_num_str_pos_y(ii_cell), "cell "+ ii_cell,...
% % %             'Units','normalized','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
% % %         switch ii_cell
% % %             case {1,2,3,4,5,6,7,8}
% % %                 cell_stats_str = {  sprintf('max=%.1fm', cell.stats.all.field_largest);...
% % %                                     sprintf('min=%.1fm', cell.stats.all.field_smallest);...
% % %                                     sprintf('ratio=%.1f', cell.stats.all.field_ratio_LS);...
% % %                                  };
% % %                 text(cell_stats_str_pos_x(ii_cell), cell_stats_str_pos_y(ii_cell)-0*0.23, cell_stats_str{1},...
% % %                     'Units','normalized','HorizontalAlignment','center','VerticalAlignment','Top','FontSize',6);
% % %                 text(cell_stats_str_pos_x(ii_cell), cell_stats_str_pos_y(ii_cell)-1*0.23, cell_stats_str{2},...
% % %                     'Units','normalized','HorizontalAlignment','center','VerticalAlignment','Top','FontSize',6);
% % %                 text(cell_stats_str_pos_x(ii_cell), cell_stats_str_pos_y(ii_cell)-2*0.23, cell_stats_str{3},...
% % %                     'Units','normalized','HorizontalAlignment','center','VerticalAlignment','Top','FontSize',6);
% % %             case 9
% % %                 cell_stats_str = {  sprintf('single field=%.1fm', cell.fields{1}.width_prc) };
% % %                 text(cell_stats_str_pos_x(ii_cell), cell_stats_str_pos_y(ii_cell)-0*0.23, cell_stats_str{1},...
% % %                     'Units','normalized','HorizontalAlignment','center','VerticalAlignment','Top','FontSize',6);
% % %         end

        % rasters
        FEs = [cell.FE];
        for ii_dir=1:2
            axes(panel_AB(ii_dataset, ii_cell, ii_dir+1));
            cla
            FE = FEs{ii_dir};
            x = [FE.spikes_pos];
            [FE.number2] = disperse(1:length(FE));
            y = arrayfun(@(FE)(FE.number2*ones(1,FE.num_spikes)),FE,'UniformOutput',0);
            y = [y{:}];
            plot(x,y,'.','Color',c{ii_dir},'MarkerSize',0.05);
            
            rescale_plot_data('x',rescale_gain_offset);
            
            box off
            h=gca;
            m = length(FE);
            h.YTick = [1 m];
            h.XLim = xlimits;
            h.YLim = [0 m+1];
            h.YRuler.TickLabelGapOffset = -1;
            switch ii_dir
                case 1
                    h.XTick = [];
                    h.YTickLabel = {'',num2str(m)};
                case 2
                    switch ii_dataset
                        case 1
                            h.XTick = 0:1:6;
                        case 2
                            h.XTick = 0:1:3;
                    end
                    h.XRuler.TickLabelGapOffset = -2;
                    h.YTickLabel = {'1',num2str(m)};
                    h.TickDir = 'out';
            end
        end
        
        % fields num (here to be above the - TODO: Move up where I plot fields
        for ii_dir=1:2
            fields = cell.fields{ii_dir};
            if isfield(fields,'in_low_speed_area')
                fields([fields.in_low_speed_area])=[];
            end
            dir_offsets =2.05+[0.18 0];
            text(1.05, dir_offsets(ii_dir), num2str(length(fields)), 'Units','normalized', 'Color', c{ii_dir},...
                'HorizontalAlignment','right', 'VerticalAlignment','middle','FontSize',6);
        end
    end
end

%% add direction arrows
arrow_x = 0.85 +[0 0.05];
arrow_y = repelem(0.95,2);
clear h
h(1)=annotation('arrow',arrow_x,      arrow_y+0.008,  'Color', prm.graphics.colors.flight_directions{1});
h(2)=annotation('arrow',flip(arrow_x),arrow_y      ,  'Color', prm.graphics.colors.flight_directions{2});
[h.HeadWidth] = disperse([5 5]);
[h.HeadLength] = disperse([5 5]);

%% add x/y labels for specific panels
for ii_dataset = 1:2
    for ii_cell = 1:5
        axes(panel_AB(ii_dataset,ii_cell, 3));
        xlabel('Position (m)', 'Units','normalized','Position',[0.5 -0.35]);
    end
end
for ii_dataset = 1:2
    ii_cell = 1;
    axes(panel_AB(ii_dataset,ii_cell, 3));
    ylabel('Lap no.',   'Units','normalized','Position',[-0.2 1]);
    axes(panel_AB(ii_dataset,ii_cell, 1));
    ylabel({'Firing rate';'(Hz)'},   'Units','normalized','Position',[-0.15 0.42]);
end
axes(panel_AB(1,1,1));
text(-0.4,1.6, 'A', 'Units','normalized','FontWeight','bold');
text(0.1,1.6, 'Bats - Small environment', 'Units','normalized','FontWeight','bold');
axes(panel_AB(2,1,1));
text(-0.4,1.6, 'B', 'Units','normalized','FontWeight','bold');
text(0.1,1.6, 'Rats - Small environment', 'Units','normalized','FontWeight','bold');

end % if plot examples ====================================================



%% ================ plot datasets histograms ========================
fprintf('======= panels CDE: plot datasets distributions ======= \n');
field_num_val = [];
field_num_grp = [];
field_size_val = [];
field_size_grp = [];
ratio_LS_val = [];
ratio_LS_grp = [];
ratio_LS_with_1s_val = [];
ratio_LS_with_1s_grp = [];
n_str_pos_x = [1 1 1.1; 1 1 1.1; 1 1 1.1];
n_str_pos_y = [1 1 1; 0.9 1 0.9; 0.9 1 0.9];
for ii_dataset = 1:length(datasets)
    %% arrange dataset struct
    name = datasets(ii_dataset).name;
    cells = datasets(ii_dataset).cells;
    signif = cat(1,cells.signif);
    signif = arrayfun(@(x)(x.TF),signif);
    signif = any(signif,2);
    cells(~signif)=[];
    signif = cat(1,cells.signif);
    signif = arrayfun(@(x)(x.TF),signif);
    stats = [cells.stats];
    stats_all = [stats.all];
    stats_dir = cat(1,stats.dir);
    fields = cat(1,cells.fields);
    fields = fields(signif);
    fields = [fields{:}];
    fields([fields.in_low_speed_area])=[];
% % %     if ~isfield(fields,'in_low_speed_area')
% % %         [fields.in_low_speed_area] = disperse(false(size(fields)));
% % %     end
% % %     if is_remove_low_speed 
% % %         fields([fields.in_low_speed_area])=[];
% % %     end
    
    %% arrange relevant data to plot
    field_num = [stats_dir(signif).field_num];
    field_size = [fields.width_prc];
    ratio_LS = [stats_all.field_ratio_LS];
    ratio_LS_with_1s = [stats_all.field_ratio_LS];
    ratio_LS_with_1s(isnan(ratio_LS_with_1s))=1;
    ratio_LS_vel = abs([stats_all.field_ratio_LS_vel]);
    ratio_LS_vel(isnan(ratio_LS_vel)) = 1;
        
    % Here we try two options to overcome the parabolic speed trajectories:
    % test 1: take only the very central part of the track
    if use_only_center & ~isnan(center_part_len(ii_dataset)) % only for small env.
        % remove fields far from the center        
        % and recalculate field num/size/ratio
        field_num = [];
        field_size = [];
        ratio_LS = [];
        ratio_LS_with_1s = [];
        ratio_LS_vel = [];
        dist2center_LS = [];
        for ii_cell = 1:length(cells)
            cell_fields = [];
            for ii_dir=1:2
                cell = cells(ii_cell);
                if ~cell.signif(ii_dir).TF
                    continue;
                end
                cell_fields_dir = cell.fields{ii_dir};
                cell_fields_dir( [cell_fields_dir.in_low_speed_area] ) = [];
                center_loc = mean(cell.stats.all.valid_speed_pos);
                cell_fields_dir( abs([cell_fields_dir.loc]-center_loc) > center_part_len(ii_dataset)/2 ) = [];
                cell_fields = [cell_fields cell_fields_dir];
                field_num = [field_num length(cell_fields_dir)];
            end
            cell_fields_sizes = [cell_fields.width_prc];
            cell_fields_vel = [cell_fields.vel];
            cell_fields_loc = [cell_fields.loc];
            field_size = [field_size cell_fields_sizes];
            switch length(cell_fields_sizes)
                case 0 % do nothing (disregard this cell)
                case 1 % add 1 only to 'with_1s' array and nan to the regular array (for code comatibility)
                    ratio_LS         = [ratio_LS            nan];
                    ratio_LS_with_1s = [ratio_LS_with_1s      1];
                    ratio_LS_vel     = [ratio_LS_vel          1];
                    dist2center_LS   = [dist2center_LS; [1 1].*abs(cell_fields_loc-center_loc)];
                otherwise
                    [field_L,IX_L] = max(cell_fields_sizes);
                    [field_S,IX_S] = min(cell_fields_sizes);
                    ratio_LS         = [ratio_LS            field_L/field_S ];
                    ratio_LS_with_1s = [ratio_LS_with_1s    field_L/field_S ];
                    ratio_LS_vel     = [ratio_LS_vel        cell_fields_vel(IX_L) / cell_fields_vel(IX_L) ];
                    dist2center_LS   = [dist2center_LS;     abs(cell_fields_loc(IX_L)-center_loc) ...
                                                            abs(cell_fields_loc(IX_S)-center_loc)    ];
                    
            end
        end
        field_num(field_num==0)=[];
    end
    % test 2: normalize by the speed ratio
    if norm_speed_ratio
        ratio_LS = ratio_LS ./ ratio_LS_vel;
        ratio_LS_with_1s = ratio_LS_with_1s ./ ratio_LS_vel;
    end
    % test 3: normalize by the distance from center ratio
    if norm_dist_from_center & ~isnan(center_part_len(ii_dataset)) % only for small env.
        ratio_LS_dist2center = dist2center_LS(:,1)' ./ dist2center_LS(:,2)';
        % here we multiply by ratio_dist2center, in contrast to normalizing
        % by the speed ratio, because the further the field from the center
        % we expect it to be smaller, so we expect ratio_LS_dist2center to
        % be values < 1
        ratio_LS = ratio_LS .* ratio_LS_dist2center;
        ratio_LS_with_1s = ratio_LS_with_1s .* ratio_LS_dist2center;
    end
    % convert cm to m
    if ii_dataset==3
        field_size = field_size / 100; % convert cm to m
    end
    
    %% store data from all datasets together
    field_num_val  = [field_num_val   field_num];
    field_size_val = [field_size_val  field_size];
    ratio_LS_val   = [ratio_LS_val    ratio_LS];
    ratio_LS_with_1s_val   = [ratio_LS_with_1s_val   ratio_LS_with_1s];
    field_num_grp  = [field_num_grp   ii_dataset.*ones(size(field_num))];
    field_size_grp = [field_size_grp  ii_dataset.*ones(size(field_size))];
    ratio_LS_grp   = [ratio_LS_grp    ii_dataset.*ones(size(ratio_LS))];
    ratio_LS_with_1s_grp   = [ratio_LS_with_1s_grp   ii_dataset.*ones(size(ratio_LS_with_1s))];
    
    %%
    datasets(ii_dataset).field_size = field_size;
    datasets(ii_dataset).fields = fields;
    datasets(ii_dataset).signif = signif;
    datasets(ii_dataset).signif = signif;
    datasets(ii_dataset).stats_all = stats_all;
    datasets(ii_dataset).stats_dir = stats_dir;
    
    %% report dataset
    fprintf('Dataset %d: %s\n', ii_dataset, datasets(ii_dataset).name{1})
    
    %% no. of fields
    x = field_num;
    axes(panel_CDE(ii_dataset,1));
    cla('reset');
    hold on
    hh=histogram(x);
    hh.FaceColor = 0.5*[1 1 1];
    hh.BinEdges = 0.5+[0:20];
    hh.Data(hh.Data > hh.BinLimits(2)) = hh.BinLimits(2);
    if hh.Values(end)
        fprintf('\tLast bin includes all cells with > %d fields\n', floor(hh.BinEdges(end)));
    end
    
    hax=gca;
    hax.YScale = 'log';
    hax.YLim(1) = 0.7;
    hax.YLim(2) = 1.15 * hax.YLim(2);
%     ha.YLim = [0.7e0 240];
    hax.XLim = [0 hh.BinLimits(2)];
    hax.XTick = [0:5:30];
    hax.YTick = [1 10 100];
    hax.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
    hax.TickDir='out';
    hax.TickLength = [0.03 0.03];
    hax.XRuler.TickLabelGapMultiplier = -0.3;
    hax.YRuler.TickLabelGapMultiplier = 0.001;
    if ii_dataset == 3
        xlabel({'No. of fields per direction'},'Units','normalized','Position',[0.5 -0.18]);
    end
    ylabel('No. of cells','Units','normalized','Position',[-0.28 0.5])
    text(n_str_pos_x(ii_dataset,1),n_str_pos_y(ii_dataset,1), "n = "+sum(~isnan(x)),...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
    
    hl=xline(nanmean(x)); hl.Color='r';
    m = hax.YLim(2) + 0.15*range(hax.YLim);
    plot(prctile(x,[25 75]), [m m], 'r-','LineWidth',1   ,'Clipping','off');
    plot(prctile(x,[50]),    m    , 'r.','MarkerSize',10 ,'Clipping','off');
    fprintf('\tNo. of fields: mean=%.1f median=%.1f IQR=%.1f-%.1f\n', nanmean(x), prctile(x,[50]), prctile(x,[25 75]) );
    
    %% Field Size
    x = field_size;
    axes(panel_CDE(ii_dataset,2));
    cla('reset');
    hold on
    h=histogram(x);
    h.FaceColor = 0.5*[1 1 1];
    h.BinEdges = 0:33;

    ylabel('counts')
    hax=gca;
    hax.YScale = 'log';
    hax.YLim(1) = 0.8;
    hax.YLim(2) = 1.15 * hax.YLim(2);
    if ii_dataset == 3
        hax.YLim(2) = 1.15 * hax.YLim(2);
    end
    hax.YTick = [1 10 100];
    hax.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
    hax.TickDir='out';
    hax.TickLength = [0.03 0.03];
    hax.XRuler.TickLabelGapMultiplier = -0.3;
    hax.YRuler.TickLabelGapMultiplier = 0.001;
    if ii_dataset == 3
        xlabel('Field size (m)')
    end
    ylabel('No. of fields','Units','normalized','Position',[-0.28 0.5])
    text(n_str_pos_x(ii_dataset,2),n_str_pos_y(ii_dataset,2), "n = "+sum(~isnan(x)),...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
    
    hl=xline(nanmean(x)); hl.Color='r';
    m = hax.YLim(2) + 0.15*range(hax.YLim);
    plot(prctile(x,[25 75]), [m m], 'r-','LineWidth',1   ,'Clipping','off');
    plot(prctile(x,[50]),    m    , 'r.','MarkerSize',10 ,'Clipping','off');
    fprintf('\t Fields sizes: mean=%.1fm median=%.1fm IQR=%.1f-%.1fm\n', nanmean(x), prctile(x,[50]), prctile(x,[25 75]) );
    [fields_size_gamma_fit] = gamfit(field_size);
    fprintf('\t smallest fields gamma distribution fit:\n')
    fprintf('\t\t Shape parameter kappa = %.2f\n', fields_size_gamma_fit(1));
    fprintf('\t\t Scale parameter theta = %.2fm\n', fields_size_gamma_fit(2));
            
    % inset - zoom in on x-axis
    if ismember(ii_dataset,[2 3])
        hax = copyobj(gca,gcf);
        axes(hax);
        hax.Position([1 2]) = hax.Position([1 2]) + [1 0.5];
        hax.Position([3 4]) = [1.2 1.3];
        delete(findobj(gca,'Type','Text'));
        hh = findobj(gca,'Type','Histogram');
        hh.BinEdges = linspace(0,max(hh.Data)*1.1,6);
        hax.FontSize=6;
        hax.XLabel.String = '';
        hax.YLabel.String = '';
    end
    
    %% ratio L/S
    x = ratio_LS;
    x2 = ratio_LS_with_1s;
    axes(panel_CDE(ii_dataset,3));
    cla('reset');
    hold on
    nBinEdges = 9;
    edges = logspace(0,log10(25),nBinEdges);
    h=histogram(x2);
    h.BinEdges = edges;
    h.FaceColor = 'k';
    h.FaceAlpha = 1;
    h=histogram(x);
    h.BinEdges = edges;
    h.FaceColor = 0.5*[1 1 1];
    h.FaceAlpha = 1;
    
    hax=gca;
    hax.XScale = 'log';
    hax.YScale = 'log';
    hax.YLim(1) = 0.7;
    hax.YLim(2) = 1.15 * hax.YLim(2);
    hax.XLim = [0 27];
    hax.YTick = [1 10 100];
    hax.XTick = [1 2 5 10 20];
    hax.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
    hax.TickDir='out';
    hax.TickLength = [0.03 0.03];
    hax.XRuler.TickLabelGapMultiplier = -0.35;
    hax.YRuler.TickLabelGapMultiplier = 0.001;
    if ii_dataset == 3
        xlabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[0.5 -0.17]);
    end
    ylabel('No. of cells','Units','normalized','Position',[-0.28 0.5])
    text(n_str_pos_x(ii_dataset,3),n_str_pos_y(ii_dataset,3), "n = "+length(x),...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
    fprintf('\tnumber of cells with single field overall (both directions): %d/%d (%.1f%%)\n',sum(isnan(x)), length(x), 100*sum(isnan(x))/length(x));
    
    hl=xline(nanmean(x2)); hl.Color='r';
    m = hax.YLim(2) + 0.15*range(hax.YLim);
    plot(prctile(x2,[25 75]), [m m], 'r-','LineWidth',1   ,'Clipping','off');
    plot(prctile(x2,[50]),    m    , 'r.','MarkerSize',10 ,'Clipping','off');
    fprintf('\tFields size ratio (L/S): mean=%.1f median=%.1f IQR=%.1f-%.1f\n', nanmean(x2), prctile(x2,[50]), prctile(x2,[25 75]) );
    
    %% ratio L/S vs ratio dist2center
    if norm_dist_from_center & ~isnan(center_part_len(ii_dataset)) % only for small env.
        x = ratio_LS_dist2center;
        y = ratio_LS ./ ratio_LS_dist2center;
        axes(norm_dist_from_center_panels(ii_dataset));
        cla('reset');
        hold on
        plot(x,y,'k.','MarkerSize',7);
        xlabel({'Distance to center ratio';'largest/smallest'});
        ylabel({'Field size ratio';'largest/smallest'});
        title(datasets(ii_dataset).name);
        hax=gca;
        hax.XScale = 'log';
        hax.YScale = 'log';
        [rho,rho_pval] = corr(x',y','rows','pairwise','type','Spearman')
        text(1,.95,sprintf('rho = %.2f',rho),'Units','normalized','FontSize',8);
        text(1,.85,sprintf('P = %.2f',rho_pval),'Units','normalized','FontSize',8);
    end
end

% add panel letter
axes(panel_CDE(1,1));
text(-0.45,1.1, 'C', 'Units','normalized','FontWeight','bold');
axes(panel_CDE(1,2));
text(-0.45,1.1, 'D', 'Units','normalized','FontWeight','bold');
axes(panel_CDE(1,3));
text(-0.45,1.1, 'E', 'Units','normalized','FontWeight','bold');

%% add dataset title
axes(panel_CDE(1,1));
text(-0.9,0.5, {'Bat';'large';'environment'}, 'Units','normalized','FontWeight','bold','FontSize',9,'HorizontalAlignment','center');
axes(panel_CDE(2,1));
text(-0.9,0.5, {'Bat';'small';'environment'}, 'Units','normalized','FontWeight','bold','FontSize',9,'HorizontalAlignment','center');
axes(panel_CDE(3,1));
text(-0.9,0.5, {'Rat';'small';'environment'}, 'Units','normalized','FontWeight','bold','FontSize',9,'HorizontalAlignment','center');

%%  ============ compare datasets distributions ============ 
FGHI_xlimits = [0.25 3.75];
asterisk_font_size = 11;
fprintf('======= panels FGHI: compare datasets distributions ======= \n');

%%  No. of fields
axes(panel_FGHI(1)); cla; hold on
text(-0.55,1., 'F', 'Units','normalized','FontWeight','bold');
n=accumarray(field_num_grp',field_num_val',[],@(x)(sum(~isnan(x))));
x=unique(field_num_grp);
y=accumarray(field_num_grp',field_num_val',[],@mean);
err=accumarray(field_num_grp',field_num_val',[],@nansem);
bar(x,y,'FaceColor',0.5*[1 1 1]);
h = errorbar(x,y,err);
h.Color = [0 0 0];                            
h.LineStyle = 'none';  
hax=gca;
hax.XLim = FGHI_xlimits;
hax.XTick=1:length(datasets);
hax.XTickLabel=[];
hax.YRuler.TickLabelGapOffset = 1;
% hax.XTickLabel=[1:length(datasets)]+":"+string({datasets.name});
% hax.XTickLabelRotation = 30;
ylabel('No. of fields per direction');
% text(x, repelem(hax.YLim(2),length(x)), "n = "+n,...
%     'HorizontalAlignment','center','VerticalAlignment','bottom');
% signif test
X1 = field_num_val(field_num_grp==1);
X2 = field_num_val(field_num_grp==2);
X3 = field_num_val(field_num_grp==3);
[~,pval12,~,tstat12] = ttest2(X1,X2,'tail','right');
[~,pval13,~,tstat13] = ttest2(X1,X3,'tail','right');
pval12_rank = ranksum(X1,X2,'tail','right');
pval13_rank = ranksum(X1,X3,'tail','right');
fprintf('Field num 1 vs. 2: ttest p=%.4g df=%d, ranksum p=%.4g\n',pval12, tstat12.df, pval12_rank);
fprintf('Field num 1 vs. 3: ttest p=%.4g df=%d, ranksum p=%.4g\n',pval13, tstat13.df, pval13_rank);
str12 = signif_astricks(pval12);
str13 = signif_astricks(pval13);
plot([1 2],5.7*[1 1],'k-')
plot([1 3],6.5*[1 1],'k-')
text(1.5,  5.9,str12,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',asterisk_font_size);
text(2,    6.7,str13,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',asterisk_font_size);

%%  fields size
axes(panel_FGHI(2)); cla; hold on
text(-0.55,1., 'G', 'Units','normalized','FontWeight','bold');
n=accumarray(field_size_grp',field_size_val',[],@(x)(sum(~isnan(x))));
x=unique(field_size_grp);
y=accumarray(field_size_grp',field_size_val',[],@mean);
err=accumarray(field_size_grp',field_size_val',[],@nansem);
bar(x,y,'FaceColor',0.5*[1 1 1]);
h = errorbar(x,y,err);    
h.Color = [0 0 0];                            
h.LineStyle = 'none';  
hax=gca;
hax.XLim = FGHI_xlimits;
hax.XTick=1:length(datasets);
hax.XTickLabel=[];
hax.YRuler.TickLabelGapOffset = 1;
% hax.XTickLabel=[1:length(datasets)]+":"+string({datasets.name});
% hax.XTickLabelRotation = 30;
ylabel('Field size (m)');
% text(x, repelem(hax.YLim(2),length(x)), "n = "+n,...
%     'HorizontalAlignment','center','VerticalAlignment','bottom');
% signif test
X1 = field_size_val(field_size_grp==1);
X2 = field_size_val(field_size_grp==2);
X3 = field_size_val(field_size_grp==3);
[~,pval12,~,tstat12] = ttest2(X1,X2,'tail','right');
[~,pval13,~,tstat13] = ttest2(X1,X3,'tail','right');
pval12_rank = ranksum(X1,X2,'tail','right');
pval13_rank = ranksum(X1,X3,'tail','right');
fprintf('Field size 1 vs. 2: ttest p=%.4g df=%d, ranksum p=%.4g\n',pval12, tstat12.df, pval12_rank);
fprintf('Field size 1 vs. 3: ttest p=%.4g df=%d, ranksum p=%.4g\n',pval13, tstat13.df, pval13_rank);
str12 = signif_astricks(pval12);
str13 = signif_astricks(pval13);
plot([1 2],6.3*[1 1],'k-')
plot([1 3],7.1*[1 1],'k-')
text(1.5,  6.5,str12,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',asterisk_font_size);
text(2,    7.3,str13,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',asterisk_font_size);

%%  field size CV
axes(panel_FGHI(3)); cla; hold on
text(-0.55,1., 'H', 'Units','normalized','FontWeight','bold');
[CV,CI,SD,SEM] = splitapply(@(X)CV_BS(X,nboot,alpha), field_size_val', field_size_grp');
x = unique(field_size_grp);
y = CV;
switch CV_err_bar
    case 'CI'
        err1 = CV-CI(:,1);
        err2 = CI(:,2)-CV;
    case 'SEM'
        err1 = SEM;
        err2 = SEM;
    case 'SD'
        err1 = SD;
        err2 = SD;
end
bar(x,y,'FaceColor',0.5*[1 1 1]);
h = errorbar(x,y,err1,err2);
h.Color = [0 0 0];                            
h.LineStyle = 'none';  
hax=gca;
hax.XLim = FGHI_xlimits;
hax.XTick=1:length(datasets);
hax.XTickLabel=[];
hax.YRuler.TickLabelGapOffset = 1;
ylabel({'Field size';'coefficient of variation'});
% signif test
z12 = (y(1)-y(2)) / sqrt(SD(1)^2+SD(2)^2);
z13 = (y(1)-y(3)) / sqrt(SD(1)^2+SD(3)^2);
pval12 = 1-normcdf(z12,0,1);
pval13 = 1-normcdf(z13,0,1);
str12 = signif_astricks(pval12);
str13 = signif_astricks(pval13);
fprintf('Field size CV 1 vs. 2 (bootstrap): z-test p=%.4g z=%.4g\n',pval12,z12);
fprintf('Field size CV 1 vs. 3 (bootstrap): z-test p=%.4g z=%.4g\n',pval13,z13);
plot([1 2],0.75*[1 1],'k-')
plot([1 3],0.85*[1 1],'k-')
text(1.5,0.77,str12,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',asterisk_font_size);
text(2,0.87,str13,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',asterisk_font_size);

%%  L/S ratio
axes(panel_FGHI(4)); cla; hold on
text(-0.55,1., 'I', 'Units','normalized','FontWeight','bold');
n=accumarray(ratio_LS_with_1s_grp',ratio_LS_with_1s_val',[],@(x)(sum(~isnan(x))));
x=unique(ratio_LS_with_1s_grp);
y=accumarray(ratio_LS_with_1s_grp',ratio_LS_with_1s_val',[],@mean);
err=accumarray(ratio_LS_with_1s_grp',ratio_LS_with_1s_val',[],@nansem);
bar(x,y,'FaceColor',0.5*[1 1 1]);
h = errorbar(x,y,err);    
h.Color = [0 0 0];                            
h.LineStyle = 'none';  
hax=gca;
hax.XLim = FGHI_xlimits;
hax.XTick=1:length(datasets);
hax.XTickLabel=string({datasets.name});
hax.XTickLabelRotation = 45;
hax.XAxis.FontWeight='bold';
hax.XRuler.TickLabelGapOffset = -2;
hax.YRuler.TickLabelGapOffset = 1;

% hax.XTick=1:length(datasets);
% hax.XTickLabel=[];
% hax.YRuler.TickLabelGapOffset = 1;

% hax.XTickLabel=[1:length(datasets)]+":"+string({datasets.name});
% hax.XTickLabelRotation = 30;
% ,'Units','normalized','Position',[0.5 -0.17]
ylabel({'Field size ratio';'largest/smallest'});
% text(x, repelem(hax.YLim(2),length(x)), "n = "+n,...
%     'HorizontalAlignment','center','VerticalAlignment','bottom');
% signif test
X1 = ratio_LS_with_1s_val(ratio_LS_with_1s_grp==1);
X2 = ratio_LS_with_1s_val(ratio_LS_with_1s_grp==2);
X3 = ratio_LS_with_1s_val(ratio_LS_with_1s_grp==3);
[~,pval12,~,tstat12] = ttest2(X1,X2,'tail','right');
[~,pval13,~,tstat13] = ttest2(X1,X3,'tail','right');
pval12_rank = ranksum(X1,X2,'tail','right');
pval13_rank = ranksum(X1,X3,'tail','right');
fprintf('L/S ratio 1 vs. 2: ttest p=%.4g df=%d, ranksum p=%.4g\n',pval12, tstat12.df, pval12_rank);
fprintf('L/S ratio 1 vs. 3: ttest p=%.4g df=%d, ranksum p=%.4g\n',pval13, tstat13.df, pval13_rank);
str12 = signif_astricks(pval12);
str13 = signif_astricks(pval13);
plot([1 2],4.5*[1 1],'k-')
plot([1 3],5.1*[1 1],'k-')
text(1.5,4.7,str12,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',asterisk_font_size);
text(2,5.3,str13,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',asterisk_font_size);

%% print/save the figure
file_out = fig_name_str;
if use_large_bins_results
    file_out = [file_out '_larger_bins'];
end
if norm_speed_ratio
    file_out = [file_out '_norm_speed_ratio'];
end
if norm_dist_from_center
    file_out = [file_out '_norm_dist2center'];
end
if use_only_center
    file_out = [file_out '_use_only_center' sprintf('_bat_%.2f_rat_%.2f',center_part_len(2),center_part_len(3)/100)];
    file_out = strrep(file_out, '.','_');
end

% file_out = sprintf('%s_A_%d_%d_%d_%d_%d_B_%d_%d_%d_%d_%d',file_out,panel_A_opt,panel_B_opt);
file_out = fullfile(res_dir, file_out);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
fprintf('Figure was successfully saved to pdf format, file name:\n\t%s\n',file_out);

%%
function [CV,CI,SD,SEM] = CV_BS(X,nboot,alpha)
    rng(0);
    BOOTSTAT = bootstrp(nboot, @(x)(std(x)/mean(x)), X);
    CV = mean(BOOTSTAT);
    CI = quantile(BOOTSTAT,[alpha 1-alpha]);
    SD = std(BOOTSTAT);
    SEM = std(BOOTSTAT)/sqrt(length(BOOTSTAT));
end

%% signif astricks string
function str = signif_astricks(pval)
    if pval>=0.05
        str = 'n.s.';
    end
    if pval<0.05
        str = '*';
    end
    if pval<0.01
        str = '**';
    end
    if pval<0.001
        str = '***';
    end
    if pval<0.0001
        str = '****';
    end
    if pval<0.00001
        str = '*****';
    end
end







%%
