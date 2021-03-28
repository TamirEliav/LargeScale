%% Sup fig: field dynamics
%%
clear 
clc

%% add data folders to path
addpath('L:\processed_data_structs');

%% define output files
res_dir =  'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S23';%'SupFig_Dynamic';
fig_caption_str = ' ';
log_name_str = [fig_name_str '_log_file' '.txt'];
log_name_str = strrep(log_name_str , ':', '-');
log_name_str = strrep(log_name_str , ' ', '_');
log_name_out = fullfile(res_dir, log_name_str);

perturbation_prc = 4;

%% open log file
% diary off
% diary(log_name_out)
% diary on
% disp('Log file');
% disp(['created: ', datestr(clock)]);
% disp('======================================================');
% disp([fig_name_str ':' fig_caption_str]);   
% disp('======================================================');
% disp('');


%% create figure
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
panel_A(1)       = axes('position', [ 2  17  8  6]);
panel_A(2)       = axes('position', [ 11.5  17  8  6]);
panel_size = [3 3];
panel_B(1)    = axes('position', [ 2  11  panel_size]);
panel_B(2)    = axes('position', [ 6  11  panel_size]);
panel_B(3)    = axes('position', [ 10  11  panel_size]);
panel_C       = axes('position', [ 15  11  panel_size]);

panel_legend  = axes('position', [ 18.2  13  2.5  1.5]);

%% load data:
cells_130m = load('cells_bat_130m.mat');
cells_130m = cells_130m.cells;
cells_130m_IDs = arrayfun(@(x)(x.cell_ID), cat(1,cells_130m.details), 'UniformOutput', false);

%% Panel A - examples:
axes(panel_A(1)); hold on;
text(-0.1125,1.13, 'A', 'Units','normalized','FontWeight','bold');

cell_examples_IX = [1 2];
cell_examples = {
            'b2382_d190808_TT15_SS04';
            'b2382_d190709_TT11_SS01';
            };
cell_examples_dir = [1 1];   

fields_colors = [     0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    0    0    1;
    1    0    0;
    0    1    0;
    1    1    0;
    0    1    1;
    1    0    1];

cell_examples = cell_examples(cell_examples_IX);
cell_examples_dir = cell_examples_dir(cell_examples_IX);
for ii_cell = 1:length(cell_examples)
    axes(panel_A(ii_cell));
    hold on
    
    example_IX = find(strcmp(cells_130m_IDs,cell_examples{ii_cell}));
    cell = cells_130m(example_IX);
    
%     cell_ID = cell_examples{ii_cell};
%     cell = cell_load_data(cell_ID,'details','fields_per_win','FE','stats');

    ii_dir = cell_examples_dir(ii_cell);
    FE = cell.FE{ii_dir};
    fields_per_win = cell.fields_per_win{ii_dir};
    
    % trajectory + spikes + fields per win:
    for ii_FE = 1:length(FE)
%         plot(FE(ii_FE).pos, ii_FE*ones([1,length(FE(ii_FE).pos)]),'.','color',[0.93,0.93,0.93],'MarkerSize',1);
        plot([FE(ii_FE).pos(1) FE(ii_FE).pos(end)], [ii_FE ii_FE] ,'-','color',[0.93,0.93,0.93],'linewidth',0.75);
        plot(FE(ii_FE).spikes_pos, ii_FE*ones([1,length(FE(ii_FE).spikes_pos)]),'.k','MarkerSize',4);    
    end
    
    n_win = (1:length(fields_per_win(1).loc)) + 0.5;
    for ii_field = 1:length(fields_per_win)
        field_c = fields_colors(mod(ii_field-1,...
            size(fields_colors,1))+1, :);
        edges = [fields_per_win(ii_field).edges];
        plot(edges',[n_win; n_win],'linewidth',1,'color',field_c);
        start_flight = fields_per_win(ii_field).start_flight;
        
        % plot start/ end of field:
        if start_flight>2
            hax = gca;
            harr = annotation('arrow', [0 1],[0 1]);
            harr.Parent = hax;
            harr.LineStyle = 'none';
            harr.HeadLength = 8;
            harr.HeadWidth = 8;
            loc = fields_per_win(ii_field).CoM;
            loc(isnan(loc)) = [];
            harr.X = loc(1) + [0 0];
            harr.Y = start_flight + [-9 0];
            harr.Color = field_c;
            harr.HeadStyle = 'plain';
        end
        end_flight = fields_per_win(ii_field).end_flight;
        if end_flight<(length(FE))
            hax = gca;
            harr = annotation('arrow', [0 1],[0 1]);
            harr.Parent = hax;
            harr.LineStyle = 'none';
            harr.HeadLength = 8;
            harr.HeadWidth = 8;
            loc = fields_per_win(ii_field).CoM;
            loc(isnan(loc)) = [];
            harr.X = loc(end) + [0 0];
            harr.Y = end_flight + [9 0];
            harr.Color = field_c;
            harr.HeadStyle = 'plain';
        end
        
        % plot segment change events:
        [~, ~, seg_change_struct] = map_detect_segments_appear_disappear(fields_per_win,FE);
            
        seg_change_flight = [seg_change_struct(ii_field).seg_appear_flight; ...
            seg_change_struct(ii_field).seg_disappear_flight];
        seg_change_flight(seg_change_flight==1) = [];
        seg_change_flight(seg_change_flight==length(FE)) = [];
        
        for ii_seg = 1:length(seg_change_flight)
            hax = gca;
            harr = annotation('arrow', [0 1],[0 1]);
            harr.Parent = hax;
            harr.LineStyle = 'none';
            harr.HeadLength = 6;
            harr.HeadWidth = 6;
            loc = fields_per_win(ii_field).CoM;
            IX_valid_loc = find(~isnan(loc));
            [~, closest_loc_IX] = min(abs(IX_valid_loc - seg_change_flight(ii_seg)));
            closest_loc_IX = IX_valid_loc(closest_loc_IX);
            harr.X = loc(closest_loc_IX) + [-9 0];
            harr.Y = seg_change_flight(ii_seg) + [0 0];
            harr.Color = field_c;
            harr.HeadStyle = 'plain';
        end  
    end
    ylabel('Flight no.','Units','normalized','Position',[-0.08 0.5]);
    xlabel('Position (m)','Units','normalized','Position',[0.5 -0.1]);
    title(['Cell ' num2str(ii_cell)]);
%     title({cell_ID, ['Iso_dis = ' num2str(cell.stats.all.IsoDist)]},'Interpreter','none');
    h = gca;
    h.TickDir='out';
    h.TickLength = [0.015 0.015];
    h.XRuler.TickLabelGapMultiplier = -0.3;
    h.YRuler.TickLabelGapMultiplier = 0.001;
end
  


%% Panel B - model results:
% load('CA3CA1_FFModel_Perturbation_FieldSegmentNumbers_1.mat') ;
load(['CA3MECCA1_FFModel_' num2str(perturbation_prc) 'PC_Perturbation_FieldSegmentNumbers_2']);
% nFieldM_0 : number of fields in original (unperturbed) maps, multi  field CA3
% nFieldS_0 : number of fields in original (unperturbed) maps, single field CA3
% nFieldP_0 : number of fields in original (unperturbed) maps, periodic MEC
% nFieldM_sSeg : number of appeared/disappeared segments (1st, 2nd column), multi  field CA3
% nFieldS_sSeg : number of appeared/disappeared segments (1st, 2nd column), single field CA3
% nFieldP_sSeg : number of appeared/disappeared segments (1st, 2nd column),  periodic MEC
axes(panel_B(1)); hold on;
text(-0.3,1.25, 'B', 'Units','normalized','FontWeight','bold');
if perturbation_prc==5
    ylimits = [0,0.4];
elseif perturbation_prc==4
    ylimits = [0,0.4];
elseif perturbation_prc==3
    ylimits = [0,0.3];
elseif perturbation_prc==2
    ylimits = [0,0.2];
elseif perturbation_prc==1
    ylimits = [0,0.15];
end
        
msz = 6 ;
n_field_max = 10;
ik = find(nFieldS_0==1) ; 
% single field CA3 maps with one field

p1 = mean(mean(nFieldS_sSeg(ik,:)==1)); 
% p1 = mean(or(nFieldS_sSeg(ik,1)==1, nFieldS_sSeg(ik,2)==1)); 
% binomial distribution parameter 

plot(1:n_field_max, (binopdf(1,1:n_field_max,p1)).^2,'-k','LineWidth',4); hold on ;
plot(1:n_field_max, binopdf(2,1:n_field_max,p1)   ,'-' ,'Color',[0.6 0.6 0.6],'LineWidth',4) ; hold on ;

% OR version:
% plot(1:n_field_max, (2*binopdf(1,1:n_field_max,p1)).^2,'-k','LineWidth',4); hold on ;
% plot(1:n_field_max, 2*binopdf(2,1:n_field_max,p1)   ,'-' ,'Color',[0.6 0.6 0.6],'LineWidth',4) ; hold on ;


for k = 1:n_field_max
    ik = find(nFieldS_0==k) ;
    p1 = mean(mean(nFieldS_sSeg(ik,:)==1)) ;
    p2 = mean(mean(nFieldS_sSeg(ik,:)==2)) ;
    % OR version:
%     p1 = mean(or(nFieldS_sSeg(ik,1)==1,nFieldS_sSeg(ik,2)==1)) ;
%     p2 = mean(or(nFieldS_sSeg(ik,1)==2,nFieldS_sSeg(ik,2)==2)) ;
%     
    plot([k k],p1^2 + 2*p1*sqrt(p1*(1-p1))*[-1 1]/sqrt(2*length(ik)),'-k','LineWidth',2) ; hold on ;
    plot(k,p1^2,'ok','MarkerSize',msz,'MarkerFaceColor','k') ; hold on ;
    if k>1
        plot([k k],p2  +sqrt(p2*(1-p2))*[-1 1]/sqrt(2*length(ik)),'-k','LineWidth',2) ; hold on ;
        plot(k,p2  ,'ok','MarkerSize',msz,'MarkerFaceColor',[0.6 0.6 0.6]) ; hold on ;
    end
end
title({'Input model:','Single-field CA3'}) ;
xlabel({'No. of fields', 'in original map'},'Units','normalized','Position',[0.5 -0.12]);
ylabel('Probability', 'Units','normalized','Position',[-0.07 0.5])
ylim(ylimits) ;
box off;
h = gca;
h.TickDir='out';
h.YTick = [0, ylimits(2)];
h.TickLength = [0.03 0.03];
h.XRuler.TickLabelGapMultiplier = -0.3;
h.YRuler.TickLabelGapMultiplier = 0.001;


axes(panel_B(2));
ik = find(nFieldM_0==1) ; 
p1 = mean(mean(nFieldM_sSeg(ik,:)==1)); 
% binomial distribution parameter 
plot(1:n_field_max, (binopdf(1,1:n_field_max,p1)).^2,'-k','LineWidth',4); hold on ;
plot(1:n_field_max, binopdf(2,1:n_field_max,p1)   ,'-' ,'Color',[0.6 0.6 0.6],'LineWidth',4) ; hold on ;

for k = 1:n_field_max
    ik = find(nFieldM_0==k) ;
    p1 = mean(mean(nFieldM_sSeg(ik,:)==1)) ;
    p2 = mean(mean(nFieldM_sSeg(ik,:)==2)) ;
    % OR version:
%     p1 = mean(or(nFieldM_sSeg(ik,1)==1,nFieldM_sSeg(ik,2)==1)) ;
%     p2 = mean(or(nFieldM_sSeg(ik,1)==2,nFieldM_sSeg(ik,2)==2)) ;
    plot([k k],p1^2 + 2*p1*sqrt(p1*(1-p1))*[-1 1]/sqrt(2*length(ik)),'-k','LineWidth',2) ; hold on ;
    plot(k,p1^2,'ok','MarkerSize',msz,'MarkerFaceColor','k') ; hold on ;
    if k>1
        plot([k k],p2  +    sqrt(p2*(1-p2))*[-1 1]/sqrt(2*length(ik)),'-k','LineWidth',2) ; hold on ;
        plot(k,p2  ,'ok','MarkerSize',msz,'MarkerFaceColor',[0.6 0.6 0.6]) ; hold on ;
    end    
end
title({'Input model:','Multi-field CA3'}) ;
xlabel({'No. of fields', 'in original map'},'Units','normalized','Position',[0.5 -0.12]);
ylabel('Probability', 'Units','normalized','Position',[-0.07 0.5])
h = gca;
ylim(ylimits) ;
box off;
h.TickDir='out';
h.YTick = [0, ylimits(2)];
h.TickLength = [0.03 0.03];
h.XRuler.TickLabelGapMultiplier = -0.3;
h.YRuler.TickLabelGapMultiplier = 0.001;


axes(panel_B(3));
ik = find(nFieldP_0==1) ; 
p1 = mean(mean(nFieldP_sSeg(ik,:)==1)); 
% binomial distribution parameter 
plot(1:n_field_max, (binopdf(1,1:n_field_max,p1)).^2,'-k','LineWidth',4); hold on ;
plot(1:n_field_max, binopdf(2,1:n_field_max,p1)   ,'-' ,'Color',[0.6 0.6 0.6],'LineWidth',4) ; hold on ;

for k = 1:n_field_max
    ik = find(nFieldP_0==k) ;
    p1 = mean(mean(nFieldP_sSeg(ik,:)==1)) ;
    p2 = mean(mean(nFieldP_sSeg(ik,:)==2)) ;
    % OR version:
%     p1 = mean(or(nFieldM_sSeg(ik,1)==1,nFieldM_sSeg(ik,2)==1)) ;
%     p2 = mean(or(nFieldM_sSeg(ik,1)==2,nFieldM_sSeg(ik,2)==2)) ;
    plot([k k],p1^2 + 2*p1*sqrt(p1*(1-p1))*[-1 1]/sqrt(2*length(ik)),'-k','LineWidth',2) ; hold on ;
    plot(k,p1^2,'ok','MarkerSize',msz,'MarkerFaceColor','k') ; hold on ;
    if k>1
        plot([k k],p2  +    sqrt(p2*(1-p2))*[-1 1]/sqrt(2*length(ik)),'-k','LineWidth',2) ; hold on ;
        plot(k,p2  ,'ok','MarkerSize',msz,'MarkerFaceColor',[0.6 0.6 0.6]) ; hold on ;
    end    
end
title({'Input model:','Periodic MEC'}) ;
xlabel({'No. of fields', 'in original map'},'Units','normalized','Position',[0.5 -0.12]);
ylabel('Probability', 'Units','normalized','Position',[-0.07 0.5])
h = gca;
ylim(ylimits) ;
box off;
h.TickDir='out';
h.YTick = [0, ylimits(2)];
h.TickLength = [0.03 0.03];
h.XRuler.TickLabelGapMultiplier = -0.3;
h.YRuler.TickLabelGapMultiplier = 0.001;


%% load population data
% prm = PARAMS_GetAll();
% cells_t = DS_get_cells_summary();
% bats = [2382, 2311];
% cells_t(~ismember(cells_t.bat, bats ),:) = [];
% cells_t([cells_t.remove]==1,:) = [];%cells chosen to remove because of repititions between days
% cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',0);
% cells = [cells{:}];
% cells = [cells.details];
% cells(~contains({cells.brain_area}, 'CA1')) = [];
% cells(ismember([cells.DayNum], [0])) = []; %take only from day1 and on data
% % cells([cells.DayNum]>2) = []; %take only from day1 and day2
% % cells([cells.DayNum]<7) = []; %take only steady state (day 7 and on)
% cells(~ismember([cells.ClusterQuality], [2])) = [];
% % remove cells based on FR:
% cells = cellfun(@(c)(cell_load_data(c,'details','meanFR')), {cells.cell_ID}, 'UniformOutput',0);
% cells = [cells{:}];
% cells_details = [cells.details];
% cells_ID = {cells_details.cell_ID};
% meanFR = [cells.meanFR];
% cells_ID([meanFR.all]>prm.inclusion.interneuron_FR_thr)=[];
% % remove non-signif cells:
% cells = cellfun(@(c)(cell_load_data(c,'details','signif')), cells_ID, 'UniformOutput',0);
% cells = [cells{:}];
% signif = cat(1,cells.signif);
% signif = arrayfun(@(x)(x.TF), signif);
% signif = any(signif,2);
% cells_ID(~signif) = [];
% 
% clear cells stats cells_details cells_t

max_num_fields_in_map = n_field_max;
count_seg_change = [];
all_maps_field_num = [];
for ii_cell = 1:length(cells_130m)
    cell = cells_130m(ii_cell);
%     cell_ID = cells_ID{ii_cell};
%     cell = cell_load_data(cell_ID, 'stats','fields_per_win','FE');
    for ii_dir = 1:2
        field_num = cell.stats.dir(ii_dir).field_num;
        field_num(field_num>max_num_fields_in_map) = max_num_fields_in_map;
        map_signif = cell.stats.dir(ii_dir).map_signif;
        if field_num>0 && map_signif %&& field_num<6
            fields_per_win = cell.fields_per_win{ii_dir};
            FE = cell.FE{ii_dir};
            [segment_appear_flight, segment_disappear_flight, ~] = map_detect_segments_appear_disappear(fields_per_win,FE);
            count_seg_change = [count_seg_change ; [segment_appear_flight' segment_disappear_flight']];
            all_maps_field_num = [all_maps_field_num ; field_num*ones([length(segment_appear_flight),1])];
            
        end
    end
end


%% Panel C: data
axes(panel_C); hold on;
cla;
text(-0.3,1.25, 'C', 'Units','normalized','FontWeight','bold');
ylimits = [0,0.08];

ik = find(all_maps_field_num==1) ;
p1 = mean(mean(count_seg_change(ik,:)==1)) ; 
% binomial distribution parameter 

h1_theory = plot(1:n_field_max, (binopdf(1,1:n_field_max,p1)).^2,'-k','LineWidth',4) ; hold on ;
h2_theory = plot(1:n_field_max, binopdf(2,1:n_field_max,p1)   ,'-' ,'Color',[0.6 0.6 0.6],'LineWidth',4) ; hold on ;

% OR version:
% h1_theory = plot(1:n_field_max, (2*binopdf(1,1:n_field_max,p1)).^2,'-k','LineWidth',4) ; hold on ;
% h2_theory = plot(1:n_field_max, 2*binopdf(2,1:n_field_max,p1)   ,'-' ,'Color',[0.6 0.6 0.6],'LineWidth',4) ; hold on ;



for k = 1:n_field_max
    ik = find(all_maps_field_num==k) ;
    p1 = mean(mean(count_seg_change(ik,:)==1)) ;
    p2 = mean(mean(count_seg_change(ik,:)==2)) ;
    % OR version:
%     p1 = mean(or(count_seg_change(ik,1)==1,count_seg_change(ik,2)==1)) ;
%     p2 = mean(or(count_seg_change(ik,1)==2,count_seg_change(ik,2)==2)) ;
    
    plot([k k],p1^2 + 2*p1*sqrt(p1*(1-p1))*[-1 1]/sqrt(2*length(ik)),'-k','LineWidth',2) ; hold on ;
    h1 = plot(k,p1^2,'ok','MarkerSize',msz,'MarkerFaceColor','k') ; hold on ;
    if k>1
        plot([k k],p2  +     sqrt(p2*(1-p2))*[-1 1]/sqrt(2*length(ik)),'-k','LineWidth',2) ; hold on ;
        h2 = plot(k,p2  ,'ok','MarkerSize',msz,'MarkerFaceColor',[0.6 0.6 0.6]) ; hold on ;
    end
end
title('Data') ;
xlabel({'No. of fields', 'in original map'},'Units','normalized','Position',[0.5 -0.12]);
ylabel('Probability', 'Units','normalized','Position',[-0.07 0.5])
ylim(ylimits) ;
box off;
h = gca;
h.TickDir='out';
d = min(0.1,ylimits(2));
h.YTick = 0:d:ylimits(2);
h.TickLength = [0.03 0.03];
h.XRuler.TickLabelGapMultiplier = -0.3;
h.YRuler.TickLabelGapMultiplier = 0.001;

%% legend
axes(panel_legend);
hold on;
cla;
plot(1, 4,'ok','MarkerSize',msz,'MarkerFaceColor','k');
line([0.95 1.05], [3 3] , 'LineWidth',4, 'color', 'k');
plot(1, 2,'ok','MarkerSize',msz,'MarkerFaceColor',[0.6 0.6 0.6]);
line([0.95 1.05], [1 1], 'LineWidth',4, 'color', [0.6 0.6 0.6]);

text(1.1, 4, 'p_1^2', 'fontsize', 7);
text(1.1, 3, 'Theoretical p_1^2', 'fontsize', 7);
text(1.1, 2, 'p_2', 'fontsize', 7);
text(1.1, 1, 'Theoretical p_2', 'fontsize', 7);

% h_leg = legend([h1 h1_theory h2 h2_theory],'p_1^2','Theoretical p_1^2','p_2','Theoretical p_2','position',[0.87 0.5 0.05 0.02]);

xlim([0.9, 1.5]);
ylim([0.5 4.5]);
set(gca,'Visible','off');

%% save figure
fig_name_out = fullfile(res_dir, sprintf('%s',fig_name_str));
% fig_name_out = fullfile(res_dir, sprintf('%s_%dPRC_new',fig_name_str,perturbation_prc));
% fig_name_out = fullfile(res_dir, sprintf('%s_ex%d_%d',fig_name_str,cell_examples_IX(1),cell_examples_IX(2)));
% fig_name_out = fullfile(res_dir, sprintf('%s__corr_%s_%d_paramset_%d',fig_name_str,corr_type,field_speed_opt,prm.parmaset));
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');
