%% Large Scale - Fig. S12 - compartmentalization

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S12';
fig_caption_str = 'compartmentalization';
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
disp([fig_name_str ': ' fig_caption_str]);
disp('======================================================');
disp('');

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
% panels_pos_x = [2 6.5 11 15.5];
% panels_pos_y = [22 18 14 10 6 2];
% panels_size = [2.5 2.5];
% for ii_model = 1:6
%     for ii_feature = 1:4
%         panels(ii_model,ii_feature) = axes('position', [panels_pos_x(ii_feature) panels_pos_y(ii_model) panels_size]);
%     end
% end

% create panels
panels_pos_x = [3 8 13]+1.5;
panels_pos_y = [20 15]-2;
panels_size = [3 3];
panels=[];
panels(1,1) = axes('position', [panels_pos_x(1) panels_pos_y(1) panels_size]);
panels(2,1) = axes('position', [panels_pos_x(1) panels_pos_y(2) panels_size]);
panels(3,1) = axes('position', [panels_pos_x(2) panels_pos_y(1) panels_size]);
panels(4,1) = axes('position', [panels_pos_x(2) panels_pos_y(2) panels_size]);
panels(5,1) = axes('position', [panels_pos_x(3) panels_pos_y(1) panels_size]);
panels(6,1) = axes('position', [panels_pos_x(3) panels_pos_y(2) panels_size]);
panel_legend = axes('position', [17.5 20  1 0.8]);

%% load data
load('L:\processed_data_structs\cells_bat_200m.mat');
% load('L:\Yonatan_theory\20200806__model_ratio_LS\SegmentSmallScaleModel_FieldStatistics.mat');
% model_data_files = {
%     "L:\Yonatan_theory\20201201__CompartmentModelResults_3versions\SegmentSmallScaleModel_6mFields_FieldStatistics.mat",...
%     "L:\Yonatan_theory\20201201__CompartmentModelResults_3versions\SegmentSmallScaleModel_6mScaledInterLandmark_FieldStatistics.mat",...
%     "L:\Yonatan_theory\20201201__CompartmentModelResults_3versions\SegmentSmallScaleModel_200mSmallestField_FieldStatistics.mat"
%     };
model_data_files = {
"L:\Yonatan_theory\20201214__CompartmentModelResults_6versions_dist2LM\SegmentSmallScaleModel_FieldStatistics.mat"
"L:\Yonatan_theory\20201214__CompartmentModelResults_6versions_dist2LM\SegmentSmallScaleModel_FieldStatistics_nField.mat"
"L:\Yonatan_theory\20201214__CompartmentModelResults_6versions_dist2LM\SegmentSmallScaleModel_InterLandmark_FieldStatistics.mat"
"L:\Yonatan_theory\20201214__CompartmentModelResults_6versions_dist2LM\SegmentSmallScaleModel_InterLandmark_FieldStatistics_nField.mat"
"L:\Yonatan_theory\20201214__CompartmentModelResults_6versions_dist2LM\SegmentSmallScaleModel_SmallestField_FieldStatistics.mat"
"L:\Yonatan_theory\20201214__CompartmentModelResults_6versions_dist2LM\SegmentSmallScaleModel_SmallestField_FieldStatistics_nField.mat"
};
model_name_str = {
    {"6m fields";"match coverage"}
    {"6m fields";"match no. fields"}
    {"6m fields (scaled)";"match coverage"}
    {"6m fields (scaled)";"match no. fields"}
    {"smallest 200m fields";"match coverage"}
    {"smallest 200m fields";"match no. fields"}
    };

%% arrange data from 200m
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

field_num = [stats_dir(signif).field_num];
field_size = [fields.width_prc];
ratio_LS = [stats_all.field_ratio_LS];
ratio_LS(isnan(ratio_LS))=[];
ratio_LS_with_1s = [stats_all.field_ratio_LS];
ratio_LS_with_1s(isnan(ratio_LS_with_1s))=1;
    
%% run over models
for ii_model = 1:length(model_data_files)

    %%
    clear FieldSize_PoissonGamma MaxMinRatio_PoissonGamma fSEG_PoissonGamma nField_PoissonGamma NearestLandmark_PoissonGamma
    load(model_data_files{ii_model});
    
%% Fields size hist
axes(panels(ii_model,1));
cla('reset')
hold on
% text(-0.4,1.15, 'A', 'Units','normalized','FontWeight','bold');

bin_size = 0.2;
x1 = field_size;
x2 = FieldSize_PoissonGamma.*bin_size;

h = histogram(x1);
h.BinWidth = 1;
% h.Normalization = 'pdf';
h.Normalization = 'probability';
h.FaceColor = 0.5*[1 1 1];

h = histogram(x2);
h.BinWidth = 0.6;
% h.Normalization = 'pdf';
h.Normalization = 'probability';
h.DisplayStyle = 'stairs';
h.EdgeColor = 'c';
h.LineWidth=1.5;

[~,pval,ks2stat] = kstest2(x1,x2);
if pval==0
    pval = realmin;
end
text(0.35,1, ['P_{KS} = ' sprintf('%.2g',pval)], 'Units','normalized','FontSize',8);

ha=gca;
ha.XLim = [0 35];
ha.YLim = [5e-4 1e0];
ha.XTick = [0:10:40];
ha.YTick = 10.^[-5:0];
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.35;
ha.YRuler.TickLabelGapMultiplier = 0.1;
ha.YScale = 'log';
xlabel('Field size (m)', 'Units','normalized','Position',[0.5 -0.13]);
ylabel('Fraction of fields')
% title(model_name_str{ii_model})

%% Fields size vs. dist2LM
% axes(panels(ii_model,4));
% cla('reset')
% hold on
% text(-0.42,1.15, 'D', 'Units','normalized','FontWeight','bold');
% 
% x = NearestLandmark_PoissonGamma * bin_size;
% y = FieldSize_PoissonGamma * bin_size;
% 
% plot(x,y,'.')
% % [rho,rhopval] = corr(x,y,'type','Spearman');
% % text(0.5,1, sprintf('%s = %.2f,  P_{%s} = %.2f', '/rho',rho,'/rho',rhopval), 'Units','normalized','FontSize',8);
% 
% ha=gca;
% % ha.XLim = [1 27];
% % ha.YLim = [5e-3 5e-1];
% % ha.XTick = [1 2 5 10 20];
% % ha.YTick = 10.^[-5:0];
% ha.TickDir='out';
% ha.TickLength = [0.03 0.03];
% ha.XRuler.TickLabelGapMultiplier = -0.35;
% ha.YRuler.TickLabelGapMultiplier = 0.1;
% xlabel({'Distance of fields';'to nearest landmark (m)'});
% ylabel('Field size (m)','Units','normalized');

%% Fields size ratio (max/min) hist
% axes(panels(ii_model,2));
% cla('reset')
% hold on
% text(-0.42,1.15, 'B', 'Units','normalized','FontWeight','bold');
% 
% x1 = ratio_LS;
% x2 = MaxMinRatio_PoissonGamma;
% 
% nBinEdges = 9;
% edges = logspace(0,log10(25),nBinEdges);
% h = histogram(x1);
% h.BinEdges = edges;
% h.Normalization = 'probability';
% h.FaceColor = 0.5*[1 1 1];
% 
% h = histogram(x2);
% h.BinWidth = 0.53;
% h.Normalization = 'probability';
% h.DisplayStyle = 'stairs';
% h.EdgeColor = 'c';
% h.LineWidth=1.5;
% 
% [~,pval,ks2stat] = kstest2(x1,x2);
% text(0.5,1, ['P_{KS} = ' sprintf('%.2g',pval)], 'Units','normalized','FontSize',8);
% 
% ha=gca;
% ha.XLim = [1 27];
% ha.YLim = [5e-3 5e-1];
% ha.XTick = [1 2 5 10 20];
% ha.YTick = 10.^[-5:0];
% ha.TickDir='out';
% ha.TickLength = [0.03 0.03];
% ha.XRuler.TickLabelGapMultiplier = -0.35;
% ha.YRuler.TickLabelGapMultiplier = 0.1;
% ha.XScale = 'log';
% ha.YScale = 'log';
% xlabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[0.5 -0.17]);
% ylabel('Fraction of cells')

%% Number of fields
% axes(panels(ii_model,3));
% cla('reset')
% hold on
% text(-0.4,1.15, 'C', 'Units','normalized','FontWeight','bold');
% 
% bin_size = 0.2;
% x1 = field_num;
% x2 = nField_PoissonGamma;
% 
% h = histogram(x1);
% h.BinWidth = 1;
% % h.Normalization = 'pdf';
% h.Normalization = 'probability';
% h.FaceColor = 0.5*[1 1 1];
% 
% h = histogram(x2);
% h.BinWidth = 1;
% % h.Normalization = 'pdf';
% h.Normalization = 'probability';
% h.DisplayStyle = 'stairs';
% h.EdgeColor = 'c';
% h.LineWidth=1.5;
% 
% [~,pval,ks2stat] = kstest2(x1,x2);
% if pval==0
%     pval = 1e-300;
% end
% text(0.35,1, ['P_{KS} = ' sprintf('%.2g',pval)], 'Units','normalized','FontSize',8);
% 
% ha=gca;
% ha.XLim = [0 35];
% ha.YLim = [5e-4 1e0];
% ha.XTick = [0:10:40];
% ha.YTick = 10.^[-5:0];
% ha.TickDir='out';
% ha.TickLength = [0.03 0.03];
% ha.XRuler.TickLabelGapMultiplier = -0.35;
% ha.YRuler.TickLabelGapMultiplier = 0.1;
% ha.YScale = 'log';
% xlabel('No. of fields', 'Units','normalized','Position',[0.5 -0.13]);
% ylabel('Fraction of cells')

    %%
%     axes(panels(ii_model,end));
%     text(1.1,0.5,model_name_str{ii_model},'Units','normalized');

end % models

%% add titles to rows/columns
axes(panels(1,1));
text(-0.85,0.5,{'Matched';'coverage'},'Units','normalized','FontSize',9,'FontWeight','bold','HorizontalAlignment','center');
axes(panels(2,1));
text(-0.85,0.5,{'Matched';'no. of fields'},'Units','normalized','FontSize',9,'FontWeight','bold','HorizontalAlignment','center');

axes(panels(3,1));
text(0.5,1.8,{'Gamma distribution fitted to:'},'Units','normalized','FontSize',9,'FontWeight','bold','HorizontalAlignment','center');

axes(panels(1,1));
text(0.5,1.4,{'Field sizes in 6m tunnel'},'Units','normalized','FontSize',9,'FontWeight','bold','HorizontalAlignment','center');
axes(panels(3,1));
text(0.5,1.4,{'Field sizes in 6m tunnel';'scaled to median';'inter-landmark interval'},'Units','normalized','FontSize',9,'FontWeight','bold','HorizontalAlignment','center');
axes(panels(5,1));
text(0.5,1.4,{'Smallest field per neuron';'in 200m tunnel'},'Units','normalized','FontSize',9,'FontWeight','bold','HorizontalAlignment','center');

%% add legend
axes(panel_legend);
cla('reset');
hold on
xlim([0 4])
ylim([0 4])
axis equal
% hp=patch([1 1 2 2], 2*[1 1 1 1]+.3*[-1 1 1 -1], 0.5*[1 1 1],'EdgeColor','k');
hp=patch([1 1 2 2], [1 2 2 1]+2, 0.5*[1 1 1],'EdgeColor','k');
hp.Clipping='off';
plot([1 2],      1*[1 1], 'c','LineWidth',2,'Clipping','off');
text(2.6, 3.6, 'Data','FontSize',7,'HorizontalAlignment','left');
text(2.6, 0.5, {'Small-scale';'compartments model'},'FontSize',7,'HorizontalAlignment','left');
hax=gca;
hax.Visible='off';


%% print/save the figure
fig_name_out = fullfile(res_dir, [fig_name_str]);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');



%%





