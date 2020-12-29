%% Large Scale - Fig. 8 - Theoretical analysis (mechanism)

%%
clear 
clc

%%
paper_fig_8_arrange_sim_data
data = paper_fig_8_arrange_real_data;


%% graphical options
multiple_panels = 1;

%% choose examples 
% cell_examples_IX = [29 198  12 121 74];
% cell_examples_IX = [29 206 1085 1094 1073];
% cell_examples_IX = [1097 1099 1100 1101 1104];
% cell_examples_IX = [1112 1135 1139 1180 1183];
% cell_examples_IX = [1278 1280 1292 1301 1306];
% cell_examples_IX = [1311 1312 1326 1344 1348];

% cell_examples_IX = [29 3988 74 121 54];
% cell_examples_IX = [29 3988 74 1100 54];
cell_examples_IX = [29 3988 74 391 54]; % chosen examples
% cell_examples_IX = [29 356 74 1100 391];
% cell_examples_IX = [29 734 1016 1768 1847];
% cell_examples_IX = [1948 1969 2000 2024 2030 ];

% good many fields / large ratio:
% [54 121 391 734 1016 1768 1847 1948 1969 2000 2024 2030 3988];


%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'Fig_8';
fig_caption_str = 'Theoretical analysis - Mechanism';
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
panel_DE_size = [2.5 2.5];
panel_F_size = [2.3 2.5];
panel_D(1) = axes('position', [ 2.0 21.5 panel_DE_size]);
panel_D(2) = axes('position', [ 5.7 21.5 panel_DE_size]);
panel_D(3) = axes('position', [ 9.4 21.5 panel_DE_size]);
panel_E(1) = axes('position', [ 2.0 17.3 panel_DE_size]);
panel_E(2) = axes('position', [ 5.7 17.3 panel_DE_size]);
panel_E(3) = axes('position', [ 9.4 17.3 panel_DE_size]);
panel_F(1) = axes('position', [ 2.0 12.5 panel_F_size.*[0.75 0.85]]);
panel_F(2) = axes('position', [ 4.4 12.5 panel_F_size.*[0.85 0.85]]);
panel_F(3) = axes('position', [ 7.0 12.5 panel_F_size.*[0.85 0.85]]);
panel_F(4) = axes('position', [ 2.4 15.0 0.35 0.85]);
panel_G(1) = axes('position', [15.0 21.1 5 3]);
panel_G(2) = axes('position', [17.2 22.4 2.5 1.7]);
panel_H(1) = axes('position', [15.0 17.7 5 2]);
panel_I(1) = axes('position', [10.5 12.5 5 3]);
panel_I(2) = axes('position', [13.2 14.0 2.5 1.5]);
panel_J(1) = axes('position', [17.2 12.5 panel_F_size.*[1.25 1]]);

% Attractor panels
panel_A = axes('position', [ 2 7.3 11 3]);
panel_B = [];
if multiple_panels
    y_pos = linspace(2,6,5)+0.5;
    for ii = 1:5
        panel_B(ii) = axes('position', [ 2 y_pos(ii) 11 0.85*mean(diff(y_pos))]);
    end
else
    panel_B = axes('position', [ 2 2 12 6]);
end
panel_B = flip(panel_B);
panel_C(1) = axes('position', [ 14.5 7.5 2.5 2.5]);
panel_C(2) = axes('position', [ 18.3 7.5 2.5 2.5]);
panel_C(3) = axes('position', [ 14.5 3.5 2.5 2.5]);

FF_offset = -10.5;
CANN_offset = +14.3;
FF_panels = [panel_D(:);panel_E(:);panel_F(:);panel_J(:);panel_G(:);panel_H(:);panel_I(:)];
CANN_panels = [panel_A(:);panel_B(:);panel_C(:)];
for ii_panel = 1:length(FF_panels)
    FF_panels(ii_panel).Position(2) = FF_panels(ii_panel).Position(2) + FF_offset;
end
for ii_panel = 1:length(CANN_panels)
    CANN_panels(ii_panel).Position(2) = CANN_panels(ii_panel).Position(2) + CANN_offset;
end

colormap gray
lw = 2;

%%
FF_title_pos = 0.56;
CANN_title_pos = 0.95;
ht=annotation('textbox', [0.3 FF_title_pos .4 0], 'String','Feedforward model', ...
    'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','none', 'FitBoxToText','off');
ht.LineStyle='none';
ht.FontSize = 10;
ht.FontWeight = 'bold';
ht=annotation('textbox', [0.3 CANN_title_pos .4 0], 'String','Attractor network model', ...
    'HorizontalAlignment','center','VerticalAlignment','middle','Interpreter','none', 'FitBoxToText','off');
ht.LineStyle='none';
ht.FontSize = 10;
ht.FontWeight = 'bold';

hl = annotation('line',[.05 .4],FF_title_pos*[1 1]);
hl = annotation('line',[.6 .95],FF_title_pos*[1 1]);
hl = annotation('line',[.05 .4],CANN_title_pos*[1 1]);
hl = annotation('line',[.6 .95],CANN_title_pos*[1 1]);


%% panels A+B CA3 - input/output cell examples
axes(panel_D(1)); cla('reset'); hold on;
text(-0.3,1.3, 'D', 'Units','normalized','FontWeight','bold');
text(1.24,1.27, 'Input from CA3', 'Units','normalized','FontWeight','bold','HorizontalAlignment','center','FontSize',8);
IX = round(linspace(1,Nmax,nplot1));
IX(5) = 888;
IX(6) = 1116;
imagesc(xFF,1:nplot1,1-fS(:,IX)') ;
ylabel('Input neuron no.','Units','normalized','Position',[-0.18 0.5]);
xlabel('Position (m)','Units','normalized','Position',[0.5 -0.07]);
title('Single-field CA3','FontWeight','bold','Units','normalized','Position',[0.5 1.05]) ;
h=gca;
h.Box = 'on';
h.XTick=[0 200];
h.YTick=[1 5 10];
h.XLim = [0 200];
h.YLim = [0.5 10.5];
h.XAxis.TickLength(1) = 0.03;
h.XRuler.TickLabelGapOffset = -1;
h.YRuler.TickLabelGapOffset = 0.5;

axes(panel_D(2)); cla; hold on;
IX = round(linspace(1,Nmax,nplot1));
IX(1) = 2;
imagesc(xFF,1:nplot1,1-fM(:,IX)') ;
xlabel('Position (m)','Units','normalized','Position',[0.5 -0.07]) ;
ylabel('Input neuron no.','Units','normalized','Position',[-0.18 0.5]);
title('Multi-field CA3','FontWeight','bold','Units','normalized','Position',[0.5 1.05]) ;
h=gca;
h.Box = 'on';
h.XTick=[0 200];
h.YTick=[1 5 10];
h.XLim = [0 200];
h.YLim = [0.5 10.5];
h.XAxis.TickLength(1) = 0.03;
h.XRuler.TickLabelGapOffset = -1;
h.YRuler.TickLabelGapOffset = 0.5;

axes(panel_E(1)); cla; hold on;
text(-0.3,1.15, 'E', 'Units','normalized','FontWeight','bold');
text(.5,1.11, 'CA1 output', 'Units','normalized','FontWeight','bold','HorizontalAlignment','center','FontSize',8);
IX = round(linspace(1,rep,nplot1));
IX([1 4 6]) = IX([1 4 6]) + [1 6 6];
imagesc(xFF,1:nplot1,1-fResS(:,IX)') ;
ylabel('Output neuron no.','Units','normalized','Position',[-0.18 0.5]);
xlabel('Position (m)','Units','normalized','Position',[0.5 -0.07]) ;
h=gca;
h.Box = 'on';
h.XTick=[0 200];
h.YTick=[1 5 10];
h.XLim = [0 200];
h.YLim = [0.5 10.5];
h.XAxis.TickLength(1) = 0.03;
h.XRuler.TickLabelGapOffset = -1;
h.YRuler.TickLabelGapOffset = 0.5;

axes(panel_E(2)); cla; hold on;
text(.5,1.11, 'CA1 output', 'Units','normalized','FontWeight','bold','HorizontalAlignment','center','FontSize',8);
IX = round(linspace(1,rep,nplot1));
IX([1 4 6]) = IX([1 4 6]) + [1 6 6];
rng(3);
IX = IX(randperm(length(IX)));
imagesc(xFF,1:nplot1,1-fResM(:,IX)') ;
xlabel('Position (m)','Units','normalized','Position',[0.5 -0.07]) ;
ylabel('Output neuron no.','Units','normalized','Position',[-0.18 0.5]);
h=gca;
h.Box = 'on';
h.XTick=[0 200];
h.YTick=[1 5 10];
h.XLim = [0 200];
h.YLim = [0.5 10.5];
h.XAxis.TickLength(1) = 0.03;
h.XRuler.TickLabelGapOffset = -1;
h.YRuler.TickLabelGapOffset = 0.5;

%% MEC input cell examples
axes(panel_D(3)); cla; hold on;
iang = 1;
IX = round(linspace(1,Nmax,nplot2));
IX(1)=IX(1)+6;
IX(3)=IX(3)+1;
IX(4)=IX(4)+7;
IX(10)=IX(10)+2;
IX(11)=IX(11)+14;
IX(15)=IX(15)+1;
IX(16)=IX(16)+2;
IX(18)=IX(18)+1;
IX(19)=IX(19)+4;
IX(20)=IX(20)-6;
imagesc(xFF,1:length(IX),1-fP(:,IX,iang)') ; hold on ;
for i = 1:naMEC-1
    plot([0 L]*ds/100,[1 1]*i*nplot2/naMEC+0.5,'-','Color',[0.6 0.6 0.6]) ;
end
for i = 1:naMEC
    y = 2.5 + (i-1)*nplot2/naMEC;
    text(210, y, "M"+i, 'FontSize',7);
end
ylabel('Input neuron no.','Units','normalized','Position',[-0.18 0.5]);
% xlabel('Position (m)','Units','normalized','Position',[0.5 -0.04]) ;
xlabel('Position (m)','Units','normalized','Position',[0.5 -0.07]) ;
title({'Input from MEC';'Periodic grid-cells'},'FontWeight','bold','Units','normalized','Position',[0.5 1.05]) ;
h=gca;
h.Box = 'on';
h.XTick=[0 200];
h.YTick=[1 5 10 15 20];
h.XLim = [0 200];
h.YLim = [0.5 20.5];
h.XAxis.TickLength(1) = 0.03;
h.XRuler.TickLabelGapOffset = -1;
h.YRuler.TickLabelGapOffset = 0.5;

%% MEC output cell examples
axes(panel_E(3)); cla; hold on;
text(.5,1.11, 'CA1 output', 'Units','normalized','FontWeight','bold','HorizontalAlignment','center','FontSize',8);
load(fullfile(sim_data_dir,['CA3MEC_FFModel_fResXw_iangle_' num2str(iang) '.mat']));
load(fullfile(sim_data_dir,['CA3MEC_FFModel_fResXs_iangle_' num2str(iang) '.mat']));
i3 = 1;
% find some periodic output examples
[~,IX] = sort(max(diff(abs(fft(fResXs(:,:,i3),[],1)).^2)),'descend');
IX = IX([2 3 5 8 51:56]);
IX([4 10]) = IX([10 4]);
% imagesc(x,1:nplot1,1-fResXs(:,round(linspace(1,rep,nplot1)),i3)') ;
imagesc(xFF,1:length(IX),1-fResXs(:,IX,i3)') ;
ylabel('Output neuron no.','Units','normalized','Position',[-0.18 0.5]);
xlabel('Position (m)','Units','normalized','Position',[0.5 -0.07]) ;
h=gca;
h.Box = 'on';
h.XTick=[0 200];
h.YTick=[1 5 10 15 20];
h.XLim = [0 200];
h.YLim = [0.5 10.5];
h.XAxis.TickLength(1) = 0.03;
h.XRuler.TickLabelGapOffset = -1;
h.YRuler.TickLabelGapOffset = 0.5;


%% add arrow from input -> output
for ii_arrow = 1:length(panel_D)  
    har=annotation('arrow');
    har.Units='centimeters';
    har.Position = [panel_D(ii_arrow).Position(1) + 0.5  * panel_D(ii_arrow).Position(3),...
                    panel_D(ii_arrow).Position(2) - 0.28 * panel_D(ii_arrow).Position(4),...
                    0 -0.5];
    har.HeadWidth = 5;
    har.HeadLength = 5;
    har.LineWidth = 1.5;
end

%% Distributions (num fields / field size / ratio)
lw = 2;
axes(panel_F(1)); cla; hold on;
xlimits = [0 25];
% xlimits = [0 33];
text(-0.45,1.62, 'F', 'Units','normalized','FontWeight','bold');
plot(xnField,pnFieldS,'-','Color',clrM(1,:),'lineWidth',lw) ; hold on ;
plot(xnField,pnFieldM,'-','Color',clrM(2,:),'lineWidth',lw) ; hold on ;
plot(xnField,squeeze(pnFieldXs(1,1,:))','-','Color',clrM(3,:),'lineWidth',lw) ; hold on ;
h=histogram(data.field_num);
h.DisplayStyle = 'stairs';
% h.BinEdges = xnField;
h.BinWidth = 1;
% h.Data(h.Data > xlimits(2)) = xlimits(2);
% h.Normalization = 'pdf';
h.Normalization = 'probability';
h.EdgeColor = 'k';
hax=gca;
hax.XScale = 'linear';
hax.YScale = 'log';
hax.XLim = xlimits+[0 1];
hax.YLim = [2e-3 0.23];
hax.XTick = [0:10:40];
hax.YTick = 10.^[-4 -3 -2 -1];
hax.XAxis.TickLength(1) = 0.025;
hax.YAxis.TickLength(1) = 0.025;
hax.YMinorTick = 'on';
hax.XRuler.TickLabelGapOffset = -1;
hax.YRuler.TickLabelGapOffset = 0;
xlabel({'No. of fields';'per direction'},'Units','normalized','Position',[0.5 -0.18]);
% ylabel({'Probability';'density function'},'Units','normalized','Position',[-0.45 0.5]);
ylabel('Probability','Units','normalized','Position',[-0.4 0.5]);

axes(panel_F(2)); cla; hold on;
xlimits = [0 35];
plot(xField,pFieldSizeS,'-','Color',clrM(1,:),'lineWidth',lw) ; hold on ;
plot(xField,pFieldSizeM,'-','Color',clrM(2,:),'lineWidth',lw) ; hold on ;
plot(xField,squeeze(pFieldSizeXs(1,1,:))','-','Color',clrM(3,:),'lineWidth',lw) ; hold on ;
h=histogram(data.field_size);
h.DisplayStyle = 'stairs';
% h.BinEdges = xField;
h.BinWidth = 1;
h.Data(h.Data > xlimits(2)) = xlimits(2);
% h.Normalization = 'pdf';
h.Normalization = 'probability';
h.EdgeColor = 'k';
hax=gca;
hax.YScale = 'log';
hax.XLim = xlimits;
hax.YLim = [3e-4 0.22];
hax.XTick = [0:10:40];
hax.YTick = 10.^[-4 -3 -2 -1];
hax.XAxis.TickLength(1) = 0.025;
hax.YAxis.TickLength(1) = 0.025;
hax.YMinorTick = 'on';
hax.XRuler.TickLabelGapOffset = -1;
hax.YRuler.TickLabelGapOffset = 0;
xlabel('Field size (m)','Units','normalized','Position',[0.5 -0.18]);
% ylabel('PDF','Units','normalized','Position',[-0.3 0.5]);

axes(panel_F(3)); cla; hold on;
% xlimits = [1 27];
xlimits = [1 100];
plot(xMaxMinRatio,pMaxMinRatioS,'-','Color',clrM(1,:),'lineWidth',lw) ; hold on ;
plot(xMaxMinRatio,pMaxMinRatioM,'-','Color',clrM(2,:),'lineWidth',lw) ; hold on ;
plot(xMaxMinRatio,squeeze(pMaxMinRatioXs(1,1,:))','-','Color',clrM(3,:),'lineWidth',lw) ; hold on ;
nBinEdges = 9;
edges = logspace(0,log10(25),nBinEdges);
h=histogram(data.ratio_LS);
h.DisplayStyle = 'stairs';
% h.DisplayStyle = 'bar';
% h.BinEdges = xnField;
h.BinEdges = edges;
h.Data(h.Data > xlimits(2)) = xlimits(2);
% h.Normalization = 'pdf';
h.Normalization = 'probability';
h.EdgeColor = 'k';
% h.FaceColor = 0.5*[1 1 1];
hax=gca;
hax.XScale = 'log';
hax.YScale = 'log';
hax.XLim = xlimits;
hax.YLim = [5e-3 0.3];
% hax.XTick = [0:10:40];
% hax.XTick = [1 2 5 10 20];
hax.XTick = [1 10 100];
hax.YTick = 10.^[-4 -3 -2 -1];
hax.XAxis.TickLength(1) = 0.025;
hax.YAxis.TickLength(1) = 0.025;
hax.YMinorTick = 'on';
hax.XRuler.TickLabelGapOffset = -1.5;
hax.YRuler.TickLabelGapOffset = 0;
xlabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[0.5 -0.2]);
% ylabel('PDF','Units','normalized','Position',[-0.3 0.5]);

% legend
axes(panel_F(4)); cla; hold on;
axis off
plot([0 1],[3 3],'-','Color',clrM(1,:),'lineWidth',lw) ; hold on ;
plot([0 1],[2 2],'-','Color',clrM(2,:),'lineWidth',lw) ; hold on ;
plot([0 1],[1 1],'-','Color',clrM(3,:),'lineWidth',lw) ; hold on ;
plot([0 1],[0 0],'-','Color','k','lineWidth',lw) ; hold on ;
text(1.2, 3, 'Single-field CA3 model', 'FontSize',7);
text(1.2, 2, 'Multi-field CA3 model', 'FontSize',7);
text(1.2, 1, 'Periodic MEC model', 'FontSize',7);
text(1.2, 0, 'Data', 'FontSize',7);

%% field dynamics
axes(panel_J(1)); cla; hold on;
text(-0.25,1.37, 'J', 'Units','normalized','FontWeight','bold');

% load data
ii_prc = 4;
filename = "fig_8_dynamics_panel_data_PRC_" + ii_prc;
load(fullfile(res_dir,filename));

% plot
x = fig_8_dynamics_panel_data.x;
y = fig_8_dynamics_panel_data.y;
err = fig_8_dynamics_panel_data.err;
hb = bar(x,y);
hb.FaceColor = 'flat';
hb.CData = [[0 0 0]; [0.6 0.6 0.6]; [0 0 0]; [0.6 0.6 0.6]; [0 0 0]; [0.6 0.6 0.6]];
% hb.CData(1:2:5,:) = 1.0.*clr(1:3,:);
% hb.CData(2:2:6,:) = 0.5.*clr(1:3,:);
he=errorbar(x,y,err);
he.Clipping = 'off';
he.CapSize = 2;
he.LineWidth = 1;
he.LineStyle = 'none';
he.Color = 'k';
ylimits = [-1 1];
h=gca;
% h.XTick = x;
% h.XTickLabels = {'MEC','model', 'Multi-field','CA3 model', 'Single-field','CA3 model'};
% h.XTickLabelRotation = 45;
% h.XTick=[1 2 3];
xtick = mean(reshape(x,2,[]));
% text(xtick, repelem(ylimits(1)-0.25*range(ylimits),length(xtick)),...
%     {{'MEC';'model'}, {'Multi-field';'CA3 model'}, {'Single-field';'CA3 model'}},...
%     'Rotation',45, 'Fontsize',6,'HorizontalAlignment','Center');

% text(xtick-0.2, repelem(ylimits(1)-0.23*range(ylimits),length(xtick)),...
%     {{'MEC'}, {'Multi-field'}, {'Single-field'}},...
%     'Rotation',45, 'Fontsize',7,'HorizontalAlignment','Center');
% text(xtick+0.2, repelem(ylimits(1)-0.32*range(ylimits),length(xtick)),...
%     {{'model'}, {'CA3 model'}, {'CA3 model'}},...
%     'Rotation',45, 'Fontsize',7,'HorizontalAlignment','Center');

text(1.1-0.2,-1.45, 'Single-field', 'Rotation',45, 'Fontsize',7,'HorizontalAlignment','Center');
text(1.8-0.2,-1.6, 'CA3 model', 'Rotation',45, 'Fontsize',7,'HorizontalAlignment','Center');
text(4.1,-1.45, 'Multi-field', 'Rotation',45, 'Fontsize',7,'HorizontalAlignment','Center');
text(4.8,-1.6, 'CA3 model', 'Rotation',45, 'Fontsize',7,'HorizontalAlignment','Center');
text(7.5+0.2,-1.3, 'MEC', 'Rotation',45, 'Fontsize',7,'HorizontalAlignment','Center');
text(8.1+0.2,-1.43, 'model', 'Rotation',45, 'Fontsize',7,'HorizontalAlignment','Center');

% add signif lines and aestricks
lw=1;
plot(xtick(1)+[-0.6 0.6], [1 1 ]*1.1, '-k', 'LineWidth',lw, 'Clipping','off');
plot(xtick(2)+[-0.6 0.6], [1 1 ]*0.5, '-k', 'LineWidth',lw, 'Clipping','off');
plot(xtick(3)+[-0.6 0.6], [1 1 ]*0.5, '-k', 'LineWidth',lw, 'Clipping','off');
plot(xtick([1 2])+[0.15 0], [1 1]*1.25, '-k', 'LineWidth',lw, 'Clipping','off');
plot(xtick([1 3])-[0.15 0], [1 1]*1.50, '-k', 'LineWidth',lw, 'Clipping','off');
plot(xtick([2 2]), [0.5 1.25], '-k', 'LineWidth',lw, 'Clipping','off');
plot(xtick([3 3]), [0.5 1.50], '-k', 'LineWidth',lw, 'Clipping','off');
plot(xtick([1 1])-0.15, [1.1 1.50], '-k', 'LineWidth',lw, 'Clipping','off');
plot(xtick([1 1])+0.15, [1.1 1.25], '-k', 'LineWidth',lw, 'Clipping','off');
text([1 1]*mean(xtick([1 2])),[1 1]*1.32, '*****','fontSize',11,'HorizontalAlignment','center');
text([1 1]*mean(xtick([1 3])),[1 1]*1.58, '*****','fontSize',11,'HorizontalAlignment','center');

h.YLim = ylimits;
h.XTick = xtick;
h.XTickLabel = [];
h.YTick = [-1 0 1];
h.TickLength = [0.01 0.01];
h.XAxis.TickLength = [1 1].*0.03;
ylabel('Corr (data, model)','Units','normalized','Position',[-0.13 0.5]);
h.XRuler.TickLabelGapMultiplier = -0.3;
h.YRuler.TickLabelGapMultiplier = 0.001;

% legend
offset_x = 0;
offset_y = -1.4;
patch([1 1 1.5 1.5]+offset_x, [0 .1 .1 0]+1.0+offset_y, 0.0*[1 1 1], 'Clipping','off');
patch([1 1 1.5 1.5]+offset_x, [0 .1 .1 0]+0.7+offset_y, 0.6*[1 1 1], 'Clipping','off');
text(1.8+offset_x,1.05+offset_y,'p_1^2', 'HorizontalAlignment','left','FontSize',7);
text(1.8+offset_x,0.75+offset_y,'p_2', 'HorizontalAlignment','left','FontSize',7);

%% MEC/CA3 spectrum
lw = 1 ;
axes(panel_G(1)); cla; hold on;
colorbar(panel_G(1), 'off');
text(-0.25,1.2, 'G', 'Units','normalized','FontWeight','bold');
for i3 = 1:nrCA3
    plot(kPow,powOutXs(:,i3,iang),'Color',clrrCA3(i3,:),'LineWidth',lw) ; hold on ;
end
plot(kPow,powOutS,'Color',clrrCA3(end,:),'LineWidth',lw) ; hold on ;
plot(kMEC,1.25e-1,'^r','MarkerFaceColor','r','MarkerSize',4) ; hold on ;
hax=gca;
hax.XLim = [0 0.25];
hax.YTick = 10.^[-1 0 1 2];
hax.YScale = 'log';
hax.TickLength(1)=0.02;
hax.XRuler.TickLabelGapOffset = -1;
xlabel('Spatial frequency (1/m)', 'Units','normalized','Position',[0.5 -0.15]);
ylabel('Power (norm.)', 'Units','normalized','Position',[-0.14 0.5]);
title('Model: Predicted spectrum of CA1 neurons','FontWeight','bold','Units','normalized','Position',[0.5 1.05]) ;
% zoom box
zoom_x = [0.035 0.08];
zoom_y = [5e-1 0.5e1];
plot(zoom_x,zoom_y([1 1]),'k-')
plot(zoom_x,zoom_y([2 2]),'k-')
plot(zoom_x([1 1]),zoom_y,'k-')
plot(zoom_x([2 2]),zoom_y,'k-')

% colorbar
colorbar_loc = hax.Position([1 2]) + hax.Position([3 4]).*[1 0.25];
colorbar_size = [0.2 0.5*hax.Position(4)];
colorbar_pos = [colorbar_loc colorbar_size];
hc=colorbar(panel_G(1),'manual');
hc.Units = 'centimeters';
hc.Position = colorbar_pos;
hc.Direction='reverse';
% hc.Ticks = [0 1];
% hc.TickLabels = {'MEC';'CA3'};
hc.Ticks = [];
text(1.02,0.808, 'MEC', 'Units','normalized','FontSize',6,'HorizontalAlignment','center');
text(1.02,0.2, 'CA3', 'Units','normalized','FontSize',6,'HorizontalAlignment','center');
colormap(panel_G(1),clrrCA3);

% zoom inset panel
axes(panel_G(2)); cla; hold on;
for i3 = 1:nrCA3
    plot(kPow,powOutXs(:,i3,iang),'Color',clrrCA3(i3,:),'LineWidth',lw) ; hold on ;
end
plot(kPow,powOutS,'Color',clrrCA3(end,:),'LineWidth',lw) ; hold on ;
hax=gca;
hax.XLim = zoom_x;
hax.YLim = zoom_y;
hax.YScale = 'log';
hax.XTick=[];
hax.YTick=[];
hax.Box='on';


%% Relative MEC/CA3 input sparsity - Spectrum peakedness
axes(panel_H(1)); cla; hold on;
text(-0.25,1.15, 'H', 'Units','normalized','FontWeight','bold');
plot(1-[rCA3 1],rMaxdpowOutXs(:,iang),'o-k','MarkerFaceColor','k','MarkerSize',4) ; hold on ;
xlabel('Relative CA3 vs. MEC input');
ylabel({'Spectral';'peak (norm.)'});
text(0,-0.3, 'CA3', 'Units','normalized','FontSize',7,'HorizontalAlignment','center');
text(1,-0.3, 'MEC', 'Units','normalized','FontSize',7,'HorizontalAlignment','center');
hax=gca;
hax.XRuler.TickLabelGapOffset = -1;
hax.TickLength(1) = 0.025;

%% Spectrum of real data maps (pooled over bats)
axes(panel_I(1)); cla; hold on;
text(-0.25,1.14, 'I', 'Units','normalized','FontWeight','bold');
shadedErrorBar(data.freq, data.maps_spec, {@mean,@nansem});
hax=gca;
hax.YScale = 'log';
hax.YMinorTick = 'on';
hax.YTick = 10.^[0 1 2];
hax.XLim = [0 0.25];
hax.TickLength(1) = 0.02;
xlabel('Spatial frequency (1/m)') ;
ylabel('Power (norm.)')
title('Data: Measured spectrum of CA1 neurons','FontWeight','bold','Units','normalized','Position',[0.5 1.05]) ;

% inset - binarized maps
axes(panel_I(2)); cla; hold on;
shadedErrorBar(data.freq, data.maps01_spec, {@mean,@nansem});
hax=gca;
hax.YScale = 'log';
hax.YMinorTick = 'on';
hax.YTick = 10.^[0 1 2];
hax.XLim = [0 0.25];
hax.TickLength(1) = 0.03;
hax.XRuler.TickLabelGapOffset = -1;
% hax.XRuler.TickLabelGapMultiplier = -0.35;
hax.YRuler.TickLabelGapOffset = -0.5;
% hax.YRuler.TickLabelGapMultiplier = 0.1;



%% ===================== Attractor netwoek model =========================
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% load data
load('L:\Misha_attractor\20200902__new_simulations\sim.mat');
load('L:\Misha_attractor\20200902__new_simulations\sim_res.mat');

%% model schematics
axes(panel_A);
cla('reset')
hold on
text(-0.1,1.05, 'A', 'Units','normalized','FontWeight','bold');
clr = [ 1 0 0;
        0 1 0;
        0 0 1;
        1 0 1;
        0 1 1;
        0.5 0.5 1;
        1 0.5 0.5;
        1 1 0;];
seg = [ 0   80;
        80  160;
        160 240;
        240 320;
        320 400;
        0   200;
        200 400;
        0   400;
        ];
bin_size = 0.5;
seg = seg .* bin_size;
dlm = 2.5;
seg_margins = [
    0 -dlm;
    dlm -dlm;
    dlm -dlm;
    dlm -dlm;
    dlm 0;
    0 -dlm;
    dlm 0;
    0 0;
];
levels = [3 3 3 3 3 2 2 1];
for ii_seg=1:size(seg,1)
    x = seg(ii_seg,:)+seg_margins(ii_seg,:);
    y = [levels([ii_seg ii_seg])];
    plot(x,y,'k','linewidth',2,'Color','k');
    text(x(1),y(1),"A"+ii_seg, 'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',8);
%     sigma = 0.05;
%     x = linspace(seg(ii,1),seg(ii,2),1000);
%     y = 0.65.*gaussmf(x,[sigma*range(seg(ii,:)) mean(seg(ii,:))])+levels(ii);
%     plot(x,y,'--','linewidth',1,'Color',clr(ii,:))
end
% cell locations on the attractors
smb = '^*dov<>phx+s';
% I flipped the order so cell 2 (green) will be on top of cell 4 (magenta)
for ii_cell = flip(1:length(cell_examples_IX))
    cell_num = cell_examples_IX(ii_cell);
    [r,c]=find(ind==cell_num);
    for ii_net = 1:length(c)
        net = c(ii_net);
        bn = th(r(ii_net),net);
        loc = bn * bin_size;
        % graphical manual jitter to avoid occlusion between cells 2 and 4:
        if ii_cell==2 && net==8
            loc = loc+1;
        end
        h=plot(loc,levels(net),smb(ii_cell),'Color',clr(ii_cell,:),'MarkerSize',4);
        h.MarkerFaceColor = clr(ii_cell,:);
    end
end
% connectivity of a single cell
cell_num = cell_examples_IX(1);
[r,c]=find(ind==cell_num);
for ii_net = 1:length(r)
    net = c(ii_net);
    bn = th(r(ii_net),net);
    loc = bn * bin_size;
    rad = 0.05 * range(seg(net,:));
    %
    x = linspace(loc,loc+rad,100);
    t = linspace(0,pi,length(x));
    y = sin(t);
    y = y.*0.25;
    y = y + levels(net);
    plot(x,y,'-','Color','r','LineWidth',1.5,'Clipping','off');
    %
    x = linspace(loc,loc-rad,100);
    t = linspace(0,pi,length(x));
    y = sin(t);
    y = y.*0.25;
    y = y + levels(net);
    plot(x,y,'-','Color','r','LineWidth',1.5,'Clipping','off');
    %
    x = linspace(loc,loc+2*rad,100);
    t = linspace(0,pi,length(x));
%     t = linspace(0,sqrt(pi),length(x));
%     t = flip(t);
%     t = t.^2;
    y = sin(t);
    y = y.*0.6;
    y = y + levels(net);
    plot(x,y,'-','Color','r','LineWidth',0.5,'Clipping','off');
    %
    x = linspace(loc,loc-2*rad,100);
    t = linspace(0,pi,length(x));
%     t = linspace(0,sqrt(pi),length(x));
%     t = flip(t);
%     t = t.^2;
    y = sin(t);
    y = y.*0.6;
    y = y + levels(net);
    plot(x,y,'-','Color','r','LineWidth',0.5,'Clipping','off');
end
hax=gca;
hax.YLim = [0 3.5];
hax.Visible = 'off';


%% cell examples 
locs = linspace(0,200,size(m,1));
run = size(m,3);
if multiple_panels

    for ii_cell = 1:length(panel_B)
        axes(panel_B(ii_cell));
        cla('reset')
        hold on
        
        plot(locs, m(:,cell_examples_IX(ii_cell),run),'LineWidth',1.5,'Color', clr(ii_cell,:));
        text(0.01,0.85,"Cell "+ii_cell,'Units','normalized','FontSize',8);
%         text(0.01,0.85,"Cell "+cell_examples_IX(ii_cell),'Units','normalized','FontSize',8);
        hax=gca;
        hax.XTick=[];
        hax.YTick=[];
    end
    axes(panel_B(1));
    text(-0.1,1.05, 'B', 'Units','normalized','FontWeight','bold');
    axes(panel_B(end));
    xlabel('Position (m)');
    xlim([0 200]);
    xticks([0:50:200]);
    axes(panel_B(3));
    ylabel('Firing rate (a.u.)');


else % single panel

    axes(panel_B);
    cla('reset')
    hold on
    text(-0.1,1.05, 'H', 'Units','normalized','FontWeight','bold');
    hl=plot(locs,m(:,cell_examples_IX,run),'LineWidth',2);
    for ii_hl = 1:length(hl)
        hl(ii_hl).Color = clr(ii_hl,:);
    end
    xlabel('Position (m)');
    xlim([0 200]);
    xticks([0:50:200]);
    ylabel('Firing rate (a.u.)');
    hax=gca;
    hax.XTick=[];
    hax.YTick=[];
    
end

%%
CANN_hist_clr = [0.55 0.5 0.45];

%% field num
axes(panel_C(1));
cla('reset')
hold on
text(-0.45,1.2, 'C', 'Units','normalized','FontWeight','bold');
h=histogram(res.fnumber);
h.Normalization = 'probability';
h.BinEdges = 0.5+[0:5];
h.FaceColor = CANN_hist_clr;
xlabel({'No. of fields per direction'},'Units','normalized','Position',[0.5 -0.18]);
ylabel('Fraction of cells','Units','normalized','Position',[-0.3 0.5])
hax=gca;
hax.YScale = 'log';
hax.YLim = [1e-3 5e-1];
% hax.XTick=[];
hax.YTick=10.^[-3:0];
hax.TickDir='out';
hax.TickLength = [0.03 0.03];
hax.XRuler.TickLabelGapOffset = -.5;

%% field size
axes(panel_C(2));
cla('reset')
hold on

h=histogram(res.fsizes_all.*bin_size);
h.Normalization = 'probability';
h.FaceColor = CANN_hist_clr;
h.BinWidth = 1;
% text(0.5,1,sprintf('bin=%.1fm',h.BinWidth),'Units','normalized','FontSize',8);

xlabel('Field size (m)','Units','normalized','Position',[0.5 -0.18]);
ylabel('Fraction of fields','Units','normalized','Position',[-0.3 0.5]);
hax=gca;
hax.YScale = 'log';
hax.XTick = [0 20 40];
hax.YTick = 10.^[-3 -2 -1 0];
hax.XLim = [0 40];
hax.YLim = [1e-3 2e-1];
hax.TickLength = [0.03 0.03];
hax.XRuler.TickLabelGapOffset = -.5;

%% field ratio L/S
axes(panel_C(3));
cla('reset')
hold on
h=histogram(res.sratio);
nbins = 40;
h.BinEdges = logspace(0,2,nbins);
h.Normalization = 'pdf';
h.FaceColor = CANN_hist_clr;
% text(0.5,1,sprintf('nbins=%d',nbins),'Units','normalized','FontSize',8);

xlabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[0.5 -0.17]);
ylabel('Fraction of cells','Units','normalized','Position',[-0.3 0.5]);
hax=gca;
hax.XScale = 'log';
hax.YScale = 'log';
% hax.XTick=[];
% hax.YTick=[];
hax.XLim = [0 20];
hax.XLim(2) = h.BinEdges(find(h.BinEdges>hax.XLim(2),1,'first'));
hax.YLim = [1e-3 3e0];
% hax.XTick = [1 2 5 10 20];
hax.YTick = 10.^[-3 -2 -1 0];
hax.TickLength = [0.03 0.03];
hax.XRuler.TickLabelGapOffset = -.5;

%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
% fig_name_out = fig_name_out + "_cells" + join("_"+cell_examples_IX,'');
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
disp('figure was successfully saved to pdf/tiff/fig formats');





%%

