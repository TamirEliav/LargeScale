%% Large Scale - Fig. S10 - decoder with shorter dt (200ms) + error vs. dt

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'Fig_S10_decoder_changing_dt';
fig_caption_str = 'Theoretical analysis - error vs. integration time (dt)';
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
panel_A_size = [3 3].*[0.95 0.85];
panel_A(1) = axes('position', [ 1.5 21 panel_A_size]);
panel_A(2) = axes('position', [ 5.0 21 panel_A_size]);
panel_A(3) = axes('position', [ 8.5 21 panel_A_size]);
panel_A(4) = axes('position', [12.0 21 panel_A_size]);
panel_A(5) = axes('position', [15.5 21 panel_A_size]);
panel_B(1)   = axes('position', [ 1.5 15 5 4]);
panel_B(2)   = axes('position', [ 7.5 15 1.5 4]);
panel_C(1,1) = axes('position', [ 1.5  9.5  5 4]);
panel_C(1,2) = axes('position', [ 2   11.5  3 2]);
panel_D(1,1) = axes('position', [ 7.7    9.5  5 4]);
panel_D(1,2) = axes('position', [ 8.2 11.5  3 2]);
panel_E(1)   = axes('position', [14    9.5  5 4]);
panel_F      = axes('position', [ 7.4 3.3 6 4]);
panel_legend = axes('position', [9.5 17 0.5 2]);

%% arrange data
dt = 0.2;
coverage = 0.20; % 20%
% coverage = 0.15; % 15%
paper_fig_4_arrange_data;
rng(0);
exampleL = 200*100; % 200*cm
[f] = MultiscalePlace_GenerateTuning_AllModels(exampleL,coverage);
pos = 1:ds:exampleL;
pos = pos / 100; % back to meter
panels_AE_dt = dt;

%% legend panel
axes(panel_legend);
cla
hold on
plot([0 1], 5*[1 1], 'Color',clr(1,:),'LineWidth',2) ; hold on ; 
plot([0 1], 4*[1 1], 'Color',clr(3,:),'LineWidth',2) ; hold on ; 
plot([0 1], 3*[1 1], 'Color',clr(2,:),'LineWidth',2) ; hold on ; 
plot([0 1], 2*[1 1], 'Color',clr(6,:),'LineWidth',2) ; hold on ; 
plot([0 1], 1*[1 1], 'Color',clr(7,:),'LineWidth',2) ; hold on ; 
text(1.2,5,'1: Single small fields', 'FontSize', 10);
text(1.2,4,'2: Single large fields', 'FontSize', 10);
text(1.2,3,'3: multiple small fields (Rich et al. 2014)', 'FontSize', 10);
text(1.2,2,'4: multi-scale (population)', 'FontSize', 10);
text(1.2,1,'5: multi-scale (single-cell)', 'FontSize', 10);
axis off
xlim([0 1])
ylim([1 5])

%% panel A - different encoding schemes
maps=[];
maps(1,:,:) = f.fA(:,1:10);
maps(1,:,[7 4 9 10]) = f.fA(:,[4 7 11 12]);
maps(2,:,:) = f.fB(:,1:10);
maps(2,:,[3 9]) = f.fB(:,[14 12]);
[~,IX] = sort(sum(f.fI,1));
maps(3,:,:) = f.fI(:,IX(50:100:950));
maps(4,:,:) = f.fF(:,950:-100:50);
maps(5,:,:) = f.fG(:,1:10);
for ii_scnr = 1:5
    axes(panel_A(ii_scnr));
    cla
    hold on
    imagesc(pos,1:10,squeeze(maps(ii_scnr,:,:))');
    cmap = gray();
    cmap = flip(cmap);
    colormap(cmap)
    xlim([0 200])
    ylim([0.5 10.5])
    xlabel('Position (m)', 'Units','normalized','Position',[0.5 -0.15]);
    h = gca;
    h.YTick = [1 5 10];
    h.XRuler.TickLabelGapOffset = -1;
    h.XRuler.TickLength = [0.03 0.03];
    title("Scheme "+ii_scnr);
end
axes(panel_A(1));
text(-0.3,1.1, 'A', 'Units','normalized','FontWeight','bold');
ylabel('Example neuron #','Units','normalized','Position',[-0.2 0.5]);

%% panel B - minimum N required for error < 1 m
axes(panel_B(1));
cla
hold on
text(-0.2,1.1, 'B', 'Units','normalized','FontWeight','bold');
text(0.5,1.1, {'Minimal no. of neurons';'required for decoding'}, 'Units','normalized','FontWeight','bold','HorizontalAlignment','center','FontSize',9);

jerr = find(ismember(errs, [2]));
plot(L/100,NerrMLA(:,jerr),'Color',clr(1,:),'LineWidth',2) ; hold on ; 
plot(L/100,NerrMLB(:,jerr),'Color',clr(3,:),'LineWidth',2) ; hold on ; 
plot(L/100,NerrMLI(:,jerr),'Color',clr(2,:),'LineWidth',2) ; hold on ;
plot(L/100,NerrMLF(:,jerr),'Color',clr(6,:),'LineWidth',2) ; hold on ; 
plot(L/100,NerrMLG(:,jerr),'Color',clr(7,:),'LineWidth',2) ; hold on ;
xlim([20 1000]) ;
ylim([10 250]) ;
xlabel('Environment size (m)', 'Units','normalized','Position',[0.5 -0.11]);
ylabel(['min N required for error < ' num2str(errs(jerr)) 'm'], 'Units','normalized','Position',[ -0.15 0.5]);
h = gca;
h.XRuler.TickLabelGapOffset = -1;
h.XRuler.TickLength = [0.02 0.02];
h.XTick = [20 200 400 600 800 1000];
% h.XScale = 'log';
% h.YScale = 'log';

%% panel B - bar plot of slopes
x = 1:5;
y = [
    lmNerrMLA{jerr}.Coefficients.Estimate(2)
    lmNerrMLB{jerr}.Coefficients.Estimate(2)
    lmNerrMLI{jerr}.Coefficients.Estimate(2)
    lmNerrMLF{jerr}.Coefficients.Estimate(2)
    lmNerrMLG{jerr}.Coefficients.Estimate(2)]';
err = [
    lmNerrMLA{jerr}.coefCI
    lmNerrMLB{jerr}.coefCI
    lmNerrMLI{jerr}.coefCI
    lmNerrMLF{jerr}.coefCI
    lmNerrMLG{jerr}.coefCI];
err  = err(2:2:end,:)';
err = err - y;

axes(panel_B(2));
cla
hold on
hb=bar(x,y);
hb.FaceColor = 'flat';
hb.CData = [clr(1,:);clr(3,:);clr(2,:);clr(6,:);clr(7,:)];
he=errorbar(x,y,err(1,:));
he.CapSize = 2;
he.LineStyle = 'none';
he.Color = 'k';
m = 1.1 * max(y+err(1,:));
m = round(m,1);
xlim([0 6]);
ylim([0 m]);
h=gca;
h.XTick = [];
h.YTick = [0 m];
ylabel('Slope (N per meter)', 'Units','normalized','Position',[-0.11 0.5]);

%% panel C - mean decoding error
jN_options = find(ismember(N, [50]));
ylimits_options = [0 20; 0 10];
for ii_N = 1:length(jN_options)
    jN = jN_options(ii_N);
    
    axes(panel_C(ii_N,1));
    cla
    hold on
    plot(L/100,meMLA(jN,:)*ds/100,'Color',clr(1,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLB(jN,:)*ds/100,'Color',clr(3,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLI(jN,:)*ds/100,'Color',clr(2,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLF(jN,:)*ds/100,'Color',clr(6,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLG(jN,:)*ds/100,'Color',clr(7,:),'LineWidth',2) ; hold on ; 
    xlim([20 1000]) ;
    xlabel('Environment size (m)', 'Units','normalized','Position',[0.5 -0.11]);
    ylabel('Mean decoding error (m)', 'Units','normalized','Position',[ -0.15 0.5]);
    h = gca;
    h.XScale = 'log';
%     h.YScale = 'log';
    h.XTick = [20 100 1000];
    h.XTickLabel = {'20';'100';'1000'};
    h.XRuler.TickLabelGapOffset = -1;
    h.XRuler.TickLength = [0.03 0.03];
%     text(0.95,1, sprintf('N = %d',N(jN)), 'Units','normalized','FontWeight','normal', 'HorizontalAlignment','right','FontSize',10);
    if ii_N == 1
        text(-0.2,1.1, 'C', 'Units','normalized','FontWeight','bold');
        text(0.5,1.1, 'Mean decoding errors', 'Units','normalized','FontWeight','bold','HorizontalAlignment','center','FontSize',9);
    end

    % zoom-in inset
    axes(panel_C(ii_N,2));
    cla
    hold on
    plot(L/100,meMLA(jN,:)*ds/100,'Color',clr(1,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLB(jN,:)*ds/100,'Color',clr(3,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLI(jN,:)*ds/100,'Color',clr(2,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLF(jN,:)*ds/100,'Color',clr(6,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLG(jN,:)*ds/100,'Color',clr(7,:),'LineWidth',2) ; hold on ; 
    xlim([20 1000]) ;
    ylim(ylimits_options(ii_N,:))
    h = gca;
    h.XScale = 'log';
%     h.YScale = 'log';
    h.XTick = [20 100 1000];
    h.XTickLabel = {'20';'100';'1000'};
    h.XRuler.TickLabelGapOffset = -1;
    h.XRuler.TickLength = [0.03 0.03];
end

%% panel D - 99th prc error
jprc = find(ismember(prc, [99]));
jN_options = find(ismember(N, [50]));
ylimits_options = [0 20; 0 10];
for ii_N = 1:length(jN_options)
    jN = jN_options(ii_N);
    
    axes(panel_D(ii_N,1));
    cla
    hold on
    plot(L/100,peMLA(jN,:,jprc)*ds/100,'Color',clr(1,:),'LineWidth',2) ; hold on ; 
    plot(L/100,peMLB(jN,:,jprc)*ds/100,'Color',clr(3,:),'LineWidth',2) ; hold on ; 
    plot(L/100,peMLI(jN,:,jprc)*ds/100,'Color',clr(2,:),'LineWidth',2) ; hold on ; 
    plot(L/100,peMLF(jN,:,jprc)*ds/100,'Color',clr(6,:),'LineWidth',2) ; hold on ; 
    plot(L/100,peMLG(jN,:,jprc)*ds/100,'Color',clr(7,:),'LineWidth',2) ; hold on ; 
%     ylim([0.1 2000]) ;
    xlim([20 1000]) ;
    xlabel('Environment size (m)', 'Units','normalized','Position',[0.5 -0.11]);
    ylabel([num2str(prc(jprc)) '% decoder error (m)'], 'Units','normalized','Position',[ -0.12 0.5]);
    h = gca;
    h.XScale = 'log';
%     h.YScale = 'log';
    h.XTick = [20 100 1000];
    h.XTickLabel = {'20';'100';'1000'};
    h.YRuler.TickLabelGapOffset = 2;
    h.XRuler.TickLabelGapOffset = -1;
    h.XRuler.TickLength = [0.03 0.03];
%     text(0.95,1, sprintf('N = %d',N(jN)), 'Units','normalized','FontWeight','normal', 'HorizontalAlignment','right','FontSize',10);
    if ii_N == 1
        text(-0.19,1.1, 'D', 'Units','normalized','FontWeight','bold');
        text(1.12,1.22, 'Catastrophic errors', 'Units','normalized','FontWeight','bold','HorizontalAlignment','center','FontSize',9);
    end
    
    % zoom-in inset
    axes(panel_D(ii_N,2));
    cla
    hold on
    plot(L/100,peMLA(jN,:,jprc)*ds/100,'Color',clr(1,:),'LineWidth',2) ; hold on ; 
    plot(L/100,peMLB(jN,:,jprc)*ds/100,'Color',clr(3,:),'LineWidth',2) ; hold on ; 
    plot(L/100,peMLI(jN,:,jprc)*ds/100,'Color',clr(2,:),'LineWidth',2) ; hold on ; 
    plot(L/100,peMLF(jN,:,jprc)*ds/100,'Color',clr(6,:),'LineWidth',2) ; hold on ; 
    plot(L/100,peMLG(jN,:,jprc)*ds/100,'Color',clr(7,:),'LineWidth',2) ; hold on ; 
%     ylim([0.1 2000]) ;
    xlim([100 1000]) ;
    h = gca;
    h.XScale = 'log';
    h.YScale = 'log';
    h.XTick = [100 1000];
    h.XTickLabel = {'100';'1000'};
    h.YTick = [1 100];
    h.YTickLabel = {'1';'100'};
    h.YRuler.TickLabelGapOffset = -0.1;
    h.XRuler.TickLabelGapOffset = -1;
    h.XRuler.TickLength = [0.03 0.03];
end

%% panel E - 99th prc error
jprc = find(ismember(prc, [99]));
jN_options = find(ismember(N, [50]));
ylimits_options = [2e-4 1; 2e-7 1];
for ii_N = 1:length(jN_options)
    jN = jN_options(ii_N);
    
    axes(panel_E(ii_N));
    cla
    hold on

    plot(L/100,pleMLA(jN,:)+(2e6)^(-1),'Color',clr(1,:),'LineWidth',2) ; hold on ; 
    plot(L/100,pleMLB(jN,:)+(2e6)^(-1),'Color',clr(3,:),'LineWidth',2) ; hold on ; 
    plot(L/100,pleMLI(jN,:)+(2e6)^(-1),'Color',clr(2,:),'LineWidth',2) ; hold on ; 
    plot(L/100,pleMLF(jN,:)+(2e6)^(-1),'Color',clr(6,:),'LineWidth',2) ; hold on ; 
    plot(L/100,pleMLG(jN,:)+(2e6)^(-1),'Color',clr(7,:),'LineWidth',2) ; hold on ; 
    xlim([20 1000]) ;
    ylim(ylimits_options(ii_N,:))
    xlabel('Environment size (m)', 'Units','normalized','Position',[0.5 -0.11]);
    ylabel('P(err>5% of environment size)', 'Units','normalized','Position',[ -0.12 0.5]);
%     text(0.95,1.03, sprintf('N = %d',N(jN)), 'Units','normalized','FontWeight','normal', 'HorizontalAlignment','right','FontSize',10);
    h = gca;
    h.XScale = 'log';
    h.YScale = 'log';
    h.XTick = [20 100 1000];
    h.XTickLabel = {'20';'100';'1000'};
    h.YRuler.TickLabelGapOffset = -0.1;
    h.XRuler.TickLabelGapOffset = -1;
    h.XRuler.TickLength = [0.03 0.03];
    if ii_N == 1
        text(-0.19,1.1, 'E', 'Units','normalized','FontWeight','bold');
    end
end

%% panel F - decoding error vs. dt
axes(panel_F);
cla
hold on
text(-0.2,1.1, 'F', 'Units','normalized','FontWeight','bold');

% load data
switch coverage
    case 0.20
        load('L:\Theory\Yonatan_code_data\Multiscale_PV_ML_decoder_5_Summary_1.mat');
    case 0.15
        load('L:\Theory\Yonatan_code_data\Multiscale_PV_ML_decoder_7_Summary_1.mat');
end
nL = length(L);
ds = 20 ;

jN = find(ismember(N, [50]));
jL = find(ismember(L, [20000]));
plot(dt.*1e3, squeeze(meMLA(jN,:,jL)), 'Color', clr(1,:), 'LineWidth', 2);
plot(dt.*1e3, squeeze(meMLB(jN,:,jL)), 'Color', clr(3,:), 'LineWidth', 2);
plot(dt.*1e3, squeeze(meMLI(jN,:,jL)), 'Color', clr(2,:), 'LineWidth', 2);
plot(dt.*1e3, squeeze(meMLF(jN,:,jL)), 'Color', clr(6,:), 'LineWidth', 2);
plot(dt.*1e3, squeeze(meMLG(jN,:,jL)), 'Color', clr(7,:), 'LineWidth', 2);

ha= gca;
ha.XLim = [20 500];
% ha.YLim = [0.8 310];
ha.XTick = [20 100:100:500];
% ha.YTick = [1 10 100];
% ha.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
ha.TickDir='out';
ha.TickLength = [0.02 0.02];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;
xlabel('Integration window (ms)')
ylabel('Mean decoding error (m)','Units','normalized','Position',[-0.12 0.5]);
ha.XScale = 'linear';
ha.YScale = 'log';

%% print/save the figure
% fig_name_out = fullfile(res_dir, fig_name_str);
fig_name_out = fullfile(res_dir, [fig_name_str '_dt_' strrep(num2str(panels_AE_dt),'.','_') '_coverage=' num2str(100*coverage)]);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');

