%% Large Scale - Fig. 4 - Theoretical analysis

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'Fig_4';
fig_caption_str = 'Theoretical analysis';
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
panel_A_size = [3 3];
panel_A(1) = axes('position', [ 1.5 21 panel_A_size]);
panel_A(2) = axes('position', [ 5.5 21 panel_A_size]);
panel_A(3) = axes('position', [ 9.5 21 panel_A_size]);
panel_A(4) = axes('position', [13.5 21 panel_A_size]);
panel_A(5) = axes('position', [17.5 21 panel_A_size]);
panel_B    = axes('position', [ 1.5 15 5 4]);
panel_C(1,1) = axes('position', [ 8.5 15 5 4]);
panel_C(1,2) = axes('position', [ 9   17 3 2]);
panel_C(2,1) = axes('position', [15   15 5 4]);
panel_C(2,2) = axes('position', [15.5   17 3 2]);
panel_D(1)   = axes('position', [ 8.5 9.5 5 4]);
panel_D(2)   = axes('position', [15   9.5 5 4]);
panel_E(1)   = axes('position', [ 8.5 4   5 4]);
panel_E(2)   = axes('position', [15   4   5 4]);
panel_legend = axes('position', [2 11.5 0.5 2]);

%% arrange data
dt = 0.5 ;
paper_fig_4_arrange_data;
rng(0);
exampleL = 200*100; % 200*cm
[f] = MultiscalePlace_GenerateTuning_AllModels(exampleL);
pos = 1:ds:exampleL;
pos = pos / 100; % back to meter

%% legend panel
axes(panel_legend);
cla
hold on
plot([0 1], 5*[1 1], 'Color',clr(1,:),'LineWidth',2) ; hold on ; 
plot([0 1], 4*[1 1], 'Color',clr(2,:),'LineWidth',2) ; hold on ; 
plot([0 1], 3*[1 1], 'Color',clr(3,:),'LineWidth',2) ; hold on ; 
plot([0 1], 2*[1 1], 'Color',clr(6,:),'LineWidth',2) ; hold on ; 
plot([0 1], 1*[1 1], 'Color',clr(7,:),'LineWidth',2) ; hold on ; 
text(1.2,5,'1: Single small fields', 'FontSize', 10);
text(1.2,4,'2: Single large fields', 'FontSize', 10);
text(1.2,3,'3: multiple small fields', 'FontSize', 10);
text(1.2,2,'4: multi-scale (population)', 'FontSize', 10);
text(1.2,1,'5: multi-scale (single-cell)', 'FontSize', 10);
axis off
xlim([0 1])
ylim([1 5])

%% panel A - different encoding schemes
sdf=[];
sdf(1,:,:) = f.fA(:,1:10);
sdf(2,:,:) = f.fB(:,1:10);
sdf(3,:,:) = f.fC(:,1:10);
sdf(4,:,:) = f.fF(:,950:-100:50);
sdf(5,:,:) = f.fG(:,1:10);
for ii_scnr = 1:5
    axes(panel_A(ii_scnr));
    cla
    hold on
    imagesc(pos,1:10,squeeze(sdf(ii_scnr,:,:))');
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
text(-0.2,1.1, 'A', 'Units','normalized','FontWeight','bold');
ylabel('Example neuron #','Units','normalized','Position',[-0.1 0.5]);

%% panel B - minimum N required for error < 1 m
axes(panel_B);
cla
hold on
text(-0.2,1.1, 'B', 'Units','normalized','FontWeight','bold');

jerr = 1;
h(1) = plot(L/100,NerrMLA(:,jerr),'Color',clr(1,:),'LineWidth',2) ; hold on ; 
h(2) = plot(L/100,NerrMLB(:,jerr),'Color',clr(2,:),'LineWidth',2) ; hold on ; 
h(3) = plot(L/100,NerrMLC(:,jerr),'Color',clr(3,:),'LineWidth',2) ; hold on ; 
h(4) = plot(L/100,NerrMLH(:,jerr),'Color',clr(8,:),'LineWidth',2) ; hold on ; 
h(5) = plot(L/100,NerrMLD(:,jerr),'Color',clr(4,:),'LineWidth',2) ; hold on ; 
h(6) = plot(L/100,NerrMLF(:,jerr),'Color',clr(6,:),'LineWidth',2) ; hold on ; 
h(7) = plot(L/100,NerrMLE(:,jerr),'Color',clr(5,:),'LineWidth',2) ; hold on ;
h(8) = plot(L/100,NerrMLG(:,jerr),'Color',clr(7,:),'LineWidth',2) ; hold on ;
xlim([20 1000]) ;
ylim([10 250]) ;
xlabel('Environment size (m)', 'Units','normalized','Position',[0.5 -0.15]);
ylabel(['min N required for error < ' num2str(errs(jerr)) '(m)'], 'Units','normalized','Position',[ -0.15 0.5]);
h = gca;
h.XRuler.TickLabelGapOffset = -1;
h.XRuler.TickLength = [0.03 0.03];
h = gca;
% h.XScale = 'log';
% h.YScale = 'log';

%% panel C - mean decoding error
jN_options = find(ismember(N, [50 100]));
ylimits_options = [0 20; 0 10];
for ii_N = 1:length(jN_options)
    jN = jN_options(ii_N);
    
    axes(panel_C(ii_N,1));
    cla
    hold on
    plot(L/100,meMLA(jN,:)*ds/100,'Color',clr(1,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLB(jN,:)*ds/100,'Color',clr(2,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLC(jN,:)*ds/100,'Color',clr(3,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLH(jN,:)*ds/100,'Color',clr(8,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLD(jN,:)*ds/100,'Color',clr(4,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLF(jN,:)*ds/100,'Color',clr(6,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLE(jN,:)*ds/100,'Color',clr(5,:),'LineWidth',2) ; hold on ;
    plot(L/100,meMLG(jN,:)*ds/100,'Color',clr(7,:),'LineWidth',2) ; hold on ; 
    xlim([20 1000]) ;
    xlabel('Environment size (m)', 'Units','normalized','Position',[0.5 -0.11]);
    ylabel('Mean decoding error (m)', 'Units','normalized','Position',[ -0.15 0.5]);
    h = gca;
    h.XScale = 'log';
%     h.YScale = 'log';
    h.XRuler.TickLabelGapOffset = -1;
    h.XRuler.TickLength = [0.03 0.03];
    title(sprintf('N = %d neurons',N(jN)))
    if ii_N == 1
        text(-0.2,1.1, 'C', 'Units','normalized','FontWeight','bold');
    end

    axes(panel_C(ii_N,2));
    cla
    hold on
    plot(L/100,meMLA(jN,:)*ds/100,'Color',clr(1,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLB(jN,:)*ds/100,'Color',clr(2,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLC(jN,:)*ds/100,'Color',clr(3,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLH(jN,:)*ds/100,'Color',clr(8,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLD(jN,:)*ds/100,'Color',clr(4,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLF(jN,:)*ds/100,'Color',clr(6,:),'LineWidth',2) ; hold on ; 
    plot(L/100,meMLE(jN,:)*ds/100,'Color',clr(5,:),'LineWidth',2) ; hold on ;
    plot(L/100,meMLG(jN,:)*ds/100,'Color',clr(7,:),'LineWidth',2) ; hold on ; 
    xlim([20 1000]) ;
    ylim(ylimits_options(ii_N,:))
    h = gca;
    h.XScale = 'log';
%     h.YScale = 'log';
    h.XRuler.TickLabelGapOffset = -1;
    h.XRuler.TickLength = [0.03 0.03];
end

%% panel D - 99th prc error
jprc = find(ismember(prc, [99]));
jN_options = find(ismember(N, [50 100]));
ylimits_options = [0 20; 0 10];
for ii_N = 1:length(jN_options)
    jN = jN_options(ii_N);
    
    axes(panel_D(ii_N));
    cla
    hold on
    plot(L/100,peMLA(jN,:,jprc)*ds/100,'Color',clr(1,:),'LineWidth',2) ; hold on ; 
    plot(L/100,peMLB(jN,:,jprc)*ds/100,'Color',clr(2,:),'LineWidth',2) ; hold on ; 
    plot(L/100,peMLC(jN,:,jprc)*ds/100,'Color',clr(3,:),'LineWidth',2) ; hold on ; 
    plot(L/100,peMLH(jN,:,jprc)*ds/100,'Color',clr(8,:),'LineWidth',2) ; hold on ; 
    plot(L/100,peMLD(jN,:,jprc)*ds/100,'Color',clr(4,:),'LineWidth',2) ; hold on ; 
    plot(L/100,peMLF(jN,:,jprc)*ds/100,'Color',clr(6,:),'LineWidth',2) ; hold on ; 
    plot(L/100,peMLE(jN,:,jprc)*ds/100,'Color',clr(5,:),'LineWidth',2) ; hold on ; 
    plot(L/100,peMLG(jN,:,jprc)*ds/100,'Color',clr(7,:),'LineWidth',2) ; hold on ; 
%     ylim([0.1 2000]) ;
    xlim([20 1000]) ;
    xlabel('Environment size (m)', 'Units','normalized','Position',[0.5 -0.15]);
    ylabel([num2str(prc(jprc)) '% decoder error (m)'], 'Units','normalized','Position',[ -0.15 0.5]);
    title(sprintf('N = %d neurons',N(jN)))
    h = gca;
    h.XScale = 'log';
%     h.YScale = 'log';
    h.XRuler.TickLabelGapOffset = -1;
    h.XRuler.TickLength = [0.03 0.03];
    if ii_N == 1
        text(-0.2,1.1, 'D', 'Units','normalized','FontWeight','bold');
    end
end

%% panel E - 99th prc error
jprc = find(ismember(prc, [99]));
jN_options = find(ismember(N, [50 100]));
ylimits_options = [2e-4 1; 2e-7 1];
for ii_N = 1:length(jN_options)
    jN = jN_options(ii_N);
    
    axes(panel_E(ii_N));
    cla
    hold on

    plot(L/100,pleMLA(jN,:)+(2e6)^(-1),'Color',clr(1,:),'LineWidth',2) ; hold on ; 
    plot(L/100,pleMLB(jN,:)+(2e6)^(-1),'Color',clr(2,:),'LineWidth',2) ; hold on ; 
    plot(L/100,pleMLC(jN,:)+(2e6)^(-1),'Color',clr(3,:),'LineWidth',2) ; hold on ; 
    plot(L/100,pleMLH(jN,:)+(2e6)^(-1),'Color',clr(8,:),'LineWidth',2) ; hold on ; 
    plot(L/100,pleMLD(jN,:)+(2e6)^(-1),'Color',clr(4,:),'LineWidth',2) ; hold on ; 
    plot(L/100,pleMLF(jN,:)+(2e6)^(-1),'Color',clr(6,:),'LineWidth',2) ; hold on ; 
    plot(L/100,pleMLE(jN,:)+(2e6)^(-1),'Color',clr(5,:),'LineWidth',2) ; hold on ; 
    plot(L/100,pleMLG(jN,:)+(2e6)^(-1),'Color',clr(7,:),'LineWidth',2) ; hold on ; 
    xlim([20 1000]) ;
    ylim(ylimits_options(ii_N,:))
    xlabel('Environment size (m)', 'Units','normalized','Position',[0.5 -0.15]);
    ylabel('P(err>5% of environment size)', 'Units','normalized','Position',[ -0.15 0.5]);
    title(sprintf('N = %d neurons',N(jN)))
    h = gca;
    h.XScale = 'log';
    h.YScale = 'log';
    h.XRuler.TickLabelGapOffset = -1;
    h.XRuler.TickLength = [0.03 0.03];
    if ii_N == 1
        text(-0.2,1.1, 'E', 'Units','normalized','FontWeight','bold');
    end
end


%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');

