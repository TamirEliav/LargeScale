%% Large Scale - Fig. S10 - PV decoder

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S10';
fig_caption_str = 'Theoretical analysis with PV decoder';
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
panel_BCDE_size = [4.6 4];
panel_A(1) = axes('position', [ 2.0 21 panel_A_size]);
panel_A(2) = axes('position', [ 5.5 21 panel_A_size]);
panel_A(3) = axes('position', [ 9.0 21 panel_A_size]);
panel_A(4) = axes('position', [12.5 21 panel_A_size]);
panel_A(5) = axes('position', [16.0 21 panel_A_size]);
panel_B(1)   = axes('position', [ 2.0 14.6 panel_BCDE_size]);
panel_B(2)   = axes('position', [ 7.6 14.6 2 4]);
panel_C(1,1) = axes('position', [ 2.0  8.7  panel_BCDE_size]);
panel_C(1,2) = axes('position', [ 2.5   10.7  2.8 2]);
panel_D(1,1) = axes('position', [ 8.1  8.7  panel_BCDE_size]);
panel_D(1,2) = axes('position', [ 8.7 10.7  2.5 2]);
panel_E(1)   = axes('position', [14.4    8.7  panel_BCDE_size]);
panel_legend = axes('position', [9.8 16.6 0.4 2]);

%% arrange data
dt = 0.5;
coverage = 0.20; % 20%
% coverage = 0.15; % 15%
paper_fig_4_arrange_data;
rng(0);
exampleL = 200*100; % 200*cm
[f] = MultiscalePlace_GenerateTuning_AllModels(exampleL,coverage);
pos = 1:ds:exampleL;
pos = pos / 100; % back to meter

%% legend panel
axes(panel_legend);
cla
hold on
plot([0 1], 5*[1 1], 'Color',clr(1,:),'LineWidth',2) ; hold on ; 
plot([0 1], 4*[1 1], 'Color',clr(3,:),'LineWidth',2) ; hold on ; 
plot([0 1], 3*[1 1], 'Color',clr(2,:),'LineWidth',2) ; hold on ; 
plot([0 1], 2*[1 1], 'Color',clr(6,:),'LineWidth',2) ; hold on ; 
plot([0 1], 1*[1 1], 'Color',clr(7,:),'LineWidth',2) ; hold on ; 
text(1.2,5,'1: Single small field', 'FontSize', 9);
text(1.2,4,'2: Single large field', 'FontSize', 9);
text(1.2,3,'3: Multiple small fields (Rich et al. 2014)', 'FontSize', 9);
text(1.2,2,'4: Multi-scale (population)', 'FontSize', 9);
text(1.2,1,'5: Multi-scale (single-cell)', 'FontSize', 9);
axis off
xlim([0 1])
ylim([1 5])

%% panel A - different encoding schemes
maps=[];
maps(1,:,:) = f.fA(:,1:10);
maps(1,:,[7 4 9 10]) = f.fA(:,[4 7 11 12]);
maps(2,:,:) = f.fB(:,1:10);
maps(2,:,[3 9 5]) = f.fB(:,[14 12 15]);
[~,IX] = sort(sum(f.fI,1));
maps(3,:,:) = f.fI(:,IX(50:100:950));
maps(3,:,10) = f.fI(:,IX(951));
rng(0);
maps(3,:,:) = maps(3,:,randperm(10));
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
    ht=title("Scheme "+ii_scnr,'Units','normalized','Position',[0.5 1.06]);
    schemes_IX = [1 3 2 6 7];
    schemes_title_color_opt = 2;
    switch schemes_title_color_opt
        case 1
            ht.Color = clr(schemes_IX(ii_scnr),:);
        case 2
            plot([0 25]+14, 11.8*[1 1], 'Color', clr(schemes_IX(ii_scnr),:), 'LineWidth', 2, 'Clipping','off');
    end
end
axes(panel_A(1));
text(-0.4,1.1, 'A', 'Units','normalized','FontWeight','bold');
text(3,1.5, 'Population Vector decoder instead of Maximum Likelihood decoder', 'Units','normalized','FontWeight','bold','FontSize',12,'HorizontalAlignment','center');
ylabel('Example neuron no.','Units','normalized','Position',[-0.2 0.5]);

%% panel B - minimum N required for error < 1 m
axes(panel_B(1));
cla
hold on
text(-0.24,1.13, 'B', 'Units','normalized','FontWeight','bold');
text(0.5,1.13, {'Minimal no. of neurons';'required for decoding'}, 'Units','normalized','FontWeight','bold','HorizontalAlignment','center','FontSize',9);

jerr = find(ismember(errs, [2]));
plot(L/100,NerrPVA(:,jerr),'Color',clr(1,:),'LineWidth',2) ; hold on ; 
plot(L/100,NerrPVB(:,jerr),'Color',clr(3,:),'LineWidth',2) ; hold on ; 
plot(L/100,NerrPVI(:,jerr),'Color',clr(2,:),'LineWidth',2) ; hold on ;
plot(L/100,NerrPVF(:,jerr),'Color',clr(6,:),'LineWidth',2) ; hold on ; 
plot(L/100,NerrPVG(:,jerr),'Color',clr(7,:),'LineWidth',2) ; hold on ;
xlim([20 1000]) ;
ylim([10 250]) ;
xlabel('Environment size (m)', 'Units','normalized','Position',[0.5 -0.11]);
ylabel(['Minimal N required for error < ' num2str(errs(jerr)) 'm'], 'Units','normalized','Position',[ -0.15 0.45]);
h = gca;
h.XRuler.TickLabelGapOffset = -1;
h.XRuler.TickLength = [0.02 0.02];
h.XTick = [20 200 400 600 800 1000];
% h.XScale = 'log';
% h.YScale = 'log';

%% panel B - bar plot of slopes
x = 1:5;
y = [
    lmNerrPVA{jerr}.Coefficients.Estimate(2)
    lmNerrPVB{jerr}.Coefficients.Estimate(2)
    lmNerrPVI{jerr}.Coefficients.Estimate(2)
    lmNerrPVF{jerr}.Coefficients.Estimate(2)
    lmNerrPVG{jerr}.Coefficients.Estimate(2)]';
CI = [
    lmNerrPVA{jerr}.coefCI
    lmNerrPVB{jerr}.coefCI
    lmNerrPVI{jerr}.coefCI
    lmNerrPVF{jerr}.coefCI
    lmNerrPVG{jerr}.coefCI];
CI  = CI(2:2:end,:)';
err = CI - y;

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

% add significance test
slopes_fit = {
    lmNerrPVA{jerr}
    lmNerrPVB{jerr}
    lmNerrPVI{jerr}
    lmNerrPVF{jerr}
    lmNerrPVG{jerr}};
signif_lines_y = [1.72 1.62 1.52 1.42];
for ii_scheme = 1:4
    % signif test
    mu1 = slopes_fit{ii_scheme}.Coefficients.Estimate(2);
    mu2 = slopes_fit{5}.Coefficients.Estimate(2);
    alpha = (1-normcdf(1,0,1));
    ci1 = coefCI(slopes_fit{ii_scheme}, alpha);
    ci2 = coefCI(slopes_fit{5}, alpha);
    sigma1 = ci1(2,2)-mu1;
    sigma2 = ci2(2,2)-mu2;
    sigma1 = sqrt(slopes_fit{ii_scheme}.CoefficientCovariance(2,2));
    sigma2 = sqrt(slopes_fit{5}.CoefficientCovariance(2,2));
    z = (mu1-mu2) / sqrt(sigma1^2+sigma2^2);
    pval = 1-normcdf(z,0,1);
    fprintf('scheme %d vs %d: pval=%g\n',ii_scheme,5,pval);
    % line + asterisks
    m = signif_lines_y(ii_scheme);
    X = [ii_scheme 5];
    Y = [1 1].*m;
    [Xf, Yf] = ds2nfu(X,Y);
    annotation('line',Xf,Yf,'LineWidth',1.1);
    if pval < 0.05
        signif_str = '*';
        signif_str_font_size = 12;
        signif_str_offset = 0;
    else
        signif_str = 'n.s.';
        signif_str_font_size = 6;
        signif_str_offset = 0.035;
    end
    text(mean(X),m+0.02+signif_str_offset, signif_str, 'HorizontalAlignment','center','FontSize',signif_str_font_size);
end


%% panel C - mean decoding error
jN_options = find(ismember(N, [50]));
ylimits_options = [0 20; 0 10];
for ii_N = 1:length(jN_options)
    jN = jN_options(ii_N);
    
    axes(panel_C(ii_N,1));
    cla
    hold on
    plot(L/100,mePVA(jN,:)*ds/100,'Color',clr(1,:),'LineWidth',2) ; hold on ; 
    plot(L/100,mePVB(jN,:)*ds/100,'Color',clr(3,:),'LineWidth',2) ; hold on ; 
    plot(L/100,mePVI(jN,:)*ds/100,'Color',clr(2,:),'LineWidth',2) ; hold on ; 
    plot(L/100,mePVF(jN,:)*ds/100,'Color',clr(6,:),'LineWidth',2) ; hold on ; 
    plot(L/100,mePVG(jN,:)*ds/100,'Color',clr(7,:),'LineWidth',2) ; hold on ; 
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
        text(-0.24,1.1, 'C', 'Units','normalized','FontWeight','bold');
        text(0.5,1.1, 'Mean decoding errors', 'Units','normalized','FontWeight','bold','HorizontalAlignment','center','FontSize',9);
    end

    % zoom-in inset
    axes(panel_C(ii_N,2));
    cla
    hold on
    plot(L/100,mePVA(jN,:)*ds/100,'Color',clr(1,:),'LineWidth',2) ; hold on ; 
    plot(L/100,mePVB(jN,:)*ds/100,'Color',clr(3,:),'LineWidth',2) ; hold on ; 
    plot(L/100,mePVI(jN,:)*ds/100,'Color',clr(2,:),'LineWidth',2) ; hold on ; 
    plot(L/100,mePVF(jN,:)*ds/100,'Color',clr(6,:),'LineWidth',2) ; hold on ; 
    plot(L/100,mePVG(jN,:)*ds/100,'Color',clr(7,:),'LineWidth',2) ; hold on ; 
    xlim([20 1000]) ;
    ylim(ylimits_options(ii_N,:))
    h = gca;
    h.XScale = 'log';
%     h.YScale = 'log';
    h.XTick = [20 100 1000];
    h.XTickLabel = {'20';'100';'1000'};
    h.XRuler.TickLabelGapOffset = -1;
    h.YRuler.TickLabelGapOffset = 0.5;
    h.XRuler.TickLength = [0.03 0.03];
    h.YRuler.TickLength = [0.015 0.015];
    h.XRuler.FontSize = 6.5;
    h.YRuler.FontSize = 6.5;
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
    plot(L/100,pePVA(jN,:,jprc)*ds/100,'Color',clr(1,:),'LineWidth',2) ; hold on ; 
    plot(L/100,pePVB(jN,:,jprc)*ds/100,'Color',clr(3,:),'LineWidth',2) ; hold on ; 
    plot(L/100,pePVI(jN,:,jprc)*ds/100,'Color',clr(2,:),'LineWidth',2) ; hold on ; 
    plot(L/100,pePVF(jN,:,jprc)*ds/100,'Color',clr(6,:),'LineWidth',2) ; hold on ; 
    plot(L/100,pePVG(jN,:,jprc)*ds/100,'Color',clr(7,:),'LineWidth',2) ; hold on ; 
%     ylim([0.1 2000]) ;
    xlim([20 1000]) ;
    xlabel('Environment size (m)', 'Units','normalized','Position',[0.5 -0.11]);
    ylabel([num2str(prc(jprc)) '% decoding error (m)'], 'Units','normalized','Position',[ -0.125 0.5]);
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
        text(0.5,1.1, 'Catastrophic errors: size', 'Units','normalized','FontWeight','bold','HorizontalAlignment','center','FontSize',9);
    end
    
    % zoom-in inset
    axes(panel_D(ii_N,2));
    cla
    hold on
    plot(L/100,pePVA(jN,:,jprc)*ds/100,'Color',clr(1,:),'LineWidth',2) ; hold on ; 
    plot(L/100,pePVB(jN,:,jprc)*ds/100,'Color',clr(3,:),'LineWidth',2) ; hold on ; 
    plot(L/100,pePVI(jN,:,jprc)*ds/100,'Color',clr(2,:),'LineWidth',2) ; hold on ; 
    plot(L/100,pePVF(jN,:,jprc)*ds/100,'Color',clr(6,:),'LineWidth',2) ; hold on ; 
    plot(L/100,pePVG(jN,:,jprc)*ds/100,'Color',clr(7,:),'LineWidth',2) ; hold on ; 
%     ylim([0.1 2000]) ;
    xlim([100 1000]) ;
    h = gca;
    h.XScale = 'log';
    h.YScale = 'log';
    h.XTick = [100 1000];
    h.XTickLabel = {'100';'1000'};
    h.YTick = [1 10 100];
    h.YTickLabel = {'1';'10';'100'};
    h.XRuler.TickLabelGapOffset = -1;
    h.YRuler.TickLabelGapOffset = 0.5;
    h.XRuler.TickLength = [0.06 0.06];
    h.YRuler.TickLength = [0.04 0.04];
    h.XRuler.FontSize = 6.5;
    h.YRuler.FontSize = 6.5;
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

    plot(L/100,plePVA(jN,:)+(2e6)^(-1),'Color',clr(1,:),'LineWidth',2) ; hold on ; 
    plot(L/100,plePVB(jN,:)+(2e6)^(-1),'Color',clr(3,:),'LineWidth',2) ; hold on ; 
    plot(L/100,plePVI(jN,:)+(2e6)^(-1),'Color',clr(2,:),'LineWidth',2) ; hold on ; 
    plot(L/100,plePVF(jN,:)+(2e6)^(-1),'Color',clr(6,:),'LineWidth',2) ; hold on ; 
    plot(L/100,plePVG(jN,:)+(2e6)^(-1),'Color',clr(7,:),'LineWidth',2) ; hold on ; 
    xlim([20 1000]) ;
    ylim(ylimits_options(ii_N,:))
    xlabel('Environment size (m)', 'Units','normalized','Position',[0.5 -0.11]);
    ylabel('Prob.(error > 5% of environment size)', 'Units','normalized','Position',[ -0.14 0.4]);
%     text(0.95,1.03, sprintf('N = %d',N(jN)), 'Units','normalized','FontWeight','normal', 'HorizontalAlignment','right','FontSize',10);
    h = gca;
    h.XScale = 'log';
    h.YScale = 'log';
    h.XTick = [20 100 1000];
    h.XTickLabel = {'20';'100';'1000'};
    h.YTick = 10.^[-3 -2 -1 0];
    h.YTickLabel = {'10^{ -3}'; '10^{ -2}'; '10^{ -1}'; '10^{ 0}'};
    ylim(10.^[-3 0])
    h.YRuler.TickLabelGapOffset = -0.1;
    h.XRuler.TickLabelGapOffset = -1;
    h.XRuler.TickLength = [0.03 0.03];
    h.YRuler.TickLength = [0.03 0.03];
    if ii_N == 1
        text(-0.20,1.1, 'E', 'Units','normalized','FontWeight','bold');
        text(0.5,1.1, 'Catastrophic errors: probability', 'Units','normalized','FontWeight','bold','HorizontalAlignment','center','FontSize',9);
    end
end


%% print/save the figure
fig_name_out = fullfile(res_dir, sprintf('%s_dt_%s_coverage=%d', ...
    fig_name_str, ...
    strrep(num2str(dt),'.','_'), ...
    100*coverage));
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');




%%


