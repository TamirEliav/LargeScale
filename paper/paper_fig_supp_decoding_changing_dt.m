%% Large Scale - fig. S18 - decoder with shorter dt (200ms) + error vs. dt

%%
% clear
clearvars -except f
clc

%% params
use_absolute_error = 1;     % 0 - 5% of env. error; 1 - 2m error

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S19';
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
panel_A_size = [3 3].*[0.85 0.80];
panel_BCDE_size = [4.6 4];
panel_A(1) = axes('position', [ 2.0 21 panel_A_size]);
panel_A(2) = axes('position', [ 5.1 21 panel_A_size]);
panel_A(3) = axes('position', [ 8.2 21 panel_A_size]);
panel_A(4) = axes('position', [11.3 21 panel_A_size]);
panel_A(5) = axes('position', [14.4 21 panel_A_size]);
panel_A(6) = axes('position', [17.5 21 panel_A_size]);
panel_B(1)   = axes('position', [ 2.0 14.6 panel_BCDE_size]);
panel_B(2)   = axes('position', [ 7.6 14.6 2.6 4]);
panel_C(1,1) = axes('position', [ 2.0  8.7  panel_BCDE_size]);
panel_C(1,2) = axes('position', [ 2.5   10.7  2.8 2]);
panel_D(1,1) = axes('position', [ 8.1  8.7  panel_BCDE_size]);
panel_D(1,2) = axes('position', [ 8.7 10.7  2.5 2]);
panel_E(1)   = axes('position', [14.4    8.7  panel_BCDE_size]);
panel_F      = axes('position', [ 4   2 6 4]);
panel_G      = axes('position', [12.5 2 6 4]);
panel_legend = axes('position', [9.8 16.6 0.4 2]);

%% color variable for different schemes
clr = [1.0 0.0 0.0 ; ...
       1.0 0.8 0.0 ; ...
       1.0 0.0 0.8 ; ...
       0.0 0.7 0.0 ; ...
       0.2 0.2 1.0 ; ...
       0.5 0.0 0.5 ] ;


%% arrange data
paper_fig_7_arrange_data;

%% generate examples
exampleL = 200*100; % 200*cm
pos = 1:ds:exampleL;
pos = pos / 100; % back to meter
if ~exist('f','var')
    rng(0);
    [f] = MultiscalePlace_GenerateTuning_AllModels(exampleL);
end

%% choose dt
panels_AE_dt = 0.2;
jdt = find(ismember(dt, panels_AE_dt));

%% panel A - different encoding schemes
maps=[];
maps(1,:,:) = f(1).maps(:,1:10);
maps(1,:,[7 4 9 10]) = f(1).maps(:,[4 7 11 12]);
maps(2,:,:) = f(2).maps(:,1:10);
maps(2,:,[3 9 5]) = f(2).maps(:,[14 12 15]);
maps(3,:,:) = f(3).maps(:,950:-100:50);
maps(3,:,1) = f(3).maps(:,955);
maps(3,:,2) = f(3).maps(:,854);
maps(3,:,10) = f(3).maps(:,52);
[~,IX] = sort(sum(f(4).maps,1));
maps(4,:,:) = f(4).maps(:,IX(50:100:950));
maps(4,:,10) = f(4).maps(:,IX(951));
rng(0);
maps(4,:,:) = maps(4,:,randperm(10));
maps(5,:,:) = f(5).maps(:,50:100:950);
maps(5,:,2) = f(5).maps(:,155);
maps(5,:,3) = f(5).maps(:,253);
maps(5,:,4) = f(5).maps(:,353);
maps(5,:,5) = f(5).maps(:,452);
maps(6,:,:) = f(6).maps(:,1:10);
maps(6,:,10) = f(6).maps(:,952);
for ii_scm = 1:6
    axes(panel_A(ii_scm));
    cla
    hold on
    imagesc(pos,1:10,squeeze(maps(ii_scm,:,:))');
    cmap = gray();
    cmap = flip(cmap);
    colormap(gca,cmap);
    xlim([0 200])
    ylim([0.5 10.5])
    xlabel('Position (m)', 'Units','normalized','Position',[0.5 -0.17]);
    h = gca;
    h.YTick = [1 5 10];
    h.XRuler.TickLabelGapOffset = -1;
    h.YRuler.TickLabelGapOffset = -0.5;
    h.XRuler.TickLength = [0.03 0.03];
    ht=title("Scheme "+ii_scm,'Units','normalized','Position',[0.5 1.06]);
    schemes_IX = [1:6];
    schemes_title_color_opt = 2;
    switch schemes_title_color_opt
        case 1
            ht.Color = clr(schemes_IX(ii_scm),:);
        case 2
            plot([0 25]+14, 11.8*[1 1], 'Color', clr(schemes_IX(ii_scm),:), 'LineWidth', 2, 'Clipping','off');
    end
end
axes(panel_A(1));
text(-0.4,1.2, 'A', 'Units','normalized','FontWeight','bold');
text(3,1.5, 'Integration time window of  {\Delta}t = 200 ms instead of 500 ms', 'Units','normalized','FontWeight','bold','FontSize',12,'HorizontalAlignment','center');
ylabel('Example neuron no.','Units','normalized','Position',[-0.18 0.5]);

%% panel B - minimum N required for error <2m (or <5%)
axes(panel_B(1));
cla
hold on
text(-0.24,1.13, 'B', 'Units','normalized','FontWeight','bold');
text(0.5,1.13, {'Minimal no. of neurons (N)';'required for decoding'}, 'Units','normalized','FontWeight','bold','HorizontalAlignment','center','FontSize',9);

if use_absolute_error
    plot(L,Nerr_ML_S1(:,jdt),'Color',clr(1,:),'LineWidth',2) ;
    plot(L,Nerr_ML_S3(:,jdt),'Color',clr(3,:),'LineWidth',2) ;
    plot(L,Nerr_ML_S2(:,jdt),'Color',clr(2,:),'LineWidth',2) ;
    plot(L,Nerr_ML_S4(:,jdt),'Color',clr(4,:),'LineWidth',2) ;
    plot(L,Nerr_ML_S5(:,jdt),'Color',clr(5,:),'LineWidth',2) ;
    plot(L,Nerr_ML_S6(:,jdt),'Color',clr(6,:),'LineWidth',2) ;
    plot(L,Nerr_ML_S5v(:,jdt),'Color',clr(5,:),'LineWidth',2,'LineStyle','--') ;
    plot(L,Nerr_ML_S6v(:,jdt),'Color',clr(6,:),'LineWidth',2,'LineStyle','--') ;
else
    plot(L,Nrerr_ML_S1(:,jdt),'Color',clr(1,:),'LineWidth',2) ;
    plot(L,Nrerr_ML_S3(:,jdt),'Color',clr(3,:),'LineWidth',2) ;
    plot(L,Nrerr_ML_S2(:,jdt),'Color',clr(2,:),'LineWidth',2) ;
    plot(L,Nrerr_ML_S4(:,jdt),'Color',clr(4,:),'LineWidth',2) ;
    plot(L,Nrerr_ML_S5(:,jdt),'Color',clr(5,:),'LineWidth',2) ;
    plot(L,Nrerr_ML_S6(:,jdt),'Color',clr(6,:),'LineWidth',2) ;
    plot(L,Nrerr_ML_S5v(:,jdt),'Color',clr(5,:),'LineWidth',2,'LineStyle','--') ;
    plot(L,Nrerr_ML_S6v(:,jdt),'Color',clr(6,:),'LineWidth',2,'LineStyle','--') ;
end
xlim([20 1000]) ;
ylim([10 max(N)]) ;
xlabel('Environment size (m)', 'Units','normalized','Position',[0.5 -0.11]);
if use_absolute_error
    ylabel(['Minimal N required for error < ' num2str(errs) 'm'], 'Units','normalized','Position',[ -0.15 0.45]);
else
    ylabel(['Minimal N required for error < ' num2str(100*rerr) '%'], 'Units','normalized','Position',[ -0.15 0.45]);
end

h = gca;
h.XRuler.TickLabelGapOffset = -1;
h.XRuler.TickLength = [0.02 0.02];
h.XTick = [20 200 400 600 800 1000];
% h.XScale = 'log';
% h.YScale = 'log';

%% panel B - bar plot of slopes
if use_absolute_error
    lms = {
        lmNerr_ML_S1{jdt}
        lmNerr_ML_S2{jdt}
        lmNerr_ML_S3{jdt}
        lmNerr_ML_S4{jdt}
        lmNerr_ML_S5{jdt}
        lmNerr_ML_S6{jdt}
        lmNerr_ML_S5v{jdt}
        lmNerr_ML_S6v{jdt}};
else
    lms = {
        lmNrerr_ML_S1{jdt}
        lmNrerr_ML_S2{jdt}
        lmNrerr_ML_S3{jdt}
        lmNrerr_ML_S4{jdt}
        lmNrerr_ML_S5{jdt}
        lmNrerr_ML_S6{jdt}
        lmNrerr_ML_S5v{jdt}
        lmNrerr_ML_S6v{jdt}};
end

switch length(lms{1}.CoefficientNames)
    case 1
        x = 1:length(lms);
        y  = cellfun(@(lm)(lm.Coefficients.Estimate), lms)';
        CI = cellfun(@(lm)(lm.coefCI),                lms,'UniformOutput',0) ;
        CI = cat(1,CI{:})';
        err = CI - y;
    case 2
        x = 1:length(lms);
        y  = cellfun(@(lm)(lm.Coefficients.Estimate(2)), lms)';
        CI = cellfun(@(lm)(lm.coefCI),                   lms,'UniformOutput',0) ;
        CI = cat(1,CI{:});
        CI = CI(2:2:end,:)';
        err = CI - y;
end
% plot
axes(panel_B(2));
cla
hold on
hb=bar(x,y);
hb.FaceColor = 'flat';
hb.CData = [clr; 1 1 1; 1 1 1];
he=errorbar(x,y,err(2,:));
he.CapSize = 2;
he.LineStyle = 'none';
he.Color = 'k';
m = 1.1 * max(y+err(2,:));
m = round(m,1);
xlim([0 length(lms)+1]);
ylim([0 m]);
h=gca;
h.XTick = [];
h.YTick = [0 m];
ylabel('Slope (N per meter)', 'Units','normalized','Position',[-0.10 0.5]);
% fill bar with dots for schemes 5v & 6v
for ii=7:8
    barw = 0.8; % bar width
    spacing_x = 0.3;
    spacing_y = 0.05;
    switch 3
        case 1
            xx = linspace(ii-barw/2,ii+barw/2,12);
            xx(1:4:end)=nan;
            xx(2:4:end)=nan;
            yy = linspace(0,y(ii),5)';
            yy = repmat(yy,1,length(xx));
        case 2
            %%
            xx = (ii-barw/2) : spacing_x : (ii+barw/2);
            yy = 0 : spacing_y : y(ii);
            xx=xx';
%             yy=yy';
            yy = repmat(yy',1,length(xx))';
            whos xx yy
        case 3
            %%
            xx = linspace(ii-barw/2,ii+barw/2,4)';
            yy = linspace(0,y(ii),5)';
            xx2 = xx + (xx(2)-xx(1))/2;
            yy2 = yy + (yy(2)-yy(1))/2;
            xx2(end)=[];
            yy2(end)=[];
            yy = repmat(yy,1,length(xx))';
            yy2 = repmat(yy2,1,length(xx2))';
    end
    plot(xx,yy,'.','Color', clr(ii-2,:),'markerSize',3);
    plot(xx2,yy2,'.','Color', clr(ii-2,:),'markerSize',3);
    plot([1 2 3],[7 7 7; 8 8 8 ;9 9 9 ; 10 10 10])
end

% add significance test
signif_lines_x = flip(linspace(5.7,6.3,5));
signif_lines_y = flip(linspace(1.8,2.35,5));
if ~use_absolute_error
    signif_lines_y = 0.5 .* signif_lines_y;
end
scm_test_IX = 6;
for ii_scm = 1:scm_test_IX-1
    % signif test
    switch length(lms{1}.CoefficientNames)
        case 1
            mu1 = lms{ii_scm}.Coefficients.Estimate;
            mu2 = lms{scm_test_IX}.Coefficients.Estimate;
            sigma1 = sqrt(lms{ii_scm}.CoefficientCovariance);
            sigma2 = sqrt(lms{scm_test_IX}.CoefficientCovariance);
        case 2
            mu1 = lms{ii_scm}.Coefficients.Estimate(2);
            mu2 = lms{scm_test_IX}.Coefficients.Estimate(2);
%             alpha = (1-normcdf(1,0,1));
%             ci1 = coefCI(lms{ii_scm}, alpha);
%             ci2 = coefCI(lms{scm_test_IX}, alpha);
            sigma1 = sqrt(lms{ii_scm}.CoefficientCovariance(2,2));
            sigma2 = sqrt(lms{scm_test_IX}.CoefficientCovariance(2,2));
    end
    z = (mu1-mu2) / sqrt(sigma1^2+sigma2^2);
    pval = 1-normcdf(z,0,1);
    fprintf('scheme %d vs %d: pval=%g\n',ii_scm,scm_test_IX,pval);
    % line + asterisks
    m = signif_lines_y(ii_scm);
    X = [ii_scm signif_lines_x(ii_scm)];
    Y = [1 1].*m;
    plot(X,Y, 'k-', 'LineWidth',0.1, 'Clipping','off');
    signif_str = signif_astricks(pval);
    if pval < 0.05
        signif_str_font_size = 12;
        signif_str_offset = 0;
    else
        signif_str_font_size = 6;
        signif_str_offset = 0.05;
    end
    text(mean(X),m+0.05+signif_str_offset, signif_str, 'HorizontalAlignment','center','FontSize',signif_str_font_size);
    X = [signif_lines_x(ii_scm) signif_lines_x(ii_scm)];
    Y = [0.5 m];
    plot(X,Y, 'k-', 'LineWidth',0.1, 'Clipping','off');
end

%% legend panel
axes(panel_legend);
cla
hold on
axis ij
for ii_scm = 1:6
    plot([0 1], ii_scm*[1 1], 'Color',clr(ii_scm,:),'LineWidth',2);
    text(1.2,   ii_scm, f(ii_scm).name_long, 'FontSize', 9);
end
dashed_line = linspace(0,1,11);
dashed_line([4 8]) = nan;
plot(dashed_line, 7*ones(size(dashed_line)), 'Color',clr(5,:),'LineWidth',2, 'Clipping','off');
plot(dashed_line, 8*ones(size(dashed_line)), 'Color',clr(6,:),'LineWidth',2, 'Clipping','off');
text(1.2,   7, f(7).name_long, 'FontSize', 9);
text(1.2,   8, f(8).name_long, 'FontSize', 9);
axis off
xlim([0 1])
ylim([1 6])

%% panel C - mean decoding error
jN_options = find(ismember(N, [50]));
ylimits_options = [0 20; 0 10];
for ii_N = 1:length(jN_options)
    jN = jN_options(ii_N);

    axes(panel_C(ii_N,1));
    cla
    hold on
    plot(L,merr_ML_S1(jN,:,jdt)*ds/100,'Color',clr(1,:),'LineWidth',2);
    plot(L,merr_ML_S2(jN,:,jdt)*ds/100,'Color',clr(2,:),'LineWidth',2);
    plot(L,merr_ML_S3(jN,:,jdt)*ds/100,'Color',clr(3,:),'LineWidth',2);
    plot(L,merr_ML_S4(jN,:,jdt)*ds/100,'Color',clr(4,:),'LineWidth',2);
    plot(L,merr_ML_S5(jN,:,jdt)*ds/100,'Color',clr(5,:),'LineWidth',2);
    plot(L,merr_ML_S6(jN,:,jdt)*ds/100,'Color',clr(6,:),'LineWidth',2);
    plot(L,merr_ML_S5v(jN,:,jdt)*ds/100,'Color',clr(5,:),'LineWidth',2,'LineStyle','--');
    plot(L,merr_ML_S6v(jN,:,jdt)*ds/100,'Color',clr(6,:),'LineWidth',2,'LineStyle','--');
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
    plot(L,merr_ML_S1(jN,:,jdt)*ds/100,'Color',clr(1,:),'LineWidth',2);
    plot(L,merr_ML_S2(jN,:,jdt)*ds/100,'Color',clr(2,:),'LineWidth',2);
    plot(L,merr_ML_S3(jN,:,jdt)*ds/100,'Color',clr(3,:),'LineWidth',2);
    plot(L,merr_ML_S4(jN,:,jdt)*ds/100,'Color',clr(4,:),'LineWidth',2);
    plot(L,merr_ML_S5(jN,:,jdt)*ds/100,'Color',clr(5,:),'LineWidth',2);
    plot(L,merr_ML_S6(jN,:,jdt)*ds/100,'Color',clr(6,:),'LineWidth',2);
    plot(L,merr_ML_S5v(jN,:,jdt)*ds/100,'Color',clr(5,:),'LineWidth',2,'LineStyle','--');
    plot(L,merr_ML_S6v(jN,:,jdt)*ds/100,'Color',clr(6,:),'LineWidth',2,'LineStyle','--');
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

%% panel D - catastrophic error - Size (99th prc error)
prc = 99;
jN_options = find(ismember(N, [50]));
ylimits_options = [0 20; 0 10];
for ii_N = 1:length(jN_options)
    jN = jN_options(ii_N);

    axes(panel_D(ii_N,1));
    cla
    hold on
    plot(L,perr_ML_S1(jN,:,jdt)*ds/100,'Color',clr(1,:),'LineWidth',2);
    plot(L,perr_ML_S2(jN,:,jdt)*ds/100,'Color',clr(2,:),'LineWidth',2);
    plot(L,perr_ML_S3(jN,:,jdt)*ds/100,'Color',clr(3,:),'LineWidth',2);
    plot(L,perr_ML_S4(jN,:,jdt)*ds/100,'Color',clr(4,:),'LineWidth',2);
    plot(L,perr_ML_S5(jN,:,jdt)*ds/100,'Color',clr(5,:),'LineWidth',2);
    plot(L,perr_ML_S6(jN,:,jdt)*ds/100,'Color',clr(6,:),'LineWidth',2);
    plot(L,perr_ML_S5v(jN,:,jdt)*ds/100,'Color',clr(5,:),'LineWidth',2,'LineStyle','--');
    plot(L,perr_ML_S6v(jN,:,jdt)*ds/100,'Color',clr(6,:),'LineWidth',2,'LineStyle','--');
%     ylim([0.1 2000]) ;
    xlim([20 1000]) ;
    xlabel('Environment size (m)', 'Units','normalized','Position',[0.5 -0.11]);
    ylabel([num2str(prc) '% decoding error (m)'], 'Units','normalized','Position',[ -0.145 0.5]);
    h = gca;
    h.XScale = 'log';
%     h.YScale = 'log';
    h.XTick = [20 100 1000];
    h.XTickLabel = {'20';'100';'1000'};
    h.YRuler.TickLabelGapOffset = 1.5;
    h.XRuler.TickLabelGapOffset = -1;
    h.XRuler.TickLength = [0.03 0.03];
%     text(0.95,1, sprintf('N = %d',N(jN)), 'Units','normalized','FontWeight','normal', 'HorizontalAlignment','right','FontSize',10);
    if ii_N == 1
        text(-0.21,1.1, 'D', 'Units','normalized','FontWeight','bold');
        text(0.5,1.1, 'Catastrophic errors: size', 'Units','normalized','FontWeight','bold','HorizontalAlignment','center','FontSize',9);
    end

    % zoom-in inset
    axes(panel_D(ii_N,2));
    cla
    hold on
    plot(L,perr_ML_S1(jN,:,jdt)*ds/100,'Color',clr(1,:),'LineWidth',2);
    plot(L,perr_ML_S2(jN,:,jdt)*ds/100,'Color',clr(2,:),'LineWidth',2);
    plot(L,perr_ML_S3(jN,:,jdt)*ds/100,'Color',clr(3,:),'LineWidth',2);
    plot(L,perr_ML_S4(jN,:,jdt)*ds/100,'Color',clr(4,:),'LineWidth',2);
    plot(L,perr_ML_S5(jN,:,jdt)*ds/100,'Color',clr(5,:),'LineWidth',2);
    plot(L,perr_ML_S6(jN,:,jdt)*ds/100,'Color',clr(6,:),'LineWidth',2);
    plot(L,perr_ML_S5v(jN,:,jdt)*ds/100,'Color',clr(5,:),'LineWidth',2,'LineStyle','--');
    plot(L,perr_ML_S6v(jN,:,jdt)*ds/100,'Color',clr(6,:),'LineWidth',2,'LineStyle','--');
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

%% panel E - catastrophic error - probability(error>5%)
jN_options = find(ismember(N, [50]));
ylimits_options = [2e-4 1; 2e-7 1];
for ii_N = 1:length(jN_options)
    jN = jN_options(ii_N);

    axes(panel_E(ii_N));
    cla
    hold on

    prob_offset = (2e6)^(-1); % at least in one simulation there was a catastrophic error
    prob_offset = 0;
    plot(L,plerr_ML_S1(jN,:,jdt)+prob_offset,'Color',clr(1,:),'LineWidth',2);
    plot(L,plerr_ML_S2(jN,:,jdt)+prob_offset,'Color',clr(2,:),'LineWidth',2);
    plot(L,plerr_ML_S3(jN,:,jdt)+prob_offset,'Color',clr(3,:),'LineWidth',2);
    plot(L,plerr_ML_S4(jN,:,jdt)+prob_offset,'Color',clr(4,:),'LineWidth',2);
    plot(L,plerr_ML_S5(jN,:,jdt)+prob_offset,'Color',clr(5,:),'LineWidth',2);
    plot(L,plerr_ML_S6(jN,:,jdt)+prob_offset,'Color',clr(6,:),'LineWidth',2);
    plot(L,plerr_ML_S5v(jN,:,jdt)+prob_offset,'Color',clr(5,:),'LineWidth',2,'LineStyle','--');
    plot(L,plerr_ML_S6v(jN,:,jdt)+prob_offset,'Color',clr(6,:),'LineWidth',2,'LineStyle','--');
    xlim([20 1000]) ;
    ylim(ylimits_options(ii_N,:))
    xlabel('Environment size (m)', 'Units','normalized','Position',[0.5 -0.11]);
    ylabel('Prob.(error > 5% of environment size)', 'Units','normalized','Position',[ -0.18 0.4]);
%     text(0.95,1.03, sprintf('N = %d',N(jN)), 'Units','normalized','FontWeight','normal', 'HorizontalAlignment','right','FontSize',10);
    h = gca;
    h.XScale = 'log';
    h.YScale = 'log';
    h.XTick = [20 100 1000];
    h.XTickLabel = {'20';'100';'1000'};
%     h.YTick = 10.^[-3 -2 -1 0];
%     h.YTickLabel = {'10^{ -3}'; '10^{ -2}'; '10^{ -1}'; '10^{ 0}'};
    h.YTick = 10.^[-4 -3 -2 -1 0];
    h.YTickLabel = {'10^{ -4}'; '10^{ -3}'; '10^{ -2}'; '10^{ -1}'; '10^{ 0}'};
    h.YLim = 10.^[-4 0];
    h.YRuler.TickLabelGapOffset = -0.1;
    h.XRuler.TickLabelGapOffset = -1;
    h.XRuler.TickLength = [0.03 0.03];
    h.YRuler.TickLength = [0.024 0.03];
    if ii_N == 1
        text(-0.24,1.1, 'E', 'Units','normalized','FontWeight','bold');
        text(0.5,1.1, 'Catastrophic errors: probability', 'Units','normalized','FontWeight','bold','HorizontalAlignment','center','FontSize',9);
    end
end


%% ========================================================================
% From here we are using different simulations resutls, with variable dt
% load data
load('L:\Yonatan_theory\20200710__simulations_updated\SummaryDecoderResults_VariableIntegrationTime_delta0.3.mat');
L = Lv*ds/100 ; % environment size variable in meters

%% panel F - decoding mean error vs. dt
axes(panel_F);
cla
hold on
text(-0.2,1.05, 'F', 'Units','normalized','FontWeight','bold');
text(1,1.25, 'Varying the integration time window  {\Delta}t', 'Units','normalized','FontWeight','bold','FontSize',12,'HorizontalAlignment','center');

jN = find(ismember(N, [50]));
jL = find(ismember(L, [200]));
plot(dt.*1e3, merr_ML_S1(:,jL,jN)*ds/100, 'Color', clr(1,:), 'LineWidth', 2);
plot(dt.*1e3, merr_ML_S2(:,jL,jN)*ds/100, 'Color', clr(2,:), 'LineWidth', 2);
plot(dt.*1e3, merr_ML_S3(:,jL,jN)*ds/100, 'Color', clr(3,:), 'LineWidth', 2);
plot(dt.*1e3, merr_ML_S4(:,jL,jN)*ds/100, 'Color', clr(4,:), 'LineWidth', 2);
plot(dt.*1e3, merr_ML_S5(:,jL,jN)*ds/100, 'Color', clr(5,:), 'LineWidth', 2);
plot(dt.*1e3, merr_ML_S6(:,jL,jN)*ds/100, 'Color', clr(6,:), 'LineWidth', 2);
plot(dt.*1e3, merr_ML_S5v(:,jL,jN)*ds/100, 'Color', clr(5,:), 'LineWidth', 2,'LineStyle','--');
plot(dt.*1e3, merr_ML_S6v(:,jL,jN)*ds/100, 'Color', clr(6,:), 'LineWidth', 2,'LineStyle','--');

ha= gca;
ha.XLim = [20 500];
ha.YLim = [0.3 100];
ha.XTick = [20 100:100:500];
ha.YTick = 10.^[0 1 2 3];
ha.YTickLabel = {'10^{ 0}'; '10^{ 1}'; '10^{ 2}'; '10^{ 3}'};
ha.TickDir='out';
ha.TickLength = [0.02 0.02];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;
xlabel('Integration time window (ms)')
ylabel('Mean decoding error (m)','Units','normalized','Position',[-0.12 0.5]);
ha.XScale = 'linear';
ha.YScale = 'log';

%% panel G - Catastrophic errors vs. dt
axes(panel_G);
cla
hold on
text(-0.25,1.05, 'G', 'Units','normalized','FontWeight','bold');

jN = find(ismember(N, [50]));
jL = find(ismember(L, [200]));
plot(dt.*1e3, plerr_ML_S1(:,jL,jN), 'Color', clr(1,:), 'LineWidth', 2);
plot(dt.*1e3, plerr_ML_S2(:,jL,jN), 'Color', clr(2,:), 'LineWidth', 2);
plot(dt.*1e3, plerr_ML_S3(:,jL,jN), 'Color', clr(3,:), 'LineWidth', 2);
plot(dt.*1e3, plerr_ML_S4(:,jL,jN), 'Color', clr(4,:), 'LineWidth', 2);
plot(dt.*1e3, plerr_ML_S5(:,jL,jN), 'Color', clr(5,:), 'LineWidth', 2);
plot(dt.*1e3, plerr_ML_S6(:,jL,jN), 'Color', clr(6,:), 'LineWidth', 2);
plot(dt.*1e3, plerr_ML_S5v(:,jL,jN), 'Color', clr(5,:), 'LineWidth', 2,'LineStyle','--');
plot(dt.*1e3, plerr_ML_S6v(:,jL,jN), 'Color', clr(6,:), 'LineWidth', 2,'LineStyle','--');

ha= gca;
ha.XLim = [20 500];
ha.YLim = [0.8e-3 1e0];
ha.XTick = [20 100:100:500];
ha.YTick = 10.^[-3 -2 -1 0];
ha.YTickLabel = {'10^{ -3}'; '10^{ -2}'; '10^{ -1}'; '10^{ 0}'};
ha.TickDir='out';
ha.TickLength = [0.02 0.02];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;
xlabel('Integration time window (ms)')
ylabel('Prob.(error > 5% of environment size)','Units','normalized','Position',[-0.14 0.4]);
ha.XScale = 'linear';
ha.YScale = 'log';

%% print/save the figure
fig_name = sprintf('%s_dt_%s', ...
    fig_name_str, ...
    strrep(num2str(panels_AE_dt),'.','_'));
file_out = fullfile(res_dir, fig_name);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');



%%

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
end


%%
