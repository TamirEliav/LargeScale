%% Large Scale - fig. S16- Theoretical analysis (decoding) - including schemes 5v,6v

%%
% clear
clearvars -except f
clc

%% params
use_absolute_error = 1;

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S16';
fig_caption_str = 'Theoretical analysis_decoding';
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
% panel_F      = axes('position', [   2    3.0  panel_BCDE_size]);
panel_F      = axes('position', [   2    3.8  3 3]);
panel_G      = axes('position', [   7    3.8  3 3]);
panel_H      = axes('position', [  14.4  3.2  panel_BCDE_size]);
panel_legend = axes('position', [9.8 16.6 0.4 2]);
panel_F_legend = axes('position', [3.2 6  0.2 0.4]);
panel_G_legend = axes('position', [9.5 6  0.2 0.75]);

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
panels_AE_dt = 0.5;
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
ylabel('Example neuron no.','Units','normalized','Position',[-0.18 0.5]);

%% panel B - minimum N required for error <2m (or <5%)
axes(panel_B(1));
cla
hold on
text(-0.24,1.13, 'B', 'Units','normalized','FontWeight','bold');
text(0.5,1.13, {'Minimal no. of neurons';'required for decoding'}, 'Units','normalized','FontWeight','bold','HorizontalAlignment','center','FontSize',9);

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
signif_lines_y = [1.72 1.62 1.52 1.42 1.32];
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
        signif_str_offset = 0.035;
    end
    text(mean(X),m+0.02+signif_str_offset, signif_str, 'HorizontalAlignment','center','FontSize',signif_str_font_size);
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
    h.YTick = 10.^[-3 -2 -1 0];
    h.YTickLabel = {'10^{ -3}'; '10^{ -2}'; '10^{ -1}'; '10^{ 0}'};
    ylim([0.7e-3 1])
    h.YRuler.TickLabelGapOffset = -0.1;
    h.XRuler.TickLabelGapOffset = -1;
    h.XRuler.TickLength = [0.03 0.03];
    h.YRuler.TickLength = [0.024 0.03];
    if ii_N == 1
        text(-0.24,1.1, 'E', 'Units','normalized','FontWeight','bold');
        text(0.5,1.1, 'Catastrophic errors: probability', 'Units','normalized','FontWeight','bold','HorizontalAlignment','center','FontSize',9);
    end
end

%% panel F - coverage exp fit
axes(panel_F);
cla
hold on
load('L:\processed_data_structs\total_area.mat');
h=histogram(total_area,'Normalization','pdf');
h.Normalization = 'pdf';
h.FaceColor = 0.5*[1 1 1];
h.FaceAlpha = 1;
cvrg=expfit(total_area);
plot(h.BinEdges,exppdf(h.BinEdges,cvrg),'LineWidth',2,'Color','k');
hax = gca;
hax.XLim = [0 160];
hax.XTick = [0:50:200];
hax.YTick = hax.YLim;
hax.YRuler.TickLabelGapOffset = -0.1;
hax.XRuler.TickLabelGapOffset = -1;
hax.XRuler.TickLength = [0.03 0.03];
hax.YRuler.TickLength = [0.024 0.03];
xlabel('Coverage (m)', 'Units','normalized','Position',[0.5 -0.14]);
ylabel({'Probability';'density function'}, 'Units','normalized','Position',[ -0.12 0.5]);
text(-0.35,1.15, 'F', 'Units','normalized','FontWeight','bold');

% add legend
axes(panel_F_legend);
cla('reset');
hold on
patch([1 1 2 2], 2*[1 1 1 1]+.3*[-1 1 1 -1], 0.5*[1 1 1],'EdgeColor','k');
plot([1 2],      1*[1 1], 'k','LineWidth',2);
text(2.6, 2, 'Data','FontSize',7,'HorizontalAlignment','left');
text(2.6, 1, 'Exponential fit','FontSize',7,'HorizontalAlignment','left');
hax=gca;
hax.Visible='off';

%% load data for panel G
% load data
load('L:\processed_data_structs\cells_bat_200m.mat');
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
% load simulations ratio LS
% load('C:\Tamir\work\Projects\LargeScale\Yonatan_theory\20200806__model_ratio_LS\Model_S6_S6v_RatioLargestSmallestField.mat');
load('L:\Yonatan_theory\20200816__model_ratio_LS_updated_remoevd_1_field_cells\Model_S6_S6v_RatioLargestSmallestField.mat')
rLS(K==1)=[];
rLSv(Kv==1)=[];

%% panel G - ratio L/S model+data
axes(panel_G);
cla
hold on
text(-0.35,1.15, 'G', 'Units','normalized','FontWeight','bold');

nBinEdges = 9;
edges = logspace(0,log10(25),nBinEdges);
% plot data
h=histogram(ratio_LS);
h.Normalization = 'probability';
h.BinEdges = edges;
h.FaceColor = 0.5*[1 1 1];
h.FaceAlpha = 1;
% plot model (S6)
h=histogram(rLS);
h.Normalization = 'probability';
h.BinEdges = edges;
h.DisplayStyle = 'stairs';
h.EdgeColor = clr(6,:);
h.FaceAlpha = 1;
h.LineWidth = 2;
% plot model (S6v)
h=histogram(rLSv);
h.Normalization = 'probability';
h.BinEdges = edges;
h.DisplayStyle = 'stairs';
h.EdgeColor = clr(6,:);
h.FaceAlpha = 1;
h.LineWidth = 2;
h.LineStyle = '-.';

hax=gca;
hax.XScale = 'log';
hax.YScale = 'log';
hax.XLim = [0 27];
% hax.YLim = [0.07  1.05 * hax.YLim(2)];
hax.YLim = [0.005 0.5];
% hax.YLim = [0.7 200];
hax.XTick = [1 2 5 10 20];
hax.YTick = hax.YLim;
hax.YTick = [0.005 0.05 0.5];
hax.TickDir='out';
hax.TickLength = [0.03 0.03];
hax.XRuler.TickLabelGapMultiplier = -0.35;
hax.YRuler.TickLabelGapMultiplier = 0.001;
xlabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[0.5 -0.17]);
ylabel('Fraction of cells','Units','normalized','Position',[-0.25 0.5])

% add legend
axes(panel_G_legend);
cla('reset');
hold on
dashed_line = linspace(1,2,11);
dashed_line([ 5 6 7]) = nan;
patch([1 1 2 2], 1*[1 1 1 1]+.3*[-1 1 1 -1], 0.5*[1 1 1],'EdgeColor','k');
plot([1 2],      2*[1 1], 'Color',clr(6,:), 'LineWidth',2);
% plot([1 2],      3*[1 1], 'Color',clr(6,:), 'LineWidth',.2);
plot(dashed_line,      3*ones(length(dashed_line)), 'Color',clr(6,:), 'LineWidth',2);
text(2.6, 1, 'Data','FontSize',7,'HorizontalAlignment','left');
text(2.6, 2, 'Model (Scheme 6)','FontSize',7,'HorizontalAlignment','left');
text(2.6, 2.5, {'Model (Scheme 6v -';'            variable coverage)'},'FontSize',7,'HorizontalAlignment','left','VerticalAlignment','top');
hax=gca;
hax.Visible='off';
hax.YDir='reverse';

%% panel H
plot_panel_different_scaling_delta(panel_H);

%% print/save the figure
fig_name_out = fullfile(res_dir, sprintf('%s_dt_%s', ...
    fig_name_str, ...
    strrep(num2str(panels_AE_dt),'.','_')));
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');


%%
function plot_panel_different_scaling_delta(panel_axes)

load('L:\Yonatan_theory\20200630__new_simulations_data+script\SummaryDecoderResults_EnvironmentSizeScaling.mat')
nalp = length(alp) ;
L = Lv*ds/100 ; % environment size variable in meters
jN = 2 ; % chooses N = 50

% color variable for S6 advantage (blue/red gradient)
nclr = 50 ;  % number of colors 
cmin = -1.5 ;
cbi  = -1   ;
c0   = 0 ;
cri  = 1 ;
cmax = 1.5 ;
cx   = linspace(cmin,cmax,nclr) ;
clrf = zeros(nclr,3) ;
iclr1 = find(cx<=cbi) ;
iclr2 = intersect(find(cx<=c0),find(cx>=cbi)) ;
iclr3 = intersect(find(cx<=cri),find(cx>=c0)) ;
iclr4 = find(cx>=cri) ;

clrf(iclr1,3) = linspace(0.5,1,length(iclr1)) ;
clrf(iclr2,1) = linspace(0  ,1,length(iclr2)) ;
clrf(iclr2,2) = linspace(0  ,1,length(iclr2)) ;
clrf(iclr2,3) = 1 ; 

clrf(iclr3,1) = 1 ; 
clrf(iclr3,2) = linspace(1,0  ,length(iclr3)) ;
clrf(iclr3,3) = linspace(1,0  ,length(iclr3)) ;
clrf(iclr4,1) = linspace(1,0.5,length(iclr4)) ;

Lint = logspace(log10(20),3,100) ; % interpolation of environment size variable
alpint = 0.1:0.01:0.9 ;            % interpolation of scaling exponent variable
[Lm,alpm] = meshgrid(L,alp) ;
[Lmint,alpmint] = meshgrid(Lint,alpint) ;
merr_ML_S1i = interp2(Lm,alpm,ones(nalp,1)*merr_ML_S1(jN,:),Lmint,alpmint) ;
merr_ML_S2i = interp2(Lm,alpm,squeeze(merr_ML_S2(jN,:,:))',Lmint,alpmint) ;
merr_ML_S3i = interp2(Lm,alpm,squeeze(merr_ML_S3(jN,:,:))',Lmint,alpmint) ;
merr_ML_S4i = interp2(Lm,alpm,squeeze(merr_ML_S4(jN,:,:))',Lmint,alpmint) ;
merr_ML_S5i = interp2(Lm,alpm,squeeze(merr_ML_S5(jN,:,:))',Lmint,alpmint) ;
merr_ML_S6i = interp2(Lm,alpm,squeeze(merr_ML_S6(jN,:,:))',Lmint,alpmint) ;

plerr_ML_S1i = interp2(Lm,alpm,ones(nalp,1)*plerr_ML_S1(jN,:),Lmint,alpmint) ;
plerr_ML_S2i = interp2(Lm,alpm,squeeze(plerr_ML_S2(jN,:,:))',Lmint,alpmint) ;
plerr_ML_S3i = interp2(Lm,alpm,squeeze(plerr_ML_S3(jN,:,:))',Lmint,alpmint) ;
plerr_ML_S4i = interp2(Lm,alpm,squeeze(plerr_ML_S4(jN,:,:))',Lmint,alpmint) ;
plerr_ML_S5i = interp2(Lm,alpm,squeeze(plerr_ML_S5(jN,:,:))',Lmint,alpmint) ;
plerr_ML_S6i = interp2(Lm,alpm,squeeze(plerr_ML_S6(jN,:,:))',Lmint,alpmint) ;

min_merr_ML_S12345i = min(merr_ML_S1i,merr_ML_S2i) ;
min_merr_ML_S12345i = min(min_merr_ML_S12345i,merr_ML_S3i) ;
min_merr_ML_S12345i = min(min_merr_ML_S12345i,merr_ML_S4i) ;
min_merr_ML_S12345i = min(min_merr_ML_S12345i,merr_ML_S5i) ;

min_plerr_ML_S12345i = min(plerr_ML_S1i,plerr_ML_S2i) ;
min_plerr_ML_S12345i = min(min_plerr_ML_S12345i,plerr_ML_S3i) ;
min_plerr_ML_S12345i = min(min_plerr_ML_S12345i,plerr_ML_S4i) ;
min_plerr_ML_S12345i = min(min_plerr_ML_S12345i,plerr_ML_S5i) ;

axes(panel_axes);
cla
imagesc(Lint,alpint,log10(plerr_ML_S6i./min_plerr_ML_S12345i),[cmin cmax]);
set(gca,'xscale','log') ;
axis square xy;
hcb = colorbar ;
hcb.Ruler.TickLabelGapOffset = 0.5;
colormap(gca,clrf)
hl = yline(0.3);
hl.Color = 'r';
hl.LineStyle = '--';
hl.LineWidth = 1.2;
hax = gca;
hax.XScale = 'log';
hax.XLim = [20 1000];
% hax.YLim = [20 1000];
hax.XTick = [100 1000];
hax.YTick = [0.1:0.1:0.9];
% hax.XTickLabel = {};
hax.YRuler.TickLabelGapOffset = -0.1;
hax.XRuler.TickLabelGapOffset = -1;
hax.XRuler.TickLength = [0.03 0.03];
hax.YRuler.TickLength = [0.024 0.03];
xlabel('Environment size (m)', 'Units','normalized','Position',[0.5 -0.16]);
ylabel('Scaling factor', 'Units','normalized','Position',[ -0.25 0.5]);
text(-0.35,1.15, 'H', 'Units','normalized','FontWeight','bold');
text(1.52,0.5, {'Log ratio prob. error';'(scheme 6 / other schemes)'}, 'Units','normalized','Rotation',-90,'FontSize',7,'HorizontalAlignment','center');

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
end


%%
