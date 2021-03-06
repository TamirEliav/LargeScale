%% Fig. 8 - Dynamics panel - perturbation analysis options suammry
%%
clear
clc
rng(0);

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'Fig_8_Dynamics_panel_options';
fig_caption_str = ' ';
log_name_str = [fig_name_str '_log_file' '.txt'];
log_name_str = strrep(log_name_str , ':', '-');
log_name_str = strrep(log_name_str , ' ', '_');
log_name_out = fullfile(res_dir, log_name_str);

%% data folder
dir_IN = 'L:\Shir\20201001__dynamics_p1_p2';

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
panel_size = [3 3];
panel_PRC(1)       = axes('position', [ 1.5  17  panel_size]);
panel_PRC(2)       = axes('position', [ 5.5  17  panel_size]);
panel_PRC(3)       = axes('position', [ 9.5  17  panel_size]);
panel_PRC(4)       = axes('position', [ 13.5  17  panel_size]);
panel_PRC(5)       = axes('position', [ 17.5  17  panel_size]);

n_field_max = 10;

%% load data
load(fullfile(dir_IN,'Data_Perturbation_FieldSegmentNumbers.mat'));
% load('D:\LargeScale_master\Code\dynamic_analysis\Data_Perturbation_FieldSegmentNumbers.mat');
% nFieldD_0 : number of fields in original (entire session) maps, real data
% nFieldD_sSeg : number of appeared/disappeared segments (1st, 2nd column), real data

p1_square_data = zeros(n_field_max,1);
p2_data = zeros(n_field_max,1);
n_samples = 10000; % number of random curves for error bars computation
p1_square_data_rnd = zeros(n_field_max, n_samples);
p2_data_rnd = zeros(n_field_max, n_samples);
for k = 1:n_field_max
    ik = find(nFieldD_0==k) ;
    p1 = mean(mean(nFieldD_sSeg(ik,:)==1)) ;
    p2 = mean(mean(nFieldD_sSeg(ik,:)==2)) ;

    p1_square_data(k) = p1^2;
    p2_data(k) = p2;

    p1_square_error_lim = 2*p1*sqrt(p1*(1-p1))/sqrt(2*length(ik));
    p2_error_lim = sqrt(p2*(1-p2))/sqrt(2*length(ik));
    p1_square_data_rnd(k,:) = p1^2 + ((rand(n_samples,1)-0.5)*2*p1_square_error_lim);
    p2_data_rnd(k,:) = p2 + ( ( rand(n_samples,1)-0.5 )*2*p2_error_lim );
end

%% panel_%PRC:
for ii_prc = 1:5
    perturbation_prc = ii_prc;
    load(fullfile(dir_IN,['CA3MECCA1_FFModel_' num2str(perturbation_prc) 'PC_Perturbation_FieldSegmentNumbers_2']));
%     load(['D:\LargeScale_master\Code\dynamic_analysis\CA3MECCA1_FFModel_' num2str(perturbation_prc) 'PC_Perturbation_FieldSegmentNumbers_2']);
    % nFieldM_0 : number of fields in original (unperturbed) maps, multi  field CA3
    % nFieldS_0 : number of fields in original (unperturbed) maps, single field CA3
    % nFieldP_0 : number of fields in original (unperturbed) maps, periodic MEC
    % nFieldM_sSeg : number of appeared/disappeared segments (1st, 2nd column), multi  field CA3
    % nFieldS_sSeg : number of appeared/disappeared segments (1st, 2nd column), single field CA3
    % nFieldP_sSeg : number of appeared/disappeared segments (1st, 2nd column),  periodic MEC

    p1_square_S = zeros(n_field_max,1);
    p2_S = zeros(n_field_max,1);
    for k = 1:n_field_max
        ik = find(nFieldS_0==k) ;
        p1 = mean(mean(nFieldS_sSeg(ik,:)==1)) ;
        p2 = mean(mean(nFieldS_sSeg(ik,:)==2)) ;

        p1_square_S(k) = p1^2;
        p2_S(k) = p2;
    end

    p1_square_M = zeros(n_field_max,1);
    p2_M = zeros(n_field_max,1);
    for k = 1:n_field_max
        ik = find(nFieldM_0==k) ;
        p1 = mean(mean(nFieldM_sSeg(ik,:)==1)) ;
        p2 = mean(mean(nFieldM_sSeg(ik,:)==2)) ;

        p1_square_M(k) = p1^2;
        p2_M(k) = p2;
    end

    p1_square_P = zeros(n_field_max,1);
    p2_P = zeros(n_field_max,1);
    for k = 1:n_field_max
        ik = find(nFieldP_0==k) ;
        p1 = mean(mean(nFieldP_sSeg(ik,:)==1)) ;
        p2 = mean(mean(nFieldP_sSeg(ik,:)==2)) ;

        p1_square_P(k) = p1^2;
        p2_P(k) = p2;
    end

    corr_p1_square_S = corr(p1_square_data, p1_square_S,'type', 'spearman');
    corr_p1_square_M = corr(p1_square_data, p1_square_M,'type', 'spearman');
    corr_p1_square_P = corr(p1_square_data, p1_square_P,'type', 'spearman');
    corr_p2_S = corr(p2_data, p2_S,'type', 'spearman');
    corr_p2_M = corr(p2_data, p2_M,'type', 'spearman');
    corr_p2_P = corr(p2_data, p2_P,'type', 'spearman');


    corr_p1_square_S_rnd = corr(p1_square_data_rnd, p1_square_S,'type', 'spearman');
    corr_p1_square_M_rnd = corr(p1_square_data_rnd, p1_square_M,'type', 'spearman');
    corr_p1_square_P_rnd = corr(p1_square_data_rnd, p1_square_P,'type', 'spearman');
    corr_p2_S_rnd = corr(p2_data_rnd, p2_S,'type', 'spearman');
    corr_p2_M_rnd = corr(p2_data_rnd, p2_M,'type', 'spearman');
    corr_p2_P_rnd = corr(p2_data_rnd, p2_P,'type', 'spearman');

    std_corr_p1_square_S = std(corr_p1_square_S_rnd);
    std_corr_p1_square_M = std(corr_p1_square_M_rnd);
    std_corr_p1_square_P = std(corr_p1_square_P_rnd);
    std_corr_p2_S = std(corr_p2_S_rnd);
    std_corr_p2_M = std(corr_p2_M_rnd);
    std_corr_p2_P = std(corr_p2_P_rnd);

    %% plot
    axes(panel_PRC(ii_prc));
    cla;
    hold on;

    x = [1,2, 4,5, 7,8 ];
    y = [corr_p1_square_S, corr_p2_S, ...
        corr_p1_square_M, corr_p2_M,  ...
        corr_p1_square_P, corr_p2_P];
    err = [ std_corr_p1_square_S std_corr_p2_S ...
            std_corr_p1_square_M std_corr_p2_M ...
            std_corr_p1_square_P std_corr_p2_P];
    fig_8_dynamics_panel_data = struct();
    fig_8_dynamics_panel_data.x = x;
    fig_8_dynamics_panel_data.y = y;
    fig_8_dynamics_panel_data.err = err;
    fig_8_dynamics_panel_data.name{1} = "CA3_single_p1";
    fig_8_dynamics_panel_data.name{2} = "CA3_single_p2";
    fig_8_dynamics_panel_data.name{3} = "CA3_multi_p1";
    fig_8_dynamics_panel_data.name{4} = "CA3_multi_p2";
    fig_8_dynamics_panel_data.name{5} = "MEC_periodic_p1";
    fig_8_dynamics_panel_data.name{6} = "MEC_periodic_p2";
    filename = "fig_8_dynamics_panel_data_PRC_" + ii_prc;
    save(fullfile(res_dir,filename),'fig_8_dynamics_panel_data');
    
    hb = bar(x,y);
    hb.FaceColor = 'flat';
    hb.CData = [[0 0 0]; [0.6 0.6 0.6]; [0 0 0]; [0.6 0.6 0.6]; [0 0 0]; [0.6 0.6 0.6]];
    he=errorbar(x,y,err);

    he.Color = 'k';
    he.CapSize = 2;
    he.LineWidth = 1;
    he.LineStyle = 'none';
    he.Clipping='off';
    
    h=gca;
    h.XTick = x;
    h.XTickLabels = {'Single-field','CA3 model', 'Multi-field','CA3 model', 'MEC','model'};
    h.XTickLabelRotation = 45;
    ylimits = [-1,1];
    h.YLim = ylimits;
    h.YTick = [-1 0 1];
    h.TickLength = [0.01 0.01];
    h.XAxis.TickLength = [0 0];
    if ii_prc==1
        ylabel('Corr(data,model)','Units','normalized','Position',[-0.2 0.5]);
    end
    h.XRuler.TickLabelGapMultiplier = -0.3;
    h.YRuler.TickLabelGapMultiplier = 0.001;
    title([num2str(perturbation_prc) '% perturbation']);
    
    %% stats
    fprintf('p_1^2\n')
    [pval zval] = my_ztest(y(1), y(3), err(1), err(3));
    fprintf('Single-field CA3 vs. Multi-field CA3:  p=%.f z=%.f\n',pval,zval);
    [pval zval] = my_ztest(y(1), y(5), err(1), err(5));
    fprintf('Single-field CA3 vs. MEC:              p=%.f z=%.f\n',pval,zval);
    fprintf('p_2\n')
    [pval zval] = my_ztest(y(2), y(4), err(2), err(4));
    fprintf('Single-field CA3 vs. Multi-field CA3:  p=%.f z=%.f\n',pval,zval);
    [pval zval] = my_ztest(y(2), y(6), err(2), err(6));
    fprintf('Single-field CA3 vs. MEC:              p=%.f z=%.f\n',pval,zval);
    
end

%%
axes(panel_PRC(end));
hax = gca;
h_app = annotation('rectangle',[0 0 0 0],'FaceColor','k');
h_app.Parent = hax;
h_app.Position = [6 0.7 0.8 0.15];
h_dis = annotation('rectangle',[0 0 0 0],'FaceColor',[0.6 0.6 0.6]);
h_dis.Parent = hax;
h_dis.Position = [6 0.4 0.8 0.15];
text(7,0.75,'P_1^2', 'HorizontalAlignment','left','FontSize',8);
text(7,0.45,'P_2', 'HorizontalAlignment','left','FontSize',8);

%%
hh = annotation('rectangle',[0 0 0 0],'FaceColor','k');
hh.Parent = panel_PRC(1);
hh.Position = [7 0.5 0.1 0.1]

%% save figure
fig_name_out = fullfile(res_dir, sprintf('%s_bar_SD_version',fig_name_str));
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
disp('figure was successfully saved to pdf/tiff/fig formats');



%%
function [pval zval] = my_ztest(mu1, mu2, sigma1, sigma2)
    zval = (mu1-mu2) / sqrt(sigma1^2+sigma2^2);
    pval = 1-normcdf(zval,0,1);
end
