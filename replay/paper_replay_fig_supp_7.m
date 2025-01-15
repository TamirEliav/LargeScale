%% Replay - Fig supp 9 - two bats same map criteria
clear 
clc
close all

%% data options 
params_opt = 11; % decoding opt 

%% plotting options

%% graphics params

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Figure_S9';
fig_caption_str = '';
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
% figure_size_cm = [21.0 29.7]; % ~A4
figure_size_cm = [21.6 27.9]; % ~US letter
fig = figure;
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
% set(gcf, 'color', 'none');
set(groot, 'defaultAxesColor','None')
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');

% create panels
panels{1}(1) = axes('position', [3 20 4 3.5]);
panels{2}(1) = axes('position', [3 15 3 3]);
panels{2}(2) = axes('position', [8 15 3 3]);

%% replay directionality bias - plot trend over exposure to enviroenment (resampled version)
load('E:\Tamir\work\PROJECTS\LargeScale\paper_replay\figures\replay_directionality.mat')
% clrs = [epoch_type_clrs,'k'];
% surprise_clipping = 50;
for ii_epoch_type = 3%1:length(epoch_types)
    axes(panels{1}(1));
    cla reset
    hold on
    yline(0,'-k');
%     c = clrs{ii_epoch_type};
    x = T.session_num_from_exposure;
    y = directionality_binom_surprise2(ii_epoch_type,:);
%     y(~ismember(T.batNum,novelty_exposure_bats)) = nan;
    y(nSeqs(ii_epoch_type,:)<minSeqsThr) = nan;
%     y(y>surprise_clipping) = surprise_clipping;
%    y(y<-surprise_clipping) = -surprise_clipping;
%    ylimits = [-1 1]*surprise_clipping;
   ylimits = [-1 1]*10;
    for ii_bat = 1:length(bats)
        bat_num = bats(ii_bat);
        c = bats_clr_map(bat_num);
        IX = T.bat_num==bat_num;
        xx = x(IX);
        yy = y(IX);
        invalid = isnan(yy);
        xx(invalid) = [];
        yy(invalid) = [];
        plot(xx,yy,'o-','color',c,'DisplayName',"bat "+bat_num,'MarkerFaceColor',c,'MarkerSize',2);
    end

%     nPointsSmooth = 3;
%     k = (nPointsSmooth-1)/2;
%     xi = [-1 1].*k + [(1-k):(max(x)+k)]';
%     xx=[];
%     yy=[];
%     for ii_xi = 1:size(xi,1)
%         TF = x>=xi(ii_xi,1) & x<=xi(ii_xi,2);
%         xx(ii_xi) = nanmedian(x(TF));
%         yy(ii_xi) = nanmedian(y(TF));
%     end
%     plot(xx,yy,'-','Color',c);
    ylim(ylimits)
    ylabel({'Replay directionality index';'(resampled data)'})
%     text(0.5,0.8,epoch_types{ii_epoch_type},'Units','normalized')
end
hax=gca;
text(.5, 2.8, {'Replay directionality';'all sessions'},'FontSize',9,'HorizontalAlignment','center','Units','normalized')
% legend
% x = [0.3 2]+32;
% y = linspace(0.35,0,3)+0.6;
% plot(x,y(1)*[1 1],'Color',clrs{1},'LineWidth',1.5,'Clipping','off')
% plot(x,y(2)*[1 1],'Color',clrs{2},'LineWidth',1.5,'Clipping','off')
% plot(x,y(3)*[1 1],'Color',clrs{3},'LineWidth',1.5,'Clipping','off')
% x = x(end)+1;
% text(x,y(1), "Sleep", 'FontSize',7)
% text(x,y(2), "Awake", 'FontSize',7)
% text(x,y(3), "Combined", 'FontSize',7)
xlabel('Session no.','Units','normalized','Position',[0.5 -0.05]);
ylim([-1 1]*7)
yticks([-5 0 5])
xticks([1 40])
hax=gca;
hax.XRuler.TickLength(1) = 0.02;
hax.YRuler.TickLength(1) = 0.035;
hax.XRuler.TickLabelGapOffset = -1;
hax.YRuler.TickLabelGapOffset = 0;
% text(5.5,-0.75,"\leftarrow"+"session #"+example_session_num_form_exposure,'FontSize',8)
h=annotation('textarrow');
h.Parent=hax;
h.X = [6 5];
h.Y = [1 1]*3.5;
h.String = {'  SAS=4'};
h.FontSize = 7; h.HeadLength = 4; h.HeadWidth = 3; h.HeadStyle = 'cback2';
h.HorizontalAlignment = 'left';
h.VerticalAlignment = 'middle';
h=annotation('textarrow');
h.Parent=hax;
h.X = 16 + [0 0.1];
h.Y = -5.65 + [0 0];
h.String = {'SAS=8   '};
h.FontSize = 7; h.HeadLength = 4; h.HeadWidth = 3; h.HeadStyle = 'cback2';
h.HorizontalAlignment = 'right';
h.VerticalAlignment = 'middle';
h=annotation('textarrow');
h.Parent=hax;
h.X = 25 + [0 1];
h.Y = -4 + [0 0];
h.String = {'SAS=27 '};
h.FontSize = 7; h.HeadLength = 4; h.HeadWidth = 3; h.HeadStyle = 'cback2';
h.HorizontalAlignment = 'right';
h.VerticalAlignment = 'middle';
h=annotation('textarrow');
h.Parent=hax;
h.X = 30 + [0 -0.1];
h.Y = -5.8 + [0 0];
h.String = {'   SAS=19'};
h.FontSize = 7; h.HeadLength = 4; h.HeadWidth = 3; h.HeadStyle = 'cback2';
h.HorizontalAlignment = 'left';
h.VerticalAlignment = 'middle';

%% Main scatter plot
axes(panels{2}(1))
hold on
% data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115.mat';
data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115 & same map.mat';
% data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115 & same map & forward.mat';
% data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115 & same map & reverse.mat';
% data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115 & replay distance_gt_5.mat';
% data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115 & replay distance_gt_10.mat';
% data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115 & forward.mat';
% data_filename = 'F:\sequences\figures\replay_vs_crossover\replay_vs_crossover_opt_11_replay position between 25-115 & reverse.mat';
data = load(data_filename);
scatter(data.x(data.TF),data.y(data.TF),5,'k','filled');
axis equal
xlim([0 135])
ylim([0 135])
xticks(linspace(0,135,4))
yticks(linspace(0,135,4))
xlabel('Previous cross-over position (m)', 'Units','normalized', 'Position',[0.5 -0.16]);
ylabel('Replay position (m)', 'Units','normalized', 'Position',[-0.2 .5]);
% text(0.8,0.3,"n = "+ data.stats.n,'Units','normalized','FontSize',9);
% text(0.6,0.2,"r = "+ sprintf('%.2g',data.stats.Pearson.r), 'Units','normalized','FontSize',7);
% text(0.6,0.1,"P = "+ sprintf('%.2g',data.stats.Pearson.p), 'Units','normalized','FontSize',7);
text(.3,1.2,"{\rho} = "+ sprintf('%.2f',data.stats.Spearman.r), 'Units','normalized','FontSize',7);
text(.3,1.1,"P = "+ sprintf('%.2g',data.stats.Spearman.p), 'Units','normalized','FontSize',7);
% text(.3,1.3,"n = "+ sprintf('%d',sum(data.TF)), 'Units','normalized','FontSize',7);
% text(0,-.4,data.msg_str, 'Units','normalized','FontSize',10);
h=refline(1,0);
h.Color = .8.*[1 1 1];
hax=gca;
hax.XRuler.TickLength(1) = 0.035;
hax.YRuler.TickLength(1) = 0.024;
hax.XRuler.TickLabelGapOffset = -.5;
hax.YRuler.TickLabelGapOffset = 1;

%% Scatter plot (control - next crossover)
axes(panels{2}(2))
cla
hold on
X = [data.seqs_all.next_co_pos];
Y = data.y;
X = X(data.TF)';
Y = Y(data.TF)';
[stats.Pearson.r, stats.Pearson.p] = corr(X,Y,'type','Pearson',rows='pairwise',tail='right');
[stats.Spearman.r, stats.Spearman.p] = corr(X,Y,'type','Spearman',rows='pairwise',tail='right');
scatter(X,Y,5,'k','filled');
axis equal
xlim([0 135])
ylim([0 135])
xticks(linspace(0,135,4))
yticks(linspace(0,135,4))
xlabel('Next cross-over position (m)', 'Units','normalized', 'Position',[0.5 -0.16]);
ylabel('Replay position (m)', 'Units','normalized', 'Position',[-0.2 .5]);
% text(0.8,0.3,"n = "+ data.stats.n,'Units','normalized','FontSize',9);
% text(0.8,0.2,"r = "+ sprintf('%.2g',stats.Pearson.r), 'Units','normalized','FontSize',7);
% text(0.8,0.1,"P = "+ sprintf('%.2g',stats.Pearson.p), 'Units','normalized','FontSize',7);
text(.3,1.2,"{\rho} = "+ sprintf('%.2f',stats.Spearman.r), 'Units','normalized','FontSize',7);
text(.3,1.1,"P = "+ sprintf('%.2f',stats.Spearman.p), 'Units','normalized','FontSize',7);
h=refline(1,0);
h.Color = .8.*[1 1 1];
hax=gca;
hax.XRuler.TickLength(1) = 0.035;
hax.YRuler.TickLength(1) = 0.024;
hax.XRuler.TickLabelGapOffset = -.5;
hax.YRuler.TickLabelGapOffset = 1;


%% add panel letters
font_size = 11;
axes(panels{1}(1))
text(-0.35,1.12, 'A', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{2}(1))
text(-0.35,1.12, 'B', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{2}(2))
text(-0.35,1.12, 'C', 'Units','normalized','FontWeight','bold','FontSize',font_size);

%%
fig_name = sprintf('%s',fig_name_str);
% fig_name = sprintf('%s_panel_B_opt_%d',fig_name,behavior_ex_opt);
% fig_name = sprintf('%s_max_tdiff_%ds',fig_name,timediff_max);
[~,data_str,~] = fileparts(data_filename);
fig_name = sprintf('%s__%s',fig_name,data_str);
file_out = fullfile(res_dir, fig_name);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
disp('figure saved!')

%%
