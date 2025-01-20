%% Replay - Fig R1 - replay directionality (contrast index version)
clear 
clc
close all

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Figure_R1';
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
panels{1}(1) = axes('position', [3 20 16 6]);

%% replay directionality bias - plot trend over exposure to enviroenment (resampled version)
load('E:\Tamir\work\PROJECTS\LargeScale\paper_replay\figures\replay_directionality.mat');
for ii_epoch_type = 3%1:length(epoch_types)
    axes(panels{1}(1));
    cla reset
    hold on
    yline(0,'-','Color',[1 1 1]*0.8,'LineWidth',0.2);
    x = T.session_num_from_exposure;
    y = directionality_contrast_index(ii_epoch_type,:);
    y(nSeqs(ii_epoch_type,:)<minSeqsThr) = nan;
   ylimits = [-1 1]*1.2;
    for ii_bat = 1:length(bats)
        bat_num = bats(ii_bat);
        c = bats_clr_map(bat_num);
        IX = T.bat_num==bat_num;
        xx = x(IX);
        yy = y(IX);
        invalid = isnan(yy);
        xx(invalid) = [];
        yy(invalid) = [];
        plot(xx,yy,'.-','color',c,'DisplayName',"bat "+bat_num,'MarkerFaceColor',c,'MarkerSize',1);
    end

    sz = nSeqs(ii_epoch_type,:);
            sz(isnan(y)) = nan;
%             sz = interp1([0 50],[0 1],sz,'linear','extrap');
            c = arrayfun(@(bat)bats_clr_map(bat),T.bat_num,'UniformOutput',0);
            c = cat(1,c{:});
            hbb=bubblechart(x,y,sz,c);
            bubblelim([min(sz) max(sz)])
            bubblesize([3 20])
            bubblelegend('No. of replays','location','eastoutside')
    ylim(ylimits)
    ylabel({'Replay directionality index';'(contrast index)'})
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
ylim([-1 1]*1.2)
yticks([-1 0 1])
xticks([1 40])
xlim([0 40.5])
hax=gca;
hax.XRuler.TickLength(1) = 0.01;
hax.YRuler.TickLength(1) = 0.01;
hax.XRuler.TickLabelGapOffset = -1;
hax.YRuler.TickLabelGapOffset = 0;
% text(5.5,-0.75,"\leftarrow"+"session #"+example_session_num_form_exposure,'FontSize',8)
h=annotation('textarrow');
h.Parent=hax;
h.X = 5.3 + [0 -1]*0.01;
h.Y = [1 1]*0.75;
% h.String = {'  SAS=4'};
h.FontSize = 7; h.HeadLength = 4; h.HeadWidth = 3; h.HeadStyle = 'cback2';
h.HorizontalAlignment = 'left';
h.VerticalAlignment = 'middle';
h=annotation('textarrow');
h.Parent=hax;
h.X = 16.3 + [0 1]*0.01;
h.Y = -1 + [0 0];
% h.String = {'SAS=8   '};
h.FontSize = 7; h.HeadLength = 4; h.HeadWidth = 3; h.HeadStyle = 'cback2';
h.HorizontalAlignment = 'right';
h.VerticalAlignment = 'middle';
h=annotation('textarrow');
h.Parent=hax;
h.X = 26.2 + [0 1]*0.01;
h.Y = -0.92 + [0 0];
% h.String = {'SAS=27 '};
h.FontSize = 7; h.HeadLength = 4; h.HeadWidth = 3; h.HeadStyle = 'cback2';
h.HorizontalAlignment = 'right';
h.VerticalAlignment = 'middle';
h=annotation('textarrow');
h.Parent=hax;
h.X = 29.9 + [0 -1]*0.01;
h.Y = -1 + [0 0];
% h.String = {'   SAS=19'};
h.FontSize = 7; h.HeadLength = 4; h.HeadWidth = 3; h.HeadStyle = 'cback2';
h.HorizontalAlignment = 'left';
h.VerticalAlignment = 'middle';

%%
fig_name = sprintf('%s',fig_name_str);
file_out = fullfile(res_dir, fig_name);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
disp('figure saved!')

%%
