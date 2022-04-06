%% PhD thesis figure 5.1 - discussion: phase-precession with different place-fields

%%
close all
clear 
clc

%% define output files
res_dir = 'E:\Tamir\PhD\Thesis\resources\ch_5_discussion';
mkdir(res_dir)
fig_name_str = 'Fig_5_1_phase_precession_multiscale';
fig_caption_str = ' ';
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
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');

% create panels
panels(1) = axes('position', [3 19   10 2]);
panels(2) = axes('position', [3 15.5 10 3]);

%
nCells = 3;
nSpikes = 1000;
fields_locs = 50.*[1 1 1];
% fields_sizes = sqrt([30 10 1]);
fields_sizes = [6 2 0.2];
cvr = -0.9;
clrs = [1 0 0; 0 1 0; 0 0 1];

rng('default');
axes(panels(2))
hold on
for ii_cell=1:nCells
    mu = [0 0];
    Sigma = [1 cvr; cvr 1];
    R = mvnrnd(mu,Sigma,nSpikes);
    R(:,1) = R(:,1).*fields_sizes(ii_cell) + + fields_locs(ii_cell);
%     R(:,1) = rescale(R(:,1),0,fields_sizes(ii_cell));
%     R(:,1) = R(:,1) + fields_locs(ii_cell)-fields_sizes(ii_cell)/2;
    R(:,2) = rescale(R(:,2),0,360);
    plot(R(:,1),[R(:,2) R(:,2)+360],'.','Color',clrs(ii_cell,:));
    hax=gca;
%     hax.XTick=[];
    hax.YLim = [0 720];
    hax.YTick = [0 360 720];
    ylabel('Phase ({\circ})')
    xlabel('Position (m)')
end
box on 

axes(panels(1))
hold on
x = linspace(30,70,1000)';
y = normpdf(x, fields_locs, fields_sizes);
y = normalize(y,1,"range");
y(y<0.001)=nan;
for ii_cell=1:nCells
%     h=area(x,y(:,ii_cell));
%     h.FaceAlpha = 0.8;
%     h.FaceColor = clrs(ii_cell,:);
    h=plot(x,y(:,ii_cell));
    h.Color = clrs(ii_cell,:);
    h.LineWidth = 2;
    hax=gca;
    hax.XTick=[];
    hax.YTick=[];
    ylabel('Firing rate')
end
legend("Cell "+[1:nCells],'Box','off')

linkaxes(panels, 'x')
xlim([30 70])

%% save fig(s)
fig_name_out = fullfile(res_dir, [fig_name_str ]);
print(fig, fig_name_out, '-dpdf', '-cmyk', '-painters');

disp('figure was successfully saved to pdf/tiff/fig formats');
diary off

%%
