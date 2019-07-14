%% Large Scale - fig 1 - Behavioral and neural recordings from bats fliying over large spatial scales.

%%
clear 
clc

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_1';
fig_caption_str = 'Behavioral and neural recordings from bats fliying over large spatial scales';
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
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none');

% create panels
panel_B_size = [9.5 1];
panel_A    = axes('position', [ 1 23.2  2 2]);
panel_B(1) = axes('position', [ 4 24  panel_B_size]);
panel_B(2) = axes('position', [ 4 23  panel_B_size]);
panel_B(3) = axes('position', [ 4 22  panel_B_size]);
panel_B(4) = axes('position', [ 4 21  panel_B_size]);
panel_C    = axes('position', [15 21  4 4]);
panel_D    = axes('position', [ 1 21  2 2]);
panel_E    = axes('position', [ 1 16 6 4]);
panel_F    = axes('position', [ 7 16 6 4]);
panel_G    = axes('position', [13 16 6 4]);
panel_H    = axes('position', [ 1 13 2 2]);
panel_I(1) = axes('position', [ 4 13 8 2]);
panel_I(2) = axes('position', [13 13 6 2]);
panel_J    = axes('position', [ 1  9 6 2.5]);
panel_K    = axes('position', [ 8  9 3 2.5]);
panel_L    = axes('position', [12  9 3 2.5]);
panel_M    = axes('position', [16  9 3 2.5]);

%
prm = PARAMS_GetAll();

%% ---------------------- spikes trace ------------------------------------
% define option for data
exp_ID = 'b0148_d170625';
exp = exp_load_data(exp_ID,'details','path');
TT = 4;
trace_duration = 0.35;
trace_ts_opt = 6;
trace_ts_list = [
    70270204396,...
    70132623245,...
    74986863162,...
    74986792310,...
    72188230643,...
    71432461878
    ];
traces_ts = trace_ts_list(trace_ts_opt) + [0 trace_duration]*1e6;
% load channels data
signal = [];
for ii_ch = 1:4
    
    file_in = fullfile(exp.path.spikes_raw, sprintf('spikes_%s_TT%d_ch%d.ncs',exp_ID,TT,ii_ch) );
    [signal(ii_ch,:), ts, fs] = Nlx_csc_read(file_in, traces_ts);
end

% plot the data
for ch = 1:4
    axes(panel_B(ch));
    plot(ts,signal(ch,:),'k','LineWidth',0.01);
    box off
    set(gca,'Visible','Off')
    xlim(ts([1 end]))
end
linkaxes(panel_B,'xy')

% add time/voltage scales
axes(panel_B(4));
xlimits = get(gca,'xlim');
ylimits = get(gca,'ylim');
scale_ms = 20;
scale_uVolt = 500;
scale_line_width = 1;
xa = xlimits(1) + [0 scale_ms*1e3];
ya = ylimits(1) + [0 scale_uVolt ];
[xaf,yaf] = ds2nfu(xa,ya);
annotation('line',xaf([1 2])+0.01,yaf([1 1])-0.02,'Linewidth',scale_line_width); % time
annotation('line',xaf([1 1])+0.01,yaf([1 2])-0.02,'Linewidth',scale_line_width)  % voltage
% h=annotation('line','Units','centimeters');
% h.X = xaf([1 2])+5;
% h.Y = yaf([1 1]);
% h

% add panel letter
axes(panel_B(1));
text(-0.05,1.1, 'B', 'Units','normalized','FontWeight','bold');

%%
clusters_colors = [
166,206,227
31,120,180
178,223,138
51,160,44
251,154,153
227,26,28
253,191,111
255,127,0
202,178,214
106,61,154
255,255,153
177,89,40];
clusters_colors = clusters_colors ./255;

%% ---------------------- spikes clusters ---------------------------------
if 1
NTT_file = 'L:\Analysis\pre_proc\0148\20170625\spikes_sorting\spikes_b0148_d170625_TT4.NTT';
limits_ts = [70019176845 70816687098];
% limits_ts = [];
if isempty(limits_ts)
    [Timestamps, CellNumbers, Samples, Header] = Nlx2MatSpike(NTT_file, [1 0 1 0 1], 1, 1, [] );
else
    [Timestamps, CellNumbers, Samples, Header] = Nlx2MatSpike(NTT_file, [1 0 1 0 1], 1, 4, limits_ts);
end
ADBitVolts = sscanf(Header{contains(Header,'ADBitVolts')},'-ADBitVolts %f %f %f %f');
Samples = Samples .* ADBitVolts' .* 1e6; % convert bits to uVolts
spikes_size = squeeze(range(Samples,1));
cells = unique(CellNumbers);
color_list = linspecer( length(cells) ,'qualitative');
rng(3)
% rng(12)
color_list = color_list(randperm(size(color_list,1)),:);
color_list(1,:) = 0.5.*[1 1  1];
% plot
axes(panel_C); hold on
ch2plot = [2 3 4];
X = spikes_size(ch2plot(1),:);
Y = spikes_size(ch2plot(2),:);
Z = spikes_size(ch2plot(3),:);
for ii_cell = 1:length(cells)
    cellNum = cells(ii_cell);
    if cellNum==0
        continue
    end
    cell_spikes_IX = find(CellNumbers == cellNum);
    x = X(cell_spikes_IX);
    y = Y(cell_spikes_IX);
    z = Z(cell_spikes_IX);
    plot3(x,y,z,'.','Color',color_list(ii_cell,:),'MarkerSize',1);
end
xlim( [0 max(X(CellNumbers ~= 0))] )
ylim( [0 max(Y(CellNumbers ~= 0))] )
zlim( [0 max(Z(CellNumbers ~= 0))] )
volt_ticks = 0:200:1000;
set(gca,'xtick',volt_ticks,'ytick',volt_ticks,'ztick',volt_ticks)
view([140 39])
xlabel(['ch' num2str(ch2plot(1)) ' ({\mu}V)'])
ylabel(['ch' num2str(ch2plot(2)) ' ({\mu}V)'])
zlabel(['ch' num2str(ch2plot(3)) ' ({\mu}V)'])
AZ_EL(1,:) = [140 39];
AZ_EL(2,:) = [488 20];
AZ_EL(3,:) = [511 -8];
viewing_option = 3;
view(AZ_EL(viewing_option,:))
end
text(-0.05,1.1, 'C', 'Units','normalized','FontWeight','bold');

%% ----------- panel D - logger reception------------------
axes(panel_D)
x = 0:100:1000;
x = x./1e3;
y = 100.*ones(size(x));
y(end-3:end) = [90 70 30 10];
plot(x,y,'.-','LineWidth',2,'MarkerSize',14);
box off
ylim([0 110]);
ha=gca;
ha.XLim = [0 1];
ha.YLim = [0 110];
ha.XTick = [0 0.5 1];
ha.YTick = [0 100];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;
xlabel('Distance (km)',  'Units','normalized','Position',[0.5  -0.15])
ylabel('Signal (%)',     'Units','normalized','Position',[-0.1 0.45])
text(-0.2,1.1, 'D', 'Units','normalized','FontWeight','bold');

%% show images
logger_image_filename = 'D:\Tamir\PROJECTS\Neurologger\miniBat\pics\minibat_good_pic.jpg';
% tunnel_view_image_file = 'L:\Videos_Photos\tunnel_area_various\20170111_170850_downsampled.jpg';
tunnel_view_image_file = 'L:\Videos_Photos\TAZOT_HAMAMA\taza3.jpg';
bespoon_image_file = 'L:\Videos_Photos\2016_ICN_poster\VirtualBox_ubuntu_21_03_2016_11_16_53 - Copy.png';
tunnel_section_file = 'L:\resources\tunnel_section.png';

axes(panel_A);
image = imread(logger_image_filename);
imshow(image);
text(-0.1,1, 'A', 'Units','normalized','FontWeight','bold');

axes(panel_E);
image = imread(tunnel_view_image_file);
imshow(image);
text(-0.1,1, 'E', 'Units','normalized','FontWeight','bold');

axes(panel_F)
image = imread(tunnel_section_file);
imshow(image);
text(-0.1,1, 'F', 'Units','normalized','FontWeight','bold');

axes(panel_G);
image = imread(bespoon_image_file);
imshow(image)
text(-0.1,1, 'G', 'Units','normalized','FontWeight','bold');

%% bespoon localization precision
axes(panel_H);
hold on
bespoon_loc_precision = load('L:\BeSpoon\testing\test_20180530__YOM_KEF_200m_static+dynamic+discretization\30-05-2018__calib_test_dynamic+non-jitter_jitter_with_kalman\data\outside_perpendicular_error.mat');
err = bespoon_loc_precision.perpendicular_error;
err = err.*100; % convert to cm
h=histogram(err);
h.NumBins = 45;
h.Normalization = 'pdf';
[muHat,sigmaHat] = normfit(err);
x = linspace(-50,50,100);
y = normpdf(x,muHat,sigmaHat);
plot(x,y,'r','LineWidth',2);
text(0.65,0.85, sprintf('\x03C3=%.1f', sigmaHat), 'Units','normalized','FontSize',8);
xlabel({'Positioning';'error (cm)'});
ylabel('Prob.');
ha = gca;
ha.XLim = [-40 40];
ha.XTick = [-40:20:40];
ha.YTick = [];
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;
text(-0.1,1.1, 'H', 'Units','normalized','FontWeight','bold');

%% behavioral trajectory is 1D (small y deviations) - example
axes(panel_I(1));
cla
hold on
panel_I_example_option = 1;
panel_I_example_list = {
'b0034_d180413';
'b2289_d180514';
'b2289_d180515';
'b2289_d180518';
'b9861_d180705';
'b9861_d180709';
};
exp = exp_load_data(panel_I_example_list{panel_I_example_option},'flight');
for ii_dir = [1 2] 
    c = prm.graphics.colors.flight_directions{ii_dir};
    ydev = exp.flight.pos_y_std(ii_dir);
    x = ydev.xy(:,1);
    y = ydev.xy(:,2);
    ymean = interp1(ydev.bin_centers, ydev.ymean, x);
    y = y-ymean;
    plot(x, y, '.', 'Color',c, 'MarkerSize',.0001);
%     h=shadedErrorBar(ydev.bin_centers, ydev.ymean, ydev.ystd, 'lineprops',{'Color',c});
%     h.patch.FaceAlpha = 0;
end
xlabel('Linearized Position (cm)','Units','normalized','Position',[0.5 -0.2]);
% ylabel('Y Position (cm)');
ylabel('Y deviation (cm)');
ha = gca;
ha.XLim = [0 200];
ha.YLim = [-1.5 1.5];
ha.XTick = [0:50:200];
ha.YTick = [-2:1:2];
ha.TickDir='out';
ha.TickLength = [0.01 0.01];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;
text(-0.1,1.1, 'I', 'Units','normalized','FontWeight','bold');

%% behavioral trajectory is 1D (small y deviations) - population
axes(panel_I(2));
cla
hold on
% use exps that we recorded cells (TODO: consider using 
cells_t = DS_get_cells_summary();
bats = [79,148,34,9861,2289];
cells_t(~ismember(cells_t.bat, bats ),:) = [];
cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells = [cells.details];
cells(~contains({cells.brain_area}, 'CA1')) = [];
cells(~ismember([cells.ClusterQuality], [2])) = [];
exp_list = unique({cells.exp_ID}'); % get the relevant exp list
exps = cellfun(@(x)(exp_load_data(x,'details','flight')), exp_list);
exps_details = [exps.details];
exps_flight = [exps.flight];
pos_y_dev_all = cat(1,exps_flight.pos_y_std);
speed_traj_all = cat(1,exps_flight.speed_traj);
ystd_median_all = arrayfun(@(x)(x.ystd_median), pos_y_dev_all);

%%
axes(panel_I(2));
cla
hold on
% arrange data
data = {};
for ii_bat = 1:length(bats)
    bat = bats(ii_bat);
    c = prm.graphics.colors.bats(bat);
    bat_exp_IX = ismember([exps_details.batNum],bat);
    bat_exp_ydev = pos_y_dev_all(bat_exp_IX,:);
    valid_pos_IX = get_data_in_ti(bat_exp_ydev(1).bin_centers, prm.fields.valid_speed_pos);
    ystd_bat_all = cat(1,bat_exp_ydev.ystd); % pooling directions together!
    ystd_bat_all = ystd_bat_all(:,valid_pos_IX); % take only valid positions
    ystd_bat_all = ystd_bat_all(:)';% pooling over all positions!
    ystd_bat_all = ystd_bat_all.*100; % convert to cm
    data{ii_bat} = ystd_bat_all;
end
% t2 = table2cell(t);
x = [data{:}];
grps = repelem(1:length(bats), cellfun(@length,data));

% plot!
cmap = prm.graphics.colors.bats;
c = arrayfun(@(x)(cmap(x)), bats,'UniformOutput',0);
% yvar_pop_plot = 'violin';
% yvar_pop_plot = 'boxplot';
yvar_pop_plot = 'median_std';
switch yvar_pop_plot 
    case 'violin'
        hv=violinplot(x,grps);
        [hv.BoxWidth] = disperse(repelem(0.02,length(hv)));
        [hv.EdgeColor ] = disperse(repelem({'none'},length(hv)));
        hs = [hv.ScatterPlot];
        [hs.SizeData] = disperse(repelem(5,length(hs)));
        hm = [hv.MedianPlot];
        [hm.SizeData] = disperse(repelem(10,length(hm)));
        ha.YLim = [0 90];
        ha.YTick = [0:30:90];
    case 'boxplot'
        boxplot(x,grps)
        box off
        ha.YLim = [0 90];
        ha.YTick = [0:30:90];
    case 'median_std'
        medians = cellfun(@median,data);
        stds    = cellfun(@std,data);
        errorbar(medians,stds,'-o','MarkerSize',3);
        ha.YLim = [0 40];
        ha.YTick = [0:20:40];
end
ha=gca;
ha.XTickLabel = bats;
ha.XLim = [0.5 length(bats)+.5];
ha.TickDir='out';
ha.TickLength = [0.01 0.01];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;

xlabel('Bats');
ylabel('Y variability (cm)');
% text(-0.15,1.1, 'I', 'Units','normalized','FontWeight','bold');

%% speed trajectory very constant along the flight - example
axes(panel_J); 
cla
hold on
panel_speed_option = 4;
panel_speed_options = {
'b0079_d160913';
'b0079_d160914';
'b0079_d160915'; % many laps with no landing in the middle but different speed between directions
'b2289_d180615'; % looks excellent but not many laps!
'b0079_d160928';
};
exp_ID = panel_speed_options{panel_speed_option};
exp = exp_load_data(exp_ID,'flight');
FE=exp.flight.FE;
x = [FE.pos];
y = [FE.vel];
y = abs(y);
ylimits = [0 9];
area([0 prm.fields.valid_speed_pos(1)],     ylimits([2 2]),'FaceColor',0.9*[1 1 1],'EdgeColor','none')
area([  prm.fields.valid_speed_pos(2) 200], ylimits([2 2]),'FaceColor',0.9*[1 1 1],'EdgeColor','none')
plot(x,y,'.','Color', 'k','MarkerSize',1);
% plot(prm.fields.valid_speed_pos([1 1]), ylimits,'--m','LineWidth',2)
% plot(prm.fields.valid_speed_pos([2 2]), ylimits,'--m','LineWidth',2)
set(gca,'xtick',0:50:200,'ytick',[-10 0 8],'xlim',[0 200])
ylim(ylimits);
xlabel('Position (m)');
ylabel('Speed (m/s)');
text(-0.1,1.1, 'J', 'Units','normalized','FontWeight','bold');

%% speed trajectory very constant along the flight - population
% KLM - TODO: use only valid flights + signif place cell in that day!
axes(panel_K);
hold on
dir_colors = prm.graphics.colors.flight_directions;
for ii_dir = 1:2
    cv = [speed_traj_all(:,ii_dir).speed_cv];
%     cv = [cv.raw_high_speed];
    cv = [cv.across_pos_high_speed];
    h=histogram(cv);
    h.BinEdges = linspace(0,0.1,20);
    h.FaceColor = dir_colors{ii_dir};
    h.Normalization = 'Count';
    h.FaceAlpha = 0.5;
end
xlabel('CV of speed');
ylabel('Counts');
text(-0.2,1.1, 'K', 'Units','normalized','FontWeight','bold');

%% number of laps
axes(panel_L);
hold on
dir_colors = prm.graphics.colors.flight_directions;
directions = [-1 1];
for ii_dir = 1:2
    direction = directions(ii_dir);
    FE_dir = cellfun(@(FE)(FE([FE.direction]==direction)), {exps_flight.FE},'UniformOutput',0);
    nFlights = cellfun(@(FE)(sum([FE.distance]>prm.flight.full_min_distance)),FE_dir);
    h=histogram(nFlights);
    h.NumBins = 10;
    h.FaceColor = dir_colors{ii_dir};
    h.Normalization = 'Count';
    h.FaceAlpha = 0.5;
end
xlabel('No. of laps');
ylabel('Counts');
text(-0.2,1.1, 'L', 'Units','normalized','FontWeight','bold');

%% Total distance 
axes(panel_M);
hold on
if 1
    %% (per direction)
    dir_colors = prm.graphics.colors.flight_directions;
    directions = [-1 1];
    for ii_dir = 1:2
        direction = directions(ii_dir);
        FE_dir = cellfun(@(FE)(FE([FE.direction]==direction)), {exps_flight.FE},'UniformOutput',0);
        total_distance = cellfun(@(FE)(sum([FE.distance])),FE_dir);
        total_distance = total_distance .*1e-3; % to km
        h=histogram(total_distance);
        h.NumBins = 15;
        h.FaceColor = dir_colors{ii_dir};
        h.Normalization = 'Count';
        h.FaceAlpha = 0.5;
    end
else
    %% Total distance (pooling directions)
    dir_colors = prm.graphics.colors.flight_directions;
    total_distance = cellfun(@(FE)(sum([FE.distance])),{exps_flight.FE});
    total_distance = total_distance .*1e-3; % to km
    h=histogram(total_distance);
    h.NumBins = 13;
    h.FaceColor = dir_colors{ii_dir};
    h.Normalization = 'Count';
    h.FaceColor = 0.5*[1 1 1];
end
xlabel('Total Distance (km)');
ylabel('Counts');
text(-0.2,1.1, 'M', 'Units','normalized','FontWeight','bold');

%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');








%%
