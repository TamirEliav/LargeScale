%% Large Scale - fig 1 - Behavioral and neural recordings from bats fliying over large spatial scales.

%%
clear 
clc

%% choose data options
panel_C_opt = 2;
panel_G_opt = 1;
panel_I_opt = 6;
panel_J_opt = 10;

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'Fig_1';
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
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');

% create panels
panel_B_size = [9.5 1];
panel_A    = axes('position', [ 1.25 23.2  2 2]);
panel_B(1) = axes('position', [ 4.5 24.5  panel_B_size]);
panel_B(2) = axes('position', [ 4.5 23.5  panel_B_size]);
panel_B(3) = axes('position', [ 4.5 22.5  panel_B_size]);
panel_B(4) = axes('position', [ 4.5 21.5  panel_B_size]);
panel_C    = axes('position', [15.2 21.0  5 5]);
panel_D    = axes('position', [ 1.1 18 8 2.3]);
panel_E    = axes('position', [ 0   12.5 9 5]);
panel_F    = axes('position', [10.5 18 2 2]);
panel_G    = axes('position', [14.1 18 6 2]);
panel_H    = axes('position', [10.5 14 2.5 2.5]);
panel_I    = axes('position', [14.0 14.25 3 2.5]);
panel_J    = axes('position', [ 1.5 10.25 6 2]);
panel_K    = axes('position', [ 9.0 10.25 2.5 2]);
panel_L    = axes('position', [13.0 10.25 2.5 2]);
panel_M    = axes('position', [17.0 10.25 2.5 2]);

%
prm = PARAMS_GetAll();

%% ---------------------- spikes trace ------------------------------------
% define option for data
exp_ID = 'b0148_d170625';
exp = exp_load_data(exp_ID,'details','path');
TT = 4;
trace_duration = 0.33;
trace_ts_opt = 6;
trace_ts_list = [
    70270204396,...
    70132623245,...
    74986863162,...
    74986792310,...
    72188230643,...
    71432461878 % in flight!
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
    plot(ts,signal(ch,:),'k','LineWidth',0.5);
    box off
    set(gca,'Visible','Off')
    xlim(ts([1 end]))
end
linkaxes(panel_B,'xy')

%% add time/voltage scales
axes(panel_B(4));
xlimits = get(gca,'xlim');
ylimits = get(gca,'ylim');
scale_ms = 20;
scale_uVolt = 500;
scale_line_width = 1;
xa = xlimits(1) + [0 scale_ms*1e3];
ya = ylimits(1) + [0 scale_uVolt ];
[xaf,yaf] = ds2nfu(xa,ya);
% xaf = xaf + 0.01;
% yaf = yaf - 0.02;
xaf = xaf - 0.035;
yaf = yaf + 0.02;
annotation('line',    xaf([1 2]),yaf([1 1]),'Linewidth',scale_line_width); % time
annotation('line',    xaf([1 1]),yaf([1 2]),'Linewidth',scale_line_width)  % voltage
h=annotation('textbox',[mean(xaf) mean(yaf)-0.003 0 0],'String',sprintf('%dms',scale_ms));
h.Text.HorizontalAlignment = 'Center';
h.Text.VerticalAlignment = 'Top';
h.Text.Position = [mean(xaf) yaf(1) 0];
h.Text.FontSize = 7;
h=annotation('textbox',[mean(xaf)-0.03 mean(yaf)+0.005 0 0],'String','');
h.Text.Rotation = 90;
% h.Text.String = sprintf('%duV',scale_uVolt);
h.Text.String = [ num2str(scale_uVolt) ' {\mu}V'];
h.Text.Position = [xaf(1) mean(yaf) 0];
h.Text.HorizontalAlignment = 'Center';
h.Text.VerticalAlignment = 'bottom';
h.Text.FontSize = 7;
h.FitBoxToText='on';
h.LineStyle='none';

% add panel letter
axes(panel_B(1));
text(-0.05,0.7, 'B', 'Units','normalized','FontWeight','bold');

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
% clusters_colors([1 2 8 9],:) = clusters_colors([7 11 9 8],:);
% clusters_colors = clusters_colors(randperm(size(clusters_colors,1)),:);

%% identify the single-units
clusters_single_unit = [1 1 1 1 1 1 1 0 1 0 0 0];
% clusters_single_unit = [1 0 0 0 0 0 0 0 0 0 0 0];

%% ---------------------- spikes clusters ---------------------------------
NTT_file = 'L:\Analysis\pre_proc\0148\20170625\spikes_sorting\spikes_b0148_d170625_TT4.NTT';
exp_ID = 'b0148_d170625';
exp=exp_load_data(exp_ID,'details','flight');
prm = PARAMS_GetAll();
FE=exp.flight.FE;
FE([FE.distance] < prm.flight.full_min_distance) = [];
behave_session_ts = exp_get_sessions_ti(exp_ID, 'Behave');
switch panel_C_opt
    case 1 
        % sleep session
        panel_C_opt_str = 'sleep_session';
        limits_ts = [70019176845 70816687098];
        [Timestamps, CellNumbers, Samples, Header] = Nlx2MatSpike(NTT_file, [1 0 1 0 1], 1, 4, limits_ts);
    case 2 
        % all spikes
        panel_C_opt_str = 'all_spikes';
        limits_ts = [];
        [Timestamps, CellNumbers, Samples, Header] = Nlx2MatSpike(NTT_file, [1 0 1 0 1], 1, 1, [] );
    case 3
        % all spikes during flight (in-air)
        panel_C_opt_str = 'flight_all';
        limits_ts = behave_session_ts;
        [Timestamps, CellNumbers, Samples, Header] = Nlx2MatSpike(NTT_file, [1 0 1 0 1], 1, 4, limits_ts);
        IX = get_data_in_ti(Timestamps, flight_ts);
        Timestamps = Timestamps(IX);
        Samples = Samples(:,:,IX);
        CellNumbers = CellNumbers(IX);
    case 4
        % spikes during flight (in-air) - 1st half
        panel_C_opt_str = 'flight_1st_half';
        limits_ts = behave_session_ts;
        [Timestamps, CellNumbers, Samples, Header] = Nlx2MatSpike(NTT_file, [1 0 1 0 1], 1, 4, limits_ts);
        FE(length(FE)/2:end) = [];
        flight_ts = [FE.start_ts; FE.end_ts]';
        IX = get_data_in_ti(Timestamps, flight_ts);
        Timestamps = Timestamps(IX);
        Samples = Samples(:,:,IX);
        CellNumbers = CellNumbers(IX);
    case 5
        % spikes during flight (in-air) - 2nd half
        panel_C_opt_str = 'flight_2nd_half';
        limits_ts = behave_session_ts;
        [Timestamps, CellNumbers, Samples, Header] = Nlx2MatSpike(NTT_file, [1 0 1 0 1], 1, 4, limits_ts);
        FE(1:length(FE)/2) = [];
        flight_ts = [FE.start_ts; FE.end_ts]';
        IX = get_data_in_ti(Timestamps, flight_ts);
        Timestamps = Timestamps(IX);
        Samples = Samples(:,:,IX);
        CellNumbers = CellNumbers(IX);
    case 6
        % all spikes during behave session (in-air+on balls)
        panel_C_opt_str = 'behave_session';
        limits_ts = behave_session_ts;
        [Timestamps, CellNumbers, Samples, Header] = Nlx2MatSpike(NTT_file, [1 0 1 0 1], 1, 4, limits_ts);
end
ADBitVolts = sscanf(Header{contains(Header,'ADBitVolts')},'-ADBitVolts %f %f %f %f');
Samples = Samples .* ADBitVolts' .* 1e6; % convert bits to uVolts
spikes_size = squeeze(range(Samples,1));
cells = unique(CellNumbers);
color_list = linspecer( length(cells) ,'qualitative');
rng(18)
color_list = color_list(randperm(size(color_list,1)),:);
color_list(1,:) = 0.5.*[1 1  1];
color_list(2,:) = [0.2 0.2 0.5];
% plot
axes(panel_C);
% figure
cla
hold on
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
    if clusters_single_unit(cellNum)
        plot3(x,y,z,'.','Color',color_list(ii_cell,:),'MarkerSize',1);
    end
%     cluCntr = prctile([x;y;z]',95);
%     text(cluCntr(1),cluCntr(2),cluCntr(3), num2str(cellNum));
%     if clusters_single_unit(cellNum)
%         plot3(cluCntr(1),cluCntr(2),cluCntr(3),'*','Color','r','MarkerSize',10);
%     end
end
%%
ha = gca;
ha.XLim = [0 max(X(CellNumbers ~= 0))];
ha.YLim = [0 max(Y(CellNumbers ~= 0))];
ha.ZLim = [0 max(Z(CellNumbers ~= 0))];
ha.XTick = floor(ha.XLim/100)*100;
ha.YTick = floor(ha.YLim/100)*100;
ha.ZTick = floor(ha.ZLim/100)*100;

h1=text([230 1610],[650 680],  [0 0], ""+ha.XRuler.TickValues,'FontSize',7);
h2=text([1560 1810],[140 683], [-20 15], ""+ha.YRuler.TickValues,'FontSize',7);
h3=text([1645 1740],[140 120], [70 960], ""+ha.ZRuler.TickValues,'FontSize',7);
ha.XRuler.TickLabels = {};
ha.YRuler.TickLabels = {};
ha.ZRuler.TickLabels = {};

% ha.XRuler.TickLabelGapMultiplier = -0.5;
% ha.YRuler.TickLabelGapMultiplier = -0.5;
% ha.ZRuler.TickLabelGapMultiplier = -0.5;
%%
AZ_EL(1,:) = [140 39];
AZ_EL(2,:) = [488 20];
AZ_EL(3,:) = [511 -8];
AZ_EL(4,:) = [876 -12];
AZ_EL(5,:) = [156 8];
AZ_EL(6,:) = [-41 -14];
AZ_EL(7,:) = [138 12];
AZ_EL(8,:) = [146 16];
viewing_option = 5;
view(AZ_EL(viewing_option,:));

hx=xlabel(['Ch' num2str(ch2plot(1)) ' ( {\mu}V)']);%,'Position',[835  250  -152]);
hy=ylabel(['Ch' num2str(ch2plot(2)) ' ( {\mu}V)']);%,'Position',[100  300  -100]);
hz=zlabel(['Ch' num2str(ch2plot(3)) ' ( {\mu}V)']);%,'Position',[1400  -44  480]);
hx.Units = 'normalized';
hy.Units = 'normalized';
hz.Units = 'normalized';
hx.Position = [0.58 -0.005 0];
hy.Position = [0.20 -0.050 0];
hz.Position = [-0.050  0.5300 0];
hx.Rotation = 2;
hy.Rotation = -24;
hz.Rotation = 90;

text(-0.215,0.95, 'C', 'Units','normalized','FontWeight','bold');

%% show images
% logger_image_filename = 'D:\Tamir\PROJECTS\Neurologger\miniBat\pics\minibat_good_pic.jpg';
logger_image_filename = 'L:\resources\minibat\minibat3.jpg';
% tunnel_view_image_file = 'L:\Videos_Photos\tunnel_area_various\20170111_170850_downsampled.jpg';
% tunnel_view_image_file = 'L:\Videos_Photos\TAZOT_HAMAMA\taza3.jpg';
tunnel_view_image_file = 'L:\Videos_Photos\TAZOT_HAMAMA\taza4.jpg';

axes(panel_A);
image = imread(logger_image_filename);
imshow(image);
text(-0.29,1.22, 'A', 'Units','normalized','FontWeight','bold');
% add scale bar
scale_mm = 10;
pixel_mm_ratio = 720/11; % 720 pixels is measured manually using ginput amd sd card width is 11mm
scale_line_width = 1.5;
scale_pixels = scale_mm * pixel_mm_ratio;
xlimits = get(gca,'xlim');
ylimits = get(gca,'ylim');
xa = xlimits(1) + [0 scale_pixels];
ya = ylimits(1) + [0 0];
[xaf,yaf] = ds2nfu(xa,ya);
xaf = xaf + 0.04;
yaf = yaf + 0.008;
annotation('line', xaf,yaf, 'Linewidth',scale_line_width);
h=annotation('textbox', [mean(xaf)-0.0005 mean(yaf)-0.008 0 0], 'String', sprintf('%dmm',scale_mm),...
    'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',8);

axes(panel_D);
image = imread(tunnel_view_image_file);
imshow(image);
axis image
% add scale bar
scale_m = 20;
pixel_m_ratio = 648/143; % manually measurement with ginput
scale_line_width = 2;
scale_pixels = scale_m * pixel_m_ratio;
xlimits = get(gca,'xlim');
ylimits = get(gca,'ylim');
xa = xlimits(2) - [0 scale_pixels];
ya = ylimits(1) + [0 0];
[xaf,yaf] = ds2nfu(xa,ya);
xaf = xaf - 0.01;
yaf = yaf + 0.008;
annotation('line', xaf,yaf, 'Linewidth',scale_line_width);
h=annotation('textbox', [mean(xaf) mean(yaf)-0.008 0 0], 'String', sprintf('%dm',scale_m),...
    'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',8);
text(-0.05,1.22, 'D', 'Units','normalized','FontWeight','bold');

%% panel I - tunnel section behavior (ZY)
axes(panel_I)
cla
axis equal
hold on
text(-0.15,1, 'I', 'Units','normalized','FontWeight','bold');
% plot tunnel section lines
plot([-1.15 -1.15],[0 1.7],'k','LineWidth',1.5);
plot([ 1.15  1.15],[0 1.7],'k','LineWidth',1.5);
plot([ 1.15 0],[1.7 2.35],'k','LineWidth',1.5);
plot([-1.15 0],[1.7 2.35],'k','LineWidth',1.5);
plot([-1.15 1.15],[0 0],'k','LineWidth',1.5);
% add Y/Z arrows
xlimits = get(gca,'xlim');
ylimits = get(gca,'ylim');
xa = -1.5 + [0 0.75];
ya = -0.25 + [0 0.75];
[xaf,yaf] = ds2nfu(xa,ya);
xaf = xaf;
yaf = yaf;
h(1)=annotation('arrow',xaf,yaf([1 1]), 'Color', 'k');
h(2)=annotation('arrow',xaf([1 1]),yaf, 'Color', 'k');
[h.HeadWidth] = disperse([5 5]);
[h.HeadLength] = disperse([5 5]);
[h.LineWidth] = disperse([1.5 1.5]);
% h=annotation('textbox', [mean(xaf) mean(yaf) 0 0]+[-0.005 -0.025 0 0], 'String', 'Y',...
%     'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',8);
% h=annotation('textbox', [mean(xaf) mean(yaf) 0 0]+[-0.03 -0.005 0 0], 'String', 'Z',...
%     'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',8);
h=annotation('textbox', [mean(xaf) mean(yaf) 0 0]+[+0.025 -0.015 0 0], 'String', 'Y',...
    'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',8);
h=annotation('textbox', [mean(xaf) mean(yaf) 0 0]+[-0.018 +0.02 0 0], 'String', 'Z',...
    'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',8);
% h.Text.Rotation = 90;
% add scale bar
scale_m = 0.5;
scale_line_width = 2;
xlimits = get(gca,'xlim');
ylimits = get(gca,'ylim');
xa = xlimits(1) + [0 scale_m];
ya = ylimits(1) + [0 0];
[xaf,yaf] = ds2nfu(xa,ya);
xaf = xaf + 0.09;
yaf = yaf - 0.01;
annotation('line', xaf,yaf, 'Linewidth',scale_line_width);
h=annotation('textbox', [mean(xaf) mean(yaf)-0.008 0 0], 'String', sprintf('%dcm',scale_m*100),...
    'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',8);

% panel_H_opt = 1;
panel_I_data_options = {
    'b0034_d180313', 1, 160;...
    'b0034_d180314', 1, 175;...
    'b0079_d160928', 1, 135;...
    'b9861_d180626', 1, 167;...
    'b9861_d180628', 1, 159;...
    'b9861_d180709', 1, 154;...
    'b9861_d180711', 1, 159;...
    };
YZ_dev_plot_opt = 'data';
% YZ_dev_plot_opt = 'illustration';
switch YZ_dev_plot_opt
    case 'data'
        exp_ID = panel_I_data_options{panel_I_opt,1};
        direction = panel_I_data_options{panel_I_opt,2};
        x0 = panel_I_data_options{panel_I_opt,3};
        load("L:\Analysis\Results\exp\pos_XYZ\"+exp_ID+"_pos_XYZ.mat");
        YZ0 = XYZ.YZ_by_dir{direction};
        YZ_mean = nanmean(YZ0,1);
%         YZ_mean(1,1,:) = 0;
        YZ0 = YZ0 - YZ_mean; % remove data mean
        YZ0 = YZ0 + [0 1.5]; % add back fixed mean
        plot(YZ0(:,1,x0),YZ0(:,2,x0),'.','Color',prm.graphics.colors.flight_directions{direction});
        set(gca,'Visible','off');
    case 'illustration'
        rng(0)
        plot(0.1*randn(1,20),1.5+0.1*randn(1,20),'.b')
        % compute points corresponding to axis-oriented ellipse
        r1 = 0.4;
        r2 = 0.4;
        xc = 0;
        yc = -1.5;
        theta = 0;
        % compute points corresponding to axis-oriented ellipse
        t = linspace(0, 2*pi, 200);
        xt = r1 * cos(t) + xc;
        yt = r2 * sin(t) + yc;
        % aply rotation by angle theta
        cot = cos(theta); sit = sin(theta);
        x = xt * cot - yt * sit;
        y = xt * sit - yt * cot;
        % draw the curbe
        plot(x, y, '-');
        text(0,0.5,'Illustration','HorizontalAlignment','center')
end
xlim([-1.5 1.5])
ylim([ 0 2.5])
xlabel('Width (m)')
ylabel('Height (m)')
ha=gca;
ha.XRuler.TickLabelGapMultiplier = -0.5;
ha.YRuler.TickLabelGapMultiplier = 0;

%% bespoon localization (anchors+tag+tunnel)
axes(panel_E);
cla
text(-0.11, 1, 'E', 'Units','normalized','FontWeight','bold');
axis equal
% axis normal
pause(eps)
hold on
exp_ID = 'b2289_d180615';
exp = exp_load_data(exp_ID, 'pos');
x = exp.pos.calib_tunnel.curvexy(:,1);
y = exp.pos.calib_tunnel.curvexy(:,2);
x = x(1:50:end);
y = y(1:50:end);
tunnel_plot = 'middle';
switch tunnel_plot
    case 'middle'
        plot(x,y, '-k','LineWidth',3);
    case 'walls';
        [joinedx, joinedy] = offsetCurve(x, y, 1.25);
        plot(joinedx, joinedy, '-k','LineWidth',0.000001)
        [joinedx, joinedy] = offsetCurve(x, y, -1.25);
        plot(joinedx, joinedy, '-k','LineWidth',0.000001);
end
set(gca,'Visible','off');
anchors_pos = [
1491.223	2464.864
1428.952	2460.655
1427.693	2377.064
1307.676	2408.918
1298.294	2472.189
1384.651	2469.469
1357.793	2487.461
1349.005	2448.650
1453.588	2514.633
1420.749	2516.462
1392.801	2512.346
1350.107	2499.272
1329.574	2484.206
1493.695	2500.963
];
plot(anchors_pos(:,1),anchors_pos(:,2),'or','MarkerSize',2,'MarkerFaceColor',[1 0 0]);
xy = [1419 2500];
for ii_anchor = 1:size(anchors_pos,1)
    anchor_pos = anchors_pos(ii_anchor,:);
    r = pdist([xy;anchor_pos]);
    circle_pos = anchor_pos - r;
    circle_pos = [circle_pos 2*r 2*r];
    h=rectangle('Position',circle_pos, 'Curvature',[1 1],'LineWidth',0.00001);
    h.LineStyle=  '-';
%     plot([anchor_pos(1) xy(1)], [anchor_pos(2) xy(2)],'-g')
end
plot(xy(1),xy(2),'.b','MarkerSize',15)
xlim([1280 1515])
ylim([2360 2541])
pause(eps)
TAZA_rotation = 9;
north_rotation = -3;
view(TAZA_rotation,90); % rotate this axis to match the TAZA rotation
% view(0,90);
pause(eps)

% add north arrow
h=annotation('arrow',[0 0],[0 0]);
h.Units = 'centimeters';
center = [7.8 15.8];
l = 1;
width  = l*cos(deg2rad(90-(TAZA_rotation+north_rotation)));
height = l*sin(deg2rad(90-(TAZA_rotation+north_rotation)));
h.Position = [center width height];
h.HeadLength = 8;
h.HeadWidth = 8;
h=text(1525,2570,'North','FontSize',9,'HorizontalAlignment','center');
h.Rotation = -(TAZA_rotation+north_rotation);

% add scale bar
scale_m = 20;
scale_line_width = 2;
xlimits = get(gca,'xlim');
ylimits = get(gca,'ylim');
xa = xlimits(1) + [0 scale_m];
ya = ylimits(1) + [0 0];
[xaf,yaf] = ds2nfu(xa,ya);
xaf = xaf + 0.06;
yaf = yaf + 0.1;
annotation('line', xaf,yaf, 'Linewidth',scale_line_width);
h=annotation('textbox', [mean(xaf) mean(yaf)-0.008 0 0], 'String', sprintf('%dm',scale_m),...
    'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',8);

% add legend
leg_pos = get(panel_E,'position');
leg_pos(1:2) = leg_pos(1:2) + [0.65 0.88].*leg_pos(3:4);
leg_pos(3:4) = [1 0.6];
leg_ax = axes('position',leg_pos);
hold on
plot(0,1,'.b','MarkerSize',10)
plot(0,2,'.r','MarkerSize',10);
% plot(0,2,'or','MarkerSize',10,'MarkerFaceColor',[1 0 0]);
text(0.18,1,'Bat','FontSize',7);
text(0.18,2,'Antenna','FontSize',7);
xlim([0 1])
ylim([0.5 2.5])
box off
set(gca,'Visible','off');
set(gca,'XTick',[],'YTick',[]);

%% bespoon localization precision
axes(panel_F);
cla
hold on
bespoon_loc_precision = load('L:\BeSpoon\testing\test_20180530__YOM_KEF_200m_static+dynamic+discretization\30-05-2018__calib_test_dynamic+non-jitter_jitter_with_kalman\data\outside_perpendicular_error.mat');
err = bespoon_loc_precision.perpendicular_error;
err = err.*100; % convert to cm
h=histogram(err);
h.NumBins = 15;
h.Normalization = 'pdf';
h.FaceColor = 0.5*[1 1 1];
h.EdgeColor = 0.5*[1 1 1];
h.LineWidth = 0.1;
[muHat,sigmaHat] = normfit(err);
x = linspace(-50,50,100);
y = normpdf(x,muHat,sigmaHat);
plot(x,y,'k','LineWidth',1.5);
text(0.65,0.85, sprintf('\x03C3=%.1fcm', sigmaHat), 'Units','normalized','FontSize',8);
xlabel({'Localization error (cm)'},'Units','normalized','Position',[0.5 -0.2]);
ylabel('Probability','Units','normalized','Position',[-0.07 0.5]);
ha = gca;
ha.XLim = [-40 40];
ha.XTick = [-40:20:40];
ha.YTick = ha.YLim;
ha.TickDir='out';
ha.TickLength = [0.03 0.03];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.25;
text(-0.45,1.15, 'F', 'Units','normalized','FontWeight','bold');

%% behavioral trajectory is 1D (small y deviations) - example
axes(panel_G);
cla
hold on
text(-0.115,1.15, 'G', 'Units','normalized','FontWeight','bold');
% panel_G_opt = 1;
panel_G_data_options = {
'b0034_d180413'; % TODO: verify this option is from a CA1 day!!!!
'b2289_d180514';
'b2289_d180515';
'b2289_d180518';
'b9861_d180705';
'b9861_d180709';
};
exp = exp_load_data(panel_G_data_options{panel_G_opt},'flight');
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
xlabel('Position (m)','Units','normalized','Position',[0.5 -0.15]);
ylabel('Y (m)','Units','normalized','Position',[-0.055 0.5]);
ha = gca;
ha.XLim = [0 200];
ha.YLim = [-1.5 1.5];
ha.XTick = [0:50:200];
ha.YTick = [-2:1:2];
ha.TickDir='out';
ha.TickLength = [0.01 0.01];
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;

% add direction arrows
arrow_x = [0.80 0.85]+0.0231;
arrow_y = repelem(0.74,2);
clear h
h(1)=annotation('arrow',arrow_x,      arrow_y+0.01, 'Color', prm.graphics.colors.flight_directions{1});
h(2)=annotation('arrow',flip(arrow_x),arrow_y     , 'Color', prm.graphics.colors.flight_directions{2});
[h.HeadWidth] = disperse([5 5]);
[h.HeadLength] = disperse([5 5]);


%% load population data
% get list of significant cells (at least in one direction)
% ---------------------------------------------------------
cells_t = DS_get_cells_summary();
% choose bats
bats = [79,148,34,9861,2289];
cells_t(~ismember(cells_t.bat, bats ),:) = [];
% choose good clusters from CA1
cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells = [cells.details];
cells(~contains({cells.brain_area}, 'CA1')) = [];
cells(~ismember([cells.ClusterQuality], [2])) = [];
% choose pyramidal cells only
cells = cellfun(@(c)(cell_load_data(c,'details','stats')), {cells.cell_ID}, 'UniformOutput',0);
cells = [cells{:}];
stats = [cells.stats];
stats = [stats.all];
cells([stats.meanFR_all]>prm.inclusion.interneuron_FR_thr)=[];
% choose only signif cells (in any direction)
cells_details = [cells.details];
cells = cellfun(@(c)(cell_load_data(c,'details','signif')), {cells_details.cell_ID}, 'UniformOutput',0);
cells = [cells{:}];
signif = cat(1,cells.signif);
signif = arrayfun(@(x)(x.TF), signif);
signif = any(signif,2);
cells(~signif)=[];

% load the final list of cells
cells_details = [cells.details];
clear cells stats signif cells_t
cells = cellfun(@(c)(cell_load_data(c,'details','stats','meanFR','stats','inclusion','signif','fields','FR_map')), {cells_details.cell_ID}, 'UniformOutput',0);
cells = [cells{:}];

% load exp data
cells_details = [cells.details];
exps_IDs = {cells_details.exp_ID};
exp_list = unique(exps_IDs);
exps = cellfun(@(x)(exp_load_data(x,'details','flight')), exp_list);
exps_details = [exps.details];
exps_flight = [exps.flight];
pos_y_dev_all = cat(1,exps_flight.pos_y_std);
speed_traj_all = cat(1,exps_flight.speed_traj);
ystd_median_all = arrayfun(@(x)(x.ystd_median), pos_y_dev_all);

%% behavioral trajectory is 1D (small y deviations) - population
axes(panel_H);
cla
text(-0.34,1.1, 'H', 'Units','normalized','FontWeight','bold');
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
yvar_pop_plot = 'ksdensity_pooled_bats';
% yvar_pop_plot = 'violin';
% yvar_pop_plot = 'boxplot';
% yvar_pop_plot = 'median_std';
% yvar_pop_plot = 'none';
switch yvar_pop_plot 
    case 'ksdensity_pooled_bats'
        for ii_dir = 1:2
            ystd_by_dir = cat(1,pos_y_dev_all(:,ii_dir).ystd); % pool over sessions
            ystd_by_dir = ystd_by_dir(:,valid_pos_IX); % take only valid position (high-speed)
            ystd_by_dir = ystd_by_dir(:); % pool over positions
            ystd_by_dir = ystd_by_dir .* 100; % convert to cm
            xi = linspace(0,50,100);
            ystd_density = ksdensity(ystd_by_dir,xi);
            plot(xi,ystd_density,'Color',prm.graphics.colors.flight_directions{ii_dir},'LineWidth',1.5);
        end
        xlabel('Y s.d. (cm)');
        ylabel('Probability','Units','normalized','Position',[-0.1 0.5]);
        ha=gca;
        ha.XLim = [0 40];
        ha.XTick = 0:10:40;
        ha.YLim = [0 0.1];
        ha.YTick = [0 0.1];
        ha.TickLength = [0.03 0.03];
        ha.XRuler.TickLabelGapMultiplier = -0.3;
        ha.YRuler.TickLabelGapMultiplier = 0.001;
    case 'violin'
        hv=violinplot(x,grps);
        [hv.BoxWidth] = disperse(repelem(0.02,length(hv)));
        [hv.EdgeColor ] = disperse(repelem({'none'},length(hv)));
        hs = [hv.ScatterPlot];
        [hs.SizeData] = disperse(repelem(5,length(hs)));
        hm = [hv.MedianPlot];
        [hm.SizeData] = disperse(repelem(10,length(hm)));
        ha=gca;
        ha.YLim = [0 90];
        ha.YTick = [0:30:90];
    case 'boxplot'
        boxplot(x,grps)
        box off
        ha=gca;
        ha.YLim = [0 90];
        ha.YTick = [0:30:90];
    case 'median_std'
        medians = cellfun(@median,data);
        stds    = cellfun(@std,data);
        errorbar(medians,stds,'-o','MarkerSize',3);
        ha=gca;
        ha.YLim = [0 40];
        ha.YTick = [0:20:40];
    case 'none'
        ha=gca;
end
switch yvar_pop_plot
    case {'ksdensity_pooled_bats'}
    otherwise
        ha.XTick = 1:length(bats);
        ha.XTickLabel = bats;
        ha.XLim = [0.5 length(bats)+.5];
        ha.TickDir='out';
        ha.TickLength = [0.01 0.01];
        ha.XRuler.TickLabelGapMultiplier = -0.3;
        ha.YRuler.TickLabelGapMultiplier = 0.001;
        xlabel('Bats');
        ylabel('Y s.d. (cm)');
end

% add direction arrows
arrow_x = 0.53 + [0 0.03] + 0.0231;
arrow_y = repelem(0.60,2);
clear h
h(1)=annotation('arrow',arrow_x,      arrow_y+0.01, 'Color', prm.graphics.colors.flight_directions{1});
h(2)=annotation('arrow',flip(arrow_x),arrow_y     , 'Color', prm.graphics.colors.flight_directions{2});
[h.HeadWidth] = disperse([5 5]);
[h.HeadLength] = disperse([5 5]);


%% speed trajectory very constant along the flight - example
axes(panel_J); 
cla
hold on
text(-0.13,1.1, 'J', 'Units','normalized','FontWeight','bold');
% panel_J_opt = 4;
panel_J_data_options = {
    'b0034_d180313';
    'b0034_d180315';
    'b0079_d161003';
    'b0148_d170613';
    'b0148_d170615';
    'b0148_d170626';
    'b0148_d170802';
    'b0148_d170807';
    'b2289_d180520';    
    'b2289_d180615';
};
exp_ID = panel_J_data_options{panel_J_opt};
exp = exp_load_data(exp_ID,'flight');
FE=exp.flight.FE;
FE([exp.flight.FE.distance]<prm.flight.full_min_distance) = [];
x = [FE.pos];
y = [FE.vel];
y = abs(y);
ylimits = [0 9];
baseval = 0.1;
area([4 prm.fields.valid_speed_pos(1)]    , ylimits([2 2]), baseval, 'FaceColor',0.8*[1 1 1],'EdgeColor','none','ShowBaseLine','off');
area([  prm.fields.valid_speed_pos(2) 194], ylimits([2 2]), baseval, 'FaceColor',0.8*[1 1 1],'EdgeColor','none','ShowBaseLine','off');
plot(x,y,'.','Color', 'k','MarkerSize',1);
% plot(prm.fields.valid_speed_pos([1 1]), ylimits,'--m','LineWidth',2)
% plot(prm.fields.valid_speed_pos([2 2]), ylimits,'--m','LineWidth',2)
set(gca,'xtick',0:50:200,'ytick',[-10 0 8],'xlim',[0 200])
set(gca,'tickdir','out','TickLength',repelem(0.01,2));
ylim(ylimits);
xlabel('Position (m)','Units','normalized','Position',[0.5 -0.25]);
ylabel('Flight speed (m/s)','Units','normalized','Position',[-0.07 0.41]);


%% panel K - speed trajectory very constant along the flight - population
axes(panel_K);
cla
hold on
text(-0.4,1.1, 'K', 'Units','normalized','FontWeight','bold');
dir_colors = prm.graphics.colors.flight_directions;
for ii_dir = 1:2
    cv = [speed_traj_all(:,ii_dir).speed_cv];
%     cv = [cv.raw_high_speed];
    cv = [cv.across_pos_high_speed];
    mean(cv)
    length(cv)
    h=histogram(cv);
    h.BinEdges = linspace(0,0.1,13);
    h.BinEdges = h.BinEdges + (ii_dir-1)*0.3*h.BinWidth;
    h.FaceColor = [1 1 1];
    h.EdgeColor = dir_colors{ii_dir};
    h.Normalization = 'Count';
    h.FaceAlpha = 0;
    h.LineWidth = 1;
end
ha=gca;
ha.XTick = [0 0.04 0.08];
ha.TickLength = [0.04 0.04];
ha.XRuler.TickLabelGapMultiplier = -0.1;
ha.YRuler.TickLabelGapMultiplier = 0.2;
xlabel('CV of speed','Units','normalized','Position',[0.5 -0.25]);
ylabel('No. of sessions','Units','normalized','Position',[-0.2 0.5]);

% add direction arrows
arrow_x = [0.48 0.50] + 0.0231;
arrow_y = repelem(0.45,2);
clear h
h(1)=annotation('arrow',arrow_x,      arrow_y+0.01, 'Color', prm.graphics.colors.flight_directions{1});
h(2)=annotation('arrow',flip(arrow_x),arrow_y     , 'Color', prm.graphics.colors.flight_directions{2});
[h.HeadWidth] = disperse([5 5]);
[h.HeadLength] = disperse([5 5]);

%% panel L - number of laps
axes(panel_L);
cla
hold on
text(-0.4,1.1, 'L', 'Units','normalized','FontWeight','bold');
dir_colors = prm.graphics.colors.flight_directions;
directions = [-1 1];
for ii_dir = 1:2
    direction = directions(ii_dir);
    FE_dir = cellfun(@(FE)(FE([FE.direction]==direction)), {exps_flight.FE},'UniformOutput',0);
    nFlights = cellfun(@(FE)(sum([FE.distance]>prm.flight.full_min_distance)),FE_dir);
    h=histogram(nFlights);
    h.NumBins = 10;
    h.BinEdges = h.BinEdges + (ii_dir-1)*0.3*h.BinWidth;
    h.FaceColor = [1 1 1];
    h.EdgeColor = dir_colors{ii_dir};
    h.Normalization = 'Count';
    h.FaceAlpha = 0;
    h.LineWidth = 1;
end
ha=gca;
ha.XLim(1) = 0;
ha.XLim(2) = 65;
ha.TickLength = [0.04 0.04];
ha.XRuler.TickLabelGapMultiplier = -0.1;
ha.YRuler.TickLabelGapMultiplier = 0.2;
xlabel('No. of flights','Units','normalized','Position',[0.5 -0.25]);
ylabel('No. of sessions','Units','normalized','Position',[-0.2 0.5]);

% add direction arrows
arrow_x = 0.67 + [0 0.02] + 0.0231;
arrow_y = repelem(0.45,2);
clear h
h(1)=annotation('arrow',arrow_x,      arrow_y+0.01, 'Color', prm.graphics.colors.flight_directions{1});
h(2)=annotation('arrow',flip(arrow_x),arrow_y     , 'Color', prm.graphics.colors.flight_directions{2});
[h.HeadWidth] = disperse([5 5]);
[h.HeadLength] = disperse([5 5]);

%% panel M - Total distance 
axes(panel_M);
cla 
hold on
text(-0.44,1.1, 'M', 'Units','normalized','FontWeight','bold');
total_dist_by_dir = 0;
if total_dist_by_dir
    %% (per direction)
    dir_colors = prm.graphics.colors.flight_directions;
    directions = [-1 1];
    for ii_dir = 1:2
        direction = directions(ii_dir);
        FE_dir = cellfun(@(FE)(FE([FE.direction]==direction)), {exps_flight.FE},'UniformOutput',0);
        total_distance = cellfun(@(FE)(sum([FE.distance])),FE_dir);
        total_distance = total_distance .*1e-3; % to km
        h=histogram(total_distance);
        h.NumBins = 10;
        h.BinEdges = h.BinEdges + (ii_dir-1)*0.3*h.BinWidth;
        h.FaceColor = 0.5.*[1 1 1];
        h.EdgeColor = dir_colors{ii_dir};
        h.Normalization = 'Count';
        h.FaceAlpha = 0;
        h.LineWidth = 1;
    end
else
    %% Total distance (pooling directions)
    dir_colors = prm.graphics.colors.flight_directions;
    total_distance = cellfun(@(FE)(sum([FE.distance])),{exps_flight.FE});
    total_distance = total_distance .*1e-3; % to km
    h=histogram(total_distance);
    fprintf('Total distance (mean) = %.1f\n', mean(total_distance));
    fprintf('Total distance (max) = %.1f\n', max(total_distance));
    length(total_distance)
    
    nBinEdges = 12;
    h.BinEdges = linspace(0,25,nBinEdges);
    h.FaceColor = 0.5.*[1 1 1];
    h.EdgeColor = 'k';
    h.Normalization = 'Count';
    h.LineWidth = 1;
    ha=gca;
    ha.XLim = [0 25];
    ha.XTick = [0:5:25];
end
ha=gca;
ha.TickLength = [0.04 0.04];
ha.XRuler.TickLabelGapMultiplier = -0.1;
ha.YRuler.TickLabelGapMultiplier = 0.2;
xlabel('Distance flown (km)','Units','normalized','Position',[0.5 -0.25]);
ylabel('No. of sessions','Units','normalized','Position',[-0.2 0.5]);

%%
data_opt_str = {
    'data options:'
    sprintf('panel C - option %d: %s               ',panel_C_opt,panel_C_opt_str);
    sprintf('panel G - option %d: %s               ',panel_G_opt,panel_G_data_options{panel_G_opt});
    sprintf('panel I - option %d: %s, dir=%d, x0=%d',panel_I_opt,panel_I_data_options{panel_I_opt,:});
    sprintf('panel J - option %d: %s               ',panel_J_opt,panel_J_data_options{panel_J_opt});
    
    };
% annotation('textbox', [0.2 0.2 0.6 0.1], 'String',data_opt_str, 'HorizontalAlignment','Left','Interpreter','none','FitBoxToText','on');

%% print/save the figure
fig_name_str = fig_name_str+"_C_opt"+panel_C_opt+panel_C_opt_str;
% fig_name_str = fig_name_str+"_G_opt"+panel_G_opt;
% fig_name_str = fig_name_str+"_I_opt"+panel_I_opt;
% fig_name_str = fig_name_str+"_J_opt"+panel_J_opt;
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');


%%






%%
