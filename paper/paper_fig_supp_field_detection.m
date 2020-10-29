%% Large Scale - Fig. supp XXX - field detection explanation

%%
clear 
clc

%% params
prm = PARAMS_GetAll();
prm.fields.valid_speed_pos = [10 187.5];
clr = prm.graphics.colors.flight_directions;

%% define output files
res_dir = hc3_get_res_dir();
res_dir = fullfile(res_dir,'paper_figures');
mkdir(res_dir)
fig_name_str = 'Fig_Sxxx_field_detection';
fig_caption_str = 'Field detection explanation';
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
panel_A_size_raster = [17 3];
panel_A_size_FR_map = [17 3];
panel_A_pos = [2 10];
panel_A = [];
for ii=1:1
    for jj=1:1
        offset_x = (ii-1)*6;
        offset_y = (jj-1)*4.3;
        offset = panel_A_pos + [offset_x offset_y];
        panel_A(ii,jj,1) = axes('position', [offset+[0 6.5 ] panel_A_size_FR_map]);
        panel_A(ii,jj,2) = axes('position', [offset+[0 3   ] panel_A_size_raster]);
        panel_A(ii,jj,3) = axes('position', [offset+[0 0   ] panel_A_size_raster]);
    end
end
panel_A = squeeze(panel_A)';

%% load data ==============================================================
load('L:\rodents_data\results\datasets\cells_bat_200m.mat');
% choose example
cell_example_opt = 1:9;
cell_examples_list = [
433;  56;  51;
609; 419; 477;
 57; 628; 337;
];
details = [cells.details];
cells_num = [details.cell_num];
IX = arrayfun(@(cell_num)(find(cell_num == cells_num)), cell_examples_list(cell_example_opt));
cell_examples = cells(IX);

%% load data ==============================================================
% cell_num = 609;
cell_num = 628;
switch 1
    case 1
        load(['L:\rodents_data\results\datasets\cell_' num2str(cell_num) '.mat']);
    case 2
        cell=cell_load_data(cell_num);
end

%% map+fields
axes(panel_A(1));
cla
hold on
maps=[cell.FR_map.all];
x = maps(1).bin_centers;
y = cat(1,maps.PSTH);
m = round(max(y(:)));
ylimits = [0 m+1];
box off
hax=gca;
hax.TickDir = 'out';
hax.XTick = [];
hax.YTick = [0 m];
hax.YLim = ylimits;
hax.XLim = [0 200];
% low-speed area
area_upperval = ylimits(1) - 1e-2*range(ylimits);
area_lowerval = ylimits(1) - 0.15*range(ylimits);
area([4 prm.fields.valid_speed_pos(1)]    , repelem(area_upperval,2), area_lowerval, 'FaceColor',0.7*[1 1 1],'EdgeColor','none','ShowBaseLine','off','Clipping','off');
area([  prm.fields.valid_speed_pos(2) 194], repelem(area_upperval,2), area_lowerval, 'FaceColor',0.7*[1 1 1],'EdgeColor','none','ShowBaseLine','off','Clipping','off');
h=plot(x,y);
[h.Color] = disperse(clr);
% fields
for ii_dir=1:2
    fields = cell.fields{ii_dir};
    if ~isempty(fields)
        fields([fields.in_low_speed_area])=[];
    end
    
    dir_offsets = [-0.065 -0.145]+0.015;
    for ii_field = 1:length(fields)
        % final detection
        x = fields(ii_field).edges_prc;
        y = dir_offsets(ii_dir)*range(hax.YLim);
        y = repelem(y,2);
        plot(x,y,'Linewidth', 2, 'Color', clr{ii_dir},'Clipping','off');
        % href
        x = fields(ii_field).edges_href;
        y = prm.fields.width_href * fields(ii_field).peak;
        y = repelem(y,2);
        plot(x,y,'Linewidth', 0.5, 'Color', clr{ii_dir},'Clipping','off');
    end
    text(1.025, dir_offsets(ii_dir), num2str(length(fields)), 'Units','normalized', 'Color', clr{ii_dir},...
        'HorizontalAlignment','right', 'VerticalAlignment','middle','FontSize',7);
end
%% rasters
FEs = [cell.FE];
for ii_dir=1:2
    axes(panel_A(ii_dir+1));
    cla
    hold on
    FE = FEs{ii_dir};
    x = [FE.spikes_pos];
    [FE.number2] = disperse(1:length(FE));
    y = arrayfun(@(FE)(FE.number2*ones(1,FE.num_spikes)),FE,'UniformOutput',0);
    y = [y{:}];
    plot(x,y,'.','Color',clr{ii_dir},'MarkerSize',5);
    box off
    hax=gca;
    m = length(FE);
    hax.YTick = [1 m];
    hax.XLim = [0 200];
    hax.YLim = [0 m+1];
    switch ii_dir
        case 1
            hax.XTick = [];
            hax.YTickLabel = {'',num2str(m)};
        case 2
            hax.XTick = 0:50:200;
            hax.XRuler.TickLabelGapOffset = -2;
            hax.YTickLabel = {'1',num2str(m)};
            hax.TickDir = 'out';
    end
    
    % field href area
    fields = cell.fields{ii_dir};
    if ~isempty(fields)
        fields([fields.in_low_speed_area])=[];
    end
    for ii_field = 1:length(fields)
        % final detection
        x = fields(ii_field).edges_href;
        y = hax.YLim([2 2]);
        area(x,y, hax.YLim(1), 'FaceColor',clr{ii_dir}, 'EdgeColor','none',...
            'ShowBaseLine','off','Clipping','off','FaceAlpha',0.1);
    end
    
end

%% add x/y labels for specific panels
axes(panel_A(1));
% text(-0.05,1.05, 'A', 'Units','normalized','FontWeight','bold');
ylabel({'Firing rate';'(Hz)'},   'Units','normalized','Position',[-0.03 0.5]);
axes(panel_A(3));
xlabel('Position (m)', 'Units','normalized','Position',[0.5 -0.2]);
ylabel('Flight no.',   'Units','normalized','Position',[-0.05 1]);

%% add direction arrows
arrow_x = 0.13 +[0 0.05];
arrow_y = repelem(0.75,2);
clear h
h(1)=annotation('arrow',arrow_x,      arrow_y+0.008,  'Color', prm.graphics.colors.flight_directions{1});
h(2)=annotation('arrow',flip(arrow_x),arrow_y      ,  'Color', prm.graphics.colors.flight_directions{2});
[h.HeadWidth] = disperse([5 5]);
[h.HeadLength] = disperse([5 5]);


%% print/save the figure
file_out = fig_name_str;
file_out = [file_out '_cell_' num2str(cell_num)];
file_out = fullfile(res_dir, file_out);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');

%%
