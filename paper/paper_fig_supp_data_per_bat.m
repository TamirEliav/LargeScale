%% Large Scale - Fig. supp XXX - data per bat

%%
clear 
clc

%% params
bat_ID_num_map = containers.Map([34 79 148 2289 9861],...
                                 1:5);

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S8';
fig_caption_str = 'Data per bat';
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
panel_ABC_size = [2.5 2.5];
panel_ABC_pos = [5 7];
panel_ABC = [];
for ii_bat=1:5
    for ii_feature=1:3
        offset_x = (ii_feature-1)*4.5;
        offset_y = (ii_bat-1)*3.25;
        offset = panel_ABC_pos + [offset_x offset_y];
        panel_ABC(ii_bat,ii_feature) = axes('position', [offset panel_ABC_size]);
    end
end
panel_ABC = flipdim(panel_ABC,1);

%% load data ==============================================================
load('L:\rodents_data\results\datasets\cells_bat_200m.mat');

% % % % % take only signif cells
% % % % signif = cat(1,cells.signif);
% % % % signif = arrayfun(@(x)(x.TF),signif);
% % % % signif = any(signif,2);
% % % % cells_all = cells(signif);
% % % % clear signif cells
        
% take all cells (later choose signif)
cells_all = cells;
clear cells
        
%%
details = [cells_all.details];
cells_bat_num = [details.bat];
bats = unique(cells_bat_num);
n_str_pos_x = [1 1 1.2; 1 1 1.1; 1 1 1.1; 1 1 1.1; 1 1 1.1];
n_str_pos_y = [1 1 1; 0.9 0.9 0.9; 0.9 0.9 0.9; 0.9 0.9 0.9; 0.9 0.9 0.9];
for ii_bat = 1:length(bats)
    
    %% arrange bat data
    bat = bats(ii_bat);
    cells_bat = cells_all(cells_bat_num == bat);
    signif_bat = cat(1,cells_bat.signif);
    signif = cat(1,cells_bat.signif);
    signif = arrayfun(@(x)(x.TF),signif);
    signif = any(signif,2);
    cells = cells_bat(signif);
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
    
    %% report per bat info
    fprintf('bat%d %d:\n',ii_bat,bat);
    fprintf('\tTotal cells\t\t\t\t\t %d\n', length(cells_bat));
    valid_flights = arrayfun(@(x)(x.has_min_flights),signif_bat);
    valid_spikes = arrayfun(@(x)(x.has_min_flights),signif_bat);
    valid_flights_spikes = valid_flights & valid_spikes;
    fprintf('\t\tvalid flights\t\t\t %d\n', sum(any(valid_flights,2)) );
    fprintf('\t\tvalid flights+spikes\t %d\n', sum(any(valid_flights_spikes,2)) );
    fprintf('\tsignif place cells count\t %d\n', length(cells));
    details_bat=[cells_bat.details];
    details=[cells.details];
    fprintf('\tsessions with CA1 cells\t\t %d\n', length(unique({details_bat.exp_ID})));
    fprintf('\t\t valid flights\t\t\t %d\n', length(unique({details_bat(any(valid_flights,2)).exp_ID})));
    fprintf('\t\t valid flights+spikes\t %d\n', length(unique({details_bat(any(valid_flights_spikes,2)).exp_ID})));
    fprintf('\tsessions with signif PC\t\t %d\n', length(unique({details.exp_ID})));
    exp_list = unique({details.exp_ID})';
    fprintf('\tfirst exp with PC\t\t\t %s\n',exp_list{1});

    %% no. of fields
    x = field_num;
    axes(panel_ABC(ii_bat,1));  cla; hold on
    hh=histogram(x);
    hh.FaceColor = 0.5*[1 1 1];
    hh.BinEdges = 0.5+[0:20];
    hh.Data(hh.Data > hh.BinLimits(2)) = hh.BinLimits(2);

    hax=gca;
    hax.YScale = 'log';
    hax.YLim(1) = 0.7;
    hax.YLim(2) = 1.15 * hax.YLim(2);
%     ha.YLim = [0.7e0 240];
    hax.XLim = [0 hh.BinLimits(2)];
    hax.XTick = [0:5:30];
    hax.YTick = [1 10 100];
    hax.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
    hax.TickDir='out';
    hax.TickLength = [0.03 0.03];
    hax.XRuler.TickLabelGapMultiplier = -0.3;
    hax.YRuler.TickLabelGapMultiplier = 0.001;
    if ii_bat == length(bats)
        xlabel({'No. of fields per direction'},'Units','normalized','Position',[0.5 -0.18]);
    end
    ylabel('No. of cells','Units','normalized','Position',[-0.28 0.5])
    n_str = "n = "+sum(~isnan(x));
    if ii_bat == 1
        n_str = n_str + " cells";
    end
    text(n_str_pos_x(ii_bat,1),n_str_pos_y(ii_bat,1), n_str,...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
    text(-0.6, 0.5, "Bat "+bat_ID_num_map(bat),...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',8,'FontWeight','bold');

    hl=xline(nanmean(x)); hl.Color='r';
    m = hax.YLim(2) + 0.15*range(hax.YLim);
    plot(prctile(x,[25 75]), [m m], 'r-','LineWidth',1   ,'Clipping','off');
    plot(prctile(x,[50]),    m    , 'r.','MarkerSize',10 ,'Clipping','off');
    
    %% Field Size
    x = field_size;
    axes(panel_ABC(ii_bat,2));  cla; hold on
    h=histogram(x);
    h.FaceColor = 0.5*[1 1 1];
    h.BinEdges = 0:33;
    
    hax=gca;
    hax.YScale = 'log';
    hax.YLim(1) = 0.8;
    hax.YLim(2) = 1.15 * hax.YLim(2);
    hax.YTick = [1 10 100];
    hax.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
    hax.TickDir='out';
    hax.TickLength = [0.03 0.03];
    hax.XRuler.TickLabelGapMultiplier = -0.3;
    hax.YRuler.TickLabelGapMultiplier = 0.001;
    if ii_bat == length(bats)
        xlabel('Field size (m)')
    end
    ylabel('No. of fields','Units','normalized','Position',[-0.28 0.5])
    n_str = "n = "+sum(~isnan(x));
    if ii_bat == 1
        n_str = n_str + " fields";
    end
    text(n_str_pos_x(ii_bat,2),n_str_pos_y(ii_bat,2), n_str,...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
    
    hl=xline(nanmean(x)); hl.Color='r';
    m = hax.YLim(2) + 0.15*range(hax.YLim);
    plot(prctile(x,[25 75]), [m m], 'r-','LineWidth',1   ,'Clipping','off');
    plot(prctile(x,[50]),    m    , 'r.','MarkerSize',10 ,'Clipping','off');
        
    %% ratio L/S
    x = ratio_LS;
    axes(panel_ABC(ii_bat,3)); cla; hold on
    nBinEdges = 9;
    edges = logspace(0,log10(25),nBinEdges);
    h=histogram(x);
    h.BinEdges = edges;
    h.FaceColor = 0.5*[1 1 1];
    h.FaceAlpha = 1;    
    
    hax=gca;
    hax.XScale = 'log';
    hax.YScale = 'log';
    hax.YLim(1) = 0.7;
    hax.YLim(2) = 1.15 * hax.YLim(2);
    hax.XLim = [0 27];
    hax.YTick = [1 10 100];
    hax.XTick = [1 2 5 10 20];
    hax.YTickLabel = {'10 ^0';'10 ^1';'10 ^2'};
    hax.TickDir='out';
    hax.TickLength = [0.03 0.03];
    hax.XRuler.TickLabelGapMultiplier = -0.35;
    hax.YRuler.TickLabelGapMultiplier = 0.001;
    if ii_bat == length(bats)
        xlabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[0.5 -0.17]);
    end
    ylabel('No. of cells','Units','normalized','Position',[-0.28 0.5])
    n_str = "n = "+length(x);
    if ii_bat == 1
        n_str = n_str + " cells";
    end
    text(n_str_pos_x(ii_bat,3),n_str_pos_y(ii_bat,3), n_str,...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',7);
    
    hl=xline(nanmean(x)); hl.Color='r';
    m = hax.YLim(2) + 0.15*range(hax.YLim);
    plot(prctile(x,[25 75]), [m m], 'r-','LineWidth',1   ,'Clipping','off');
    plot(prctile(x,[50]),    m    , 'r.','MarkerSize',10 ,'Clipping','off');
    
end

% add panel letters
axes(panel_ABC(1,1));
text(-0.3,1.2, 'A', 'Units','normalized','FontWeight','bold');
axes(panel_ABC(1,2));
text(-0.3,1.2, 'B', 'Units','normalized','FontWeight','bold');
axes(panel_ABC(1,3));
text(-0.3,1.2, 'C', 'Units','normalized','FontWeight','bold');

%% print/save the figure
file_out = fig_name_str;
file_out = fullfile(res_dir, file_out);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');

%%
