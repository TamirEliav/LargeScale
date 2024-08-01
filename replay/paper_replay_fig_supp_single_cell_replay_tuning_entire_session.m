%% Replay - Fig supp XXX - single cell replay tuning
clear 
clc
close all

%% data options 
cat_IX = 3; % sleep+rest pooled
ii_cat = 3;
params_opt = 11;
exp_ID = 'b0184_d191205'; % best one, same session from main fig 2

%% graphics params
lw = 1.5;
sym = '|';
sym_line_width = 1;

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Extended_Data_Fig_XXX_single cell_replay_tuning_entire_session';
fig_caption_str = 'single_cell_replay_tuning_entire_session';
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
clear panels
h = [1.5 1.5 8 8];
w = 8.5;
x = 2.5;
y = [2 4.5 7 16]+1;
panels(1,1) = axes('Units','centimeters','Position',[x y(1) w h(1)]);
panels(1,2) = axes('Units','centimeters','Position',[x y(2) w h(2)]);
panels(1,3) = axes('Units','centimeters','Position',[x y(3) w h(3)]);
panels(1,4) = axes('Units','centimeters','Position',[x y(4) w h(4)]);
x = x + w + 1;
panels(2,1) = axes('Units','centimeters','Position',[x y(1) w h(1)]);
panels(2,2) = axes('Units','centimeters','Position',[x y(2) w h(2)]);
panels(2,3) = axes('Units','centimeters','Position',[x y(3) w h(3)]);
panels(2,4) = axes('Units','centimeters','Position',[x y(4) w h(4)]);


%% load cells data
cells_t = DS_get_cells_summary();
cells_exp_ID = cellfun(@(c)DS_get_exp_ID_from_cell_ID(c),cells_t.cell_ID,'UniformOutput',false);
TF = string(cells_exp_ID) == string(exp_ID);
cells_t(~TF,:)=[];
% whos cells_t
cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',1);
details = [cells.details];
details(~contains({details.brain_area}, {'CA1','CA3'})) = [];
details(~ismember([details.ClusterQuality], [2])) = [];
cells_t = cells_t({details.cell_ID},:);
cells = cellfun(@(c)cell_load_data(c,'details','signif','inclusion'),cells_t.cell_ID);
% whos cells_t cells
clear cells_t
if isempty(cells)
    return
end
inclusion = cat(1,cells.inclusion);
cells(~[inclusion(:,1).pyr])=[];
details = [cells.details];
if isempty(details)
    return
end
cells = cellfun(@(c)cell_load_data(c,'details','signif','inclusion','fields','FR_map','replay_FR_map','spikes','FE'),{details.cell_ID});
inclusion = cat(1,cells.inclusion);
signif = cat(1,cells.signif);


%% plot entire session
for ii_dir = 1:2
    inclusion_dir = inclusion(:,ii_dir);
    signif_dir = signif(:,ii_dir);
    TF = [inclusion_dir.TF];
    TF = TF & [signif_dir.SI_thr_shuffle];
    TF = TF & [signif_dir.SI_thr_signif];
    cells_dir = cells(TF);
    nCellsDir = length(cells_dir);
    if nCellsDir==0
        continue
    end
    clrs = distinguishable_colors(nCellsDir);
    
    axes(panels(ii_dir,1));
    cla reset
    hold on
%     hax=gca; hax.ColorOrder = clrs;
    for ii_cell = 1:nCellsDir
        cell = cells_dir(ii_cell);
        x = cell.FR_map(ii_dir).all.bin_centers;
        y = cell.FR_map(ii_dir).all.PSTH;
        c = clrs(ii_cell,:);
        plot(x,y,'LineWidth',lw,'Color',c)
    end
    xlabel('Position (m)')
    if ii_dir == 1
        ylabel({'Flight';'FR map';'(Hz)'})
    end
    
    axes(panels(ii_dir,2));
    cla reset
    hold on
%     hax=gca; hax.ColorOrder = clrs;
    for ii_cell = 1:nCellsDir
        cell = cells_dir(ii_cell);
        x = cell.replay_FR_map.replay_PSTH_all(ii_cat,ii_dir).bin_centers;
        y = cell.replay_FR_map.replay_PSTH_all(ii_cat,ii_dir).PSTH;
        c = clrs(ii_cell,:);
        plot(x,y,'LineWidth',lw,'Color',c)
    end
    if ii_dir == 1
        ylabel({'Replay';'FR map';'(Hz)'})
    end
    
    axes(panels(ii_dir,3));
    cla reset
    hold on
    FE = cells_dir(1).FE{ii_dir};
    FE_start_pos = arrayfun(@(fe)fe.pos([1]),FE,'UniformOutput',true);
    FE_end_pos = arrayfun(@(fe)fe.pos([end]),FE,'UniformOutput',true);
    x = [FE_start_pos; FE_end_pos];
    y = repmat(1:length(FE),2,1);
    plot(x,y,'Color',0.5*[1 1 1]);
%     hax=gca; hax.ColorOrder = clrs;
    for ii_cell = 1:nCellsDir
        cell = cells_dir(ii_cell);
        FE = cell.FE{ii_dir};
        x = [FE.spikes_pos];
        y = arrayfun(@(ii,n)ii.*ones(1,n) , 1:length(FE), [FE.num_spikes], 'UniformOutput', false);
        y = [y{:}];
        c = clrs(ii_cell,:);
        plot(x,y,sym,'LineWidth',sym_line_width,'Color',c)
    end
    if ii_dir == 1
        ylabel('Flight no.')
    end

    axes(panels(ii_dir,4));
    cla reset
    hold on
    events = cells_dir(1).replay_FR_map.replay_PSTH_all(ii_cat,ii_dir).events;
    if ~isempty(events)
        seqs = [events.seq_model];
        x = [[seqs.start_pos]; [seqs.end_pos]];
        y = repmat(1:length(seqs),2,1);
        plot(x,y,'Color',0.5*[1 1 1]);
    %     hax=gca; hax.ColorOrder = clrs;
        for ii_cell = 1:nCellsDir
            cell = cells_dir(ii_cell);
            replay = cell.replay_FR_map.replay_PSTH_all(ii_cat,ii_dir);
            x = replay.spikes_pos;
            y = replay.spikes_replay_num;
            c = clrs(ii_cell,:);
            plot(x,y,sym,'LineWidth',lw,'Color',c)
        end
        ylim([0 length(seqs)+1])
    end
    if ii_dir == 1
        ylabel('Replay no.')
    end
    dir_arrow_str_map = containers.Map([1 2],{'\rightarrow','\leftarrow'});
    map_dir_str = sprintf('Flight direction %d %s',ii_dir,dir_arrow_str_map(ii_dir));
    title(map_dir_str)
end
linkaxes(panels,'x')

%% add panel letters
font_size = 11;
axes(panels(1,1))
text(-0.12,1.15, 'd', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels(1,2))
text(-0.12,1.15, 'c', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels(1,3))
text(-0.1,1.05, 'b', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels(1,4))
text(-0.1,1.05, 'a', 'Units','normalized','FontWeight','bold','FontSize',font_size);

%%
fig_name = sprintf('%s',fig_name_str);
file_out = fullfile(res_dir, fig_name);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
disp('figure saved!')

%%
