function exp_explore_units_FR(NTT_file)

%% check input
if nargin==0
    [NTT_filename NTT_path] = uigetfile(fullfile('L:\Analysis\pre_proc','*.NTT'), 'choose NTT file');
    NTT_file = fullfile(NTT_path, NTT_filename);
end

%% get exp_ID from NTT filename
[~, name, ~] = fileparts(NTT_file);
exp_ID = regexp(name, '(b[\d]+_d[\d]+)','tokens');
exp_ID = exp_ID{1}{1};

%% load exp
exp = exp_load_data(exp_ID,'details','path', 'pos', 'flight');

%% read units from NTT file
[Timestamps, CellNumbers, Samples, Header] = ...
     Nlx2MatSpike(NTT_file, [1 0 1 0 1], 1, 1, [] );
% parse header and convert bits to uVolts
ADBitVolts = sscanf(Header{contains(Header,'ADBitVolts')},'-ADBitVolts %f %f %f %f');
Samples = Samples .* ADBitVolts' .* 1e6;
CellNumList = unique(CellNumbers);
CellNumList = setdiff(CellNumList, 0); % remove the un-sorted unit
nCells = length(CellNumList);

%% create FE struct for each unit
cells = {};
for ii_cell = 1:nCells
    %% get cell spikes data
    cell = struct();
    cellNum = CellNumList(ii_cell);
    IX = find(CellNumbers==cellNum);
    cell.spikes.ts = Timestamps(IX);
    cell.spikes.waveforms = Samples(:,:,IX);
    cell.spikes.pos = interp1(exp.pos.proc_1D.ts, exp.pos.proc_1D.pos,       cell.spikes.ts, 'linear');
    cell.spikes.vel = interp1(exp.pos.proc_1D.ts, exp.pos.proc_1D.vel_csaps, cell.spikes.ts, 'linear');
    cell.cellNum = cellNum;
    
    %% get relevant position data
    FE_by_dir = {};
    directions = [1 -1];
    for ii_dir = 1:length(directions)
        %% get FE struct
        FE_dir_IX = find([exp.flight.FE.direction]==directions(ii_dir));
        FE = exp.flight.FE(FE_dir_IX);
        % remove flight epochs without spikes
        FE_remove_1 = find([FE.end_ts] < cell.spikes.ts(1));
        FE_remove_2 = find([FE.start_ts] > cell.spikes.ts(end));
        FE([FE_remove_1 FE_remove_2]) = [];

        %% add spikes data to FE struct
        ti = [FE.start_ts; FE.end_ts]';
        [~, spikes_IX_per_ti] = get_data_in_ti(cell.spikes.ts, ti);
        spikes_ts_by_epoch    = cellfun(@(x)(cell.spikes.ts(x)),            spikes_IX_per_ti, 'UniformOutput',false);
        spikes_wvfrm_by_epoch = cellfun(@(x)(cell.spikes.waveforms(:,:,x)), spikes_IX_per_ti, 'UniformOutput',false);
        spikes_pos_by_epoch   = cellfun(@(x)(cell.spikes.pos(x)),           spikes_IX_per_ti, 'UniformOutput',false);
        spikes_vel_by_epoch   = cellfun(@(x)(cell.spikes.vel(x)),           spikes_IX_per_ti, 'UniformOutput',false);
        num_spikes_by_epoch = cellfun(@length, spikes_ts_by_epoch);
        [FE(:).spikes_ts] = disperse(spikes_ts_by_epoch);
        [FE(:).spikes_pos] = disperse(spikes_pos_by_epoch);
        [FE(:).spikes_vel] = disperse(spikes_vel_by_epoch);
        [FE(:).spikes_wvfrm] = disperse(spikes_wvfrm_by_epoch);
        [FE(:).num_spikes] = disperse(num_spikes_by_epoch);
        FE_by_dir{ii_dir} = FE;
    end
    cell.FE = FE_by_dir;
    cells{ii_cell} = cell;
end
cells = [cells{:}];

%% calc FR maps for all units
for ii_cell = 1:nCells
    for ii_dir = 1:2
        FE = cells(ii_cell).FE{ii_dir};
        cells(ii_cell).FR_map(ii_dir) = FE_compute_PSTH(FE);
    end
end

%% calc correlations between all maps (of different units)
% first, arrange all maps together in one matrix
cells.FR_map;
FR_maps_all = arrayfun(@(x)([x.FR_map(1).PSTH x.FR_map(2).PSTH]), cells, 'UniformOutput', 0);
FR_maps_all = cat(1,FR_maps_all{:});
[FR_maps_corr,~] = corr(FR_maps_all', 'rows','pairwise');

%% create figures output dir
[pathstr, name, ext] = fileparts(NTT_file);
dir_out = fullfile(pathstr,name);
mkdir(dir_out);

%% fig 1 - spikes and FR maps from all units together!
figure('Units','normalized','Position',[0 0 1 1]);
pnl = panel();
pnl.pack('h',[50 50],'v',[35 65])
pnl.margin = [15 25 15 15];
pnl.title(NTT_file)

ti = exp_get_sessions_ti(exp_ID, 'Behave');
directions = [1 -1];
directions_str = {'--->>>';'<<<---'};
% colors from https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
units_colors = [...
[230, 25, 75];...
[60, 180, 75];...
[255, 225, 25];...
[0, 130, 200];...
[245, 130, 48];...
[145, 30, 180];...
[70, 240, 240];...
[240, 50, 230];...
[210, 245, 60];...
[250, 190, 190];...
[0, 128, 128];...
[230, 190, 255];...
[170, 110, 40];...
[255, 250, 200];...
[128, 0, 0];...
[170, 255, 195];...
[128, 128, 0];...
[255, 215, 180];...


[210, 245, 60];...
[250, 190, 190];...
[0, 128, 128];...
[230, 190, 255];...
[170, 110, 40];...
[255, 250, 200];...
[128, 0, 0];...
[170, 255, 195];...
[128, 128, 0];...
[255, 215, 180];...
];
units_colors = units_colors ./ 255;
% units_colors = repelem(units_colors,10,1);
obj_handles_units = {};
obj_handles_pos = [];
for ii_dir = 1:2
    % FR maps
    pnl(ii_dir,1).select();hold on
    set(gca,'Color','k');
    for ii_cell = 1:length(cells)
        cell = cells(ii_cell);
        h = plot(cell.FR_map(ii_dir).bin_centers, cell.FR_map(ii_dir).PSTH,...
            'Color', units_colors(ii_cell,:), 'LineWidth',1.5);
        obj_handles_units{1,ii_dir,ii_cell} = h;
        if ~isempty(h)
            h.DisplayName = sprintf('unit %d (%s)', cell.cellNum, 'a'+cell.cellNum-1);
        end
    end
    xlabel('Position (m)')
    ylabel('Firing rate (Hz)')
    h=title(directions_str{ii_dir}); h.FontSize = 12;
    
    % trajectory + spikes
    pnl(ii_dir,2).select();hold on
    set(gca,'Color','k');
    exp_FE = exp.flight.FE([exp.flight.FE.direction]==directions(ii_dir));
    obj_handles_pos(ii_dir) = plot( [exp_FE.pos], [exp_FE.ts], '.', 'Color', 0.9.*[1 1 1]);
    for ii_cell = 1:length(cells)
        cell = cells(ii_cell);
        cell_FE = cell.FE{ii_dir};
        h = plot([cell_FE.spikes_pos],[cell_FE.spikes_ts], '.', 'Color', units_colors(ii_cell,:));
        obj_handles_units{2,ii_dir,ii_cell} = h;
        if ~isempty(h)
            h.DisplayName = sprintf('unit %d (%s)', cell.cellNum, 'a'+cell.cellNum-1);
        end
    end
    ylim(ti)
    rescale_plot_data('y',[1e-6/60 ti(1)])
    xlabel('Position (m)')
    ylabel('Time (min)')
    h=title(directions_str{ii_dir}); h.FontSize = 12;
end

% link xlim
linkaxes(pnl.de.axis, 'x');

% link visability property
plotbrowser('on')
Link = {};
for ii_cell = 1:size(obj_handles_units,3)
    handles = [obj_handles_units{:,:,ii_cell}];
    Link{ii_cell} = linkprop(handles, 'Visible');
end
setappdata(gcf, 'LinkUnits', Link);
Link = linkprop(obj_handles_pos, 'Visible');
setappdata(gcf, 'LinkPosition', Link);

% save figure
file_out = fullfile(dir_out, 'maps_all_units');
saveas(gcf, file_out, 'fig')
saveas(gcf, file_out, 'tif')

%% fig 2 - correlation matrix (of FR maps)
figure
% replace main diagonal with nans
M = FR_maps_corr;
mask = eye(size(M));
mask(mask==1) = nan;
mask(mask==0) = 1;
M = M.*mask;
h=imagesc(M);
set(h,'AlphaData',~isnan(M));
colorbar
% colormap jet
title('maps correlation (between all units)')
xlabel('units')
ylabel('units')
set(gca,'xtick',CellNumList, 'ytick',CellNumList);
file_out = fullfile(dir_out, 'maps_corr');
saveas(gcf, file_out, 'tif')

%%
% close all

end






