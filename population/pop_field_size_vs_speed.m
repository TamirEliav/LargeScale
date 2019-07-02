%%
clear
clc

%%
dir_out = 'L:\Analysis\Results\population\field_size_vs_speed';
mkdir(dir_out);

%% load cells
cells_t = DS_get_cells_summary();
bats = [79,148,34,9861,2289];
cells_t(~ismember(cells_t.bat, bats ),:) = [];
prm = PARAMS_GetAll();

%% filter cells - brain region / sorting quality
cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells = [cells.details];
cells(~contains({cells.brain_area}, 'CA1')) = [];
cells(~ismember([cells.ClusterQuality], [2])) = [];
% cells(cellfun(@isempty, {cells.stable_ts})) = [];
cells_t = cells_t({cells.cell_ID},:);

%% filter cells - mean FR
cells = cellfun(@(c)(cell_load_data(c,'stats')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells = [cells.stats];
cells = [cells.all];
cells_t([cells.meanFR_all]>prm.inclusion.interneuron_FR_thr,:) = [];

%% disp final cells table
clear cells
cells_t
whos cells_t

%% load cells data
cells = cellfun(@(c)(cell_load_data(c,'fields','details','stats','signif')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells_bat = arrayfun(@(x)(x.bat), [cells.details]);

%% compare averaged field size between directions
figure('Units','normalized','Position',[0 0 1 1]);
pnl=panel();
pnl.pack('h',2);
pnl(2).pack(3,2);
pnl.margin = 20;
pnl.de.margin = 10;
h=pnl.title('comparing averaged field size between directions');
h.Position = [0.5 1.05];
h.FontSize = 12;
avg_field_size_all = [];
lm_h=[];
for ii_bat = 1:length(bats)
    bat = bats(ii_bat);
    c = prm.graphics.colors.bats(bat);
    bat_cells_TF = ismember(cells_bat,bat);
    avg_field_size_dir1 = arrayfun(@(x)(x.dir(1).field_size_mean), [cells(bat_cells_TF).stats]);
    avg_field_size_dir2 = arrayfun(@(x)(x.dir(2).field_size_mean), [cells(bat_cells_TF).stats]);
    IX1 = cellfun(@(x)(x(1).TF), {cells(bat_cells_TF).signif});
    IX2 = cellfun(@(x)(x(2).TF), {cells(bat_cells_TF).signif});
    IX = IX1 & IX2; % take significant cells in both directions!
    
    x = avg_field_size_dir1(IX)';
    y = avg_field_size_dir2(IX)';
    avg_field_size_all = [avg_field_size_all;[x y]];
    % statistics
    lm = fitlm(x,y);
    [RHO,RHO_PVAL] = corr(x,y);
    [H,P,CI,STATS] = ttest(x,y);
    
    pnl(1).select(); hold on
    lm_h(ii_bat) = plot(x, y, '.', 'Color', c);
    plot(x,predict(lm),'Color',c)
    
    [I J] = ind2sub([3,2],ii_bat);
    pnl(2,I,J).select(); hold on
    plot([x y]','.-', 'Color', c);
%     text(0.9,0.9, sprintf('%.2f',P),'Units','normalized');
    text(0.8,0.9, {'t-test';"p="+P},'Units','normalized','HorizontalAlignment','left');
    xlim([0 3]);
    set(gca,'xtick',[1 2],'xticklabel',"dir"+[1 2])
    title("bat "+bat);
end
pnl(1).select(); hold on
% lsline
lm = fitlm(avg_field_size_all(:,1), avg_field_size_all(:,2));
lm_h(end+1)=plot(avg_field_size_all(:,1),predict(lm),'Color','k');
axis equal
xlabel('Avg. field size dir1 (m)')
ylabel('Avg. field size dir2 (m)')
title('All bats')
xlimits = get(gca,'xlim');
ylimits = get(gca,'ylim');
xlimits(1) = 0;
ylimits(1) = 0;
xlimits(2) = max([xlimits ylimits]);
ylimits = xlimits;
xlim(xlimits);
ylim(ylimits);
h=refline(1,0);
h.Color = 0.8*[1 1 1];
legend(lm_h,[""+bats 'all'],'Location','northwest')
% save fig
fig_name = 'compare_directions_field_size';
fig_name = fullfile(dir_out, fig_name);
saveas(gcf, fig_name, 'tif');
saveas(gcf, fig_name, 'fig');

%% paired test figure


%% plot field size vs. velocity
figure('Units','normalized','Position',[0 0 1 1]);
pnl=panel();
pnl.pack('h',2);
pnl(2).pack(3,2);
pnl.margin = 30;
% pnl(1).margin = 45;
% pnl(2).margin = [];
pnl(2).de.margin = 15;
h=pnl.title('field size vs. velocity');
h.Position = [0.5 1.05];
h.FontSize = 12;
avg_field_size_all = [];
lm_h=[];
% dir_sym = {'>';'<'};
dir_sym = {'.';'.'};
hl = [];
fields_size_all = {};
fields_vel_all = {};
for ii_bat = 1:length(bats)
    bat = bats(ii_bat);
    bat_cells_TF = ismember(cells_bat,bat);
    bat_cells = cells(bat_cells_TF);
    for ii_dir = 1:2
        fields_size = [];
        fields_vel = [];
        for ii_cell = 1:length(bat_cells)
            cell = bat_cells(ii_cell);
            if ~isempty(cell.fields{ii_dir})
                IX = ~[cell.fields{ii_dir}.in_low_speed_area];
                cell_field_size = [cell.fields{ii_dir}.width_prc];
                cell_field_vel = [cell.fields{ii_dir}.vel];
                cell_field_size = cell_field_size(IX);
                cell_field_vel  = cell_field_vel (IX);
                fields_size = [fields_size cell_field_size];
                fields_vel = [fields_vel cell_field_vel];
            end
        end
        fields_size_all{ii_bat,ii_dir} = fields_size;
        fields_vel_all{ii_bat,ii_dir} = fields_vel;
        
        pnl(1).select(); hold on
        c = prm.graphics.colors.bats(bat);
        hl(ii_bat,ii_dir) = plot(fields_vel, fields_size, dir_sym{ii_dir}, 'Color', c);
        
        [I J] = ind2sub([3,2],ii_bat);
        pnl(2,I,J).select(); hold on
        c = prm.graphics.colors.flight_directions{ii_dir};
        h=bar(ii_dir, mean(fields_size));
        h.FaceColor = c;
        h.FaceAlpha = 0.5;
        h=errorbar(ii_dir, mean(fields_size), std(fields_size)/sqrt(length(fields_size)),'Color',c);
        h.LineWidth = 2;
    end

    [I J] = ind2sub([3,2],ii_bat);
    pnl(2,I,J).select(); hold on
    [H,P,CI,STATS] = ttest2(fields_size_all{ii_bat,1}, fields_size_all{ii_bat,2});
    text(0.05,0.95, {'t-test';"p="+P},'Units','normalized','HorizontalAlignment','left', 'VerticalAlignment','bottom');
    set(gca, 'xtick',[1 2], 'xticklabel', "dir"+[1 2]);
    xlim([0 3])
    ylabel('field size (m)')
    title("bat "+bat)
end
pnl(1).select(); hold on
title('teting effects of speed on field size')
xlabel('field velocity (m/s)')
ylabel('field size (m)')
legend(hl(:,1),""+bats,'Location','northwest')
% save fig
fig_name = 'field_size_vs_velocity_ttest';
fig_name = fullfile(dir_out, fig_name);
saveas(gcf, fig_name, 'tif');
saveas(gcf, fig_name, 'fig');


%% plot field size vs. velocity
figure('Units','normalized','Position',[0 0 1 1]);
pnl=panel();
pnl.pack('v',2);
% pnl(2).pack(3,2);
pnl(2).pack('h',length(bats));
pnl.margin = 25;
% pnl(1).margin = 45;
% pnl(2).margin = [];
pnl(2).de.margin = 15;
h=pnl.title('field size vs. velocity');
h.Position = [0.5 1.05];
h.FontSize = 12;
avg_field_size_all = [];
lm_h=[];
% dir_sym = {'>';'<'};
dir_sym = {'.';'.'};
hl = [];
fields_size_all = {};
fields_vel_all = {};
for ii_bat = 1:length(bats)
    bat = bats(ii_bat);
    bat_cells_TF = ismember(cells_bat,bat);
    bat_cells = cells(bat_cells_TF);
    for ii_dir = 1:2
        fields_size = [];
        fields_vel = [];
        for ii_cell = 1:length(bat_cells)
            cell = bat_cells(ii_cell);
            if ~isempty(cell.fields{ii_dir})
                IX = ~[cell.fields{ii_dir}.in_low_speed_area];
                cell_field_size = [cell.fields{ii_dir}.width_prc];
                cell_field_vel = [cell.fields{ii_dir}.vel];
                cell_field_size = cell_field_size(IX);
                cell_field_vel  = cell_field_vel (IX);
                fields_size = [fields_size cell_field_size];
                fields_vel = [fields_vel cell_field_vel];
            end
        end
        fields_size_all{ii_bat,ii_dir} = fields_size;
        fields_vel_all{ii_bat,ii_dir} = fields_vel;
        
        pnl(1).select(); hold on
        c = prm.graphics.colors.bats(bat);
        hl(ii_bat,ii_dir) = plot(fields_vel, fields_size, dir_sym{ii_dir}, 'Color', c);
    end

%     [I J] = ind2sub([3,2],ii_bat);
%     pnl(2,I,J).select(); hold on
    pnl(2,ii_bat).select(); hold on

    data = struct();
    data.dir1 = fields_size_all{ii_bat,1};
    data.dir2 = fields_size_all{ii_bat,2};
    rng(0);
    hv = violinplot(data);
    [hv.BoxWidth] = disperse(repelem(0.03,length(hv)));
    [hv.ViolinColor] = disperse(prm.graphics.colors.flight_directions);
    P = ranksum(data.dir1,data.dir2);
%     P = ranksum(data.dir1,data.dir2,'tail','right');
    text(0.05,0.95, {'ranksum';"p="+P},'Units','normalized','HorizontalAlignment','left', 'VerticalAlignment','bottom');
    
    ylabel('field size (m)')
    title("bat "+bat)
end
pnl(1).select(); hold on
title('teting effects of speed on field size')
xlabel('field velocity (m/s)')
ylabel('field size (m)')
legend(hl(:,1),""+bats,'Location','northwest')
% save fig
fig_name = 'field_size_vs_velocity_ranksum_two_tailed';
% fig_name = 'field_size_vs_velocity_ranksum_one_tailed';
fig_name = fullfile(dir_out, fig_name);
saveas(gcf, fig_name, 'tif');
saveas(gcf, fig_name, 'fig');

%%










%%
