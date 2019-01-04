function cell_plot_time_AC(cell_ID)

%% load cell/exp data
cell = cell_load_data(cell_ID,'details','time_AC');
% exp = exp_load_data(cell.details.exp_ID,'details');
prm = PARAMS_GetAll();
dir_colors = prm.graphics.colors.flight_directions;

%% create figure
fig_size = [30 20];
figure('Units','centimeters','Position',[5 5 fig_size], 'PaperPosition',[0 0 fig_size]);
pnl = panel();
pnl.pack('v',2);

%% plot (flight + by directions)
pnl(1).select(); hold on
plot_AC(cell.time_AC.in_flight, {'color','k', 'LineWidth',1.5});
for ii_dir = 1:2
    c = dir_colors{ii_dir};
    plot_AC(cell.time_AC.by_dir(ii_dir), {'color',c, 'LineWidth',1});
end
legend({'in-flight';'dir1';'dir2'})

%% plot (flight by directions)
pnl(2).select(); hold on
if ~(isempty(cell.time_AC.in_field) | isempty(cell.time_AC.by_fields) )
    plot_AC(cell.time_AC.in_field, {'color','k', 'LineWidth',1.5});
    fields_color = get(groot,'defaultAxesColorOrder');
    for ii_field = 1:length(cell.time_AC.by_fields)
        c_IX = mod( ii_field-1, size(fields_color,1)-1 )+1;
        c = fields_color(c_IX,:);
        plot_AC(cell.time_AC.by_fields(ii_field) , {'color',c, 'LineWidth',0.1});
    end
    legend_str = cellfun(@(x)(sprintf('#%d',x)),num2cell(1:length(cell.time_AC.by_fields)),'UniformOutput',false);
    legend_str = {'all' legend_str{:}};
    leg_max_row = 8;
    leg_ncol = ceil(length(legend_str)/leg_max_row);
    hleg = columnlegend(leg_ncol, legend_str, 'location', 'NorthEast','boxoff');
    hax = gca;
    hleg.Position([1:2]) = hax.Position([1:2])+[0.8 0.2];
    hleg.Position([3:4]) = [0.03 0.2];
end

end


function plot_AC(AC, line_prop)
    plot(AC.lags, AC.c/max(AC.c), line_prop{:});
end