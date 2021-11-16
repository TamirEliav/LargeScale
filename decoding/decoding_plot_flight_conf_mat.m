function decoding_plot_flight_conf_mat(exp_ID, params_opt)

%% load data
decode_dir_IN = 'F:\sequences\decoded\flight';
out_dir = 'F:\sequences\decoded_figs\flight\conf_mat';
mkdir(out_dir)

%% load data
exp = exp_load_data(exp_ID, 'details','path','flight','pos');
decode_filename = fullfile(decode_dir_IN, exp_ID, sprintf('%s_flight_opt_%d.nc',exp_ID,params_opt) );
decode = decoding_read_decoded_file(decode_filename);

%% arrange data to compare
t = decode.time;
pos_real = interp1(exp.pos.proc_1D.ts, exp.pos.proc_1D.pos, t);
direction_real = interp1(exp.pos.proc_1D.ts, sign(exp.pos.proc_1D.vel_csaps), t);
pos_predict = decode.MAP_pos;
direction_predict = decode.MAP_direction;

%% plotting params
cmap_str = 'bone';
bin_edges = decode.pos;
bin_centers = edges2centers(bin_edges);
bin_size = median(diff(bin_edges));
pos_ticks = 0:20:200;

%% plot figure
hf = figure;
hf.WindowState = 'maximized';
pnl = panel();
pnl.pack('h',[70 30]);
pnl(1).pack(3,3);
pnl(2).pack('v',[30 30 30]);
pnl(2,1).pack('h',2);
pnl(2,3).pack('h',2);
% pnl.select('all');
% pnl.identify()
pnl.margin = 15;
pnl.de.margin = 5;
pnl(1).margin = 30;
pnl(2).de.margin = 15;
pnl(1,1,1).margin = 5;
pnl(1,1,2).margin = 5;
pnl(1,1,3).margin = 5;
pnl(1,1,1).marginleft = -10;
pnl(1,1,2).marginleft = -10;
pnl(1,1,3).marginleft = -10;
pnl(1,2,1).margin = 15;
pnl(1,2,2).margin = 15;
pnl(1,2,3).margin = 15;
pnl(1,3,1).margin = 15;
pnl(1,3,2).margin = 15;
pnl(1,3,3).margin = 15;
pnl(2,1,1).margin = 25;

directions_to_use_opts = {[-1;1],[1],[-1]};
for ii_dir_pooling_opt = 1:length(directions_to_use_opts)
    directions_to_use = directions_to_use_opts{ii_dir_pooling_opt};
    IX = any(direction_real == directions_to_use,1)';
    x = pos_real(IX);
    y = pos_predict(IX);
    N = histcounts2(x,y,bin_edges,bin_edges)';
    N_norm_by_real = N ./ sum(N,1);
    N_norm_by_predict = N ./ sum(N,2);
    
    pnl(1,1,ii_dir_pooling_opt).select();
    title("directions: "+strjoin(string(directions_to_use'),','));
    plot(x,y,'.');
    set_pos_prediction_axis(bin_edges,pos_ticks);
    
    pnl(1,2,ii_dir_pooling_opt).select();
    imagesc(bin_centers,bin_centers,N_norm_by_real);
    colormap(flip(colormap(cmap_str)));
    colorbar('Location','eastoutside')
    set_pos_prediction_axis(bin_edges,pos_ticks);
    
    pnl(1,3,ii_dir_pooling_opt).select();
    imagesc(bin_centers,bin_centers,N_norm_by_predict);
    colormap(flip(colormap(cmap_str)));
    colorbar('Location','eastoutside')
    set_pos_prediction_axis(bin_edges,pos_ticks);
end

pnl(2,1,1).select();
switch 1
    case 1
        [cm,gorder] = confusionmat(direction_real,direction_predict);
        cm = cm';
        cm = cm./sum(cm);
        cm = 100.*cm;
        h=heatmap(gorder,gorder,cm);
        h.CellLabelFormat='%0.2g%%';
        xlabel('Real direction');
        ylabel('Predicted direction');
    case 2
        h = confusionchart(direction_real,direction_predict);
        h.Normalization = 'row-normalized';
        h.Title = 'Direction confusion matrix';
        h.YLabel = 'Real direction';
        h.XLabel = 'Predicted direction';
end

pnl(2,1,2).select();
pos_err = abs(pos_real-pos_predict);
h=histogram(pos_err);
h.BinWidth = bin_size/2;
h.Normalization = 'probability';
h.LineStyle = '-';
hl=xline(median(pos_err));
hl.Color = 'r';
hl.Label = sprintf('median=%.3gm',hl.Value);
xlim([0 10]);
xlabel('Position error (m)');
ylabel('Probability');
params_str = {
    "params opt: " + params_opt;
    "pos_bin_size: "+decode.params.pos_bin_size + "m";
    "pos_std: "+decode.params.pos_std + "m";
    "mark_std: "+decode.params.mark_std + "uV";
    "state_decay_timescale: " + decode.params.state_decay_timescale + "s";
    };
text(.97,.97,params_str,...
    'horizontalalignment','right',...
    'verticalalignment','top',...
    'units','normalized',...
    'interpreter','none');

% median error given real/predicted location
pnl(2,2).select();
hold on
[~,bin_edges,BIN] = histcounts(pos_real, 'BinWidth', bin_size);
bin_centers = edges2centers(bin_edges);
G_real = findgroups(BIN);
G_pos_real = bin_centers(unique(BIN));
pos_error_median_by_real_pos = splitapply(@median, pos_err, G_real);
[~,~,BIN] = histcounts(pos_predict, bin_edges);
G_predict = findgroups(BIN);
G_pos_predict = bin_centers(unique(BIN));
pos_error_median_by_predicted_pos = splitapply(@median, pos_err, G_predict);
plot(G_pos_real, pos_error_median_by_real_pos,         '-', 'LineWidth',2);
plot(G_pos_predict, pos_error_median_by_predicted_pos, ':', 'LineWidth',2);
set(gca,'YScale','log')
ylabel('median error (m)')
xlabel('Position (m)')
legend('@real','@predicted','Location',[.9 .35 .01 .01],'units','normalized');

pnl(2,3,1).select();
x = pos_error_median_by_real_pos;
y = interp1(G_pos_predict, pos_error_median_by_predicted_pos, G_pos_real);
% plot(x,y,'.k')
lm = fitlm(x,y,'RobustOpts','on');
plot(lm);
lm_coef = lm.Coefficients;
cc=corrcoef(x,y);
% title("r="+cc(1,2)+" p="+pval(1,2));
title(sprintf('R^2=%.2g',lm.Rsquared.Ordinary));
legend off
hax=gca;
hax.TickLength(1) = 0.05;
hax.XScale = 'log';
hax.YScale = 'log';
xlabel('Median error @ real (m)')
ylabel('Median error @ predict (m)')

h=pnl.title(exp_ID);
h.Interpreter='none';
h.Position(2) = 1.04;
h.FontSize=14;

fig_filename = fullfile(out_dir, sprintf('%s_flight_decoding_opt_%d',exp_ID,params_opt));
saveas(hf, fig_filename, 'jpg');
close(hf);

%% save flight decoding results to mat file
res = struct();
res.params = decode.params;
res.params_opt = params_opt;
res.pos = decode.pos;
res.pos_real = pos_real;
res.pos_predict = pos_predict;
res.direction_real = direction_real;
res.direction_predict = direction_predict;
res.direction_confusion_matrix = cm;
res.direction_mean_acc = mean(diag(cm));
res.pos_err = pos_err;
res.pos_err_mean = mean(pos_err);
res.pos_err_median = median(pos_err);
res_filename = fullfile(out_dir, sprintf('%s_flight_decoding_opt_%d',exp_ID,params_opt));
save(res_filename, 'res');

end
%%
function set_pos_prediction_axis(bin_edges,pos_ticks)
    axis xy
    axis square
    axis equal
    xlim(bin_edges([1 end]))
    ylim(bin_edges([1 end]))
    xticks(pos_ticks)
    yticks(pos_ticks)
    xlabel('Real position (m)')
    ylabel('Predicted position (m)')
    hax=gca;
    hax.TickDir = 'out';
    hax.XRuler.TickLabelGapOffset = -1;
    hax.YRuler.TickLabelGapOffset = -1;
end

