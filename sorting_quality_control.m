%%
cells = cellfun(@(c)(cell_load_data(c,'details','stats','cluster_quality')), {cells_details.cell_ID}, 'UniformOutput',0);
cells = [cells{:}];
CQ = [cells.cluster_quality];
stats = [cells.stats];
stats_all = [stats.all];

%%
figure('Units','normalized','Position',[0 0 1 1])
pnl = panel();
pnl.pack(2,3);
pnl.margin = 25;
h=pnl.title('Sorting quality control for field ratio')
h.FontSize = 15;
h.Position = [0.5 1.05];

pnl(1,1).select();
hold on
histogram([CQ.Isolation_dis_median],[0:5:100 inf])
histogram([CQ.Isolation_dis_mean],[0:5:100 inf])
xlabel('Isolation index')
ylabel('Counts')
legend({'mean';'median'})
title('by segments')

% subplot(232)
pnl(2,1).select();
hold on
plot([stats_all.IsoDist], [CQ.Isolation_dis_mean], '.')
% plot([stats_all.IsoDist], [CQ.Isolation_dis_median], '.')
refline(1,0)
axis equal
xlim([0 100])
ylim([0 100])
xlabel('Isolation index (all)')
ylabel('Isolation index (by segments)')

% subplot(233)
pnl(1,2).select();
hold on
plot([stats_all.IsoDist], [stats_all.field_ratio_LS], '.k')
% xlim([100 200])
xlabel('Isolation index (all)')
ylabel({'Field ratio';'largest/smallest'})
title('pooled all recordings')

% subplot(234)
pnl(2,2).select();
hold on
plot([CQ.Isolation_dis_mean], [stats_all.field_ratio_LS], '.')
plot([CQ.Isolation_dis_median], [stats_all.field_ratio_LS], '.')
% xlim([0 100])
xlabel('Isolation index (all)')
ylabel({'Field ratio';'largest/smallest'})
legend({'mean';'median'})
title('by segments')

% subplot(235)
pnl(1,3).select();
hold on
plot([stats_all.IsoDist], [stats_all.field_ratio_LS], '.k')
xlim([0 100])
xlabel('Isolation index (all)')
ylabel({'Field ratio';'largest/smallest'})
title('pooled all recordings')

% subplot(236)
pnl(2,3).select();
hold on
plot([CQ.Isolation_dis_mean], [stats_all.field_ratio_LS], '.')
plot([CQ.Isolation_dis_median], [stats_all.field_ratio_LS], '.')
xlim([0 100])
xlabel('Isolation index (all)')
ylabel({'Field ratio';'largest/smallest'})
legend({'mean';'median'})
title('by segments')

