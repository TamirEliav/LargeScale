function POS_test_join_short_flights_with_gap(exp_ID)

%%
% exp_ID = 'b0034_d180228';

%% load data
exp = exp_load_data(exp_ID, 'details', 'pos', 'flight');
prm = PARAMS_GetAll();
FE = exp.flight.FE;

%% params
gap_vel_thr = 4;
full_flight_distance = 100;

%% 
[FE.full] = disperse([FE.distance] > full_flight_distance);
% find two consequitive flights in the same direction
IX = find(diff([FE.direction])==0);
FE_candidate_merge = [FE(IX);FE(IX+1)];
% arrayfun(@(x)(x.distance),FE_same_dir);
gap_start_pos = arrayfun(@(FE)(FE.pos(end)), FE_candidate_merge(1,:));
gap_end_pos = arrayfun(@(FE)(FE.pos(1)), FE_candidate_merge(2,:));
gap_start_ts = arrayfun(@(FE)(FE.ts(end)), FE_candidate_merge(1,:));
gap_end_ts = arrayfun(@(FE)(FE.ts(1)), FE_candidate_merge(2,:));
gap_distance = gap_end_pos - gap_start_pos;
gap_time = gap_end_ts - gap_start_ts;
gap_time = gap_time .* 1e-6;
% gap_vel = 1e6 * (gap_end_pos-gap_start_pos) ./ (gap_end_ts-gap_start_ts)
gap_vel = gap_distance ./ gap_time;
valid_gaps = abs(gap_vel)>gap_vel_thr & sign(gap_vel)==[FE(IX).direction];
FE_merged = FE;
[FE_merged(IX).full]   = disperse([FE_merged(IX  ).full] | valid_gaps);
[FE_merged(IX+1).full] = disperse([FE_merged(IX+1).full] | valid_gaps);
[FE_merged.changed] = disperse([FE_merged.full] & ~[FE.full]);

%%
figure('Units','normalized','Position',[0 0 1 1]);
pnl = panel();
pnl.pack('h',[80 20]);
pnl(1).pack('h',2);
pnl(2).pack('v',4);
pnl.margin = [30 30 30 15];
pnl.de.margin = 20;
h=pnl.title(exp.details.exp_ID);
% h=suptitle(exp.details.exp_ID);
h.Interpreter = 'none';
h.FontSize = 16;
h.Position = [0.5 1.03];

pnl(1,1).select();
title('Before merging gaps')
xlabel('position')
ylabel('time')
hold on
plot(exp.pos.proc_1D.pos, exp.pos.proc_1D.ts, '.k')
FE_plot = FE;
for ii_flight = 1:length(FE_plot)
    flight = FE_plot(ii_flight);
    if flight.full
        c = 'b';
    else
        c = 'r';
    end
    plot(flight.pos, flight.ts, '.', 'Color',c)
end

pnl(1,2).select();
title('After merging gaps')
xlabel('position')
ylabel('time')
hold on
plot(exp.pos.proc_1D.pos, exp.pos.proc_1D.ts, '.k')
FE_plot = FE_merged;
for ii_flight = 1:length(FE_plot)
    flight = FE_plot(ii_flight);
    if flight.full
        c = 'b';
    else
        c = 'r';
    end
    if flight.changed
        c = 'g';
    end
    plot(flight.pos, flight.ts, '.', 'Color',c)
end

% hax=findobj(gcf,'type','axes');
linkaxes(pnl(1).de.axis,'xy')


pnl(2,1).select();
hold on
histogram( abs(gap_vel), length(gap_vel) )
xline( gap_vel_thr,'r');
% xline(-gap_vel_thr,'r');
title('newly added flight segement')
xlabel('Distance')
ylabel('Counts')

pnl(2,2).select();
histogram([FE([FE_merged.changed]).distance])
title('newly added flight segement')
xlabel('Distance')
ylabel('Counts')

pnl(2,3).select();
total_dist_before = sum([FE([FE.full]).distance]);
total_dist_after  = sum([FE_merged([FE_merged.full]).distance]);
total_dist = [total_dist_before total_dist_after] .* 1e-3;
dist_added = diff(total_dist);
dist_added_prc = 100 * diff(total_dist) / total_dist(1);
bar(total_dist);
title('Total distance')
set(gca,'Xticklabel',{'before','after'})
ylabel('Distance (km)')
h=text(1,0.9, sprintf('+%.1f%%',dist_added_prc), 'Units','normalized','HorizontalAlignment','right','Color',[0 0.5 0]);
h=text(1,0.8, sprintf('(%.1fkm)',dist_added), 'Units','normalized','HorizontalAlignment','right','Color',[0 0.5 0]);

%% save figure
fig_filename = fullfile('L:\Analysis\Results\exp\flight_merge_gap_test', [exp_ID '_exp_flight_merge_gap_test']);
saveas(gcf,fig_filename,'tif')
saveas(gcf,fig_filename,'fig')








%%
