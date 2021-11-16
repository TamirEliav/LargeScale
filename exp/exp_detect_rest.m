function exp_detect_rest(exp_ID)

%% load exp data
exp = exp_load_data(exp_ID,'pos','LM');
prm = PARAMS_GetAll();

%% arrange relevant data
session_ti = exp_get_sessions_ti(exp_ID,'Behave');
IX = get_data_in_ti(exp.pos.proc_1D.ts, session_ti);
pos.fs = exp.pos.proc_1D.fs;
pos.ts = exp.pos.proc_1D.ts(IX);
pos.pos = exp.pos.proc_1D.pos(IX);
% pos.pos = exp.pos.proc_1D.pos_csaps(IX);
pos.vel = exp.pos.proc_1D.vel_csaps(IX);
nanvel = isnan(pos.vel);
nanpos = isnan(pos.pos);
pos.vel(nanvel) = interp1(pos.ts(~nanvel), pos.vel(~nanvel), pos.ts(nanvel));
pos.pos(nanpos) = interp1(pos.ts(~nanpos), pos.pos(~nanpos), pos.ts(nanpos));

%% params
prm.balls.balls_loc_detect_speed_thr = 0.5;
prm.balls.rest_speed_thr = 1;
prm.balls.rest_speed_thr_edges = 0.2; % TODO: impolement!
prm.balls.min_duration = 0; % in seconds
prm.balls.merge_thr = 0.5; % in seconds
prm.balls.max_dist_from_ball = 0.5; % in meters

%% extract balls exact location from the data
IX = abs(pos.vel) < prm.balls.balls_loc_detect_speed_thr;
pos_on_balls = pos.pos(IX)';
% time_on_balls = pos.ts(IX)';
% vel_on_balls = pos.vel(IX)';
num_balls = 2;
cidx = kmeans(pos_on_balls, num_balls, 'Replicates',100);
balls_loc = splitapply(@median, pos_on_balls, cidx);
balls_loc = sort(balls_loc);
[pos.dist_from_ball, pos.nearest_ball_num] = min(abs(pos.pos-balls_loc));

%% find rest epochs
xthr_IX = find( abs(pos.vel)<prm.balls.rest_speed_thr & ...
                pos.dist_from_ball < prm.balls.max_dist_from_ball );
xthr_ts = pos.ts(xthr_IX);
start_IX = [xthr_IX(1) xthr_IX(find(diff(xthr_ts) > prm.balls.merge_thr*1e6)+1 )              ];
end_IX   = [           xthr_IX(find(diff(xthr_ts) > prm.balls.merge_thr*1e6)   )  xthr_IX(end)];
xthr = false(1,length(pos.ts));
for ii = 1:length(start_IX)
    xthr(start_IX(ii):end_IX(ii)) = true;
end
cc = bwconncomp(xthr);

% create rest events struct
events=struct();
events.duration = 1/pos.fs .* cellfun(@length,cc.PixelIdxList); % in seconds
events.start_IX = cellfun(@min, cc.PixelIdxList);
events.end_IX = cellfun(@max, cc.PixelIdxList);
events.start_ts = pos.ts(events.start_IX);
events.end_ts = pos.ts(events.end_IX);
g = bwlabel(xthr);
g(g==0)=nan;
events.mean_vel = splitapply(@mean,pos.vel,g);
events.ball_num = splitapply(@mode,pos.nearest_ball_num,g);
events.valid = events.duration > prm.balls.min_duration;
events = soa2aos(events);
events(~[events.valid])=[];
[events.ts] = disperse(cellfun(@(IX)(pos.ts(IX)),cc.PixelIdxList,'UniformOutput',false));
[events.pos] = disperse(cellfun(@(IX)(pos.pos(IX)),cc.PixelIdxList,'UniformOutput',false));
[events.vel] = disperse(cellfun(@(IX)(pos.pos(IX)),cc.PixelIdxList,'UniformOutput',false));

%%
% figure
% hold on
% plot(pos.ts(:), pos.pos(:),'.k')
% plot(pos.ts(xthr), pos.pos(xthr),'.r')
% for ii_event = 1:length(events)
%     event = events(ii_event);
%     x1 = event.start_ts;
%     x2 = event.end_ts;
%     y1 = 0;
%     y2 = 200;
%     h=fill([x1 x2 x2 x1],[y1 y1 y2 y2],'red');
%     h.FaceAlpha=0.3;
%     h.EdgeColor='none';
% end
% rescale_plot_data('x',[1e-6 pos.ts(1)])

%%
% figure
% ecdf([events.duration],'Function','survivor');
% xlim([0 10])
% hax=gca;
% hax.XScale = 'log';

%% save rest on the balls results to mat file
rest = struct();
rest.params = prm.balls;
rest.events = events;
rest.ti = [events.start_ts; events.end_ts]';
rest.fs = pos.fs;
rest.balls_loc = balls_loc;
file_name = fullfile('L:\Analysis\Results\exp\rest',[exp_ID '_exp_rest']);
save(file_name,'rest');


return


%%




















%%
figure
hold on
plot(pos.ts, pos.pos,'k.');
scatter(time_on_balls, pos_on_balls , 5, cidx, '.')
yline(balls_loc(1),'r')
yline(balls_loc(2),'r')

%%
p = 1e-100;
w = real(abs(pos.vel)<0.02);

pos.pos_csaps_balls_pp = csaps(1:length(pos.ts), pos.pos', p, [], w);
pos.pos_csaps_balls = fnval(pos.pos_csaps_balls_pp, 1:length(pos.ts) );

% pos.pos_csaps_balls_pp = csaps(pos.ts, pos.pos', p);
% pos.pos_csaps_balls = fnval(pos.pos_csaps_balls_pp, pos.ts );

pos.vel_csaps_balls = fnval( fnder(pos.pos_csaps_balls_pp), 1:length(pos.ts)) .* prm.pos.resample_fs;

%
figure
hold on
plot(pos.ts, pos.pos_csaps_balls,'.r')
plot(time_on_balls, pos_on_balls,'.k')
% ylim(balls_loc(2)+[-1 1].*0.5)
% xlim([27411154261.163448333740234375, 27446126290.86737060546875])
yline(balls_loc(1),'c')
yline(balls_loc(2),'c')

yyaxis right
plot(pos.ts, abs(pos.vel_csaps_balls),'.m')
% plot(pos.ts, abs(pos.vel),'.m')

%%
ker_STD = 1; % in seceonds
win = round(5 * ker_STD * exp.pos.proc_1D.fs);
pos.pos_smooth = smoothdata(pos.pos, 'gaussian', win);
pos.vel_smooth = [0 diff(pos.pos_smooth)];
figure
hold on
plot(pos.ts, pos.pos_smooth,'.r')
plot(time_on_balls, pos_on_balls,'.k')
yline(balls_loc(1),'c')
yline(balls_loc(2),'c')
xlim([27411154261.163448333740234375, 27446126290.86737060546875])
ylim(balls_loc(2)+[-1 1].*0.5)
yyaxis right
plot(pos.ts, abs(pos.vel_smooth),'.g')

%% 1. extract basic epochs (low speed thr crossing)
xthr_IX = find( abs(pos.vel) < prm.balls.speed_thr );
start_IX = [xthr_IX(1) xthr_IX( find(diff(xthr_IX)>1)+1 )               ];
end_IX   = [           xthr_IX( find(diff(xthr_IX)>1)   )  xthr_IX(end) ];
start_ts = pos.ts(start_IX)';
end_ts = pos.ts(end_IX)';
duration = (end_ts - start_ts) * 1e-6;
direction = sign(pos.pos(end_IX) - pos.pos(start_IX))';
distance = abs(pos.pos(end_IX) - pos.pos(start_IX))';

% create struct array of flight epochs
FE = repelem(struct,length(start_IX));
[FE(:).start_IX] = disperse(start_IX);
[FE(:).end_IX] = disperse(end_IX);
[FE(:).start_ts] = disperse(start_ts);
[FE(:).end_ts] = disperse(end_ts);
[FE(:).duration] = disperse(duration);
[FE(:).direction] = disperse(direction);
[FE(:).distance] = disperse(distance);

%% add ts/position/velocity of all samples per epoch
ti = [FE.start_ts;FE.end_ts]';
[~, IX_per_ti] = get_data_in_ti(pos.ts,ti);
ts_by_epoch = cellfun(@(x)(pos.ts(x)), IX_per_ti, 'UniformOutput',false);
pos_by_epoch = cellfun(@(x)(pos.pos(x)), IX_per_ti, 'UniformOutput',false);
vel_by_epoch = cellfun(@(x)(pos.vel(x)), IX_per_ti, 'UniformOutput',false);

[FE(:).ts] = disperse(ts_by_epoch);
[FE(:).pos] = disperse(pos_by_epoch);
[FE(:).vel] = disperse(vel_by_epoch);

%% 2. remove epochs without high speed
xthr_IX = find( abs(pos.vel) > prm.flight.speed_high_thr );
xthr_ts = pos.ts(xthr_IX);
ti = [FE.start_ts;FE.end_ts]';
[~, IX_per_ti] = get_data_in_ti(xthr_ts, ti);
[FE(:).duration_high_speed] = disperse( cellfun(@length, IX_per_ti) ./ pos.fs );
IX = find( [FE.duration_high_speed] < prm.flight.high_speed_min_duration );
FE(IX) = [];

%% TODO: remove epochs with extreme speed (or that should be already done in the position processing?!)

%% 3. divide to full/partial epochs
% TODO: decide if to do that here or when calculating the FR map...

%% 4. assign epoch numbers
[FE.number] = disperse(1:length(FE));


%% create flight struct
flight.FE = FE;
flight.speed_low_thr = prm.flight.speed_low_thr;
flight.speed_high_thr = prm.flight.speed_high_thr;
flight.high_speed_min_duration = prm.flight.high_speed_min_duration;

%% save flight struct
file_name = fullfile('L:\Analysis\Results\exp\flight',[exp_ID '_exp_flight']);
save(file_name,'flight');
    

end



