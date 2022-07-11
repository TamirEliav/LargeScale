function exp_detect_rest(exp_ID)

%% load exp data
exp = exp_load_data(exp_ID,'pos','LM');
prm = PARAMS_GetAll();
dir_out = 'L:\Analysis\Results\exp\rest';

%% arrange relevant data
session_ti = exp_get_sessions_ti(exp_ID,'Behave','Behave_6m','Light1','Light2');
% session_ti = exp_get_sessions_ti(exp_ID,'Behave_6m','Light1','Light2');
session_ti(any(isnan(session_ti),2),:) = []; % remove nan in case of missing sessions
IX = get_data_in_ti(exp.pos.proc_1D.ts, session_ti);
pos.fs = exp.pos.proc_1D.fs;
pos.ts = exp.pos.proc_1D.ts(IX);
pos.pos_original = exp.pos.proc_1D.pos(IX);
pos_smoothing_opt = 5;
csaps_p = 1e-10;
% csaps_p = 1;
switch pos_smoothing_opt
    case 1 % raw
        pos.pos = exp.pos.proc_1D.pos(IX);
        pos.vel = exp.pos.proc_1D.vel(IX);
    case 2 % original csaps used for all analysis
        pos.pos = exp.pos.proc_1D.pos_csaps(IX);
        pos.vel = exp.pos.proc_1D.vel_csaps(IX);
    case 3 % different csaps params
        pos.proc_1D.vel = [0 diff(exp.pos.proc_1D.pos).*prm.pos.resample_fs];
        pos_csaps_pp = csaps(1:length(exp.pos.proc_1D.ts), exp.pos.proc_1D.pos', csaps_p);
        pos.pos = fnval(pos_csaps_pp, 1:length(exp.pos.proc_1D.ts) );
        pos.vel = fnval( fnder(pos_csaps_pp), 1:length(exp.pos.proc_1D.ts)) .* prm.pos.resample_fs;
        pos.pos(isnan(exp.pos.proc_1D.pos)) = nan;
        pos.vel(isnan(exp.pos.proc_1D.pos)) = nan;
        pos.pos = pos.pos(IX);
        pos.vel = pos.vel(IX);
    case 4 % low pass
        t = exp.pos.proc_1D.ts;
        x = exp.pos.proc_1D.pos;
        nanx = isnan(x);
        x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));
        cutoff_freq = 1e-10;
        order = 1000;
        Wn   = cutoff_freq ./ (exp.pos.proc_1D.fs/2);
        b = fir1(order, Wn,'low');
        a = 1;
        xfilt = filtfilt(b,a,x);
        figure
        hold on
        plot(t,x,'.')
        plot(t,xfilt,'.')
%         xlim([59142191680.91274261474609375, 59279743710.0490570068359375])
%         ylim([3 10])
        pos.pos = xfilt(IX);
        pos.vel = exp.pos.proc_1D.vel(IX);
    case 5 
        t = exp.pos.proc_1D.ts;
        x = exp.pos.proc_1D.pos;
        nanx = isnan(x);
        x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));
%         csaps_p = 1e-5;
        win_s = 1;
        win = round(win_s * exp.pos.proc_1D.fs);
        xfilt = smoothdata(x,2, "movmedian",win);
        vel = [0 diff(xfilt).*prm.pos.resample_fs];
%         pos_csaps_pp = csaps(1:length(t), xfilt', csaps_p);
%         xfilt = fnval(pos_csaps_pp, 1:length(exp.pos.proc_1D.ts) );
%         vel = fnval( fnder(pos_csaps_pp), 1:length(exp.pos.proc_1D.ts)) .* prm.pos.resample_fs;
%         win_s = .1;
%         win = round(win_s * exp.pos.proc_1D.fs);
%         vel = smoothdata(vel,2, "movmedian",win);

        pos.pos = xfilt(IX);
        pos.vel = vel(IX);
%         figure
%         hold on
%         plot(t,x,'.')
%         plot(t,xfilt,'.')
%         xlim([59142191680.91274261474609375, 59279743710.0490570068359375])
%         ylim([3 10])
end
nanpos = isnan(pos.pos);
nanvel = isnan(pos.vel);
pos.pos(nanpos) = interp1(pos.ts(~nanpos), pos.pos(~nanpos), pos.ts(nanpos));
pos.vel(nanvel) = interp1(pos.ts(~nanvel), pos.vel(~nanvel), pos.ts(nanvel));
% sum(nanpos)
% sum(nanvel)

%%
% figure;
% plot(pos.pos, pos.vel,'.-');
% % xlim([4 10])
% title([pos_smoothing_opt;csaps_p])


%% params
prm.balls.balls_loc_detect_speed_thr = 0.2; % m/s
prm.balls.rest_speed_thr = 1; % m/s
prm.balls.rest_speed_thr_edges = 0.1; % m/s
prm.balls.rest_speed_thr_edges_duration = 0.1; % s
prm.balls.min_duration = 1; % in seconds
prm.balls.merge_thr = 1; % in seconds
prm.balls.max_dist_from_ball = 2; % in meters. TODO: fit some curve to the detected rest epochs and realign the ball locations PER rest epoch!

%% extract balls exact location from the data
IX = abs(pos.vel) < prm.balls.balls_loc_detect_speed_thr;
pos_on_balls = pos.pos(IX)';
% time_on_balls = pos.ts(IX)';
% vel_on_balls = pos.vel(IX)';
num_balls = 2;
cidx = kmeans(pos_on_balls, num_balls, 'Replicates',100);
balls_loc = splitapply(@median, pos_on_balls, cidx);
balls_loc = sort(balls_loc)';
[pos.dist_from_ball, pos.nearest_ball_num] = min(abs(pos.pos-balls_loc'));

%% find rest epochs
xthr_IX = find( abs(pos.vel) < prm.balls.rest_speed_thr & ...
                pos.dist_from_ball < prm.balls.max_dist_from_ball );
xthr_ts = pos.ts(xthr_IX);
start_IX = [xthr_IX(1) xthr_IX(find(diff(xthr_ts) > prm.balls.merge_thr*1e6)+1 )              ];
end_IX   = [           xthr_IX(find(diff(xthr_ts) > prm.balls.merge_thr*1e6)   )  xthr_IX(end)];
xthr = false(1,length(pos.ts));
for ii = 1:length(start_IX)
    xthr(start_IX(ii):end_IX(ii)) = true;
end
% the problem that we workaround here is that the two balls are so
% close to each other such that 2m from any ball is almost anywhere.
% Thus, we cut adjacent rest epochs from different ball. 
% This is to overcome erronous detection of rest epochs that merged two 
% epochs from different balls into one epoch.
xthr(diff(pos.nearest_ball_num)~=0) = false;
cc = bwconncomp(xthr);

% cut edges with more strict thr
xthr2 = false(size(xthr));
for ii = 1:cc.NumObjects
    IX = cc.PixelIdxList{ii};
    nPointsMargin = round(prm.balls.rest_speed_thr_edges_duration * pos.fs);
    new_start_IX = find(abs(pos.vel(IX))<prm.balls.rest_speed_thr_edges, nPointsMargin, "first");
    new_end_IX = find(abs(pos.vel(IX))<prm.balls.rest_speed_thr_edges, nPointsMargin, "last");
    if length(new_start_IX)<nPointsMargin |...
       length(new_end_IX)<nPointsMargin | ...
       new_start_IX(end)>new_end_IX(1)
        continue;
    end
    new_IX = IX(new_start_IX(end):new_end_IX(1));
    xthr2(new_IX) = true;
end
xthr = xthr2;
cc = bwconncomp(xthr);
g = bwlabel(xthr);
g(g==0)=nan;

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
invalid = ~[events.valid];
events(invalid)=[];
xthr(ismember(g,find(invalid)))=false;
cc = bwconncomp(xthr);
[events.ts] = disperse(cellfun(@(IX)(pos.ts(IX)),cc.PixelIdxList,'UniformOutput',false));
[events.pos] = disperse(cellfun(@(IX)(pos.pos(IX)),cc.PixelIdxList,'UniformOutput',false));
[events.pos_original] = disperse(cellfun(@(IX)(pos.pos_original(IX)),cc.PixelIdxList,'UniformOutput',false));
[events.vel] = disperse(cellfun(@(IX)(pos.vel(IX)),cc.PixelIdxList,'UniformOutput',false));

%%
fig=figure;
fig.WindowState = 'maximized';
hax=gca;
hax.ColorOrder = [1 0 0; 0 0 1];
hold on
plot(pos.ts,pos.pos_original,'.k')
for ii_event = 1:length(events)
    event = events(ii_event);
    plot(event.ts, event.pos_original,'.')
end
rescale_plot_data('x',[1e-6 pos.ts(1)])
xlabel('Time (s)')
ylabel('Position (m)')
title({'rest on the balls detection';exp_ID},'Interpreter','none')
file_name = fullfile(dir_out ,[exp_ID '_exp_rest_detection']);
saveas(fig,file_name,'fig')
saveas(fig,file_name,'jpg')

%% save rest on the balls results to mat file
rest = struct();
rest.params = prm.balls;
rest.events = events;
rest.ti = [events.start_ts; events.end_ts]';
rest.fs = pos.fs;
rest.balls_loc = balls_loc;
rest.ts = pos.ts;
rest.vel_smooth = pos.vel;
file_name = fullfile(dir_out ,[exp_ID '_exp_rest']);
save(file_name,'rest');

%%











%%