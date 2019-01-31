function exp_create_position(exp_ID)

%% load exp details
exp = exp_load_data(exp_ID,'details','path');
prm = PARAMS_GetAll();

%% extract/load bespoon position data
load(exp.path.bsp_tag_file);
pos.raw = bsp_pos;
pos.raw.fs = 1e6/median(diff(pos.raw.ts_nlg_usec));

%% Linearize (project to the tunnel midline curve)
load(exp.path.calib_tunnel_file);
pos.calib_tunnel = calib_tunnel;
pos.proj.pos = POS_calc_linearized(pos.raw.pos,calib_tunnel);
pos.proj.ts = pos.raw.ts_nlg_usec;
pos.proj.speed = 1e6.*[0 ;sqrt(sum(diff(pos.proj.pos,1,1).^2,2)) ./ diff(pos.proj.ts)];

%% identify outliers
% by distance from the tunnel
pos.proj.outliers_pos_IX = find( abs(pos.proj.pos(:,2)) > prm.pos.outliers_pos_dist_from_midline );
% by high momentray velocity
pos.proj.outliers_speed_IX = find( pos.proj.speed > prm.pos.outliers_speed);
pos.proj.outliers_IX = union(pos.proj.outliers_pos_IX, pos.proj.outliers_speed_IX);
pos.proj.outliers_pos_dist_from_midline = prm.pos.outliers_pos_dist_from_midline;
pos.proj.outliers_speed = prm.pos.outliers_speed;

%% create proccessed 1D data from linearized X position
pos.proc_1D.pos = pos.proj.pos(:,1)';
pos.proc_1D.ts = pos.proj.ts';
pos.proc_1D.pos(pos.proj.outliers_IX) = [];
pos.proc_1D.ts(pos.proj.outliers_IX) = [];

%% Fill holes (interp/exterp)
[pos.proc_1D.pos, pos.proc_1D.ts] = POS_fill_holes(pos.proc_1D.pos, pos.proc_1D.ts);

%% UP-sample
Ts = 1/prm.pos.resample_fs;
ts_upsampled = pos.proc_1D.ts(1) : Ts*1e6 : pos.proc_1D.ts(end);
pos.proc_1D.pos = interp1(pos.proc_1D.ts, pos.proc_1D.pos, ts_upsampled);
pos.proc_1D.ts = ts_upsampled;
pos.proc_1D.fs = prm.pos.resample_fs;

%% Calc velocity (using csaps smoothing)
pos.proc_1D.vel = [0 diff(pos.proc_1D.pos).*prm.pos.resample_fs];
pos.proc_1D.pos_csaps_pp = csaps(1:length(pos.proc_1D.ts), pos.proc_1D.pos', prm.pos.csaps_p);
pos.proc_1D.pos_csaps_p = prm.pos.csaps_p;
pos.proc_1D.pos_csaps = fnval(pos.proc_1D.pos_csaps_pp, 1:length(pos.proc_1D.ts) );
pos.proc_1D.vel_csaps = fnval( fnder(pos.proc_1D.pos_csaps_pp), 1:length(pos.proc_1D.ts)) .* prm.pos.resample_fs;
pos.proc_1D.pos_csaps(isnan(pos.proc_1D.pos)) = nan;
pos.proc_1D.vel_csaps(isnan(pos.proc_1D.pos)) = nan;


%% save data to file
filename = fullfile('L:\Analysis\Results\exp\pos',[exp_ID '_exp_pos']);
save(filename, 'pos');


end