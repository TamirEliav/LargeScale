function POS_pre_process(exp_ID)

%% read exp info
[exp_path exp_info] = DS_get_path(exp_ID);

%% tunnel calibration + LM positions
load(exp_path.calib_tunnel_file );
landmarks = readtable(exp_path.calib_landmarks_file);
landmarks.Y = feval(calib_tunnel.spline_fitresult,landmarks.X);
mapxy = [landmarks.X landmarks.Y];
[~,~,t] = distance2curve(calib_tunnel.curvexy,mapxy,'linear');
landmarks.pos_linearized = t .* calib_tunnel.tunnel_length;

%% load bsp data
load( fullfile(exp_path.bsp, ['bsp_pos_tag_' num2str(exp_info.bsp_tag_ID) '.mat'] ) );

%% project data to tunnel midline (calib)
mapxy = [bsp_pos.pos(:,1) bsp_pos.pos(:,2)];
[xy,distance,t] = distance2curve(calib_tunnel.curvexy,mapxy,'linear');
bsp_pos.pos_linearized = t .* calib_tunnel.tunnel_length;
bsp_pos.pos_linearized_dist = distance;
bsp_pos.pos_xy_calib = xy;

%% Use csaps to smooth the position data and derivate the velocity 
pos_csaps_smoothing = 0.0005;
fs = 1e9 / median( diff( bsp_pos.ts_ns ) );

pos_csaps_pp = csaps(1:length(bsp_pos.ts_ns), bsp_pos.pos', pos_csaps_smoothing);
pos_csaps = fnval(pos_csaps_pp, 1:length(bsp_pos.ts_ns))';
vel_csaps = fnval( fnder(pos_csaps_pp), 1:length(bsp_pos.ts_ns)) * fs';

pos_linearized_csaps_pp = csaps(1:length(bsp_pos.ts_ns), bsp_pos.pos_linearized', pos_csaps_smoothing );
pos_linearized_csaps = fnval(pos_linearized_csaps_pp, 1:length(bsp_pos.ts_ns) );
vel_linearized_csaps = fnval( fnder(pos_linearized_csaps_pp), 1:length(bsp_pos.ts_ns)) * fs;

bsp_pos.csaps_pp = pos_csaps_pp;
bsp_pos.pos_csaps_smoothing = pos_csaps_smoothing;
bsp_pos.pos_csaps = pos_csaps;
bsp_pos.vel_csaps = vel_csaps;
bsp_pos.pos_linearized_csaps_pp = pos_linearized_csaps_pp;
bsp_pos.pos_linearized_csaps = pos_linearized_csaps;
bsp_pos.vel_linearized_csaps = vel_linearized_csaps;

%% calc velocity
% TODO: make sure to consider that there are position samples that differ
% in their ISI (if there is no measurement in bsp, it simply does not
% record anything... consider adding nans in those places!!
bsp_pos.fs = 1e9/median(diff(bsp_pos.ts_ns));
bsp_pos.pos_linearized_vel = diff(bsp_pos.pos_linearized) .*1e9./diff(bsp_pos.ts_ns);
bsp_pos.pos_linearized_vel = [0; bsp_pos.pos_linearized_vel];

vel_estim_smooth_win_sec = 1.5;
vel_estim_smooth_th = 0.5;
bsp_pos.pos_linearized_vel_smooth = smooth(bsp_pos.pos_linearized_vel,bsp_pos.fs*vel_estim_smooth_win_sec);
bsp_pos.pos_linearized_dir = sign( bsp_pos.pos_linearized_vel_smooth );
IX = find( abs(bsp_pos.pos_linearized_vel_smooth ) < vel_estim_smooth_th );
bsp_pos.pos_linearized_dir( IX ) = 0;

%% TODO: add some plots (in a different function)

%% save pos data
pos = bsp_pos;
pos.exp_ID = exp_ID;
mkdir(exp_path.position);
file_name = fullfile(exp_path.position,['pos_' exp_ID]);
save(file_name,'pos');







end
