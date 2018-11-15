function pos_linearized = POS_calc_linearized(pos,calib)

%% for now ignore the calib and use a fix calibration (TODO: use the input calib)
% load('L:\DATA\9892_Rocky\calib\calib_tunnel.mat');
% calib = calib_tunnel;

%% calc linearized
[~,~,t] = distance2curve(calib.curvexy, pos, 'linear');
pos_linearized = t .* calib.tunnel_length;


end