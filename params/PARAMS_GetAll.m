function p = PARAMS_GetAll()

% units are in (unless otherwise indicated):
% - meters for position
% - meters/second for velocity
% - uV for voltage
% radians for phase

p.arena.balls_detect_area = [6 7; 191 192];
p.arena.balls_detect_vel_thr = 0.1;

% p.pos.pos_smooth
% p.pos.pos_max_interp
% p.pos.vel_min
% p.pos.vel_max
% p.pos.resample_fs = 100;

p.flight.vel_min = 1;
p.rest.vel_thr = 0.1;
p.rest.duration_min = 2;
p.rest.merge_thr = 0.5;
p.rest.balls_margins = 0.2; % max distance to balls

% p.spikes.thr_std
% p.spikes.thr_uV
% p.spikes.lockout_time_usec;
% p.spikes.noise_coincidence_win_msec

p.oscillations.time_AC.bin_size = 0.010;
p.oscillations.time_AC.win_size = 2;
p.oscillations.phase_AC.bin_size = 2*pi/15;
p.oscillations.phase_AC.win_size = 4 * 2*pi;
p.oscillations.phase_hist.bin_size = 2*pi/20;
p.oscillations.STA.win_size = 2;

p.LFP.ripples.freq_band = [120 180];
p.LFP.ripples.thr_std = 4;
p.LFP.theta.band = [5 8];
p.LFP.theta.resample_fs = 1000;
p.LFP.flight_rhythm.band = [3 10];
p.LFP.flight_rhythm.resample_fs = 1000;

% p.FR_map.map_1D.bin_size = 0.1;
% p.FR_map.map_1D.ker_SD = 0.2;
% p.FR_map.map_1D.ker_type = 'gaussian';

p.Ipos.time_bin_size = 0.050;
p.Ipos.pos_bin_size = 1;
p.Ipos.pos_bin_limits = [0 200];

p.graphics.colors.flight_directions = {...
    [0         0.4470    0.7410];...
    [0.8500    0.3250    0.0980]};

end

