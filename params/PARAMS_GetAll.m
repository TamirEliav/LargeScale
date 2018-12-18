function prm = PARAMS_GetAll()

% units are in (unless otherwise indicated):
% - meters for position
% - meters/second for velocity
% - uV for voltage
% - Hz for sample rate
% - Hz for firing rate
% - radians for phase
% - percentage in 0-100

% prm.arena.balls_detect_area = [6 7; 191 192];
% prm.arena.balls_detect_vel_thr = 0.1;

% prm.pos.pos_smooth
% prm.pos.pos_max_interp
% prm.pos.vel_min
% prm.pos.vel_max
prm.pos.outliers_pos_dist_from_midline = 2;
prm.pos.outliers_speed = 20;
prm.pos.resample_fs = 100; % note that the csaps smoothing parameter is dependent 
% on the sample rate, so if you change this, you should also change the 
% smoothing param (there is no clear linear correspondense between the fs
% and p).
prm.pos.csaps_p = 1e-5;

prm.flight.speed_high_thr = 4;
prm.flight.speed_low_thr = 1;
prm.flight.high_speed_min_duration = 1;
prm.flight.full_min_distance = 100;

% prm.rest.merge_thr = 0.5;
% prm.rest.balls_margins = 0.2; % max distance from balls

% prm.spikes.thr_std
% prm.spikes.thr_uV
% prm.spikes.lockout_time_usec;
% prm.spikes.noise_coincidence_win_msec

prm.RecStability.BinSize = 4; % in minutes

% prm.oscillations.time_AC.bin_size = 0.010;
% prm.oscillations.time_AC.win_size = 2;
% prm.oscillations.phase_AC.bin_size = 2*pi/15;
% prm.oscillations.phase_AC.win_size = 4 * 2*pi;
% prm.oscillations.phase_hist.bin_size = 2*pi/20;
% prm.oscillations.STA.win_size = 2;

% prm.LFP.ripples.freq_band = [120 180];
% prm.LFP.ripples.thr_std = 4;
% prm.LFP.theta.band = [5 8];
% prm.LFP.theta.resample_fs = 1000;
% prm.LFP.flight_rhythm.band = [3 10];
% prm.LFP.flight_rhythm.resample_fs = 1000;

prm.FR_map.bin_size = 0.2;
prm.FR_map.ker_SD = 0.5;
prm.FR_map.min_time_spent_per_meter = 0.75;
prm.FR_map.ker_type = 'gaussian';
prm.FR_map.shuffles_num = 1000;
prm.FR_map.shuffles_max_shift = 30;

prm.fields.FR_thr = 1;
prm.fields.overlap_href = 0.5; % href=horizontal reference
prm.fields.width_href = 0.2;
prm.fields.width_prc = [5 95];
prm.fields.min_spikes = 10;
prm.fields.min_flights_with_spikes = 5;
prm.fields.min_flights_with_spikes_prc = 0.2;
prm.fields.local_shuffle.margin = 0.5; % relative to field width
prm.fields.local_shuffle.n_shuffles = 1000;
prm.fields.local_shuffle.max_shift = 30;
prm.fields.local_shuffle.signif_SI_prc = 95;

prm.Ipos.time_bin_size = 0.050;
prm.Ipos.pos_bin_size = 1;
prm.Ipos.pos_bin_limits = [0 200];

prm.signif.SI_thr = 0.5;
prm.signif.SI_thr_shuffle = 99;
prm.signif.odd_even_FR_map_corr_thr = 0.5;


prm.graphics.colors.flight_directions = {...
    [0         0.4470    0.7410];...
    [0.8500    0.3250    0.0980]};

end

