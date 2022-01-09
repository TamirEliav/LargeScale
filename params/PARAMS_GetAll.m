function prm = PARAMS_GetAll()

%% global param set definition
global paramset_global_var;
if isempty(paramset_global_var)
    paramset_global_var = 0;
end

%% units
% units are in (unless otherwise indicated):
% - meters for position
% - meters/second for velocity
% - uV for voltage
% - Hz for sample rate
% - Hz for firing rate
% - radians for phase
% - percentage in 0-100

%% params
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

prm.oscillations.time_AC.bin_size = 0.010;
prm.oscillations.time_AC.win_size = 2;
% prm.oscillations.phase_AC.bin_size = 2*pi/15;
% prm.oscillations.phase_AC.win_size = 4 * 2*pi;
% prm.oscillations.phase_hist.bin_size = 2*pi/20;
% prm.oscillations.STA.win_size = 2;

% prm.ripples.smooth_ker = 4; % ms
% prm.ripples.high_thr_std = 3; % no minimum duration
% prm.ripples.low_thr_std = 3; % for event edges
% prm.ripples.merge_thr_msec = 25;
% prm.ripples.min_width_msec = 30;
% prm.ripples.ripple_gamma_power_ratio_thr = 0; % zero means technically not applied

prm.ripples.smooth_ker = 4; % ms
prm.ripples.high_thr_std = 2.5; % no minimum duration
prm.ripples.low_thr_std = 1; % for event edges
prm.ripples.merge_thr_msec = 25;
prm.ripples.min_width_msec = 25;
prm.ripples.ripple_gamma_power_ratio_thr = 0; % zero means technically not applied

prm.MUA.bin_size = 1; % ms
prm.MUA.smooth_ker = 15; % ms
prm.MUA.high_thr_std = 2.5; % no minimum duration
prm.MUA.low_thr_std = 1;  % for event edges

prm.FR_map.bin_size = 0.2;
prm.FR_map.bin_limits = [0 200];
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
prm.fields.valid_speed_pos = [10 187.5];

prm.Ipos.time_bin_size = 0.050;
prm.Ipos.pos_bin_size = 1;
prm.Ipos.pos_bin_limits = [0 200];

prm.signif.SI_thr = 0.25;
prm.signif.SI_thr_shuffle = 99;
prm.signif.odd_even_FR_map_corr_thr = 0.5;

prm.inclusion.min_full_flights = 10; % >=
prm.inclusion.min_spikes_air   = 50; % >=
prm.inclusion.interneuron_FR_thr = 5; % 5 Hz for mean FR over entire session

prm.graphics.colors.flight_directions = {...
    [0         0.4470    0.7410];...
    [0.8500    0.3250    0.0980]};
bats_colors_mapping = {
    34,   [0.3020    0.6863    0.290];
    79,   [0.8941    0.1020    0.109];
    148,  [0.2157    0.4941    0.721];
    2289, [0.5961    0.3059    0.639];
    9861, [1.0000    0.4980    0    ];
    };
M = containers.Map([bats_colors_mapping{:,1}],...
                    bats_colors_mapping(:,2));
prm.graphics.colors.bats = M; % this is a colormap

%% decoding params
prm.decode.err_thr_prc = 5;
prm.decode.inc_criteria.err_prob = 0.4;
prm.decode.inc_criteria.max_predict_err_prob = 0.05;

%% different param sets
prm.parmaset = paramset_global_var;
switch paramset_global_var
    case 1
        prm.fields.width_href = 0.1;
    case 2
        prm.fields.width_href = 0.3;
    case 3
        prm.fields.width_prc = [0  100];
    case 4
        prm.fields.width_prc = [10  90];
    case 5
        prm.fields.width_prc = [25  75];
    case 6
        prm.fields.min_flights_with_spikes_prc = 0.25;
    case 7
        prm.fields.min_flights_with_spikes_prc = 0.30;
    case 8
        prm.fields.width_prc = [20 80];
end

%%

                
end












