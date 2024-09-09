%%
clear
clc

%% create all figures!
% main
paper_replay_fig_1          % fig 1 - replay examples
paper_replay_fig_2_new      % fig 2 - single units activity in replay
paper_replay_fig_2          % fig 3 - population
                            % fig 4 - shir (short tunnel)
paper_replay_fig_3          % fig 5 - 2 bats (behavioral relevance)

% supp
paper_replay_fig_supp_1     % EDF1 - many examples
paper_replay_fig_supp_2     % EDF2 - replay decoding deconstructed
paper_replay_fig_supp_3     % EDF3 - replay vs flight speed + stats
paper_replay_fig_supp_4     % EDF4 - per bat results
                            % EDF5 - shir long replays
paper_replay_fig_supp_5     % EDF6 - novelty from day 1
paper_replay_fig_supp_MUA_FR_maps % EDF7
paper_replay_fig_supp_6     % EDF8 - reverse/forward + future+past + takeoff/landing/midair + ball1/ball2

% paper_replay_fig_supp_replay_directionality % EDF9 --> move to main fig 5
paper_replay_fig_supp_7     % EDF9 - 2bats