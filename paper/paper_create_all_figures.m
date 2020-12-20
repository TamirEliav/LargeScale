%%
clear
clc

%%
delay_start = 30;
delay_between_figs = 1;

%% list of scripts
fig_scripts = {
    1, 'paper_fig_1';...
    1, 'paper_fig_2';...
    1, 'paper_fig_3';...
    1, 'paper_fig_4';...
%     'paper_fig_5';... % shir
%     'paper_fig_6';... % shir
    1, 'paper_fig_7';...
    1, 'paper_fig_8';...
    
    1, 'paper_fig_anatomy';...                     % fig S1
%     '';...                                    % fig S2 (in illustrator)
    1, 'paper_fig_supp_1D_behavior';...
    1, 'paper_fig_supp_field_detection';...
    1, 'paper_fig_supp_in_field_spikes';...        % fig S5
    1, 'paper_supp_fig_compare_arms';...
    1, 'paper_fig_supp_decoding_gamma_fit';...
    1, 'paper_fig_supp_data_per_bat';...
    1, 'paper_supp_fig_compare_directions';...
    1, 'paper_supp_fig_spike_sorting_control';...  % fig S10
    1, 'paper_fig_supp_different_paramsets';...
    1, 'paper_fig_supp_compartmentalization';...
    1, 'paper_fig_supp_over_representation';...
    1, 'paper_fig_supp_inter_LM_distance';...
%     'lab_vs_wild';...                         % fig S15 (by shir)
    1, 'paper_fig_supp_decoding_variable_coverage';...
    1, 'paper_fig_supp_decoding_PV';...
    1, 'paper_fig_supp_decoding_changing_dt';...
    1, 'paper_fig_supp_attractor';...
    1, 'paper_fig_supp_mechanism';...              % fig S20
%     'field_dynamics';...                      % fig S21 (by shir)       
    };

%% run scripts!
fprintf('you have %d seconds to exit the remote control screen... go!', delay_start)
pause(delay_start)
for ii_fig = 1:size(fig_scripts,1)
    
    if ~fig_scripts{ii_fig,1}
        continue;
    end
    
    script_filename = fig_scripts{ii_fig,2};
    disp(script_filename)
    create_figure(script_filename)
    
    pause(delay_between_figs)
    close all
    pause(delay_between_figs)
end


%% create figure function
% we have this function because we use scripts to create the figures, and
% they start with clear, so we need to avoid clearing the general workspace
% (with the scripts list in it...)
function create_figure(script_filename)
    eval(script_filename)
end






%%





