%% full tunnel flights  vs short replays
close all
clear
clc

%%
out_dir = 'E:\Tamir\work\PROJECTS\LargeScale\paper_replay\figures\cell_revision';

%%
load('L:\processed_data_structs\replay_events.mat');

%%
FE_dist_norm_all = {};
nUturnsPerSession = [];
nFlightsPerSession = [];
for ii_exp = 1:length(exp_list)
    exp_ID = exp_list{ii_exp};
    exp = exp_load_data(exp_ID,'flight','uturns','rest');
    L = diff(exp.rest.balls_loc);
    FE_dist = [exp.flight.FE.distance];
    FE_dist_norm = FE_dist / L;
    FE_dist_norm(FE_dist_norm>1)=1;
    FE_dist_norm_all{ii_exp} = FE_dist_norm;
    TF = exp.uturns.pos_norm > 0.05 & exp.uturns.pos_norm < 0.95;
%     nUturnsPerSession(ii_exp) = sum(TF);
    nUturnsPerSession(ii_exp) = length(exp.uturns.pos);
    nFlightsPerSession(ii_exp) = length(FE_dist_norm);
end

%%
events = [events_all_per_session{:}];
seqs = [events.seq_model];

%%
fig1 = figure(Units="centimeters",Position=[5 5 30 20]);
tiledlayout("flow")
nexttile
bins = linspace(0,1,100);
histogram([FE_dist_norm_all{:}],bins,'normalization','probability','DisplayStyle','stairs','DisplayName','flights');
hold on
histogram([seqs.distance_norm],bins,'normalization','probability','DisplayStyle','stairs','DisplayName','replays')
legend(Location="north")
xlabel('Distance (norm.)')
ylabel('probability')
nexttile
histogram(nUturnsPerSession,10)
xlabel('No. of uturns')
ylabel('Counts (sessions)')
nexttile
histogram(nUturnsPerSession./nFlightsPerSession,10,'Normalization','probability')
xlim([0 1])
xlabel('No. of u-turns per flight')
ylabel('Fraction of sessions')
filename = fullfile(out_dir, "flights are continous but replays are fragmented.pdf");
exportgraphics(fig1,filename)

%%











%%















%%