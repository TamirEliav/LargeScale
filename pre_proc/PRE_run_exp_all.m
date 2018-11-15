% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Run all exp from excel
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%%
clear all
clc

%% load exp data (excel file)
exp_data_file = 'L:\Analysis\Code\inclusion_lists\recordings_summary.xlsx';
exp_t = readtable(exp_data_file, 'Sheet', 'Experiments', 'ReadRowNames',1)

%% choose subset of exp
exp_t = exp_t(exp_t.batNum==148, :)

%% run over the list
for ii_exp = 1:size(exp_t,1)
    exp_ID = exp_t.Properties.RowNames{ii_exp}
%     PRE_run_exp( exp_ID );
    PRE_calc_flight_rhythm( exp_ID );
end
