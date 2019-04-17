%% cell analysis pipline
clear
clc
%% open log file
log_name_str = ['cell_analysis ' datestr(clock, 'yyyy-mm-dd HH-MM-SS') '.txt'];
log_name_out = fullfile('L:\Analysis\Results\pipelines', log_name_str );
% diary off; diary(log_name_out); diary on
% TODO: save script code to log

%% 
cell_list = {
'b0148_d170711_TT4_SS01',...
'b0148_d170711_TT4_SS02',...
'b0148_d170711_TT4_SS03',...
'b0148_d170711_TT4_SS04',...
'b0148_d170711_TT4_SS05',...
'b0148_d170711_TT1_SS01',...
};

%% load cells summary and choose cells
cells_t = DS_get_cells_summary();
% cells_t(~strcmp(cells_t.brain_area, 'CA1'),:)=[];
% cells_t(~ismember(cells_t.bat, [79,148,34,9861] ),:) = [];
% cells_t(~ismember(cells_t.bat, [34] ),:) = [];
% cells_t(~ismember(cells_t.bat, [148] ),:) = [];
% cells_t(cells_t.date~='10/03/2018', :)=[];
cells_t(~contains(cells_t.cell_ID, cell_list),:) = [];
cells_t

%%
for ii_cell = 1:height(cells_t)
    %%
    cell_ID = cells_t.cell_ID{ii_cell};
    fprintf('cell %d/%d %s\n', ii_cell, height(cells_t), cell_ID);
    %%
try
    tic
    cell_create_details(cell_ID);
    cell_create_spikes_data(cell_ID);
    cell_calc_time_stability(cell_ID);
    cell_create_flight_data(cell_ID);
    cell_calc_FR_map(cell_ID);
    cell_calc_FR_map_shuffles(cell_ID);
    cell_calc_Ipos(cell_ID);
    cell_calc_fields(cell_ID);
    cell_calc_significant(cell_ID);
    cell_calc_mean_FR(cell_ID);
    cell_calc_stats(cell_ID);
    cell_plot_map_fields(cell_ID);
    
    cell_calc_time_AC(cell_ID);
    cell_plot_time_AC(cell_ID);
    toc
    
catch err
    getReport(err)
end

% pause
    close all
    
end


%% close log file
diary off



