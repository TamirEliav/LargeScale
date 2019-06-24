%% cell analysis pipline
clear
clc
%% open log file
log_name_str = ['cell_analysis ' datestr(clock, 'yyyy-mm-dd HH-MM-SS') '.txt'];
log_name_out = fullfile('L:\Analysis\Results\pipelines', log_name_str );
diary off; diary(log_name_out); diary on
% log script code
disp('-------------------------------------------------------------------')
p = mfilename('fullpath')
code_txt = fileread([p '.m'])
disp('-------------------------------------------------------------------')

%% 
cell_list = {
'b0148_d170703_TT4_SS01'...
'b0148_d170703_TT4_SS02'...
'b0148_d170703_TT4_SS03'...
'b0148_d170703_TT4_SS04'...
};

%% load cells summary and choose cells
cells_t = DS_get_cells_summary();
% cells_t(~strcmp(cells_t.brain_area, 'CA1'),:)=[];
cells_t(~ismember(cells_t.bat, [79,148,34,9861,2289] ),:) = [];
% cells_t(~ismember(cells_t.bat, [34] ),:) = [];
% cells_t(~ismember(cells_t.bat, [148] ),:) = [];
% cells_t(cells_t.date~='10/03/2018', :)=[];
% cells_t(~contains(cells_t.cell_ID, cell_list),:) = [];

%% filter cells - brain region / sorting quality
cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells = [cells.details];
cells(~contains({cells.brain_area}, 'CA1')) = [];
cells(~ismember([cells.ClusterQuality], [2])) = [];
% cells(cellfun(@isempty, {cells.stable_ts})) = [];
cells_t = cells_t({cells.cell_ID},:);

%% filter cells - mean FR
cells = cellfun(@(c)(cell_load_data(c,'stats')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells = [cells.stats];
cells = [cells.all];
cells_t([cells.meanFR_all]>5,:) = [];

%% disp final cells table
cells_t 

%% run over cells
err_list = {};
for ii_cell = 1:height(cells_t)
    %%
    cell_ID = cells_t.cell_ID{ii_cell};
    fprintf('cell %d/%d %s\n', ii_cell, height(cells_t), cell_ID);
        
    %%
try
    tic
    cell_create_details(cell_ID);
%     cell_create_spikes_data(cell_ID);
%     
%     cell_calc_time_stability(cell_ID);
%     cell_create_flight_data(cell_ID);
%     cell_calc_FR_map(cell_ID);
%     cell_calc_FR_map_shuffles(cell_ID);
%     cell_calc_Ipos(cell_ID);
% % %     cell_calc_fields(cell_ID);
% % %     cell_calc_significant(cell_ID);
% % %     cell_calc_mean_FR(cell_ID);
% % %     cell_calc_stats(cell_ID);
% % %     cell_calc_inclusion(cell_ID)
%     cell_calc_time_AC(cell_ID);
    
% % %     cell_plot_map_fields(cell_ID);
%     cell_plot_time_AC(cell_ID);
    toc
    
catch err
    fprintf('error in cell: %s\n', cell_ID);
    getReport(err)
    err_list{ii_cell}.cell_ID = cell_ID;
    err_list{ii_cell}.err = err;
end

%     pause
    close all
    
end
err_list = [err_list{:}];

%% disp all cells with errors
disp('------------------------------------------------')
disp('------------------------------------------------')
disp('------------------------------------------------')
disp('List of errors')
for ii_cell = 1:length(err_list)
    disp(err_list(ii_cell).cell_ID)
    getReport(err_list(ii_cell).err)
end
disp('------------------------------------------------')
fprintf('%d/%d cells had error!\nSee details above\n', length(err_list), height(cells_t));
if ~isempty(err_list)
    {err_list.cell_ID}'
end


%% close log file
diary off



