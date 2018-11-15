%%
clear 
clc

%%
dir_IN = 'L:\Analysis\from_Shir_20181011\Wild_cells';
dir_OUT = 'L:\Analysis\from_Shir_20181011\Wild_cells_rearranged';
mkdir(dir_OUT);
files = dir([dir_IN '\spikes*bat*'])

%%
for ii_file = 1:length(files)
    
    %% load cell data
    ii_file
    filename = fullfile(files(ii_file).folder, files(ii_file).name);
    load(filename);
    cell_data = cell_struct;
    clear cell_struct
    
    %% re arrange data struct
    if isfield(cell_data.data, 'direction1')
        cell_data.data = [cell_data.data.direction1 cell_data.data.direction2];
    end
    if isfield(cell_data.fields, 'direction1')
        cell_data.fields = [cell_data.fields.direction1 cell_data.fields.direction2];
    end
%     if isfield(cell_data.shuffle, 'direction1')
%         cell_data.shuffle = [cell_data.shuffle.direction1 cell_data.shuffle.direction2];
%     end
	
    cell_data.details.bat = cell_data.bat;
    cell_data.details.day = cell_data.day;
    cell_data.details.expname = cell_data.expname;
    cell_data.details.TT = cell_data.TT;
    cell_data.details.SS = cell_data.SS;
    cell_data.details.BspTag = cell_data.general.BspTag;
    cell_data.details.path = cell_data.path;
    cell_data.details.cell_ID = sprintf('b%s_d%s_TT%s_SS%s',...
                                cell_data.bat,...
                                cell_data.day(3:end),...
                                cell_data.TT,...
                                cell_data.SS );
    
    %% save in struct format!
    fileout = fullfile(dir_OUT, files(ii_file).name);
    save(fileout, '-struct', 'cell_data');
    
end


%% rename files ...
clear 
clc
dir_IN = 'L:\Analysis\from_Shir_20181011\Wild_cells_rearranged';
files = dir([dir_IN '\spikes*bat*'])
for ii_file = 1%2:length(files)
    %% load cell data
    ii_file
    old_file_name = fullfile(files(ii_file).folder, files(ii_file).name)
    cell_data = load(old_file_name, 'details')
    new_file_name = fullfile(files(ii_file).folder, ['cell_data_' cell_data.details.cell_ID])
    movefile(old_file_name, new_file_name);
end











%%
