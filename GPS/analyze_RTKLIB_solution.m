

%%
clear
clc

%%
main_dir = 'L:\GPS\test\20180809\data';
rover_files_PREFIX = 'rover';
base_files_PREFIX = 'base';
solution_files_PREFIX = 'solutionENU';

%%
files = dir( fullfile(main_dir,['*.txt']) );
fileinfo_cell = regexpi({files.name}, '^(?<prefix>[a-z]+)_(?<date>\d+)_(?<time>\d+)', 'names','lineanchors');
fileinfo = vertcat(fileinfo_cell{:});
rover_files_IX = find(~cellfun(@isempty,strfind({fileinfo.prefix}, rover_files_PREFIX)));
base_files_IX = find(~cellfun(@isempty,strfind({fileinfo.prefix}, base_files_PREFIX)));
solution_files_IX = find(~cellfun(@isempty,strfind({fileinfo.prefix}, solution_files_PREFIX)));
% base_files_IX = strfind([fileinfo.prefix], base_files_PREFIX);
% solution_files_IX = strfind([fileinfo.prefix], solution_files_PREFIX);

%%
n = 10;
durations = zeros(3,n);
for ii = 1:n
    tic
    sdf1 = importdata(filename, ' ');
    durations(1,ii) = toc;

    tic
    sdf2=readtable(filename);
    durations(2,ii) = toc;
    
    tic
    sdf2=textread(filename);
    durations(3,ii) = toc;
    
end
mean(durations,2)

%% extract data
clear data

% The first file contain the variable names
ii_solution_file = 1;
ii_file = solution_files_IX(ii_solution_file);
filename = fullfile(files(ii_file).folder, files(ii_file).name);
T = readtable(filename);
VariableNames = T.Properties.VariableNames;

% run over all files
for ii_solution_file = 1:50%length(solution_files_IX)
    %%
    ii_file = solution_files_IX(ii_solution_file);
    filename = fullfile(files(ii_file).folder, files(ii_file).name)
    
    
    %%
    continue
    switch ii_solution_file
        case 1
            % only the first file contain the variable names
            T = readtable(filename);
            VariableNames = T.Properties.VariableNames;
        otherwise
            if length(T.Properties.VariableNames) ~= length(VariableNames)
                continue;
            end
            T.Properties.VariableNames = VariableNames;
            data = [data; T];
    end
    
end

%% analyze/plot
clear ts;
for ii = 1:size(data,1)
    ts(ii) = datetime([data.x_{ii} ' ' data.GPST{ii}], 'InputFormat', 'yyyy/MM/dd HH:mm:ss.SSS');
end





%%













%%
