function PARAMS_SetParamset(paramset)

%% set the global variable
global paramset_global_var;
paramset_global_var = paramset;

%% create symbolic links for results folders
% 0 - shared field across paramsets
% 1 - specific field per paramset
field_options = {
    'details',          0;
    'spikes',           0;
    'RecStability',     0;
    'FE',               0;
    'FR_map',           0;
    'Ipos',             0;
    'fields',           1;
    'shuffle',          0;
    'signif',           1;
    'meanFR',           0;
    'time_AC',          0;
    'stats',            1;
    'inclusion',        0;
    'cluster_control',  0;
    'cluster_quality',  0;
    
    'figures',          1;
    'choose_examples',  1;
    };

% main folders
cells_dir = 'L:\Analysis\Results\cells';
cells_paramsets_dir = ['L:\Analysis\Results\cells_paramset_' num2str(paramset)];
cells_default_dir = ['L:\Analysis\Results\cells_paramset_' num2str(0)];

% run over results folders
for ii_fields = 1:length(field_options)
    
    % folders to link
    link_dir = fullfile(cells_dir, field_options{ii_fields,1} );
    if field_options{ii_fields,2}
        target_dir = fullfile(cells_paramsets_dir, field_options{ii_fields,1});
    else
        target_dir = fullfile(cells_default_dir, field_options{ii_fields,1});
    end
    
    % first, remove the link folder
    delete(link_dir)
    
    % now create the symbolic link!
    command = sprintf('mklink /D "%s" "%s"', link_dir , target_dir );
    [status,cmdout] = system(command,'-echo');
    if status
        fprintf('mklink error %d: %s', status,cmdout);
        error('Paramset mklink error!')
    end
    
end


end










%%
