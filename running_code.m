%%
clear
clc

%% git commmands
% git status
% git add -A
% git commit -a -m "comment"
% git push
% git pull
% git fetch -v --dry-run
% git remote update

%% merge session ts in recording_summary to one single column
clear
clc
exp_t = DS_get_exp_summary();
sessions_names = {'Sleep1';'Behave';'Sleep2'};
% exp_t.sessions_names = {};
T = table();
for ii_exp = 1:height(exp_t)
    exp = exp_t(ii_exp,:);
    
    sessions_ts = [
    exp.session_sleep1_start_ts exp.session_sleep1_end_ts;    
    exp.session_behave_start_ts exp.session_behave_end_ts;
    exp.session_sleep2_start_ts exp.session_sleep2_end_ts;
    ];
    T.sessions_ts{ii_exp} = mat2str(sessions_ts);
    T.sessions_names{ii_exp} = sessions_names;
end
writetable(T,'testTable','FileType','spreadsheet')

%%
clear
clc
exp_t = DS_get_exp_summary();

%% 19/11/2108 disp events (to fill the session ts in exp_summary)
clc
exp_ID = 'b0079_d160902';

[exp_path exp_info] = DS_get_path(exp_ID);

if ~exist(exp_path.nlx,'dir')
    disp('create EVENT files')
    Nlg2Nlx(exp_path.raw);
end
PRE_event_disp(exp_ID);

filename_template = sprintf('L:\\DATA\\%04d_%s\\%s\\nlg\\LOG*.txt',...
    exp_t.batNum(exp_ID),...
    exp_t.batName{exp_ID},...
    datestr(exp_t.date(exp_ID),'yyyymmdd'));
file = dir(filename_template);
fid = fopen(fullfile(file.folder,file.name),'rt');
while ~feof(fid)
    tline = fgetl(fid);
    if contains(tline, {'***','Sent string for remote log'})
        disp(tline)
    end
end


%% 20/11/2018 rename sorted NTT files from Shir
clear
clc
bat = '9861';
dir_main = ['L:\Analysis\pre_proc\SpikeSorting\bat' bat '\2*'];
folders = dir(dir_main);
for ii_folder = 1:length(folders)
    dir_files = fullfile(folders(ii_folder).folder, folders(ii_folder).name,'spikes_NTT');
    files = dir(dir_files);
    files = files(~[files.isdir]);
    file_tmpl = ['spikes_b' bat '_d' folders(ii_folder).name(3:end) '_'];
    for ii_file = 1:length(files)
        filename_orig = fullfile(files(ii_file).folder, files(ii_file).name);
        filename_new = strrep(filename_orig, 'spikes__', file_tmpl);
        [status,msg,msgID] = movefile(filename_orig,filename_new)
    end
end


%% 20/11/2018 plot all the cells that passed all criteria EXCEPT for IsoDist>10
clear
clc
load('L:\Analysis\Results_OLD_20181118\posters\SfN_2018\population\inclusion_list_IsoDist_10.mat')
load('L:\Analysis\Results_OLD_20181118\posters\SfN_2018\population\pop_data.mat')
T=pop_data.T;

SI_thr = 0.5;
IsoDist_thr = 10;
mean_FR_thr = 5;
corr_even_odd_thr = 0.5;
non_repeating_cells_IX = repmat(cellfun(@isempty, T.same_cell),1,2);
dCA1_cells_IX = repmat(strcmp(T.brain_area,'dCA1'),1,2);
good_isolated_cells = repmat(pop_data.IsoDist' > IsoDist_thr,1,2);
pyramidal_cells = repmat(pop_data.meanFR < mean_FR_thr,2,1)';
% pyramidal_cells = pop_data.meanPSTH < mean_FR_thr;
SI_thr_cross_cells = pop_data.SI_bit_spikes > SI_thr;
SI_signif_cell = boolean(pop_data.signif);
stable_cells = pop_data.corr_even_odd > corr_even_odd_thr;

inc_list_no_IsoDist = non_repeating_cells_IX & ...
        dCA1_cells_IX & ...
        pyramidal_cells & ...
        stable_cells & ...
        SI_thr_cross_cells & ...
        SI_signif_cell;
IX = find(any(inc_list_no_IsoDist') & ~any(incl_list_cell_dir'))

dir_IN = 'L:\Analysis\from_Shir_20181120\cells_0148';
dir_OUT = 'L:\Analysis\Results\21081120_cells_without_IsoDist_thr';
% plot cell data for all of those cells
for ii_cell = 1:length(IX)
    %%
    ii_cell 
    cell_ID = pop_data.cell_ID{IX(ii_cell)}
    
    %% load cell data 
%     cell_load_data(cell_ID, 'details');
    filename = fullfile('L:\Analysis\from_Shir_20181011\Wild_cells_rearranged',['cell_data_' cell_ID]);
    cell_data = load(filename,'details','data','session','timeStability','clusterQuality');
    cell_ID_shir = [cell_data.details.day '_bat' cell_data.details.bat '_TT' cell_data.details.TT '_SS_' cell_data.details.SS];
    
    %% look for already existing plot of this cell (from shir)
    file_IN = fullfile(dir_IN,['spikes_' cell_ID_shir '.tif']);
    file_OUT = fullfile(dir_OUT, ['cell_' cell_ID '.tif']);
    if ~isempty(dir(file_IN))
        % copy the figure
        [status,msg,msgID] = copyfile(file_IN, file_OUT);
    else 
        %% create the figure and save it
        fig_size = [20 20];
        figure('Units','centimeters', 'PaperPosition',[0 0 fig_size], 'Position',[2 2 fig_size]);
        pnl=panel();
        pnl.pack('h',[80 20])
        pnl(1).pack('v',[20 30 30 20])
%         pnl(2).pack('v',[30 70])
        pnl.margin = 20;
        pnl.de.margin = 15;
        h=pnl.title(cell_ID); set(h,'interpreter','none','fontsize',14,'position',[0.5 1.05]);
        prm = PARAMS_GetAll();
        t0 = cell_data.session.behavioral.start;
        for ii_dir = 1:2
            dir_color = prm.graphics.colors.flight_directions{ii_dir};
            pnl(1,1).select(); hold on
            plot(cell_data.data(ii_dir).bin_centers, cell_data.data(ii_dir).PSTH,...
                'color',dir_color, 'LineWidth', 1.5);
            xlabel('Position (m)')
            ylabel('Firing rate (Hz)')
            
            pnl(1,ii_dir+1).select(); hold on
            pos = [cell_data.data(ii_dir).flights.pos];
            pos_ts = [cell_data.data(ii_dir).flights.ts_nlg_usec];
            spikes_pos = [cell_data.data(ii_dir).flights.spike_pos];
            spikes_ts = [cell_data.data(ii_dir).flights.spike_ts];
            plot(pos, pos_ts, '.', 'color', 0.9.*[1 1 1], 'MarkerSize',0.1);
            plot(spikes_pos, spikes_ts, '.', 'color', dir_color);
            rescale_plot_data('y',[1e-6/60,t0]);
            xlabel('Position (m)')
            ylabel('Time (min)')
        end
        pnl(1,4).select(); hold on
        binned_FR = [cell_data.timeStability.pre_sleep_fr cell_data.timeStability.behave_fr' cell_data.timeStability.post_sleep_fr];
        t = cell_data.timeStability.time_bins_behave_min';
        t = [t(1)-20 t t(end)+20];
        bar(t,binned_FR);
        xlabel('Time (min)')
        ylabel('Firing rate (Hz)')
        
        pnl(2).select(); hold on
        axis off
        text(0,1,{...
            sprintf('Isolation Distance: %g', cell_data.clusterQuality.Isolation_dis),...
            sprintf('L-ratio: %g', cell_data.clusterQuality.L_Ratio),...
            },'Units','normalized')
        
        % save figure
        file_OUT = fullfile(dir_OUT, ['cell_' cell_ID]);
        saveas(gcf,file_OUT,'tif');
        close(gcf)
    end
    

end

%% 21/11/2018 rename sorted NTT files from Shir (remove '_')
clear
clc
bat = '0148';
dir_main = ['L:\Analysis\pre_proc\SpikeSorting\' bat '\2*'];
folders = dir(dir_main);
for ii_folder = 1:length(folders)
    dir_files = fullfile(folders(ii_folder).folder, folders(ii_folder).name,'spikes_NTT');
    files = dir(dir_files);
    files = files(~[files.isdir]);
    file_tmpl = ['spikes_b' bat '_d' folders(ii_folder).name(3:end) '_'];
    for ii_file = 1:length(files)
        filename_orig = fullfile(files(ii_file).folder, files(ii_file).name);
%         filename_new = strrep(filename_orig, '_.', '.');
        filename_new = regexprep(filename_orig, ['spikes_\d+_bat' bat '_'], file_tmpl);
        [status,msg,msgID] = movefile(filename_orig,filename_new)
    end
end

%%




%% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % Internal functions section  % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
















%%





