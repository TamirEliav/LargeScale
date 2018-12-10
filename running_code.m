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

%% 22/11/2018 calc y data inposition projection (linearization)
IX = 10000:20000;
IX = 20000:40000;
plot(pos.raw.pos(IX,1),pos.raw.pos(IX,2),'.-')
hold on
plot(calib_tunnel.spline_fitresult)

%%
figure
hold on
IX = pos.proj.dir==1;
plot(pos.proj.pos(IX,1), pos.proj.pos(IX,2),'b.');
IX = pos.proj.dir==-1;
plot(pos.proj.pos(IX,1), pos.proj.pos(IX,2),'r.');
xlim([0 200])
ylim([-1 1])

%%
figure
hold on
plot(pos.raw.pos(:,1),pos.raw.pos(:,2),'.b')
plot(pos.calib_tunnel.spline_fitresult)

%% 26/11/2018 structure array stuff
clear
clc
subs = [1 nan 1 1 1 2 2 3  3  3]';
val  = [8 50  8 7 8 4 3 10 9 11]';
IX = ~isnan(subs);
A = accumarray(subs(IX),val(IX),[],@mean)
A = accumarray(subs(IX),val(IX),[],@median)
A = accumarray(subs(IX),val(IX),[],@range)

%% 27/11/2018 compare csaps using the indices vs. the time as X input (and how this relates to the smoothing factor)
% csaps_p1 = 1e-6;
csaps_p1 = 1e-8;
% csaps_p2 = 1-1e-1000;
interval = 1;
fs = pos.proc_1D.fs / interval;
t = pos.proc_1D.ts(1:interval:end);
x = pos.proc_1D.pos(1:interval:end);
pp1 = csaps(1:length(t), x, csaps_p1);
% pp2 = csaps(t, x, csaps_p2);
xpp1 = fnval(pp1, 1:length(t) );
% xpp2 = fnval(pp2, t );
vel = diff(x) .* fs;
vel_pp1 = fnval(fnder(pp1), 1:length(t)).* fs;
% vel_pp2 = fnval(fnder(pp2), t) .* 1e6;

figure
hold on
plot(t(2:end),vel,'.k')
plot(t,vel_pp1,'.b')
% plot(t,vel_pp2,'.r')
% legend({'raw';'pp1';'pp2'})
rescale_plot_data('x',[1e-6/60 t(1)])
xlim([42 45])
ylim([-15 15])
title(sprintf('fs=%d, p=%d',fs,csaps_p1))

%% 27/11/2018 compare usage of my new util function get_data_in_ti
clc
ti = 0:100:100000;
ti = [ti' ti'+30];
t = rand(1,10000);
t = normalize(t,'range', [min(ti(:)) max(ti(:))]);
tic
% [IX1, IX_per_ti] = get_data_in_ti(t,ti,1);
IX1 = get_data_in_ti(t,ti,1);
toc
tic
[IX2, IX_per_ti] = get_data_in_ti(t,ti,2);
toc
whos ti t IX1 IX2
length(IX1) / length(t)
length(IX2) / length(t)

sum(sort(IX1)~=sort(IX2))

IX_per_ti;
cellfun(@(x)(t(x)), IX_per_ti, 'UniformOutput',false);

%% 28/11/2018 test shir's fill holes function
clear
clc
% cell_ID = 'b0034_d180312_TT4_SS01';
% cell = cell_load_data(cell_ID);
% exp_ID = cell.details.exp_ID;
exp_ID = 'b0148_d170608';
% exp_create_position(exp_ID);
exp = exp_load_data(exp_ID);

figure
hold on
plot(exp.pos.proc_1D.ts,exp.pos.proc_1D.pos,'.k')
fs_old = 1e6/median(diff(exp.pos.raw.ts_nlg_usec));
fs_new = exp.pos.proc_1D.fs;
ib = find_nearest_point(exp.pos.proc_1D.ts, exp.pos.raw.ts_nlg_usec);
dist_from_orig_ts = abs(exp.pos.proc_1D.ts - exp.pos.raw.ts_nlg_usec(ib)');
IX = dist_from_orig_ts >(1e6/fs_old);
plot(exp.pos.proc_1D.ts(IX),exp.pos.proc_1D.pos(IX),'.r')

%%
sdf = repelem(struct,10);
[sdf(:).a] = disperse([1:10:100;(1:10:100)+3])

%%
clc
sdf = {[1 2 3], [11 12]};
sdf2 = cellfun(@(x)(ones(size(x))), sdf, 'UniformOutput', false)
% [sdf{:}]

%% 04/12/2018 - compare pos vs. csaps
exp_ID = 'b9861_d180711'
exp=exp_load_data(exp_ID);
err = exp.pos.proc_1D.pos - exp.pos.proc_1D.pos_csaps;
plot(exp.pos.proc_1D.vel_csaps,err,'.')

%%
cells_list = {'b0034_d180312_TT4_SS01';'b0034_d180312_TT4_SS03'};
cells = cellfun(@(x)(cell_load_data(x,'details','FE','fields')), cells_list, 'UniformOutput', false);
cells = [cells{:}];

%%
% edges = [1 2; 3 4; 5 6]'
edges = {[1 2],[3 4],[5 6]}
sdf = repelem(struct(),3);
[sdf(:).width] = disperse([1 2 3]);
[sdf(:).edges] = disperse(edges);
% sdf.width
sdf.edges

%%
edges = [1 2; 3 4; 5 6]
edges = num2cell(edges,2)

%%
cells=cellfun(@(x)(cell_load_data(x,'details','fields','FR_map','Ipos','RecStability','signif')), cells_t.cell_ID,'UniformOutput', false);
cells = [cells{:}];
num_fields = cellfun(@(x)([length(x{1}) length(x{2})]),{cells.fields}, 'UniformOutput', false);
num_fields = cat(1,num_fields{:});
figure
plot(num_fields(:,1), num_fields(:,2), '.')

%%
signif = cellfun(@(x)([x(:).TF]),{cells.signif}, 'UniformOutput', false);
signif = cat(1,signif{:});

%%
FR_maps = cellfun(@(x)([x(:).all]),{cells.FR_map}, 'UniformOutput', false);
FR_maps = cat(1,FR_maps{:});
SI_bits_spike  = arrayfun(@(x)(x.SI_bits_spike), FR_maps)
plot(SI_bits_spike(:,1),SI_bits_spike(:,2),'.')
axis equal
refline(1,0)
xlim([0 8])
ylim([0 8])
xlabel('SI direction 1 (bits/spike)')
ylabel('SI direction 2 (bits/spike)')
title('comparing spatial info between flight directions')
saveas(gcf,'L:/Analysis/Results/comparing spatial info between flight directions','tif');

%%
fields1 = cellfun(@(x)(x{1}),{cells.fields}, 'UniformOutput', false)
cellfun(@(x)( sum([x.width_prc]) ), cells(:).fields{1}, 'UniformOutput', false);
% sum([sdf.width_prc])
% fields_total_area = sum(arrayfun(@(x)(x.width_prc), sdf))

%% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % Internal functions section  % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
















%%





