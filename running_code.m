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
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%
cells=cellfun(@(x)(cell_load_data(x,'details','fields','FR_map','Ipos','RecStability','signif','stats')), cells_t.cell_ID,'UniformOutput', false);
cells = [cells{:}];

%% remove IN
mean_FR_flight = arrayfun(@(x)(x.all.meanFR_flight),[cells.stats]);
mean_FR_all = arrayfun(@(x)(x.all.meanFR_all),[cells.stats]);

figure
subplot(1,2,1)
hold on
bins = [0:0.2:10];
histogram(mean_FR_flight,bins)
histogram(mean_FR_all,bins)
legend({'flight';'all'})
subplot(1,2,2)
plot(repelem([1;2],1,length(mean_FR_flight)), [mean_FR_flight;mean_FR_all],'.-')
h=gca;
h.TickDir='both';
h.XTick=[1 2];
h.XTickLabel={'flight';'all'};

%% take only signif cells (in both directions)
signif = cellfun(@(x)([x.TF]),{cells.signif},'UniformOutput',0);
signif = cat(1,signif{:});
signif = all(signif,2);
cells = cells(signif);

%% get num fields
num_fields = cellfun(@(x)([length(x{1}) length(x{2})]),{cells.fields}, 'UniformOutput', false);
num_fields = cat(1,num_fields{:});

%% plot num fields DIR1 vs. DIR2
figure
plot_std = 0.1;
scatter(num_fields(:,1)+plot_std.*randn(size(num_fields,1),1),...
        num_fields(:,2)+plot_std.*randn(size(num_fields,1),1),...
        10, arrayfun(@(x)(find(unique(bats) == x)), bats),'filled')
axis equal
refline(1,0)
xlim([0 20])
ylim([0 20])
xlabel('#fields direction 1')
ylabel('#fields direction 2')
title('comparing number of fields between flight directions')
saveas(gcf,'L:/Analysis/Results/comparing number of fields between flight directions','tif');

%% plot SI DIR1 vs. DIR2
FR_maps = cellfun(@(x)([x(:).all]),{cells.FR_map}, 'UniformOutput', false);
FR_maps = cat(1,FR_maps{:});
SI_bits_spike  = arrayfun(@(x)(x.SI_bits_spike), FR_maps);
scatter(SI_bits_spike(:,1),SI_bits_spike(:,2),...
        10, arrayfun(@(x)(find(unique(bats) == x)), bats),'filled')
axis equal
xlim([0 5])
ylim([0 5])
refline(1,0)
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
clc
pnl = panel();
pnl.pack(2,3)
sdf1=pnl.ch
sdf2=pnl.de
pnl.select('all')
pnl.identify()
cellfun(@(x)(x.axis),sdf2,'UniformOutput',false)
sdf2(cellfun(@(x)(isempty(x.axis)),sdf2)) = []
% arrayfun(@(ax)({ax.Position(3:4) = [0.2 0.2]}), pnl.de.axis)
% for ii_ax = 1:length(pnl.de.axis)
%     ax = pnl.de.axis(ii_ax)
%     ax.Position(3:4) = [0.2 0.2]
% end

%%
exp_ID ='b0034_d180312';
% exp=exp_load_data(exp_ID);
prm = PARAMS_GetAll();
pos = exp.pos.proj.pos;
ts = exp.pos.proj.ts;
pos(exp.pos.proj.outliers_IX,:)=[];
ts(exp.pos.proj.outliers_IX)=[];
ti = [[exp.flight.FE.start_ts];[exp.flight.FE.end_ts]]';
IX = get_data_in_ti(ts',ti);
ts = ts(IX);
pos = pos(IX,:);
directions = interp1(exp.pos.proc_1D.ts, sign(exp.pos.proc_1D.vel_csaps), ts);
dir_colors = prm.graphics.colors.flight_directions;
dir_sign = [1 -1];
figure
hold on
bin_edges = 0:0.5:200;
bin_centers = edges2centers(bin_edges);
for ii_dir = 1:2
    IX = directions == dir_sign(ii_dir);
    plot(pos(IX,1),pos(IX,2),'.', 'color', dir_colors{ii_dir},'MarkerSize',1);
    [~,~,BIN] = histcounts(pos(IX,1),bin_edges);
    ymean = accumarray(BIN,pos(IX,2),[],@mean);
    ystd = accumarray(BIN,pos(IX,2),[],@std);
    
    x = bin_centers;
    y = nan(size(bin_centers));
    err = nan(size(bin_centers));
    y(unique(BIN)) = ymean(unique(BIN));
    err(unique(BIN)) = ystd(unique(BIN));
    
    yyaxis left
    hold on
    h=shadedErrorBar(x, y, err,'lineprops',{'color', dir_colors{ii_dir}});
    
    yyaxis right
    hold on
    plot(x, err, '.-', 'color', dir_colors{ii_dir});
end

%%
figure
t = 1:100;
yyaxis left
plot(t, sin(t))
rescale_plot_data('x',[1e-1 50])

yyaxis right
plot(t, cos(t))

rescale_plot_data('x',[1e-1 50])


%%
cells = cellfun(@(x)(cell_load_data(x,'details','spikes')), cells_t.cell_ID, 'UniformOutput', false);
cells = [cells{:}];
%%
CQ = arrayfun(@(x)(x.ClusterQuality),[cells.details]);
CQ = CQ==1;
IsoDist = arrayfun(@(x)(x.Isolation_dis),[cells.spikes]);
figure
hold on
data = {IsoDist(~CQ),IsoDist(CQ)};
plotSpread(data)
violin(data)
h=gca;
h.XTick = [1 2];
h.XTickLabel = {'bad','good'};
title('bat 34 IsoDist comparing bad/good manualy quality')
ylabel('IsoDist')


%% 24/12/2018 - debug sync problems (bsp<->nlg)
clear 
clc
exp_ID = 'b0034_d180313';
exp=exp_load_data(exp_ID);

% BSP TTL
load( fullfile(exp.path.bsp, 'bsp_TTL.mat') );
bsp_TTL_ts_msec = round(1e-6.*bsp_TTL_ts_ns');
bsp_TTL_intervals = diff(bsp_TTL_ts_msec);
bsp_TTL_intervals_inc = diff(bsp_TTL_intervals);
% NLG TTL
nlg_TTL_file_name = fullfile(exp.path.nlx, 'EVENTS__Digital in.nev');
FieldSelection = [1 0 0 0 0];
ExtractHeader = 0;
ExtractMode = 1;
ModeArray = [];
nlg_TTL_ts_usec = Nlx2MatEV( nlg_TTL_file_name ,FieldSelection,ExtractHeader,ExtractMode,ModeArray);
nlg_TTL_ts_msec = nlg_TTL_ts_usec*1e-3;
nlg_TTL_intervals = diff(nlg_TTL_ts_msec);
nlg_TTL_intervals_inc = diff(nlg_TTL_intervals);

whos bsp_TTL_ts_msec nlg_TTL_ts_msec

time_conv_p_msec = sync_TTL_polyfit(round(bsp_TTL_ts_msec), round(nlg_TTL_ts_msec), 2, 100);


%%
thr = 5;
x = diff(bsp_TTL_ts_msec);
y = diff(nlg_TTL_ts_msec);
[dist,ix,iy] = dtw(x,y);
pairs = [x(ix);y(iy)];
rsdl = diff(pairs);
IX = find( abs(rsdl) < thr );
pairs_good = pairs(:,IX);

if length(unique(pairs(1,IX))) ~= length(IX) || ...
   length(unique(pairs(2,IX))) ~= length(IX)
   error()
end

mathing_TTL_bsp = bsp_TTL_ts_msec(union(ix(IX), ix(IX)+1));
mathing_TTL_nlg = nlg_TTL_ts_msec(union(iy(IX), iy(IX)+1));

%%
a=unique(ix(IX))
b=unique(iy(IX))
if length(a)~=length(b)
    error
%     TODO: choose the closest ones to the bsp!! (because bsp is the
%     master)
end
pairs_x_IX = union(a,a+1);
pairs_y_IX = union(b,b+1);
whos x y ix iy a b pairs_x_IX pairs_y_IX IX

%%
figure
subplot(2,1,1); hold on
plot(rsdl,   '.-b')
plot(IX,rsdl(IX),'or')
plot(get(gca,'xlim'),  [thr thr],'m--')
plot(get(gca,'xlim'), -[thr thr],'m--')
xlabel('interval pair index')
ylabel('interval pair diff (msec)')
% zoom-in 
subplot(2,1,2); hold on
plot(rsdl,   '.-b')
plot(IX,rsdl(IX),'or')
plot(get(gca,'xlim'),  [thr thr],'m--')
plot(get(gca,'xlim'), -[thr thr],'m--')
ylim([-100 100])
xlabel('interval pair index')
ylabel('interval pair diff (msec)')
title('zoom-in')
suptitle('finding matching pairs of intervals');


%% 27/12/2018 plexon header stuff...
clear
% clc
% file_IN = 'L:\Analysis\pre_proc\SpikeSorting\0034\20180310\spikes_NTT\spikes_b0034_d180310_TT3.NTT';
% file_OUT = 'L:\Analysis\pre_proc\SpikeSorting\0034\20180310\spikes_NTT\spikes_b0034_d180310_TT3_plexon.NTT';
% file_IN = 'D:\Tamir\PROJECTS\Plexon\Data\spikes_b0148_d170608_TT4_.NTT';
% file_OUT = 'D:\Tamir\PROJECTS\Plexon\Data\spikes_b0148_d170608_TT4__new_header_++.NTT';
% file_IN = 'D:\Tamir\PROJECTS\Plexon\Data\spikes_b0034_d180310_TT1_shir.NTT';
% file_OUT = 'D:\Tamir\PROJECTS\Plexon\Data\spikes_b0034_d180310_TT1__plexon.NTT';

header_file = 'plexon_header_NTT.txt';
header_new = textread(header_file, '%s', 'delimiter', '\n', 'whitespace', '');

% [Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] =...
%     Nlx2MatSpike(file_IN , [1 1 1 1 1], 1, 1, [] );
% whos

% ScNumbers(:) = 1;
% CellNumbers(:) = 0;
% Mat2NlxSpike(file_OUT, 0, 1, [], [1 1 1 1 1 1], ...
%     Timestamps, ScNumbers, CellNumbers, Features, Samples, header_new);
% Mat2NlxSpike(file_OUT, 0, 1, [], [1 0 0 1 1 1], ...
%     Timestamps, Features, Samples, header_new);

%% cont - create demo data
clear
clc
file_OUT = 'D:\Tamir\PROJECTS\Plexon\Data\test_20181227.NTT';
header_file = 'plexon_header_NTT.txt';
Header = textread(header_file, '%s', 'delimiter', '\n', 'whitespace', '');

n = 1e4;
Timestamps = linspace(0,1e6*60*60,n);
ScNumbers = zeros(1,n);
CellNumbers = zeros(1,n);
% Features = 
spike_shape = normpdf(1:32,8,2)';
spike_shape = spike_shape ./ max(spike_shape);
spike_shape = spike_shape .* 250;
% spike_shape = spike_shape + randn(size(spike_shape));
Samples = repmat(spike_shape,1,4,n);
Samples = Samples + randn(size(Samples)).*10;

AppendToFileFlag = 0;
ExportMode = 1;
ExportModeVector = [];
FieldSelectionFlags = [1 1 1 0 1 1];
Mat2NlxSpike( file_OUT, AppendToFileFlag, ExportMode, ExportModeVector,...
              FieldSelectionFlags, Timestamps, ScNumbers, CellNumbers, Samples, Header);

%% 30/12/2018 - nlx header stuff...
clear
clc
file_IN = 'D:\Tamir\PROJECTS\Plexon\Data\bat6255_Day120326_1_PreSleep_TT1.Ntt';
file_OUT = 'D:\Tamir\PROJECTS\Plexon\Data\bat6255_Day120326_1_PreSleep_TT1_copy.Ntt';
[Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] =...
    Nlx2MatSpike(file_IN, [1 1 1 1 1], 1, 1, [] );

ADMaxValue = 32767;
ADC=0.000000009155552760375940;
figure
plot(mean(Samples,[2 3]).*ADC*1e6)

% ADC=0.000000009155552760375940;
% ADC=0.000000009160000000000000;

% ADC_new = ADC;
InputRange_new = max(Samples(:)) *ADC*1e6;
% ADC_new = 1e-6;
ADC_new = InputRange_new / ADMaxValue / 1e6;

Samples = Samples .* ADC ./ ADC_new;
Samples = round(Samples);

% InputRange = ADMaxValue * ADC_new * 1e6;
% InputRange = round(InputRange);

ADC_str = sprintf('%.24f',ADC_new);
InputRange_str = sprintf('%g',InputRange_new);
Header{16} = sprintf('-ADBitVolts %s %s %s %s', ADC_str, ADC_str, ADC_str, ADC_str);
Header{21} = sprintf('-InputRange %s %s %s %s', InputRange_str, InputRange_str, InputRange_str, InputRange_str);

figure
plot(mean(Samples,[2 3]).*ADC_new*1e6)

Mat2NlxSpike(file_OUT, 0, 1, [], [1 1 1 1 1 1], ...
    Timestamps, ScNumbers, CellNumbers, Features, Samples, Header);
whos

% close all

%%
file_IN = 'L:\Analysis\pre_proc\SpikeSorting\0034\20180312\spikes_NTT\spikes_b0034_d180312_TT4.NTT';
nlx_change_header(file_IN, []);

%% 01/01/2019 - change header for NTTs I got from Shir
clear 
clc
dir_IN = 'L:\Analysis\pre_proc\SpikeSorting_Shir';
dir_OUT = 'L:\Analysis\pre_proc\SpikeSorting';
files = dir( fullfile(dir_IN,'**','*.NTT') )
for ii_file = 1:length(files)
    file = files(ii_file);
    file_IN = fullfile(file.folder,file.name)
    file_OUT = strrep(file_IN, dir_IN, dir_OUT);
    nlx_change_header(file_IN, file_OUT);
end


%% 03/01/2019 - plexon CSC spikes detection test
% I did the following steps:
% 1. loaded raw ncs data to plexon
% 2. saved it as PLX file
% 3. convert 4 channels to 1 tetrode 
% 4. filter
% 5. thr detection (only positive)
% 6. saved only the spikes plx channel
% 7. exported the waveforms to matlab (only sorted waveforms), with the
% following data columns:
%      col. 1: unit
%      col. 2: ts
%      col. 3-end: wvfrm
% now I want to load it and save it as ntt file
clear 
clc
% load data
file_IN = 'D:\__TEMP\plexon_CSC_test\TT4.mat';
file_OUT = strrep(file_IN,'.mat','.NTT');
load(file_IN);
TT4 = CSC13;
% arrange data
% remove invalidated waveforms
TT4( TT4(:,1)==-1 , :) = [];
CellNumbers = TT4(:,1)';
Timestamps = TT4(:,2)' .* 1e6;
Samples = TT4(:,4:end)';
Samples = reshape(Samples,32,4,length(CellNumbers));
% plot( mean(Samples,[2 3]) )
% histogram(categorical(CellNumbers))
% convert to ntt
header_file = 'Nlx_header_NTT.txt';
header_new = textread(header_file, '%s', 'delimiter', '\n', 'whitespace', '');
ADMaxValue = 32767;
InputRange = max(Samples(:));
ADC = InputRange / ADMaxValue / 1e6;
Samples = Samples ./ ADC ./ 1e6;
ADC_str = sprintf('%.24f',ADC);
InputRange_str = sprintf('%g',InputRange);
ADC_str_IX = contains(header_new, 'ADBitVolts');
InputRange_str_IX = contains(header_new, 'InputRange');
header_new{ADC_str_IX} = sprintf('-ADBitVolts %s %s %s %s', ADC_str, ADC_str, ADC_str, ADC_str);
header_new{InputRange_str_IX} = sprintf('-InputRange %s %s %s %s', InputRange_str, InputRange_str, InputRange_str, InputRange_str);
% save ntt file
Mat2NlxSpike(file_OUT, 0, 1, [], [1 0 1 0 1 1], ...
    Timestamps, CellNumbers, Samples, header_new);






%% try to load spikes from PLX file!
clear 
clc
filename = 'D:\__TEMP\plexon_CSC_test\TT4-01.plx';
channel = 'CSC13';
unit = 2;
% [n, npw, ts, wave] = plx_waves_v(filename, channel, unit);
% wave = wave.*1e3;
% plot(mean(wave))

% [OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreTresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information(filename)

% [adfreq, n, ts, fn, ad] = plx_ad_v(filename, channel);
[n,names] = plx_chan_names(filename)
% [tscounts, wfcounts, evcounts, slowcounts] = plx_info(filename,1)

whos


%% develop fast coincident detector (vectorized)
clear 
clc
x = zeros(4,25);
x(:,5)=1;
x([1 2 3],10)=1;
x([4],11)=1;
x(1,20)=1;
x(2,21)=1;
x(3,22)=1;
x(4,23)=1;
x

n = 3;
win = rectwin(n);
CD_filtfilt = filtfilt(win, 1, sum(x)) ./ size(x,1) ./ n
% CD_conv1 = conv(sum(x),win,'full')
CD_conv2 = conv(sum(x),win,'same')
% CD_conv3 = conv(sum(x),win,'valid')
% T = table(x',CD_filtfilt',CD_conv1',CD_conv2',CD_conv3')
T = table(x',CD_filtfilt',CD_conv2')

%% test new spike detection
clear 
clc
exp_ID = 'b0034_d180306';
exp = exp_load_data(exp_ID,'details','path');
dir_IN = exp.path.spikes_raw;
dir_OUT = exp.path.spikes_detection;
dir_OUT = [dir_OUT '_old_code'];
% dir_IN = 'L:\Analysis\pre_proc\0034\20180306\spikes_raw';
% dir_OUT = 'L:\Analysis\pre_proc\0034\20180306\spikes_detection_test2';
% dir_OUT = 'L:\Analysis\pre_proc\0034\20180315\spikes_detection_test_runtime4_sleep_original';
% dir_OUT = 'L:\Analysis\pre_proc\0034\20180315\spikes_detection_test_runtime5_sleep_new_bug';
% dir_OUT = 'L:\Analysis\pre_proc\0034\20180315\spikes_detection_test_runtime6_sleep_new_with_abs_in_lib';
% dir_OUT = 'L:\Analysis\pre_proc\0034\20180315\spikes_detection_test_runtime7_all_new_no_abs_lib_yes_abs_align';
% dir_OUT = 'L:\Analysis\pre_proc\0034\20180315\spikes_detection_test_runtime8_all_new_yes_abs_lib_yes_abs_align';
% dir_OUT = 'L:\Analysis\pre_proc\0034\20180315\spikes_detection_test_runtime9_updated_code_sleep_only';
% dir_OUT = 'L:\Analysis\pre_proc\0034\20180315\spikes_detection_test_runtime10_Sleep1_report_thr_std_median';

params.thr_uV = [repelem(30,4);...
                 repelem(25,4);...
                 repelem(30,4);...
                 repelem(35,4)];
% params.thr_uV = 50;
params.ref_ch = [2 3];
% params.t_start_end = [61004856668 61526309828];
params.t_start_end = [61639062973 62365439740];
params.t_start_end = [];
% params.t_start_end = exp_get_sessions_ti(exp_ID,'Sleep1');
% params.t_start_end = mean(params.t_start_end) + 1e6*60*2.*[-1 1];
params.use_neg_thr = 0;
params.TT_to_use = [1 2 3 4];
params.merge_thr_crs_width = 4;
params.lib_spike_shapes = 'library_of_acceptable_spike_shapes_new.mat';
params.lib_corr_thr = 0.9;
params.active_TT_channels = ones(4,4);
params.CD_detect_win_len = 32;
params.CD_invalid_win_len = 32*2;
params.CD_n_TT_thr = 3;
params.CD_n_ch_thr = 9;
% params.CD_n_TT_thr  = length(params.TT_to_use);
% params.CD_n_ch_thr = 0.5 * sum(params.active_TT_channels(:)); % at least on half of the channels
params.is_save_artifacts = 1;
% run!
Nlx_detect_spikes_CSC(dir_IN,dir_OUT,params)
% Nlx_detect_spikes_CSC3(dir_IN,dir_OUT,params)


%% convert plx to ntt
clear
clc
% file_IN = 'L:\Analysis\pre_proc\0034\20180315\spikes_detection\spikes_b0034_d180315_TT4-02.plx';
file_IN = 'L:\Analysis\pre_proc\0034\20180306\spikes_detection_new_code\spikes_b0034_d180306_TT4_.plx';

% file_OUT = 'D:\__TEMP\spike_detection\spikes_detection\spikes_b0034_d180314_TT1-11.NTT';
plexon_plx2ntt(file_IN, [])
% plexon_plx2ntt(file_IN, file_OUT)

%% test convolution idea to make coincidence detection faster
sdf = zeros(length(Timestamps_accepted_spikes_TT), n_samples);
for ii_TT = 1:length(Timestamps_accepted_spikes_TT)
    spikes_ts_IX = ismember(timestamps,Timestamps_accepted_spikes_TT{ii_TT});
	sdf(ii_TT,spikes_ts_IX) = 1;
end
CD_win_len = 32;
win = rectwin(CD_win_len);
CD_conv2 = conv(sum(sdf),win,'same');
plot(CD_conv2)

%% play with spaese matrix
clear
clc
nSamples = 32*1e3*60*60*2.5;
nSpikes = 1e5;
nTT = 16;
for ii = 1:nTT
    spikes_IX(ii,:) = randsample(nSamples,nSpikes);
end

%%
tic
for ii = 1:nTT
    A = sparse(ii, spikes_IX(ii,:), 1,nTT,nSamples);
end
toc

tic
B = zeros(nTT,nSamples);
for ii = 1:nTT
    B(ii, spikes_IX(ii,:)) = 1;
end
toc

%% play with findpeaks
clear 
clc
nSamples = 32*1e3*60*60*2.5;
nSpikes = 1e5;
% nSamples = 100;
% nSpikes = 20;
x = zeros(nSamples,1);
x( randsample(nSamples,nSpikes) ) = rand(nSpikes,1)+1;

%%
tic
sdf=findpeaks(x,'MinPeakDistance',5);
toc
tic
sdf=findpeaks(x,'MinPeakDistance',5, 'MinPeakHeight', 1);
toc

%% 22/01/2019 - test if accum array is faster than for loop
clear
clc
n = 1e5;
k = 32;
x = zeros(1,k*n);
x(5:k:length(x)) = 150;
subs = ceil((1:length(x))./k);
tic
accumarray(subs',x',[],@max);
toc
tic
accumarray(subs',x',[],@my_func);
toc
plot(subs)

%%
clc
iter = 10000;
tic
for ii=1:iter
    a = max(x);
end
toc

tic
for ii=1:iter
    [a,b] = max(x);
end
toc

%% test std and median runtime for large arrays
fs = 32000;
x = randn(1,fs*60*60*2);
tic; std(x); toc
tic; median(abs(x)); toc
tic; mean(abs(x)); toc

%% plot median/std of neural raw channels for bats across days
clear;clc
% load exp summary and choose exps
exp_t = DS_get_exp_summary();
bat_num = 34;
exp_t(~contains(exp_t.recordingArena, '200m'),:) = [];
exp_t(exp_t.position_data_exist==0,:) = [];
exp_t(exp_t.neural_data_exist==0,:) = [];
exp_t(~ismember(exp_t.batNum, [bat_num] ),:) = [];
exp_t

exps = cellfun(@(x)(exp_load_data(x,'details','csc_raw_stats')), exp_t.exp_ID, 'UniformOutput', false);
exps = exps(cellfun(@(x)(isfield(x,'csc_raw_stats')),exps));
exps = [exps{:}];
dates = cellfun(@(x)(x.date), {exps.details});
stats = [exps.csc_raw_stats];
abs_median = cat(3,stats.csc_abs_median);
reref_abs_median = cat(3,stats.csc_reref_abs_median);

% plot to figure
figure('Units','normalized','Position',[0 0 1 1])
pnl = panel();
pnl.pack('v',2,'h',4);
pnl.margin=25;
h=pnl.title(sprintf('Raw neural data noise levels, bat %04d',bat_num));
h.FontSize = 20;
h.Position = [0.5 1.05];
h=pnl(1).title('No re-ref');
h.FontSize = 16;
for TT=1:4
    pnl(1,TT).select();
    title(sprintf('TT %d',TT))
    plot(dates, squeeze(abs_median(TT,:,:)),'.-', 'LineWidth',2);
    xlabel('Date')
    ylabel('median(abs(raw data)) [uVolt]')
    ylim([0 10])
    legend({'ch1','ch2','ch3','ch4'})
end
h=pnl(2).title('with re-ref');
h.FontSize = 16;
for TT=1:4
    pnl(2,TT).select();
    title(sprintf('TT %d',TT))
    plot(dates, squeeze(reref_abs_median (TT,:,:)),'.-', 'LineWidth',2);
    xlabel('Date')
    ylabel('median(abs(raw data)) [uVolt]')
    ylim([0 10])
    legend({'ch1','ch2','ch3','ch4'})
end
% save figure
fig_filename = fullfile('L:\Analysis\Results', sprintf('csc_raw_stats_bat_%04d',bat_num));
saveas(gcf,fig_filename, 'tif')

%%
clear;clc
limits_ts = [60389322391 60496267561];
file_in = 'L:\DATA\0034_Ace\20180228\nlx\CSC0.ncs';
file_out = 'L:\DATA\0034_Ace\20180228\nlx\CSC0__.ncs';
[signal, ts, fs, params] = Nlx_csc_read(file_in, limits_ts);
signal = signal(1:2:end);
ts = ts(1:2:end);
fs = fs/2;
plot(ts,signal)
nlx_csc_write(file_out , signal, ts, fs)

%% nested diary
diary off
diary('diary1'); diary on
disp('external')
diary('diary2'); diary on
disp('internal')
diary off
diary off
% can't use nested diary (over-write!)

%% plx2ntt 
PLX_filename = 'D:\__TEMP\plx2ntt\spikes_b0034_d180304_TT4.plx';
NTT_filename_out = 'D:\__TEMP\plx2ntt\spikes_b0034_d180304_TT4+.NTT';
NTT_filename_header = 'D:\__TEMP\plx2ntt\spikes_b0034_d180304_TT4.NTT';
plexon_plx2ntt(PLX_filename, NTT_filename_out, NTT_filename_header)

%%
clear;clc
rng(0);
sdf = randn(1,1e6);
% sdf1 = sdf + 1;
% sdf2 = sdf + 0.5;
sdf1 = randn(1,1e6) + 2;
sdf2 = randn(1,1e6) + 1;
sdf3 = sdf1-sdf2;
median(abs(sdf1))
median(abs(sdf2))
median(abs(sdf3))

%% try filter design
clear; clc
filter_params.passband  = [600 6000];
% filter_params.type = 'butter';
% filter_params.type = 'fir1';
filter_params.type = 'firpm';
for order = [8]
    filter_params.order = order;
    t_start_end = [39968960821 40164431783];
    file_IN = 'D:\__TEMP\filtering\CSC8.ncs';
    file_OUT = ['D:\__TEMP\filtering\' sprintf('CSC8_filtered_%d-%d_%s_order_%d.ncs',filter_params.passband, filter_params.type,filter_params.order)];
    Nlx_filter_CSC2(file_IN, file_OUT, t_start_end, filter_params);
end


%% play with linking show/hide properties across axes
clear;clc
figure
h(1) = subplot(1,2,1);
hold on
h1=plot(1:10, 'k')
h2=plot(2:11, 'r')
h(2) = subplot(1,2,2);
hold on
h3=plot(1:10, 'k')
h4=plot(2:11, 'r')

linkprop([h1 h3], 'Visible')
% h1.Visible = 'off'
plotbrowser('on')


%% 06/03/2019 - test artifacts in the logger from bespoon RF (-> it was flash write artifacts!!!!!)
clear;clc
exp_ID = 'b9861_d180603';
exp=exp_load_data(exp_ID,'details','path','pos');
fs_bsp = exp.pos.raw.fs;
ts_bsp = 1/fs_bsp;
M_avg = [];
M_abs_avg = [];
M_std = [];
for TT=2
    for ch=3
        % load raw data
        fprintf('TT%dch%d\n',TT,ch)
        filename = sprintf('spikes_%s_TT%d_ch%d.ncs',exp_ID, TT, ch);
        filename = fullfile(exp.path.spikes_raw, filename);
        [signal, ts, fs, params] = Nlx_csc_read(filename, []);
        bsp_ts = exp.pos.raw.ts_nlg_usec;
%         ti = exp_get_sessions_ti(exp_ID,'Sleep1');
%         [~,~,bsp_ts,~] = get_data_in_ti(bsp_ts, ti);
        trigger_IX = find_nearest_point(bsp_ts, ts);
%         trigger_IX = trigger_IX(1:2:end);
%         n = 6000;
%         trigger_IX = find_nearest_point(linspace(ti(1),ti(2),n), ts);
%         test_fs = n/(diff(ti)*1e-6)
        win_len = [0 1/fs_bsp];
        window_num_samples = round(win_len.*fs);
        [trig_signal] = trigger_signal_by_IX(signal, trigger_IX, window_num_samples);
        M_avg(TT,ch,:) = mean(trig_signal,1);
        M_abs_avg(TT,ch,:) = mean(abs(trig_signal),1);
        M_std(TT,ch,:) = std(trig_signal,1,1);
    end
end

%%
opt = 1;
switch opt
    case 1
        M = M_avg;
        thr = 0.23;
        M_str = 'averaged triggered on bsp time';
    case 2
        M = M_abs_avg;
        M_str = 'averaged absolute triggered on bsp time';
    case 3
        M = M_std;
        thr = 100;
        M_str = 'std triggered on bsp time';
end
for TT=2
    figure
    for ch=3
        subplot(2,2,ch);
        title(sprintf('ch %d',ch))
        t = linspace(win_len(1),win_len(2),size(M_avg,3));
        signal_STA = squeeze(M(TT,ch,:));
        signal_STA = smooth(signal_STA,10);
        findpeaks(signal_STA,t, 'MinPeakHeight', thr);
        [pks,locs,w,p] = findpeaks(signal_STA,t, 'MinPeakHeight', thr);
        xlabel('Time lag from Bespoon acquisition (sec)')
        ylabel('Voltage (uVolt)')
        h=suptitle({'neural channel (filtered>600Hz)';M_str;exp_ID});
        h.Interpreter = 'none';
    end
end

%%
maxlag = 2*round(fs/fs_bsp);
[c,lags] = xcorr(abs(signal), maxlag);
lags_in_sec = lags./fs;

%%
plot(lags_in_sec,c)
xlabel('Time (s)')

%%
filename_out = 'D:\__TEMP\New folder (6)\bsp_ts.ntt';
Timestamps = bsp_ts';
Samples = triang(32).*200-100;
Samples = repelem(Samples,1, 4, length(Timestamps));
% Samples = ones(32,4,length(Timestamps));
CellNumbers = zeros(size(Timestamps));
nlx_ntt_write(filename_out, Timestamps, Samples, CellNumbers, fs);

%% 19/03/2019 - explore same interneuron recorded for two days
clear
clc
cells(1) = load('L:\Analysis\Results\cells\FR_map\b2289_d180514_TT1_SS01_cell_FR_map.mat')
cells(2) = load('L:\Analysis\Results\cells\FR_map\b2289_d180515_TT1_SS01_cell_FR_map.mat')
map1 = [cells(1).FR_map(1).all.PSTH cells(1).FR_map(2).all.PSTH];
map2 = [cells(2).FR_map(1).all.PSTH cells(2).FR_map(2).all.PSTH];
corr(map1',map2','rows','pairwise')

dx = 0.2;
fs = 1/dx;
window = round(25*fs);
noverlap = round(window/2);
figure
hold on
% pnl = panel();
% pnl.pack(2,2);
for ii_cell = 1:2
    for ii_dir = 1:2
        map = cells(ii_cell).FR_map(ii_dir).all.PSTH;
        map(isnan((map))) = 0;
%         pnl(ii_cell,ii_dir).select();
        [Pxx,F] = pwelch(map, window, noverlap,window,fs);
        plot(F,log(abs(Pxx)))
    end
end


%%
clear
clc
index = @(a,b)(abs(a-b)./(a+b));
x = 0.5:0.1:20;
y = 0.5:0.1:20;
y=1;
[X,Y] = meshgrid(x,y);
Z = index(X,Y);
figure
mesh(X,Y,Z)


%%
clear
clc
file_IN = 'D:\__TEMP\New folder (18)\CSC14.ncs';
ts_limits = [74468891905 74470238429];
[signal, ts, fs, params] = Nlx_csc_read(file_IN, [ts_limits]);

% filter_params.type = 'highpassfir1';
% filter_params.passband  = [600];
% filter_params.order = 12;
% Wn   = filter_params.passband / (fs/2);
% b = fir1(filter_params.order, Wn,'high');
% % b = fir1(filter_params.order, Wn,'bandpass');
% a = 1;

stopfreq = 400;
passfreq = 700;
custom_filt_params = {'highpassfir',...
    'StopbandFrequency', stopfreq,...
    'PassbandFrequency', passfreq,...
    'StopbandAttenuation', 60,...
    'PassbandRipple', 1,...
    'SampleRate', fs};
custom_filt = designfilt(custom_filt_params{:});
order = filtord(custom_filt);
signal1 = filtfilt(custom_filt,signal);
file_OUT1 = ['D:\__TEMP\New folder (18)\' sprintf('CSC14_filtfilt_custom_order_%d_%d_%d.ncs', order, stopfreq, passfreq)];
nlx_csc_write(file_OUT1, signal1, ts, fs, params.header);

% file_OUT1 = ['D:\__TEMP\New folder (18)\' sprintf('CSC14_filtfilt_order_%d.ncs',filter_params.order)];
% signal1 = filtfilt(b,a,signal);
% nlx_csc_write(file_OUT1, signal1, ts, fs, params.header);

% file_OUT2 = ['D:\__TEMP\New folder (18)\' sprintf('CSC14_filter_order_%d.ncs',filter_params.order)];
% signal2 = filter(b,a,signal);
% nlx_csc_write(file_OUT2, signal2, ts, fs, params.header);
            
%% playing with brush
clear
clc
figure
hold on
h1 = plot(1:10,'r');
h2 = plot(11:20,'b');
%%
clc
h1.BrushData
h2.BrushData
h=gca

%% 29/04/2019 - make sure all cell details are valid
clear; clc
cells_t = DS_get_cells_summary();
for ii_cell = 1:height(cells_t)
    %%
    cell = cells_t(ii_cell,:);
    cell_ID1 = cell.cell_ID{1};
    cell_ID2 = sprintf('b%04d_d%s_TT%d_SS%02d', ...
        cell.bat,...
        datestr(cell.date,'yymmdd'),...
        cell.TT,...
        cell.unit);
    if ~strcmp(cell_ID1 , cell_ID2 )
        fprintf('problem with %s (%d)\n', cell_ID1 , ii_cell)
    end
end
disp('finished')


%% 30/04/2019 - play with conv/imfilter (test nan and edge effects)
clear
clc
% rng(0);
x = randn(1,100)+10;
x(50:60) = nan;
x(1:10) = x(1:10)+10;
x = [nan(1,10) x];
nanx = isnan(x);
IX    = 1:numel(x);
x(nanx) = interp1(IX(~nanx), x(~nanx), IX(nanx), 'nearest', 'extrap');
% x(1) = 10;
% win = gausswin(10,5);
win = fspecial('gaussian',20,1);
x2 = imfilter(x,win,'same','conv','symmetric');
% x2 = imfilter(x,win,'same','conv',4);

x(nanx) = nan;
x2(nanx) = nan;

% figure
cla
hold on
plot(x)
plot(x2)
xlim([1 length(x)])
% ylim([0 30])

%% 01/05/2019 - look at time spent across the population
clear
clc
cells_t = DS_get_cells_summary();
cells_t(~ismember(cells_t.bat, [79,9861,2289] ),:) = [];
cells = cellfun(@(c)(cell_load_data(c,'FR_map','details')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
FR_maps = [cells.FR_map];
M = [];
for ii_cell = 1:length(cells)
    M(ii_cell,1,:) = cells(ii_cell).FR_map(1).all.time_spent;
    M(ii_cell,2,:) = cells(ii_cell).FR_map(2).all.time_spent;
end

figure
x = cells(1).FR_map(1).all.bin_centers;
subplot(2,1,1); hold on
plot(x, squeeze(M(:,1,:))')
plot(x, squeeze(mean(M(:,1,:),1))','k','LineWidth',2)
subplot(2,1,2); hold on
plot(x,squeeze(M(:,2,:))');
plot(x, squeeze(mean(M(:,2,:),1))','k','LineWidth',2)

%% 01/05/2019 - test parfor behavior
clear
clc
parfor ii=1:5
    func1(ii)
%     disp(ii)
%     figure(ii)
%     plot(0:1,0:1)
%     text(0.5,0.5,num2str(ii));
%     pause(ii)
    close all
end

% function func1(ii)
% %     ii = rand(1,1);
%     disp(ii)
%     figure
%     plot(0:1,0:1)
%     text(0.5,0.5,num2str(ii));
%     pause(ii)
%     saveas(gcf,sprintf('fig%d',ii),'tif');
% end

%% 19/05/2019 - copy good cells to new folder
dir_IN = 'L:\Analysis\Results\cells\figures';
dir_OUT = 'L:\Analysis\Results\cells\figures\good_cluster_CA1';
tmpl_list = cells_t.cell_ID;
ext = 'tif';
util_copy_files_by_template(dir_IN, dir_OUT, tmpl_list, ext);

%% play with dictionary
c = containers.Map;
c(79) = [1 0 1];
c(34) = [1 0 ];
keys(c)
values(c)

%%
bats_colors_mapping = {
    79,     'r';
    148,    'g';
    34,     [1 0 0];
    2289,   'c';
    9861,   'm';
    }
M = containers.Map([bats_colors_mapping{:,1}],...
                    bats_colors_mapping(:,2));
M(34)
M(79)

%%
clear;clc
prm = PARAMS_GetAll();
bats = keys(prm.graphics.colors.bats);
figure
hold on
for ii_bat = 1:length(keys(prm.graphics.colors.bats))
    bat = bats{ii_bat};
    c = prm.graphics.colors.bats(bat)
    plot(randn(1,10),'Color', c,'LineWidth',2);
end
legend("bat "+bats);

%% check even/odd corr values
cells = cellfun(@(x)(cell_load_data(x,'details','FR_map')), cells_t.cell_ID, 'UniformOutput', false);
cells = [cells{:}];
sdf=[cells.FR_map];
sdf2=[sdf.corr_odd_even];
ccc = [sdf2.rho];
histogram(ccc)
sum(isnan(ccc))

%% paper fig 1 panel H - YZ small variabiltiy
exp_ID = 'b0034_d180313';
X0 = 1:200;
x0 = 164;
exp = exp_load_data(exp_ID,'details','pos','flight');
prm = PARAMS_GetAll();
FE = exp.flight.FE([exp.flight.FE.distance]>prm.flight.full_min_distance);
for ii_fe = 1:length(FE)
    FE(ii_fe).Y = interp1(exp.pos.proj.ts, exp.pos.proj.pos(:,2), FE(ii_fe).ts);
    FE(ii_fe).Z = interp1(exp.pos.raw.ts_nlg_usec, exp.pos.raw.pos(:,3), FE(ii_fe).ts);
end
FE_dir = FE([FE.direction]==1);
YZ0 = [];
for ii_flight = 1:length(FE_dir)
    flight = FE_dir(ii_flight);
    x = flight.pos;
    [C,ia,idx] = unique(x,'stable');
    y = accumarray(idx,flight.Y,[],@median);
    z = accumarray(idx,flight.Z,[],@median);
    Y0 = interp1(C,y,X0,'linear');
    Z0 = interp1(C,z,X0,'linear');
    YZ0(ii_flight,1,:) = Y0;
    YZ0(ii_flight,2,:) = Z0;
end
YZ0 = YZ0 - nanmean(YZ0,1);

%%
YZ_std = squeeze(nanstd(YZ0,1));
figure
yyaxis left
plot(X0,YZ_std(1,:))
ylabel('Y std (m)')
yyaxis right
plot(X0,YZ_std(2,:))
ylabel('Z std (m)')

%%
clear
clc
options = {
    'b0034_d180313', 1, 160;...
    'b0034_d180314', 1, 175;...
    'b0034_d180315', 1, 162;...
    'b0079_d160909', 1,  49;...
    'b0079_d160911', 1,  52;...
    'b0079_d160913', 1,  57;...
    'b0079_d160914', 1,  54;...
    'b0079_d160916', 1, 133;...
    'b0079_d160918', 1, 147;...
    'b0079_d160919', 1, 144;...
    'b0079_d160920', 1, 144;...
    'b0079_d160921', 1, 147;...
    'b0079_d160925', 1, 140;...
    'b0079_d160926', 1,  84;...
    'b0079_d160928', 1, 135;...
    'b0079_d160930', 1,  88;...
    'b0079_d161003', 1, 143;...
    'b0079_d161005', 1, 139;...
    'b2289_d180514', 1, 160;...
    'b2289_d180525', 1, 162;...
    'b2289_d180528', 1, 164;...
    'b2289_d180531', 1, 165;...
    'b9861_d180526', 1, 156;...
    'b9861_d180601', 1, 163;...
    'b9861_d180606', 1,  22;...
    'b9861_d180608', 1, 160;...
    'b9861_d180609', 1, 160;...
    'b9861_d180610', 1, 160;...
    'b9861_d180612', 1, 164;...
    'b9861_d180613', 1, 160;...
    'b9861_d180614', 1, 168;...
    'b9861_d180615', 1, 164;...
    'b9861_d180616', 1, 162;...
    'b9861_d180617', 1, 168;...
    'b9861_d180618', 1, 159;...
    'b9861_d180619', 1, 179;...
    'b9861_d180620', 1, 158;...
    'b9861_d180621', 1, 168;...
    'b9861_d180622', 1, 158;...
    'b9861_d180623', 1, 158;...
    'b9861_d180624', 1, 167;...
    'b9861_d180626', 1, 167;...
    'b9861_d180627', 1, 162;...
    'b9861_d180628', 1, 159;...
    'b9861_d180709', 1, 154;...
    'b9861_d180710', 1, 171;...
    'b9861_d180711', 1, 164;...
    };
for option = 1:length(options)
    %%
    exp_ID = options{option,1};
    direction = options{option,2};
    x0 = options{option,3};
    load("L:\Analysis\Results\exp\pos_XYZ\"+exp_ID+"_pos_XYZ.mat");
    YZ0 = XYZ.YZ_by_dir{direction};
    YZ0 = YZ0 - nanmean(YZ0,1);
    YZ0 = YZ0 + [0 1.5];
    figure
    hold on
    plot(YZ0(:,1,x0),YZ0(:,2,x0),'.');
    plot([-1.25 -1.25],[0 1.7],'k')
    plot([ 1.25  1.25],[0 1.7],'k')
    plot([ 1.25 0],[1.7 2.35],'k')
    plot([-1.25 0],[1.7 2.35],'k')
    xlim([-1.5 1.5])
    ylim([-0 2.5])
    xlabel('Y position (m)');
    ylabel('Z position (m)');
    n=sum(~any(isnan(YZ0(:,:,x0))'));
    text(0.9,0.9,"n="+n, 'Units','normalized','FontSize',14);
    title(exp_ID + " dir"+direction + " @x="+x0,'Interpreter','none');
    fig_filename = exp_ID + "_dir"+direction + "_x="+x0;
    fig_filename = fullfile('L:\Analysis\Results\exp\pos_XYZ\good_examples',fig_filename);
    saveas(gcf,fig_filename,'tif')
    close all
end

%% 
figure
set(groot, 'defaultAxesTickDir', 'out');
% set(groot,  'defaultAxesTickDirMode', 'manual');
% set(0,'defaultAxesTickDir','out');
plot(1:10)



%% fig 1 - find panel J example (speed trajectory)
prm = PARAMS_GetAll();
for ii_exp = 1:height(exp_t)
    %% get exp data
    ii_exp 
    exp_ID = exp_t.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'details','flight');
    FE = exp.flight.FE;
    FE = FE([FE.distance]>prm.flight.full_min_distance);
    
    %% plot
    figure('Units','centimeters','Position',[5 5 15 20])
    pnl=panel();
    pnl.pack('v',3);
    pnl.margin = [15 15 5 10 ];
    pnl.de.margin = 10;
    h=pnl.title(exp.details.exp_ID);
    h.Interpreter = 'none';
    h.FontSize = 16;
    h.Position = [0.5 1];
    
    pnl(1).select();
    hold on
    directions = [1 -1];
    for ii_dir = 1:2
        FE_dir = FE([FE.direction]==directions(ii_dir));
        plot([FE_dir.pos],[FE_dir.vel],'.', 'Color',prm.graphics.colors.flight_directions{ii_dir});
    end
    ylabel('velocity (m/s)')
    text(1,1.05,"n="+length(FE),'Units','normalized','HorizontalAlignment','right','FontSize',14);
    
    pnl(2).select();
    hold on
    plot([FE_dir.pos],abs([FE_dir.vel]),'.k','MarkerSize',1);
    ylabel('Speed (m/s)')
    ylim([0 10])
    
    pnl(3).select();
    hold on
    for ii_dir = [1 2] 
        c = prm.graphics.colors.flight_directions{ii_dir};
        ydev = exp.flight.pos_y_std(ii_dir);
        x = ydev.xy(:,1);
        y = ydev.xy(:,2);
        ymean = interp1(ydev.bin_centers, ydev.ymean, x);
        y = y-ymean;
        plot(x, y, '.', 'Color',c, 'MarkerSize',.0001);
    end
    xlabel('X position (m)')
    ylabel('Y (m)')
    
    figname = fullfile('L:\paper_figures\speed_trajectory', exp.details.exp_ID+"_speed_trajectory");
    saveas(gcf, figname, 'tif');
    close all
end

%%
figure
subplot(121)
plot(abs(field_vel_all) , field_size_all ,'.k')
xlabel('speed (m/s)')
ylabel('Size (m)')
title('speed at peak of place field')
subplot(122)
plot(abs(field_vel2_all) , field_size_all ,'.k')
xlabel('speed (m/s)')
ylabel('Size (m)')
title('speed averaged over in-field spikes')
linkaxes(findobj(gcf,'Type','Axe'))





%%

%% multi-scale: control for speed changes !
clc
cells_signif = cat(1,cells.signif);
cells_signif = arrayfun(@(x)(x.TF), cells_signif);
signif_IX = any(cells_signif,2);
stats = [cells(signif_IX).stats];
stats_all = [stats.all];

figure
clear h

subplot(131)
hold on
axis equal
x = abs([stats_all.field_smallest_vel]);
y = abs([stats_all.field_largest_vel]);
h(1) = plot(x,y, '.k');
refline(1,0)
xlabel('speed at smallest field (m/s)')
ylabel('speed at largest field (m/s)')
p_ranksum=ranksum(x,y);
text(0.1,0.95,sprintf('p ranksum=%.2f',p_ranksum),'Units','normalized')

subplot(132)
hold on
x = abs([stats_all.field_ratio_LS_vel]);
y = [stats_all.field_ratio_LS];
h(2) = plot(x,y, '.k');
xlabel('speed ratio')
ylabel('size ratio')
lm = fitlm(x,y)
[r,p] = corr(x',y','rows','pairwise')
text(0.1,0.95,sprintf('r=%.2f,p=%.3f',r,p),'Units','normalized')
text(0.1,0.90,sprintf('R^2=%.3f',r^2),'Units','normalized')

subplot(133)
hold on
x = abs([stats_all.field_ratio_LS_vel]);
y = [stats_all.field_ratio_LS];
h(3) = plot(x,y, '.k');
xlabel('speed ratio')
ylabel('size ratio (log)')
set(gca,'yscale','log')
set(gca,'ytick',[1 2 3 5 10 15 20])

hlink = linkprop(h,'BrushData');
setappdata(gcf, 'brush_data_link', hlink);




%% 14/11/2019
fields_size = [];
fields_peak_FR = [];
for ii_dir = 1:2
    for ii_cell = 1:length(cells)
        cell = cells(ii_cell);
        if ~cell.signif(ii_dir).TF % check signif per direction
            continue;
        end
        fields = cell.fields{ii_dir};
        fields([fields.in_low_speed_area]) = []; % remove fields in low speed area
        fields_size = [fields_size fields.width_prc];
        fields_peak_FR = [fields_peak_FR fields.peak];
    end
end
% fields_peak_FR = fields_size .* 2 + randn(size(fields_peak_FR));
N = length(fields_size);

figure

subplot(221)
plot(fields_size, fields_peak_FR, '.')
[r,p] = corr(fields_size', fields_peak_FR','type','Pearson');
text(0.98,0.95,sprintf('r=%.2f',r),'Units','normalized','HorizontalAlignment','right')
text(0.98,0.90,sprintf('p=%.2g',p),'Units','normalized','HorizontalAlignment','right')
xlabel('fields size (m)')
ylabel('fields peak FR (Hz)')
title('linear-linear')

subplot(222)
plot(fields_size, fields_peak_FR, '.')
ha=gca;
ha.YScale = 'log';
xlabel('fields size (m)')
ylabel('fields peak FR (Hz)')
title('log-linear')
[r,p] = corr(fields_size', log(fields_peak_FR'),'type','Pearson');
text(0.98,0.95,sprintf('r=%.2f',r),'Units','normalized','HorizontalAlignment','right')
text(0.98,0.90,sprintf('p=%.2g',p),'Units','normalized','HorizontalAlignment','right')

subplot(223)
plot(fields_size, fields_peak_FR, '.')
ha=gca;
ha.YScale = 'log';
ha.XScale = 'log';
xlabel('fields size (m)')
ylabel('fields peak FR (Hz)')
title('log-log')
[r,p] = corr(log(fields_size'), log(fields_peak_FR'),'type','Pearson');
text(0.98,0.95,sprintf('r=%.2f',r),'Units','normalized','HorizontalAlignment','right')
text(0.98,0.90,sprintf('p=%.2g',p),'Units','normalized','HorizontalAlignment','right')

subplot(224)
hold on
[~,IX1] = sort(fields_size);
[~,IX2] = sort(fields_peak_FR);
rank1 = 1:N;
rank1(IX1) = rank1;
rank2 = 1:N;
rank2(IX2) = rank2;
plot(rank1, rank2, '.', 'MarkerSize',4)
ksdensity([rank1;rank2]','PlotFcn','contour')
[r,p] = corr(rank1', rank2','type','Pearson');
text(0.98,0.95,sprintf('r=%.2f',r),'Units','normalized','HorizontalAlignment','right')
text(0.98,0.90,sprintf('p=%.2g',p),'Units','normalized','HorizontalAlignment','right')
xlabel('fields size (rank)')
ylabel('fields peak FR (rank)')
title('by rank')

h=suptitle('Field peak FR vs. size correlations');
h.FontSize = 16;

%%
L = 50000;
[f] = MultiscalePlace_GenerateTuning_AllModels(L);
%%
figure
sdf=f.fI;
[~,IX] = sort(sum(sdf));
ds = 20;
pos = 1:ds:L;
pos = pos / 100; % back to meter
subplot(121)
imagesc(pos, 1:1000, 1-sdf(:,IX)');
colormap bone
subplot(122)
mean(sum(sdf)./5)
histogram(sum(sdf)./5)

%%
clear
clc
file_in = 'L:\DATA\9861_Somo\20180704\nlx\CSC5.ncs';
limits_ts = [28121638315 28985664907];
[signal, ts, fs, params] = Nlx_csc_read(file_in, []);
limits_ts = ts(1) + round([0 0.5*60*1e6]);
[signal, ts, fs, params] = Nlx_csc_read(file_in, limits_ts);

% set up the parameters of the filter
ops.fshigh = 600;
ops.fs = fs;
[b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high');
nCh = 16;
signal = repmat(signal, nCh, 1);

dataRAW = gpuArray(signal); % move int16 data to GPU
dataRAW = dataRAW';
dataRAW = single(dataRAW); % convert to float32 so GPU operations are fast
% dataRAW = dataRAW(:, chanMap); % subsample only good channels

% subtract the mean from each channel
dataRAW = dataRAW - mean(dataRAW, 1); % subtract mean of each channel

% next four lines should be equivalent to filtfilt (which cannot be used because it requires float64)
tic
datr = filter(b1, a1, dataRAW); % causal forward filter
datr = flipud(datr); % reverse time
datr = filter(b1, a1, datr); % causal forward filter again
datr = flipud(datr); % reverse time back
toc

%%
clear
clc
switch 2
    case 1
        % file_in = 'L:\DATA\9861_Somo\20180704\nlx\CSC5.ncs';
        % [signal, ts, fs, params] = Nlx_csc_read(file_in, []);
        nCh = 16;
        signal = repmat(signal, nCh, 1);
    case 2
        fs = 32000;
        dt = 1/fs;
        T = 1*60*60;    % 1 hours
%         T = 10*60;      % 10 minutes
        t = 0:dt:T;
        L = length(t);
        nCh = 1;
%         signal = randn(nCh,L);
        f1=10;
        f2=1000;
        signal = sin(2*pi*f1.*t) + sin(2*pi*f2.*t);
        signal = repmat(signal, nCh, 1);
end

% set up the parameters of the filter
ops.fshigh = 600;
ops.fs = fs;
[b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high');

end_loop = 0;
chunk_size = round(2*60*fs);
nDataChunks = ceil(L/chunk_size)
chunk_start_IX = 1;
% ch_chunks = num2cell(1:nCh);  % channel by channel
ch_chunks = {1:nCh};            % all channels together
signal_filtered = zeros(size(signal));
timing_all = tic;
for ii_ch_chunk = 1:length(ch_chunks)
    ch_chunk_IX = ch_chunks{ii_ch_chunk};
    ii_data_chunk = 1;
    while ~end_loop
        fprintf('data chunk %d: %d/n', ii_data_chunk);
        chunk_end_IX = chunk_start_IX + chunk_size-1;
        if chunk_end_IX > L
            chunk_end_IX = L;
            end_loop = 1;
        end
        data_chunk_IX = chunk_start_IX:chunk_end_IX;

        dataRAW = gpuArray(signal(ch_chunk_IX,data_chunk_IX)); % move int16 data to GPU
        dataRAW = dataRAW';
        dataRAW = single(dataRAW); % convert to float32 so GPU operations are fast
        % dataRAW = dataRAW(:, chanMap); % subsample only good channels

        % subtract the mean from each channel
    %     dataRAW = dataRAW - mean(dataRAW, 1); % subtract mean of each channel

        % next four lines should be equivalent to filtfilt (which cannot be used because it requires float64)

        tic
        datr = filter(b1, a1, dataRAW); % causal forward filter
        datr = flipud(datr); % reverse time
        datr = filter(b1, a1, datr); % causal forward filter again
        datr = flipud(datr); % reverse time back
        signal_filtered(ch_chunk_IX,data_chunk_IX) = gather(datr)';
        toc
        
        chunk_start_IX = chunk_start_IX + chunk_size;
        ii_data_chunk = ii_data_chunk + 1;
        disp('...')
    end
end
toc(timing_all)

%%
figure
hold on
plot(signal(1,:));
plot(signal_filtered(1,:));

%%
% 64ch 0.5min 1.7322sec
% 16ch 0.5min 0.9354sec

%% 26/12/2019 - all kind of correlations
figure
hold on
h=[];
h(1)=plot(LS_field_ratio_all, LS_field_size(1,:),'*');
h(2)=plot(LS_field_ratio_all, LS_field_size(2,:),'+');
linkprop(h, 'BrushData');
xlabel('LS field ratio')
ylabel('LS field size (m)')
[r1,p1]=corr(LS_field_ratio_all', LS_field_size(1,:)', 'rows','pairwise');
[r2,p2]=corr(LS_field_ratio_all', LS_field_size(2,:)', 'rows','pairwise');
% text(0.7,0.9, sprintf('S: r=%.2f, p=%.2g',r1,p1), 'Units','normalized');
% text(0.7,0.8, sprintf('L: r=%.2f, p=%.2g',r2,p2), 'Units','normalized');
legend({sprintf('S: r=%.2f, p=%.2g',r1,p1);...
        sprintf('L: r=%.2f, p=%.2g',r2,p2)}, 'Location','northwest')
h=gca;
h.XScale = 'log';
h.YScale = 'log';

%%
figure
plot(LS_field_size(1,:), LS_field_size(2,:),'.')
lsline
axis equal

%%
smallest_field_peak_FR = nan(1,length(cells));
for ii_cell = 1:length(cells)
    cell = cells(ii_cell);
    % pooled stats - check at least one direction is signif
    if any([cell.signif.TF])
        LS_field_ratio_all(ii_cell) = cell.stats.all.field_ratio_LS;
    end
    % per dir stats - check signif per direction
    for ii_dir = 1:2
        if cell.signif(ii_dir).TF 
            LS_field_ratio_dir(ii_dir,ii_cell) = cell.stats.dir(ii_dir).field_ratio_LS;
        end
    end
end

%%
figure
subplot(2,2,1)
plot(in_field_spikes_prc(:),total_area(:),'.')
subplot(2,2,3)
sdf = repelem(LS_field_ratio_all,2,1)';
plot(in_field_spikes_prc(:),sdf(:),'.')
subplot(2,2,4)
plot(in_field_spikes_prc(:),LS_field_ratio_dir(:),'.')

%% 29/12/2019 - simulate perfect gaussian field 
% and see how many spikes are counted as in-field (with the same detections
% process and params)
clear
clc
T = 60*60; % in sec
dt=0.01;
v = 8;
dx = v*dt
x = 0:dx:200;
field_pos = 100;
field_size_std = 5;
FR_map = gaussmf(x,[field_size_std field_pos]);
FR_map_peak = 20;
FR_map = FR_map .* FR_map_peak;
field_href = 0.2;
n = round(T/dt)
pos_IX = randi([1 length(x)], 1,n);
pos = x(pos_IX);
flight = randi([1 40], 1,n);
spikes_rate = poissrnd(FR_map(pos_IX).*dt);
unique(spikes_rate) % dt is so small -> binary case (spike/no spike)
spikes_pos = pos;
spikes_flight = flight;
spikes_pos(spikes_rate==0) = [];
spikes_flight(spikes_rate==0) = [];
left_border_IX = find( FR_map < field_href*FR_map_peak & x < field_pos, 1,'last');
right_border_IX = find( FR_map < field_href*FR_map_peak & x > field_pos, 1,'first');
field_border_IX = [left_border_IX right_border_IX];
field_border_pos = x(field_border_IX);
in_field_spikes = spikes_pos > field_border_pos(1) & spikes_pos < field_border_pos(2);
in_field_spikes_prc = sum(in_field_spikes) / length(in_field_spikes)
out_of_field_spikes_prc = 1-in_field_spikes_prc
% 
h=[];
figure
h(1)=subplot(2,1,1);
hold on
plot(x,FR_map,'k');
xline(x(field_border_IX(1)));
xline(x(field_border_IX(2)));
yline(field_href*FR_map_peak);
plot(x(field_border_IX), FR_map(field_border_IX), 'ob');
h(2)=subplot(2,1,2);
hold on
% plot(pos,flight,'.k')
plot(spikes_pos,spikes_flight,'.r')
plot(spikes_pos(~in_field_spikes),spikes_flight(~in_field_spikes),'.b')
linkaxes(h,'x')

text(1,0.8,sprintf('in field spikes=%.1f%%',100*in_field_spikes_prc),'Units','normalized','HorizontalAlignment','right')

%% 30/12/2019 - test symlink (for the different paramsets)
clear
clc
dir_link     = 'D:\__TEMP\symlink_test\cells\fields';
dir_original = 'D:\__TEMP\symlink_test\cells_paramset_1\fields';
command = sprintf('mklink /D "%s" "%s"', dir_link, dir_original);
delete(dir_link)
[status,cmdout] = system(command,'-echo')

%% 30/12/2019 - Run different paramsets
clear
clc
for paramset = 0:8
    PARAMS_SetParamset(paramset);
    PIP_cell_analyses;
end

%% 06/01/2020 - check how many cells have PF ONLY near balls, and thus is not considered to be signif
flags = zeros(length(cells),2);
for ii_cell = 1:length(cells)
    cell = cells(ii_cell);
    for ii_dir = 1:2
        if ~cell.signif(ii_dir).TF
            fields = cell.fields{ii_dir};
            if ~isempty(fields)
                if all([fields.in_low_speed_area])
                    flags(ii_cell,ii_dir) = 1;
                end
            end
        end
    end
end
sum(flags)
[rows,cols]=find(flags);
arrayfun(@(x)(x.cell_ID),[cells(rows).details], 'UniformOutput',0)'
cols    

signif_TF = arrayfun(@(x)(x.TF), signif);
cells_excluded = ~any(signif_TF,2) & any(flags,2);
sum(cells_excluded)
arrayfun(@(x)(x.cell_ID),[cells(cells_excluded).details], 'UniformOutput',0)'

%%
nFE1 = nan(length(cells1),2);
nFE2 = nan(length(cells2),2);
for ii_cell = 1:length(cell1)
    cell1 = cell_load_data(cells1(ii_cell).details.cell_ID, 'FE');
    cell2 = cell_load_data(cells2(ii_cell).details.cell_ID, 'FE');
    nFE1(ii_cell,:) = cellfun(@length, cell1.FE);
    nFE2(ii_cell,:) = cellfun(@length, cell2.FE);
end
plot(nFE1, nFE2, '.')

%% 06/10/2020 - check what happened to Fig. 2K after introducing paramsets
clc
cells1 = pop_data(1).cells; % paramset 0
cells2 = pop_data(2).cells; % paramset 100
signif1 = any(arrayfun(@(x)(x.TF), cat(1,cells1.signif)),2);
signif2 = any(arrayfun(@(x)(x.TF), cat(1,cells2.signif)),2);
if any(signif1 ~= signif2)
    error('different signif cells between paramsets')
end
cells1(~signif1)=[];
cells2(~signif2)=[];
stats1 = [cells1.stats];
stats2 = [cells2.stats];
stats1 = [stats1.all];
stats2 = [stats2.all];
nFEs1 = cellfun(@length, cat(1,cells1.FE));
nFEs2 = cellfun(@length, cat(1,cells2.FE));

figure
subplot(2,3,1)
x1 = [stats1.field_num];
x2 = [stats2.field_num];
stem(x2-x1);
ndiff = length(find(x2~=x1));
xlabel('cells');
ylabel('{\Delta} no. of fields');
title(sprintf('#changes=%d',ndiff))

subplot(2,3,2); hold on
axis equal
x1 = [stats1.field_ratio_LS];
x2 = [stats2.field_ratio_LS];
plot(x1,x2, 'ok')
plot(x1(x1~=x2),x2(x1~=x2), '*r')
ndiff = length(find(x2~=x1 & ~isnan(x1) & ~isnan(x2)));
xlabel('field LS ratio NEW');
ylabel('field LS ratio OLD');
title(sprintf('#changes=%d',ndiff))
xlim([min([x1,x2]) max([x1,x2])])
ylim([min([x1,x2]) max([x1,x2])])
refline(1,0)

subplot(2,3,3); hold on
axis equal
x1 = abs([stats1.field_ratio_LS_vel]);
x2 = abs([stats2.field_ratio_LS_vel]);
plot(x1,x2, '.k')
plot(x1(x1~=x2),x2(x1~=x2), '.r')
ndiff = length(find(x2~=x1 & ~isnan(x1) & ~isnan(x2)));
xlabel('field LS speed ratio NEW');
ylabel('field LS speed ratio OLD');
title(sprintf('#changes=%d',ndiff))
xlim([min([x1,x2]) max([x1,x2])])
ylim([min([x1,x2]) max([x1,x2])])
refline(1,0)

subplot(2,3,4); hold on
axis equal
refline(1,0)
x1 = abs(speed_ratio_max_min1);
x2 = abs([stats2.field_ratio_LS_vel]);
plot(x1,x2, '.k')
plot(x1(x1~=x2),x2(x1~=x2), '.r')
ndiff = length(find(x2~=x1 & ~isnan(x1) & ~isnan(x2)));
xlabel('field LS speed ratio NEW (calced with OLD bug)');
ylabel('field LS speed ratio OLD');
title(sprintf('#changes=%d',ndiff))
xlim([min([x1,x2]) max([x1,x2])])
ylim([min([x1,x2]) max([x1,x2])])
refline(1,0)

subplot(2,3,5); hold on
axis equal
refline(1,0)
x1 = abs(speed_ratio_max_min1);
x2 = abs(speed_ratio_max_min2);
plot(x1,x2, '.k')
plot(x1(x1~=x2),x2(x1~=x2), '.r')
ndiff = length(find(x2~=x1 & ~isnan(x1) & ~isnan(x2)));
xlabel('field LS speed ratio NEW (calced with OLD bug)');
ylabel('field LS speed ratio OLD (calced with OLD bug)');
title(sprintf('#changes=%d',ndiff))
xlim([min([x1,x2]) max([x1,x2])])
ylim([min([x1,x2]) max([x1,x2])])
refline(1,0)

subplot(2,3,6); hold on
axis equal
x1 = nFEs1;
x2 = nFEs2;
plot(x1(:),x2(:), '.k')
plot(x1(x1~=x2),x2(x1~=x2), '.r')
ndiff = length(find(any(x2~=x1 & ~isnan(x1) & ~isnan(x2))));
xlabel('no. of flights NEW');
ylabel('no. of flights OLD');
title(sprintf('#changes=%d',ndiff))
xlim([min([x1(:);x2(:)]) max([x1(:);x2(:)])])
ylim([min([x1(:);x2(:)]) max([x1(:);x2(:)])])
refline(1,0)

saveas(gcf, 'L:\paper_figures\20200106_compare_BEFORE_AFTER_field_in_low_speed_correction\backwards_comparison', 'fig')
saveas(gcf, 'L:\paper_figures\20200106_compare_BEFORE_AFTER_field_in_low_speed_correction\backwards_comparison', 'tif')

%%
clc
ii_cell = 10
cells1(ii_cell).stats.all.field_ratio_LS_vel
cells2(ii_cell).stats.all.field_ratio_LS_vel

[cells1(ii_cell).fields{1}.vel]
[cells2(ii_cell).fields{1}.vel]

%%
cells = cells2;
speed_ratio_max_min = nan(size(cells));
for ii_cell = 1:length(cells)
    cell = cells(ii_cell);
    fields = cell.fields;
    fields_all = [];
    for ii_dir = 1:2
        if isempty(fields{ii_dir})
            continue
        end
        fields_to_add = fields{ii_dir};
        % workaround to solve the problem that sometimes I don't have the
        % field 'overlap_edges'... maybe change that in 'cell_calc_fields'...
        if isfield(fields_to_add,'overlap_edges')
            fields_to_add = rmfield(fields_to_add,'overlap_edges');
        end
        fields_all = [fields_all fields_to_add];
    end
    fields_all([fields_all.in_low_speed_area])=[];
    speed_ratio_max_min(ii_cell) = max([fields_all.vel]) / min([fields_all.vel]);
end

%% Yohai code example
% spikes_ts; % spikes time in ms
% FR_bin_size = 1; % in ms
% ker_SD = 15; % in ms
% sleep_ts; % start/end time of single sleep session in ms
% FR_bin_edges = sleep_ts(1) : FR_bin_size : sleep_ts(2);
% N=histcounts(spikes_ts, FR_bin_edges);
% FR_sleep = N * (1e3/FR_bin_size); % in Hz (Actually it is not really neccessary to convert to Hz because we use z-score)
% % smooth
% hsize = 1 + (5*ker_SD/FR_bin_size);
% hsize = round(hsize);
% alpha = hsize*FR_bin_size/(2*ker_SD); %the equation is:
% % alpha = hsize/2*sigma  where sigma is in number of samples - so the sigma_PF_smoothing
% % (m) is devided by bin size (m) in order to get the sigma in number of bins unit.
% ker = gausswin(hsize,alpha)'./(sqrt(2*pi)*ker_SD);
% FR_sleep_smoothed = imfilter(FR_sleep,ker,'same','conv','symmetric');
% FR_sleep_zscored = zscore(FR_sleep_smoothed);
% findpeaks(FR_sleep_zscored, 'MinPeakHeight', 3);

%%
cell_ID = 'b0148_d170625_TT4_SS02';
exp_ID = cell.details.exp_ID;
exp=exp_load_data(exp_ID);

%%
figure
spikes_size = squeeze(range(cell.spikes.waveforms,1));
% x = squeeze(cell.spikes.waveforms(8,2,:));
% y = squeeze(cell.spikes.waveforms(8,3,:));
% z = squeeze(cell.spikes.waveforms(8,4,:));
x = spikes_size(2,:);
y = spikes_size(3,:);
z = spikes_size(4,:);
plot3(x,y,z,'.')

%%
IX = find(z>500);
length(IX)
bad_spikes_ts = cell.spikes.ts(IX);
[IX, IX_per_ti, t2, t2_per_ti] = get_data_in_ti(bad_spikes_ts, exp.details.session_ts)
FE_dir1_ts = [cell.FE{1}.start_ts; cell.FE{1}.end_ts]';
FE_dir2_ts = [cell.FE{2}.start_ts; cell.FE{2}.end_ts]';
IX1 = get_data_in_ti(bad_spikes_ts, FE_dir1_ts);
IX2 = get_data_in_ti(bad_spikes_ts, FE_dir2_ts);

%% check what nan do to histogram plot
figure
x = [1 2 3 4 5 nan nan nan];
subplot(211)
histogram(x,'Normalization','pdf')
subplot(212)
histogram(x(~isnan(x)),'Normalization','pdf')


%% check de facto numbers...
stats= [cells.stats];
stats_dir = cat(1,stats.dir);
sdf=cat(1,stats_dir(signif_cells_IX).spikes_num_air);

%% #spikes in-air
sdf=arrayfun(@(x)(x.spikes_num_air), stats_dir);
IX = sdf<50 & ~signif_cells_IX;

%% #spikes in-air (only valid speed zone)
nSpikeHighSpeed = zeros(length(cells),2);
nSpikeInAir = zeros(length(cells),2);
for ii_cell = 1:length(cells)
    for ii_dir = 1:2
        FE = cells(ii_cell).FE{ii_dir};
        nSpikeHighSpeed(ii_cell,ii_dir) = length(get_data_in_ti( [FE.spikes_pos], prm.fields.valid_speed_pos));
        nSpikeInAir(ii_cell,ii_dir) = length( [FE.spikes_pos] );
    end
end
nSpikeHighSpeed(~signif_cells_IX) = nan;
nSpikeInAir(~signif_cells_IX) = nan;

%% Run code with the different paramsets
paramsets = 1:9;
for ii_paramset = 1:length(paramsets)
    ii_paramset
    paramsets = 1:9;
    paramset = paramsets(ii_paramset);
    fprintf('loading paramset %d results...\n',paramset)
    PARAMS_SetParamset(paramset);
    PIP_cell_analyses
    close all
end
% go back to default paramset (0)
disp('Go back to default paramset (0)')
PARAMS_SetParamset(0);

%% 3D plot of tunnel calib
load('L:\DATA\9861_Somo\calib\20180606_calib_middle_line_top_14_anchors_active\bsp\client\bsp_pos_tag_708.mat');
figure
plot3(bsp_pos.pos(:,1),bsp_pos.pos(:,2),bsp_pos.pos(:,3),'.')
% axis equal

%% 27/02/2020 - calc arms lengths and angle between them
exp_ID = 'b2289_d180615';
exp = exp_load_data(exp_ID, 'pos','LM');
turnpoint_LM = exp.LM(contains({exp.LM.name},{'turn-point'}));
% turnpoint_LM.pos
ball1_LM = exp.LM(2);
ball2_LM = exp.LM(end);

%%
figure
hold on
plot(exp.pos.calib_tunnel.curvexy(:,1), exp.pos.calib_tunnel.curvexy(:,2), '.')
plot(turnpoint_LM.pos_X, turnpoint_LM.pos_Y, 'or')
plot(turnpoint_LM.pos_X+[-50 0], turnpoint_LM.pos_Y+[-50 0], '--r')
axis equal
lm1=fitlm(xy1(:,1),xy1(:,2));
lm2=fitlm(xy2(:,1),xy2(:,2));
plot(lm1)
plot(lm2)

%% calc tunnel total length
xy = exp.pos.calib_tunnel.curvexy;
tunnel_length = sum( sqrt(sum(diff(xy).^2,2)) )
arm1_length   = sum( sqrt(sum(diff(xy1).^2,2)) )
arm2_length   = sum( sqrt(sum(diff(xy2).^2,2)) )
m1 = lm1.Coefficients.Estimate(2);
m2 = lm2.Coefficients.Estimate(2);
phi = rad2deg( atan(m1) - atan(m2) )

%% check
figure
hold on
h=histogram(LS_field_ratio_all);
h.DisplayStyle = 'stairs';
h=histogram(LS_field_ratio_all2);
h.DisplayStyle = 'stairs';
legend({'original';'corrected'})

%%
clc
LS_field_ratio_all(LS_field_ratio_all ~= LS_field_ratio_all2)
LS_field_ratio_all2(LS_field_ratio_all ~= LS_field_ratio_all2)
sum(~isnan(LS_field_ratio_all(LS_field_ratio_all ~= LS_field_ratio_all2)))


%% 
cells = cellfun(@(c)(cell_load_data(c,'details','stats','meanFR','signif','fields','FR_map','cluster_quality')), cells_ID, 'UniformOutput',0);
signif = arrayfun(@(x)(x.TF), cat(1,cells.signif));

%%
corr_even_odd = arrayfun(@(c)([c.stats.dir.corr_odd_even]), cells,'UniformOutput',0);
corr_even_odd = cat(1,corr_even_odd{:});
corr_begin_end = arrayfun(@(c)([c.stats.dir.corr_begin_end]), cells,'UniformOutput',0);
corr_begin_end = cat(1,corr_begin_end{:});
corr_bins = linspace(-1,1,50);
figure
subplot(211)
hold on
% histogram(corr_even_odd(signif),corr_bins )
histogram(corr_even_odd(signif(:,1)),corr_bins )
histogram(corr_even_odd(signif(:,2)),corr_bins )
title('corr even odd')
subplot(212)
hold on
% histogram(corr_begin_end(signif),corr_bins )
histogram(corr_begin_end(signif(:,1)),corr_bins )
histogram(corr_begin_end(signif(:,2)),corr_bins )
title('corr begin end')

%% 07/04/2020 - copy cells figure for cells with large fields (> 10 m)
clear
clc
load('L:\Analysis\Results\cells_paramset_0\cells_data.mat')
%% find relevant cells
field_size_thr = 10;
signif=cat(1,cells.signif);
signif=arrayfun(@(x)(x.TF),signif);
cells2copy = zeros(1,length(cells));
for ii_cell = 1:length(cells)
    cell = cells(ii_cell);
    for ii_dir = 1:2
        if cell.signif(ii_dir).TF
            fields = cell.fields{ii_dir};
            if any([fields.width_prc] > field_size_thr )
                cells2copy(ii_cell) = 1;
            end
        end
    end
end
cells_details = [cells.details];
cells_ID = {cells_details.cell_ID}';
cells2copy_ID = cells_ID( find(cells2copy ));
%% copy
dir_IN = 'L:\Analysis\Results\cells_paramset_0\figures';
dir_OUT = '';
util_copy_files_by_template(dir_IN, dir_OUT, cells2copy_ID, 'tif')

%% 29/11/2020 - debugging the gamma fit for compartmentalization model
figure
hold on
xxx = linspace(0,30,100);
plot(xxx, gampdf(xxx, 3.84, 0.85) ); % smallest
% plot(xxx, gampdf(xxx, 3.56, 1.3 ) ); % 6m - Yonatan's truncation @ 60cm (also maybe 3.56 should be 3.16 as in 200m)
% plot(xxx, gampdf(xxx, 3.75, 0.39 ) ); % 6m - normal gamfit
plot(xxx, gampdf(xxx, 3.56, 0.37 ) ); % 6m - Yonatan correction
plot(xxx, gampdf(xxx, 3.16, 1.8 ) ); % 200m
legend({'Smallest per cell in 200m','6m','200m'})
xlabel('Field size (m)')
ylabel('pdf')

%% fitting fields counts distribution with a poisson
figure
hold on
h=histogram(x);
h.Normalization = 'pdf';
xxx=1:10;
plot(xxx,poisspdf(xxx,poissfit(x)),'.-r')
% plot(xxx,poisspdf(xxx,1.5),'.-c')
xlabel('no. of fields')
ylabel('pdf')
legend('data','fit')

%%
ccc = [];
for ii_cell=1:length(cells)
    cell = cells(ii_cell);
    for ii_dir = 1:2
        if ~cell.signif(ii_dir).TF
            continue;
        end
        sdf = cell.FR_map(ii_dir).corr_odd_even;
        ccc=[ccc sdf.rho];
    end
end
figure
h=histogram(ccc);
h.BinEdges = linspace(-1,1,10000);
h.Normalization='cumcount'
xlim([-0.2 1])
ylim([0 25])
xlim([-1 1])
ylim([0 350])
xlabel('correlation')
ylabel('Cummulative number of cells')

%%
details=[cells.details];
IX1=r([stats_dir.corr_odd_even] < 0.2);
IX2=c([stats_dir.corr_odd_even] < 0.2);
cell2plot = sort({details(IX1).cell_ID}')
[details(IX1).cell_num]
IX2

%%
for ii_cell = 1:length(cell2plot)
    cell_ID = cell2plot{ii_cell};
    cell_plot_map_fields(cell_ID)
end

%% debug the low correlations in good maps!
cell_num = 473;
ii_dir = 1;
ii_cell = find([details.cell_num] == cell_num );
figure
hold on
cell = cells(ii_cell);
xxx = cell.FR_map(ii_dir).odd.bin_centers;
PSTH1 = cell.FR_map(ii_dir).odd.PSTH;
PSTH2 = cell.FR_map(ii_dir).even.PSTH;
ker = fspecial('gaussian',[1 10],4);
corr(PSTH1',PSTH2','rows','complete')
PSTH1 = imfilter(PSTH1,ker);
PSTH2 = imfilter(PSTH2,ker);
corr(PSTH1',PSTH2','rows','complete')
plot(xxx, PSTH1)
plot(xxx, PSTH2)


%%
sdf=arrayfun(@(cell)([cell.signif.has_min_spikes]), cells , 'UniformOutput', false);
sdf = cat(1,sdf{:});
active_cells = any(sdf');
sdf=arrayfun(@(cell)([cell.signif.TF]), cells , 'UniformOutput', false);
sdf = cat(1,sdf{:});
signif_cells = any(sdf');
active_non_signif_cells = active_cells & ~signif_cells;
fprintf('num active cells = %d\n',sum(active_cells))
fprintf('num signif cells = %d\n',sum(signif_cells))
active_non_signif_cells_ID = arrayfun(@(cell)([cell.details.cell_ID]), cells(active_non_signif_cells), 'UniformOutput', false);
fprintf('%s\n',active_non_signif_cells_ID{:})


%% fig S6 batch
clearvars -except data
% clear
clc
for grp=-1:-1:-4
    fprintf('grp %d\n',grp);
    paper_supp_fig_many_FR_map_examples_2
end






%%






%%
