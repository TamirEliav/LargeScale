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



%%



