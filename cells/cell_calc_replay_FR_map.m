function cell_calc_replay_FR_map(cell_ID,opts)
arguments
    cell_ID
    opts.plot_figs = true
end

%%
% clear
% clc
% close all
% cell_ID = 1017; %% good!
% cell_ID = 1074; %% good!
% cell_ID = 1138; %% very good!
% cell_ID = 1140; %% very good!
% cell_ID = 1141; %% very good!
% cell_ID = 1208; %% good!
% cell_ID = 1289; %% good!
% cell_ID = 1297; %% very good! (multiple fields!)
% cell_ID = 1299; % good
% cell_ID = 1302; % good - multiple fileds
% cell_ID = 1313; % good
% cell_ID = 1323; % good?
% cell_ID = 1337; % good?
% cell_ID = 1342; % good?

%%
cell_analysis_name = 'replay_FR_map';
dir_out = fullfile('L:\Analysis\Results\cells\',cell_analysis_name);
mkdir(dir_out);

%% load cell/exp data
cell = cell_load_data(cell_ID,'details','spikes','FR_map','fields');
exp_ID = cell.details.exp_ID;
exp = exp_load_data(exp_ID,'details','rest');

%% load replay data
epoch_types = {'sleep','rest'};
params_opt = 11;
event_type = 'posterior';
events_all = builtin('cell',[1 length(epoch_types)]);
for ii_epoch_type = 1:length(epoch_types)
    epoch_type = epoch_types{ii_epoch_type};
    [events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
    [~, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
    events(~TF) = [];
    if isempty(events)
        continue;
    end
    if isfield(events,'rest_ball_num')
        [events.rest_ball_loc] = disperse(exp.rest.balls_loc([events.rest_ball_num]));
    else
        [events.rest_ball_num] = disperse(nan(size(events)));
        [events.rest_ball_loc] = disperse(nan(size(events)));
    end
    events_all{ii_epoch_type} = events;
end

%% arrange replay events into categories
events_cat_names = {'sleep','rest','pooled','sleep1','sleep2','rest_ball_1','rest_ball_2'};
nCat = length(events_cat_names);
events_all_per_cat = builtin('cell',[1 length(events_cat_names)]);
events_all_per_cat{1} = events_all{1};
events_all_per_cat{2} = events_all{2};
events_all_per_cat{3} = [events_all{:}];
if ~isempty(events_all{1})
    events_all_per_cat{4} = events_all{1}([events_all{1}.epoch_num]==1);
    events_all_per_cat{5} = events_all{1}([events_all{1}.epoch_num]==2);
end
if ~isempty(events_all{2})
    events_all_per_cat{6} = events_all{2}([events_all{2}.rest_ball_num]==1);
    events_all_per_cat{7} = events_all{2}([events_all{2}.rest_ball_num]==2);
end

%% params
prm = PARAMS_GetAll();
pos_fs = 1/0.002; % decoder uses 2 ms time bins
pos_fs = 5000; % upsample do handle the high speed of replay
pos_dt = 1/pos_fs;
bin_size = prm.FR_map.bin_size;
bin_limits = prm.FR_map.bin_limits;
bin_edges = bin_limits(1):bin_size:bin_limits(end);
bin_centers = (bin_edges(1:end-1)+bin_edges(2:end))./2;
% min_time_spent_per_meter = prm.FR_map.min_time_spent_per_meter/50;
% min_time_spent_per_bin = min_time_spent_per_meter .* bin_size;
min_time_spent_per_meter = 0;
min_time_spent_per_bin = 0;
% min_time_spent_per_bin = 1/pos_fs;
% min_time_spent_per_bin = min_time_spent_per_bin*3;
min_time_spent_per_meter = min_time_spent_per_bin / bin_size;
ker_SD = prm.FR_map.ker_SD;

%%
dir2state_mapping = containers.Map([1 2],[4 1]);
nDir = length(dir2state_mapping.keys);
replay_PSTH_all = builtin('cell',[nCat,nDir]);
clear replay_PSTH_all
for ii_event_cat = 1:nCat
    events_cat = events_all_per_cat{ii_event_cat};
    for ii_dir = 1:nDir
        if ~isempty(events_cat)
            state_num = dir2state_mapping(ii_dir);
            IX = [events_cat.state_num]==state_num;
            events = events_cat(IX);
        else
            events = [];
        end
        if isempty(events)
            seqs = [];
            pos = [];
            spikes_pos = [];
            spikes_replay_num = [];
        else
            seqs = [events.seq_model];
            seqs_ti = [seqs.start_ts;seqs.end_ts];
            seqs_xi = [seqs.start_pos;seqs.end_pos];
            k = floor( [seqs.duration].*pos_fs );
            kmax = max(k);
            t = 0:pos_dt:(pos_dt*(kmax-1));
            t = t.*1e6;
            M = t' + seqs_ti(1,:);
            mask = [1:kmax]'.*ones(1,length(seqs));
            mask = logical(mask <= k);
            pos_ts = M(mask)';
            [~,~,spikes_ts,~] = get_data_in_ti(cell.spikes.ts, seqs_ti');
            pos = interp1(seqs_ti(:), seqs_xi(:), pos_ts, "linear","extrap");
            spikes_pos = interp1(seqs_ti(:), seqs_xi(:), spikes_ts, "linear","extrap");
            sdf = repmat(1:size(seqs_ti,2),[2 1]);
            spikes_replay_num = interp1(seqs_ti(:), sdf(:), spikes_ts,'previous');
        end

        [PSTH,spike_density,time_spent] = computePSTH(pos,pos_fs,spikes_pos,bin_edges,min_time_spent_per_bin,ker_SD);
        [SI_bits_spike, SI_bits_sec] = computeSI(PSTH,time_spent);
        sparsity = cumputeSparsity(PSTH,time_spent);
        [rho_pearson,pval_pearson] = corr(cell.FR_map(ii_dir).all.PSTH', PSTH', "rows","pairwise",'type','Pearson');
        [rho_spearman,pval_spearman] = corr(cell.FR_map(ii_dir).all.PSTH', PSTH', "rows","pairwise",'type','Spearman');

        valid_coverage_TF = ~isnan(PSTH);

        % arrange into a struct
        replay_PSTH = struct();
        replay_PSTH.PSTH = PSTH;
        replay_PSTH.spike_density = spike_density;
        replay_PSTH.time_spent = time_spent;
        replay_PSTH.min_time_spent_per_meter = min_time_spent_per_meter;
        replay_PSTH.bin_size = bin_size;
        replay_PSTH.bin_edges = bin_edges;
        replay_PSTH.bin_centers = bin_centers;
        replay_PSTH.SI_bits_spike = SI_bits_spike;
        replay_PSTH.SI_bits_sec = SI_bits_sec;
        replay_PSTH.sparsity  = sparsity;
        replay_PSTH.cat_name = events_cat_names{ii_event_cat};
        replay_PSTH.map_dir = ii_dir;
        replay_PSTH.nSpikes = length(spikes_pos);
        replay_PSTH.nReplays = length(events);
        replay_PSTH.replayTotalDuration = length(pos)/pos_fs;
        replay_PSTH.valid_env_coverage_m   = sum(valid_coverage_TF)*bin_size;
        replay_PSTH.valid_env_coverage_prc = sum(valid_coverage_TF)*bin_size / diff(exp.rest.balls_loc);
        replay_PSTH.events = events;
        replay_PSTH.spikes_pos = spikes_pos;
        replay_PSTH.spikes_replay_num = spikes_replay_num;
        replay_PSTH.maps_corr_rho_pearson = rho_pearson;
        replay_PSTH.maps_corr_pval_pearson = pval_pearson;
        replay_PSTH.maps_corr_rho_spearman= rho_spearman;
        replay_PSTH.maps_corr_pval_spearman= pval_spearman;
        
        fields = cell.fields{ii_dir};
        if ~isempty(fields)
            fields_xi = cat(1,fields.edges_prc);
            n_valid_bins_in_fields = length(get_data_in_ti(bin_centers(valid_coverage_TF), fields_xi));
            fields_total_area = sum(diff(fields_xi,1,2));
            replay_PSTH.valid_fields_coverage_m   = n_valid_bins_in_fields * bin_size;
            replay_PSTH.valid_fields_coverage_prc = n_valid_bins_in_fields * bin_size / fields_total_area;
        else
            replay_PSTH.valid_fields_coverage_m   = 0;
            replay_PSTH.valid_fields_coverage_prc = 0;
        end
        replay_PSTH_all(ii_event_cat,ii_dir) = replay_PSTH;
    end
end

%%
if opts.plot_figs
cat_to_plot_options = {[3];[1:3];[1:7]};
cat_to_plot_options_names = {'pooled';'sleep_rest';'extended'};
for ii_cat_to_plot = 1:length(cat_to_plot_options)
    cat_to_plot = cat_to_plot_options{ii_cat_to_plot};
    dir_clrs = {'b','r'};
    fig = figure;
    fig.WindowState = 'maximized';
    ax_posx = linspace(.05,1,length(cat_to_plot)+1);
    ax_posy = linspace(.06,.95,5+1);
    w = mean(diff(ax_posx))*.9;
    h = mean(diff(ax_posy))*.9;
    ax_posx(end)=[]; ax_posy(end)=[];
    clear panels
    for c = 1:length(ax_posx)
        for r = 1:length(ax_posy)
            panels(r,c) = axes('units','normalized','position', [ax_posx(c) ax_posy(r) w h]);
            hold on
        end
    end
    
    for col = 1:length(cat_to_plot)
        ii_cat = cat_to_plot(col);
        axes(panels(end,col))
        title(events_cat_names{ii_cat},'Interpreter','none')
        for ii_dir = 1:nDir
            axes(panels(1,col))
            plot(cell.FR_map(ii_dir).all.bin_centers, cell.FR_map(ii_dir).all.PSTH,'-','Color',dir_clrs{ii_dir});
            xlabel('Position (m)')
            if col == 1
                ylabel('Flight FR map (Hz)')
            end
    
            replay_PSTH = replay_PSTH_all(ii_cat,ii_dir);
            events = replay_PSTH_all(ii_cat,ii_dir).events;
    
            axes(panels(2,col))
            plot(replay_PSTH.bin_centers, replay_PSTH.PSTH,'-','Color',dir_clrs{ii_dir});
            if col == 1
                ylabel('Replay FR map (Hz)')
            end
    
            axes(panels(3,col))
            plot(replay_PSTH.bin_centers, replay_PSTH.time_spent,'-','Color',dir_clrs{ii_dir});
            yline(min_time_spent_per_meter,'--')
            if col == 1
                ylabel('Replay time Spent (s)')
            end
    
            axes(panels(3+ii_dir,col))
            if ~isempty(events)
                seqs = [events.seq_model];
                seqs_xi = [seqs.start_pos; seqs.end_pos];
                plot(seqs_xi, repmat(1:size(seqs_xi,2),[2 1]), 'Color',.5*[1 1 1]);
                plot(seqs_xi(1,:), 1:size(seqs_xi,2), '.', 'Color',.5*[1 1 1]);
                plot(replay_PSTH.spikes_pos, replay_PSTH.spikes_replay_num, '.', 'Color',dir_clrs{ii_dir})
                if col == 1
                    ylabel("#replay, dir "+ii_dir)
                end
            end
            
        end
    end
    linkaxes(panels(1,:),'y')
    linkaxes(panels(2,:),'y')
    linkaxes(panels(3,:),'y')
    linkaxes(panels,'x')
    hSG = sgtitle(sprintf('%s (%d)',cell.details.cell_ID,cell.details.cell_num),'interpreter','none');
    pause(0.001)
    hNC=hSG.NodeChildren.Children(2);
    hNC.Position(2)=0.93;
    pause(0.001)

    filename = sprintf('%s(%d)_%s_%s',cell.details.cell_ID,cell.details.cell_num,cell_analysis_name,cat_to_plot_options_names{ii_cat_to_plot});
    file_OUT = fullfile(dir_out,filename);
    saveas(fig,file_OUT,'jpg')
    
end
end

%% save data to file
replay_FR_map = struct();
replay_FR_map.replay_PSTH_all = replay_PSTH_all;
replay_FR_map.events_cat_names = events_cat_names;
filename = fullfile(dir_out, sprintf('%s_cell_%s',cell.details.cell_ID,cell_analysis_name));
save(filename, 'replay_FR_map');




%%
