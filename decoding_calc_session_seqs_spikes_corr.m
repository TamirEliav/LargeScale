function decoding_calc_session_seqs_spikes_corr(exp_ID)
% arguments
%     %% 
%     exp_ID = 'b9861_d180526'
%     epoch_type {mustBeMember(epoch_type,{'sleep','rest','flight'})} = 'sleep'
%     params_opt = 11;
%     event_type {mustBeMember(event_type,{'PE','posterior','ripples','MUA'})} = 'posterior'
% end

%%
% exp_ID = 'b2382_d190623';
% exp_ID = 'b2382_d190627';
% exp_ID = 'b0184_d191129';
% exp_ID = 'b0184_d191130';
% exp_ID = 'b0184_d191201';
% exp_ID = 'b0184_d191202';
% exp_ID = 'b0184_d191203';
% exp_ID = 'b0184_d191204';
% exp_ID = 'b0184_d191205'; % dir 1 replay 28!!!!!!!!!!!!

%%
lw = 3;
sym = '|';
ii_cat = 3;

%% load decoding results
params_opt = 11;

%% load replay data
% [events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
% [~, TF] = decoding_apply_seq_inclusion_criteria([events.seq_model]);
% events(~TF) = [];
% if isempty(events)
%     return;
% end
% seqs = [events.seq_model];

%% load cells data
cells_t = DS_get_cells_summary();
cells_exp_ID = cellfun(@(c)DS_get_exp_ID_from_cell_ID(c),cells_t.cell_ID,'UniformOutput',false);
TF = string(cells_exp_ID) == string(exp_ID);
cells_t(~TF,:)=[];
% whos cells_t
cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',1);
details = [cells.details];
details(~contains({details.brain_area}, {'CA1','CA3'})) = [];
details(~ismember([details.ClusterQuality], [2])) = [];
cells_t = cells_t({details.cell_ID},:);
cells = cellfun(@(c)cell_load_data(c,'details','signif','inclusion'),cells_t.cell_ID);
% whos cells_t cells
clear cells_t
if isempty(cells)
    return
end
inclusion = cat(1,cells.inclusion);
cells(~[inclusion(:,1).pyr])=[];
details = [cells.details];
if isempty(details)
    return
end
cells = cellfun(@(c)cell_load_data(c,'details','signif','inclusion','fields','FR_map','replay_FR_map','spikes','FE'),{details.cell_ID});
inclusion = cat(1,cells.inclusion);
signif = cat(1,cells.signif);

%% plot entire session
fig1=figure;
fig1.WindowState = 'maximized';
clear panels
panels(1,1) = axes('Units','normalized','Position',[0.05 0.03 0.45 0.05]);
panels(1,2) = axes('Units','normalized','Position',[0.05 0.10 0.45 0.05]);
panels(1,3) = axes('Units','normalized','Position',[0.05 0.16 0.45 0.3]);
panels(1,4) = axes('Units','normalized','Position',[0.05 0.47 0.45 0.50]);
panels(2,1) = axes('Units','normalized','Position',[0.53 0.03 0.45 0.05]);
panels(2,2) = axes('Units','normalized','Position',[0.53 0.10 0.45 0.05]);
panels(2,3) = axes('Units','normalized','Position',[0.53 0.16 0.45 0.3]);
panels(2,4) = axes('Units','normalized','Position',[0.53 0.47 0.45 0.50]);

for ii_dir = 1:2
    inclusion_dir = inclusion(:,ii_dir);
    signif_dir = signif(:,ii_dir);
    TF = [inclusion_dir.TF];
    TF = TF & [signif_dir.SI_thr_shuffle];
    TF = TF & [signif_dir.SI_thr_signif];
    cells_dir = cells(TF);
    nCellsDir = length(cells_dir);
    if nCellsDir==0
        continue
    end
    clrs = distinguishable_colors(nCellsDir);
    
    axes(panels(ii_dir,1));
    cla reset
    hold on
%     hax=gca; hax.ColorOrder = clrs;
    for ii_cell = 1:nCellsDir
        cell = cells_dir(ii_cell);
        x = cell.FR_map(ii_dir).all.bin_centers;
        y = cell.FR_map(ii_dir).all.PSTH;
        c = clrs(ii_cell,:);
        plot(x,y,'LineWidth',lw,'Color',c)
    end
    ylabel({'Flight';'FR map';'(Hz)'})
    
    axes(panels(ii_dir,2));
    cla reset
    hold on
%     hax=gca; hax.ColorOrder = clrs;
    for ii_cell = 1:nCellsDir
        cell = cells_dir(ii_cell);
        x = cell.replay_FR_map.replay_PSTH_all(ii_cat,ii_dir).bin_centers;
        y = cell.replay_FR_map.replay_PSTH_all(ii_cat,ii_dir).PSTH;
        c = clrs(ii_cell,:);
        plot(x,y,'LineWidth',lw,'Color',c)
    end
    ylabel({'Replay';'FR map';'(Hz)'})
    
    axes(panels(ii_dir,3));
    cla reset
    hold on
    FE = cells_dir(1).FE{ii_dir};
    FE_start_pos = arrayfun(@(fe)fe.pos([1]),FE,'UniformOutput',true);
    FE_end_pos = arrayfun(@(fe)fe.pos([end]),FE,'UniformOutput',true);
    x = [FE_start_pos; FE_end_pos];
    y = repmat(1:length(FE),2,1);
    plot(x,y,'Color',0.5*[1 1 1]);
%     hax=gca; hax.ColorOrder = clrs;
    for ii_cell = 1:nCellsDir
        cell = cells_dir(ii_cell);
        FE = cell.FE{ii_dir};
        x = [FE.spikes_pos];
        y = arrayfun(@(ii,n)ii.*ones(1,n) , 1:length(FE), [FE.num_spikes], 'UniformOutput', false);
        y = [y{:}];
        c = clrs(ii_cell,:);
        plot(x,y,sym,'LineWidth',lw,'Color',c)
    end
    ylabel('Flight no.')

    axes(panels(ii_dir,4));
    cla reset
    hold on
    events = cells_dir(1).replay_FR_map.replay_PSTH_all(ii_cat,ii_dir).events;
    if ~isempty(events)
        seqs = [events.seq_model];
        x = [[seqs.start_pos]; [seqs.end_pos]];
        y = repmat(1:length(seqs),2,1);
        plot(x,y,'Color',0.5*[1 1 1]);
    %     hax=gca; hax.ColorOrder = clrs;
        for ii_cell = 1:nCellsDir
            cell = cells_dir(ii_cell);
            replay = cell.replay_FR_map.replay_PSTH_all(ii_cat,ii_dir);
            x = replay.spikes_pos;
            y = replay.spikes_replay_num;
            c = clrs(ii_cell,:);
            plot(x,y,sym,'LineWidth',lw,'Color',c)
        end
        ylim([0 length(seqs)+1])
    end
    ylabel('Replay no.')
    
%     details_dir = [cells_dir.details];
%     cells_dir = cellfun(@(c)cell_load_data(c,'details','signif','inclusion','fields','FR_map','replay_FR_map'),{details_dir.cell_ID});
%     for ii_cell = 1:length(cells_dir)
%         [cells_dir(ii_cell).fields{ii_dir}.ii_cell] = deal(ii_cell);
%     end
    
%     fields_dir = arrayfun(@(c)c.fields{ii_dir},cells_dir,'UniformOutput',false);
%     fields_dir = [fields_dir{:}];
%     nFieldsDir = length([fields_dir{:}]);
%     replay_dir = arrayfun(@(cell)cell.replay_FR_map.replay_PSTH_all(ii_cat,ii_dir),cells_dir);
    events = cells_dir(1).replay_FR_map.replay_PSTH_all(ii_cat,ii_dir).events;
%     events = replay_dir(1).events;
    if isempty(events)
        continue
    end
    seqs = [events.seq_model];
%     for ii_seq = 1:length(seqs)
%         %%
%         seq = seqs(ii_seq);
%         X = [];
%         Y = [];
%         figure
%         hold on
%         plot([seq.start_ts seq.end_ts], [seq.start_pos seq.end_pos], '-r') ;
%         for ii_cell = 1:nCellsDir
%             cell = cells_dir(ii_cell);
%             cell.replay_FR_map.replay_PSTH_all(ii_cat,ii_dir).spike
%         end
%     end
end
linkaxes(panels,'x')
zoom
htb = annotation('textbox',[0.48 0.99 0.1 0.01]);
htb.String = exp_ID;
htb.Interpreter = 'none';
htb.FitBoxToText = 'on';

dir_out = 'F:\sequences\session_replay_spikes';
mkdir(dir_out)
filename = sprintf('%s_session_replays_spikes',exp_ID);
file_out = fullfile(dir_out,filename);
saveas(fig1,file_out,'fig')
saveas(fig1,file_out,'jpg')



%% plot individual replays
fig2 = figure;
fig2.Units = 'centimeters';
fig2.Position = [15 1 14 25];
clear panels
w = 0.8;
% h = 0.185;
h = 0.235;
x = 0.1;
% y = [0.05 0.28 0.5 0.75];
y = cumsum(h.*ones(1,4))-h+0.04;
y = flip(y);
y = [y(1) y];
for ii_y = 1:length(y)
    panels(ii_y) = axes('Units','normalized','Position',[x y(ii_y) w h]);
end

for ii_dir = 1:2
    inclusion_dir = inclusion(:,ii_dir);
    signif_dir = signif(:,ii_dir);
    TF = [inclusion_dir.TF];
    TF = TF & [signif_dir.SI_thr_shuffle];
    TF = TF & [signif_dir.SI_thr_signif];
    cells_dir = cells(TF);
    nCellsDir = length(cells_dir);
    clrs = distinguishable_colors(nCellsDir);
    if nCellsDir == 0
        continue;
    end
    events = cells_dir(1).replay_FR_map.replay_PSTH_all(ii_cat,ii_dir).events;
    if isempty(events)
        continue
    end
    seqs = [events.seq_model];
    
    for ii_seq = 1:length(seqs)
%         ii_seq
        seq = seqs(ii_seq);
        event = events(ii_seq);

        seq_ts = [seq.start_ts seq.end_ts];
        seq_pos = [seq.start_pos seq.end_pos];
        seq_pos_limits = sort(seq_pos,'ascend');

        XX = {};
        YY = {};
        YY2 = {};
        ZZ = {};
        TFTF = {};
        nFieldsInvolved = 0;
        allFieldLocs = [];
        fieldsCellnum = [];
        for ii_cell = 1:nCellsDir
            cell = cells_dir(ii_cell);
            fields = cell.fields{ii_dir};
            if isempty(fields)
                continue;
            end
            fields_edges = cat(1,fields.edges_prc);
            fields_TF = false(size(fields));
            for ii_field = 1:length(fields)
                field_int = fixed.Interval( fields_edges(ii_field,1),fields_edges(ii_field,1) );
                replay_int = fixed.Interval(seq_pos_limits(1),seq_pos_limits(2));
                fields_TF(ii_field) = overlaps(field_int,replay_int);
            end
            fields = fields(fields_TF);
            if isempty(fields)
                continue;
            end
            [~,~,x,~] = get_data_in_ti(cell.spikes.ts, seq_ts);
            if isempty(x)
                continue;
            end
            y = [fields.loc];
            y2 = [1:length(fields)]+nFieldsInvolved;
            nFieldsInvolved = nFieldsInvolved + length(fields);
            z = [fields.width_prc];
            allFieldLocs = [allFieldLocs y];
            fieldsCellnum = [fieldsCellnum ii_cell*ones(1,length(fields))];
            expected_spikes_loc = interp1(seq_ts,seq_pos, x, 'linear');
            
            [X,Y] = meshgrid(x, y);
            [X,Y2] = meshgrid(x, y2);
            [X,Z] = meshgrid(x, z);
            TF = false(size(Y));
            [~,nearest_field_per_spike] = min(abs(Y-expected_spikes_loc),[],1);
            IND = sub2ind(size(TF),nearest_field_per_spike,1:length(x));
            TF(IND) = true;
            
            X = X(:);
            Y = Y(:);
            Y2 = Y2(:);
            Z = Z(:);
            TF = TF(:);
            XX{ii_cell} = X;
            YY{ii_cell} = Y;
            YY2{ii_cell} = Y2;
            ZZ{ii_cell} = Z;
            TFTF{ii_cell} = TF;
%             plot(X, Y,sym, 'Color',clrs(ii_cell,:),'LineWidth',lw)
        end
        if nFieldsInvolved<4
            continue;
        end
        XXX = cat(1,XX{:});
        YYY = cat(1,YY{:});
        YYY2 = cat(1,YY2{:});
        ZZZ = cat(1,ZZ{:});
        TFTFTF = logical(cat(1,TFTF{:}));
        mark_height = interp1([min(ZZZ) max(ZZZ)],[20 50],ZZZ,'linear');

        if isnan(event.rest_ball_num)
            if ~exist('decode_sleep','var')
                decode_sleep = decoding_load_data(exp_ID, 'sleep', params_opt);
            end
            decode = decode_sleep;
        else
            if ~exist('decode_rest','var')
                decode_rest = decoding_load_data(exp_ID, 'rest', params_opt);
            end
            decode = decode_rest;
        end
        seq_ts_with_margin = seq_ts + [-1 1].*1e6*0.1;
        decode_IX = find(decode.time > seq_ts_with_margin(1) & decode.time < seq_ts_with_margin(2));
        decode_posterior = squeeze(decode.posterior(:,event.state_num,decode_IX));
        decode_time = decode.time(decode_IX);
        decode_pos = decode.pos;

        linkaxes(panels,'off');
        
        axes(panels(1))
        cla reset
        hold on
        imagesc(decode_time, decode_pos, decode_posterior);
        axis tight
        hax=gca;
        hax.Colormap = flip(gray);
        rescale_plot_data('x',[1e-6 seq_ts(1)])
        ylimits = seq_pos_limits + [-1 1].*0.1*range(seq_pos_limits);
        ylim(ylimits);
        ylabel('Position (m)')
        title(sprintf('%s dir:%d seq:%d',exp_ID,ii_dir,ii_seq),'Interpreter','none');

        axes(panels(2))
        cla reset
        hold on
        plot(seq_ts, seq_pos, '-r');
        hax=gca;
        hax.Colormap = clrs;
        hsc = scatter(XXX(~TFTFTF),YYY(~TFTFTF),[],fieldsCellnum(YYY2(~TFTFTF))','|','LineWidth',0.5,'SizeData',mark_height(~TFTFTF));
        hsc = scatter(XXX(TFTFTF),YYY(TFTFTF),[],fieldsCellnum(YYY2(TFTFTF))','|','LineWidth',2.5,'SizeData',mark_height(TFTFTF));
        ylim(ylimits);
        rescale_plot_data('x',[1e-6 seq_ts(1)])
        [rho_pearson, pval_pearson] = corr(XXX,YYY,"type","pearson");
        [rho_spearman, pval_spearman] = corr(XXX,YYY,"type","Spearman");
        text(0.95,0.95,sprintf('r=%.2g',rho_pearson),'Units','normalized');
        text(0.95,0.9,sprintf('P=%.2g',pval_pearson),'Units','normalized');
        text(0.95,0.80,sprintf('\\rho=%.2g',rho_spearman),'Units','normalized');
        text(0.95,0.75,sprintf('P=%.2g',pval_spearman),'Units','normalized');
        text(0.95,0.60,sprintf('nFields=%d',nFieldsInvolved),'Units','normalized');
        ylabel('Position (m)')
        panels(2).Color = 'None';

        axes(panels(3))
        cla reset
        hold on
        [allFieldLocs_sorted, IX_fields_locs_sorted] = sort(allFieldLocs,'ascend');
        YYY3 = interp1(IX_fields_locs_sorted, 1:nFieldsInvolved, YYY2, 'linear');
        hax=gca;
        hax.Colormap = clrs;
        hsc = scatter(XXX(~TFTFTF),YYY3(~TFTFTF),[],fieldsCellnum(YYY2(~TFTFTF))','|','LineWidth',0.5,'SizeData',mark_height(~TFTFTF));
        hsc = scatter(XXX(TFTFTF),YYY3(TFTFTF),[],fieldsCellnum(YYY2(TFTFTF))','|','LineWidth',2.5,'SizeData',mark_height(TFTFTF));
        ylim([1 nFieldsInvolved] + [-1 1].*0.1*(nFieldsInvolved-1))
        yticks(1:nFieldsInvolved)
        rescale_plot_data('x',[1e-6 seq_ts(1)])
        [rho_pearson, pval_pearson] = corr(XXX,YYY3,"type","pearson");
        [rho_spearman, pval_spearman] = corr(XXX,YYY3,"type","Spearman");
        text(0.95,0.95,sprintf('r=%.2g',rho_pearson),'Units','normalized');
        text(0.95,0.9,sprintf('P=%.2g',pval_pearson),'Units','normalized');
        text(0.95,0.80,sprintf('\\rho=%.2g',rho_spearman),'Units','normalized');
        text(0.95,0.75,sprintf('P=%.2g',pval_spearman),'Units','normalized');
        ylabel('field no.')

        % firing pattern matrix
        pos_fs = 5000;
        pos_dt = 1/pos_fs;
        pos_t = seq_ts(1) : pos_dt*1e6 : seq_ts(end);
        pos = interp1(seq_ts, seq_pos, pos_t, 'linear');
        bin_size = 0.1;
        ker_SD = 0.5;
        bin_edges = seq_pos_limits(1) : bin_size : seq_pos_limits(end);
        bin_centers = edges2centers(bin_edges);
        bin_centers_time = interp1(seq_pos, seq_ts, bin_centers, 'linear');
        bin_size_time_ms = median(diff(bin_centers_time))*1e-3;
        min_time_spent = 0;
        FR_replay = zeros(nCellsDir, length(bin_centers));
        FR_flight = zeros(nCellsDir, length(bin_centers));
        for ii_cell = 1:nCellsDir
            cell = cells_dir(ii_cell);
            [~,~,spikes_ts,~] = get_data_in_ti(cell.spikes.ts, seq_ts);
            spikes_pos = interp1(seq_ts,seq_pos, spikes_ts,'linear');
            FR_replay(ii_cell,:) = computePSTH(pos,pos_fs,spikes_pos,bin_edges,min_time_spent,ker_SD);
            FR_flight(ii_cell,:) = interp1(cell.FR_map(ii_dir).all.bin_centers, cell.FR_map(ii_dir).all.PSTH, bin_centers, 'linear');
        end
        FR_replay = FR_replay ./ max(FR_replay,[],2);
        FR_flight = FR_flight ./ max(FR_flight,[],2);
        [rho_spearman, pval_spearman] = corr(FR_replay(:),FR_flight(:),'type','Spearman','rows','pairwise');
        [rho_pearson, pval_pearson] = corr(FR_replay(:),FR_flight(:),'type','Pearson','rows','pairwise');
        % cell-shuffle
        nShuffle = 1000;
        rho_pearson_shuffles = zeros(1,nShuffle);
        rho_spearman_shuffles = zeros(1,nShuffle);
        for ii_shuffle = 1:nShuffle
            shuffle_IX = randperm(nCellsDir,nCellsDir);
            FR_flight_shuffle = FR_flight(shuffle_IX,:);
            rho_spearman_shuffles(ii_shuffle) = corr(FR_replay(:),FR_flight_shuffle(:),'type','Spearman','rows','pairwise');
            rho_pearson_shuffles(ii_shuffle) = corr(FR_replay(:),FR_flight_shuffle(:),'type','Pearson','rows','pairwise');
        end
        rho_pearson_shuffles_z = (rho_pearson - mean(rho_pearson_shuffles)) / std(rho_pearson_shuffles);
        rho_spearman_shuffles_z = (rho_spearman - mean(rho_spearman_shuffles)) / std(rho_spearman_shuffles);
       
        axes(panels(4))
        cla reset
        hold on
        imagesc(bin_centers_time, 1:nCellsDir, FR_replay);
        scatter(bin_centers_time(1)-0.02*range(bin_centers_time), 1:nCellsDir,50,clrs,'|','LineWidth',10)
%         xlim(seq_ts + [-1 1].*0.05*range(seq_ts))
        rescale_plot_data('x',[1e-6 seq_ts(1)])
        yticks(1:nCellsDir)
        text(0.97,0.85,sprintf('r=%.g',rho_pearson),'Units','normalized');
        text(0.97,0.8,sprintf('P=%.g',pval_pearson),'Units','normalized');
        text(0.97,0.70,sprintf('\\rho=%.g',rho_spearman),'Units','normalized');
        text(0.97,0.65,sprintf('P=%.g',pval_spearman),'Units','normalized');
        text(0.97,0.55,'shuffle','Units','normalized');
        text(0.97,0.45,sprintf('r=%.2g (z)',rho_pearson_shuffles_z),'Units','normalized');
        text(0.97,0.35,sprintf('\\rho=%.2g (z)',rho_spearman_shuffles_z),'Units','normalized');
        ylabel('cell no.')
        title('Replay firing rate','Units','normalized','Position',[0.5 0.92])

        axes(panels(5))
        cla reset
        hold on
        imagesc(bin_centers_time, 1:nCellsDir, FR_flight)
        scatter(bin_centers_time(1)-0.02*range(bin_centers_time), 1:nCellsDir,50,clrs,'|','LineWidth',10)
        rescale_plot_data('x',[1e-6 seq_ts(1)])
        xlabel('Time (s)')
        ylabel('cell no.')
        yticks(1:nCellsDir)
        title('Flight firing rate','Units','normalized','Position',[0.5 0.92])

        linkaxes(panels,'x');
        xlim([0 seq.duration]+ [-1 1].*0.05*seq.duration);

        dir_out = 'F:\sequences\session_replay_spikes\individual_replays';
        mkdir(dir_out)
        filename = sprintf('%s_dir%d_replay_%d_spikes',exp_ID,ii_dir,ii_seq);
        file_out = fullfile(dir_out,filename);
        saveas(fig2,file_out,'jpg')

    end
end


%%