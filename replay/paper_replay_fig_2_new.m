%% Replay - Fig 2( new) - replay at single-units level
%%
clear 
clc
close all

%% plotting options
ii_cat = 3; % pooled
params_opt = 11;
seq_examples_list = [ ...FR_peaks_loc
    struct('exp_ID','b0184_d191205', 'ii_dir',1, 'ii_seq', 18, 'ii_FE', 9);...% FE=2/5/9/18/26/27/31
    struct('exp_ID','b0184_d191205', 'ii_dir',1, 'ii_seq', 28, 'ii_FE', 30);...% FE=20/30/31
    ];
norm_FR_matrix_per_cell = 1;

%% graphics params

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Fig_2_new';
fig_caption_str = ' ';
log_name_str = [fig_name_str '_log_file' '.txt'];
log_name_str = strrep(log_name_str , ':', '-');
log_name_str = strrep(log_name_str , ' ', '_');
log_name_out = fullfile(res_dir, log_name_str);

%% open log file
diary off
diary(log_name_out)
diary on
disp('Log file');
disp(['created: ', datestr(clock)]);
disp('======================================================');
disp([fig_name_str ':' fig_caption_str]);   
disp('======================================================');
disp('');

%% create figure
% figure_size_cm = [21.0 29.7]; % ~A4
figure_size_cm = [21.6 27.9]; % ~US letter
close all
fig = figure;
% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',7);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'DefaultAxesUnits','centimeters');
set(gcf,'PaperType','usletter')
% set(gcf,'PaperType','<custom>');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 figure_size_cm]);
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]); % position on screen...
set(gcf, 'Renderer', 'painters');
% set(gcf, 'color', 'none');
set(groot, 'defaultAxesColor','None')
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');

% create panels
clear panels

w = 5;
h = 3;
x = 3;
y = 10+[0 4 8 12 12];
y = flip(y);
panels{1}(1,1) = axes('position', [x y(1) w h]);
panels{1}(1,2) = axes('position', [x y(2) w h]);
panels{1}(1,3) = axes('position', [x y(3) w h]);
panels{1}(1,4) = axes('position', [x y(4) w h]);
panels{1}(1,5) = axes('position', [x y(5) w h]);
x = x + 8;
panels{1}(2,1) = axes('position', [x y(1) w h]);
panels{1}(2,2) = axes('position', [x y(2) w h]);
panels{1}(2,3) = axes('position', [x y(3) w h]);
panels{1}(2,4) = axes('position', [x y(4) w h]);
panels{1}(2,5) = axes('position', [x y(5) w h]);

x = 3;
y = 5;
panels{2}(1) = axes('position', [x y w h]);
panels{2}(2) = axes('position', [x+0.2 y+2 .5 .5]);

total_offset = [0 0];
for ii = 1:length(panels)
    subpanels = panels{ii};
    subpanels = subpanels(:);
    for jj = 1:length(subpanels)
        subpanels(jj).Position([1 2]) = subpanels(jj).Position([1 2]) + total_offset;
    end
end

%%
for ii_ex = 1:length(seq_examples_list)
    
    %% load cells data
    ex = seq_examples_list(ii_ex);
    exp_ID = ex.exp_ID;
    ii_dir = ex.ii_dir;
    ii_seq = ex.ii_seq;
    ii_FE = ex.ii_FE;
    cells_t = DS_get_cells_summary();
    cells_exp_ID = cellfun(@(c)DS_get_exp_ID_from_cell_ID(c),cells_t.cell_ID,'UniformOutput',false);
    TF = string(cells_exp_ID) == string(exp_ID);
    cells_t(~TF,:)=[];
    cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',1);
    details = [cells.details];
    details(~contains({details.brain_area}, {'CA1','CA3'})) = [];
    details(~ismember([details.ClusterQuality], [2])) = [];
    cells_t = cells_t({details.cell_ID},:);
    cells = cellfun(@(c)cell_load_data(c,'details','signif','inclusion'),cells_t.cell_ID);
    clear cells_t
    inclusion = cat(1,cells.inclusion);
    cells(~[inclusion(:,1).pyr])=[];
    details = [cells.details];
    cells = cellfun(@(c)cell_load_data(c,'details','signif','inclusion','fields','FR_map','replay_FR_map','spikes','FE'),{details.cell_ID});
    inclusion = cat(1,cells.inclusion);
    signif = cat(1,cells.signif);

    %% get relevant cells for map direction
    inclusion_dir = inclusion(:,ii_dir);
    signif_dir = signif(:,ii_dir);
    TF = [inclusion_dir.TF];
    TF = TF & [signif_dir.SI_thr_shuffle];
    TF = TF & [signif_dir.SI_thr_signif];
    cells_dir = cells(TF);
    nCellsDir = length(cells_dir);
    clrs = distinguishable_colors(nCellsDir);
    events = cells_dir(1).replay_FR_map.replay_PSTH_all(ii_cat,ii_dir).events;
    seqs = [events.seq_model];

    %% get the relevant seq
    seq = seqs(ii_seq);
    event = events(ii_seq);
    seq_ts = [seq.start_ts seq.end_ts];
    seq_pos = [seq.start_pos seq.end_pos];
    seq_pos_limits = sort(seq_pos,'ascend');

    %% arrange data for plotting
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
    end
    XXX = cat(1,XX{:});
    YYY = cat(1,YY{:});
    YYY2 = cat(1,YY2{:});
    ZZZ = cat(1,ZZ{:});
    TFTFTF = logical(cat(1,TFTF{:}));
    mark_height = interp1([min(ZZZ) max(ZZZ)],[20 50],ZZZ,'linear');

    if isnan(event.rest_ball_num)
        epoch_type = 'sleep';
    else
        epoch_type = 'rest';
    end
    fprintf('epoch type = %s (#%d)\n',epoch_type,event.epoch_num)
    decode = decoding_load_data(exp_ID, epoch_type, params_opt);
    seq_ts_with_margin = seq_ts + [-1 1].*1e6*0.1;
    decode_IX = find(decode.time > seq_ts_with_margin(1) & decode.time < seq_ts_with_margin(2));
    decode_posterior = squeeze(decode.posterior(:,event.state_num,decode_IX));
    decode_time = decode.time(decode_IX);
    decode_pos = decode.pos;

    %% plot seq spikes
    % this panel is made of two panels: 1) posterior imagesc and 2) spikes
    axes(panels{1}(ii_ex,1))
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

    axes(panels{1}(ii_ex,2))
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
    xoffset = 1.05;
    text(xoffset,0.95,sprintf('r=%.2g',rho_pearson),'Units','normalized');
    text(xoffset,0.83,sprintf('P=%.2g',pval_pearson),'Units','normalized');
    text(xoffset,0.60,sprintf('\\rho=%.2g',rho_spearman),'Units','normalized');
    text(xoffset,0.48,sprintf('P=%.2g',pval_spearman),'Units','normalized');
    text(xoffset,0.30,sprintf('nFields=%d',nFieldsInvolved),'Units','normalized');
    ylabel('Position (m)')
    hax.Color = 'None';

    %% plot flight spikes (single flight example)
    axes(panels{1}(ii_ex,3))
    cla reset
    hold on
    for ii_cell = 1:nCellsDir
        cell = cells_dir(ii_cell);
        FE = cell.FE{ii_dir}(ii_FE);
        if ii_cell==1
            lm = fitlm([FE.pos],[FE.ts]);
            xlimits = feval(lm,seq_pos);
%             xlimits = interp1([FE.pos],[FE.ts],seq_pos,'linear','extrap');
            time_convert_src = xlimits;
            time_convert_dst = seq_ts;
            xlimits = seq_ts;
            x = interp1(time_convert_src,time_convert_dst,[FE.ts], 'linear','extrap');
            y = [FE.pos];
            plot(x,y,'k-');
        end
        x = interp1(time_convert_src,time_convert_dst, [FE.spikes_ts], 'linear','extrap');
        y = [FE.spikes_pos];
        plot(x,y, '|', 'Color',clrs(ii_cell,:))
    end
    xlim(xlimits)
    ylim(ylimits)
    rescale_plot_data('x',[1e-6 xlimits(1)])
    text(0.05,0.85,"flight #"+ii_FE,'HorizontalAlignment','Left','units','normalized');

    %% calc firing pattern matrix
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
    if norm_FR_matrix_per_cell
        FR_replay = FR_replay ./ max(FR_replay,[],2);
        FR_flight = FR_flight ./ max(FR_flight,[],2);
    end
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

    %% sort cells for FR plots, according to FR peaks in flight
    [~,FR_peaks_loc] = max(FR_flight,[],2);
    [~, FR_sort_cells_IX] = sort(FR_peaks_loc);

    %% plot replay firing rate pattern
    axes(panels{1}(ii_ex,4))
    cla reset
    hold on
    imagesc(bin_centers_time, 1:nCellsDir, FR_replay(FR_sort_cells_IX,:));
    scatter(bin_centers_time(1)-0.02*range(bin_centers_time), 1:nCellsDir,50,clrs(FR_sort_cells_IX,:),'|','LineWidth',10,'clipping','off')
%         xlim(seq_ts + [-1 1].*0.05*range(seq_ts))
    rescale_plot_data('x',[1e-6 seq_ts(1)])
    yticks(1:nCellsDir)
    xoffset = 1.05;
    text(xoffset,0.90,sprintf('r=%.g',rho_pearson),'Units','normalized');
    text(xoffset,0.80,sprintf('P=%.g',pval_pearson),'Units','normalized');
    text(xoffset,0.65,sprintf('\\rho=%.g',rho_spearman),'Units','normalized');
    text(xoffset,0.55,sprintf('P=%.g',pval_spearman),'Units','normalized');
    text(xoffset,0.30,'shuffle','Units','normalized');
    text(xoffset,0.18,sprintf('r=%.2g (z)',rho_pearson_shuffles_z),'Units','normalized');
    text(xoffset,0.05,sprintf('\\rho=%.2g (z)',rho_spearman_shuffles_z),'Units','normalized');
    ylabel('cell no.')
    title('Replay firing rate','Units','normalized','Position',[0.5 0.92])

    %% plot flight firing rate pattern
    axes(panels{1}(ii_ex,5))
    cla reset
    hold on
    imagesc(bin_centers_time, 1:nCellsDir, FR_flight(FR_sort_cells_IX,:))
    scatter(bin_centers_time(1)-0.02*range(bin_centers_time), 1:nCellsDir,50,clrs(FR_sort_cells_IX,:),'|','LineWidth',10,'clipping','off')
    rescale_plot_data('x',[1e-6 seq_ts(1)])
    xlabel('Time (s)')
    ylabel('cell no.')
    yticks(1:nCellsDir)
    title('Flight firing rate','Units','normalized','Position',[0.5 0.92])

    %% link axes
    linkaxes(panels{1}(ii_ex,[1:5]),'x');
    xlim([0 seq.duration]+ [-1 1].*0.05*seq.duration);
end

%% single unit FR corr hist
axes(panels{2}(1))
cla reset
hold on
FR_corr = load("E:\Tamir\work\PROJECTS\LargeScale\paper_replay\figures\single_unit_replay_tuning_corr.mat");
histogram(cat(1,FR_corr.ccc_data{:}),'FaceColor',0.5*[1 1 1],'normalization','pdf','NumBins',21,'EdgeColor','None');
histogram(cat(1,FR_corr.ccc_shuffle{:}),'EdgeColor','k','DisplayStyle','Stairs','normalization','pdf','LineWidth',1.5);
xlim([-1 1])
xlabel({'Correlation of firing rate maps';'during flight vs. during replay'})
ylabel('Probability density')

%% add legend
axes(panels{2}(2))
cla reset
hold on
axis off
axis equal
patch([0 0 1 1],[0 1 1 0], 0.5*[1 1 1],'EdgeColor','None');
plot([0 1],[1 1].*2, 'k-','LineWidth',1.5)
text(1.5,0.5, 'shuffle','FontSize',7)
text(1.5,2, 'Data','FontSize',7)

%% add panel letters
font_size = 11;
axes(panels{1}(1,2))
text(-0.2,1.12, 'a', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{1}(1,3))
text(-0.2,1.12, 'b', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{1}(1,4))
text(-0.2,1.12, 'c', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{1}(1,5))
text(-0.2,1.12, 'd', 'Units','normalized','FontWeight','bold','FontSize',font_size);

axes(panels{2}(1))
text(-0.2,1.12, 'e', 'Units','normalized','FontWeight','bold','FontSize',font_size);

%%
fig_name = sprintf('%s_norm_FR_per_cell=%d_ex_1_FE_%d_ex_2_FE_%d', ...
    fig_name_str, ...
    norm_FR_matrix_per_cell,...
    seq_examples_list(1).ii_FE,...
    seq_examples_list(2).ii_FE)
file_out = fullfile(res_dir, fig_name);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
disp('figure saved!')

%%
