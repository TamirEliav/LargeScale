%% Replay - Fig 2( new) - replay at single-units level
%%
clear 
clc
close all

%% plotting options
ii_cat = 3; % pooled
params_opt = 11;
seq_examples_list = [ ...FR_peaks_loc
    struct('exp_ID','b0184_d191205', 'ii_dir',1, 'ii_seq', 18, 'ii_FE', 9);...%  FE=2/5/9/18/26/27/31,  best=9
    struct('exp_ID','b0184_d191205', 'ii_dir',1, 'ii_seq', 28, 'ii_FE', 31);...% FE=20/30/31,           best=30, nachums choice 31
    ];
norm_FR_matrix_per_cell = 1;
undetected_field_cell_ID = 'b0184_d191205_TT12_SS02';
undetected_field_map_ii_dir = 1;
undetected_field_loc = 36.2;

% tuning_corr_cell_examples_list = [...
%     struct('cell_num',2070, 'ii_dir', 2),... % 1540
%     struct('cell_num',2046, 'ii_dir', 1),... % 1017
%     struct('cell_num',1297, 'ii_dir', 1),... %  499
%     struct('cell_num', 796, 'ii_dir', 1),... %  211
%     struct('cell_num',2044, 'ii_dir', 1),... % 1490
%     struct('cell_num',2055, 'ii_dir', 1),... % 1020
% struct('cell_num', 999, 'ii_dir', 2),... %  1309
% struct('cell_num',2181, 'ii_dir', 2),... %  1560
% struct('cell_num',1387, 'ii_dir', 2),... %  1379
% struct('cell_num',2044, 'ii_dir', 1),... %  1016
% struct('cell_num', 525, 'ii_dir', 1),... %    72
% struct('cell_num',1607, 'ii_dir', 1),... %   746
%     struct('cell_num',2070, 'ii_dir', 1),... %  1022
%     struct('cell_num', 752, 'ii_dir', 1),... %   178
%     struct('cell_num',1450, 'ii_dir', 2),... %  1395
%     struct('cell_num',1406, 'ii_dir', 2),... %  1383
%     struct('cell_num',1984, 'ii_dir', 1),... %   853
%     struct('cell_num',2055, 'ii_dir', 1),... %  1020
% struct('cell_num',1066, 'ii_dir', 2),... %  1329
% struct('cell_num',2107, 'ii_dir', 1),... %  1059
% struct('cell_num',2046, 'ii_dir', 1),... %  1017
% struct('cell_num',1832, 'ii_dir', 1),... %   830
% struct('cell_num',2104, 'ii_dir', 2),... %  1552
% struct('cell_num',1922, 'ii_dir', 2),... %  1490
%     struct('cell_num',1074, 'ii_dir', 2),... %  1331
%     struct('cell_num',2148, 'ii_dir', 1),... %  1064
%     struct('cell_num',2070, 'ii_dir', 2),... %  1540
%     struct('cell_num', 796, 'ii_dir', 1),... %   211
%     struct('cell_num', 830, 'ii_dir', 1),... %   216
%     struct('cell_num',1680, 'ii_dir', 1),... %   766
% struct('cell_num',1738, 'ii_dir', 2),... %  1444
% struct('cell_num',2117, 'ii_dir', 2),... %  1554
% struct('cell_num',1297, 'ii_dir', 1),... %   499
% struct('cell_num',1738, 'ii_dir', 2),... %  1444
% struct('cell_num',2117, 'ii_dir', 2),... %  1554
% struct('cell_num',1297, 'ii_dir', 1),... %   499
% ];
% panel_d_opt = 6;
% tuning_corr_cell_examples_list = tuning_corr_cell_examples_list([1:6]+(panel_d_opt-1)*6);
tuning_corr_cell_examples_list = [...
    struct('cell_num',2046, 'ii_dir', 1),...
    struct('cell_num',2044, 'ii_dir', 1),...
    struct('cell_num',1984, 'ii_dir', 1),...
    struct('cell_num', 796, 'ii_dir', 1),...
    struct('cell_num',2070, 'ii_dir', 2),...
    struct('cell_num',1738, 'ii_dir', 2),...
    ];

%% graphics params
spikes_tick_height = 10;

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Fig_2';
fig_caption_str = 'single unit activity during replays';
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

panels{1}(1) = axes('position', [2 23 13 2.5]);
x = [2 7 12];
y = 17.5;
panels{2}(1,1) = axes('position', [x(1) y 4 4]);
panels{2}(1,2) = axes('position', [x(2) y 4 4]);
panels{2}(1,3) = axes('position', [x(3) y 3 4]);
y = 12.5;
panels{2}(2,1) = axes('position', [x(1) y 4 4]);
panels{2}(2,2) = axes('position', [x(2) y 4 4]);
panels{2}(2,3) = axes('position', [x(3) y 3 4]);
x = [2 9];
y = [3 1.5 0]+7;
w = 6;
h = 1.0;
for ii_x = 1:length(x)
    for ii_y = 1:length(y)
        panels{3}(ii_x,ii_y) = axes('position', [x(ii_x) y(ii_y) w h]);
    end
end
panels{3} = panels{3}(:);
x = 2;
y = 3;
panels{4}(1) = axes('position', [x y 5 3]);
panels{4}(2) = axes('position', [x+0.3 y+2.2 .5 .5]);

% w = 5;
% h = 3;
% x = 3;
% y = 6+[0 4 8 12 12];
% y = flip(y);
% panels{1}(1,1) = axes('position', [x y(1) w h]);
% panels{1}(1,2) = axes('position', [x y(2) w h]);
% panels{1}(1,3) = axes('position', [x y(3) w h]);
% panels{1}(1,4) = axes('position', [x y(4) w h]);
% panels{1}(1,5) = axes('position', [x y(5) w h]);
% x = x + 8;
% panels{1}(2,1) = axes('position', [x y(1) w h]);
% panels{1}(2,2) = axes('position', [x y(2) w h]);
% panels{1}(2,3) = axes('position', [x y(3) w h]);
% panels{1}(2,4) = axes('position', [x y(4) w h]);
% panels{1}(2,5) = axes('position', [x y(5) w h]);
% 
% x = 3;
% y = 5;
% panels{2}(1) = axes('position', [x y w h]);
% panels{2}(2) = axes('position', [x+0.2 y+2 .5 .5]);

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
    events = cells_dir(1).replay_FR_map.replay_PSTH_all(ii_cat,ii_dir).events;
    seqs = [events.seq_model];

    %% manually add the undetected field 
    for ii_cell = 1:nCellsDir
        if strcmp(cells_dir(ii_cell).details.cell_ID , undetected_field_cell_ID)
            fields = cells_dir(ii_cell).fields{undetected_field_map_ii_dir};
            new_field = fields(1);
            new_field.loc = undetected_field_loc;
            fields(end+1) = new_field;
            [~,sorted_field_locs] = sort([fields.loc],'ascend');
            fields = fields(sorted_field_locs);
            cells_dir(ii_cell).fields{undetected_field_map_ii_dir} = fields;
        end
    end

    %% sort cells for FR plots, according to first field
    first_field_locs_all_cells = [];
    xlimits = [7 75];
    for ii_cell = 1:nCellsDir
        cell = cells_dir(ii_cell);
        cell_fields = cell.fields{ii_dir};
        field_locs = [cell_fields.loc];
        field_locs(field_locs<xlimits(1) | field_locs>xlimits(2)) =  [];
        first_field_locs_all_cells(ii_cell) = field_locs(1);
    end
    [~, sort_cells_IX] = sort(first_field_locs_all_cells,'ascend');
    cells_dir = cells_dir(sort_cells_IX);

    %% get the relevant seq
    seq = seqs(ii_seq);
    event = events(ii_seq);
    seq_ts = [seq.start_ts seq.end_ts];
    seq_pos = [seq.start_pos seq.end_pos];
    seq_pos_limits = sort(seq_pos,'ascend');

    %% choose colors for cells
    clrs = distinguishable_colors(nCellsDir);
%     rng(1);
%     clrs_order = randperm(nCellsDir);
    clrs_order = [1 3 2 4 5];
%     clrs_order = 1:5;
    clrs = clrs(clrs_order,:);
    clrs(3,:) = [1 0.5 0];

    %% plot cells FR map
    xlimits = [7 75];
    if ii_ex == 1
        axes(panels{1}(1))
        cla reset
        hold on
        for ii_cell = length(cells_dir):-1:1
            cell = cells_dir(ii_cell);
            FR = cell.FR_map(ii_dir).all;
            x = FR.bin_centers;
            y = FR.PSTH;
            IX = x>xlimits(1) & x<xlimits(2);
            y = y./max(y(IX))*1.2;
            y = y + ii_cell*0.2;
            plot(x,y, 'Color', clrs(ii_cell,:),'LineWidth',2);
            text(xlimits(1)-2, ii_cell*0.2+0.06, "Cell "+ii_cell, 'FontSize',6, 'HorizontalAlignment','right','Units','data','Color',clrs(ii_cell,:),'FontWeight','bold')
        end
        xlim(xlimits)
        xlabel('Position (m)','Units','normalized','Position',[0.5 -0.17])
        ylabel('firing rate (norm.)')
        hax=gca;
        hax.XRuler.TickLength(1) = 0.009;
        hax.YRuler.TickLength(1) = 0.007;
        hax.XRuler.TickLabelGapOffset = -2;
        hax.YRuler.TickLabelGapOffset = 1;
        hax.YAxis.Visible = 'off';
    end

    %% plot cells FR map (zoom in to seq_pos)
    axes(panels{2}(ii_ex,3))
    cla reset
    hold on
    for ii_cell = length(cells_dir):-1:1
        cell = cells_dir(ii_cell);
        FR = cell.FR_map(ii_dir).all;
        x = FR.bin_centers;
        y = FR.PSTH;
        IX = x>seq_pos(1) & x<seq_pos(2);
        y = y./max(y(IX))*1;
        y = y + ii_cell*0.2;
        plot(y,x, 'Color', clrs(ii_cell,:),'LineWidth',.65);
    end
    ylim(seq_pos)
    xlabel('firing rate (norm.)')
%     ylabel('Position (m)')
    hax=gca;
    hax.XRuler.TickLength(1) = 0.02;
    hax.YRuler.TickLength(1) = 0.01;
    hax.XRuler.TickLabelGapOffset = -2;
    hax.YRuler.TickLabelGapOffset = 1;
    hax.XAxis.Visible = 'off';

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
%     decode = decoding_load_data(exp_ID, epoch_type, params_opt);
%     seq_ts_with_margin = seq_ts + [-1 1].*1e6*0.1;
%     decode_IX = find(decode.time > seq_ts_with_margin(1) & decode.time < seq_ts_with_margin(2));
%     decode_posterior = squeeze(decode.posterior(:,event.state_num,decode_IX));
%     decode_time = decode.time(decode_IX);
%     decode_pos = decode.pos;

    %% plot seq spikes
%     % this panel is made of two panels: 1) posterior imagesc and 2) spikes
%     axes(panels{2}(ii_ex,1))
%     cla reset
%     hold on
%     imagesc(decode_time, decode_pos, decode_posterior);
%     axis tight
%     hax=gca;
%     hax.Colormap = flip(gray);
%     rescale_plot_data('x',[1e-6 seq_ts(1)])
%     ylimits = seq_pos_limits + [-1 1].*0.1*range(seq_pos_limits);
%     ylim(ylimits);
%     ylabel('Position (m)')
%     title(sprintf('%s dir:%d seq:%d',exp_ID,ii_dir,ii_seq),'Interpreter','none');

    axes(panels{2}(ii_ex,2))
    cla reset
    hold on
    plot(seq_ts, seq_pos, '-r');
    hax=gca;
    hax.Colormap = clrs;
%     hsc = scatter(XXX(~TFTFTF),YYY(~TFTFTF),[],fieldsCellnum(YYY2(~TFTFTF))','|','LineWidth',0.5,'SizeData',mark_height(~TFTFTF));
%     hsc = scatter(XXX(TFTFTF),YYY(TFTFTF),[],fieldsCellnum(YYY2(TFTFTF))','|','LineWidth',2.5,'SizeData',mark_height(TFTFTF));
    hsc = scatter(XXX(~TFTFTF),YYY(~TFTFTF),[],fieldsCellnum(YYY2(~TFTFTF))','|','LineWidth',0.5,'SizeData',spikes_tick_height);
    x = XXX(TFTFTF);
    y = YYY(TFTFTF);
    c = fieldsCellnum(YYY2(TFTFTF))';
    rng_opts = [0 27];
    rng(rng_opts(ii_ex)); % 11 12 16 27 47
    rand_order_IX = randperm(length(x));
    x = x(rand_order_IX);
    y = y(rand_order_IX);
    c = c(rand_order_IX);
    hsc = scatter(x,y,[],c,'|','LineWidth',2,'SizeData',spikes_tick_height);
    ylimits = seq_pos_limits + [-1 1].*0.1*range(seq_pos_limits);
    ylim(ylimits);
    xlim(seq_ts)
    rescale_plot_data('x',[1e-6 seq_ts(1)])
    [rho_pearson, pval_pearson] = corr(XXX,YYY,"type","pearson");
    [rho_spearman, pval_spearman] = corr(XXX,YYY,"type","Spearman");
%     xoffset = 1.05;
%     text(xoffset,0.95,sprintf('r=%.2g',rho_pearson),'Units','normalized');
%     text(xoffset,0.83,sprintf('P=%.2g',pval_pearson),'Units','normalized');
%     text(xoffset,0.60,sprintf('\\rho=%.2g',rho_spearman),'Units','normalized');
%     text(xoffset,0.48,sprintf('P=%.2g',pval_spearman),'Units','normalized');
%     text(xoffset,0.30,sprintf('nFields=%d',nFieldsInvolved),'Units','normalized');
    xlabel('Time (s)','units','normalized','position',[0.5 -0.09]);
    if ii_ex == 1
        title('Spiking during replay');
    end
%     ylabel('Position (m)')
%     hax.Color = 'None';
    hax.XRuler.TickLength(1) = 0.02;
    hax.YRuler.TickLength(1) = 0.01;
    hax.XRuler.TickLabelGapOffset = -2;
    hax.YRuler.TickLabelGapOffset = 1;

    %% plot flight spikes (single flight example)
    axes(panels{2}(ii_ex,1))
    cla reset
    hold on
    x_all = [];
    y_all = [];
    c_all = [];
    g_all = [];
    for ii_cell = 1:nCellsDir
        cell = cells_dir(ii_cell);
        FE = cell.FE{ii_dir}(ii_FE);
        if ii_cell==1
%             lm = fitlm([FE.pos],[FE.ts]);
%             xlimits = feval(lm,seq_pos);
% %             xlimits = interp1([FE.pos],[FE.ts],seq_pos,'linear','extrap');
%             time_convert_src = xlimits;
%             time_convert_dst = seq_ts;
%             xlimits = seq_ts;
%             x = interp1(time_convert_src,time_convert_dst,[FE.ts], 'linear','extrap');
            x = [FE.ts];
            y = [FE.pos];
            plot(x,y,'k-','LineWidth',0.5,'Color',0.5*[1 1 1]);
            lm = fitlm(y,x);
            xlimits = feval(lm, ylimits);
        end
%         x = interp1(time_convert_src,time_convert_dst, [FE.spikes_ts], 'linear','extrap');
        x = [FE.spikes_ts]';
        y = [FE.spikes_pos]';
        c = clrs(ii_cell,:);
        c = repmat(c,length(x),1);
        g = ii_cell.*ones(size(x));
        x_all = [x_all; x];
        y_all = [y_all; y];
        c_all = [c_all; c];
        g_all = [g_all; g];
    end
    rng(4);
    rand_order_IX = randperm(length(x_all));
    x_all = x_all(rand_order_IX);
    y_all = y_all(rand_order_IX);
    c_all = c_all(rand_order_IX,:);
    g_all = g_all(rand_order_IX,:);
    yoffset = zeros(size(g_all));
    yoffset(g_all==5) = +2;
    yoffset(g_all==4) = +1;
    yoffset(g_all==3) = 0;
    yoffset(g_all==2) = -1;
    yoffset(g_all==1) = -2;
    yoffset_gain = 0.01.*range(ylimits);
    yoffset_gain = yoffset_gain*2; % opt 2
    yoffset = yoffset.*yoffset_gain;
    scatter(x_all,y_all+yoffset,[], c_all, '|', 'LineWidth',.1, 'SizeData',spikes_tick_height);
    xlim(xlimits)
    ylim(ylimits)
    rescale_plot_data('x',[1e-6 xlimits(1)])
    text(0.05,0.85,"Flight #"+ii_FE,'HorizontalAlignment','Left','units','normalized','FontSize',7);
    xlabel('Time (s)','units','normalized','position',[0.5 -0.09]);
    ylabel('Position (m)')
    if ii_ex == 1
        title('Spiking during flight');
    end
    hax=gca;
    hax.XRuler.TickLength(1) = 0.02;
    hax.YRuler.TickLength(1) = 0.01;
    hax.XRuler.TickLabelGapOffset = -1.5;
    hax.YRuler.TickLabelGapOffset = 1;

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

    %% plot replay firing rate pattern
%     axes(panels{1}(ii_ex,4))
%     cla reset
%     hold on
%     imagesc(bin_centers_time, 1:nCellsDir, FR_replay(sort_cells_IX,:));
%     scatter(bin_centers_time(1)-0.02*range(bin_centers_time), 1:nCellsDir,50,clrs(sort_cells_IX,:),'|','LineWidth',10,'clipping','off')
% %         xlim(seq_ts + [-1 1].*0.05*range(seq_ts))
%     rescale_plot_data('x',[1e-6 seq_ts(1)])
%     yticks(1:nCellsDir)
%     xoffset = 1.05;
%     text(xoffset,0.90,sprintf('r=%.g',rho_pearson),'Units','normalized');
%     text(xoffset,0.80,sprintf('P=%.g',pval_pearson),'Units','normalized');
%     text(xoffset,0.65,sprintf('\\rho=%.g',rho_spearman),'Units','normalized');
%     text(xoffset,0.55,sprintf('P=%.g',pval_spearman),'Units','normalized');
%     text(xoffset,0.30,'shuffle','Units','normalized');
%     text(xoffset,0.18,sprintf('r=%.2g (z)',rho_pearson_shuffles_z),'Units','normalized');
%     text(xoffset,0.05,sprintf('\\rho=%.2g (z)',rho_spearman_shuffles_z),'Units','normalized');
%     ylabel('cell no.')
%     title('Replay firing rate','Units','normalized','Position',[0.5 0.92])
% 
%     %% plot flight firing rate pattern
%     axes(panels{1}(ii_ex,5))
%     cla reset
%     hold on
%     imagesc(bin_centers_time, 1:nCellsDir, FR_flight(sort_cells_IX,:))
%     scatter(bin_centers_time(1)-0.02*range(bin_centers_time), 1:nCellsDir,50,clrs(sort_cells_IX,:),'|','LineWidth',10,'clipping','off')
%     rescale_plot_data('x',[1e-6 seq_ts(1)])
%     xlabel('Time (s)')
%     ylabel('cell no.')
%     yticks(1:nCellsDir)
%     title('Flight firing rate','Units','normalized','Position',[0.5 0.92])

    %% link axes
    linkaxes(panels{2}(ii_ex,:),'y');
%     xlim([0 seq.duration]+ [-1 1].*0.05*seq.duration);
    ylim(seq_pos)
end

%% plot tuning curves in flight vs in replay
for ii_cell = 1:length(tuning_corr_cell_examples_list)
    axes(panels{3}(ii_cell))
    cla reset
    hold on
    cell_num = tuning_corr_cell_examples_list(ii_cell).cell_num;
    ii_dir = tuning_corr_cell_examples_list(ii_cell).ii_dir;
    cell = cell_load_data(cell_num, 'FR_map','replay_FR_map','details');
    exp=exp_load_data(cell.details.exp_ID,'details');
    x = cell.FR_map(ii_dir).all.bin_centers;
    y = cell.FR_map(ii_dir).all.PSTH;
    y = y./max(y);
    plot(x,y,'k-','LineWidth',2.4);
    x = cell.replay_FR_map.replay_PSTH_all(ii_cat,ii_dir).bin_centers;
    y = cell.replay_FR_map.replay_PSTH_all(ii_cat,ii_dir).PSTH;
    y = y./max(y);
    plot(x,y,'r-','LineWidth',1);
%     title("cell "+cell_num+" dir"+ii_dir)
    hax=gca;
    hax.XRuler.TickLength(1) = 0.02;
    hax.YRuler.TickLength(1) = 0.007;
    hax.XRuler.TickLabelGapOffset = -2;
    hax.YRuler.TickLabelGapOffset = 1;
    switch exp.details.recordingArena
        case '120m'
            xlim([0 135])
        case '200m'
            xlim([0 200])
    end
    ylim([0 1])
    xticks(0:50:200)
    yticks([0 1])
%     text(0.01,1.2,sprintf("cell %d dir %d",cell_num,ii_dir),'units','normalized','FontSize',6.5);
    if ii_cell==5 || ii_cell==6
        xlabel('Position (m)','Units','normalized','position',[0.5 -0.43])
    end
    if ii_cell == 5
        ylabel('Firing rate (norm.)','Units','normalized','position',[-0.1 2])
    end
    if ii_cell == 4
        x = 100+[0 6];
        y = [0.8 1.05];
        plot(x,y([2 2]), 'k-','LineWidth', 3,'Clipping','Off');
        plot(x,y([1 1]), 'r-','LineWidth', 1,'Clipping','Off');
        text(x(2)+2, y(2), 'Flight','HorizontalAlignment','left','FontSize',7,'VerticalAlignment','middle')
        text(x(2)+2, y(1), 'Replay','HorizontalAlignment','left','FontSize',7,'VerticalAlignment','middle')
    end
end

%% single unit FR corr hist
axes(panels{4}(1))
cla reset
hold on
FR_corr = load("E:\Tamir\work\PROJECTS\LargeScale\paper_replay\figures\single_unit_replay_tuning_corr.mat");
x = cat(1,FR_corr.ccc_data{:});
x_shuffle = cat(1,FR_corr.ccc_shuffle{:});
histogram(x,'FaceColor',0.15*[1 1 1],'normalization','pdf','NumBins',21,'EdgeColor','k','LineWidth',1);
histogram(x_shuffle,'EdgeColor','k','DisplayStyle','Stairs','normalization','pdf','LineWidth',1.5);
[~,P_KS, KS_stats] = kstest2(x, x_shuffle);
[~,P_ttest, ~, ttest_stats] = ttest2(x, x_shuffle);
fprintf('FR corr distribution, data vs shuffle, ks-test: p=%g\n',P_KS);
fprintf('FR corr distribution, data vs shuffle, t-test: p=%g\n',P_ttest);
text(0.53,0.88,['{\itP}_{KS} = ' sprintf('%.2g',P_KS)],'units','normalized','FontSize',7,'HorizontalAlignment','center');
text(0.53,0.95,['{\itP}_{t} = ' sprintf('%.2g',P_ttest)],'units','normalized','FontSize',7,'HorizontalAlignment','center');
xlim([-1 1])
xlabel({'Correlation of firing rate maps (r)';'during flight vs. during replay'})
ylabel('Probability density')
hax=gca;
hax.XRuler.TickLength(1) = 0.02;
hax.YRuler.TickLength(1) = 0.007;
hax.XRuler.TickLabelGapOffset = -2;
hax.YRuler.TickLabelGapOffset = 1;

%% add legend
axes(panels{4}(2))
cla reset
hold on
axis off
axis equal
patch([0 0 1 1],[0 1 1 0], 0.15*[1 1 1],'EdgeColor','k');
plot([0 1],[1 1].*2, 'k-','LineWidth',1.5)
text(1.5,0.5, 'Data','FontSize',7)
text(1.5,2, 'Shuffle','FontSize',7)
axis ij

%% add panel letters
font_size = 11;
axes(panels{1}(1))
text(-0.08,1.0, 'a', 'Units','normalized','FontWeight','bold','FontSize',font_size);

axes(panels{2}(1,1))
text(-0.26,1.1, 'b', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{2}(2,1))
text(-0.26,1.1, 'c', 'Units','normalized','FontWeight','bold','FontSize',font_size);

axes(panels{3}(1,1))
text(-0.17,1.12, 'd', 'Units','normalized','FontWeight','bold','FontSize',font_size);

axes(panels{4}(1))
text(-0.2,1.1, 'e', 'Units','normalized','FontWeight','bold','FontSize',font_size);

%%
fig_name = sprintf('%s_ex_1_FE_%d_ex_2_FE_%d', ...
    fig_name_str, ...
    seq_examples_list(1).ii_FE,...
    seq_examples_list(2).ii_FE)
% fig_name = sprintf('%s_panel_d_opt_%d',fig_name,panel_d_opt)
file_out = fullfile(res_dir, fig_name);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
disp('figure saved!')

%%
