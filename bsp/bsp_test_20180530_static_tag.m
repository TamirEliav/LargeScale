%% brspoon static tag test (30/05/2018)

%%
clear 
clc

%% data params
% tags from highest to lowest: 1751 -> 68cm -> 2659 -> 59cm -> 4013 ***
% all height meausrements will be from lowest tag 4013 ***
tag_IDs = [1751 2659 4013];
tag_offsets = [0 59 68]./100; % in meters
tag_offsets = cumsum(tag_offsets);
tag_offsets = flip(tag_offsets);
tag_colors = [...
    1 0 0;...
    1 1 0;...
    0 1 0;...
    ];

%% sync nlg comments with bsp data
main_dir = 'L:\BeSpoon\testing\test_20180530__YOM_KEF_200m_static+dynamic+discretization\30-05-2018__calib_test_static';
nlg_dir = fullfile(main_dir, 'nlg');
nlx_dir = fullfile(main_dir, 'nlx');
bsp_dir = fullfile(main_dir, 'bsp');
sync_dir = fullfile(main_dir, 'sync');
results_dir = fullfile(main_dir, 'results');
mkdir(results_dir)
% Nlg2Nlx(main_dir);
% PRE_sync_bsp_to_nlg(bsp_dir, nlx_dir, sync_dir);

%% load recording comments
T = readtable('L:\BeSpoon\testing\test_20180530__YOM_KEF_200m_static+dynamic+discretization\30-05-2018__calib_test_static\static_test.xlsx');

%% load bsp data
clear tags_data
for ii_tag = 1:length(tag_IDs)
    bsp_data_file = fullfile(bsp_dir, sprintf('bsp_pos_tag_%d',tag_IDs(ii_tag)));
    load( bsp_data_file )
    tags_data{ii_tag} = bsp_pos;
end

%% divide data to segments on the points and 
tag_pos_at_points = {};
tag_pos_at_points_mean = [];
tag_pos_at_points_std = [];
tag_pos_at_points_err = [];
points_locs = [T.loc_X T.loc_Y T.loc_Z];
clear raw_err
is_raw_err_abs = 0;
if is_raw_err_abs
    raw_err_str = 'abs_error';
else
    raw_err_str = 'error';
end
for ii_tag = 1:length(tag_IDs)
    for ii_point = 1:length(T.point)
        IX = find( ...
            tags_data{ii_tag}.ts_nlg_usec > T.t_start(ii_point)*1e3 &...
            tags_data{ii_tag}.ts_nlg_usec < T.t_end(ii_point)*1e3    ...
            );
        tag_pos_at_points{ii_tag, ii_point} = tags_data{ii_tag}.pos(IX,:);
        tag_pos_at_points_mean(ii_tag, ii_point,:) = mean(tag_pos_at_points{ii_tag, ii_point});
        tag_pos_at_points_std(ii_tag, ii_point,:) = std(tag_pos_at_points{ii_tag, ii_point});
        ref_loc = points_locs(ii_point,:);
        ref_loc(3) = ref_loc(3) + T.offset_z(ii_point) + tag_offsets(ii_tag);
        tag_pos_at_points_err(ii_tag, ii_point,:) = squeeze(tag_pos_at_points_mean(ii_tag, ii_point,:)) - ref_loc';
        
        npoints = size(tag_pos_at_points{ii_tag, ii_point},1);
        raw_err{ii_tag, ii_point} = tag_pos_at_points{ii_tag, ii_point} - repmat(ref_loc,npoints,1);
        if is_raw_err_abs
            raw_err{ii_tag, ii_point} = abs(raw_err{ii_tag, ii_point});
        end
    end
end
raw_err_pooled = cat(1,raw_err{:});

%% plot x/y position
figure
p = panel();
plot( tag_pos_at_points_mean(1,:,1), tag_pos_at_points_mean(1,:,2),'.-')
axis equal

%% plot Relative heights
figure
hold on
ii_dim = 3;
x = (1:length(T.point))';
y = squeeze(tag_pos_at_points_mean(:,:,ii_dim))';
y = y - repmat(y(:,3),1,3); % remove the height of the lower tag
h=bar(x,y);
for ii_tag=1:3
    h(ii_tag).FaceColor = tag_colors(ii_tag,:);
    href = refline(0, tag_offsets(ii_tag));
    href.Color = tag_colors(ii_tag,:);
    href.LineWidth = 1.2;
end
title('Relative heights')
xlabel('#point')
ylabel('Relative height to lower tag (m)')
legend({...
    'upper tag';'middle tag';'lower tag';...
    'upper tag (real)';'middle tag (real)';'lower tag (real)';...
    });

filename = fullfile(results_dir, 'relative_heights');
saveas(gcf, filename, 'tif')
saveas(gcf, filename, 'fig')

%% plot Accuracy + precision (bars+errorbars)
figure
set(gcf, 'units', 'normalized', 'position', [0 0 1 1])
p = panel();
p.pack('v',3)
p.de.margin = 5;
p.margin = [30 25 10 15];
dim_strs = 'XYZ';
for ii_dim = 1:3
    p(ii_dim).select();
    barvalues = squeeze(tag_pos_at_points_err(:,:,ii_dim))';
    errors = squeeze(tag_pos_at_points_std(:,:,ii_dim))';
    width = [];
    groupnames = [];
    bw_title = [];
    bw_xlabel = [];
    bw_ylabel = [];
    bw_colormap = tag_colors;
    gridstatus = [];
    bw_legend = [];
    error_sides = [];
    legend_type = [];
    handles = barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides, legend_type);
    ylabel( sprintf('%s (m)',dim_strs(ii_dim)) );
%     h.Position = [-0.05 0.5];
%     x = (1:length(T.point))';
%     y = squeeze(tag_pos_at_points_err(:,:,ii_dim))';
%     err = squeeze(tag_pos_at_points_std(:,:,ii_dim))';
%     bar(x,y);
%     hold on 
%     errorbar(x,y,err)
end
xlabel('#point');
h=p.title('Accuracy (bars) + precision (errorbars)');
h.Position = [0.5 1.01 1];
h.FontSize = 24;

filename = fullfile(results_dir, 'Accuracy+precision');
saveas(gcf, filename, 'tif')
saveas(gcf, filename, 'fig')


%% plot accuracy alone as bars
figure
set(gcf, 'units', 'normalized', 'position', [0 0 1 1])
p = panel();
p.pack('v',3)
p.de.margin = 5;
p.margin = [30 25 10 15];
dim_strs = 'XYZ';
for ii_dim = 1:3
    p(ii_dim).select();
    x = (1:length(T.point))';
    y = squeeze(tag_pos_at_points_err(:,:,ii_dim))';
    h=bar(x,y);
    for ii_h = 1:length(h)
        h(ii_h).FaceColor = tag_colors(ii_h,:);
    end
    ylabel( sprintf('%s (m)',dim_strs(ii_dim)) );
end
xlabel('#point');
h=p.title('Accuracy (bars)');
h.Position = [0.5 1.01 1];
h.FontSize = 24;

filename = fullfile(results_dir, 'Accuracy');
saveas(gcf, filename, 'tif')
saveas(gcf, filename, 'fig')

%% plot precision alone as bars
figure
set(gcf, 'units', 'normalized', 'position', [0 0 1 1])
p = panel();
p.pack('v',3)
p.de.margin = 5;
p.margin = [30 25 10 15];
dim_strs = 'XYZ';
for ii_dim = 1:3
    p(ii_dim).select();
    x = (1:length(T.point))';
    y = squeeze(tag_pos_at_points_std(:,:,ii_dim))';
    h=bar(x,y);
    for ii_h = 1:length(h)
        h(ii_h).FaceColor = tag_colors(ii_h,:);
    end
    ylabel( sprintf('%s (m)',dim_strs(ii_dim)) );
end
xlabel('#point');
h=p.title('Precision (bars)');
h.Position = [0.5 1.01 1];
h.FontSize = 24;

filename = fullfile(results_dir, 'Precision');
saveas(gcf, filename, 'tif')
saveas(gcf, filename, 'fig')


%% scatter of accuracy in X vs. Y vs. Z
figure

err_x = tag_pos_at_points_err(:,:,1);
err_y = tag_pos_at_points_err(:,:,2);
err_z = tag_pos_at_points_err(:,:,3);

subplot(2,2,1)
plot(err_x(:), err_z(:), '.')
xlabel('X Accuracy (m)');
ylabel('Z Accuracy (m)');

subplot(2,2,2)
plot(err_y(:), err_z(:), '.')
xlabel('Y Accuracy (m)');
ylabel('Z Accuracy (m)');

subplot(2,2,3)
plot(err_x(:)+err_y(:), err_z(:), '.')
xlabel('X+Y Accuracy (m)');
ylabel('Z Accuracy (m)');

subplot(2,2,4)
plot(err_x(:), err_y(:), '.')
xlabel('X Accuracy (m)');
ylabel('Y Accuracy (m)');

filename = fullfile(results_dir, 'scatter_accuracy_X_vs_Y_vs_Z');
saveas(gcf, filename, 'tif')
saveas(gcf, filename, 'fig')

%% precision vs. accuracy
figure
hold on
dim_sym = 'o*d';
for ii_dim = 1:3
    accuracy = abs(tag_pos_at_points_err(:,:,ii_dim));
    precision = tag_pos_at_points_std(:,:,ii_dim);
    plot( accuracy(:), precision(:), dim_sym(ii_dim));
end
xlabel('Accuracy (m)')
ylabel('Precision (m)')
title('precision vs. accuracy')

filename = fullfile(results_dir, 'precision_vs_accuracy');
saveas(gcf, filename, 'tif')
saveas(gcf, filename, 'fig')


%% raw error correlations (X vs. Y vs. Z)
figure
set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1])

subplot(2,2,1)
plot(raw_err_pooled(:,1), raw_err_pooled(:,2),'.')
title('Y vs. X')
xlabel('err (m)')
ylabel('err (m)')

subplot(2,2,2)
plot(raw_err_pooled(:,1), raw_err_pooled(:,3),'.')
title('Z vs. X')
xlabel('err (m)')
ylabel('err (m)')

subplot(2,2,3)
plot(raw_err_pooled(:,2), raw_err_pooled(:,3),'.')
title('Z vs. Y')
xlabel('err (m)')
ylabel('err (m)')

subplot(2,2,4)
plot(sqrt(sum((raw_err_pooled(:,1:2)).^2,2)), raw_err_pooled(:,3),'.')
title('Z vs. X/Y')
xlabel('err (m)')
ylabel('err (m)')

suptitle({'Raw measurements error correlations';'pooled acroos locations/tags'});

filename = ['Raw_measurements_' raw_err_str  '_correlations_pooled'];
fileout = fullfile(results_dir , filename);
saveas(gcf, fileout, 'tif')
saveas(gcf, fileout, 'fig')

%% raw error correlations (X vs. Y vs. Z) - by locations X tags
figure
for ii_loc = 1:size(raw_err,2)
    clf
    p=panel();
    set(p.figure, 'Units', 'normalized', 'Position', [0 0 1 1])
    p.margin = [25 35 25 35];
    p.pack('h',4,'v',3)
    p.de.margin = 15;
    p.fontsize=16
    
    axis_str = 'axis square';
    
    for ii_tag = 1:size(raw_err,1)
        
        err = raw_err{ii_tag, ii_loc};
        
        p(1,ii_tag).select()
        plot(err(:,1), err(:,2), '.', 'Color', tag_colors(ii_tag,:))
        eval(axis_str)
%         xlabel('err (m)')
%         ylabel('err (m)')
        
        p(2,ii_tag).select()
        plot(err(:,1), err(:,3), '.', 'Color', tag_colors(ii_tag,:))
        eval(axis_str)
%         xlabel('err (m)')
%         ylabel('err (m)')
        
        p(3,ii_tag).select()
        plot(err(:,2), err(:,3), '.', 'Color', tag_colors(ii_tag,:))
        eval(axis_str)
%         xlabel('err (m)')
%         ylabel('err (m)')
        
        p(4,ii_tag).select()
        plot( sqrt(sum(err(:,1:2).^2,2)) , err(:,3), '.', 'Color', tag_colors(ii_tag,:))
        eval(axis_str)
%         xlabel('err (m)')
%         ylabel('err (m)')
    end
    h = p.title({'Raw measurements error correlations';sprintf('location #%d', ii_loc)}); pause(1); h.Position = [0.5 1.08];
    h = p.xlabel('Error (m)'); pause(0.1); h.Units = 'normalized'; h.Position = [0.5 -0.05];
    h = p.ylabel('Error (m)'); pause(0.1); h.Units = 'normalized'; h.Position = [-0.02 0.5];
    h = p(1).xlabel('Y vs. X'); pause(0.1); h.Units = 'normalized'; h.FontWeight = 'bold'; h.Position = [0.5 1.04];
    h = p(2).xlabel('Z vs. X'); pause(0.1); h.Units = 'normalized'; h.FontWeight = 'bold'; h.Position = [0.5 1.04];
    h = p(3).xlabel('Z vs. Y'); pause(0.1); h.Units = 'normalized'; h.FontWeight = 'bold'; h.Position = [0.5 1.04];
    h = p(4).xlabel('Z vs. X/Y'); pause(0.1); h.Units = 'normalized'; h.FontWeight = 'bold'; h.Position = [0.5 1.04];
    h = p(1,1).ylabel('Upper tag'); h.Units = 'normalized'; h.FontWeight = 'bold'; h.Position = [-0.25 0.5]; h.Color = tag_colors(1,:);
    h = p(1,2).ylabel('Middle tag'); h.Units = 'normalized'; h.FontWeight = 'bold'; h.Position = [-0.25 0.5]; h.Color = tag_colors(2,:);
    h = p(1,3).ylabel('Lower tag'); h.Units = 'normalized'; h.FontWeight = 'bold'; h.Position = [-0.25 0.5]; h.Color = tag_colors(3,:);
    filename = sprintf(['Raw_measurements_' raw_err_str '_correlations_loc#%d'],ii_loc);
    fileout = fullfile(results_dir , filename);
    saveas(gcf, fileout, 'tif')
    saveas(gcf, fileout, 'fig')
end

%%







