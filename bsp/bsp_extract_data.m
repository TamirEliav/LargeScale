function bsp_extract_data(main_dir)


%%
% clear all
% clc

%% read files
% main_dir = 'D:\Tamir\DATA\9892_Rocky\20160306\Bsp';
% main_dir = 'P:\BeSPoon\testing\20160316__timestamps_bug_version_test2_LIORA\bsp\dist';
dist_csv_file_names = dir(fullfile(main_dir,['*_distances.csv']));
TTL_ts_ns_all = [];
bsp_data = struct('tag_ID',[], 'pos', [], 'ts_ns', []);
for ii_file=1:length(dist_csv_file_names)
    file_name = fullfile(main_dir,dist_csv_file_names(ii_file).name)
    
    % replace ';' with ','
    fid  = fopen(file_name,'r');
    f=fread(fid,'*char')';
    fclose(fid);
    f = strrep(f,';',',');
    file_name = strrep(file_name, 'distances.csv', 'distances_COMMA.csv');
    fid  = fopen(file_name,'w');
    fprintf(fid,'%s',f);
    fclose(fid);
    
    [NUM,TXT,RAW]=xlsread(file_name);
    
    % read header
    num_tags = RAW{2,2};
    for ii_tag=1:num_tags
        clear tag_id ts_ns X Y Z pos
        tag_id = RAW{3,(ii_tag-1)*3+2};
        tag_id = tag_id(6:end);
        tag_id = str2num(tag_id);
        ts_ns = NUM(4:end,1);
        X = NUM(4:end,(ii_tag-1)*3+2);
        Y = NUM(4:end,(ii_tag-1)*3+3);
        Z = NUM(4:end,(ii_tag-1)*3+4);
        pos = [X Y Z];
        nan_IX = find(isnan(X));
        pos(nan_IX,:) = [];
        ts_ns(nan_IX) = [];
        
        % insert data to relevant tag (if tag entry does not exist, create new one)
        if ismember(tag_id, [bsp_data.tag_ID])
            tag_IX = find([bsp_data.tag_ID] == tag_id);
        else
            tag_IX = length([bsp_data.tag_ID])+1;
            bsp_data(tag_IX).tag_ID = tag_id;
        end
        bsp_data(tag_IX).pos = [bsp_data(tag_IX).pos; pos];
        bsp_data(tag_IX).ts_ns = [bsp_data(tag_IX).ts_ns; ts_ns];
    end
    
    % add TTL (if exist)
    TTL_col = 1+num_tags*3+1;
    if size(NUM,2)>= TTL_col
        TTL_IX = find(NUM(4:end, TTL_col) == 1);
        ts_ns = NUM(4:end,1);
        TTL_ts_ns = ts_ns(TTL_IX);
        TTL_ts_ns_all = [TTL_ts_ns_all; TTL_ts_ns];
    else
        warning('No TTL data in current file, CHECK WHY!!!!')
    end
    
end

%% save as mat file
bsp_TTL_ts_ns = TTL_ts_ns_all;
save(fullfile(main_dir,'bsp_data'), 'bsp_data');
save(fullfile(main_dir,'bsp_TTL'), 'bsp_TTL_ts_ns');
for ii_tag = 1:length([bsp_data.tag_ID])
    tag_id = bsp_data(ii_tag).tag_ID;
    bsp_pos = bsp_data(ii_tag);
    save(fullfile(main_dir,['bsp_pos_tag_' num2str(tag_id)]), 'bsp_pos');
end

%% plot some dummy localization figures (TODO: move this to separate code, should be in the PRE-processing code)
for ii_tag = 1:length([bsp_data.tag_ID])

    bsp_pos = bsp_data(ii_tag);

    %%
    pos_diff = diff( [[0 0] ;bsp_pos.pos(:,1:2)] );
    time_diff_sec = diff( [0; bsp_pos.ts_ns]).*1e-9;
    vel_m_sec = pos_diff./repmat(time_diff_sec,1,2);
    abs_vel_m_sec = sqrt(    vel_m_sec(:,1).^2   +  vel_m_sec(:,2).^2  );
    abs_vel_m_sec(abs_vel_m_sec==inf) = nan;
    % % % 
    % % % bad_IX = find(isnan(abs_vel_m_sec) | abs_vel_m_sec>15);
    % % % good_IX = setdiff(1:length(ts_ns_all), bad_IX);
    % % % good_IX = 1:length(ts_ns_all);
    
    %%
    ax_h = [];
    figure
    figure_size_cm = [40 20];
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 figure_size_cm]);
    set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[5 5 0 0]); % position on screen...
    p = panel();
    p.pack('h',[0.8 0.2]);
    p(1).pack('v',[0.25 0.75]);
    p(2).pack('v',[0.5 0.5]);
    p.select('all')

    % 2D coverage (X vs Y)
    p(1,1).select(); hold on;
    ax_h(1) = gca;
    plot(bsp_pos.pos(:,1), bsp_pos.pos(:,2) );
    xlabel('X(m)')
    ylabel('Y(m)')
%     axis equal
    xlim([1315 1500])
    ylim([2460 2510])
    title('2D position')

	% 1D coverage (X vs time)
    p(1,2).select(); hold on;
    ax_h(2) = gca;
    plot(bsp_pos.pos(:,1), bsp_pos.ts_ns*1e-9/60, '.')
    xlim([1300 1550])
    xlabel('X(m)')
    ylabel('Time(minutes)')
    linkaxes(ax_h,'x')

    % cummulative distance
    p(2,1).select(); hold on;
    cum_dist = cumsum(abs_vel_m_sec.* time_diff_sec) ;
    cum_dist = cum_dist * 1e-3;
    cum_dist = cum_dist - cum_dist(1);
    t = bsp_pos.ts_ns;
    t = t - t(1);
    t = t * 1e-9 / 60;
    plot(t, cum_dist,'.');
    xlabel('Time (minutes)');
    ylabel('Cummulative distance (km)');
    
    % velocity histogram
    p(2,2).select(); hold on;
    hist(abs_vel_m_sec, [0:0.25:20]);
    xlabel('Velocity (m/sec)');
    ylabel('Counts');
    xlim([0 20])
    
    saveas(gcf, fullfile(main_dir,['bsp_pos_tag_' num2str(bsp_pos.tag_ID)]), 'jpeg')
    saveas(gcf, fullfile(main_dir,['bsp_pos_tag_' num2str(bsp_pos.tag_ID)]), 'fig')
end



%%
% [n x]=hist(abs_vel_m_sec,[0:0.5:30]);
% bar(x,n);xlim([0 20])
% IX = find(n>200);
% median(diff(x(IX)))
% fs = 1e9/nanmedian(diff(ts_ns_all))
% Ts = 1/fs












