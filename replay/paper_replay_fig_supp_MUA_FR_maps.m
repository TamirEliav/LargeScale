%% Replay - Fig supp 8 - Firing rate maps
clear 
clc
close all

%% data options 

%% plotting options

%% graphics params

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Extended_Data_Fig_8';
fig_caption_str = 'MUA_FR_maps';
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
panels{1}(1,1) = axes('position', [2.5 20 17 4]);
panels{1}(1,2) = axes('position', [2.5 15 17 4]);
panels{1}(2,1) = axes('position', [2.5 9 17 4]);
panels{1}(2,2) = axes('position', [2.5 4 17 4]);

%% load data
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
FR_maps = [];
for ii_exp = 1:length(exp_list)
    exp_ID = exp_list{ii_exp};
    exp = exp_load_data(exp_ID,'MUA_FR_map');
    FR_maps(ii_exp,:,:) = exp.MUA_FR_map.maps;
end

%%
dir_arrow_str_map = containers.Map([1 2],{'\rightarrow','\leftarrow'});

%%
for ii_dir = 1:2
    map_dir_str = sprintf('Flight direction %d %s',ii_dir,dir_arrow_str_map(ii_dir));

    axes(panels{1}(1,ii_dir));
    cla reset
    hold on
    x = linspace(0,1,size(FR_maps,3));
    y = squeeze(FR_maps(:,ii_dir,:))';
    plot(x,y);
    plot(x,mean(y,2),'-k','LineWidth',2)
    xlabel('Position (norm.)')
    ylabel('MUA firing rate (Hz)')
    text(0.1,0.8,map_dir_str,'units','normalized');

    axes(panels{1}(2,ii_dir));
    cla reset
    hold on
    x = linspace(0,1,size(FR_maps,3));
    y = squeeze(FR_maps(:,ii_dir,:))';
    y = normalize(y,'zscore');
    plot(x,y)
    plot(x,mean(y,2),'-k','LineWidth',2)
    xlabel('Position (norm.)')
    ylabel('MUA firing rate (z)')
    text(0.1,0.95,map_dir_str,'units','normalized');
end

%% add panel letters
font_size = 11;
axes(panels{1}(1,1))
text(-0.08,1.15, 'a', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{1}(2,1))
text(-0.08,1.15, 'b', 'Units','normalized','FontWeight','bold','FontSize',font_size);

%%
fig_name = sprintf('%s',fig_name_str);
file_out = fullfile(res_dir, fig_name);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
disp('figure saved!')


