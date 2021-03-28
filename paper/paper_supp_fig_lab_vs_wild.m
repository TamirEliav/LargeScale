%% Large Scale - fig supp - Lab born vs Wild bats
%%
clear 
clc

%% plotting options
field_speed_opt = 1;
% corr_type = 'pearson';
corr_type = 'spearman';
match_TT_pos = true;

% examples_option = 5; %1:5
hist_normalization_option = 'probability'; %count
hist_Yscale_option = 'log';

color_Lab = [0    0.75    0];%[0.4940    0.1840    0.5560];
color_Wild = [0.3 0.25 0.2];%color_wild = [0.65 0.6 0.65];
color_wild_text = [0.45 0.4 0.45];
Lab_alpha = 0.35;
Wild_alpha = 0.6;

%% add data folders to path
addpath('L:\processed_data_structs');

%% define output files
res_dir =  'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S16';%'Sup_lab_and_wild';
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
figure ;
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
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');

pause(0.2); % workaround to solve matlab automatically changing the axes positions...

% create panels

panel_A(1)    = axes('position', [ 2  18  6  6]);
panel_A(2)    = axes('position', [ 8.7  18  6.1  6]);
% panel_B    = axes('position', [ 10  19  3  3]);
% panel_C    = axes('position', [ 15  19  3  3]);
% panel_D    = axes('position', [ 2  14  6  3]);
panel_B    = axes('position', [ 2.5  14  3  3]);
panel_C    = axes('position', [ 7.5  14  3  3]);
panel_D    = axes('position', [ 12.5  14  6  3]);

%% mamtak room image
image_filename = 'L:\Resources\IMG_3678_brightness35.tif';
axes(panel_A(1));
image1 = imread(image_filename);
imshow(image1(:,:,1:3));
text(-0.17,1.1, 'A', 'Units','normalized','FontWeight','bold');

image_filename = 'L:\Resources\Bat_neurobiology_lab_Moment(11)_cropped_brightness150.tif';
axes(panel_A(2));
image2 = imread(image_filename);
imshow(image2(:,:,1:3));

%% load data

prm = PARAMS_GetAll();
% cells_t = DS_get_cells_summary();
% bats = [9845, 102, 194];
% cells_t(~ismember(cells_t.bat, bats ),:) = [];
% cells_t([cells_t.remove]==1,:) = [];%cells chosen to remove because of repititions between days
% cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',0);
% cells = [cells{:}];
% cells = [cells.details];
% cells(~contains({cells.brain_area}, 'CA1')) = [];
% cells(~ismember([cells.ClusterQuality], [2])) = [];
% % remove cells based on FR:
% cells = cellfun(@(c)(cell_load_data(c,'details','meanFR')), {cells.cell_ID}, 'UniformOutput',0);
% cells = [cells{:}];
% cells_details = [cells.details];
% cells_ID = {cells_details.cell_ID};
% meanFR = [cells.meanFR];
% cells_ID([meanFR.all]>prm.inclusion.interneuron_FR_thr)=[];
% % remove non-signif cells:
% cells = cellfun(@(c)(cell_load_data(c,'details','signif')), cells_ID, 'UniformOutput',0);
% cells = [cells{:}];
% signif = cat(1,cells.signif);
% signif = arrayfun(@(x)(x.TF), signif);
% signif = any(signif,2);
% cells_ID(~signif) = [];
% % load exps for significant cells only:
% cells = cellfun(@(c)(cell_load_data(c,'details')), cells_ID, 'UniformOutput',0);
% cells = [cells{:}];
% cells_details = [cells.details];
% exps_ID = {cells_details.exp_ID};
% exps_ID = unique(exps_ID);
% exps = cellfun(@(c)(exp_load_data(c,'details','flight','balls')), exps_ID, 'UniformOutput',0);
% exps = [exps{:}];
% save('Exps_lab.mat','exps');

exps = load('Exps_lab.mat');
exps_lab = exps.exps;
exps = load('exps_including_FE.mat');
exps_wild = exps.exps;


%% Panel B - total distance
axes(panel_B); hold on;
text(-0.5,1.2, 'B', 'Units','normalized','FontWeight','bold');

distance_exp_lab = zeros([1,length(exps_lab)]);
for ii_exp = 1:length(exps_lab)
    distance_exp_lab(ii_exp) = sum([exps_lab(ii_exp).flight.FE(:).distance]);
    total_time_exp_lab(ii_exp) = (exps_lab(ii_exp).flight.FE(end).end_ts ...
        - exps_lab(ii_exp).flight.FE(1).start_ts) * (10^-6)/3600; %in hours;
end

distance_exp_wild = zeros([1,length(exps_wild)]);
for ii_exp = 1:length(exps_wild)
    distance_exp_wild(ii_exp) = sum([exps_wild(ii_exp).flight.FE(:).distance]);
    total_time_exp_wild(ii_exp) = (exps_wild(ii_exp).flight.FE(end).end_ts ...
        - exps_wild(ii_exp).flight.FE(1).start_ts) * (10^-6)/3600; %in hours
end
    
x = [1,2];
dis_per_hr_lab = (distance_exp_lab./total_time_exp_lab)/1000;
dis_per_hr_wild = (distance_exp_wild./total_time_exp_wild)/1000;

y = [mean(dis_per_hr_lab),mean(dis_per_hr_wild)];
hb = bar(x,y);
hb.FaceColor = 'flat';
hb.CData = [color_Lab; 0.5*[1 1 1]];
sem_lab = std(dis_per_hr_lab)/sqrt(length(dis_per_hr_lab));
sem_wild = std(dis_per_hr_wild)/sqrt(length(dis_per_hr_wild));
he=errorbar(x,y,[sem_lab sem_wild]);
he.CapSize = 2;
he.LineStyle = 'none';
he.Color = 'k';
he.LineWidth = 0.7;
h=gca;
h.XTick = [1 2];
h.XTickLabels = {'Lab','Wild'};
h.TickLength = [0.03 0.03];    
ylim([0,15]);
ylabel({'Total flight distance','per 1 hour session (km)'});

[~,P_ttest,~,test_stats] = ttest2(dis_per_hr_lab, dis_per_hr_wild);
% text(1.2,1, sprintf('P_{T} = %.02f',P_ttest),'Units','normalized','FontSize',7,'HorizontalAlignment','right');
if P_ttest>=0.05
    signif_str = 'n.s.';
    signif_str_font_size = 7;
%     annotation('line',Xf,Yf,'LineWidth',1.1);
    line([1,2],[13,13],'color','k','LineWidth',1.1);
    text(1.5,14, signif_str, 'HorizontalAlignment','center','FontSize',signif_str_font_size);
end

%% Panel C - Flight speed
axes(panel_C); hold on;
text(-0.4,1.2, 'C', 'Units','normalized','FontWeight','bold');

pos_lab = [];
vel_lab = [];
valid_pos_lab = [];
valid_vel_lab = [];
vel_lab_per_exp = zeros([1,length(exps_lab)]);
vel_lab_per_exp_entireTunnel = zeros([1,length(exps_lab)]);

for ii_exp = 1:length(exps_lab)
    FE = exps_lab(ii_exp).flight.FE;
    FE([FE.distance]<prm.flight.full_min_distance) = [];
    pos_tmp = [FE.pos];
    vel_tmp = [FE.vel];
    pos_lab = [pos_lab [FE.pos]];
    vel_lab = [vel_lab [FE.vel]];
    balls = exps_lab(ii_exp).balls;
    valid_area = balls + prm.fields.valid_speed_pos;
    valid_ind = and(pos_tmp>valid_area(1), pos_tmp<valid_area(2));
    valid_pos_lab = [valid_pos_lab pos_tmp(valid_ind)];
    valid_vel_lab = [valid_vel_lab vel_tmp(valid_ind)];
    
    vel_lab_per_exp(ii_exp) = mean(abs(vel_tmp(valid_ind)));
    vel_lab_per_exp_entireTunnel(ii_exp) = mean(abs([FE.vel]));
%     vel_lab_per_exp(ii_exp) = median(abs([FE.vel]));
end
vel_lab = abs(vel_lab);
valid_vel_lab = abs(valid_vel_lab);

pos_wild = [];
vel_wild = [];
valid_pos_wild = [];
valid_vel_wild = [];
vel_wild_per_exp = zeros([1,length(exps_wild)]);
 vel_wild_per_exp_entireTunnel = zeros([1,length(exps_wild)]);
for ii_exp = 1:length(exps_wild)
    FE = exps_wild(ii_exp).flight.FE;
    FE([FE.distance]<prm.flight.full_min_distance) = [];
    pos_tmp = [FE.pos];
    vel_tmp = [FE.vel];
    pos_wild = [pos_wild [FE.pos]];
    vel_wild = [vel_wild [FE.vel]];
    valid_area = [10 187.5];
    valid_ind = and(pos_tmp>valid_area(1), pos_tmp<valid_area(2));
    valid_pos_wild = [valid_pos_wild pos_tmp(valid_ind)];
    valid_vel_wild = [valid_vel_wild vel_tmp(valid_ind)];
    vel_wild_per_exp(ii_exp) = mean(abs(vel_tmp(valid_ind)));
    vel_wild_per_exp_entireTunnel(ii_exp) = mean(abs([FE.vel]));
    
%     vel_wild_per_exp(ii_exp) = median(abs([FE.vel]));
end
vel_wild = abs(vel_wild);
valid_vel_wild = abs(valid_vel_wild);

y = [mean(vel_lab_per_exp),mean(vel_wild_per_exp)];
hb = bar(x,y);
hb.FaceColor = 'flat';
hb.CData = [color_Lab; 0.5*[1 1 1]];
sem_lab = std(vel_lab_per_exp)/sqrt(length(vel_lab_per_exp));
sem_wild = std(vel_wild_per_exp)/sqrt(length(vel_wild_per_exp));
he=errorbar(x,y,[sem_lab sem_wild]);
he.CapSize = 2;
he.LineStyle = 'none';
he.Color = 'k';
he.LineWidth = 0.7;
h=gca;
h.XTick = [1 2];
h.XTickLabels = {'Lab','Wild'};
h.TickLength = [0.03 0.03];    
ylim([0,10]);
ylabel('Flight speed (m/s)'); 

[~,P_ttest,~,test_stats] = ttest2(vel_lab_per_exp, vel_wild_per_exp);
% text(1.2,0.9, sprintf('P = %0.00e',P_ttest),'Units','normalized','FontSize',7,'HorizontalAlignment','right');
if P_ttest<10^-6
    signif_str = '*****';
    signif_str_font_size = 12;
%     annotation('line',Xf,Yf,'LineWidth',1.1);
    line([1,2],[8.75,8.75],'color','k','LineWidth',1.1);
    text(1.5,9.1, signif_str, 'HorizontalAlignment','center','FontSize',signif_str_font_size);
end

%% Panel D - Flight speed vs position
axes(panel_D)
cla
hold on
text(-0.17,1.2, 'D', 'Units','normalized','FontWeight','bold');

bin_size = 1; %m
bin_edges = 0:bin_size:200;
bin_centers = bin_edges(1:end-1) + bin_size/2;

pos_wild_binned = discretize(pos_wild,bin_edges);
[pos_wild_groups, group_wild_bin] = findgroups(pos_wild_binned);
vel_group_stats_wild = splitapply(@(x) [mean(x) std(x) (std(x)/sqrt(length(x)))], vel_wild', pos_wild_groups');
% vel_group_stats_wild = splitapply(@(x) [median(x) prctile(x,25) prctile(x,75)], vel_wild', pos_wild_groups');
IX_begining = group_wild_bin<6;
vel_group_stats_wild(IX_begining,:) = [];
group_wild_bin(IX_begining) = [];
plot_stdshade(bin_centers(group_wild_bin), vel_group_stats_wild(:,1)', vel_group_stats_wild(:,2)' , 0.3, [0.5 0.5 0.5])
% plot_QIshade(bin_centers(group_wild_bin),vel_group_stats_wild(:,1)',vel_group_stats_wild(:,2)',vel_group_stats_wild(:,3)',0.3, [0.5 0.5 0.5])

pos_lab_binned = discretize(pos_lab,bin_edges);
[pos_lab_groups, group_lab_bin] = findgroups(pos_lab_binned);
vel_group_stats_lab = splitapply(@(x) [mean(x) std(x) (std(x)/sqrt(length(x)))], vel_lab', pos_lab_groups');
% vel_group_stats_lab = splitapply(@(x) [median(x) prctile(x,25) prctile(x,75)], vel_lab', pos_lab_groups');
IX_begining = group_lab_bin<6;
vel_group_stats_lab(IX_begining,:) = [];
group_lab_bin(IX_begining) = [];
plot_stdshade(bin_centers(group_lab_bin), vel_group_stats_lab(:,1)', vel_group_stats_lab(:,2)' , 0.3, color_Lab)
% plot_QIshade(bin_centers(group_lab_bin),vel_group_stats_lab(:,1)',vel_group_stats_lab(:,2)',vel_group_stats_lab(:,3)',0.3, color_Lab)

ylabel('Flight speed (m/s)'); 
xlabel('Position (m)'); 


%%
fig_name_out = fullfile(res_dir, fig_name_str);

% fig_name_out = fullfile(res_dir, sprintf('%s__corr_%s_%d_paramset_%d',fig_name_str,corr_type,field_speed_opt,prm.parmaset));
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');

diary off
    
    
    
    
    
    
    
    
    
    
    
    