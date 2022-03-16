%% PhD thesis figure 4.6 - replay - over-representations

%%
close all
clear 
clc

%% plot options
hist_style = 'stairs';
% hist_style = 'bar';
switch hist_style 
    case 'stairs'
        hist_color_prop = 'EdgeColor';
    case 'bar'
        hist_color_prop = 'FaceColor';
end

%% options
% bats_to_include = [34 148 9861 2289];
% bats_to_include = [194 184 2382];
% bats_to_include = 34;
% bats_to_include = 148;
% bats_to_include = 9861;
% bats_to_include = 2289;
% bats_to_include = 194;
% bats_to_include = 184;
% bats_to_include = 2382;
%     bat_num    GroupCount
%     _______    __________
% 
%      2289           1    
%        34           2    
%      9861           4    
%       148          15    
%       184          17    
%       194          17    
%      2382          26    

%% define output files
res_dir = 'E:\Tamir\PhD\Thesis\resources\ch_4_seq';
mkdir(res_dir)
fig_name_str = 'Fig_4_6_replay_over_representations';
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
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');

% create panels
panels_size = [8 3];
panels_A(1,1) = axes('position', [2  19 panels_size]);
panels_A(2,1) = axes('position', [2  14.5 panels_size]);
panels_A(1,2) = axes('position', [12  19 panels_size]);
panels_A(2,2) = axes('position', [12  14.5 panels_size]);
panels_B(1) = axes('position', [2 9.5 3 3]);
panels_C(1) = axes('position', [7.5 9.5 panels_size]);

%% choose bats / sessions
[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
clear exp_list
groupsummary(T,'bat_num')
if exist('bats_to_include','var')
    T = groupfilter(T,"bat_num",@(x)ismember(x,bats_to_include),'bat_num');
end
bats = unique(T.bat_num)

%% load data
events = {};
clear MUA_maps;
for ii_exp = 1:height(T)
    % load exp data
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'details','MUA_FR_map');
    MUA_maps(ii_exp) = exp.MUA_FR_map;
    epoch_type = 'sleep';
    params_opt = 11;
    [events_session, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
    events{ii_exp} = events_session;
end
T.nEvents = cellfun(@length,events)'; % note this is without filtering sequence by features!
sortrows( groupsummary(T,'bat_num',["median","mean","max","sum"],"nEvents"),"sum_nEvents", 'descend')

%% pool events 
events = [events{:}];

%% apply inclusion criteria 
seqs = [events.seq_model];
[seqs, TF] = decoding_apply_seq_inclusion_criteria(seqs);
events(~TF)=[];

%% group by forward/reverse 
g = categorical([seqs.forward],[true false],{'Forward','Reverse'});

%% plot replay start/end position histograms (per forward/reverse replays)
nbins = 100;
cmap = brewermap(100,'RdBu');
clrs = cmap([0.9 0.1].*size(cmap,1),:);
directions = [1 -1];
for ii_dir = 1:length(directions)
    seqs_dir_TF = [seqs.state_direction]==directions(ii_dir);

    axes(panels_A(1,ii_dir))
    cla
    hold on
    x = [seqs.start_pos_norm];
    histogram(x(g=='Forward' & seqs_dir_TF),linspace(0,1,nbins+1),'Normalization','count',hist_color_prop,clrs(1,:),'DisplayStyle',hist_style);
    histogram(x(g=='Reverse' & seqs_dir_TF),linspace(0,1,nbins+1),'Normalization','count',hist_color_prop,clrs(2,:),'DisplayStyle',hist_style);
%     h=legend('Forward','Reverse');
%     h.Box='off';
    xlabel('Start position of replay (norm.)')
    ylabel('Counts')
    hax=gca;
    hax.XRuler.TickLabelGapOffset = -1;

    axes(panels_A(2,ii_dir))
    cla
    hold on
    x = [seqs.end_pos_norm];
    histogram(x(g=='Forward' & seqs_dir_TF),linspace(0,1,nbins+1),'Normalization','count',hist_color_prop,clrs(1,:),'DisplayStyle',hist_style);
    histogram(x(g=='Reverse' & seqs_dir_TF),linspace(0,1,nbins+1),'Normalization','count',hist_color_prop,clrs(2,:),'DisplayStyle',hist_style);
    xlabel('End position of replay (norm.)')
    ylabel('Counts')
    hax=gca;
    hax.XRuler.TickLabelGapOffset = -1;
end

% set ylimits
panels_A(1,1).YLim = [0 80];
panels_A(2,1).YLim = [0 150];
panels_A(1,2).YLim = [0 150];
panels_A(2,2).YLim = [0 220];

% linkaxes(panels_A(:),'y')

axes(panels_A(1,1))
text(-0.15,1.15, 'A', 'Units','normalized','FontWeight','bold');

for ii=1:numel(panels_A)
    pnl = panels_A(ii); axes(pnl);
    text(0,-0.27*range(pnl.YLim),{'Resting';'ball 1'},'HorizontalAlignment','center','FontSize',8);
    text(1,-0.27*range(pnl.YLim),{'Resting';'ball 2'},'HorizontalAlignment','center','FontSize',8);
end

% pointers to takeoff/landing in the histogram (dir 1)
TL_fnt_sz = 7;
h=annotation('textarrow');
h.Parent=panels_A(1,1);
h.X = 0.015*[1 1];
h.Y = 75+[0 -7];
h.String = ' Takeoffs';
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels_A(1,1);
h.X = 0.045*[1 1];
h.Y = 30+[0 -7];
h.String = {'           Reverse','           takeoffs'};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels_A(1,1);
h.X = 0.92.*[1 1];
h.Y = 55+[0 -7];
h.String = {'Landings        '};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels_A(1,1);
h.X = 0.985*[1 1];
h.Y = 65+[0 -7];
h.String = {'     Reverse';'     landings'};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels_A(2,1);
h.X = 0.015*[1 1];
h.Y = 65+[0 -12];
h.String = {'  Reverse','  takeoffs'};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels_A(2,1);
h.X = 0.1*[1 1];
h.Y = 40+[0 -12];
h.String = {'        Takeoffs'};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels_A(2,1);
h.X = 0.985*[1 1];
h.Y = 130+[0 -12];
h.String = 'Landings';
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels_A(2,1);
h.X = 0.95*[1 1];
h.Y = 60+[0 -12];
h.String = {'Reverse          ','landings          '};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

% pointers to takeoff/landing in the histogram (dir 2)
h=annotation('textarrow');
h.Parent=panels_A(1,2);
h.X = 0.015*[1 1];
h.Y = 143+[0 -12];
h.String = {'     Reverse','     landings'};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels_A(1,2);
h.X = 0.07*[1 1];
h.Y = 85+[0 -12];
h.String = '        Landings';
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels_A(1,2);
h.X = 0.99*[1 1];
h.Y = 72+[0 -12];
h.String = '      Takeoffs';
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
% h.Parent=panels_A(1,2);
h.Units='normalized';
h.X = 0.877+[0 0.01];
h.Y = 0.752+[0 -0.005];
h.String = {'Reverse','takeoffs'};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels_A(2,2);
h.X = 0.015*[1 1];
h.Y = 205+[0 -20];
h.String = '      Landings';
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels_A(2,2);
h.X = 0.07*[1 1];
h.Y = 85+[0 -20];
h.String = {'       Reverse','       landings'};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels_A(2,2);
h.X = 0.985*[1 1];
h.Y = 100+[0 -20];
h.String = 'Takeoffs';
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
% h.Parent=panels_A(1,2);
h.Units='normalized';
h.X = 0.873+[0 0.01];
h.Y = 0.567+[0 -0.005];
h.String = {'Takeoffs'};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

%% legend (dir 1)
pnl = panels_A(1,1);
axes(pnl);
h=annotation('textarrow');
h.Units='centimeters';
h.X = pnl.Position(1) + pnl.Position(3).*[0.6 0.8];
h.Y = pnl.Position(2) + pnl.Position(4).*[1.1 1.1];
h.Text.HorizontalAlignment = 'right';
h.String = 'For map of direction  '; h.FontSize = 9;

h=annotation('textarrow');
h.Units='centimeters';
h.X = pnl.Position(1) + pnl.Position(3).*[0.6 0.7];
h.Y = pnl.Position(2) + pnl.Position(4).*[0.9 0.9];
h.Color = clrs(1,:); h.String = 'Forward  '; h.FontSize = 9;
h.HeadLength = 5; h.HeadWidth = 8;

h=annotation('textarrow');
h.Units='centimeters';
h.X = pnl.Position(1) + pnl.Position(3).*[0.7 0.6];
h.Y = pnl.Position(2) + pnl.Position(4).*[0.73 0.73];
h.Color = clrs(2,:);
h.Text.HorizontalAlignment = 'right';
h.String = 'Reverse           '; h.FontSize = 9;
h.HeadLength = 5; h.HeadWidth = 8;

%% legend (dir 2)
pnl = panels_A(1,2);
axes(pnl);
h=annotation('textarrow');
h.Units='centimeters';
h.X = pnl.Position(1) + pnl.Position(3).*[0.8 0.6];
h.Y = pnl.Position(2) + pnl.Position(4).*[1.1 1.1];
h.Text.HorizontalAlignment = 'right';
h.String = 'For map of direction                    '; h.FontSize = 9;

h=annotation('textarrow');
h.Units='centimeters';
h.X = pnl.Position(1) + pnl.Position(3).*[0.7 0.6];
h.Y = pnl.Position(2) + pnl.Position(4).*[0.9 0.9];
h.Text.HorizontalAlignment = 'right';
h.Color = clrs(1,:); h.String = 'Forward           '; h.FontSize = 9;
h.HeadLength = 5; h.HeadWidth = 8;

h=annotation('textarrow');
h.Units='centimeters';
h.X = pnl.Position(1) + pnl.Position(3).*[0.6 0.7];
h.Y = pnl.Position(2) + pnl.Position(4).*[0.73 0.73];
h.Color = clrs(2,:);
h.String = 'Reverse  '; h.FontSize = 9;
h.HeadLength = 5; h.HeadWidth = 8;

%% calc num takeoff/landing replays
% thr = how far (normalized position from the ball is defined takeoff/landing)
% many thresholds
thrs = linspace(0,0.5,100); 
takeoff_counts = zeros(size(thrs));
landing_counts = zeros(size(thrs));
for ii_thr = 1:length(thrs)
    thr = thrs(ii_thr);
    g = classify_replay_landing_takeoff_other(seqs,thr);
    takeoff_counts(ii_thr) = sum(g=='Takeoff');
    landing_counts(ii_thr) = sum(g=='Landing');
end
takeoff_fractions = takeoff_counts ./ length(seqs);
landing_fractions = landing_counts ./ length(seqs);

%% plot (seperate fig) takeoff/landing fractions vs thrs
lw = 2;
figure
subplot(211)
hold on
plot(thrs, landing_fractions,'-g','LineWidth',lw);
plot(thrs, takeoff_fractions,'-m','LineWidth',lw);
xlabel('Distance from the ball (norm.)')
ylabel('Fraction of replays')
legend("Landing","Takeoff",'Location','northwest')
subplot(212)
hold on
plot(thrs, landing_fractions./takeoff_fractions,'-k','LineWidth',lw);
xlabel('Distance from the ball (norm.)')
ylabel('ratio : Landing / Takeoff')
yline(1)

%% now choose one thr
thr = 0.05;
g = classify_replay_landing_takeoff_other(seqs,thr);
nTakeoff = sum(g=='Takeoff');
nLanding = sum(g=='Landing');
nTakeoffLanding = nTakeoff + nLanding;
binom_pval = myBinomTest(nLanding,nTakeoffLanding,0.5,'two');
disp('Panel B:');
disp('Comparing Landing vs. takeoff')
fprintf('Landing: n=%d, Takeoff: n=%d, ratio=%.2g\n',nLanding,nTakeoff, nLanding/nTakeoff);
fprintf('Binomial test (two-sided), pal=%.2g\n',binom_pval);

%% plot bar graph of comparing number of takeoff/landing replays
axes(panels_B(1))
cla
hold on
h=histogram( removecats(g(~(g=='Other')),'Other') ,'Normalization','probability');
h.BarWidth = 0.7;
h.FaceColor = 0.5*[1 1 1];
h.DisplayOrder = 'ascend';
text(1,h.Values(1), "n="+nTakeoff, 'HorizontalAlignment','center','VerticalAlignment','bottom', FontSize=8);
text(2,h.Values(2), "n="+nLanding, 'HorizontalAlignment','center','VerticalAlignment','bottom', FontSize=8);
ylabel('Fraction')
text(-0.4,1.1, 'B', 'Units','normalized','FontWeight','bold');
plot([1 1],[0.5 1],'-k','Clipping','off');
plot([2 2],[0.9 1],'-k','Clipping','off');
plot([1 2],[1 1],'-k','Clipping','off');
% text(1.5,1.05,'***',FontSize=12,HorizontalAlignment='center');
text(1.5,1.11,sprintf('p < 10^{-%d}',ceil(log10(binom_pval))),FontSize=8,HorizontalAlignment='center');

%% MUA firing rate maps over-representations (equal for takeoff vs. landing)
maps = cat(3,MUA_maps.maps);
maps = nanzscore(maps,0,2);
maps_avg = mean(maps,3);
x = linspace(0,1,size(maps,2));
figure
subplot(211)
plot(x,maps_avg')
xlabel('Position in the tunnel (norm)');
ylabel('MUA firing rate (z)')
subplot(212)
y = [maps_avg(1,:);flip(maps_avg(2,:))];
plot(x, mean(y), 'k','LineWidth',2);
% shadedErrorBar(x, mean(y), nansem(y));
xticks([0 1])
xticklabels(["Takeoff","Landing"])
ylabel('MUA firing rate (z)')

%%
pnl = panels_C(1);
axes(pnl);
cla
hold on
plot(x, mean(y), 'k','LineWidth',2);
xticks([0 1])
xticklabels(["Takeoff","Landing"])
% annotation('arrow',0.37+[0 0.28], 0.3975.*[1 1]);
h=annotation('arrow');
h.Units = 'centimeters';
h.X = pnl.Position(1) + pnl.Position(3).*[0.15 0.85];
h.Y = pnl.Position(2) + pnl.Position(4).*[1 1].*-0.13;
xlabel('Flight trajectory')
ylabel({'MUA firing rate during flight';'population average (z)'})
box on 
text(-0.2,1.1, 'C', 'Units','normalized','FontWeight','bold');

%% save fig(s)
if exist('bats_to_include','var')
    bats_str = ['_bats_' char(strjoin(""+bats,'_'))];
else
    bats_str = '_bats_all';
end

fig_name_out = fullfile(res_dir, [fig_name_str bats_str]);
print(fig, fig_name_out, '-dpdf', '-cmyk', '-painters');

disp('figure was successfully saved to pdf/tiff/fig formats');
diary off


%%
function g = classify_replay_landing_takeoff_other(seqs,thr)
%     near_balls_TF = [seqs.start_pos_norm] < thr | [seqs.start_pos_norm] > (1-thr) | ...
%                     [seqs.end_pos_norm] < thr   | [seqs.end_pos_norm] > (1-thr) ;
    
    landing_TF = (~[seqs.forward] & (  ([seqs.start_pos_norm] < thr     & [seqs.state_direction]==-1) | ...
                                       ([seqs.start_pos_norm] > (1-thr) & [seqs.state_direction]== 1) ) )...
                 | ...
                 ( [seqs.forward] & (  ([seqs.end_pos_norm]   < thr     & [seqs.state_direction]==-1) | ...
                                       ([seqs.end_pos_norm]   > (1-thr) & [seqs.state_direction]== 1) ) );
    
    takeoff_TF = (~[seqs.forward] & (  ([seqs.end_pos_norm] < thr     & [seqs.state_direction]== 1) | ...
                                       ([seqs.end_pos_norm] > (1-thr) & [seqs.state_direction]==-1) ) )...
                 | ...
                 ( [seqs.forward] & (  ([seqs.start_pos_norm]   < thr     & [seqs.state_direction]== 1) | ...
                                       ([seqs.start_pos_norm]   > (1-thr) & [seqs.state_direction]==-1) ) );
    g = zeros(size(seqs));
    g(landing_TF) = 1;
    g(takeoff_TF) = 2;
    g = categorical(g,[0 1 2],{'Other','Landing','Takeoff'});
end



