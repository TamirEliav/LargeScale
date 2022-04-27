%% PhD thesis figure 4.7 - replay - over-representations during (rest)

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

%% define output files
res_dir = 'E:\Tamir\PhD\Thesis\resources\ch_4_seq';
mkdir(res_dir)
fig_name_str = 'Fig_4_7_rest_replay_over_representations';
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
panels_A(1,1,1) = axes('position', [2  21 panels_size]);
panels_A(1,2,1) = axes('position', [2  16.5 panels_size]);
panels_A(1,1,2) = axes('position', [12  21 panels_size]);
panels_A(1,2,2) = axes('position', [12  16.5 panels_size]);
panels_A(2,1,1) = axes('position', [2  12 panels_size]);
panels_A(2,2,1) = axes('position', [2  7.5 panels_size]);
panels_A(2,1,2) = axes('position', [12  12 panels_size]);
panels_A(2,2,2) = axes('position', [12  7.5 panels_size]);

panels_B(1) = axes('position', [2  2.5 3 3]);
panels_C(1) = axes('position', [7  2.5 3 3]);
panels_D(1) = axes('position', [12 2.5 3 3]);

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
for ii_exp = 1:height(T)
    % load exp data
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'details','rest');
    epoch_type = 'rest';
    params_opt = 11;
    [events_session, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
    % assign resting ball num
    X = [exp.rest.events.start_ts];
    V = [exp.rest.events.ball_num];
    Xq = [events_session.peak_ts];
    [events_session.ball_num] = disperse( interp1(X, V, Xq, 'previous','extrap') );
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

%% classify seq as forward / reverse 
gForRev = categorical([seqs.forward],[true false],{'Forward','Reverse'})';
nFor = sum(gForRev=='Forward');
nRev = sum(gForRev=='Reverse');
ForRev_binom_pval = myBinomTest(nFor,length(gForRev),0.5,'two');
disp('Comparing (immediate) Forward vs. Reverse')
fprintf('Forward: n=%d, Reverse: n=%d, ratio=%.2g\n',nFor,nRev, nFor/nRev);
fprintf('Binomial test (two-sided), pal=%.2g\n', ForRev_binom_pval);

%% classify seq as Takeoff / Landing (or non of them)
thr = 0.05;
gTakeLand = classify_replay_landing_takeoff_other(seqs,thr);
nTakeoff = sum(gTakeLand=='Takeoff');
nLanding = sum(gTakeLand=='Landing');
nTakeoffLanding = nTakeoff + nLanding;
TakeLand_binom_pval = myBinomTest(nLanding,nTakeoffLanding,0.5,'two');
disp('Comparing Landing vs. takeoff')
fprintf('Landing: n=%d, Takeoff: n=%d, ratio=%.2g\n',nLanding,nTakeoff, nLanding/nTakeoff);
fprintf('Binomial test (two-sided), pal=%.2g\n',TakeLand_binom_pval);

%% classify seq as (immediate) past / future
gPastFuture = categorical( ([events.ball_num] == 1 & [seqs.state_direction] == 1) | ...
                           ([events.ball_num] == 2 & [seqs.state_direction] == -1), ...
    [true false],{'Future','Past'})';
nPast = sum(gPastFuture =='Past');
nFuture = sum(gPastFuture =='Future');
PastFuture_binom_pval = myBinomTest(nPast,length(gPastFuture),0.5,'two');
disp('Comparing (immediate) Past vs. Future')
fprintf('Past: n=%d, Future: n=%d, ratio=%.2g\n',nPast,nFuture, nPast/nFuture);
fprintf('Binomial test (two-sided), pal=%.2g\n',PastFuture_binom_pval);

%%
tbl_seq_cls = table(gForRev,gTakeLand,gPastFuture);
fig2=figure;
fig2.WindowState = 'maximized';
tiledlayout(3,3)
% nexttile
% heatmap(tbl_seq_cls,'gTakeLand','gForRev');
% nexttile
% heatmap(tbl_seq_cls,'gForRev','gPastFuture')
% nexttile
% heatmap(tbl_seq_cls,'gTakeLand','gPastFuture')
nexttile
histogram(gForRev); ylabel('Counts')
nexttile
histogram(gTakeLand); ylabel('Counts')
nexttile
histogram(gPastFuture); ylabel('Counts')
nexttile
histogram(gTakeLand.*gForRev,'DisplayOrder','descend'); ylabel('Counts')
nexttile
histogram(gTakeLand.*gPastFuture,'DisplayOrder','descend'); ylabel('Counts')
nexttile
histogram(gForRev.*gPastFuture,'DisplayOrder','descend'); ylabel('Counts')
nexttile
histogram(gTakeLand.*gForRev.*gPastFuture,'DisplayOrder','descend'); ylabel('Counts')
sgtitle('Rest replay classifications');

%% plot replay start/end position histograms (per forward/reverse replays)
nbins = 100;
cmap = brewermap(100,'RdBu');
clrs = cmap([0.9 0.1].*size(cmap,1),:);
directions = [1 -1];
rest_balls = [1 2];
for ii_ball = 1:2
    for ii_dir = 1:length(directions)
        seqs_dir_TF = [seqs.state_direction]==directions(ii_dir);
        seqs_ball_TF = [events.ball_num]==ii_ball;
    
        axes(panels_A(ii_ball,1,ii_dir))
        cla
        hold on
        x = [seqs.start_pos_norm];
        histogram(x(gForRev=='Forward' & seqs_dir_TF' & seqs_ball_TF'),linspace(0,1,nbins+1),'Normalization','count',hist_color_prop,clrs(1,:),'DisplayStyle',hist_style);
        histogram(x(gForRev=='Reverse' & seqs_dir_TF' & seqs_ball_TF'),linspace(0,1,nbins+1),'Normalization','count',hist_color_prop,clrs(2,:),'DisplayStyle',hist_style);
        xlabel('Start position of replay (norm.)')
        ylabel('Counts')
        hax=gca;
        hax.XRuler.TickLabelGapOffset = -1;
    
        axes(panels_A(ii_ball,2,ii_dir))
        cla
        hold on
        x = [seqs.end_pos_norm];
        histogram(x(gForRev=='Forward' & seqs_dir_TF' & seqs_ball_TF'),linspace(0,1,nbins+1),'Normalization','count',hist_color_prop,clrs(1,:),'DisplayStyle',hist_style);
        histogram(x(gForRev=='Reverse' & seqs_dir_TF' & seqs_ball_TF'),linspace(0,1,nbins+1),'Normalization','count',hist_color_prop,clrs(2,:),'DisplayStyle',hist_style);
        xlabel('End position of replay (norm.)')
        ylabel('Counts')
        hax=gca;
        hax.XRuler.TickLabelGapOffset = -1;
    end
end
% set ylimits
% panels_A(1,1).YLim = [0 80];
% panels_A(2,1).YLim = [0 150];
% panels_A(1,2).YLim = [0 150];
% panels_A(2,2).YLim = [0 220];

% linkaxes(panels_A(:),'y')

axes(panels_A(1,1))
text(-0.15,1.15, 'A', 'Units','normalized','FontWeight','bold');

for ii=1:numel(panels_A)
    pnl = panels_A(ii); axes(pnl);
    text(0,-0.27*range(pnl.YLim),{'Resting';'ball 1'},'HorizontalAlignment','center','FontSize',8);
    text(1,-0.27*range(pnl.YLim),{'Resting';'ball 2'},'HorizontalAlignment','center','FontSize',8);
end

text(panels_A(1,1,2),1.1,-0.5,'Replay while resting on ball 1','Units','normalized','Rotation',-90,'HorizontalAlignment','center','FontWeight','bold')
text(panels_A(2,1,2),1.1,-0.5,'Replay while resting on ball 2','Units','normalized','Rotation',-90,'HorizontalAlignment','center','FontWeight','bold')

%% legend (dir 1)
pnl = panels_A(1,1,1);
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
pnl = panels_A(1,1,2);
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
    takeoff_counts(ii_thr) = sum(g =='Takeoff');
    landing_counts(ii_thr) = sum(g =='Landing');
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

%% comparing takeoff vs landing (bar plot)
axes(panels_B(1))
cla
hold on
gTakeLand2 = gTakeLand;
h=histogram( removecats(gTakeLand2(~(gTakeLand2=='Other')),'Other') ,'Normalization','probability');
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
text(1.5,1.11,sprintf('p < 10^{-%d}',ceil(log10(TakeLand_binom_pval))),FontSize=8,HorizontalAlignment='center');

%% comparing Forward vs Reverse (bar plot)
axes(panels_C(1))
cla
hold on
h=histogram( gForRev ,'Normalization','probability');
h.BarWidth = 0.7;
h.FaceColor = 0.5*[1 1 1];
h.DisplayOrder = 'ascend';
for ii_cat = 1:length(h.Categories)
    text(ii_cat,h.Values(ii_cat), "n="+sum(gForRev ==h.Categories{ii_cat}), 'HorizontalAlignment','center','VerticalAlignment','bottom', FontSize=8);
end
ylabel('Fraction')
text(-0.4,1.1, 'C', 'Units','normalized','FontWeight','bold');
plot([1 1],[0.9 1],'-k','Clipping','off');
plot([2 2],[0.9 1],'-k','Clipping','off');
plot([1 2],[1 1],'-k','Clipping','off');
text(1.5,1.11,sprintf('p < 10^{-%d}',ceil(log10(ForRev_binom_pval))),FontSize=8,HorizontalAlignment='center');

%% comparing Past vs Future (bar plot)
axes(panels_D(1))
cla
hold on
h=histogram( gPastFuture ,'Normalization','probability');
h.BarWidth = 0.7;
h.FaceColor = 0.5*[1 1 1];
h.DisplayOrder = 'ascend';
for ii_cat = 1:length(h.Categories)
    text(ii_cat,h.Values(ii_cat), "n="+sum(gPastFuture==h.Categories{ii_cat}), 'HorizontalAlignment','center','VerticalAlignment','bottom', FontSize=8);
end
ylabel('Fraction')
text(-0.4,1.1, 'D', 'Units','normalized','FontWeight','bold');
plot([1 1],[0.9 1],'-k','Clipping','off');
plot([2 2],[0.9 1],'-k','Clipping','off');
plot([1 2],[1 1],'-k','Clipping','off');
text(1.5,1.11,sprintf('p < 10^{-%d}',ceil(log10(PastFuture_binom_pval))),FontSize=8,HorizontalAlignment='center');

%% save fig(s)
if exist('bats_to_include','var')
    bats_str = ['_bats_' char(strjoin(""+bats,'_'))];
else
    bats_str = '_bats_all';
end

fig_name_out = fullfile(res_dir, [fig_name_str bats_str]);
print(fig, fig_name_out, '-dpdf', '-cmyk', '-painters');

fig_name_out = fullfile(res_dir, [fig_name_str '_seq_classifications_' bats_str]);
saveas(fig2, fig_name_out, 'jpg');

disp('figure was successfully saved to pdf/tiff/fig formats');
diary off


%%




