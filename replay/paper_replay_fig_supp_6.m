%% Replay - Fig supp 8 - reverse/forward + future+past + takeoff/landing/midair + ball1/ball2
%%
clear 
clc
close all

%% data options 
params_opt = 11; % decoding opt 
% params_opt = 21; % decoding opt (random walk = fixed speed)

%% plotting options

%% graphics params
TakeLandCategoriesOrder = {
    'Takeoff'
    'Mid-air'
    'Landing'};
ForRevCategoriesOrder = {
    'Forward'
    'Reverse'};
PastFutureCategoriesOrder = {
    'Past'
    'Future'};

hist_style = 'stairs';
% hist_style = 'bar';
switch hist_style 
    case 'stairs'
        hist_color_prop = 'EdgeColor';
    case 'bar'
        hist_color_prop = 'FaceColor';
end


%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Extended_Data_Fig_8';
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
panels{1}(1) = axes('position', [2  19 4 3]);
panels{2}(1) = axes('position', [8  19 3 3]);
panels{3}(1) = axes('position', [13 19 3 3]);
panels{4}(1) = axes('position', [2 14.5 4 3]);
panels{4}(2) = axes('position', [2 14.5 .5 2]+[.3 1.4 0 0]);
panels{5}(1) = axes('position', [8 14.5 4 3]);

w = 8.5;
h = 3.5;
x = [2 12];
y = [8 2.2];
panels{6}(1,1) = axes('position', [x(1) y(1) w h]);
panels{6}(1,2) = axes('position', [x(2) y(1) w h]);
panels{6}(2,1) = axes('position', [x(1) y(2) w h]);
panels{6}(2,2) = axes('position', [x(2) y(2) w h]);

sdf = cellfun(@(x)x(:)', panels, 'UniformOutput', false);
sdf = [sdf{:}];
total_offset = [0 2.5];
for ii = 1:length(sdf)
    sdf(ii).Position([1 2]) = sdf(ii).Position([1 2]) + total_offset;
end

%% load data
if ~exist('events','var')

[exp_list,T] = decoding_get_inclusion_list();
T = T(exp_list,:);
clear exp_list
groupsummary(T,'bat_num')
if exist('bats_to_include','var')
    T = groupfilter(T,"bat_num",@(x)ismember(x,bats_to_include),'bat_num');
end
bats = unique(T.bat_num)

events = {};
for ii_exp = 1:height(T)
    % load exp data
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'details','rest');
    epoch_type = 'rest';
    [events_session, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
    [~, TF] = decoding_apply_seq_inclusion_criteria([events_session.seq_model]);
    events_session(~TF) = [];
    if isempty(events_session)
        continue;
    end
    % assign resting ball num
    X = [exp.rest.events.start_ts];
    V = [exp.rest.events.ball_num];
    Xq = [events_session.peak_ts];
    [events_session.ball_num] = disperse( interp1(X, V, Xq, 'previous','extrap') );
    events{ii_exp} = events_session;
end
T.nEvents = cellfun(@length,events)';
sortrows( groupsummary(T,'bat_num',["median","mean","max","sum"],"nEvents"),"sum_nEvents", 'descend')
events = [events{:}];
seqs = [events.seq_model];
end

%% classify seq as forward / reverse 
gForRev = categorical([seqs.forward],[true false],{'Forward','Reverse'})';
gForRev = categorical(gForRev,ForRevCategoriesOrder);
nFor = sum(gForRev=='Forward');
nRev = sum(gForRev=='Reverse');
ForRev_binom_pval = myBinomTest(nFor,length(gForRev),0.5,'two');
disp('Comparing Forward vs. Reverse')
fprintf('Forward: n=%d, Reverse: n=%d, ratio=%.2g\n',nFor,nRev, nFor/nRev);
fprintf('Binomial test (two-sided), p=%.2g\n', ForRev_binom_pval);

%% classify seq as Takeoff / Landing (or non of them)
thr = 0.05;
gTakeLand = classify_replay_landing_takeoff_other(seqs,thr);
gTakeLand = categorical(gTakeLand,TakeLandCategoriesOrder);
nTakeoff = sum(gTakeLand=='Takeoff');
nLanding = sum(gTakeLand=='Landing');
nTakeoffLanding = nTakeoff + nLanding;
TakeLand_binom_pval = myBinomTest(nLanding,nTakeoffLanding,0.5,'two');
disp('Comparing Landing vs. takeoff')
fprintf('Landing: n=%d, Takeoff: n=%d, ratio=%.2g\n',nLanding,nTakeoff, nLanding/nTakeoff);
fprintf('Binomial test (two-sided), p=%.2g\n',TakeLand_binom_pval);

%% classify seq as (immediate) past / future
gPastFuture = categorical( ([events.ball_num] == 1 & [seqs.state_direction] == 1) | ...
                           ([events.ball_num] == 2 & [seqs.state_direction] == -1), ...
    [true false],{'Future','Past'})';
gPastFuture = categorical(gPastFuture,PastFutureCategoriesOrder);
nPast = sum(gPastFuture =='Past');
nFuture = sum(gPastFuture =='Future');
PastFuture_binom_pval = myBinomTest(nPast,length(gPastFuture),0.5,'two');
disp('Comparing (immediate) Past vs. Future')
fprintf('Past: n=%d, Future: n=%d, ratio=%.2g\n',nPast,nFuture, nPast/nFuture);
fprintf('Binomial test (two-sided), p=%.2g\n',PastFuture_binom_pval);


%% takeoff vs landing
axes(panels{1}(1))
cla
hold on
G = categorical(categories(gTakeLand));
gcounts = tabulate(gTakeLand);
counts = [gcounts{:,2}];
proportions = counts./sum(counts);
expectedProportion = [.05 .9 0.05];
p = myBinomTest(counts,sum(counts),expectedProportion,'two').*length(counts);
p(p==0)=realmin;
fprintf('\nboth balls: takeoff vs landing\n')
for ii = 1:length(p)
    fprintf('%s: P binom (two-sided) = %d\n',TakeLandCategoriesOrder{ii},p(ii));
end

xG = 1:length(G);
Y = proportions;
h=bar(xG,Y);
h.BarWidth = 0.75;
h.FaceColor = 0.5*[1 1 1];
hax=gca;
% hax.XTickLabel = TakeLandCategoriesOrder;
for ii = 1:length(G)
    x = ii+[-1 1].*h.BarWidth/2;
    y = expectedProportion(ii).*sum(Y);
    line(x, [y, y], 'Color', 'r', 'LineWidth', 2);
    x = ii;
    y = Y(ii)+0.001*max(Y);
%     text(x,y, "{\itP} < 10^{"+ceil(log10(p(ii)))+"}", horizontalAlignment = 'center',VerticalAlignment='bottom',FontSize=6);
    signifStr = genSignifStrAstricks(p(ii));
    str = {signifStr;"{\itn} = " + counts(ii)};
    text(x,y, str, horizontalAlignment = 'center',VerticalAlignment='bottom',FontSize=8);
end
% xlabel('Replay speed (m/s)', 'Units','normalized', 'Position',[0.5 -0.12]);
ylabel('Fraction', 'Units','normalized', 'Position',[-0.18 .5]);
hax=gca;
hax.TickLength(1) = [0.025];
hax.XRuler.TickLabelGapOffset = -1;
hax.YRuler.TickLabelGapOffset = 0;
hax.XLim = [min(xG) max(xG)] + [-1 1].*0.5;
hax.XTick = xG;
hax.XTickLabel = [];
for ii=1:length(TakeLandCategoriesOrder)
    text(ii, hax.YLim(1)-0.15*range(hax.YLim), {TakeLandCategoriesOrder{ii},'zone'},'HorizontalAlignment','center','FontSize',7)
end

%% takeoff vs landing - separated by balls
% for ii_ball = 1:2
% axes(panels{4}(ii_ball))
% cla
% hold on
% 
% events_ball_IX = [events.ball_num] == ii_ball;
% gTakeLand_by_ball = gTakeLand(events_ball_IX);
% G = categorical(categories(gTakeLand_by_ball));
% gcounts = tabulate(gTakeLand_by_ball);
% counts = [gcounts{:,2}];
% proportions = counts./sum(counts);
% expectedProportion = [.05 .9 0.05];
% expectedCounts = expectedProportion .* sum(counts);
% p = myBinomTest(counts,sum(counts),expectedProportion,'two').*length(counts);
% p(p==0)=realmin;
% xG = 1:length(G);
% 
% Y = proportions;
% h=bar(xG,Y);
% h.BarWidth = 0.75;
% h.FaceColor = 0.5*[1 1 1];
% hax=gca;
% hax.XTick = xG;
% hax.XTickLabel = TakeLandCategoriesOrder;
% for ii = 1:length(G)
%     x = ii+[-1 1].*h.BarWidth/2;
%     y = expectedProportion(ii).*sum(Y);
%     line(x, [y, y], 'Color', 'r', 'LineWidth', 2);
%     x = ii;
%     y = Y(ii)+0.001*max(Y);
% %     text(x,y, "{\itP} < 10^{"+ceil(log10(p(ii)))+"}", horizontalAlignment = 'center',VerticalAlignment='bottom',FontSize=6);
%     signifStr = genSignifStrAstricks(p(ii));
%     str = {signifStr;"{\itn} = " + counts(ii)};
%     text(x,y, str, horizontalAlignment = 'center',VerticalAlignment='bottom',FontSize=8);
% end
% % xlabel('Replay speed (m/s)', 'Units','normalized', 'Position',[0.5 -0.12]);
% ylabel('Fraction', 'Units','normalized', 'Position',[-0.2 .5]);
% hax=gca;
% hax.TickLength(1) = [0.025];
% hax.XRuler.TickLabelGapOffset = -1;
% hax.YRuler.TickLabelGapOffset = -1;
% hax.XLim = [min(xG) max(xG)] + [-1 1].*0.5;
% 
% end


%% Forward vs Reverse
axes(panels{3}(1))
cla
hold on
G = categorical(categories(gForRev));
gcounts = tabulate(gForRev);
counts = [gcounts{:,2}];
proportions = counts./sum(counts);
expectedProportion = [.5 .5];
expectedCounts = expectedProportion .* sum(counts);
p = myBinomTest(counts,sum(counts),expectedProportion,'two').*length(counts);
p(p==0)=realmin;
fprintf('\nboth balls: Forward vs Reverse\n')
for ii = 1:length(p)
    fprintf('%s: P binom (two-sided) = %d\n',ForRevCategoriesOrder{ii},p(ii));
end

xG = 1:length(G);
Y = proportions;
h=bar(xG,Y);
h.BarWidth = 0.75;
h.FaceColor = 0.5*[1 1 1];
hax=gca;
hax.XTick = xG;
hax.XTickLabel = ForRevCategoriesOrder;
for ii = 1:length(G)
    x = ii+[-1 1].*h.BarWidth/2;
    y = expectedProportion(ii).*sum(Y);
    line(x, [y, y], 'Color', 'r', 'LineWidth', 2);
    x = ii;
    y = Y(ii)+0.03*max(Y);
%     text(x,y, "{\itP} < 10^{"+ceil(log10(p(ii)))+"}", horizontalAlignment = 'center',VerticalAlignment='bottom',FontSize=6);
    signifStr = genSignifStrAstricks(p(ii));
    str = {signifStr;"{\itn} = " + counts(ii)};
    text(x,y, str, horizontalAlignment = 'center',VerticalAlignment='bottom',FontSize=8);
end
% xlabel('Replay speed (m/s)', 'Units','normalized', 'Position',[0.5 -0.12]);
ylabel('Fraction', 'Units','normalized', 'Position',[-0.2 .5]);
hax=gca;
hax.TickLength(1) = [0.025];
hax.XRuler.TickLabelGapOffset = -1;
hax.YRuler.TickLabelGapOffset = 0;
hax.XLim = [min(xG) max(xG)] + [-1 1].*0.5;
hax.YLim = [0 0.7];


%% Forward vs Reverse - separated by balls
% for ii_ball = 1:2
% axes(panels{5}(ii_ball))
% cla
% hold on
% 
% events_ball_IX = [events.ball_num] == ii_ball;
% gForRev_by_ball = gForRev(events_ball_IX);
% G = categorical(categories(gForRev_by_ball));
% gcounts = tabulate(gForRev_by_ball);
% counts = [gcounts{:,2}];
% proportions = counts./sum(counts);
% expectedProportion = [.5 .5];
% expectedCounts = expectedProportion .* sum(counts);
% p = myBinomTest(counts,sum(counts),expectedProportion,'two').*length(counts);
% p(p==0)=realmin;
% xG = 1:length(G);
% 
% Y = proportions;
% h=bar(xG,Y);
% h.BarWidth = 0.75;
% h.FaceColor = 0.5*[1 1 1];
% hax=gca;
% hax.XTick = xG;
% hax.XTickLabel = ForRevCategoriesOrder;
% for ii = 1:length(G)
%     x = ii+[-1 1].*h.BarWidth/2;
%     y = expectedProportion(ii).*sum(Y);
%     line(x, [y, y], 'Color', 'r', 'LineWidth', 2);
%     x = ii;
%     y = Y(ii)+0.001*max(Y);
% %     text(x,y, "{\itP} < 10^{"+ceil(log10(p(ii)))+"}", horizontalAlignment = 'center',VerticalAlignment='bottom',FontSize=6);
%     signifStr = genSignifStrAstricks(p(ii));
%     str = {signifStr;"{\itn} = " + counts(ii)};
%     text(x,y, str, horizontalAlignment = 'center',VerticalAlignment='bottom',FontSize=8);
% end
% % xlabel('Replay speed (m/s)', 'Units','normalized', 'Position',[0.5 -0.12]);
% ylabel('Fraction', 'Units','normalized', 'Position',[-0.2 .5]);
% hax=gca;
% hax.TickLength(1) = [0.025];
% hax.XRuler.TickLabelGapOffset = -1;
% hax.YRuler.TickLabelGapOffset = -1;
% hax.XLim = [min(xG) max(xG)] + [-1 1].*0.5;
% hax.YLim = [0 0.7];
% 
% end



%% Future vs Past
axes(panels{2}(1))
cla
hold on
G = categorical(categories(gPastFuture));
gcounts = tabulate(gPastFuture);
counts = [gcounts{:,2}];
proportions = counts./sum(counts);
expectedProportion = [.5 .5];
expectedCounts = expectedProportion .* sum(counts);
p = myBinomTest(counts,sum(counts),expectedProportion,'two').*length(counts);
p(p==0)=realmin;
fprintf('\nboth balls: Future vs Past\n')
for ii = 1:length(p)
    fprintf('%s: P binom (two-sided) = %d\n',PastFutureCategoriesOrder{ii},p(ii));
end

xG = 1:length(G);
Y = proportions;
h=bar(xG,Y);
h.BarWidth = 0.75;
h.FaceColor = 0.5*[1 1 1];
hax=gca;
hax.XTick = xG;
hax.XTickLabel = PastFutureCategoriesOrder;
for ii = 1:length(G)
    x = ii+[-1 1].*h.BarWidth/2;
    y = expectedProportion(ii).*sum(Y);
    line(x, [y, y], 'Color', 'r', 'LineWidth', 2);
    x = ii;
    y = Y(ii)+0.1*max(Y);
%     text(x,y, "{\itP} < 10^{"+ceil(log10(p(ii)))+"}", horizontalAlignment = 'center',VerticalAlignment='bottom',FontSize=6);
    signifStr = genSignifStrAstricks(p(ii));
    str = {signifStr;"{\itn} = " + counts(ii)};
    text(x,y, str, horizontalAlignment = 'center',VerticalAlignment='bottom',FontSize=8);
end
% xlabel('Replay speed (m/s)', 'Units','normalized', 'Position',[0.5 -0.12]);
ylabel('Fraction', 'Units','normalized', 'Position',[-0.2 .5]);
hax=gca;
hax.TickLength(1) = [0.025];
hax.XRuler.TickLabelGapOffset = -1;
hax.YRuler.TickLabelGapOffset = 0;
hax.XLim = [min(xG) max(xG)] + [-1 1].*0.5;
hax.YLim = [0 0.7];

%% Future vs Past - separated by balls
% for ii_ball = 1:2
% axes(panels{6}(ii_ball))
% cla
% hold on
% 
% events_ball_IX = [events.ball_num] == ii_ball;
% gPastFuture_by_ball = gPastFuture(events_ball_IX);
% G = categorical(categories(gPastFuture_by_ball));
% gcounts = tabulate(gPastFuture_by_ball);
% counts = [gcounts{:,2}];
% proportions = counts./sum(counts);
% expectedProportion = [.5 .5];
% expectedCounts = expectedProportion .* sum(counts);
% p = myBinomTest(counts,sum(counts),expectedProportion,'two').*length(counts);
% p(p==0)=realmin;
% xG = 1:length(G);
% 
% Y = proportions;
% h=bar(xG,Y);
% h.BarWidth = 0.75;
% h.FaceColor = 0.5*[1 1 1];
% hax=gca;
% hax.XTick = xG;
% hax.XTickLabel = PastFutureCategoriesOrder;
% for ii = 1:length(G)
%     x = ii+[-1 1].*h.BarWidth/2;
%     y = expectedProportion(ii).*sum(Y);
%     line(x, [y, y], 'Color', 'r', 'LineWidth', 2);
%     x = ii;
%     y = Y(ii)+0.001*max(Y);
% %     text(x,y, "{\itP} < 10^{"+ceil(log10(p(ii)))+"}", horizontalAlignment = 'center',VerticalAlignment='bottom',FontSize=6);
%     signifStr = genSignifStrAstricks(p(ii));
%     str = {signifStr;"{\itn} = " + counts(ii)};
%     text(x,y, str, horizontalAlignment = 'center',VerticalAlignment='bottom',FontSize=8);
% end
% % xlabel('Replay speed (m/s)', 'Units','normalized', 'Position',[0.5 -0.12]);
% ylabel('Fraction', 'Units','normalized', 'Position',[-0.2 .5]);
% % text(1.6,0.5,"Ball "+ii_ball,'Units','normalized','FontSize',10,'HorizontalAlignment','center');
% hax=gca;
% hax.TickLength(1) = [0.025];
% hax.XRuler.TickLabelGapOffset = -1;
% hax.YRuler.TickLabelGapOffset = -1;
% hax.XLim = [min(xG) max(xG)] + [-1 1].*0.5;
% hax.YLim = [0 0.7];
% 
% end









%% Takeoff/Landing X Past/Future X Forward/Reverse
axes(panels{4}(1))
cla reset
hold on

% from fig 2  (start) ---------------
TakeLand_thr = 0.05;
g = gTakeLand .* gPastFuture;
[N,G] = histcounts(g);
N = reshape(N,2,3)';
G = reshape(G,2,3)';
takeoff_IX = find(all(contains(G,'Takeoff'),2));
midair_IX = find(all(contains(G,'Mid-air'),2));
landing_IX = find(all(contains(G,'Landing'),2));
new_order = [takeoff_IX midair_IX landing_IX];
N = N(new_order,:);
G = G(new_order,:);
N([1 3],:) = N([1 3],:)./TakeLand_thr;
N([2],:) = N([2],:)./(1-2*TakeLand_thr);
N = N ./ sum(N,"all");
% from fig 2  (end) ---------------

g = gTakeLand .* gPastFuture .*gForRev;
[N,G] = histcounts(g);
N = reshape(N,2,2,3);
G = reshape(G,2,2,3);
N = permute(N,[3 2 1]);
G = permute(G,[3 2 1]);
new_order = [takeoff_IX midair_IX landing_IX];  
N = N(new_order,:,:);
G = G(new_order,:,:);

N([1 3],:,:) = N([1 3],:,:)./TakeLand_thr;
N([2],  :,:) = N([2],  :,:)./(1-2*TakeLand_thr);
N = N ./ sum(N,"all");

w = 0.28;
x = linspace(0,1,3);
hs = (x(2)-x(1))*0.35;
clear hb1 hb2
hb1 = bar(x-hs/2,   squeeze(N(:,1,:)), w, 'stacked');
hb2 = bar(x+hs/2,   squeeze(N(:,2,:)), w, 'stacked');
hb1(1).FaceColor = 'flat';
hb1(2).FaceColor = 'flat';
hb2(1).FaceColor = 'flat';
hb2(2).FaceColor = 'flat';
hb1(1).CData = 0.00*[1 1 1];
hb1(2).CData = 0.35*[1 1 1];
hb2(1).CData = 0.85*[1 1 1];
hb2(2).CData = 1.00*[1 1 1];
% hb1(1).FaceColor = 0.0*[1 1 1];
% hb1(2).FaceColor = 0.0*[1 1 1];
% hb2(1).FaceColor = 1.0*[1 1 1];
% hb2(2).FaceColor = 1.0*[1 1 1];
% plot(hb1(1).XData+w/4.*[-1 1]',hb1(1).YData([1 1;2 2;3 3]'),'w','LineWidth',1.1);
% plot(hb2(1).XData+w/4.*[-1 1]',hb2(1).YData([1 1;2 2;3 3]'),'k','LineWidth',1.1);
% hb1(1).EdgeColor = 'w';
% hb1(2).EdgeColor = 'w';
% hb2(1).EdgeColor = 'k';
% hb2(2).EdgeColor = 'k';
xlim([0 1]+w.*[-1 1])
xticklabels()
hax=gca;
hax.XTickLabel = [];
for ii=1:length(TakeLandCategoriesOrder)
    text(x(ii), hax.YLim(1)-0.15*range(hax.YLim), {TakeLandCategoriesOrder{ii},'zone'},'HorizontalAlignment','center','FontSize',7)
end

ylabel('Fraction')

%% add legend 
clc
axes(panels{4}(2))
cla reset
hold on
axis off
axis ij

V = [.2 0; .2 .7; .8 .7; .8 0];
F = [1 2 3 4 1];
h=patch('Faces',F,'Vertices',V+[0 1]*1); h.FaceColor=[1 1 1]*.1;
h=patch('Faces',F,'Vertices',V+[0 1]*2); h.FaceColor=[1 1 1]*.3;
h=patch('Faces',F,'Vertices',V+[0 1]*3); h.FaceColor=[1 1 1]*.7;
h=patch('Faces',F,'Vertices',V+[0 1]*4); h.FaceColor=[1 1 1]*.9;

% x = [0.2 0.2 .8 .8];
% y = [0 .7 .7 0];
% h=patch('XData',x, 'YData', y+1, 'FaceColor','Flat')
% patch(x,y+1, 'r', 'FaceColor','Flat', 'CData',0.35*[1 1 1])
% patch(x,y+2, 'r', 'FaceColor','Flat', 'CData',0.35*[1 1 1])
% patch(x,y+3, 'r', 'FaceColor','Flat', 'CData',0.35*[1 1 1])
% patch(x,y+4, 'r', 'FaceColor','Flat', 'CData',0.35*[1 1 1])
% patch(x,y+2, 0.*[1 1 1])
% patch(x,y+3, 0.85*[1 1 1])
% patch(x,y+4, 1.0*[1 1 1])
x = 1.1;
y = [1:4]+0.35;
font_size = 7;
text(x,y(1),'Past forward','FontSize',font_size,'HorizontalAlignment','left','VerticalAlignment','middle')
text(x,y(2),'Past reverse','FontSize',font_size,'HorizontalAlignment','left','VerticalAlignment','middle')
text(x,y(3),'Future forward','FontSize',font_size,'HorizontalAlignment','left','VerticalAlignment','middle')
text(x,y(4),'Future reverse','FontSize',font_size,'HorizontalAlignment','left','VerticalAlignment','middle')
xlim([0 1])
ylim([0 5])

%% replay distance from current location
axes(panels{5}(1))
cla reset
hold on
x = abs([seqs.middle_pos_norm] - interp1([1 2],[0 1],[events.rest_ball_num]));
histogram(x,'BinWidth',0.05,'BinLimits',[0 1],'DisplayStyle','stairs','Normalization','probability','EdgeColor','k');
% gTakeLand = reordercats(gTakeLand,["Takeoff","Mid-air","Landing"]);
% G = categories(gTakeLand);
% for ii_g = 1:length(G)
%     gcurr = G{ii_g};
%     IX = gTakeLand==gcurr;
% %     IX = IX & gForRev=="Forward";
% %     IX = IX & gForRev=="Reverse";
%     seqs_g = seqs(IX);
%     events_g = events(IX);
%     x = abs([seqs_g.middle_pos_norm] - interp1([1 2],[0 1],[events_g.rest_ball_num]));
%     histogram(x,'BinWidth',0.05,'BinLimits',[0 1],'DisplayStyle','stairs','DisplayName',gcurr,'Normalization','probability');
% end
xlim([0 1])
xticks([0 1])
% legend('Box','off')
xlabel({'Distance of awake replay';'from current position (norm.)'},'units','normalized','position',[0.5 -.05])
ylabel('Probability')








%% ========================================================================

%% plot replay start/end position histograms (per forward/reverse replays)
nbins = 100;
cmap = brewermap(100,'RdBu');
clrs = cmap([0.9 0.1].*size(cmap,1),:);
directions = [1 -1];
g = gForRev';
for ii_dir = 1:length(directions)
    seqs_dir_TF = [seqs.state_direction]==directions(ii_dir);

    axes(panels{6}(1,ii_dir))
    cla reset
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

    axes(panels{6}(2,ii_dir))
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
panels{6}(1,1).YLim = [0 80];
panels{6}(2,1).YLim = [0 180];
panels{6}(1,2).YLim = [0 120];
panels{6}(2,2).YLim = [0 150];

for ii=1:numel(panels{6})
    pnl = panels{6}(ii); axes(pnl);
    text(0,-0.24*range(pnl.YLim),{'Resting';'platform 1'},'HorizontalAlignment','center','FontSize',8);
    text(1,-0.24*range(pnl.YLim),{'Resting';'platform 2'},'HorizontalAlignment','center','FontSize',8);
end

% pointers to takeoff/landing in the histogram (dir 1)
TL_fnt_sz = 7;
h=annotation('textarrow');
h.Parent=panels{6}(1,1);
h.X = 0.011*[1 1];
h.Y = 45+[0 -7];
h.String = '    Takeoffs';
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels{6}(1,1);
h.X = 0.052*[1 1];
h.Y = 22+[0 -7];
h.String = {'           Reverse','           takeoffs'};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels{6}(1,1);
h.X = 0.92.*[1 1];
h.Y = 75+[0 -7];
h.String = {'Landings          '};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels{6}(1,1);
h.X = 0.985*[1 1];
h.Y = 70+[0 -7];
h.String = {'            Reverse';'            landings'};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels{6}(2,1);
h.X = 0.005*[1 1];
h.Y = 61+[0 -30];
h.String = {'     Reverse','     takeoffs'};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels{6}(2,1);
h.X = 0.05*[1 1];
h.Y = 35+[0 -15];
h.String = {'        Takeoffs'};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels{6}(2,1);
h.X = 0.985*[1 1];
h.Y = 180+[0 -12];
h.String = 'Landings';
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels{6}(2,1);
h.X = 0.93*[1 1];
h.Y = 60+[0 -12];
h.String = {'Reverse          ','landings          '};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

% pointers to takeoff/landing in the histogram (dir 2)
h=annotation('textarrow');
h.Parent=panels{6}(1,2);
h.X = 0.015*[1 1];
h.Y = 110+[0 -12];
h.String = {'     Reverse','     landings'};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels{6}(1,2);
h.X = 0.06*[1 1];
h.Y = 65+[0 -12];
h.String = '          Landings';
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels{6}(1,2);
h.X = 0.99*[1 1];
h.Y = 70+[0 -20];
h.String = '       Takeoffs';
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels{6}(1,2);
h.Units='normalized';
h.X = 0.95+[0 0];
h.Y = 47+[0 -15];
h.String = {'Reverse          ','takeoffs          '};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels{6}(2,2);
h.X = 0.015*[1 1];
h.Y = 140+[0 -15];
h.String = '    Landings';
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels{6}(2,2);
h.X = 0.07*[1 1];
h.Y = 65+[0 -15];
h.String = {'       Reverse','       landings'};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels{6}(2,2);
h.X = 0.99*[1 1];
h.Y = 60+[0 -20];
h.String = {'       Reverse','       takeoffs'};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

h=annotation('textarrow');
h.Parent=panels{6}(2,2);
h.Units='normalized';
h.X = 0.93+[0 0];
h.Y = 50+[0 -15];
h.String = {'Takeoffs         '};
h.FontSize = TL_fnt_sz; h.HeadLength = 4; h.HeadWidth = 4;

%% legend (dir 1)
pnl = panels{6}(1,1);
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
pnl = panels{6}(1,2);
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


%% ========================================================================



%% add panel letters
font_size = 11;
axes(panels{1}(1))
text(-0.25,1.1, 'a', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{2}(1))
text(-0.4,1.1, 'b', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{3}(1))
text(-0.3,1.1, 'c', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{4}(1))
text(-0.25,1.1, 'd', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{5}(1))
text(-0.28,1.1, 'e', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{6}(1,1))
text(-0.115,1.1, 'f', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{6}(2,1))
text(-0.115,1.1, 'g', 'Units','normalized','FontWeight','bold','FontSize',font_size);

%% add line titles
% axes(panels{1}(1))
% text(-0.7,0.5,{"Both landing","balls"},'Units','normalized','FontSize',10,'HorizontalAlignment','center');
% axes(panels{4}(1))
% text(-0.7,0.5,{"Ball 1"},'Units','normalized','FontSize',10,'HorizontalAlignment','center');
% axes(panels{4}(2))
% text(-0.7,0.5,{"Ball 2"},'Units','normalized','FontSize',10,'HorizontalAlignment','center');

%%
fig_name = sprintf('%s_decoding_opt_%d',fig_name_str, params_opt);
file_out = fullfile(res_dir, fig_name);
print(gcf, file_out, '-dpdf', '-cmyk', '-painters');
disp('figure saved!')


%%
function str = genSignifStrAstricks(pval)
if pval < 0.001
    str = '***';
elseif pval < 0.01
    str = '**';
elseif pval < 0.05
    str = '*';
else
    str = 'n.s.';
end
end