%% Replay - Fig supp 4 - reverse/forward + future+past + takeoff/landing/midair + ball1/ball2
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

%% define output files
res_dir =  'L:\paper_replay\figures';
mkdir(res_dir)
fig_name_str = 'Fig_supp_4';
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
panels{1}(1) = axes('position', [3  19 4 4]);
panels{2}(1) = axes('position', [9  19 3 4]);
panels{3}(1) = axes('position', [14 19 3 4]);
panels{4}(1) = axes('position', [3  13 4 4]);
panels{5}(1) = axes('position', [9  13 3 4]);
panels{6}(1) = axes('position', [14 13 3 4]);
panels{4}(2) = axes('position', [3   7 4 4]);
panels{5}(2) = axes('position', [9   7 3 4]);
panels{6}(2) = axes('position', [14  7 3 4]);
sdf=[panels{:}]
for ii = 1:length(sdf)
    sdf(ii).Position(1) = sdf(ii).Position(1) + 2;
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
disp('Comparing (immediate) Forward vs. Reverse')
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
p = myBinomTest(counts,sum(counts),expectedProportion,'one').*length(counts);
p(p==0)=realmin;
xG = 1:length(G);

Y = proportions;
h=bar(xG,Y);
h.BarWidth = 0.75;
h.FaceColor = 0.5*[1 1 1];
hax=gca;
hax.XTick = xG;
hax.XTickLabel = TakeLandCategoriesOrder;
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
ylabel('Fraction', 'Units','normalized', 'Position',[-0.2 .5]);
hax=gca;
hax.TickLength(1) = [0.025];
hax.XRuler.TickLabelGapOffset = -1;
hax.YRuler.TickLabelGapOffset = -1;
hax.XLim = [min(xG) max(xG)] + [-1 1].*0.5;


%% takeoff vs landing - separated by balls
for ii_ball = 1:2
axes(panels{4}(ii_ball))
cla
hold on

events_ball_IX = [events.ball_num] == ii_ball;
gTakeLand_by_ball = gTakeLand(events_ball_IX);
G = categorical(categories(gTakeLand_by_ball));
gcounts = tabulate(gTakeLand_by_ball);
counts = [gcounts{:,2}];
proportions = counts./sum(counts);
expectedProportion = [.05 .9 0.05];
expectedCounts = expectedProportion .* sum(counts);
p = myBinomTest(counts,sum(counts),expectedProportion,'one').*length(counts);
p(p==0)=realmin;
xG = 1:length(G);

Y = proportions;
h=bar(xG,Y);
h.BarWidth = 0.75;
h.FaceColor = 0.5*[1 1 1];
hax=gca;
hax.XTick = xG;
hax.XTickLabel = TakeLandCategoriesOrder;
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
ylabel('Fraction', 'Units','normalized', 'Position',[-0.2 .5]);
hax=gca;
hax.TickLength(1) = [0.025];
hax.XRuler.TickLabelGapOffset = -1;
hax.YRuler.TickLabelGapOffset = -1;
hax.XLim = [min(xG) max(xG)] + [-1 1].*0.5;

end


%% Forward vs Reverse
axes(panels{2}(1))
cla
hold on
G = categorical(categories(gForRev));
gcounts = tabulate(gForRev);
counts = [gcounts{:,2}];
proportions = counts./sum(counts);
expectedProportion = [.5 .5];
expectedCounts = expectedProportion .* sum(counts);
p = myBinomTest(counts,sum(counts),expectedProportion,'one').*length(counts);
p(p==0)=realmin;
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
    y = Y(ii)+0.001*max(Y);
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
hax.YRuler.TickLabelGapOffset = -1;
hax.XLim = [min(xG) max(xG)] + [-1 1].*0.5;
hax.YLim = [0 0.7];



%% Forward vs Reverse - separated by balls
for ii_ball = 1:2
axes(panels{5}(ii_ball))
cla
hold on

events_ball_IX = [events.ball_num] == ii_ball;
gForRev_by_ball = gForRev(events_ball_IX);
G = categorical(categories(gForRev_by_ball));
gcounts = tabulate(gForRev_by_ball);
counts = [gcounts{:,2}];
proportions = counts./sum(counts);
expectedProportion = [.5 .5];
expectedCounts = expectedProportion .* sum(counts);
p = myBinomTest(counts,sum(counts),expectedProportion,'one').*length(counts);
p(p==0)=realmin;
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
    y = Y(ii)+0.001*max(Y);
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
hax.YRuler.TickLabelGapOffset = -1;
hax.XLim = [min(xG) max(xG)] + [-1 1].*0.5;
hax.YLim = [0 0.7];

end



%% Future vs Past
axes(panels{3}(1))
cla
hold on
G = categorical(categories(gPastFuture));
gcounts = tabulate(gPastFuture);
counts = [gcounts{:,2}];
proportions = counts./sum(counts);
expectedProportion = [.5 .5];
expectedCounts = expectedProportion .* sum(counts);
p = myBinomTest(counts,sum(counts),expectedProportion,'one').*length(counts);
p(p==0)=realmin;
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
    y = Y(ii)+0.001*max(Y);
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
hax.YRuler.TickLabelGapOffset = -1;
hax.XLim = [min(xG) max(xG)] + [-1 1].*0.5;
hax.YLim = [0 0.7];

%% Future vs Past - separated by balls
for ii_ball = 1:2
axes(panels{6}(ii_ball))
cla
hold on

events_ball_IX = [events.ball_num] == ii_ball;
gPastFuture_by_ball = gPastFuture(events_ball_IX);
G = categorical(categories(gPastFuture_by_ball));
gcounts = tabulate(gPastFuture_by_ball);
counts = [gcounts{:,2}];
proportions = counts./sum(counts);
expectedProportion = [.5 .5];
expectedCounts = expectedProportion .* sum(counts);
p = myBinomTest(counts,sum(counts),expectedProportion,'one').*length(counts);
p(p==0)=realmin;
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
    y = Y(ii)+0.001*max(Y);
%     text(x,y, "{\itP} < 10^{"+ceil(log10(p(ii)))+"}", horizontalAlignment = 'center',VerticalAlignment='bottom',FontSize=6);
    signifStr = genSignifStrAstricks(p(ii));
    str = {signifStr;"{\itn} = " + counts(ii)};
    text(x,y, str, horizontalAlignment = 'center',VerticalAlignment='bottom',FontSize=8);
end
% xlabel('Replay speed (m/s)', 'Units','normalized', 'Position',[0.5 -0.12]);
ylabel('Fraction', 'Units','normalized', 'Position',[-0.2 .5]);
% text(1.6,0.5,"Ball "+ii_ball,'Units','normalized','FontSize',10,'HorizontalAlignment','center');
hax=gca;
hax.TickLength(1) = [0.025];
hax.XRuler.TickLabelGapOffset = -1;
hax.YRuler.TickLabelGapOffset = -1;
hax.XLim = [min(xG) max(xG)] + [-1 1].*0.5;
hax.YLim = [0 0.7];

end

%% add panel letters
font_size = 11;
axes(panels{1}(1))
text(-0.3,1.1, 'a', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{2}(1))
text(-0.3,1.1, 'b', 'Units','normalized','FontWeight','bold','FontSize',font_size);
axes(panels{3}(1))
text(-0.3,1.1, 'c', 'Units','normalized','FontWeight','bold','FontSize',font_size);

%% add line titles
axes(panels{1}(1))
text(-0.7,0.5,{"Both landing","balls"},'Units','normalized','FontSize',10,'HorizontalAlignment','center');
axes(panels{4}(1))
text(-0.7,0.5,{"Ball 1"},'Units','normalized','FontSize',10,'HorizontalAlignment','center');
axes(panels{4}(2))
text(-0.7,0.5,{"Ball 2"},'Units','normalized','FontSize',10,'HorizontalAlignment','center');

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