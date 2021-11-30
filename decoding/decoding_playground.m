
%%
prob=decode.posterior_pos;
sdf=prob(:,1000);
sdf2=sort(sdf,'descend');
% sdf2=sort(sdf,'ascend');
sdf3=cumsum(sdf2);
plot(sdf3,'o-')
thr = 0.8;
find(sdf3>thr,1,'first')

%% load exp data
exp_ID = 'b9861_d180526';
epoch_type  = 'sleep';
params_opt = 11;
event_type = 'posterior';
exp = exp_load_data(exp_ID, 'details','path','rest','ripples','MUA','PE','pos','flight');
events = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
decode = decoding_load_data(exp_ID, epoch_type, params_opt);

%% load all cells data
cells_all = load('L:\processed_data_structs\cells_bat_200m.mat');

%% arrange exp_ID cells
cells = cells_all.cells;
details = [cells.details];
cells = cells(contains({details.exp_ID},exp_ID));
signif = cat(1,cells.signif);
ii_dir = 2;
cells = cells([signif(:,ii_dir).TF]);
details = [cells.details];
fields = cat(1,cells.fields);
fields = fields(:,ii_dir);
n = cumsum(cellfun(@length,fields));
grps = interp1(n,1:length(n),1:n(end), 'next','extrap');
fields = [fields{:}];
[fields.cell_num] = disperse(grps);
[~,sorted_IX] = sort([fields.loc],'ascend');
fields = fields(sorted_IX);
cells_spikes = cellfun(@(cell_ID)(cell_load_data(cell_ID,'spikes')),{details.cell_ID});
[cells.spikes] = disperse([cells_spikes.spikes]);

%%
event_num = 49;
event = events(event_num);
IX = [event.start_IX:event.end_IX];
prob = decode.posterior_pos(:,IX);
time = decode.time(IX);
ti = [event.start_ts event.end_ts];
figure
hold on
imagesc('CData',prob, 'XData',time, 'YData',decode.pos);
axis xy 
for ii_cell = 1:length(cells)
    cell = cells(ii_cell);
    [~,~,spikes_ts] = get_data_in_ti([cell.spikes.ts], ti);
    x = spikes_ts;
    fields_IX = find([fields.cell_num] == ii_cell);
%     y = fields_IX;
    y = [fields(fields_IX).loc];
    y = repelem(y,length(x),1);
    plot(x,y,'|k');
end

%% quantifying decoding confidence
% event_num = 49;
% event_num = 50;
% % event_num = 57;
event = events(event_num);
IX = [event.start_IX:event.end_IX];

% event1_num = 49;
% event2_num = 50;
% event1 = events(event1_num);
% event2 = events(event2_num);
% IX1 = [event1.start_IX:event1.end_IX];
% IX2 = [event2.start_IX:event2.end_IX];
% IX = [event1.start_IX:event2.end_IX];
% events12_lag = event2.start_IX - event1.start_IX;

ii_state = event.state_num;
prob=squeeze(decode.posterior(:,ii_state,:));
sdf=prob(:,IX);
sdf2=sort(sdf,1,'descend');
% sdf2=sort(sdf,'ascend');
sdf3=cumsum(sdf2,1);
% plot(sdf3,'o-');
thr = 0.0;
nbins = size(sdf,1);
m = zeros(1,size(sdf3,2));
m2 = zeros(1,size(sdf3,2));
sparsity = zeros(1,size(sdf3,2));
test_cdf = cdf('Discrete Uniform',1:nbins,nbins);
for ii_bin = 1:size(sdf3,2)
    m(ii_bin) = find(sdf3(:,ii_bin)>thr,1,'first');
    m2(ii_bin) = max(abs(sdf3(:,ii_bin)'-test_cdf));
    sparsity(ii_bin) = cumputeSparsity(sdf3(:,ii_bin));
end
% normalize m and m2
m = interp1([1 size(sdf,1)], [1 0], m, 'linear');
m2_min = max(abs(test_cdf-test_cdf)); % simply zero...
m2_max = max(abs(ones(1,nbins)-test_cdf));
m3 = interp1(linspace(m2_min,m2_max,nbins), linspace(0,1,nbins), m2, 'linear');
m4 = max(abs(sdf3'-test_cdf),[],2); % use this!
% plot
figure
hax=[];
hax(1)=subplot(311);
yyaxis left
hold on
plot(m)
plot(m2)
plot(m3)
% plot(m4)
% plot(sparsity)
legend({sprintf('%d%% HPD (norm.)',thr*100);'ksstat';'ksstat (norm.)';'sparsity'})
yyaxis right
plot(sparsity)
title("ksstat")
hax(2)=subplot(312);
hold on
imagesc(sdf)
plot(event.seq.xx, event.seq.yy,'r-')
% plot(event1.seq.xx, event1.seq.yy,'r-')
% plot(event2.seq.xx+events12_lag, event2.seq.yy,'r-')
cmap = bone;
cmap = flipud(cmap);
colormap(cmap);
axis xy
hax(3)=subplot(313);
hold on
M = squeeze(decode.likelihood(:,ii_state,IX));
% M = squeeze(sum(decode.likelihood(:,:,IX),2));
M = M./sum(M);
% M = smoothdata(M,1,'movmean',10);
% M = smoothdata(M,2,'movmean',50);
% M = smoothdata(M,1,'gaussian',10);
% M = smoothdata(M,2,'gaussian',50);
imagesc(M);
plot(event.seq.xx, event.seq.yy,'r-')
% plot(event1.seq.xx, event1.seq.yy,'r-')
% plot(event2.seq.xx+events12_lag, event2.seq.yy,'r-')
colormap(cmap);
axis xy
climits = prctile(M(:),[0 100]+[1 -1].*20);
set(gca,'CLim', climits)
linkaxes(hax,'x');

%%
M2 = M;
M2 = smoothdata(M2,1,'movmean',10);
M2 = smoothdata(M2,2,'movmean',10);

calc_plot_radon(M);
calc_plot_radon(sdf);
calc_plot_radon(M2);



%%
figure
hold on
[R, Xp] = radon(M,theta);
% [R, Xp] = radon(smoothdata(M,'movmean',5),theta);
[~,max_IX] = max(R,[],'all','linear');
[r,c]=ind2sub(size(R),max_IX);
imagesc('CData',R, 'XData',theta, 'YData',Xp);
plot(theta(c),Xp(r),'or')

figure
hold on
[R, Xp] = radon(sdf,theta);
[~,max_IX] = max(R,[],'all','linear');
[r,c]=ind2sub(size(R),max_IX);
imagesc('CData',R, 'XData',theta, 'YData',Xp);
plot(theta(c),Xp(r),'or')

%%
figure
hold on
plot(m,m2,'.')
plot(m,m3,'.')
xlim([0.8 1])
ylim([0.8 1])
axis equal
refline(1,0)

%%
fig=figure;
fig.WindowState = 'maximized';
pnl = panel();
pnl.pack('v',[.3 .7])
pnl(1).pack('h',[.85 .05 .05]);
pnl(2).pack('h',[.85 .05 .05]);
pnl.de.margin = 10;
pnl(2,1).select();
hold on
splitapply(@(IX,x)(plot(IX,x,'.')),1:length(decode.time), decode.MAP_pos, decode.MAP_state_IX)
gaps_IX = find(diff(decode.time) > median(diff(decode.time)));
for ii_gap = 1:length(gaps_IX)
    hl = xline(gaps_IX(ii_gap),'g','time gap');
    hl.LineWidth = 2;
end
rescale_plot_data('x',[1/decode.Fs 0]);
[hl, objhl] = legend(decode.state,'Interpreter','none','Location','northoutside');
hl.Position = [0.87 0.87 0.1 0.1];
icons = findobj(objhl,'Marker','.');
set(icons,'MarkerSize',20);
nItems = length(hl.PlotChildren);
% legend_marker_size = 20;
% [objhl(nItems+[1:nItems].*2).MarkerSize] = disperse(legend_marker_size.*ones(1,nItems));
arrayfun(@(y,str)(yline(y,'-',str,'Color',0.5.*[1 1 1],'LineWidth',0.5,'LabelVerticalAlignment','middle')),[exp.LM.pos_proj], string({exp.LM.name}))
xlabel('Time (s)')
ylabel('Position (m)')
pnl(2,2).select();
hold on
histogram(decode.MAP_pos,       'BinWidth',decode.params.pos_bin_size, 'Orientation','horizontal','Normalization','probability');
arrayfun(@(y,str)(yline(y,'Color',0.5.*[1 1 1],'LineWidth',0.5)),[exp.LM.pos_proj], string({exp.LM.name}))
pnl(2,3).select();
hold on
histogram(decode_flight.MAP_pos,'BinWidth',decode.params.pos_bin_size, 'Orientation','horizontal','Normalization','probability');
arrayfun(@(y,str)(yline(y,'Color',0.5.*[1 1 1],'LineWidth',0.5)),[exp.LM.pos_proj], string({exp.LM.name}))
linkaxes(pnl(2).de.axis,'y')
pnl(1,1).select();
% vel = gradient(decode.MAP_pos) .* decode.Fs;
% vel = smoothdata(vel,'gaussian',10);
MAP_pos_csaps_smoothing = 0.0005;
MAP_pos_csaps_smoothing = 0.00005;
MAP_pos_csaps_smoothing = 0.000005;
MAP_pos_csaps_pp = csaps(1:length(decode.time), decode.MAP_pos, MAP_pos_csaps_smoothing);
MAP_vel_csaps = fnval( fnder(MAP_pos_csaps_pp), 1:length(decode.time)) * decode.Fs;
hold on
% splitapply(@(IX,x)(plot(IX,x,'.')),1:length(decode.time), MAP_vel_csaps , decode.MAP_state_IX)
MUA_ts_IX = interp1(decode.time, 1:length(decode.time), exp.MUA.t, 'nearest');
ripples_ts_IX = interp1(decode.time, 1:length(decode.time), exp.ripples.t, 'nearest');
PE_events_ts_IX = interp1(decode.time, 1:length(decode.time), [exp.PE.thr.peak_ts], 'nearest');
plot(MUA_ts_IX, exp.MUA.zFR);
% plot(ripples_ts_IX , exp.ripples.zpripple_all);
plot(PE_events_ts_IX, [exp.PE.thr.peak_zFR],'*r');
legend("Firing rate","Population Events")
rescale_plot_data('x',[1/decode.Fs 0]);
ylim([0 10])
xlabel('Time (s)')
% ylabel('Velocity (m/s)')
ylabel('Firing rate (z)')
linkaxes([pnl(1,1).axis pnl(2,1).axis],'x')
% xlim([355 359])
% pnl(1,2).select();
% hax=gca;
% v=violinplot(MAP_vel_csaps, decode.MAP_state_IX);
% xticklabels(decode.state);
% hax.XTickLabelRotation = -45;
% hax.XAxis.TickLabelInterpreter = 'none';
% [v.ViolinColor] = disperse( hax.ColorOrder')
% linkaxes([pnl(1,1).axis pnl(1,2).axis],'y')
% ylim([-1 1].*300);
% % ylim([]);
% title(MAP_pos_csaps_smoothing)

%%
figure
plot(randn(10,6),'.')
[hl objhl ] = legend("data "+[1:6]);
[objhl(6+[1:6].*2).MarkerSize] = disperse(10.*ones(1,6));

%%
I = randn(30,200);
[R xp] = radon(I);
whos R xp I

%%
start_pos_opts = linspace(min(decode.pos),max(decode.pos),1000);
vel_opts = linspace(-1,1,10000).*8*200;
median(diff(vel_opts));
radius = 5;
for ii_start_pos = 1:length(start_pos_opts)
    for ii_vel = 1:length(vel_opts)
        
    end
end


%%
I = zeros(100,100);
I(25:75, 25:75) = 1;
theta = 0:180;
[R,xp] = radon(I,theta);
[~,max_IX] = max(R,[],'all','linear')
[r c] = ind2sub(size(R),max_IX);

figure
subplot(211)
imagesc(I)
refline()
subplot(212)
imagesc(R)
xlabel('\theta (degrees)')
ylabel('x''')
colormap(gca,hot), colorbar
iptsetpref('ImshowAxesVisible','on')


%% sequence line fitting (radon)
theta = linspace(0,180,1000);
[R xp] = radon(prob,theta);
clear x y
y(1) = floor((size(prob,1)+1)./2);
x(1) = floor((size(prob,2)+1)./2);
[Y,I] = max(R,[],1);
[~,pk] = max(Y);
for ii_theta = 1:length(theta)
    ii_theta = pk;
    angle = theta(ii_theta);
    offset = xp(I(ii_theta));
    y(2) = y(1) + offset*sin(deg2rad(-angle));
    x(2) = x(1) + offset*cos(deg2rad(-angle));
    coeffs1 = polyfit(x, y, 1);
    xx = 1:size(prob,2);
    yy = (-1/coeffs1(1))*(xx - x(1)) + y(1) - offset;
    coeffs2 = polyfit(xx, yy, 1);
    slope = coeffs2(1);
%     for ii_ts = 1:size(prob,2)
%         radius = 5;
%         pos_IX = 
%         prob
%     end
end

%%
figure
subplot(211)
hold on
imagesc(prob)
plot(x(1),y(1),'.r')
plot(x(2),y(2),'*r')
plot(xx,yy,'r-')
axis ij
% axis equal
subplot(212)
imagesc(R,'XData',theta,'YData',xp);
% axis xy
% axis equal
xlabel('Theta')
ylabel('xp')


%%
function calc_plot_radon(M)
fig=figure;
fig.WindowState = 'maximized';
hold on
theta = linspace(0,180,1000);
[R, Xp] = radon(M,theta);
[~,max_IX] = max(R,[],'all','linear');
[r,c]=ind2sub(size(R),max_IX);
imagesc('CData',R, 'XData',theta, 'YData',Xp);
plot(theta(c),Xp(r),'or')
end








%%
