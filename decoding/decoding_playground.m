
%%
prob=decode.posterior_pos;
sdf=prob(:,1000);
sdf2=sort(sdf,'descend');
% sdf2=sort(sdf,'ascend');
sdf3=cumsum(sdf2);
plot(sdf3,'o-')
thr = 0.8;
find(sdf3>thr,1,'first')

%%
event_num = 49;
% event_num = 57;
event = events_all(event_num);
IX = [event.start_IX:event.end_IX];
prob=decode.posterior_pos;
sdf=prob(:,IX);
sdf2=sort(sdf,1,'descend');
% sdf2=sort(sdf,'ascend');
sdf3=cumsum(sdf2,1);
% plot(sdf3,'o-');
thr = 0.95;
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
m2_max = max(abs(ones(1,nbins)-test_cdf)); % simply zero...
m3 = interp1(linspace(m2_min,m2_max,nbins), linspace(0,1,nbins), m2, 'linear');
m4 = max(abs(sdf3'-test_cdf),[],2); % use this!
% plot
figure
hax=[];
hax(1)=subplot(211);
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
hax(2)=subplot(212);
imagesc(sdf)
cmap = bone;
cmap = flipud(cmap);
colormap(cmap);
axis xy
linkaxes(hax,'x');
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
pnl(1).pack('h',[.85 .1]);
pnl(2).pack('h',[.85 .1]);
pnl.de.margin = 20;
pnl(2,1).select();
hold on
splitapply(@(IX,x)(plot(IX,x,'.')),1:length(decode.time), decode.MAP_pos, decode.MAP_state_IX)
gaps_IX = find(diff(decode.time) > median(diff(decode.time)));
for ii_gap = 1:length(gaps_IX)
    hl = xline(gaps_IX(ii_gap),'g','time gap');
    hl.LineWidth = 2;
end
rescale_plot_data('x',[1/decode.Fs 0]);
hl=legend(decode.state,'Interpreter','none','Location','northoutside');
hl.Position = [0.87 0.87 0.1 0.1];
xlabel('Time (s)')
ylabel('Position (m)')
pnl(2,2).select();
histogram(decode.MAP_pos,'BinWidth',decode.params.pos_bin_size, 'Orientation','horizontal')
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
splitapply(@(IX,x)(plot(IX,x,'.')),1:length(decode.time), MAP_vel_csaps , decode.MAP_state_IX)
rescale_plot_data('x',[1/decode.Fs 0]);
xlabel('Time (s)')
ylabel('Velocity (m/s)')
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


%%
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









%%
