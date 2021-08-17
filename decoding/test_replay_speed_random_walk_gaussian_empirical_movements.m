%%
clear
clc

%%
dir_OUT = 'D:\sequences\seq_uri_eden\results';

%% rat
real_world_speed = 17;
variance = 6;
bin_size = 3;
dt = 0.002;
fs = 1/dt;
x = linspace(0,10,100);
% x = 0:bin_size:100;
v = x./dt;
mu = 0;
sigma = sqrt(variance);
p1 = normcdf( x,mu,sigma);
p2 = normcdf(-x,mu,sigma);
p = p1 - p2;
median_model_speed = interp1(p,v,0.5);
replay_speed_factor = median_model_speed / real_world_speed;

track_length = 300;
T = 60*30;
t = 0:dt:T;
pos = t.*real_world_speed;
pos = mod(pos,track_length);
pos = [pos flip(pos)];
t = [t t+T+dt];

figure
subplot(211)
yyaxis left
plot(x,p,'.-')
ylabel('cdf')
yyaxis right
plot(x,1-p,'.-')
ylabel('1-cdf')
xlabel('distance after 1 timebin (cm)')
subplot(212)
yyaxis left
plot(v./100,p,'.-')
ylabel('cdf')
yyaxis right
plot(v./100,1-p,'.-')
ylabel('1-cdf')
xlabel('speed (m/s)')
text(0.5,0.4, ...
    {sprintf('median real speed = %.2gcm/s',real_world_speed);
     sprintf('median model speed = %.2gm/s',median_model_speed/100);
     sprintf('speed factor = %.2g',replay_speed_factor);},...
    'Units','normalized','Interpreter','none')
sgtitle('Rat model speed')
figname = fullfile(dir_OUT, "rat_model_speed_cdf");
saveas(gcf,figname,'jpg')

%%
[N,XEDGES,YEDGES] = histcounts2(pos(1:end-1),pos(2:end),'BinWidth',bin_size,'Normalization','count',...
    'XBinLimits',[100 200],'YBinLimits',[100 200]);
BIN_SIZE = mean(diff(XEDGES));
N = N ./ sum(N,1);
% replay_speed_factor = 550;
N2 = N^replay_speed_factor;
N2 = N2 ./ sum(N2,1);
N_diag_proj_prob = arrayfun(@(k)(mean(diag(N,k))),1-size(N,2):size(N,1)-1);
N2_diag_proj_prob = arrayfun(@(k)(mean(diag(N2,k))),1-size(N,2):size(N,1)-1);
N_diag_proj_prob = N_diag_proj_prob ./ sum(N_diag_proj_prob);
N2_diag_proj_prob = N2_diag_proj_prob ./ sum(N2_diag_proj_prob);
N_diag_proj_prob = N_diag_proj_prob ./ BIN_SIZE;
N2_diag_proj_prob = N2_diag_proj_prob ./ BIN_SIZE;
diag_proj_lag = [-length(N):length(N)] .* BIN_SIZE;
diag_proj_lag([1 end])=[];
diag_model = normpdf(diag_proj_lag,0,sigma);
hf=figure;
hf.WindowState='maximized';
subplot(221)
imagesc(XEDGES,YEDGES,N)
colormap(flip(bone))
xlabel('position@t_i (m)')
ylabel('position@t_{i+1} (m)')
title('empirical')
subplot(222)
imagesc(XEDGES,YEDGES,N2)
title('empirical spedup')
xlabel('position@t_i (m)')
ylabel('position@t_{i+1} (m)')
subplot(212)
hold on
plot(diag_proj_lag, N_diag_proj_prob)
plot(diag_proj_lag, N2_diag_proj_prob)
plot(diag_proj_lag, diag_model)
title('Diagonal projection')
xlabel('Position lag (cm)')
ylabel('pdf')
legend('empirical','empirical_spedup','model','interpreter','none')
% xlim([-1 1].*30)
ylim([0 .4])
sum([N_diag_proj_prob;N2_diag_proj_prob;diag_model],2)
diag_proj_lag([1 end])
text(0.8,0.4, "bin_size="+bin_size,'Units','normalized','Interpreter','none')
figname = fullfile(dir_OUT, "rat_model_movement_dynamics_bin_size_"+bin_size+"cm");
saveas(gcf,figname,'jpg')

