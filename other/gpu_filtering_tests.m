%% 
clear
clc

%%
fs = 2000;
dt = 1/fs;
T = 2*60*60;
% T = 10*60;
% T=1;
t = 0:dt:T;
f = [1:16]';
noise_std = 0.1;
x = cos(2*pi*f.*t);
x = x + noise_std*randn(size(x));

%% filter params
passband = [6 10];
Wn   = passband / (fs/2);
b = fir1(100, Wn,'bandpass');
a = 1;

%% cpu filtering
tic 
x_filtfilt = filtfilt(b,a,x')';
toc

%% cpu filtering (parfor)
tic 
x_filtfilt2 = zeros(size(x));
parfor ii=1:size(x,1)
    x_filtfilt2(ii,:) = filtfilt(b,a,x(ii,:));
end
toc

%% cpu filtering (filter-filter)
tic 
x_filter = filter(b,a,x')';
x_filterfilter = filter(b,a,fliplr(x_filter)')';
toc

%% gpu filtering
tic
xgpu = gpuArray(x);
for ch = 1:size(xgpu,1)
    xgpu(ch,:) = filter(b,a,xgpu(ch,:));
    xgpu(ch,:) = filter(b,a,fliplr(xgpu(ch,:)));
end
toc

%%
figure
hold on
plot(t,x(8,:))
plot(t,x(16,:))
% plot(t,x_filtfilt(8,:))
% plot(t,x_filtfilt(16,:))
plot(t,x_filterfilter(8,:))
plot(t,x_filterfilter(16,:))
legend()

%%
figure
hold on
plot(t,x(8,:))
plot(t,x(16,:))
plot(t,xgpu(8,:))
plot(t,xgpu(16,:))


%%
figure
hold on
plot(t,x)
plot(t,x_filt)
plot(t,x_filtfilt)
plot(t,x_filterfilter)
legend('x','filt','filtfilt','filterfilter')