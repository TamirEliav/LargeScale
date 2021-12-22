%%
clear
clc

%%
exp_ID='b0184_d191208';
epoch_type ='flight';
params_opt = 4;
exp = exp_load_data(exp_ID,'details','path','pos','flight');
decode = decoding_load_data(exp_ID, epoch_type, params_opt);

%%
t = decode.time;
pos_real = interp1(exp.pos.proc_1D.ts, exp.pos.proc_1D.pos, t, 'linear');
pos_estimate = decode.MAP_pos;
a = -decode.params.pos_bin_size/2;
b =  decode.params.pos_bin_size/2;
noise = a+(b-a).*rand(size(pos_estimate));
pos_estimate = pos_estimate + noise;
% pos_estimate = smoothdata(pos_estimate,'gaussian',50);
pos_diff = pos_estimate-pos_real;
pos_diff_smooth = smoothdata(pos_diff,'gaussian',0.5*decode.Fs);
g1 = interp1(decode.pos,1:length(decode.pos), pos_real, 'nearest','extrap');
g2 = interp1(decode.pos,1:length(decode.pos), pos_estimate, 'nearest','extrap');
g = zeros(size(t));
g(diff(t)>2e6/decode.Fs) = 1;
g = cumsum(g)+1;

%%
figure
subplot(211)
hold on
plot(t,pos_estimate,'r.','MarkerSize',5);
plot(t,pos_real,'k.','MarkerSize',5);
rescale_plot_data('x',[1e-6 t(1)]);
subplot(212)
hold on
plot(t,pos_diff,'.')
plot(t,pos_diff_smooth,'r')
ylim([-5 5])
rescale_plot_data('x',[1e-6 t(1)]);
linkaxes(findall(gcf,'type','axes'),'x')
zoom on

%%
figure
hold on
x = decode.pos;
y = splitapply(@median, pos_diff, g1);
err = splitapply(@nansem, pos_diff, g1);
splitapply(@(x,y)(plot(x,y,'.')), pos_real,pos_diff,g);
shadedErrorBar(x,y,err)
ylim([-5 5])
zoom on

%%








%%
