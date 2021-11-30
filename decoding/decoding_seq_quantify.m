function decoding_seq_quantify(exp_ID, epoch_type, params_opt, event_type)
arguments
    %% 
    exp_ID = 'b9861_d180526'
    epoch_type {mustBeMember(epoch_type,{'sleep','rest','flight'})} = 'sleep'
    params_opt = 11;
    event_type {mustBeMember(event_type,{'PE','posterior','ripples','MUA'})} = 'posterior'
end

%% IN/OUT folders
dir_IN = 'F:\sequences\decoded';
dir_OUT = 'F:\sequences\events_quantification';
mkdir(dir_OUT);

%% load data
exp = exp_load_data(exp_ID, 'details','path','rest','ripples','MUA','PE','pos','flight');
events = decoding_load_events(exp_ID, epoch_type, params_opt, event_type);
decode_filename = fullfile(dir_IN, epoch_type, exp_ID, sprintf('%s_%s_opt_%d.nc',exp_ID,epoch_type,params_opt));
decode = decoding_read_decoded_file(decode_filename);

%% params & consts
dt = 1/decode.Fs;
dx = decode.params.pos_bin_size;
radius = 5;
% RW_speed = median(abs([exp.flight.speed_traj.vel_median_median]));
RW_speed = median(abs([exp.flight.FE.vel]));

%% calc
clear res_all
for ii_event = 1:length(events)
    %%
    event = events(ii_event);
    IX = [event.start_IX:event.end_IX];
    posterior = decode.posterior(:,:,IX);
    prob = squeeze(posterior(:,event.state_num,:));
    res = event_fit_line(prob,dt,dx,radius);
    start_IX = find(~isnan(res.xx),1,'first');
    end_IX = find(~isnan(res.xx),1,'last');
    res.start_pos = interp1(1:size(prob,1), decode.pos, res.yy(start_IX),'linear','extrap');
    res.end_pos = interp1(1:size(prob,1), decode.pos, res.yy(end_IX),'linear','extrap');
    res.middle_pos = mean([res.start_pos res.end_pos]);
    res.direction = sign(res.end_pos-res.start_pos);
    if contains(decode.state(event.state_num), 'Outbound')
        res.state_direction = 1;
    else
        res.state_direction = -1;
    end
    res.forward = res.direction == res.state_direction;
    res.compression = res.speed / RW_speed;
    res.start_ts = decode.time(IX(start_IX));
    res.end_ts = decode.time(IX(end_IX));
    [KS, HPD, sparsity] = decoding_calc_prob_confidence(prob(:,~isnan(res.xx)));
%     [KS, HPD, sparsity] = decoding_calc_prob_confidence(prob);
    res.confidence_KS = mean(KS);
    res.confidence_HPD = mean(HPD);
    res.confidence_sparsity = mean(sparsity);
    
    res_all(ii_event) = res;
end
% add quantification to events struct
if exist('res_all','var')
    [events.seq] = disperse(res_all);
else
    [events.seq] = struct();
end

%% save all events to mat file
params = struct();
params.decode = decode.params;
params.dx = dx;
params.dt = dt;
params.radius = radius;
params.RW_speed = RW_speed;
filename = fullfile(dir_OUT, sprintf('%s_events_%s_dec_prm_%d_%s',exp_ID,epoch_type,params_opt,event_type));
save(filename, 'events','params');

%% plot scatter matrix figure
% % % % T=struct2table(res_all);
% % % % % T1 = T(:,["score","compression"]);
% % % % % g = kmeans(table2array(T1),2);
% % % % T2 = T(:,["score","range","duration","start_pos","end_pos","compression","confidence_KS"]);
% % % % % gplotmatrix(table2array(T2),[],g,[],[],[],[],[],T2.Properties.VariableNames);
% % % % fig = figure;
% % % % fig.WindowState = 'maximized';
% % % % [h,ax,bigax]=gplotmatrix(table2array(T2),[],[],[],[],[],[],[],T2.Properties.VariableNames);
% % % % set([ax.XLabel ax.YLabel],'Interpreter','none');
% % % % saveas(gcf, filename, 'jpg');

end


function res = event_fit_line(prob,dt,dx,radius)
%% sequence line fitting (radon)
% based (with some modifications) on code from here:
% https://github.com/DavidTingley/RnR_methods/blob/9fc48e7e90043b02012adf76605df35d2c8526bc/utilities/Pr2Radon.m
radius_bins = radius/dx;
theta = linspace(0,180,1000);
[R, xp] = radon(prob,theta);
clear x y
y(1) = floor((size(prob,1)+1)./2);
x(1) = floor((size(prob,2)+1)./2);
[Y,I] = max(R,[],1);
[~,pk] = max(Y);
slopes = 1./tan(deg2rad(theta));
% vels = slopes .* dx/dt;
scores = zeros(size(theta));
xx = 1:size(prob,2);
pos_IX = 1:size(prob,1);
switch 2
    case 1
        theta_IX_to_run_over = 1:length(theta);
    case 2
        % this option is (1) faster and (2) more robust against short line
        % fitting (note that we set nans wherever the line is outside of
        % the image!)
        theta_IX_to_run_over = pk;
end
for ii_theta = theta_IX_to_run_over 
%     angle = theta(ii_theta);
    slope = slopes(ii_theta);
    offset = xp(I(ii_theta));
%     y(2) = y(1) + offset*sin(deg2rad(-angle));
%     x(2) = x(1) + offset*cos(deg2rad(-angle));
%     slope = 1/tan(deg2rad(angle));
    yy = slope*(xx - x(1)) + y(1) - offset;
    mask = pos_IX' < (yy+radius_bins) & pos_IX' > (yy-radius_bins);
    mask = double(mask);
    mask(:,~any(mask)) = nan;
    scores(ii_theta) = nanmean(sum(prob.*mask,1));
end
[~, max_score_IX] = max(scores);
score = scores(max_score_IX);
slope = slopes(max_score_IX);
offset = xp(I(max_score_IX));
yy = slope.*(xx - x(1)) + y(1) - offset;
invalid_yy = (yy<0.5) | (yy>size(prob,1)+0.5);
xx(invalid_yy)=nan;
yy(invalid_yy)=nan;

res = struct();
res.slope = slope;
res.offset = offset; % offset from the center pixel (see doc radon)
res.score = score;
res.xx = xx;
res.yy = yy;
res.distance = range(yy)*dx;
res.duration = range(xx)*dt;
res.speed = res.distance / res.duration;

%%
% % % figure
% % % hold on
% % % % xxx = theta;
% % % % xxx = slopes;
% % % % xxx = [res.slope];
% % % % xxx = xxx .* dx./dt;
% % % xxx = vels;
% % % plot(xxx, [res.score]./max([res.score]))
% % % plot(xxx, Y./max(Y))
% % % xlim([-1 1].*500)

%%
% figure
% subplot(211)
% hold on
% imagesc(prob)
% plot(x(1),y(1),'.r')
% % plot(x(2),y(2),'*r')
% plot(xx,yy,'r-')
% axis ij
% % axis equal
% subplot(212)
% imagesc(R,'XData',theta,'YData',xp);
% % axis xy
% % axis equal
% xlabel('Theta')
% ylabel('xp')

end






%%
