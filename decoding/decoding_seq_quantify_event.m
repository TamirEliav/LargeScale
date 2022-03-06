function res = decoding_seq_quantify_event(prob,opts)
arguments 
    prob
    opts.dx = 1
    opts.dt = 1
    opts.radius = 1
    opts.pos
    opts.time
    opts.state = 'N/A'
    opts.RW_speed = 1
    opts.time_units_factor = 1;
end

%%
% if isfield(opts,'pos');  pos  = opts.pos;  else pos  = [1:size(prob,1)] .* opts.dx; end
% if isfield(opts,'time'); time = opts.time; else time = [1:size(prob,2)] .* opts.dt; end

%%
if isfield(opts,'pos');  opts.dx  = median(diff(opts.pos));  else opts.pos  = [1:size(prob,1)] .* opts.dx; end
if isfield(opts,'time');  opts.dt  = median(diff(opts.time));  else opts.time = [1:size(prob,2)] .* opts.dt; end

%% sequence line fitting (radon)
% based (with some modifications) on code from here:
% https://github.com/DavidTingley/RnR_methods/blob/9fc48e7e90043b02012adf76605df35d2c8526bc/utilities/Pr2Radon.m
radius_bins = opts.radius/opts.dx;
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
    case 3
        % TODO: consider adding a third option (as done in Denovelis 2021):
        % convolve the prob matrix with a box kernel (sum up over
        % neighoboring positions) and then apply Radon. it will:
        % 1) smooth the data which will give a better estimate of the
        % best fitted line.
        % 2) calculates the replay score immediately instead of doing that
        % with the for-loop below
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
res.theta = theta(max_score_IX);
res.score = score;
res.xx = xx;
res.yy = yy;
res.distance = range(yy)*opts.dx;
res.duration = range(xx)*opts.dt * opts.time_units_factor;
res.speed = res.distance / res.duration;

start_IX = find(~isnan(res.xx),1,'first');
end_IX = find(~isnan(res.xx),1,'last');
res.start_pos = interp1(1:size(prob,1), opts.pos, res.yy(start_IX),'linear','extrap');
res.end_pos = interp1(1:size(prob,1), opts.pos, res.yy(end_IX),'linear','extrap');
res.middle_pos = mean([res.start_pos res.end_pos]);
res.direction = sign(res.end_pos-res.start_pos);
if contains(opts.state, 'Outbound'); res.state_direction =  1; end
if contains(opts.state, 'Inbound');  res.state_direction = -1; end
if ~isfield(res,'state_direction'); res.state_direction = nan; end % TODO: we should get the "state direction" also in the bayseian decoder case!
res.forward = res.direction == res.state_direction;
res.compression = res.speed / opts.RW_speed;
res.start_ts = opts.time(start_IX);
res.end_ts = opts.time(end_IX);
[KS, HPD, sparsity] = decoding_calc_prob_confidence(prob(:,~isnan(res.xx)));
%     [KS, HPD, sparsity] = decoding_calc_prob_confidence(prob);
res.confidence_KS = mean(KS);
res.confidence_HPD = mean(HPD);
res.confidence_sparsity = mean(sparsity);

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
