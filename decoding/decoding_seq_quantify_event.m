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
nPos = size(prob,1);
nTime = size(prob,2);
if isfield(opts,'pos');  opts.dx  = median(diff(opts.pos));  else opts.pos  = [1:nPos] .* opts.dx; end
if isfield(opts,'time');  opts.dt  = median(diff(opts.time));  else opts.time = [1:nTime] .* opts.dt; end
opts.dt = opts.dt * opts.time_units_factor;

%% sequence line fitting (radon with convolution)
L = round(opts.radius / opts.dx);
L = 2*L+1;
win = ones(L,1);
prob2 = imfilter(prob,win,0,'same','conv');
thetas = linspace(0,180,1000);
slopes = 1./tan(deg2rad(thetas));
[R, xp] = radon(prob2,thetas);
[Rmax, max_IX] = max(R,[],'all');
[max_IX_offset, max_IX_theta] = ind2sub(size(R),max_IX);
theta = thetas(max_IX_theta);
score = Rmax / nTime; % TODO: consider dividing by the number of VALID time bins (sum(~isnan(xx)))
offset = xp(max_IX_offset);
slope = slopes(max_IX_theta);
vel = slope * opts.dx / opts.dt;
compression = abs(vel) / opts.RW_speed;
xc = floor((nTime+1)./2);
yc = floor((nPos +1)./2);
xx = 1:nTime;
yy = slope.*(xx - xc) + yc - offset;
invalid_yy = (yy<0.5) | (yy>nPos+0.5);
xx(invalid_yy)=nan;
yy(invalid_yy)=nan;

%% sequence line fitting (radon)
% based (with some modifications) on code from here:
% https://github.com/DavidTingley/RnR_methods/blob/9fc48e7e90043b02012adf76605df35d2c8526bc/utilities/Pr2Radon.m
% % % radius_bins = opts.radius/opts.dx;
% % % thetas = linspace(0,180,1000);
% % % [R, xp] = radon(prob,thetas);
% % % clear x y
% % % y(1) = floor((nPos +1)./2);
% % % x(1) = floor((nTime+1)./2);
% % % [Y,I] = max(R,[],1);
% % % [~,pk] = max(Y);
% % % slopes = 1./tan(deg2rad(thetas));
% % % % vels = slopes .* dx/dt;
% % % scores = zeros(size(thetas));
% % % xx = 1:nTime;
% % % pos_IX = 1:nPos;
% % % switch 2
% % %     case 1
% % %         theta_IX_to_run_over = 1:length(thetas);
% % %     case 2
% % %         % this option is (1) faster and (2) more robust against short line
% % %         % fitting (note that we set nans wherever the line is outside of
% % %         % the image!)
% % %         theta_IX_to_run_over = pk;
% % %     case 3
% % %         % TODO: consider adding a third option (as done in Denovelis 2021):
% % %         % convolve the prob matrix with a box kernel (sum up over
% % %         % neighoboring positions) and then apply Radon. it will:
% % %         % 1) smooth the data which will give a better estimate of the
% % %         % best fitted line.
% % %         % 2) calculates the replay score immediately instead of doing that
% % %         % with the for-loop below
% % % end
% % % for ii_theta = theta_IX_to_run_over 
% % % %     angle = theta(ii_theta);
% % %     slope = slopes(ii_theta);
% % %     offset = xp(I(ii_theta));
% % % %     y(2) = y(1) + offset*sin(deg2rad(-angle));
% % % %     x(2) = x(1) + offset*cos(deg2rad(-angle));
% % % %     slope = 1/tan(deg2rad(angle));
% % %     yy = slope*(xx - x(1)) + y(1) - offset;
% % %     mask = pos_IX' < (yy+radius_bins) & pos_IX' > (yy-radius_bins);
% % %     mask = double(mask);
% % %     mask(:,~any(mask)) = nan;
% % %     scores(ii_theta) = nanmean(sum(prob.*mask,1));
% % % end
% % % [~, max_score_IX] = max(scores);
% % % score = scores(max_score_IX);
% % % slope = slopes(max_score_IX);
% % % offset = xp(I(max_score_IX));
% % % yy = slope.*(xx - x(1)) + y(1) - offset;
% % % invalid_yy = (yy<0.5) | (yy>nPos+0.5);
% % % xx(invalid_yy)=nan;
% % % yy(invalid_yy)=nan;

%%
res = struct();
res.slope = slope;
res.offset = offset; % offset from the center pixel (see doc radon)
res.theta = theta;
res.score = score;
res.compression = compression;              % this is defined directly by the fitted line
res.xx = xx;
res.yy = yy;
res.distance = range(yy)*opts.dx;           % this is defined by the valid pixels 
res.duration = range(xx)*opts.dt;           % this is defined by the valid pixels 
res.speed = res.distance / res.duration;    % this is defined by the valid pixels 
res.vel = vel;                              % this is defined directly by the fitted line

start_IX = find(~isnan(res.xx),1,'first');
end_IX = find(~isnan(res.xx),1,'last');
res.start_pos = interp1(1:nPos, opts.pos, res.yy(start_IX),'linear','extrap');
res.end_pos = interp1(1:nPos, opts.pos, res.yy(end_IX),'linear','extrap');
res.middle_pos = mean([res.start_pos res.end_pos]);
pos_edges = centers2edges(opts.pos);
res.env_edges = pos_edges([1 end]);
res.env_size = nPos * opts.dx;
res.env_size = range(pos_edges);
res.start_pos_norm = interp1(res.env_edges, [0 1], res.start_pos, 'linear','extrap');
res.end_pos_norm = interp1(res.env_edges, [0 1], res.end_pos, 'linear','extrap');
res.middle_pos_norm = interp1(res.env_edges, [0 1], res.middle_pos, 'linear','extrap');
res.distance_norm = range([res.start_pos_norm res.end_pos_norm]);
res.start_ts = opts.time(start_IX);
res.end_ts = opts.time(end_IX);
res.direction = sign(res.end_pos-res.start_pos);
if contains(opts.state, 'Outbound'); res.state_direction =  1; end
if contains(opts.state, 'Inbound');  res.state_direction = -1; end
if ~isfield(res,'state_direction'); res.state_direction = nan; end % TODO: we should get the "state direction" also in the bayseian decoder case!
res.forward = res.direction == res.state_direction;
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
