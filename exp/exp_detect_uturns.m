function exp_detect_uturns(exp_ID)

%% load exp data
exp = exp_load_data(exp_ID,'pos','rest');
prm = PARAMS_GetAll();
dir_out = 'L:\Analysis\Results\exp\uturns';

%% params
min_dist_from_balls = 2; % m
min_time_around_uturns = 1; % sec
min_vel_continuous = 1; % m/sec

%% arrange relevant data
session_ti = exp_get_sessions_ti(exp_ID,'Behave','Behave_6m','Light1','Light2');
% session_ti = exp_get_sessions_ti(exp_ID,'Behave_6m','Light1','Light2');
session_ti(any(isnan(session_ti),2),:) = []; % remove nan in case of missing sessions
IX = get_data_in_ti(exp.pos.proc_1D.ts, session_ti);
pos.ts = exp.pos.proc_1D.ts(IX);
pos.pos = exp.pos.proc_1D.pos(IX);
pos.vel = exp.pos.proc_1D.vel_csaps(IX);
pos.fs = exp.pos.proc_1D.fs;

%% detect U-turns away from the balls:
U_turn_IX = find(abs(diff(sign(pos.vel)))==2); % vel flips sign
U_turn_IX(pos.pos(U_turn_IX)<(exp.rest.balls_loc(1)+min_dist_from_balls)) = [];
U_turn_IX(pos.pos(U_turn_IX)>(exp.rest.balls_loc(2)-min_dist_from_balls)) = [];

%% require continous movement before and after the u_turn:
samples_around = 1:(min_time_around_uturns*pos.fs);

% remove U-turns at the edges of the data:
U_turn_IX(U_turn_IX>(length(pos.vel)-max(samples_around))) = [];
U_turn_IX(U_turn_IX<max(samples_around)) = [];

cont_dir_before = zeros(1,length(U_turn_IX));
cont_dir_after = zeros(1,length(U_turn_IX));
high_vel_before = zeros(1,length(U_turn_IX));
high_vel_after = zeros(1,length(U_turn_IX));
uturn_directions = zeros(1,length(U_turn_IX)); % based on the direction BEFORE the uturn
for ii_uturn = 1:length(U_turn_IX)
    % condition on the same velocity sign for t_around_u time before and 
    % after the U-turn:
    cont_dir_before(ii_uturn) = sum(abs(diff(sign(pos.vel(U_turn_IX(ii_uturn)-samples_around)))))==0;
    cont_dir_after(ii_uturn) = sum(abs(diff(sign(pos.vel(U_turn_IX(ii_uturn)+samples_around)))))==0;
    % condition on abs high velocity for t_around_u time before and after 
    % the U-turn (not landing on the floor times):
    high_vel_before(ii_uturn) = median(abs(pos.vel(U_turn_IX(ii_uturn)-samples_around)))>=min_vel_continuous;
    high_vel_after(ii_uturn) = median(abs(pos.vel(U_turn_IX(ii_uturn)+samples_around)))>=min_vel_continuous;
    vel_sign_before_uturn = sign(median(pos.vel(U_turn_IX(ii_uturn)-samples_around)));
    switch vel_sign_before_uturn
        case 1
            uturn_directions(ii_uturn) = 1;
        case -1
            uturn_directions(ii_uturn) = 2;
    end
end
cont_before = and(cont_dir_before,high_vel_before);
cont_after = and(cont_dir_after,high_vel_after);
TF = and(cont_before,cont_after);

U_turn_IX(~TF) = [];
uturn_directions(~TF) = [];

%% 
uturns = struct();
uturns.ts = pos.ts(U_turn_IX);
uturns.pos = pos.pos(U_turn_IX);
uturns.dir = uturn_directions;

%% plot 
fig=figure;
fig.WindowState = 'maximized';
hax=gca;
hax.ColorOrder = [1 0 0; 0 0 1];
hold on
plot(pos.ts, pos.pos,'.k','MarkerSize',1)
scatter(uturns.ts, uturns.pos,50,uturns.dir,'filled');
cmap = [1 0 0; 0 0 1];
colormap(cmap)
rescale_plot_data('x',[1e-6 pos.ts(1)])
xlabel('Time (s)')
ylabel('Position (m)')
title({'uturns detection';exp_ID},'Interpreter','none')
file_name = fullfile(dir_out, [exp_ID '_exp_uturns_detection']);
saveas(fig,file_name,'fig')
saveas(fig,file_name,'jpg')

%% save rest on the balls results to mat file
file_name = fullfile(dir_out ,[exp_ID '_exp_uturns']);
save(file_name,'uturns');

%%











%%