%% speed data for Yonatan
clear
clc

%% load data
load('L:\processed_data_structs\exps_including_FE.mat');
prm = PARAMS_GetAll();
option = 2;

%%
speed_data = struct();

%%
% figure
for ii_exp = 1:length(exps)
    exp = exps(ii_exp);
    speed_data(ii_exp).bat_num = exp.details.batNum;
    FE = exp.flight.FE;
    dirs = [1 -1];
    for ii_dir = 1:length(dirs)
        switch option
            case 1
                speed_traj = exp.flight.speed_traj(ii_dir);
                IX = get_data_in_ti(speed_traj.bins_centers,prm.fields.valid_speed_pos);
                vels = speed_traj.vel_median(IX);
            case 2
                FE_dir_IX = [FE.direction] == dirs(ii_dir);
                FE_dir = FE(FE_dir_IX);
                FE_dir([FE_dir.distance] < 100) = [];
                vels = [FE_dir.vel];
        end
        speed_data(ii_exp).vel_mean(ii_dir) = mean(vels);
        speed_data(ii_exp).vel_std(ii_dir) = std(vels);
        speed_data(ii_exp).direction(ii_dir) = ii_dir;
    end
end

%%
save('L:\DATA_for_people\Yonatan\speed_data','speed_data');

%%
clear
clc
load('L:\DATA_for_people\Yonatan\speed_data.mat');

speeds = abs(cat(1,speed_data.vel_mean));
bats = [speed_data.bat_num];
[G grps]= findgroups(bats);
y = splitapply(@mean,speeds,G');
% err = splitapply(@std,speeds,G');
err = splitapply(@nansem,speeds,G');
figure
hold on
x = 1:length(grps);
x=repmat(x',1,size(y,2))
bar(x,y);
% errorbar(x,y,err)
xticks(1:length(grps))
xticklabels(grps)
xlabel('bat')
ylabel('speed (m/s)')
legend('dir1' , 'dir2','Location','bestoutside')
saveas(gcf,'L:\DATA_for_people\Yonatan\speed_data','jpg');

%%








%%
