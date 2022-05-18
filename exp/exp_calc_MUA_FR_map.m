function exp_calc_MUA_FR_map(exp_ID)

%% load data
exp = exp_load_data(exp_ID,'details','MUA','flight','flight_6m','rest');

%%
if isfield(exp,'flight_6m')
    FE = [exp.flight_6m.FE];
else
    FE = [exp.flight.FE];
    FE([FE.distance]<100)=[];
end

%% arrange data
MUA = exp.MUA;
FR_during_FE = interp1([MUA.t],[MUA.FR],[FE.ts],'linear','extrap');
pos = [FE.pos];
% pos_norm = interp1(exp.rest.balls_loc',[0 1],pos,'linear','extrap');
gdir = sign([FE.vel]);

%% calc maps
nBins = 1000;
% EDGES = linspace(0,1,nBins+1);
EDGES = linspace(exp.rest.balls_loc(1), exp.rest.balls_loc(2), nBins+1);
CENTERS = edges2centers(EDGES);
maps = [];
for direction = [1 -1]
    dir_IX = gdir==direction;
    pos_dir = pos(dir_IX);
    vals = FR_during_FE(dir_IX);
    bins = discretize(pos_dir,EDGES);
    invalid_IX = isnan(bins);
    bins(invalid_IX)=[];
    vals(invalid_IX)=[];
    maps(end+1,:) = accumarray(bins', vals', [length(CENTERS) 1], @mean, nan);
end

%% plot
% plot(CENTERS, MUA_FR_map','linewidth',2)

%% save res
MUA_FR_map = struct();
MUA_FR_map.bin_centers = CENTERS;
MUA_FR_map.maps = maps;
fn = 'MUA_FR_map';
dir_out = fullfile('L:\Analysis\Results\exp\', fn);
mkdir(dir_out);
file_name = fullfile(dir_out ,[exp_ID '_exp_' fn]);
save(file_name,'MUA_FR_map');
