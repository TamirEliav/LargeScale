function data = paper_fig_8_arrange_real_data()

%% load data
load('L:\processed_data_structs\cells_bat_200m.mat')

%% work only with signif cells
signif = cat(1,cells.signif);
signif = arrayfun(@(x)(x.TF),signif);
signif = any(signif,2);
cells(~signif)=[];

%% Basic population stats (fields num/size/ratio) =========================
signif = cat(1,cells.signif);
signif = arrayfun(@(x)(x.TF),signif);
stats = [cells.stats];
stats_all = [stats.all];
stats_dir = cat(1,stats.dir);
fields = cat(1,cells.fields);
fields = fields(signif);
fields = [fields{:}];
fields([fields.in_low_speed_area])=[];

field_num = [stats_dir(signif).field_num];
field_size = [fields.width_prc];
ratio_LS = [stats_all.field_ratio_LS];
%     ratio_LS(isnan(ratio_LS))=[];
ratio_LS_with_1s = [stats_all.field_ratio_LS];
ratio_LS_with_1s(isnan(ratio_LS_with_1s))=1;



%% spatail periodicity (fft) ==============================================

%% define locations
entire_tunnel_edges = [10 187.5]; % valid speed pos
long_arm_edges     = [10 139.5]; % valid speed pos, but only the long arm part

%% get valid maps
signif = cat(1,cells.signif);
signif = arrayfun(@(x)(x.TF), signif);
x = cells(1).FR_map(1).all.bin_centers;
fs = 1 / cells(1).FR_map(1).all.bin_size;

% full FR PSTH maps
maps = cat(1,cells.FR_map);
maps = maps(signif);
maps = [maps.all];
maps = cat(1,maps.PSTH);
maps_long = maps;
% binarized maps
maps01 = nan([size(signif),length(x)]);
for ii_cell = 1:length(cells)
    for ii_dir = 1:2
        if signif(ii_cell,ii_dir)
            fields = cells(ii_cell).fields{ii_dir};
            fields_edges = cat(1,fields.edges_prc);
            IX = get_data_in_ti(x,fields_edges);
            maps01(ii_cell,ii_dir,:) = 0;
            maps01(ii_cell,ii_dir,IX) = 1;
        end
    end
end
TF=repmat(signif,1,1,length(x));
maps01 = maps01(TF);
maps01 = reshape(maps01,[],length(x));

% select locations
entire_tunnel_IX = x>entire_tunnel_edges(1) & x<entire_tunnel_edges(2);
long_arm_IX = x>long_arm_edges(1) & x<long_arm_edges(2);
maps(:,~entire_tunnel_IX)=[];
maps01(:,~entire_tunnel_IX)=[];
maps_long(:,~long_arm_IX)=[];
% get bat identity per map
details = [cells.details];
bats = [details.bat];
bats = repelem(bats',1,2);
maps_bat = bats(signif);

%% calc maps fft
n = size(maps,2);
n_long = size(maps_long,2);
maps_spec = fft(maps,[],2);
maps_spec = fftshift(maps_spec);
maps_spec = abs(maps_spec).^2/n;
maps_spec = maps_spec ./ mean(maps_spec,2);
maps01_spec = fft(maps01,[],2);
maps01_spec = fftshift(maps01_spec);
maps01_spec = abs(maps01_spec).^2/n;
maps01_spec = maps01_spec ./ mean(maps01_spec,2);
maps_long_spec = fft(maps_long,[],2);
maps_long_spec = fftshift(maps_long_spec);
maps_long_spec = abs(maps_long_spec).^2/n;
maps_long_spec = maps_long_spec ./ mean(maps_long_spec,2);
freq = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
freq_long = (-n_long/2:n_long/2-1)*(fs/n_long); % zero-centered frequency range

%% arrange ouput struct

% field related
data.field_num = field_num;
data.field_size = field_size;
data.ratio_LS = ratio_LS;
data.ratio_LS_with_1s = ratio_LS_with_1s;

% fft related 
data.freq = freq;
data.freq_long = freq_long;
data.maps_spec = maps_spec;
data.maps01_spec = maps01_spec;
data.maps_long_spec = maps_long_spec;
data.maps_bat = maps_bat;


%%







%%
