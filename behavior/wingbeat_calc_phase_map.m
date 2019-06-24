function wingbeat_calc_phase_map(exp_ID)

%% load exp data
exp = exp_load_data(exp_ID,'flight');
prm = PARAMS_GetAll();

%% load data

% wingbeat
wingbeat_file = fullfile('L:\Analysis\Results\wingbeats',[exp_ID '_wingbeats.mat' ]);
if ~exist(wingbeat_file,'file')
    wingbeats_detect(exp_ID)
end
load(wingbeat_file)

% landmarks
LM_file = 'L:\DATA\0148_Boson\calibrations\landmarks_copy_from_Mario.xlsx';
LM = readtable(LM_file);
% calc LM linearized position
load(exp_path.calib_tunnel_file );
mapxy = [LM.Y LM.X];
[~,~,t] = distance2curve(calib_tunnel.curvexy,mapxy,'linear');
LM.pos_linearized = t .* calib_tunnel.tunnel_length;

%% remove some flights
if 0
IX = 60:length(flights.direction);
flights.direction(IX) = [];
flights.durations(IX) = [];
flights.end_IX(IX) = [];
flights.end_ts(IX) = [];
flights.start_IX(IX) = [];
flights.start_ts(IX) = [];
end 

%%
wingbeat.pos = interp1(pos.ts_nlg_usec, pos.pos_linearized, wingbeat.ts);

% take wingbeat for each flight
pos_bin_size = 0.05;
pos_bin_edges = 0:pos_bin_size:200;
pos_bin_centers = (pos_bin_edges(1:end-1) + pos_bin_edges(2:end))/2;
wingbeat_phase_map = zeros( length(flights.durations), length(pos_bin_centers) );
wingbeat_wavelnegth_map = zeros( length(flights.durations), length(pos_bin_centers) );
speed_map = zeros( length(flights.durations), length(pos_bin_centers) );
clear wingbeats_by_flight
for ii_flight = 1:length(flights.durations)
    wingbeat_IX = find(  wingbeat.ts > flights.start_ts(ii_flight) &...
                wingbeat.ts < flights.end_ts(ii_flight));
    wingbeats_by_flight(ii_flight).wingbeat_IX = wingbeat_IX;
    wingbeats_by_flight(ii_flight).wingbeat_ts = wingbeat.ts(wingbeat_IX);
    wingbeats_by_flight(ii_flight).wingbeat_pos = wingbeat.pos(wingbeat_IX);
    [N,~,bins] = histcounts(wingbeat.pos(wingbeat_IX), pos_bin_edges);
    wingbeat_pos_bin_IX = find(N);
    wingbeat_phase_accum = interp1(wingbeat_pos_bin_IX, [1:length(wingbeat_pos_bin_IX)].*2*pi, 1:length(N),'linear',0);
    wingbeat_phase_map(ii_flight,:) = mod(wingbeat_phase_accum,2*pi);
%     wingbeat_wavelnegth_map(ii_flight,:) = [0 pos_bin_size*2*pi./diff(wingbeat_phase_accum)];
    wingbeat_wavelnegth_map(ii_flight,:) =interp1(bins(2:end),abs(diff(bins)*pos_bin_size), 1:length(N), 'nearest','extrap');
    flight_pos_IX = flights.start_IX(ii_flight) : flights.end_IX(ii_flight);
    speed_map(ii_flight,:) = interp1(   pos.pos_linearized_csaps(flight_pos_IX),...
                                        pos.vel_linearized_csaps(flight_pos_IX),...
                                        pos_bin_centers,'nearest','extrap');
end

clear flights_dir_IX
flights_dir_IX{1} = find(flights.direction==1);
flights_dir_IX{2} = find(flights.direction==-1);

% % % % delays = zeros(size(wingbeat_phase_map,1),size(wingbeat_phase_map,1));
% % % % phase_diff = zeros(size(wingbeat_phase_map,1),size(wingbeat_phase_map,1), size(wingbeat_phase_map,2));
% % % % for ii_flight1 = 1:size(wingbeat_phase_map,1)
% % % %     for ii_flight2 = 1:size(wingbeat_phase_map,1)
% % % %         delays(ii_flight1,ii_flight2) = finddelay(wingbeat_phase_map(ii_flight1,:),wingbeat_phase_map(ii_flight2,:),20);
% % % %         phase_diff(ii_flight1,ii_flight2,:) = circ_dist( wingbeat_phase_map(ii_flight1,:), wingbeat_phase_map(ii_flight2,:));
% % % %     end
% % % % end
% % % % delays = diag(delays,1);
% % % % % plot(delays,'.-')

%%
pvals=[];
A_r = [];
for ii_dir = 1:2
    A_r(ii_dir,:) = circ_r(wingbeat_phase_map(flights_dir_IX{ii_dir},:));
    for ii_bin = 1:length(pos_bin_centers)
        pvals(ii_dir,ii_bin) = circ_rtest( wingbeat_phase_map(flights_dir_IX{ii_dir},ii_bin));
    end
end

%%
% dir_IN = 'L:\Analysis\Results\Wingbeats';
% file_IN = fullfile(dir_IN,[exp_ID '_wingbeat_maps' ]);
% load(file_IN)

%%
figure
p=panel();
fig_size = [45 20];
set(gcf, 'Units','Centimeters', 'Position',[2 2 fig_size], 'PaperPosition',[1 1 fig_size]);
p.pack('h',2,'v',{0.1 0.1 0.1 0.25 0.25 0.25});
p.margintop = 15;
p.marginleft = 27;
p.marginbottom = 24;
p.ch.margin = [17 10 10 10];
% p(1,2).margintop = 30;
% p(2,2).margintop = 30;
% p(,2).margintop = 30;
h=p.title(exp_ID);
set(h,'position',[0.5 1.025],'FontWeight','Bold','FontSize',12,'Interpreter','none');
direction_str = {'>>>>>>>>>>';'<<<<<<<<<<'};
wavelength_lim = [0.85 1.45];
phase_lim = [0 2*pi];
speed_lim = [5.5 7.5];
pos_lim = [0 200];
pval_lim = [0 0.52];
plot_LM = 0;
for ii_dir = 1:2
    
    p(ii_dir,1).select(); hold on
    bin_IX = find( pvals(ii_dir,:) < 0.05 );
    ylimits = pval_lim;
    plot(pos_bin_centers, A_r(ii_dir,:),'.')
    plot(pos_bin_centers(bin_IX),A_r(ii_dir,bin_IX),'.r')
    ylim(ylimits)
    xlim(pos_lim)
    ylabel('reighley')
    title(['direction: ' direction_str{ii_dir}])
    if plot_LM
        for ii_LM = 1:size(LM,1)
            x = LM.pos_linearized(ii_LM);
            plot([x x],ylimits,'--m');
            text(x,ylimits(1)-0.2*diff(ylimits),LM.name(ii_LM),'rotation',-70,'interpreter','none')
        end
    end
    
    p(ii_dir,2).select(); hold on
    x = pos_bin_centers;
    y = trimmean(wingbeat_wavelnegth_map(flights_dir_IX{ii_dir},:),5);
    err = std(wingbeat_wavelnegth_map(flights_dir_IX{ii_dir},:),0,1) ./ sqrt(length(flights_dir_IX{ii_dir}));
    shadedErrorBar(x,y,err);
    ylim(wavelength_lim)
    ylabel({'wave-';'length (m)'} )
    xlim(pos_lim)
    
    p(ii_dir,3).select(); hold on
    x = pos_bin_centers;
    y = trimmean(abs(speed_map(flights_dir_IX{ii_dir},:)),5);
    err = std(speed_map(flights_dir_IX{ii_dir},:),0,1) ./ sqrt(length(flights_dir_IX{ii_dir}));
    shadedErrorBar(x,y,err);
    ylim(speed_lim)
    ylabel('Speed (m/s)')
    
    p(ii_dir,4).select()
    imagesc(pos_bin_centers, 1:length(flights_dir_IX{ii_dir}),wingbeat_phase_map(flights_dir_IX{ii_dir},:))
    xlim(pos_lim)
    ylim([0 length(flights_dir_IX{ii_dir})] )
    set(gca,'clim',phase_lim)
    ylabel('#flight')
    title('wingbeat phase map')
    ax_pos = get(gca,'position');
    colorbar_ax_pos = [0.025 ax_pos(2) 0.01 0.1];
    h=colorbar('location','WestOutside', 'position', colorbar_ax_pos);
    h.TickLength = 0.025;
    h.TickDirection = 'out';
    h.Ticks = phase_lim;
    h.TickLabels = {'0';'2{\pi}'};
    h.Label.String = 'phase ({\circ})';
    h.Label.Position = [-0.5 mean(phase_lim)];
    
    p(ii_dir,5).select()
    imagesc(pos_bin_centers, 1:length(flights_dir_IX{ii_dir}),wingbeat_wavelnegth_map(flights_dir_IX{ii_dir},:))
    xlim(pos_lim)
    ylim([0 length(flights_dir_IX{ii_dir})] )
    set(gca,'clim',wavelength_lim)
    ylabel('#flight')
    title('wingbeat wavelength map')
    ax_pos = get(gca,'position');
    colorbar_ax_pos = [0.025 ax_pos(2) 0.01 0.1];
    h=colorbar('location','WestOutside', 'position', colorbar_ax_pos);
    h.TickLength = 0.025;
    h.TickDirection = 'out';
    h.Ticks = wavelength_lim;
    h.Label.String = {'wave-';'length (m)'};
    h.Label.Position = [-0.3 mean(wavelength_lim)];
%     h.Label.FontSize = 6;
%     h.Label.String = {'wave-';'length(m)'};
%     h.Label.Units = 'Normalized';
%     h.Label.Position = [0.5 2];
%     h.Label.Rotation = 0;

    p(ii_dir,6).select()
    imagesc(pos_bin_centers, 1:length(flights_dir_IX{ii_dir}),abs(speed_map(flights_dir_IX{ii_dir},:)));
    xlim(pos_lim)
    ylim([0 length(flights_dir_IX{ii_dir})] )
    set(gca,'clim',speed_lim)
    xlabel('Position (m)')
    ylabel('#flight')
    title('bat speed map')
    ax_pos = get(gca,'position');
    colorbar_ax_pos = [0.025 ax_pos(2) 0.01 0.1];
    h=colorbar('location','WestOutside', 'position', colorbar_ax_pos);
    h.TickLength = 0.025;
    h.TickDirection = 'out';
    h.Ticks = speed_lim;
    h.Label.String = {'Speed';'(m/s)'};
    h.Label.Position = [-0.3 mean(speed_lim)];
end

set(gcf, 'Units','Centimeters', 'Position',[2 2 fig_size], 'PaperPosition',[1 1 fig_size]);

dir_OUT = 'L:\Analysis\Results\wing_beat_artifacts_position_locking\';
saveas(gcf,fullfile(dir_OUT,['wingbeat_phase_reighley_' exp_ID]),'tif')
% saveas(gcf,fullfile(dir_OUT,['wingbeat_phase_reighley_' exp_ID]),'pdf')

%%
