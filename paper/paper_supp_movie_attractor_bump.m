%%
clear
clc

%% load data
load('L:\Misha_attractor\20200902__new_simulations\sim.mat')

%%
fig=figure;
fig.WindowState = 'maximized';
fig.Units='normalized';
pnls=[];
pnls(1) = axes('position', [0.05 0.65+0.05 0.15 0.25]);
pnls(2) = axes('position', [0.23 0.65+0.05 0.15 0.25]);
pnls(3) = axes('position', [0.41 0.65+0.05 0.15 0.25]);
pnls(4) = axes('position', [0.59 0.65+0.05 0.15 0.25]);
pnls(5) = axes('position', [0.77 0.65+0.05 0.15 0.25]);
pnls(6) = axes('position', [0.05 0.35+0.05 0.42 0.25]);
pnls(7) = axes('position', [0.50 0.35+0.05 0.42 0.25]);
pnls(8) = axes('position', [0.05 0.05+0.05 0.87 0.25]);

%%
seg = [ 0 80;
        80 160;
        160 240;
        240 320;
        320 400;
        0 200;
        200 400;
        0 400;];
bin_size = 0.5;
seg = seg .* bin_size;

%%
run = 10;
nloops = size(m,1);
% nloops = 10;
ylimits = [0 0.3];
ylimits = [0 0.15];
M(nloops) = struct('cdata',[],'colormap',[]);
for bin = 1:nloops
    for ii_net=1:size(ind,2)
        axes(pnls(ii_net));
        plot(th(:,ii_net), m(bin,ind(:,ii_net),run));
        cla('reset');
        hold on
        plot(th(:,ii_net).*bin_size, m(bin,ind(:,ii_net),run),'Color',0*[1 1 1]);
        if ii_net == size(ind,2)
            hl = xline(bin*bin_size);
            hl.LineWidth = 1.5;
            hl.Color = 'r';
            xlabel('Position (m)');
        end
        if ii_net == 6
            ylabel('Firing rate (a.u.)');
        end
        hax=gca;
        hax.XLim = seg(ii_net,:);
        hax.YLim = ylimits;
        hax.YTick = [];
        hax.TickLength(1) = 0.005 / max(hax.Position([3 4]));
        box off
        switch ii_net
            case {1,2,3,4,5}
                hax.XTick = seg(ii_net,1) : 20 : seg(ii_net,2);
            case {6,7}
                hax.XTick = seg(ii_net,1) : 25 : seg(ii_net,2);
            case 8
                hax.XTick = seg(ii_net,1) : 50 : seg(ii_net,2);
        end
    end
    drawnow
    M(bin) = getframe(gcf);
end

%% create AVI file
Q = 10;

res_dir = hc3_get_res_dir();
res_dir = fullfile(res_dir,'paper_figures');
mkdir(res_dir)
filename_str = "movie_run_"+run+"_Q_"+Q;
filename_out = fullfile(res_dir, filename_str);

v = VideoWriter(filename_out, 'MPEG-4');
v.Quality=Q;
open(v)
writeVideo(v,M);
close(v)

%%




