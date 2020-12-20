%% Large Scale - Fig. S19 - Attractor network model

%%
clear 
clc

%% load data
load('L:\Misha_attractor\20200902__new_simulations\sim.mat');
load('L:\Misha_attractor\20200902__new_simulations\sim_res.mat');

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S19';
fig_caption_str = 'attractor_model';
log_name_str = [fig_name_str '_log_file' '.txt'];
log_name_str = strrep(log_name_str , ':', '-');
log_name_str = strrep(log_name_str , ' ', '_');
log_name_out = fullfile(res_dir, log_name_str);

%% open log file
diary off
diary(log_name_out)
diary on
disp('Log file');
disp(['created: ', datestr(clock)]);
disp('======================================================');
disp([fig_name_str ': ' fig_caption_str]);
disp('======================================================');
disp('');

%% create figure
% =========================================================================
% figure_size_cm = [21.0 29.7]; % ~A4
figure_size_cm = [21.6 27.9]; % ~US letter
figure ;
% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',7);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'DefaultAxesUnits','centimeters');
set(gcf,'PaperType','usletter')
% set(gcf,'PaperType','<custom>');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 figure_size_cm]);
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]); % position on screen...
set(gcf, 'Renderer', 'painters');
set(groot, 'defaultAxesTickDir', 'out');
set(groot,  'defaultAxesTickDirMode', 'manual');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none', 'FitBoxToText','on');
pause(0.2); % workaround to solve matlab automatically changing the axes positions...

% create panels
pnls(1,1) = axes('Units','normalized', 'position', [0.05 0.85 0.17 0.07]);
pnls(1,2) = axes('Units','normalized', 'position', [0.23 0.85 0.17 0.07]);
pnls(1,3) = axes('Units','normalized', 'position', [0.41 0.85 0.17 0.07]);
pnls(1,4) = axes('Units','normalized', 'position', [0.59 0.85 0.17 0.07]);
pnls(1,5) = axes('Units','normalized', 'position', [0.77 0.85 0.17 0.07]);
pnls(1,6) = axes('Units','normalized', 'position', [0.05 0.75 0.44 0.07]);
pnls(1,7) = axes('Units','normalized', 'position', [0.50 0.75 0.44 0.07]);
pnls(1,8) = axes('Units','normalized', 'position', [0.05 0.65 0.89 0.07]);

pnls(2,1) = axes('Units','normalized', 'position', [0.05 0.85-0.33 0.17 0.07]);
pnls(2,2) = axes('Units','normalized', 'position', [0.23 0.85-0.33 0.17 0.07]);
pnls(2,3) = axes('Units','normalized', 'position', [0.41 0.85-0.33 0.17 0.07]);
pnls(2,4) = axes('Units','normalized', 'position', [0.59 0.85-0.33 0.17 0.07]);
pnls(2,5) = axes('Units','normalized', 'position', [0.77 0.85-0.33 0.17 0.07]);
pnls(2,6) = axes('Units','normalized', 'position', [0.05 0.75-0.33 0.44 0.07]);
pnls(2,7) = axes('Units','normalized', 'position', [0.50 0.75-0.33 0.44 0.07]);
pnls(2,8) = axes('Units','normalized', 'position', [0.05 0.65-0.33 0.89 0.07]);

panel_C(1) = axes('position', [ 3  3 2.5 3]);
panel_C(2) = axes('position', [ 7  3 2.5 3]);
% panel_C(3) = axes('position', [ 11  11.5 2.5 3]);
    

%% attractor snapshots
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
run = 10;
pos_snapshots = [100 300];
ylimits = [0 0.3];
% ylimits = [0 0.15];
for ii_bin = 1:length(pos_snapshots)
    bin = pos_snapshots(ii_bin);
    for ii_net=1:size(ind,2)
        axes(pnls(ii_bin,ii_net));
        cla('reset');
        hold on
        plot(th(:,ii_net).*bin_size, m(bin,ind(:,ii_net),run),'Color',0*[1 1 1]);
        hl = xline(bin*bin_size);
        hl.LineWidth = 1.5;
        hl.Color = 'r';
        if ii_net == size(ind,2)
            xlabel('Position (m)','Units','normalized','Position',[0.5 -0.25]);
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
                if bin*bin_size > hax.XLim(1) & bin*bin_size < hax.XLim(2)
                    text(bin*bin_size, ylimits(2)+0.15*range(ylimits),...
                        'Input position','HorizontalAlignment','center','Color','r');
                end
            case {6,7}
                hax.XTick = seg(ii_net,1) : 25 : seg(ii_net,2);
            case 8
                hax.XTick = seg(ii_net,1) : 50 : seg(ii_net,2);
        end
        text(hax.XLim(1)+1, 0.9*range(hax.YLim), "A"+ii_net, 'HorizontalAlignment','left','FontSize', 8);
        ticks = hax.XTick;
        ticklabels = hax.XTickLabel;
        hax.XTickLabel=[];
        offsets = zeros(size(ticks));
        offsets([1 end]) = [1 -1];
        offsets(ticks>=100) = offsets(ticks>=100).*2;
        text(ticks+offsets,repelem(-0.18*range(ylimits),length(ticks)),ticklabels, 'FontSize',7.5,'HorizontalAlignment','center');
    end
end

%%
axes(pnls(1,3));
text(0.5,1.35, 'Continuous attractors: Bumps of activity', 'Units','normalized','FontWeight','bold','HorizontalAlignment','center');
axes(pnls(1,1));
text(-0.15,1.2, 'A', 'Units','normalized','FontWeight','bold');
axes(pnls(2,1));
text(-0.15,1.2, 'B', 'Units','normalized','FontWeight','bold');

%%
axes(panel_C(1));
cla('reset')
hold on
text(-0.45, 1.25, 'C', 'Units','normalized','FontWeight','bold');
text( 0.07, 1.25, 'Control model: Independent attractors', 'Units','normalized','FontWeight','bold');
h=histogram(res0.fnumber);
h.Normalization = 'probability';
h.BinEdges = 0.5+[0:5];
% h.DisplayStyle = 'stairs';
h.LineWidth = 2;
h.EdgeColor = 'k';
h.FaceColor = 0.5*[1 1 1];
xlabel('No. of fields')
ylabel('Fraction of cells');
hax=gca;
hax.YScale = 'linear';
% hax.YLim = 10.^[-3 0];
% hax.XTick=[];
% hax.YTick=10.^[-3:0];
hax.TickLength = [0.03 0.03];
hax.XRuler.TickLabelGapOffset = -.5;

%%
axes(panel_C(2));
cla('reset')
hold on
seg_size = range(seg,2);
x = res0.fsizes_all ./ max(seg_size) .* bin_size;
x = x .* seg_size';
h=histogram(x(:));
% h.NumBins = 200;
h.BinWidth = .01;
h.Normalization = 'probability';
% h.DisplayStyle = 'stairs';
h.LineWidth = 0.1;
h.EdgeColor = 'k';
h.FaceColor = 0.5*[1 1 1];
xlabel('Field size (m)')
ylabel('Fraction of fields')
hax=gca;
% hax.YScale = 'log';
hax.YScale = 'linear';
hax.XTick = [0 20 40];
% hax.YTick=[];
hax.XLim = [0 45];
% hax.YLim = [1e-3 2e-1];
hax.TickLength = [0.03 0.03];
hax.XRuler.TickLabelGapOffset = -.5;

% 35 17 6

%%
% axes(panel_C(3));
% cla('reset')
% hold on
% % text(-0.25,1.2, 'E', 'Units','normalized','FontWeight','bold');
% h=histogram(res0.sratio);
% h.Normalization = 'pdf';
% h.DisplayStyle = 'stairs';
% h.LineWidth = 2;
% h.EdgeColor = 'c';
% xlabel('Fields ratio')
% ylabel('PDF')
% xlim([0 20])
% hax=gca;
% hax.XScale = 'log';
% hax.YScale = 'log';
% % hax.XTick=[];
% % hax.YTick=[];
% hax.XLim = [0 20];
% hax.YLim = [1e-3 1e0];




%% print/save the figure
fig_name_out = fullfile(res_dir, [fig_name_str]);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');



%%





