%%
clear
clc

%% create figure
fig_size = [15 7];
figure('Units','centimeters','Position',[5 5 fig_size], 'PaperPosition', [0 0 fig_size])
pnl = panel();
pnl.pack('h',2);
pnl(1).pack('v',3);
pnl(2).pack('v',[20 60 20]); % workaround...
pnl.margin = [5 5 5 10];
pnl.de.margin = [20 5 0 5];
LINE_WIDTH = 1.5;

%% create data
x = 0:0.1:200;
CA3_maps(1,:) = gaussmf(x, [10 40]);
CA3_maps(2,:) = gaussmf(x, [2  80]);
CA3_maps(3,:) = gaussmf(x, [5  150]);

%% CA3
CA3_colors = {'r';'g';'b'}
for ii = 1:3
    pnl(1,ii).select(); hold on;
    axis off
    plot(x, CA3_maps(ii,:),'Color',CA3_colors{ii}, 'LineWidth', LINE_WIDTH);
    text(0.85,0.7,sprintf('cell %d',ii),'Units','normalized');
end

%% CA1
pnl(2,2).select(); hold on;
axis off
plot(x, sum(CA3_maps,1), 'Color','k', 'LineWidth', LINE_WIDTH);

%%
h=pnl(1).title('CA3'); h.FontSize=14;
h=pnl(2).title('CA1'); h.FontSize=14;

%% add arrows
annotation(gcf, 'arrow', [0.47 0.58], [0.72 0.55]);
annotation(gcf, 'arrow', [0.47 0.58], [0.50 0.50]);
annotation(gcf, 'arrow', [0.47 0.58], [0.28 0.45]);

%% save figure
dir_out = 'L:\Analysis\Results\midterm';
mkdir(dir_out);
filename = fullfile(dir_out, 'fig7')
saveas(gcf,filename,'fig')
saveas(gcf,filename,'tif')
saveas(gcf,filename,'pdf')



%%
