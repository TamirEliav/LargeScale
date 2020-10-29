%% Large Scale - Fig. 8 - Attractor network model

%%
clear 
clc

%% load data
load('L:\Misha_attractor\20200902__new_simulations\sim.mat');
load('L:\Misha_attractor\20200902__new_simulations\sim_res.mat');

%% options
multiple_panels = 1;

%% choose examples 
% 12 29 48 74 98 107 121 198 206 230 283 296 308 327 400
cell_examples_IX = [29 198 12 737 74];
connectivity_cell_IX = 29;

%% define output files
res_dir = hc3_get_res_dir();
res_dir = fullfile(res_dir,'paper_figures');
mkdir(res_dir)
fig_name_str = 'fig_8_attractor_model';
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
panel_A = axes('position', [ 3 20 11 4]);
if multiple_panels
    panel_B(1) = axes('position', [ 3 16 11 0.8]);
    panel_B(2) = axes('position', [ 3 17 11 0.8]);
    panel_B(3) = axes('position', [ 3 18 11 0.8]);
    panel_B(4) = axes('position', [ 3 19 11 0.8]);
    panel_B(5) = axes('position', [ 3 20 11 0.8]);
else
    panel_B = axes('position', [ 3 16 11 4]);
end
panel_C(1) = axes('position', [ 3   11.5 2.5 3]);
panel_C(2) = axes('position', [ 7   11.5 2.5 3]);
panel_C(3) = axes('position', [ 11  11.5 2.5 3]);
    

%% model schematics
axes(panel_A);
cla('reset')
hold on
text(-0.1,1.05, 'A', 'Units','normalized','FontWeight','bold');
clr = [ 1 0 0;
        0 1 0;
        0 0 1;
        1 0 1;
        0 1 1;
        0.5 0.5 1;
        1 0.5 0.5;
        1 1 0;];
seg = [ 0   80;
        80  160;
        160 240;
        240 320;
        320 400;
        0   200;
        200 400;
        0   400;
        ];
bin_size = 0.5;
seg = seg .* bin_size;
dlm = 2.5;
seg_margins = [
    0 -dlm;
    dlm -dlm;
    dlm -dlm;
    dlm -dlm;
    dlm 0;
    0 -dlm;
    dlm 0;
    0 0;
];
levels = [3 3 3 3 3 2 2 1];
for ii_seg=1:size(seg,1)
    plot(seg(ii_seg,:)+seg_margins(ii_seg,:),[levels([ii_seg ii_seg])],'k','linewidth',2,'Color','k');
%     sigma = 0.05;
%     x = linspace(seg(ii,1),seg(ii,2),1000);
%     y = 0.65.*gaussmf(x,[sigma*range(seg(ii,:)) mean(seg(ii,:))])+levels(ii);
%     plot(x,y,'--','linewidth',1,'Color',clr(ii,:))
end
% cell locations on the attractors
smb = 'o*d^v<>phx+s';
for ii_cell = 1:length(cell_examples_IX)
    cell_num = cell_examples_IX(ii_cell);
    [r,c]=find(ind==cell_num);
    for ii_net = 1:length(c)
        net = c(ii_net);
        bn = th(r(ii_net),net);
        loc = bn * bin_size;
        h=plot(loc,levels(net),smb(ii_cell),'Color',clr(ii_cell,:));
        h.MarkerFaceColor = clr(ii_cell,:);
    end
end
% connectivity of a single cell
cell_num = connectivity_cell_IX;
[r,c]=find(ind==cell_num);
for ii_net = 1:length(r)
    net = c(ii_net);
    bn = th(r(ii_net),net);
    loc = bn * bin_size;
    rad = 0.05 * range(seg(net,:));
    %
    x = linspace(loc,loc+rad,100);
    t = linspace(0,pi,length(x));
    y = sin(t);
    y = y.*0.2;
    y = y + levels(net);
    plot(x,y,'-','Color','r','LineWidth',1.5);
    %
    x = linspace(loc,loc-rad,100);
    t = linspace(0,pi,length(x));
    y = sin(t);
    y = y.*0.2;
    y = y + levels(net);
    plot(x,y,'-','Color','r','LineWidth',1.5);
    %
    x = linspace(loc,loc+2*rad,100);
    t = linspace(0,pi,length(x));
%     t = linspace(0,sqrt(pi),length(x));
%     t = flip(t);
%     t = t.^2;
    y = sin(t);
    y = y.*0.5;
    y = y + levels(net);
    plot(x,y,'-','Color','r','LineWidth',0.5);
    %
    x = linspace(loc,loc-2*rad,100);
    t = linspace(0,pi,length(x));
%     t = linspace(0,sqrt(pi),length(x));
%     t = flip(t);
%     t = t.^2;
    y = sin(t);
    y = y.*0.5;
    y = y + levels(net);
    plot(x,y,'-','Color','r','LineWidth',0.5);
end
hax=gca;
hax.YLim = [0 3.5];
hax.Visible = 'off';


%% cell examples 
locs = linspace(0,200,size(m,1));
if multiple_panels

    for ii_cell = 1:length(panel_B)
        axes(panel_B(ii_cell));
        cla('reset')
        hold on
        run=10;
        plot(locs, m(:,cell_examples_IX(ii_cell),run),'LineWidth',2,'Color', clr(ii_cell,:));
        hax=gca;
        hax.XTick=[];
        hax.YTick=[];
    end
    axes(panel_B(end));
    text(-0.1,1.05, 'B', 'Units','normalized','FontWeight','bold');
    axes(panel_B(1));
    xlabel('Position (m)');
    xlim([0 200]);
    xticks([0:50:200]);
    axes(panel_B(3));
    ylabel('Firing rate (a.u.)');


else % single panel

    axes(panel_B);
    cla('reset')
    hold on
    text(-0.1,1.05, 'B', 'Units','normalized','FontWeight','bold');
    run=10;
    hl=plot(locs,m(:,cell_examples_IX,run),'LineWidth',2);
    for ii_hl = 1:length(hl)
        hl(ii_hl).Color = clr(ii_hl,:);
    end
    xlabel('Position (m)');
    xlim([0 200]);
    xticks([0:50:200]);
    ylabel('Firing rate (a.u.)');
    hax=gca;
    hax.XTick=[];
    hax.YTick=[];


end

%%
axes(panel_C(1));
cla('reset')
hold on
text(-0.45,1.2, 'C', 'Units','normalized','FontWeight','bold');
h=histogram(res.fnumber);
h.Normalization = 'pdf';
h.FaceColor = 0.5*[1 1 1];
xlabel({'No. of fields per direction'},'Units','normalized','Position',[0.5 -0.18]);
ylabel('PDF')
hax=gca;
hax.YScale = 'log';
hax.YLim = 10.^[-3 0];
% hax.XTick=[];
hax.YTick=10.^[-3:0];

%%
axes(panel_C(2));
cla('reset')
hold on
% text(-0.25,1.2, 'D', 'Units','normalized','FontWeight','bold');

h=histogram(res.fsizes_all.*bin_size);
h.Normalization = 'pdf';
h.FaceColor = 0.5*[1 1 1];
h.EdgeColor = 'none';

h=histogram(res.fsizes_all.*bin_size);
h.Normalization = 'pdf';
h.DisplayStyle = 'stairs';
h.EdgeColor = 'k';
xlabel('Field size (cm)','Units','normalized','Position',[0.5 -0.18]);
ylabel('PDF')
hax=gca;
hax.YScale = 'log';
hax.XTick=[0 20 40];
% hax.YTick=[];
hax.XLim = [0 45];
hax.YLim = [1e-3 2e-1];

%%
axes(panel_C(3));
cla('reset')
hold on
% text(-0.25,1.2, 'E', 'Units','normalized','FontWeight','bold');
h=histogram(res.sratio);
h.Normalization = 'pdf';
h.FaceColor = 0.5*[1 1 1];
xlabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[0.5 -0.17]);
ylabel('PDF')
xlim([0 20])
hax=gca;
hax.XScale = 'log';
hax.YScale = 'log';
% hax.XTick=[];
% hax.YTick=[];
hax.XLim = [0 20];
hax.YLim = [1e-3 1e0];


%% print/save the figure
fig_name_out = fullfile(res_dir, [fig_name_str]);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');


%%





