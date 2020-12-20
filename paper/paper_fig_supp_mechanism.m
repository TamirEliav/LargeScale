%% Large Scale - Fig. S20 - Theoretical analysis (mechanism)

%%
clear 
clc

%%
paper_fig_8_arrange_sim_data
data = paper_fig_8_arrange_real_data();

%% define output files
res_dir = 'L:\paper_figures';
mkdir(res_dir)
fig_name_str = 'fig_S20';
fig_caption_str = 'Theoretical analysis - Mechanism';
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
disp([fig_name_str ':' fig_caption_str]);
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
panel_A_size = [2.5 2.5];
panel_D_size = [2.5 3];
panel_BC_size = [3 3];
panel_A(1) = axes('position', [ 2.0 22 panel_A_size]);
panel_A(2) = axes('position', [ 5.5 22 panel_A_size]);
panel_A(3) = axes('position', [ 9.0 22 panel_A_size]);
panel_A(4) = axes('position', [12.5 22 panel_A_size]);
panel_A(5) = axes('position', [16.0 22 panel_A_size]);
panel_BC(1) = axes('position', [ 2.0 11.8+4.5 panel_BC_size]);
panel_BC(2) = axes('position', [ 2.0  7.0+4.5 panel_BC_size]);
panel_D(1) = axes('position', [ 2.0 17.5-11 panel_D_size]);
panel_D(2) = axes('position', [ 5.7 17.5-11 panel_D_size]);
panel_D(3) = axes('position', [ 9.0 17.5-11 panel_D_size]);
panel_D(4) = axes('position', [ 12 19.5-11 0.5 0.65]);
panel_E(1) = axes('position', [ 7.0 11.8+4.5 5 3]);
panel_E(2) = axes('position', [ 9.2 13.1+4.5 2.5 1.7]);
panel_F(1) = axes('position', [ 7.0 7+4 4 3]);
panel_G(1) = axes('position', [ 7.0 1.7 5 3]);
for ii_bat = 1:5
    panel_size = [5 2.5];
%     panel_size_inset = [2.5 1.25];
    panel_size_inset = [1.4 1.25];
%     inset_offset = [2.6 1.5];
    inset_offset1 = [2.00 1.55];
    inset_offset2 = [3.60 1.55];
    spacing = 1;
    panel_pos_x = 14.2;
    panel_pos_y = 2.35 + (ii_bat-1)*(panel_size(2)+spacing);
    panel_pos = [panel_pos_x panel_pos_y];
%     panel_pos_inset = panel_pos + inset_offset;
    panel_pos_inset1 = panel_pos + inset_offset1;
    panel_pos_inset2 = panel_pos + inset_offset2;
    panel_H(ii_bat,1) = axes('position', [panel_pos         panel_size]);
    panel_H(ii_bat,2) = axes('position', [panel_pos_inset1   panel_size_inset]);
    panel_H(ii_bat,3) = axes('position', [panel_pos_inset2   panel_size_inset]);
end
panel_H = flipud(panel_H);

%%
colormap gray
lw = 2;

%% ========================== SUPP panels =================================
nplot1 = 10 ; % number of examples to plot
iang = 20 ;   % slice angle to plot in supplementary (iang = 20 ==> 19 degrees) 
i3   = 1 ;    % i3 = 1 means purely MEC input

%% 2D grid
f2 = zeros(L,L,naMEC) ;
for i = 1:naMEC
    f2(:,:,i) = reshape(f2MEC(xplot,aMEC(i)),L,L) ;
end

for ia = 1:naMEC
    axes(panel_A(ia)); cla; hold on;
    imagesc(xFF,xFF,1-f2(:,:,ia)) ; hold on ;
    h=gca;
    h.XLim=[0 50];
    h.YLim=[0 50];
    h.XTick = [0:10:50];
    h.YTick = [0:10:50];
    h.XAxis.TickLength(1) = 0.025;
    h.YAxis.TickLength(1) = 0.025;
    h.XRuler.TickLabelGapOffset = -1;
    h.YRuler.TickLabelGapOffset = -0.5;
    h.Box = 'on';
    xlabel('X (m)', 'Units','normalized','Position',[0.5 -0.2]);
    if ia==1
        ylabel('Y (m)', 'Units','normalized','Position',[-0.25 0.5]);
    end
    text(0.5,1.15, "Module "+ia, 'Units','normalized','HorizontalAlignment','center');
    for iang = [1 iang] 
        plot([0 xFF(end)*cos(thMEC2(iang))],3+[0 xFF(end)*sin(thMEC2(iang))],'-r','LineWidth',lw/5) ; hold on ;
        axis square xy ; xlim([0 50]) ; ylim([0 50]) ;
        
    end
end
axes(panel_A(1));
text(-0.4,1.25, 'A', 'Units','normalized','FontWeight','bold');


%% models stats - field number/size/ratio
lw=2;

axes(panel_D(1)); cla; hold on;
text(-0.55,1.15, 'D', 'Units','normalized','FontWeight','bold');
plot(xnField,squeeze(pnFieldXs(i3,:,:))','-','Color',clrM(3,:),'lineWidth',lw) ; hold on ;
pnFieldS(pnFieldS==0)=eps;
% pnFieldS(pnFieldS==eps)=0;
plot(xnField,pnFieldS,'-','Color',clrM(1,:),'lineWidth',lw) ; hold on ;
plot(xnField,pnFieldM,'-','Color',clrM(2,:),'lineWidth',lw) ; hold on ;
xlabel({'No. of fields per direction'},'Units','normalized','Position',[0.5 -0.18]);
ylabel({'Probability';'density function'},'Units','normalized','Position',[-0.3 0.5]);
hax=gca;
hax.YScale='log';
hax.XLim = [0 25];
hax.YLim = [1e-3 1.5e-1];
hax.XTick = 0:10:40;
hax.YTick = 10.^[-4:0];

axes(panel_D(2)); cla; hold on;
plot(xField,squeeze(pFieldSizeXs(i3,:,:))','-','Color',clrM(3,:),'lineWidth',lw) ; hold on ;
plot(xField,pFieldSizeS,'-','Color',clrM(1,:),'lineWidth',lw) ; hold on ;
plot(xField,pFieldSizeM,'-','Color',clrM(2,:),'lineWidth',lw) ; hold on ;
xlabel('Field size (m)','Units','normalized','Position',[0.5 -0.18]);
% ylabel('PDF','Units','normalized','Position',[-0.3 0.5]);
hax=gca;
hax.YScale='log';
hax.XLim = [0 30];
hax.YLim = [2e-4 0.2];
hax.XTick = 0:10:40;
hax.YTick = 10.^[-4:0];

axes(panel_D(3)); cla; hold on;
plot(xMaxMinRatio,squeeze(pMaxMinRatioXs(i3,:,:))','-','Color',clrM(3,:),'lineWidth',lw) ; hold on ;
plot(xMaxMinRatio,pMaxMinRatioS,'-','Color',clrM(1,:),'lineWidth',lw) ; hold on ;
plot(xMaxMinRatio,pMaxMinRatioM,'-','Color',clrM(2,:),'lineWidth',lw) ; hold on ;
xlabel({'Field size ratio';'largest/smallest'},'Units','normalized','Position',[0.5 -0.17]);
% ylabel('PDF','Units','normalized','Position',[-0.25 0.5]);
hax=gca;
hax.YScale='log';
hax.XLim = [0 30];
hax.YLim = [5e-3 0.12];
hax.XTick = 0:10:40;
hax.YTick = 10.^[-4:0];

% legend
axes(panel_D(4)); cla; hold on;
axis off
plot([0 1],[1 1],'-','Color',clrM(2,:),'lineWidth',lw) ; hold on ;
plot([0 1],[2 2],'-','Color',clrM(1,:),'lineWidth',lw) ; hold on ;
plot([0 1],[3 3],'-','Color',clrM(3,:),'lineWidth',lw) ; hold on ;
% text(1.2, 2, 'Single-field CA3', 'FontSize',7);
% text(1.2, 1, 'Multi-field CA3', 'FontSize',7);
% text(1.2, 3, 'MEC, all slice angles', 'FontSize',7);
text(-0.2, 2, 'Single-field CA3', 'FontSize',7,'HorizontalAlignment','right');
text(-0.2, 1, 'Multi-field CA3', 'FontSize',7,'HorizontalAlignment','right');
text(-0.2, 3, 'MEC, all slice angles', 'FontSize',7,'HorizontalAlignment','right');

%% MEC input/output cell examples
axes(panel_BC(1)); cla; hold on;
text(-0.4,1.25, 'B', 'Units','normalized','FontWeight','bold');
imagesc(xFF,1:nplot2,1-fP(:,round(linspace(1,Nmax,nplot2)),iang)') ; hold on ;
for i = 1:naMEC-1
    plot([0 L]*ds/100,[1 1]*i*nplot2/naMEC+0.5,'-','Color',[0.6 0.6 0.6]) ;
end
for i = 1:naMEC
    y = 2.5 + (i-1)*nplot2/naMEC;
    text(210, y, "M"+i, 'FontSize',7);
end
ylabel('Input neuron no.') ;
xlabel('Position (m)','Units','normalized','Position',[0.5 -0.05]) ;
% title('Periodic grid-cells MEC','FontWeight','bold','Units','normalized','Position',[0.5 1.05]) ;
% title('1D slice through 2D grid','FontWeight','bold','Units','normalized','Position',[0.5 1.05]) ;
title({'Input from MEC';'1D slice through 2D grid'},'FontWeight','bold','Units','normalized','Position',[0.5 1.05]) ;
h=gca;
h.Box = 'on';
h.XTick=[0 200];
h.YTick=[1 5 10 15 20];
h.XLim = [0 200];
h.YLim = [0.5 20.5];
h.XAxis.TickLength(1) = 0.03;
h.YRuler.TickLabelGapOffset = 0.5;

axes(panel_BC(2)); cla; hold on;
text(-0.4,1.15, 'C', 'Units','normalized','FontWeight','bold');
i3 = 1;
load(fullfile(sim_data_dir,['CA3MEC_FFModel_fResXw_iangle_' num2str(iang) '.mat']));
load(fullfile(sim_data_dir,['CA3MEC_FFModel_fResXs_iangle_' num2str(iang) '.mat']));
% find some periodic output examples
[~,IX] = sort(max(diff(abs(fft(fResXs(:,:,i3),[],1)).^2)),'descend');
IX = IX([9 7 51:58]);
IX([5 10]) = IX([10 5]);
imagesc(xFF,1:nplot1,1-fResXs(:,IX,i3)') ;
box on
ylabel('Output neuron no.') ;
xlabel('Position (m)','Units','normalized','Position',[0.5 -0.07]) ;
text(0.5,1.1, 'CA1 output', 'Units','normalized','FontWeight','bold','HorizontalAlignment','center','FontSize',8);
h=gca;
h.Box = 'on';
h.XTick=[0 200];
h.YTick=[1 5 10];
h.XLim = [0 200];
h.YLim = [0.5 10.5];
h.XAxis.TickLength(1) = 0.03;
h.YRuler.TickLabelGapOffset = 0.5;

%% add arrow from input -> output
har=annotation('arrow');
har.Units='centimeters';
har.Position = [3.5 15.6 0 -0.5];
har.HeadWidth = 5;
har.HeadLength = 5;
har.LineWidth = 1.5;

%% MEC/CA3 spectrum
axes(panel_E(1)); cla; hold on;
colorbar(panel_E(1), 'off');
text(-0.25,1.25, 'E', 'Units','normalized','FontWeight','bold');
lw = 1 ;
for i3 = 1:nrCA3
    plot(kPow,powOutXs(:,i3,iang),'Color',clrrCA3(i3,:),'LineWidth',lw) ; hold on ;
end
plot(kPow,powOutS,'Color',clrrCA3(end,:),'LineWidth',lw) ; hold on ;
plot(kMEC,1.25e-1,'^r','MarkerFaceColor','r','MarkerSize',4) ; hold on ;
hax=gca;
hax.XLim = [0 0.25];
hax.YTick = 10.^[-1 0 1 2];
hax.YScale = 'log';
hax.TickLength(1)=0.02;
xlabel('Spatial frequency (1/m)') ;
ylabel('Power (norm.)') ;
title('Model: Predicted spectrum of CA1 neurons','FontWeight','bold','Units','normalized','Position',[0.5 1.15]) ;
% zoom box
zoom_x = [0.03 0.09];
zoom_y = [6.5e-1 0.7e1];
plot(zoom_x,zoom_y([1 1]),'k-')
plot(zoom_x,zoom_y([2 2]),'k-')
plot(zoom_x([1 1]),zoom_y,'k-')
plot(zoom_x([2 2]),zoom_y,'k-')

% colorbar
colorbar_loc = hax.Position([1 2]) + hax.Position([3 4]).*[1 0.25];
colorbar_size = [0.2 0.5*hax.Position(4)];
colorbar_pos = [colorbar_loc colorbar_size];
hc=colorbar(panel_E(1),'manual');
hc.Units = 'centimeters';
hc.Position = colorbar_pos;
hc.Direction='reverse';
% hc.Ticks = [0 1];
% hc.TickLabels = {'MEC';'CA3'};
hc.Ticks = [];
text(1.02,0.808, 'MEC', 'Units','normalized','FontSize',6,'HorizontalAlignment','center');
text(1.02,0.2, 'CA3', 'Units','normalized','FontSize',6,'HorizontalAlignment','center');
colormap(panel_E(1),clrrCA3);

% zoom inset panel
axes(panel_E(2)); cla; hold on;
for i3 = 1:nrCA3
    plot(kPow,powOutXs(:,i3,iang),'Color',clrrCA3(i3,:),'LineWidth',lw) ; hold on ;
end
plot(kPow,powOutS,'Color',clrrCA3(end,:),'LineWidth',lw) ; hold on ;
hax=gca;
hax.XLim = zoom_x;
hax.YLim = zoom_y;
hax.YScale = 'log';
hax.XTick=[];
hax.YTick=[];
hax.Box='on';

%% spectrum at all grid slice angles (with 50-50 CA3-MEC input)
axes(panel_F(1)); cla; hold on;
text(-0.2,1.2, 'F', 'Units','normalized','FontWeight','bold');
i3 = 11 ; % color plot of spectrum at all slice angels shown for half/half MEC/CA3 input
xlimits = [0 0.2];
% IX = kPow < xlimits(2);
IX = kPow > xlimits(1) & kPow < xlimits(2);
M = squeeze(powOutXs(IX,i3,:))';
% M(:,1) =  nan;
% M(:,1) =  M(:,2);
% M = log10(M);
% imagesc(kPow,0:29, M, [-1 3.6]) ; axis xy ;  hold on ;
imagesc(kPow(IX),0:29, M) ; axis xy ;  hold on ;
plot(1./(aMEC/2.5),0,'^r','MarkerFaceColor','r','MarkerSize',5) ; hold on ;
hax=gca;
hax.ColorScale='log';
% hax.XLim = xlimits;
% hax.XLim = [min(kPow(IX)) max(kPow(IX))];
hax.XLim = [min(kPow(IX))/2 xlimits(2)];
xlabel('Spatial frequency (1/m)') ;
ylabel('Slice angle (\circ)') ;
title({'CA1 spectra for';'equal CA3 and MEC input'},'FontWeight','bold') ;

% colorbar
colorbar_loc = hax.Position([1 2]) + hax.Position([3 4]).*[1.07 0.25];
colorbar_size = [0.2 0.5*hax.Position(4)];
colorbar_pos = [colorbar_loc colorbar_size];
hc=colorbar(panel_F(1),'manual');
hc.Units = 'centimeters';
hc.Position = colorbar_pos;
% IX = kPow<xlimits(2);
% climits = [min(M(:,IX),[],'all') max(M(:,IX),[],'all')];
hc.TickDirection = 'out';
hc.TickLength = 0.08;
hc.TickLength = 0.1;
% hc.Ticks = -5:5;
% hc.Ticks = hc.Limits;
% hc.Ticks = 0:5:40;
hc.Ticks = 10.^[-5:5];
% hc.Ticks = [];
hc.TickLabels = [];
text(1.1,0.2, 'min', 'Units','normalized','FontSize',6,'HorizontalAlignment','center');
text(1.1,0.808, 'max', 'Units','normalized','FontSize',6,'HorizontalAlignment','center');
text(1.24,0.5, 'Power', 'Units','normalized','FontSize',7,'HorizontalAlignment','center','Rotation',-90);
colormap(panel_F(1),gray);

%% Relative MEC/CA3 input sparsity - Spectrum peakedness
axes(panel_G(1)); cla; hold on;
text(-0.25,1.1, 'G', 'Units','normalized','FontWeight','bold');
plot(1-[rCA3 1],rMaxdpowOutXs,'Color',[0.6 0.6 0.6]) ; hold on ;
plot(1-[rCA3 1],rMaxdpowOutXs(:,iang),'o-k','MarkerFaceColor',0.6*[1 1 1],'MarkerSize',4) ; hold on ;
plot(1-[rCA3 1],mean(rMaxdpowOutXs,2),'c','LineWidth',3) ; hold on ;
xlabel('Relative CA3 vs. MEC input');
ylabel('Spectral peak (norm.)');
text(0,-0.2, 'CA3', 'Units','normalized','FontSize',7,'HorizontalAlignment','center');
text(1,-0.2, 'MEC', 'Units','normalized','FontSize',7,'HorizontalAlignment','center');hax=gca;
hax.XRuler.TickLabelGapOffset = -1;
hax.TickLength(1) = 0.025;

%% Spectrum of real data maps (pooled over bats)
long_arm_clr = [0.4 0.5 0.6];
bats=unique(data.maps_bat);
for ii_bat = 1:length(bats)
    
    bat = bats(ii_bat);
    maps_IX = data.maps_bat==bat;
    
    axes(panel_H(ii_bat,1)); cla; hold on;
    if ii_bat==1
        text(-0.22,1.45, 'H', 'Units','normalized','FontWeight','bold');
        title({'Data: Measured spectrum of';'CA1 neurons in each bat'},'FontWeight','bold','Units','normalized','Position',[0.5 1.2]) ;
    end
    shadedErrorBar(data.freq, data.maps_spec(maps_IX,:), {@mean,@nansem});
    hax=gca;
    hax.YScale = 'log';
    hax.XLim = [0 0.25];
    hax.YTick = 10.^[-5:5];
    hax.TickLength(1) = 0.02;
    hax.XRuler.TickLabelGapOffset = -0;
    if ii_bat==length(bats)
        xlabel('Spatial frequency (1/m)') ;
    end
    ylabel('Power (norm.)')
    text(0.05,0.9, "Bat "+ii_bat, 'Units','normalized','FontSize',8,'FontWeight','bold',...
        'HorizontalAlignment','left','VerticalAlignment','middle');
    text(0.05,0.2, "n = "+sum(maps_IX), 'Units','normalized','FontSize',8,'FontWeight','normal',...
        'HorizontalAlignment','left','VerticalAlignment','middle');

    % inset1 - binarized maps
    axes(panel_H(ii_bat,2)); cla; hold on;
    shadedErrorBar(data.freq, data.maps01_spec(maps_IX,:), {@mean,@nansem});
    hax=gca;
    hax.YScale = 'log';
    hax.XLim = [0 0.25];
    hax.YTick = 10.^[-5:5];
    hax.YTickLabel = [];
    hax.TickLength(1) = 0.015 * 2;
    hax.XRuler.TickLabelGapOffset = -3;
    hax.YRuler.TickLabelGapOffset = -0.5;
    
    % inset2 - only long arm
    axes(panel_H(ii_bat,3)); cla; hold on;
    shadedErrorBar(data.freq_long, data.maps_long_spec(maps_IX,:), {@mean,@nansem});
    hax=gca;
    hax.YScale = 'log';
    hax.XLim = [0 0.25];
    hax.YTick = 10.^[-5:5];
    hax.YTickLabel = [];
    hax.TickLength(1) = 0.015 * 2;
    hax.XRuler.TickLabelGapOffset = -3;
    hax.YRuler.TickLabelGapOffset = -0.5;
end


%%

%% print/save the figure
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
disp('figure was successfully saved to pdf/tiff/fig formats');
