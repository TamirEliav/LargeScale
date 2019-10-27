%%
clear
clc

%%
res_dir = 'L:\Analysis\Results\posters_presentations\SfN_2019';
mkdir(res_dir)
fig_name_str = 'SfN_2019_Behavior_constant_speed_flights';

%% create figure
figure_size_cm = [21.6 27.9]; % ~US letter
figure ;
% Some WYSIWYG options:
set(gcf,'DefaultAxesFontSize',8);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'DefaultAxesUnits','centimeters');
set(gcf,'PaperType','usletter')
% set(gcf,'PaperType','<custom>');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 figure_size_cm]);
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]); % position on screen...
set(gcf, 'Renderer', 'painters');
set(groot, 'defaultAxesTickDir', 'out');
annotation('textbox', [0.5 1 0 0], 'String',fig_name_str, 'HorizontalAlignment','center','Interpreter','none');

% panel_A = axes('position', [5  5 10  3.5]);
panel_B = axes('position', [5 10 5.371 6.288]);

%% choose data
prm = PARAMS_GetAll();
panel_J_opt = 11;
panel_J_data_options = {
    'b0034_d180313';
    'b0034_d180315';
    'b0079_d161003';
    'b0148_d170613';
    'b0148_d170615';
    'b0148_d170626';
    'b0148_d170802';
    'b0148_d170807';
    'b2289_d180520';    
    'b2289_d180615';
    'b0079_d160928';
};
exp_ID = panel_J_data_options{panel_J_opt};
exp = exp_load_data(exp_ID,'details','flight');
FE=exp.flight.FE;
FE([exp.flight.FE.distance]<prm.flight.full_min_distance) = [];

annotation('textbox', [0.5 0.9 0 0], 'String',exp_ID, 'HorizontalAlignment','center','Interpreter','none');

%% speed
% axes(panel_A);
% cla
% hold on
% x = [FE.pos];
% y = [FE.vel];
% y = abs(y);
% ylimits = [0 9];
% baseval = 0.1;
% area([1 prm.fields.valid_speed_pos(1)]    , ylimits([2 2]), baseval, 'FaceColor',0.8*[1 1 1],'EdgeColor','none','ShowBaseLine','off');
% area([  prm.fields.valid_speed_pos(2) 199], ylimits([2 2]), baseval, 'FaceColor',0.8*[1 1 1],'EdgeColor','none','ShowBaseLine','off');
% plot(x,y,'.','Color', 'k','MarkerSize',1);
% set(gca,'xtick',0:50:200,'ytick',[-10 0 8],'xlim',[0 200])
% set(gca,'tickdir','out','TickLength',repelem(0.01,2));
% ylim(ylimits);
% xlabel('Position (m)','Units','normalized','Position',[0.5 -0.25]);
% ylabel('Speed (m/s)','Units','normalized','Position',[-0.05 0.45]);
% ha=gca;
% ha.XRuler.TickLabelGapMultiplier = -0.3;
% ha.YRuler.TickLabelGapMultiplier = 0.001;


%% flights
axes(panel_B);
cla
hold on
switch 3
    case 1
        x = [FE.pos];
        y = [FE.ts];
        subsample = 10;
        x = x(1:subsample:end);
        y = y(1:subsample:end);
        plot(x,y,'.k','MarkerSize',1);
        ti = exp_get_sessions_ti(exp.details.exp_ID, 'Behave');
        t0=ti(1);
        rescale_plot_data('y',[1e-6/60,t0]);
        xlabel('Position (m)','Units','normalized','Position',[0.5 -0.05]);
        ylabel('Time (min)','Units','normalized','Position',[-0.08 0.45]);
    case 2
        x = [FE.pos];
        [FE.number2] = disperse(1:length(FE));
        y = arrayfun(@(FE)(FE.number2*ones(size(FE.pos))),FE,'UniformOutput',0);
        y = [y{:}];
        subsample = 10;
        x = x(1:subsample:end);
        y = y(1:subsample:end);
        plot(x,y,'.k','MarkerSize',1);
        xlabel('Position (m)','Units','normalized','Position',[0.5 -0.05]);
        ylabel('Flight no.','Units','normalized','Position',[-0.08 0.45]);
    case 3
        [FE.number2] = disperse(1:length(FE));
        x = arrayfun(@(FE)([FE.pos([1 end])]),FE,'UniformOutput',0);
        x = cat(1,x{:});
        y = arrayfun(@(FE)(FE.number2*[1 1]),FE,'UniformOutput',0);
        y = cat(1,y{:});
        plot(x',y','-k');
        xlabel('Position (m)','Units','normalized','Position',[0.5 -0.05]);
        ylabel('Flight no.','Units','normalized','Position',[-0.08 0.45]);
end
ha=gca;
ha.XTick = 0:50:200;
ha.XLim = [0 200];
ha.YLim(1) = 0;
% ha.YLim(2) = 100;
ha.TickDir = 'out';
ha.TickLength = repelem(0.01,2);
ha.XRuler.TickLabelGapMultiplier = -0.3;
ha.YRuler.TickLabelGapMultiplier = 0.001;

%% save figure
fig_name_str = fig_name_str+"_opt"+panel_J_opt+"_"+exp_ID;
% fig_name_str = fig_name_str + "_subsample_"+subsample;
fig_name_out = fullfile(res_dir, fig_name_str);
print(gcf, fig_name_out, '-dpdf', '-cmyk', '-painters');
% print(gcf, fig_name_out, '-dtiff', '-cmyk', '-painters');
% saveas(gcf , fig_name_out, 'fig');
disp('figure was successfully saved to pdf/tiff/fig formats');

% close(gcf)








%%
