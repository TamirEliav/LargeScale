function decoding_plot_session_seqs(exp_ID, epoch_type, params_opt, event_type)
switch epoch_type
    case 'rest'
        plot_rest_session(exp_ID, epoch_type, params_opt, event_type);
    case 'sleep'
        plot_sleep_session(exp_ID, epoch_type, params_opt, event_type);
end
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function plot_rest_session(exp_ID, epoch_type, params_opt, event_type)
%% load data
% exp = exp_load_data(exp_ID, 'details','rest','pos','LM','MUA_FR_map','ripples','MUA');
exp = exp_load_data(exp_ID, 'details','rest','pos','LM','MUA_FR_map');
events= decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
seqs = [events.seq_model];
[seqs, TF] = decoding_apply_seq_inclusion_criteria(seqs);
events(~TF) = [];

%% epochs ts/name
epochs_ti = exp.rest.ti;
epochs_rest_ball_num = [exp.rest.events.ball_num];
n_epochs = size(epochs_ti,1);
epochs_durations = range(epochs_ti,2) * 1e-6; % in seconds

%% arrange data for plotting
xlimits = exp_get_sessions_ti(exp_ID,'Behave','Behave_6m');
xlimits = [min(xlimits,[],'all') max(xlimits,[],'all')];
IX = get_data_in_ti(exp.pos.proc_1D.ts,xlimits);
t = exp.pos.proc_1D.ts(IX);
pos = exp.pos.proc_1D.pos(IX);
% nanpos = isnan(pos);
% pos(nanpos) = interp1(t(~nanpos), pos(~nanpos), t(nanpos),'linear','extrap');
ds = 10;
t = t(1:ds:end);
pos = pos(1:ds:end);

%%
fig=figure;
fig.WindowState = 'Maximized';
ax=[];

% ax(1)=axes('Units','normalized','Position',[.03 .86 .85 .08]);
% hold on
% ripples_t_rest_TF = any(exp.ripples.t>epochs_ti(:,1) & exp.ripples.t<epochs_ti(:,2),1);
% mua_t_rest_TF = any(exp.MUA.t>epochs_ti(:,1) & exp.MUA.t<epochs_ti(:,2),1);
% % ripples_t_rest_TF = get_data_in_ti(exp.ripples.t,epochs_ti);
% % mua_t_rest_IX = get_data_in_ti(exp.ripples.t,epochs_ti);
% exp.ripples.zpripple(~ripples_t_rest_TF)=nan;
% exp.MUA.zFR(~mua_t_rest_TF)=nan;
% plot(exp.ripples.t, exp.ripples.zpripple);
% plot(exp.MUA.t, exp.MUA.zFR);
% ylim([-1 5])
% xlim(xlimits)
% rescale_plot_data('x',[1e-6 xlimits(1)]);

ax(1)=axes('Units','normalized','Position',[.03 .05 .85 .85]);
hold on
plot(t,pos, 'Color',0.7*[1 1 1])
hax=gca;
clrs = hax.ColorOrder([1 4],:);
for ii_seq = 1:length(seqs)
    seq = seqs(ii_seq);
    x = [seq.start_ts seq.end_ts];
    y = [seq.start_pos seq.end_pos];
    c = 'r';
   switch seq.state_direction
       case 1
           c = clrs(2,:);
       case -1
           c = clrs(1,:);
   end
   if seq.forward
       ls = '-';
   else
       ls = ':';
   end
    plot(x,y,'Color',c,'LineWidth',1.5,'LineStyle',ls);
    plot(x(1),y(1),'o','Color',c,'MarkerSize',4,'MarkerFaceColor',c);
    plot(x(2),y(2),'o','Color',c,'MarkerSize',4);
end
% plot cross-overs
if isfield(exp.pos.proc_1D,'co')
    x = exp.pos.proc_1D.co.ts;
    y = exp.pos.proc_1D.co.pos;
    c = clrs(exp.pos.proc_1D.co.direction,:);
    scatter(x,y,50,c,'x');
end
xlim(xlimits)
ylimits = hax.YLim;
rests_ball_num = [exp.rest.events.ball_num];
for ii_epoch = 1:n_epochs
    switch rests_ball_num(ii_epoch)
        case 1
            c = 'g';
        case 2
            c = 'y';
    end
    area(epochs_ti(ii_epoch,:),ylimits([2 2]),'FaceAlpha',.05,'FaceColor',c,'EdgeColor','none');
end
arrayfun(@(y,str)(yline(y,'-',str,'Color',0.5.*[1 1 1],'LineWidth',0.5,'LabelVerticalAlignment','middle')),[exp.LM.pos_proj], string({exp.LM.name}))
xlabel('Time (s)')
ylabel('Position (m)')
rescale_plot_data('x',[1e-6 xlimits(1)]);

% marginal
ax(2)=axes('Units','normalized','Position',[.9 .05 .09 .85]);
hold on
seqs_edges = [seqs.start_pos; seqs.end_pos];
seqs_edges = [min(seqs_edges ); max(seqs_edges )]';
x = linspace(exp.rest.balls_loc(1),exp.rest.balls_loc(2),200);
[~,~,~,x_per_seq] = get_data_in_ti(x,seqs_edges);
clear hh
hh(1)=histogram([x_per_seq{[seqs.state_direction]==1}],x,'Orientation','horizontal','DisplayStyle','stairs','normalization','count','LineWidth',2,'EdgeColor',clrs(2,:));
hold on
hh(2)=histogram([x_per_seq{[seqs.state_direction]==-1}],x,'Orientation','horizontal','DisplayStyle','stairs','normalization','count','LineWidth',2,'EdgeColor',clrs(1,:));
hax=gca;
m = hax.XLim(2) / max(exp.MUA_FR_map.maps,[],"all");
plot(exp.MUA_FR_map.maps(1,:).*m, exp.MUA_FR_map.bin_centers, ':', 'Color',clrs(2,:),'LineWidth',1);
plot(exp.MUA_FR_map.maps(2,:).*m, exp.MUA_FR_map.bin_centers, ':', 'Color',clrs(1,:),'LineWidth',1);
arrayfun(@(y,str)(yline(y,'-',str,'Color',0.5.*[1 1 1],'LineWidth',0.5,'LabelVerticalAlignment','middle')),[exp.LM.pos_proj], string({exp.LM.name}))
xlabel('Counts')
ylabel('Position (m)')
linkaxes(ax([1 2]),'y');
ylim(exp.rest.balls_loc+3.*[-1 1])

% add legend
axes('Units','normalized','Position',[.8 .935 .01 .05]);
hold on
plot([0 0],[0 1],'Color',clrs(1,:),'LineWidth',2)
plot([1 1],[0 1],'Color',clrs(2,:),'LineWidth',2)
plot(0,0,'v','Color',clrs(1,:),'MarkerFaceColor',clrs(1,:))
plot(1,1,'^','Color',clrs(2,:),'MarkerFaceColor',clrs(2,:))
axis off

axes('Units','normalized','Position',[.82 .935 .02 .05]);
hold on
plot([0 1],[0 0],'k:','LineWidth',2)
plot([0 1],[1 1],'k-','LineWidth',2)
plot(1,2,'ko')
plot(1,3,'ko','MarkerFaceColor','k')
text(1.3,0,'Reverse','VerticalAlignment','middle')
text(1.3,1,'Forward','VerticalAlignment','middle')
text(1.3,2,'End','VerticalAlignment','middle')
text(1.3,3,'Start','VerticalAlignment','middle')
axis off

axes('Units','normalized','Position',[.91 .93 .01 .02]);
hold on
plot([0 1],[1 1],'-k','LineWidth',2)
plot([0 1],[0 0],':k','LineWidth',1)
text(1.3,1,'replay (counts)','VerticalAlignment','middle');
text(1.3,0,'MUA FR map (norm.)','VerticalAlignment','middle');
axis off

h=sgtitle({sprintf('%s session replay sequences - %s',epoch_type,exp_ID)},'interpreter','none');

%% save fig
fig_name = sprintf('%s_session_seqs_%s_opt_%d',exp_ID, epoch_type, params_opt);
dir_OUT = 'F:\sequences\decoded_figs\seqs_entire_session';
mkdir(dir_OUT);
filename = fullfile(dir_OUT, fig_name);
saveas(fig, filename , 'fig');
saveas(fig, filename , 'jpg');

end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function plot_sleep_session(exp_ID, epoch_type, params_opt, event_type)
%% load data
exp = exp_load_data(exp_ID, 'details','LM','rest','MUA_FR_map');
events= decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
seqs = [events.seq_model];
[seqs, TF] = decoding_apply_seq_inclusion_criteria(seqs);
events(~TF) = [];

%% epochs ts/name
sleep_session_TF = contains(exp.details.session_names, 'sleep','IgnoreCase',true);
epochs_ti = exp.details.session_ts(sleep_session_TF,:);
epochs_names = exp.details.session_names(sleep_session_TF);
n_epochs = size(epochs_ti,1);
epochs_durations = range(epochs_ti,2) * 1e-6; % in seconds

%%
fig=figure;
fig.WindowState = 'Maximized';
pnl = panel();
pnl.pack('h',[.9 .1],'v',n_epochs)
pnl.margintop = 20;
pnl.de.margintop = 5;
for ii_epoch = 1:n_epochs
    seqs_epoch_IX = [events.epoch_num] == ii_epoch;
    seqs_epoch = seqs(seqs_epoch_IX);
    pnl(1,ii_epoch).select()
    hold on
    hax=gca;
    clrs = hax.ColorOrder([1 4],:);
    for ii_seq = 1:length(seqs_epoch)
        seq = seqs_epoch(ii_seq);
        x = [seq.start_ts seq.end_ts];
        y = [seq.start_pos seq.end_pos];
        c = 'r';
       switch seq.state_direction
           case 1
               c = clrs(2,:);
           case -1
               c = clrs(1,:);
       end
       if seq.forward
           ls = '-';
       else
           ls = ':';
       end
        plot(x,y,'Color',c,'LineWidth',1.5,'LineStyle',ls);
        plot(x(1),y(1),'o','Color',c,'MarkerSize',4,'MarkerFaceColor',c);
        plot(x(2),y(2),'o','Color',c,'MarkerSize',4);
    end
    xlimits = epochs_ti(ii_epoch,:);
    xlim(xlimits)
    arrayfun(@(y,str)(yline(y,'-',str,'Color',0.5.*[1 1 1],'LineWidth',0.5,'LabelVerticalAlignment','middle')),[exp.LM.pos_proj], string({exp.LM.name}))
    xlabel('Time (s)')
    ylabel('Position (m)')
    title(epochs_names{ii_epoch},'FontWeight','bold');
    rescale_plot_data('x',[1e-6 xlimits(1)]);

    % marginal
    pnl(2,ii_epoch).select()
    seqs_edges = [seqs_epoch.start_pos; seqs_epoch.end_pos];
    seqs_edges = [min(seqs_edges ); max(seqs_edges )]';
    x = linspace(exp.rest.balls_loc(1),exp.rest.balls_loc(2),200);
    [~,~,~,x_per_seq] = get_data_in_ti(x,seqs_edges);
    histogram([x_per_seq{[seqs_epoch.state_direction]==1}],x,'Orientation','horizontal','DisplayStyle','stairs','normalization','count','LineWidth',2,'EdgeColor',clrs(2,:));
    hold on
    histogram([x_per_seq{[seqs_epoch.state_direction]==-1}],x,'Orientation','horizontal','DisplayStyle','stairs','normalization','count','LineWidth',2,'EdgeColor',clrs(1,:));
    hax=gca;
    m = hax.XLim(2) / max(exp.MUA_FR_map.maps,[],"all");
    plot(exp.MUA_FR_map.maps(1,:).*m, exp.MUA_FR_map.bin_centers, ':', 'Color',clrs(2,:),'LineWidth',1);
    plot(exp.MUA_FR_map.maps(2,:).*m, exp.MUA_FR_map.bin_centers, ':', 'Color',clrs(1,:),'LineWidth',1);
    arrayfun(@(y,str)(yline(y,'-',str,'Color',0.5.*[1 1 1],'LineWidth',0.5,'LabelVerticalAlignment','middle')),[exp.LM.pos_proj], string({exp.LM.name}))
    xlabel('Counts')
    ylabel('Position (m)')
end
linkaxes(pnl.de.axis,'y');
ylim(exp.rest.balls_loc+3.*[-1 1])

% add legend
axes('Units','normalized','Position',[.8 .935 .01 .05]);
hold on
plot([0 0],[0 1],'Color',clrs(1,:),'LineWidth',2)
plot([1 1],[0 1],'Color',clrs(2,:),'LineWidth',2)
plot(0,0,'v','Color',clrs(1,:),'MarkerFaceColor',clrs(1,:))
plot(1,1,'^','Color',clrs(2,:),'MarkerFaceColor',clrs(2,:))
axis off

axes('Units','normalized','Position',[.82 .935 .02 .05]);
hold on
plot([0 1],[0 0],'k:','LineWidth',2)
plot([0 1],[1 1],'k-','LineWidth',2)
plot(1,2,'ko')
plot(1,3,'ko','MarkerFaceColor','k')
text(1.3,0,'Reverse','VerticalAlignment','middle')
text(1.3,1,'Forward','VerticalAlignment','middle')
text(1.3,2,'End','VerticalAlignment','middle')
text(1.3,3,'Start','VerticalAlignment','middle')
axis off

axes('Units','normalized','Position',[.91 .93 .01 .02]);
hold on
plot([0 1],[1 1],'-k','LineWidth',2)
plot([0 1],[0 0],':k','LineWidth',1)
text(1.3,1,'replay (counts)','VerticalAlignment','middle');
text(1.3,0,'MUA FR map (norm.)','VerticalAlignment','middle');
axis off

h= pnl.title( sprintf('%s session replay sequences - %s',epoch_type,exp_ID) );
h.Interpreter = 'none';
h.Position = [.5 1.05];
h.FontSize = 14;
h.FontWeight = 'bold';
% h=sgtitle({sprintf('%s session replay sequences - %s',epoch_type,exp_ID)},'interpreter','none');

%% save fig
fig_name = sprintf('%s_session_seqs_%s_opt_%d',exp_ID, epoch_type, params_opt);
dir_OUT = 'F:\sequences\decoded_figs\seqs_entire_session';
mkdir(dir_OUT);
filename = fullfile(dir_OUT, fig_name);
saveas(fig, filename , 'fig');
saveas(fig, filename , 'jpg');

end

