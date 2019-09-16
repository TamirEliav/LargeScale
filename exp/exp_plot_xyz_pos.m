function exp_plot_xyz_pos(exp_ID)

%%
% clear
% clc
% exp_ID = 'b0034_d180313';

%% load exp data
exp = exp_load_data(exp_ID,'details','pos','flight');
prm = PARAMS_GetAll();

%% take only full flights
FE = exp.flight.FE([exp.flight.FE.distance]>prm.flight.full_min_distance);

%% calc Y/Z for FE data
for ii_fe = 1:length(FE)
    FE(ii_fe).Y = interp1(exp.pos.proj.ts, exp.pos.proj.pos(:,2), FE(ii_fe).ts);
    FE(ii_fe).Z = interp1(exp.pos.raw.ts_nlg_usec, exp.pos.raw.pos(:,3), FE(ii_fe).ts);
end

%% calc Y/Z over the tunnel locations
X = 1:200;
YZ_by_dir = {};
directions = [1 -1];
YZ_std = nan(2,2,length(X));
for ii_dir = 1:length(directions)
    dir = directions(ii_dir);
    FE_dir = FE([FE.direction]==dir);
    YZ = nan(length(FE_dir),2,length(X));
    for ii_flight = 1:length(FE_dir)
        flight = FE_dir(ii_flight);
        [C,ia,idx] = unique(flight.pos,'stable');
        y = accumarray(idx,flight.Y,[],@median);
        z = accumarray(idx,flight.Z,[],@median);
        Y = interp1(C,y,X,'linear');
        Z = interp1(C,z,X,'linear');
        YZ(ii_flight,1,:) = Y;
        YZ(ii_flight,2,:) = Z;
        YZ_by_dir{ii_dir} = YZ;
        YZ_std(ii_dir, :,:) = nanstd(YZ);
    end
end

%% find lowest variability in both Y&Z
invalid_pos_IX = ~(X>prm.fields.valid_speed_pos(1) & X<prm.fields.valid_speed_pos(end));
sdf = YZ_std;
sdf(:,:,invalid_pos_IX)  = nan;
sdf = sum(sdf,2);
[~,low_std_IX] = min(sdf,[],3);

%%
figure('Units','normalized','Position',[0 0 1 1]);
pnl=panel();
pnl.pack(3,2);
pnl.margin = [20 30 20 20];
pnl.de.margin = 30;
h=pnl.title(exp.details.exp_ID);
h.Position = [0.5 1.03];
h.FontSize = 16;
h.Interpreter = 'none';
for ii_dir = 1:length(directions)
    
    dir = directions(ii_dir);
    FE_dir = FE([FE.direction]==dir);
    
    pnl(1,ii_dir).select();
    plot([FE_dir.pos],[FE_dir.Y], '.', 'Color', prm.graphics.colors.flight_directions{ii_dir});
    xlabel('X position (m)')
    ylabel('Y position (m)')
    title("dir "+ii_dir)
    text(0.99,1,"n="+length(FE_dir),'Units','normalized','HorizontalAlignment','right','FontSize',14);
    xline(X(low_std_IX(ii_dir)),'k',X(low_std_IX(ii_dir)),'LabelHorizontalAlignment','center','LabelVerticalAlignment','top');
    
    pnl(2,ii_dir).select();
    plot([FE_dir.pos],[FE_dir.Z], '.', 'Color', prm.graphics.colors.flight_directions{ii_dir});
    xlabel('X position (m)')
    ylabel('Z position (m)')
    title("dir "+ii_dir);
    text(0.99,1,"n="+length(FE_dir),'Units','normalized','HorizontalAlignment','right','FontSize',14);
    xline(X(low_std_IX(ii_dir)),'k',X(low_std_IX(ii_dir)),'LabelHorizontalAlignment','center','LabelVerticalAlignment','top');
    
    pnl(3,ii_dir).select();
    yyaxis left
    plot(X,squeeze(YZ_std(ii_dir,1,:)),'k')
    set(gca,'YColor', 'k');
    ylabel('Y std (m)')
    yyaxis right
    plot(X,squeeze(YZ_std(ii_dir,2,:)),'m')
    set(gca,'YColor', 'm');
    ylimits = get(gca,'ylim');
    ylimits(2) = min(1,ylimits(2));
    set(gca,'ylim',ylimits);
    ylabel('Z std (m)')
    xlabel('X position (m)')
    title("dir "+ii_dir);
    text(0.99,1,"n="+length(FE_dir),'Units','normalized','HorizontalAlignment','right','FontSize',14);
    xline(X(low_std_IX(ii_dir)),'k',X(low_std_IX(ii_dir)),'LabelHorizontalAlignment','center','LabelVerticalAlignment','top');
end
% h=findobj(gcf,'type','axe');
linkaxes(pnl.de.axis,'x')

%% save figure
fig_filename = fullfile('L:\Analysis\Results\exp\pos_XYZ', [exp_ID '_pos_XYZ']);
saveas(gcf,fig_filename,'tif')
% saveas(gcf,fig_filename,'fig')

%% save data
XYZ = struct();
XYZ.YZ_by_dir = YZ_by_dir;
XYZ.YZ_std = YZ_std;
XYZ.low_std_IX = low_std_IX;
file_name = fullfile('L:\Analysis\Results\exp\pos_XYZ', [exp_ID '_pos_XYZ']);
save(file_name,'XYZ');


%%



