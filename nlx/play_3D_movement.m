function play_3D_movement(ts,pos,speed,chunk)
% ts        1XN vector, in seconds
% pos       3xN matix, in meters
% speed     use 1 to play in orignial speed
% chunk     number of points to plot together (should be used when setting
%           speed to high values.

%%
N = length(ts);
if chunk==0
    chunk = N;
end

%% create figure
figure
hold on
colorbar
colormap cool
axis equal
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
xlimits = [min(pos(1,:))  max(pos(1,:))];
ylimits = [min(pos(2,:))  max(pos(2,:))];
zlimits = [min(pos(3,:))  max(pos(3,:))];
margins = 0.1;
xlimits = xlimits + margins*[-1 +1]*range(diff(xlimits));
ylimits = ylimits + margins*[-1 +1]*range(diff(ylimits));
zlimits = zlimits + margins*[-1 +1]*range(diff(zlimits));
xlim(xlimits)
ylim(ylimits)
zlim(zlimits)

%% plot!
pause_time = [0 diff(ts)];
pause_time = pause_time./speed;
S = 10.*ones(1,N);
C = ts;
for ii_point = 1:chunk:N
    IX_points_chunk_start = ii_point;
    IX_points_chunk_end = min((ii_point+chunk-1),N);
    IX_points = IX_points_chunk_start:IX_points_chunk_end;
    scatter3(pos(1,IX_points), pos(2,IX_points), pos(3,IX_points), S(IX_points), C(IX_points));
    refreshdata;
    drawnow; 
    chunk_pause_time = sum(pause_time(IX_points));
    pause(chunk_pause_time);
end








%%