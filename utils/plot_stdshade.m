function plot_stdshade(x_data,y_data,std_values,alpha,acolor)
% plot mean and sem/std is shown as shading.
% - acolor defines the used color (default is red) 
% - alpha defines transparency of the shading (default is no shading and black mean line)

fill([x_data fliplr(x_data)],[y_data+std_values fliplr(y_data-std_values)],...
    acolor, 'FaceAlpha', alpha, 'EdgeAlpha', alpha, 'linestyle','none');
hold on;
plot(x_data, y_data, 'color', acolor, 'linewidth' , 1.5); %% change color or linewidth to adjust mean line

end

