function plot_boxplot(data1,data2,ax,color1,color2)
box_width = 0.3; line_width = 1; wisker_width = 0.2;

median1 = median(data1);
median2 = median(data2);
quarntile1_25 = prctile(data1,25);
quarntile2_25 = prctile(data2,25);
quarntile1_75 = prctile(data1,75);
quarntile2_75 = prctile(data2,75);
prctile1_90 = prctile(data1,90);
prctile2_90 = prctile(data2,90);
prctile1_10 = prctile(data1,10);
prctile2_10 = prctile(data2,10);

line(ax,1+(box_width/2)*[-1,+1],median1*[1,1],'color',color1,'linewidth',line_width);
line(ax,1+(box_width/2)*[-1,+1],quarntile1_25*[1,1],'color',color1,'linewidth',line_width);
line(ax,1+(box_width/2)*[-1,+1],quarntile1_75*[1,1],'color',color1,'linewidth',line_width);
line(ax,1+(wisker_width/2)*[-1,+1],prctile1_90*[1,1],'color',color1,'linewidth',line_width);
line(ax,1+(wisker_width/2)*[-1,+1],prctile1_10*[1,1],'color',color1,'linewidth',line_width);
line(ax,[1,1],[quarntile1_75,prctile1_90],'color',color1,'linewidth',line_width);
line(ax,[1,1],[prctile1_10,quarntile1_25],'color',color1,'linewidth',line_width);
line(ax,(1-(box_width/2))*[1,1],[quarntile1_25,quarntile1_75],'color',color1,'linewidth',line_width);
line(ax,(1+(box_width/2))*[1,1],[quarntile1_25,quarntile1_75],'color',color1,'linewidth',line_width);

line(ax,2+(box_width/2)*[-1,+1],median2*[1,1],'color',color2,'linewidth',line_width);
line(ax,2+(box_width/2)*[-1,+1],quarntile2_25*[1,1],'color',color2,'linewidth',line_width);
line(ax,2+(box_width/2)*[-1,+1],quarntile2_75*[1,1],'color',color2,'linewidth',line_width);
line(ax,2+(wisker_width/2)*[-1,+1],prctile2_90*[1,1],'color',color2,'linewidth',line_width);
line(ax,2+(wisker_width/2)*[-1,+1],prctile2_10*[1,1],'color',color2,'linewidth',line_width);
line(ax,[2,2],[quarntile2_75,prctile2_90],'color',color2,'linewidth',line_width);
line(ax,[2,2],[prctile2_10,quarntile2_25],'color',color2,'linewidth',line_width);
line(ax,(2-(box_width/2))*[1,1],[quarntile2_25,quarntile2_75],'color',color2,'linewidth',line_width);
line(ax,(2+(box_width/2))*[1,1],[quarntile2_25,quarntile2_75],'color',color2,'linewidth',line_width);

xlim([0.5,2.5]);

end