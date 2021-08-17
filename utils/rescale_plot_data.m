function rescale_plot_data(varargin)
% usage examples:
% ---------------
% rescale only y axes (of gca), convert to minutes, starting from t0
% ts in usec with some t0 offset
% plot(pos,pos_ts,'k'); 
% plot(spikes_pos,spikes_ts,'r')
% rescale_plot_data('y',[1e-6/60,t0]);
% 
% specify different axes (not gca)
% rescale_plot_data(AX, 'x', [gain_x offset_x], 'y', [gain_y offset_y]); %

%% parse variables
if nargin == 0
    return
end

ii_argin = 1;

%% check if user specified a different axis
if isgraphics(varargin(ii_argin),'axes')
    h=varargin(1);
    ii_argin = ii_argin+1;
else
    h= gca;
end

while ii_argin<nargin
    % check which axes to change (X/Y/Z)
    axes_str = upper(varargin{ii_argin});
    gain   = varargin{ii_argin+1}(1);
    offset = varargin{ii_argin+1}(2);
    
    switch axes_str
        case 'X'
            h.XLim = my_rescale(h.XLim, offset, gain);
            for ii_child = 1:length(h.Children)
                hchild = h.Children(ii_child);
                switch class(hchild)
                    case {'matlab.graphics.chart.primitive.Line',...
                          'matlab.graphics.chart.primitive.Area',...
                          'matlab.graphics.chart.primitive.Bar'}
                        hchild.XData = my_rescale(hchild.XData, offset, gain);
                    case 'matlab.graphics.primitive.Text'
                        hchild.Position(1) = my_rescale(hchild.Position(1), offset, gain);
                    case 'matlab.graphics.primitive.Patch'
                        hchild.Vertices(:,1) = my_rescale(hchild.Vertices(:,1), offset, gain);
                end
            end
        case 'Y'
            h.YLim = my_rescale(h.YLim, offset, gain);
            for ii_child = 1:length(h.Children)
                hchild = h.Children(ii_child);
                switch class(hchild)
                    case {'matlab.graphics.chart.primitive.Line',...
                          'matlab.graphics.chart.primitive.Area',...
                          'matlab.graphics.chart.primitive.Bar'}
                        hchild.YData = my_rescale(hchild.YData, offset, gain);
                    case 'matlab.graphics.primitive.Text'
                        hchild.Position(2) = my_rescale(hchild.Position(2), offset, gain);
                    case 'matlab.graphics.primitive.Patch'
                        hchild.Vertices(:,2) = my_rescale(hchild.Vertices(:,2), offset, gain);
                end
            end
        case 'Z'
            h.ZLim = my_rescale(h.ZLim, offset, gain);
            for ii_child = 1:length(h.Children)
                hchild = h.Children(ii_child);
                switch class(hchild)
                    case {'matlab.graphics.chart.primitive.Line',...
                          'matlab.graphics.chart.primitive.Area',...
                          'matlab.graphics.chart.primitive.Bar'}
                        hchild.ZData = my_rescale(hchild.ZData, offset, gain);
                    case 'matlab.graphics.primitive.Text'
                        hchild.Position(3) = my_rescale(hchild.Position(3), offset, gain);
                end
            end
    end
    ii_argin = ii_argin+2;
end

end

function data = my_rescale(data,offset, gain)
    data = (data-offset).*gain;
end