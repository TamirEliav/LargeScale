function MyBrushFunc(varargin)
    %% get selected data (with brush) - from position plots
    hfigdata = getappdata(gcf);
    hlines = [hfigdata.obj_handles_units{2,:,:}];
    brush_ts = [];
    for ii = 1:length(hlines)
        IX = find(hlines(ii).BrushData);
        brush_ts = [brush_ts hlines(ii).UserData(IX)];
    end
    
    %% set selected data (with brush) - for spikes voltage plots
    hfigdata = getappdata(gcf);
    hlines = [hfigdata.obj_handles_units{3,:,:}];
    for ii = 1:length(hlines)
        spikes_ts = hlines(ii).UserData;
        data_to_brush = ismember(spikes_ts, brush_ts);
        hlines(ii).BrushData = double(data_to_brush);
    end
end
