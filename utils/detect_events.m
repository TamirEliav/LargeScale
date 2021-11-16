function [events, g, opts] = detect_events(x,opts)
    arguments
        x (1,:)
        opts.Fs = 1
        opts.thr = 3
        opts.thr_edges = 1
        opts.merge_thr = nan % no merging by default, units same as Fs
        opts.min_width = nan
        opts.plot = false
    end

    %% thr crossing (+merging)
    if isnan(opts.merge_thr)
        % no merging
        xthr1 = x > opts.thr;
    else
        % merging 
        xthr1 = false(1,length(x));
        xthr_IX = find(x > opts.thr);
        if ~isempty(xthr_IX)
            start_IX = [xthr_IX(1) xthr_IX(find(diff(xthr_IX)./opts.Fs > opts.merge_thr)+1 )              ];
            end_IX   = [           xthr_IX(find(diff(xthr_IX)./opts.Fs > opts.merge_thr)   )  xthr_IX(end)];

            for ii = 1:length(start_IX)
                xthr1(start_IX(ii):end_IX(ii)) = true;
            end
        end
    end
    
    %% extend thr crossing to to a lower edge_thr
    xthr2 = x > opts.thr_edges;
    xthr2(xthr1) = true; % in case we used merging
    xthr12 = extend_intervals(xthr1,xthr2);
    
    %% extract events
    cc = bwconncomp(xthr12);
    % create events struct
    events=struct();
    events.duration = 1/opts.Fs .* cellfun(@length,cc.PixelIdxList);
    events.start_IX = cellfun(@min, cc.PixelIdxList);
    events.end_IX = cellfun(@max, cc.PixelIdxList);
%     events.start_ts = ts(events.start_IX);
%     events.end_ts = ts(events.end_IX);
    g = bwlabel(xthr12);
    g(g==0)=nan;
    [~,max_IX] = splitapply(@max,x,g); % index relative to event
    events.peak_IX = cellfun(@(IX,max_IX)(IX(max_IX)), cc.PixelIdxList, num2cell(max_IX));
%     events.peak_ts = ts(events.peak_IX);
    events.peak_val = x(events.peak_IX);

    events = soa2aos(events);
    if ~isnan(opts.min_width)
        invalid = [events.duration] < opts.min_width;
        events(invalid) = [];
    end
   
    %% 
    if opts.plot || nargout==0
        figure
        hold on
        plot(x,'b')
        yline(opts.thr,'r','thr')
        yline(opts.thr_edges,'r','thr edges')
        xthr12_vals = nan(size(x));
        xthr12_vals(xthr12) = x(xthr12);
        plot(xthr12_vals,'r')
        plot([events.peak_IX], [events.peak_val], '*r');
    end
end







