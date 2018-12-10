function fields = fields_limit_width_by_neighbors(FR_map, fields, prm)

if prm.fields.overlap_href == prm.fields.width_href
    return
end
if length(fields) <= 1
    return
end

%% sort by location
[~,IX] = sort([fields.loc],'ascend');
fields = fields(IX);

%%
fields_edges = cat(1,fields.edges_href);
checkEdges = (fields_edges(2:end,1) - fields_edges(1:end-1,1))<0;
ind = find(checkEdges);
%field i left border minus field i-1 right border. if it's  smalller
%than 0 then, the field i overlap with field i-1.
% TODO: there is an option for bug/undefined condition here! because:
% (1) we only compare neighboring fields
% (2) we don't consider the option of field inside another field (like we did in the overlap function...)
for ii_ind = 1:length(ind)
    ii_field = ind(ii_ind);
    win_edges = [fields(ii_field).loc fields(ii_field+1).loc];
    PSTH_win_IX = get_data_in_ti(FR_map.bin_centers, win_edges);
    PSTH_win = nan(size(FR_map.PSTH));
    PSTH_win(PSTH_win_IX) = FR_map.PSTH(PSTH_win_IX);
    [~,valley_IX] = min(PSTH_win);
    valley_loc = FR_map.bin_centers(valley_IX);
    
    if valley_loc < fields(ii_field).edges_href(2)
        fields(ii_field).edges_href(2) = valley_loc ;
        fields(ii_field).width_href = range(fields(ii_field).edges_href);
    end
    if valley_loc > fields(ii_field+1).edges_href(1)
        fields(ii_field+1).edges_href(1) = valley_loc;
        fields(ii_field+1).width_href = range(fields(ii_field+1).edges_href);
    end
end




end