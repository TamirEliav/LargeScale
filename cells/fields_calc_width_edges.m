function [widths, edges] = fields_calc_width_edges(FR_map, fields, ref_height_for_width)
widths = [];
edges = [];
for ii_field = 1:length(fields)
    
    %% arrange data
    PSTH = FR_map.PSTH;
    field = fields(ii_field);
    peak = field.peak;
    loc_ind = field.loc_IX;

    %% divide to before and after the field peak
    PSTH_right = PSTH(loc_ind+1:end);
    PSTH_left = PSTH(1:loc_ind-1);

    %% right edge
    ind = find(PSTH_right - (peak*ref_height_for_width) <= 0);
    if ~isempty(ind)
        ind_R = loc_ind + ind(1); % the right border is the first time the PSTH-half height is negative.
        loc_R = interp1([PSTH(ind_R-1) PSTH(ind_R)], [(ind_R-1) ind_R],...
            peak*ref_height_for_width,'linear');
    else
        ind = find(~isnan(PSTH_right));
        if isempty(ind)%the case in which the peak is at the edge
            ind_R = loc_ind + 1;
        else
            ind_R = loc_ind + ind(end);
        end
        loc_R = ind_R;
    end

    %% left edge
    ind = find(PSTH_left - (peak*ref_height_for_width) <= 0);
    if ~isempty(ind)
        ind_L = ind(end);
        ind_L_next = ind_L+1;
        while isnan(PSTH(ind_L_next)) && ind_L_next<=loc_ind
            ind_L_next = ind_L_next+1;
        end
        loc_L = interp1([PSTH(ind_L) PSTH(ind_L_next)], [ind_L (ind_L+1)],...
            peak*ref_height_for_width,'linear');
    else
        ind = find(~isnan(PSTH_left));
        if isempty(ind)%the case in which the peak is at the edge
            ind_L = length(PSTH_left);
        else
            ind_L = ind(1);
        end
        loc_L = ind_L;
    end

    field_edges_bin = [loc_L loc_R];
    field_size_bin = range(field_edges_bin);

    widths(ii_field) = field_size_bin .* FR_map.bin_size;
    edges(ii_field,:) = interp1(1:length(FR_map.bin_centers), FR_map.bin_centers, field_edges_bin);
end

end



