function fields = fields_remove_overlaps(FR_map, fields, prm)

    if length(fields)<=1
        return;
    end

    %% sort by peak FR
    [~,ind_sorted] = sort( [fields.peak] ,'descend');
    fields = fields(ind_sorted);
    
    %% field edges for overlap
    [~, edges] = fields_calc_width_edges(FR_map, fields, prm.fields.overlap_href);
    [fields(:).overlap_edges] = disperse(edges);
    
    %%
    ii_field = 2;
    while ii_field <= length(fields)
        
        %%
        edges = [fields(ii_field).overlap_edges];
        edges_higher = cat(1,fields(1:ii_field-1).edges_href);
        % case 1: the field intersect with other fields
        overlapping_fields_ind1 = ...
            or( and( edges(1)>edges_higher(:,1), edges(1)<edges_higher(:,2) ),...
                and( edges(2)>edges_higher(:,1), edges(2)<edges_higher(:,2) )    );
        % case 2: the field include other fields
        overlapping_fields_ind2 = ...
            or( and( edges_higher(:,1)>edges(1), edges_higher(:,1)<edges(2) ),...
                and( edges_higher(:,2)>edges(1), edges_higher(:,2)<edges(2) )     );
        is_overlap = any(or(overlapping_fields_ind1,overlapping_fields_ind2));
        
        if is_overlap
            %if overlapping fields, remove current field (no need to increase ii_field!)
            fields(ii_field) = [];
        else
            % if not, continue to next fields
            ii_field = ii_field + 1;
        end
    end
end
