function [segment_appear_flight, segment_disappear_flight, fields_per_win] = map_detect_segments_appear_disappear(fields_per_win,FE,options)
% calc probability of appearing and disappearing field's segments:
% plan - 
% identify changes in field size which are larger than a fixed threshold.
% then check if these changes are larger than a value relative to the mean
% field size before the change (50% of field size before change?). then
% chack if the changes is consistent, i.e. stable over the next time window
% of at least n windows (3 - 5 windows) which are not nans.

min_change_size = 1; %1m
min_n_windows = 5;
min_prc_change = 50; % 50% of previousley mean field size (not including the change)
min_large_change = 5; %if a change is larger than this number, it is 
%considered a change even if it doesn't reach the precentage criterion 
large_gaps_th = 5;

detect_field_size_change = false;
detect_field_fr_change = false;
detect_field_loc_change = false;
detect_large_gaps = false;
if nargin == 2 %no option input
    %default detection:
    detect_field_size_change = true;
    detect_field_fr_change = false;
    detect_field_loc_change = true;
    detect_large_gaps = true;
else
    if any(strcmp(options, 'size'))
        detect_field_size_change = true;
    end
    if any(strcmp(options, 'fr'))
        detect_field_fr_change = true;
    end
    if any(strcmp(options, 'loc'))
        detect_field_loc_change = true;
    end
    if any(strcmp(options, 'gap'))
        detect_large_gaps = true;
    end
end
    
n_flight_in_session = length(FE);

segment_appear_flight = zeros([1,n_flight_in_session]); %appear events
segment_disappear_flight = zeros([1,n_flight_in_session]); %disappear events
    
for ii_field = 1:length(fields_per_win)
    IX_size_change_appear = [];
    IX_size_change_disappear = [];
    IX_FR_change_appear = [];
    IX_FR_change_disappear = [];
    IX_loc_change_appear = [];
    IX_loc_change_disappear = [];
    IX_seg_appear = [];
    IX_seg_disappear = [];
    
    %% detect abrupt field size changes:
    if detect_field_size_change
        field_size = fields_per_win(ii_field).field_size;
        field_size_spikes = fields_per_win(ii_field).field_size_spikes;
        % ignore value before start of field and after end of field:
        if fields_per_win(ii_field).start_flight>1
            field_size_spikes(1:fields_per_win(ii_field).start_flight-1) = nan;
        end
        if fields_per_win(ii_field).end_flight<length(field_size_spikes)
            field_size_spikes(fields_per_win(ii_field).end_flight:end) = nan;
        end
        if sum(~isnan(field_size_spikes))>min_n_windows
            [IX_change, part_mean, part_change, part_change_type] = ...
                get_parts(field_size_spikes, min_n_windows, min_change_size);
            % remove relative small changes:
            [IX_size_change,IX_size_change_type] = remove_relative_small_changes(...
                IX_change, part_mean, part_change, part_change_type, min_prc_change, min_change_size, min_large_change);
            
            IX_size_change_appear = IX_size_change(IX_size_change_type>0);
            IX_size_change_disappear = IX_size_change(IX_size_change_type<0);
            num_spikes_per_flight = fields_per_win(ii_field).num_spikes_per_flight;
            [IX_size_change_disappear] = get_last_window_before_loc_change(IX_size_change_disappear, field_size, num_spikes_per_flight);
        end
    end
    
    %% detect abrupt FR changes:
    if detect_field_fr_change
        field_size = fields_per_win(ii_field).field_size;
        %spike count is taken here as all the spikes occuring inside
        %the field's edges (constant value 25-75%), and not only for
        %detected windows:
        spikes_count = fields_per_win(ii_field).num_spikes_constant_edges;
        % ignore value before start of field and after end of field:
        if fields_per_win(ii_field).start_flight>1
            spikes_count(1:fields_per_win(ii_field).start_flight-1) = nan;
        end
        if fields_per_win(ii_field).end_flight<length(field_size)
            spikes_count(fields_per_win(ii_field).end_flight:end) = nan;
        end
        % ignore zeros spikes flights:
        spikes_count(spikes_count==0) = nan;
        
        if sum(~isnan(spikes_count))>min_n_windows
            mean_field_size = nanmean(field_size);
            spikes_count_per_m = spikes_count./mean_field_size;
            [IX_change, part_mean, part_change, part_change_type] = ...
                get_parts(spikes_count_per_m, min_n_windows, min_change_size);
            % remove relative small changes:
            [IX_FR_change,IX_FR_change_type] = remove_relative_small_changes(...
                IX_change, part_mean, part_change, part_change_type, min_prc_change, min_change_size, min_large_change);
            IX_FR_change_appear = IX_FR_change(IX_FR_change_type>0);
            IX_FR_change_disappear = IX_FR_change(IX_FR_change_type<0);
        end
    end
    
    %% detect abrupt position changes:
    if detect_field_loc_change
        field_size = fields_per_win(ii_field).field_size;
        field_loc = fields_per_win(ii_field).CoM;
        [IX_change, part_mean, part_change, part_change_type] = ...
            get_parts(field_loc, min_n_windows, min_change_size);
        % remove relative small changes:
        IX_part = [1; IX_change; length(field_size)];
        part_mean_size = [];
        for ii_part = 1:length(part_mean)
            part_mean_size(ii_part) = nanmean(field_size(IX_part(ii_part):IX_part(ii_part+1)));
        end
        
        [IX_loc_change_appear,~] = remove_relative_small_changes(...
            IX_change, part_mean_size, part_change, part_change_type, min_prc_change, min_change_size, min_large_change);
        
        % remove loc change events based on overlap of fields before and
        % after the change:
        IX_loc_change = IX_loc_change_appear;
        IX_not_size_change = [];
        IX_not_loc_change = [];
        for ii_change = 1:length(IX_loc_change)
            edges_before = fields_per_win(ii_field).edges(1:IX_loc_change(ii_change)-1,:);
            edges_before(isnan(edges_before(:,1)),:) = [];
            edges_before = edges_before(end,:);
            
            edges_after = fields_per_win(ii_field).edges(IX_loc_change(ii_change):end,:);
            edges_after(isnan(edges_after(:,1)),:) = [];
            edges_after = edges_after(1,:);
            
            overlap = min(edges_before(2),edges_after(2)) - max(edges_before(1),edges_after(1));
            overlap(overlap<0) = 0;
            overlap_prc = 100 * overlap/min(diff(edges_before),diff(edges_after));
            if overlap_prc>=50
                IX_not_loc_change = [IX_not_loc_change ii_change];
            else
                IX_not_size_change = [IX_not_size_change IX_loc_change(ii_change)];
            end
        end
        
        IX_loc_change_appear(IX_not_loc_change) = [];
        num_spikes_per_flight = fields_per_win(ii_field).num_spikes_per_flight;
        [IX_loc_change_appear] = get_closest_flight_with_spikes(IX_loc_change_appear, field_loc, num_spikes_per_flight);
        [IX_loc_change_disappear] = get_last_window_before_loc_change(IX_loc_change_appear, field_loc, num_spikes_per_flight);
        
        % remove duplicates of events:
        IX_size_change_appear = remove_events_duplications(IX_size_change_appear,IX_not_size_change);
        IX_size_change_disappear = remove_events_duplications(IX_size_change_disappear,IX_not_size_change);
        
        IX_not_FR_change = unique([IX_loc_change_appear; IX_loc_change_disappear])';
        IX_not_FR_change = unique([IX_not_FR_change IX_not_FR_change-1 IX_not_FR_change+1]);
        IX_FR_change_appear = remove_events_duplications(IX_FR_change_appear,IX_not_FR_change);
        IX_FR_change_disappear = remove_events_duplications(IX_FR_change_disappear,IX_not_FR_change);
    end
    
    %% detect large gaps:
    if detect_large_gaps
        num_spikes_per_flight = fields_per_win(ii_field).num_spikes_per_flight;
        IX_flights_w_spikes = find(num_spikes_per_flight>0);
        diff_IX_flights_w_spikes = diff(IX_flights_w_spikes);
        IX_gaps_start = IX_flights_w_spikes(diff_IX_flights_w_spikes>=large_gaps_th)+1;
        IX_gaps_end = IX_flights_w_spikes(find(diff_IX_flights_w_spikes>=large_gaps_th)+1);
        
        gap_to_remove = [];
        for ii_gap = 1:length(IX_gaps_start)
            zero_spikes_flights = sum(...
                num_spikes_per_flight(IX_gaps_start(ii_gap):(IX_gaps_end(ii_gap)-1))==0);
            if zero_spikes_flights<large_gaps_th
                gap_to_remove = [gap_to_remove ii_gap];
            end
        end
        IX_gaps_start(gap_to_remove) = [];
        IX_gaps_end(gap_to_remove) = [];
        
        % refine gap definition with spikes within range and not only the
        % detected spikes that were added to certain window detection:
        gap_to_remove = [];
        edges = fields_per_win(ii_field).edges;
        IX_valid_edges = find(~isnan(edges(:,1)));
        for ii_gap = 1:length(IX_gaps_start)
            [~, closest_IX] = min(abs(IX_valid_edges - IX_gaps_start(ii_gap)));
            closest_IX = IX_valid_edges(closest_IX);
            edges_start_gap = edges(closest_IX,:);
            [~, closest_IX] = min(abs(IX_valid_edges - IX_gaps_end(ii_gap)));
            closest_IX = IX_valid_edges(closest_IX);
            edges_end_gap = edges(closest_IX,:);
            edges_gap = [min(edges_start_gap(1),edges_end_gap(1)) ...
                max(edges_start_gap(2), edges_end_gap(2))];
            
            FE_no_spikes_IX = cellfun(@(x)(isempty(get_data_in_ti(x,edges_gap))),{FE.spikes_pos});
            FE_no_spikes_IX = FE_no_spikes_IX((IX_gaps_start(ii_gap)-1):IX_gaps_end(ii_gap));
            
            IX_spikes_flights_in_gap = [1 find(~FE_no_spikes_IX) length(FE_no_spikes_IX)];
            gaps_in_gap = diff(IX_spikes_flights_in_gap)-1;
            
            if max(gaps_in_gap)<large_gaps_th % spikes_flights_in_gap
                gap_to_remove = [gap_to_remove ii_gap];
            end
        end
        IX_gaps_start(gap_to_remove) = [];
        IX_gaps_end(gap_to_remove) = [];
        IX_seg_appear = IX_gaps_end;
        IX_seg_disappear = IX_gaps_start;
        
        % remove event duplications:
        IX_not_size_change = unique([IX_gaps_start; IX_gaps_end])';
        IX_not_size_change = unique([IX_not_size_change IX_not_size_change-1 IX_not_size_change+1]);
        IX_size_change_appear = remove_events_duplications(IX_size_change_appear,IX_not_size_change);
        IX_size_change_disappear = remove_events_duplications(IX_size_change_disappear,IX_not_size_change);
        IX_loc_change_appear = remove_events_duplications(IX_loc_change_appear,IX_not_size_change);
        IX_loc_change_disappear = remove_events_duplications(IX_loc_change_disappear,IX_not_size_change);
        IX_FR_change_appear = remove_events_duplications(IX_FR_change_appear,IX_not_size_change);
        IX_FR_change_disappear = remove_events_duplications(IX_FR_change_disappear,IX_not_size_change);
    end
    
    %% merge changes types
    IX_change_appear = unique([IX_size_change_appear; IX_loc_change_appear; IX_seg_appear; IX_FR_change_appear]);
    IX_change_disappear = unique([IX_size_change_disappear; IX_loc_change_disappear; IX_seg_disappear; IX_FR_change_disappear]);
    IX_all_changes = unique([IX_change_appear; IX_change_disappear]);
    segment_appear_flight(IX_change_appear) = segment_appear_flight(IX_change_appear) +1;
    segment_disappear_flight(IX_change_disappear) = segment_disappear_flight(IX_change_disappear) +1;
    
    % save detected segments changes
    fields_per_win(ii_field).seg_appear_flight = IX_change_appear;
    fields_per_win(ii_field).seg_disappear_flight = IX_change_disappear;
    fields_per_win(ii_field).seg_all_changes = IX_all_changes;
    
    
    %add flight start and end:
    fields_start_flight = fields_per_win(ii_field).start_flight;
    segment_appear_flight(fields_start_flight) = segment_appear_flight(fields_start_flight) +1;
    fields_end_flight = fields_per_win(ii_field).end_flight;
    segment_disappear_flight(fields_end_flight) = segment_disappear_flight(fields_end_flight) +1;
    
end

%% discard flights in which a segment could not appear or disappear:
prm = PARAMS_GetAll();
min_n_flights_per_field = ceil(max(prm.fields.min_flights_with_spikes, ...
    prm.fields.min_flights_with_spikes_prc * n_flight_in_session));
% counting all the flights in which a field potentially could
% appear: all the flights in the session minus the first 2 flights and
% the last (5-1) flights. in these flights if a field was created it
% could not be counted as a new field.
% similarly for the case of a disappearing field, we will remove the
% last 1 flight and the first 5 flights (it can disappear on the 6th
% flight and on the one before last flights).

remove_flightsIX_appear = [1 2 (((-min_n_flights_per_field+2):0)+n_flight_in_session)];
remove_flightsIX_disappear = [(1:min_n_flights_per_field) n_flight_in_session];
segment_appear_flight(remove_flightsIX_appear) = [];
segment_disappear_flight(remove_flightsIX_disappear) = [];

end



%% 
function [IX_change, part_mean, part_change, part_change_type] = get_parts(field_size, min_n_windows, min_change_size)
% field_size is a seuqence of field sizes per window
% get parts devides the field size sequence into parts that has different
% mean field size.

[TF,mean_s] = ischange(field_size,'MaxNumChanges',5);
IX_change = find(TF);
part_start = [1; IX_change];
part_end = [IX_change-1; length(TF)];
part_len = part_end - part_start +1;
part_mean = zeros([length(part_len),1]);
part_len_noNans = zeros([length(part_len),1]);
for ii_part = 1:length(part_len)
    part_fields = field_size(part_start(ii_part):part_end(ii_part));
    part_mean(ii_part) = nanmean(part_fields);
    part_fields(isnan(part_fields)) = [];
    part_len_noNans(ii_part) = length(part_fields);
end
part_change = [abs(diff(part_mean)) ; nan];
part_change_type = [diff(part_mean) ; nan];

% remove 1 flight outliers:
to_remove = [];
to_nan = [];
for ii_part = 1:length(part_len)
    if (part_len_noNans(ii_part)==1)
        if ii_part==1 || ii_part==length(part_len)
            continue;
        end
        cond = min(part_change(ii_part),part_change(ii_part-1)) >= ...
            0.5*max([part_mean(ii_part), part_mean(ii_part+1), part_mean(ii_part-1)]);
        if cond
            to_remove = [to_remove ii_part];
            to_nan = [to_nan part_start(ii_part):part_end(ii_part)];
            continue;
        end
    end
end
part_start(to_remove) = [];
part_end(to_remove) = [];
part_mean(to_remove) = [];
part_change = [abs(diff(part_mean)) ; nan];
part_change_type = [diff(part_mean) ; nan];
part_len_noNans(to_remove) = [];
field_size(to_nan) = nan;

% remove short and small change parts:
ii_part = 1;
while ii_part<=length(part_start)
    if length(part_start)==1
        break
    end
    if (part_len_noNans(ii_part)<min_n_windows) || (part_change(ii_part)<min_change_size)
        if ii_part==1 || part_change(ii_part)<part_change(ii_part-1)
            part_start(ii_part+1) = part_start(ii_part);
            part_fields = field_size(part_start(ii_part+1):part_end(ii_part+1));
            part_mean(ii_part+1) = nanmean(part_fields);
            part_fields(isnan(part_fields)) = [];
            part_len_noNans(ii_part+1) = length(part_fields);
        else
            part_end(ii_part-1) = part_end(ii_part);
            part_fields = field_size(part_start(ii_part-1):part_end(ii_part-1));
            part_mean(ii_part-1) = nanmean(part_fields);
            part_fields(isnan(part_fields)) = [];
            part_len_noNans(ii_part-1) = length(part_fields);
        end
        part_start(ii_part) = [];
        part_end(ii_part) = [];
        part_mean(ii_part) = [];
        part_change = [abs(diff(part_mean)) ; nan];
        part_change_type = [diff(part_mean) ; nan];
        part_len_noNans(ii_part) = [];
        
        continue;
    end
    ii_part = ii_part+1;
end

IX_change = part_start;
IX_change(1) = [];
        
end

%% 
function [IX_change, IX_change_type] = remove_relative_small_changes(IX_change, part_mean, part_change, part_change_type, min_prc_change, min_change_size, min_large_change)

IX_remove = [];
for ii_part = 1:length(IX_change)
    prc_th = min_prc_change * max(part_mean(ii_part),part_mean(ii_part+1)) / 100;
    prc_th = max(prc_th, min_change_size); %must be sill larger than 1m
    % criterion to pass: change need to be larger than the 50% precentile
    % OR larger than 5m in order to be counted as a valid change. if none of
    % these 2 conditions is reached - remove the change:
    if ~or(part_change(ii_part)>prc_th, part_change(ii_part)>min_large_change)
        IX_remove = [IX_remove ii_part];
    end
end
IX_change(IX_remove) = [];
part_change_type(IX_remove) = [];
part_change_type(isnan(part_change_type)) = [];
IX_change_type = sign(part_change_type);

if isempty(IX_change) %just fix mess-ups with 0xN empty arrays
    IX_change = [];
    part_change_type = [];
end
end

%%
function [IX_loc_change_appear] = get_closest_flight_with_spikes(IX_loc_change_appear, field_loc, num_spikes_per_flight)

IX_flights_w_spikes = flip(find(num_spikes_per_flight>0));
%flip so next when finding the minimum there will be higher priority for
%the flights with spikes occuring after the loc change
for ii_change = 1:length(IX_loc_change_appear)
    
    [~, closest_IX] = min(abs(IX_flights_w_spikes - IX_loc_change_appear(ii_change)));
    closest_IX = IX_flights_w_spikes(closest_IX);
    IX_loc_change_appear(ii_change) = closest_IX;
end


end
%%
function [IX_loc_change_disappear] = get_last_window_before_loc_change(IX_loc_change_appear, field_loc, num_spikes_per_flight)

IX_loc_change_disappear = [];
for ii_change = 1:length(IX_loc_change_appear)
    loc_tmp = field_loc(1:(IX_loc_change_appear(ii_change)-1));
    ind_tmp = find(~isnan(loc_tmp));
    
    win_tmp = num_spikes_per_flight(1:ind_tmp(end));
    ind_tmp = find(win_tmp>0);
    IX_loc_change_disappear(ii_change) = ind_tmp(end)+1;
    %maybe think of adding here acondition for nan flights that didn't pass
    %through the field.
end
IX_loc_change_disappear = reshape(IX_loc_change_disappear,[max(size(IX_loc_change_disappear)),1]);

end

%%
function new_IX = remove_events_duplications(IX,IX_to_remove)

if ~isempty(IX_to_remove) && ~isempty(IX)
    IX(sum(IX==IX_to_remove,2)>0) = [];
end

new_IX = IX;
end