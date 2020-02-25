function fields = fields_remove_unstable(fields, prm)

if isempty(fields)
    return
end

unstable_fields_IX = [];
for ii_field = 1:length(fields)
    num_flights_thr = max([prm.fields.min_flights_with_spikes, ...
                           prm.fields.min_flights_with_spikes_prc * fields(ii_field).FE_field_pass_num]);
    
    switch prm.parmaset 
        case 9
            % in this case, we count flights with at least the mean number
            % of spikes.
            num_spikes_expected_per_flight = nanmean(fields(ii_field).num_spikes_per_flight);
            num_flights_with_expected_spikes = sum(fields(ii_field).num_spikes_per_flight > num_spikes_expected_per_flight);
            if num_flights_with_expected_spikes < num_flights_thr
                unstable_fields_IX = [unstable_fields_IX ii_field];
            end
        otherwise
            % normally, we count flights with at least one spike
            if fields(ii_field).num_flights_with_spikes < num_flights_thr
                unstable_fields_IX = [unstable_fields_IX ii_field];
            end
   end

end

fields(unstable_fields_IX) = [];

end
