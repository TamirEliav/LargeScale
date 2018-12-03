function fields = fields_remove_unstable(fields, prm)

if isempty(fields)==0
    return
end

prm.fields.min_flights_with_spikes = 5;
prm.fields.min_flights_with_spikes_prc = 0.2;
total_num_flights = length([fields(1).num_spikes_per_flight]);
num_flights_thr = max(prm.fields.min_flights_with_spikes, ...
                       prm.fields.min_flights_with_spikes_prc * total_num_flights);

unstable_fields = [fields.num_flights_with_spikes] < num_flights_thr;
fields(unstable_fields) = [];

end
    