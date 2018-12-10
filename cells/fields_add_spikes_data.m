function fields = fields_add_spikes_data(FE,FR_map,fields,prm)

for ii_field = 1:length(fields)
    field_edges = fields(ii_field).edges_href;
    IX = get_data_in_ti([FE.spikes_pos], field_edges);
    spikes_ts  = [FE.spikes_ts];
    spikes_pos = [FE.spikes_pos];
    spikes_vel = [FE.spikes_vel];
%     spikes_wvfrm = cat(3,FE.spikes_wvfrm);
    fields(ii_field).spikes_ts = spikes_ts(IX);
    fields(ii_field).spikes_pos = spikes_pos(IX);
    fields(ii_field).spikes_vel = spikes_vel(IX);
%     fields(ii_field).spikes_wvfrm = spikes_wvfrm(:,:,IX);
    % calc field size based on spikes
    fields(ii_field).edges_prc = prctile(fields(ii_field).spikes_pos, prm.fields.width_prc);
    fields(ii_field).width_prc = range(fields(ii_field).edges_prc);
    fields(ii_field).width_std = std(fields(ii_field).spikes_pos);
    % get number of spikes in field per flight
    % first, check if the bat passed the field in that flight! (if not set nan)
    FE_no_field_pass_IX = cellfun(@(x)(isempty(get_data_in_ti(x,field_edges))),{FE.pos});
    FE_field_pass_IX = ~FE_no_field_pass_IX;
    fields(ii_field).FE_field_pass_IX = FE_field_pass_IX;
    fields(ii_field).FE_field_pass_num = sum(FE_field_pass_IX);
    spikes_IX = cellfun(@(x)(get_data_in_ti(x,field_edges)),{FE.spikes_pos},'UniformOutput',false);
    fields(ii_field).num_spikes_per_flight = cellfun(@length, spikes_IX);
    fields(ii_field).num_spikes_per_flight(FE_no_field_pass_IX) = nan;
    fields(ii_field).num_spikes_per_flight_cv = nanstd(fields(ii_field).num_spikes_per_flight) / ...
                                                nanmean(fields(ii_field).num_spikes_per_flight);
    fields(ii_field).num_flights_with_spikes     = sum(fields(ii_field).num_spikes_per_flight>0);
    fields(ii_field).num_flights_with_spikes_prc = sum(fields(ii_field).num_spikes_per_flight>0) / ...
                                                   fields(ii_field).FE_field_pass_num;
end


end