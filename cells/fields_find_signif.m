function fields = fields_find_signif(FE, fields, prm)

%%
[fields(:).signif] = disperse( repelem(false,length(fields)) ); % pre-define field 'signif' to allow parfor
if length(fields)==0
    return;
end
    

%%
parfor ii_field = 1:length(fields)
    
    %% create FE struct specifically for the flights passed through the field
    FE_field = repelem(struct(),length(FE)); % create struct
    xi = fields(ii_field).edges_href;
    xi = xi + [-1 1] .* prm.fields.local_shuffle.margin .* range(xi); % take margins around the field
    [pos_IX_per_flight,~,pos_per_flight,~] = cellfun(@get_data_in_ti, {FE.pos}, repelem({xi},length(FE)) , 'UniformOutput', false);
    ts_per_flight = cellfun(@(x,IX)(x(IX)), {FE.ts}, pos_IX_per_flight, 'UniformOutput', false);
    [spikes_pos_IX_per_flight,~,spikes_pos_per_flight,~] = cellfun(@get_data_in_ti, {FE.spikes_pos}, repelem({xi},length(FE)) , 'UniformOutput', false);
    spikes_ts_per_flight = cellfun(@(x,IX)(x(IX)), {FE.spikes_ts}, spikes_pos_IX_per_flight, 'UniformOutput', false);
    [FE_field.pos] = disperse(pos_per_flight);
    [FE_field.ts] = disperse(ts_per_flight);
    [FE_field.spikes_pos] = disperse(spikes_pos_per_flight);
    [FE_field.spikes_ts] = disperse(spikes_ts_per_flight);
    FE_field(cellfun(@isempty,{FE_field.pos})) = []; % remove flights that didn't pass the field!
    [FE_field.start_ts] = disperse( cellfun(@(x)(x(1)), {FE_field.ts}, 'UniformOutput',false) );
    [FE_field.end_ts] = disperse( cellfun(@(x)(x(end)), {FE_field.ts}, 'UniformOutput',false) );

    %% calc FR specifically for FE_field
    FE_PSTH = FE_compute_PSTH(FE_field);
    
    %% run shuffling 
    FE_PSTH_shuffle = FE_compute_PSTH_shuffle(FE_field,...
                                              prm.fields.local_shuffle.n_shuffles,...
                                              prm.fields.local_shuffle.max_shift);
    SI_shuffle = [FE_PSTH_shuffle.FE_PSTH.SI_bits_spike];
    fields(ii_field).signif = FE_PSTH.SI_bits_spike > prctile(SI_shuffle,prm.fields.local_shuffle.signif_SI_prc);
   
end
% remove non-significant fields!
fields(~[fields.signif])  = [];

end




