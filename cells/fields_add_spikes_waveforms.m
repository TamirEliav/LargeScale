function cell = fields_add_spikes_waveforms(cell)

%%
for ii_dir = 1:2
    for ii_field = 1:length(cell.fields{ii_dir})
        %%
        field_spikes_ts = cell.fields{ii_dir}(ii_field).spikes_ts;
        [~,IX] = ismember(field_spikes_ts , cell.spikes.ts);
        cell.fields{ii_dir}(ii_field).spikes_wvfrm = cell.spikes.waveforms(:,:,IX);
    end
end




end