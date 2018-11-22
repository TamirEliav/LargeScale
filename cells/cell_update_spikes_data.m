function cell_update_spikes_data(cell_ID)

%% load spikes ts

% save spikes data in mat file
filename_cell_spikes = ['L:\Analysis\Results\cells\spikes\' cell_ID 'cell_spikes'];
save(filename_cell_spikes, 'spikes')

end