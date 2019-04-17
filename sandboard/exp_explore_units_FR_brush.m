function exp_explore_units_FR_brush()

%% get selected data (with brush)
hfigdata = getappdata(gcf);
hlines = [hfigdata.obj_handles_units{2,:,:}];
brush_ts = [];
for ii = 1:length(hlines)
    IX = find(hlines(ii).BrushData);
    brush_ts = [brush_ts hlines(ii).UserData(IX)];
end

%% read units from NTT file
NTT_IN = hfigdata.NTT_file;
[Timestamps, CellNumbers, Samples, Header] = ...
     Nlx2MatSpike(NTT_IN, [1 0 1 0 1], 1, 1, [] );

 %% create new unit with the selected spikes
% find the selected spikes by ts
IX = find(ismember(Timestamps, brush_ts));
new_unit_num = max(unique(CellNumbers))+1;
CellNumbers(IX) = new_unit_num;
 
%% save new NTT
[FILEPATH,NAME,EXT] = fileparts(NTT_IN);
NAME = [NAME '___selected_' string(datetime('now', 'Format', 'yyyyMMdd_HHmmss')).char];
% NAME = [NAME '___selected'];
NTT_OUT_filename = fullfile(FILEPATH, [NAME EXT]);
Mat2NlxSpike(NTT_OUT_filename, 0, 1, [], [1 0 1 0 1 1], ...
    Timestamps, CellNumbers, Samples, Header);


end
