function decode = decoding_read_decoded_file(filename)

%% read decoded file
h = h5info(filename);
decode = struct();
decode.pos = h5read(filename,'/position')';
decode.time = h5read(filename,'/time')';
decode.state = h5read(filename,'/state')';
decode.state = deblank(string(cell2mat(decode.state)))';
decode.posterior = h5read(filename,'/acausal_posterior');
if contains('likelihood',{h.Datasets.Name})
    decode.likelihood = h5read(filename,'/likelihood');
    decode.likelihood = decode.likelihood ./ sum(decode.likelihood,[1 2]);
end
decode.params = attr2struct(h.Attributes);

%% process decoded results
decode.posterior_state = squeeze(sum(decode.posterior,1));
decode.posterior_pos = squeeze(sum(decode.posterior,2));
[~,decode.MAP_pos_IX] = max(decode.posterior_pos,[],1);
[~,decode.MAP_state_IX] = max(decode.posterior_state,[],1);
decode.MAP_pos = decode.pos(decode.MAP_pos_IX);
direction_by_state = zeros(1,length(decode.state));
direction_by_state(contains(decode.state,'Outbound')) = 1;
direction_by_state(contains(decode.state,'Inbound')) = -1;
decode.MAP_direction = direction_by_state(decode.MAP_state_IX);
decode.Fs = 1e6/median(diff(decode.time));

end

%%
function attrs_st = attr2struct(attrs)

    %%
    values = {attrs.Value};
    names = {attrs.Name};
    names = regexprep(names,'^_','','emptymatch');
    attrs_st = cell2struct(values,names,2);

end
