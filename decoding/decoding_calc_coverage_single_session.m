function [coverage, n_seqs] = decoding_calc_coverage_single_session(exp_ID, epoch_type, params_opt,nbins)

coverage = zeros(2,nbins);
n_seqs = zeros(2,1);

[events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, 'posterior');
if isempty(events)
    return;
end
seqs = [events.seq_model];
[seqs, TF] = decoding_apply_seq_inclusion_criteria(seqs);
events(~TF)=[];
if isempty(seqs)
    return;
end
seqs_edges = [seqs.start_pos_norm; seqs.end_pos_norm];
seqs_edges = [min(seqs_edges ); max(seqs_edges )]';

xbins = linspace(0,1,nbins);
directions = [1 -1];
for ii_dir = 1:2
    direction = directions(ii_dir);
    coverage(ii_dir,:) = sum( xbins>seqs_edges([seqs.state_direction]==direction,1) & ...
                              xbins<seqs_edges([seqs.state_direction]==direction,2),1);
    n_seqs(ii_dir) = sum([seqs.state_direction]==direction);
end
