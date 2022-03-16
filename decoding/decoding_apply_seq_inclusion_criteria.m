function [seqs, TF] = decoding_apply_seq_inclusion_criteria(seqs)

TF = true(size(seqs));

TF([seqs.score]<0.5)=false;
TF([seqs.distance]<3)=false;

seqs(~TF) = [];