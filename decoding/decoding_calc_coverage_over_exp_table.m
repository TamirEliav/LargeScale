function coverage = decoding_calc_coverage_over_exp_table(T,params_opt,nbins)

coverage_all = zeros(height(T),2,2,nbins); % [exp] X [sleep/rest] X [2 directions] X [space]
ccc_all = zeros(height(T),2); % [exp] X [2 directions]
n_seqs_all = zeros(height(T),2,2); % [exp] X [sleep/rest] X [2 directions]
for ii_exp = 1:height(T)
    % load exp data
    exp_ID = T.exp_ID{ii_exp};
    exp = exp_load_data(exp_ID,'details');

    epoch_type = 'sleep';
    [coverage_all(ii_exp,1,:,:), n_seqs_all(ii_exp,1,:)] = decoding_calc_coverage_single_session(exp_ID, epoch_type, params_opt,nbins);
    epoch_type = 'rest';
    [coverage_all(ii_exp,2,:,:),  n_seqs_all(ii_exp,2,:)]= decoding_calc_coverage_single_session(exp_ID, epoch_type, params_opt,nbins);
    ccc1 = corr(squeeze(coverage_all(ii_exp,:,1,:))');
    ccc2 = corr(squeeze(coverage_all(ii_exp,:,2,:))');
    ccc_all(ii_exp,1) = ccc1(2);
    ccc_all(ii_exp,2) = ccc2(2);
end

coverage = struct();
coverage.T = T;
coverage.coverage_all = coverage_all;
coverage.ccc_all = ccc_all;
coverage.n_seqs_all = n_seqs_all;
coverage.params_opt = params_opt;