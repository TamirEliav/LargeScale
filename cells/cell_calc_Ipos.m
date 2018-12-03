function cell_calc_Ipos(cell_ID)

%% load cell/exp data
cell = cell_load_data(cell_ID,'details','FE');
exp = exp_load_data(cell.details.exp_ID,'details');

%% params
prm = PARAMS_GetAll();
pos_bin_size = prm.Ipos.pos_bin_size;
time_bin_size = prm.Ipos.time_bin_size;
pos_bin_limits = prm.Ipos.pos_bin_limits;

%%
for ii_dir = 1:2
    
    %% arrange data for current direction
    spikes_ts = [cell.FE{ii_dir}.spikes_ts];
    pos = [cell.FE{ii_dir}.pos];
    pos_ts = [cell.FE{ii_dir}.ts];
    pos_invalid_IX = find(isnan(pos_ts));
    pos(pos_invalid_IX) = [];
    pos_ts(pos_invalid_IX) = [];
    
    %% TODO: check why I did that shift (I think it is to avoid lareg memory allocation)
    ti = exp_get_sessions_ti(cell.details.exp_ID,'Behave');
    time_zero = ti(1) .* 1e-6;
    spikes_ts = spikes_ts.*1e-6 - time_zero;
    pos_ts = pos_ts.*1e-6 - time_zero;
    
    %% discretization/binning
    pos_bin_edges = pos_bin_limits(1):pos_bin_size:pos_bin_limits(2);
    time_bin_edges = pos_ts(1):time_bin_size:pos_ts(end);
    pos_bins_centers = (pos_bin_edges(1:end-1) + pos_bin_edges(2:end)) / 2;
    nbins_pos = length(pos_bin_edges)-1;
    pos_pos_discretize = discretize(pos,pos_bin_edges, 1:length(pos_bin_edges)-1);
    pos_time_discretize = discretize(pos_ts, time_bin_edges, 1:length(time_bin_edges)-1);
    pos_time_discretize(isnan(pos_time_discretize)) = max(pos_time_discretize);
    X = accumarray(pos_time_discretize', pos_pos_discretize',[],@mode,nan)';
    n_pos_bins_per_time_bin = accumarray(pos_time_discretize', pos_pos_discretize',[], @(x)(length(unique(x))), nan)';
    K = histcounts(spikes_ts, time_bin_edges);
% %     K = [0 abs(diff(K))]; %% TODO: TEMP!!!! (looking at the spikes count (temporal) CHANGE!)

    invalid_time_bins_IX = find(isnan(X));
    valid_time_bins_IX = find(~isnan(X));
    X(invalid_time_bins_IX) = [];
    K(invalid_time_bins_IX) = [];
    n_pos_bins_per_time_bin(invalid_time_bins_IX) = [];

    %% calc basic probabilities
    kmax = max(K);
    nspikes_bin_edges = -0.5:1:(kmax+0.5);
    nbins_nspikes = length(nspikes_bin_edges)-1;
    % #spikes probability in a general temporal bin
    P_k = histcounts(K, nspikes_bin_edges,'Normalization','probability');
    % #spikes probability in a temporal bin, given specific position
    P_k_xi = zeros(nbins_pos, nbins_nspikes);
    for ii_x = 1:nbins_pos
        IX = find( X == ii_x );
        P_k_xi(ii_x,:) = histcounts(K(IX), nspikes_bin_edges,'Normalization','probability');
    end
    
    %% calc local spatial information
    Ipos_data = zeros(1, nbins_pos);
    for ii_x = 1:nbins_pos
        sdf = 0;
        for ii_k = 1:nbins_nspikes
            sdf = sdf + P_k_xi(ii_x, ii_k) * ( log2( P_k_xi(ii_x,ii_k)+eps) - log2(P_k(ii_k)+eps) );
        end
        Ipos_data(ii_x) = sdf;
    end

    %% arrange results
    Ipos.params.pos_bin_size = pos_bin_size;
    Ipos.params.time_bin_size = time_bin_size;
    Ipos.data(ii_dir).Ipos = Ipos_data;
    Ipos.data(ii_dir).pos_bin_edges = pos_bin_edges;
    Ipos.data(ii_dir).pos_bins_centers = pos_bins_centers;
    Ipos.data(ii_dir).time_bin_edges = time_bin_edges;
    Ipos.data(ii_dir).X = X;
    Ipos.data(ii_dir).K = K;
    Ipos.data(ii_dir).P_k = P_k;
    Ipos.data(ii_dir).P_k_xi = P_k_xi;
    Ipos.data(ii_dir).n_pos_bins_per_time_bin = n_pos_bins_per_time_bin;
end

%% save data to file
filename = fullfile('L:\Analysis\Results\cells\Ipos',[cell_ID '_cell_Ipos']);
save(filename, 'Ipos');









end







