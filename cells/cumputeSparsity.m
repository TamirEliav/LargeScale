function sparsity = cumputeSparsity(PSTH,time_spent)

ind_notNaN = ~isnan(PSTH); 
r_i = PSTH(ind_notNaN);
p_i = time_spent(ind_notNaN)./sum(time_spent(ind_notNaN));
sparsity = sum( p_i .* r_i )^2 / sum( p_i .* ( r_i .^2 ) ) ; % Sparsity

end