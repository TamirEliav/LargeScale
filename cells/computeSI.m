function [SI_bits_spike, SI_bits_sec] = computeSI(PSTH,time_spent)

ind_notNaN = ~isnan(PSTH); 
prob = time_spent(ind_notNaN)./sum(time_spent(ind_notNaN));
r_i = PSTH(ind_notNaN);
meanFR = sum(r_i .* prob);
SI = sum((prob.*r_i).*log2((r_i+eps)/meanFR));

SI_bits_spike = SI/meanFR;
SI_bits_sec = SI;

end