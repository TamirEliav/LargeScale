function [p, p_inverse internals] = sync_TTL_polyfit(TTL_ts_1, TTL_ts_2, N_smoothing, max_diff)
% TTL ts in msec !!!
% TTL 1 and 2 can have different number of TTLs
% p is the linear polyfit that converts time_1 to time_2
% p_inverse is the linear polyfit that converts time_2 to time_1
% TODO: plot more details in the out figures so user can estimate easily if
% the sync is OK (maybe create a automatic WARNING if the values are not good...)


%%
if isempty(N_smoothing)
    N_smoothing = 5;
end
figure('units','normalized', 'position', [0 0 1 1])

%% create time-series vectors
TTL_vec_1 = zeros(1, TTL_ts_1(end)-TTL_ts_1(1)+1);
TTL_vec_2 = zeros(1, TTL_ts_2(end)-TTL_ts_2(1)+1);
TTL_vec_1( ceil(TTL_ts_1-TTL_ts_1(1)) + 1 ) = 1;
TTL_vec_2( ceil(TTL_ts_2-TTL_ts_2(1)) + 1 ) = 1;
% % % eliminate the edge effects (there is always a big peak at the edge...)
% % zero_padd_duration = 60*1e3;
% % zero_padd = zeros(1,zero_padd_duration);
% % TTL_vec_1 = [zero_padd TTL_vec_1 zero_padd];
% % TTL_vec_2 = [zero_padd TTL_vec_2 zero_padd];
if N_smoothing > 1
    TTL_vec_1 = filtfilt(hamming(N_smoothing),1,TTL_vec_1);
    TTL_vec_2 = filtfilt(hamming(N_smoothing),1,TTL_vec_2);
end



%% cross correlate (find shift)
[c lags] = xcorr(TTL_vec_1, TTL_vec_2);
[~, IX] = max(c);
shift_msec = -lags(IX) + TTL_ts_2(1) - TTL_ts_1(1);
TTL_ts_1_shifted = TTL_ts_1 + shift_msec;

[PKS LOCS] = findpeaks(c,'MINPEAKHEIGHT',max(c)/100);
[PKS_sorted PKS_sorted_IX] = sort(PKS, 'descend');
LOCS_sorted = LOCS(PKS_sorted_IX);
sidelobe_ratio = PKS_sorted(1) / PKS_sorted(2);
sidelobe_dB = 20*log10(sidelobe_ratio);
sidelobe_delay = abs(LOCS_sorted(1) - LOCS_sorted(2));

subplot(1,3,1)
hold on
plot(lags,c)
plot(lags(LOCS),c(LOCS),'or')
text(lags(LOCS_sorted(2))+200,c(LOCS_sorted(2)), ...
    {['1^{st} sidelobe: ' num2str(sidelobe_dB) 'dB'];...
     [num2str(sidelobe_delay) 'ms lag']},...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
 
%% find best TTL pairs
temp_diff = [];
IX_closest_TTL_2 = [];
number_of_pairs = min( length(TTL_ts_1), length(TTL_ts_2) );
for ii_TTL = 1:length(TTL_ts_1_shifted)
    [temp_diff(ii_TTL), IX_closest_TTL_2(ii_TTL)] = min( abs(TTL_ts_1_shifted(ii_TTL) - TTL_ts_2) );
end
[temp_diff_sorted,  IX_TTL_1_temp_diff_sorted] = sort(temp_diff,'ascend');
IX_matching_TTL_1 = IX_TTL_1_temp_diff_sorted(1:number_of_pairs);
IX_matching_TTL_2 = IX_closest_TTL_2(IX_matching_TTL_1);
IX_matching_TTL_1 = sort(IX_matching_TTL_1);
IX_matching_TTL_2 = sort(IX_matching_TTL_2);

%% remove pairs with temporal diff > max_diff
X = TTL_ts_1_shifted(IX_matching_TTL_1);
Y = TTL_ts_2(IX_matching_TTL_2);
bad_pairs_IX = find( abs(X-Y) > max_diff);
IX_matching_TTL_1(bad_pairs_IX) = [];
IX_matching_TTL_2(bad_pairs_IX) = [];

%% 1st fit (before removing bad points)
X = TTL_ts_1(IX_matching_TTL_1);
Y = TTL_ts_2(IX_matching_TTL_2);
p = polyfit(X,Y,1);
p_inverse = polyfit(Y,X,1);

subplot(1,3,2)
hold on
plot(X,Y,'o')
plot(X,polyval(p,X),'r')
legend({'data points';'polyval'})
title(['1^{st} fit (before removing bad points), n=' num2str(length(X))])

%% remove pairs with temporal diff > max_diff
X = TTL_ts_1(IX_matching_TTL_1);
Y = TTL_ts_2(IX_matching_TTL_2);
bad_pairs_IX = find( abs(polyval(p,X)-Y) > max_diff);
IX_matching_TTL_1(bad_pairs_IX) = [];
IX_matching_TTL_2(bad_pairs_IX) = [];

%% 2nd fit (after removing bad points)
X = TTL_ts_1(IX_matching_TTL_1);
Y = TTL_ts_2(IX_matching_TTL_2);
p = polyfit(X,Y,1);
p_inverse = polyfit(Y,X,1);

subplot(1,3,3)
hold on
plot(X,Y,'o')
plot(X,polyval(p,X),'r')
legend({'data points';'polyval'})
title(['2^{nd} fit (after removing bad points), n=' num2str(length(X))])

%%
internals.shift_msec = shift_msec;
internals.number_of_pairs = number_of_pairs;
internals.c = c;
internals.lags = lags;


