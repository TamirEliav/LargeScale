function ripples_xcorr(exp_ID)

%% get exp info
exp = exp_load_data(exp_ID,'details','path','ripples');
prm = PARAMS_GetAll();

%% params
max_tdiff = 50e3;

%%
ripples_ts = {};
nTT = length(exp.ripples.by_TT);
for TT = 1:nTT
    ripples = exp.ripples.by_TT{TT};
    if isempty(ripples)
        continue;
    end
    ripples_ts{TT} = [ripples.peak_ts];
end
% add ripples ts detected by all TTs
ripples = exp.ripples.all;
ripples_ts{nTT+1} = [ripples.peak_ts];

%%
N = cellfun(@length, ripples_ts);
N_joint = nan(length(ripples_ts),length(ripples_ts));
for ii = 1:length(ripples_ts)
    for jj = 1:length(ripples_ts)
        %%
        ts1 = ripples_ts{ii};
        ts2 = ripples_ts{jj};
        if isempty(ts1) || isempty(ts2)
            continue;
        end
        tdiff=abs(ts1-ts2');
        N_joint(ii,jj) = sum(any(tdiff<max_tdiff,2));
    end
end

%%
figure('Units','centimeters','Position',[5 5 25 12]);
subplot(121)
h=heatmap(N_joint);
h.Title='Counts';
h=gca;
h.XDisplayLabels{end} = 'all';
h.YDisplayLabels{end} = 'all';
xlabel('Reference detection TT')
ylabel('Detection TT')
subplot(122)
h=heatmap(round(100.*N_joint./N));
h.Title='Percentage';
h=gca;
h.XDisplayLabels{end} = 'all';
h.YDisplayLabels{end} = 'all';
xlabel('Reference detection TT')
ylabel('Detection TT')
h=sgtitle({exp_ID,'Cross-tetordes ripples detection'});
h.Interpreter='none';
% h.Colormap = flipud(h.Colormap)
% colormap bone
TT_pos_str = "TT"+[1:nTT]+": "+exp.details.TT_loc;
annotation('textbox',[.9 .3 .4 .5],'String',TT_pos_str,'FitBoxToText','on')

%% save figure
fig_filename = fullfile('L:\Analysis\Results\exp\ripples', [exp_ID '_ripples_ts_corr_matrix']);
saveas(gcf,fig_filename,'jpeg')

%%







% % %% bin
% % bin_size = 10e3; % in usec
% % win = 1e6;
% % win = 500e3;
% % tstart = min([ripples_ts{:}]) - win;
% % tend = max([ripples_ts{:}])   + win;
% % edges = tstart:bin_size:tend;
% % N=[];
% % for ii = 1:length(ripples_ts)
% %     [N(ii,:),edges] = histcounts(ripples_ts{ii},edges);
% % end
% % 
% % %% smooth N to overcome small time diffs
% % N = smoothdata(N,2,'movmean',5);
% % 
% % %% cross-correlation
% % [ccc,lags] = xcorr(N', round(win/bin_size));
% % lags_t = lags .* bin_size * 1e-3;
% % ccc = reshape(ccc,size(ccc,1), size(N,1), size(N,1));
% % 
% % %% plot
% % figure
% % pnl=panel();
% % pnl.pack(size(ccc,2),size(ccc,3));
% % 
% % for ii = 1:size(ccc,2)
% %     for jj = 1:size(ccc,3)
% %         pnl(ii,jj).select()
% %         plot(lags_t,ccc(:,ii,jj))
% %         yticks([0 max(ccc(:,ii,jj))])
% %     end
% % end






%%











%%
