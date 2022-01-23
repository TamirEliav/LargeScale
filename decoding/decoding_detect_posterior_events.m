function decoding_detect_posterior_events(decode)
%%
exp_ID = decode.exp_ID;
epoch_type = decode.epoch_type;
params_opt = decode.params_opt;

%% load data
fprintf('Detecting posterior events for data:\n')
fprintf('exp: %s, epoch type: %s, decoding paramset: %d\n', exp_ID, epoch_type, params_opt);
exp = exp_load_data(exp_ID, 'details','path','rest');
dir_IN = 'F:\sequences\decoded';
dir_OUT = 'F:\sequences\posterior_events';
figs_dir = fullfile(dir_OUT,'figs');
mkdir(dir_OUT);
mkdir(figs_dir);

%% get movement states
mvmnt_states_IX = find(contains(decode.state, "empirical_movement"));
events_all=struct();
for ii_state = 1:length(mvmnt_states_IX)
    %%
    state_num = mvmnt_states_IX(ii_state);
    state_prob = decode.posterior_state(state_num,:);
    state_str = decode.state(state_num);
    [events,g,opts] = detect_events(state_prob,...
                            "Fs",decode.Fs,...
                            "thr",0.8,...
                            "thr_edges",0.5,...
                            "merge_thr",0.1,...
                            "min_width",0.05,...
                            "plot",true);
	details_str = {
        "#events = "+length(events);...
        "thr = " + opts.thr;...
        "thr edges = " + opts.thr_edges;...
        "merge thr = " + opts.merge_thr;...
        "min width = " + opts.min_width + "s";...
        "";...
        sprintf('Decoding params opt: %d',params_opt),
        sprintf('bin size: %.2gm',decode.params.pos_bin_size),
        sprintf('replay speed: x%d',decode.params.replay_speed),
        sprintf('state decay timescale: %.3g s',decode.params.state_decay_timescale),
        };
    text(1, 1, details_str ,'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top');
    gaps_IX = find(diff(decode.time) > median(diff(decode.time)));
    for ii_gap = 1:length(gaps_IX)
        hl = xline(gaps_IX(ii_gap),'g','time gap');
        hl.LineWidth = 2;
    end
    rescale_plot_data('x',[1/decode.Fs 0]);
    ylim([0 1]);
    xlabel('Time (s)');
    ylabel('State probability');
    title(  {exp_ID;...
            'State high probability events detection';...
            sprintf('%s (%s)',state_str,epoch_type)},...
          'Interpreter','none');
    fig=gcf;
    fig.WindowState = 'maximized';
    mkdir(figs_dir, 'detect')
    filename = fullfile(figs_dir, 'detect', sprintf('%s_posterior_events_detect_%s_dec_prm_%d_%s',exp_ID,epoch_type,params_opt,state_str));
    saveas(gcf, filename, 'jpg');
    close(gcf);
    
    %% calc some events features
    sym_index = @(x)(corr(x',flip(x')));
    [events.peak_ts] = disperse(decode.time([events.peak_IX]));
    [events.start_ts] = disperse(decode.time([events.start_IX]));
    [events.end_ts] = disperse(decode.time([events.end_IX]));
	[events.mean_prob] = disperse(splitapply(@mean,state_prob,g));
    [events.symmetry] = disperse(splitapply(sym_index,state_prob,g));
    [~,max_IX] = splitapply(@max,state_prob,g);
    lengths  = splitapply(@length,state_prob,g);
    [events.peak_relative_loc] = disperse(max_IX./lengths);

    %% arrange results in a struct
    events_all(ii_state).state_str = state_str;
    events_all(ii_state).state_num = state_num;
    events_all(ii_state).events = events;
    events_all(ii_state).opts = opts;
    
end

%% save all events to mat file
filename = fullfile(dir_OUT, sprintf('%s_posterior_events_%s_dec_prm_%d',exp_ID,epoch_type,params_opt));
save(filename, 'events_all');


end

