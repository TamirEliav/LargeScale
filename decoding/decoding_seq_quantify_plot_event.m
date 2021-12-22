function decoding_seq_quantify_plot_event(exp_ID, epoch_type, params_opt, event_type,event_num,margin,clim_prc)
arguments
    %% 
    exp_ID = 'b9861_d180526'
    epoch_type {mustBeMember(epoch_type,{'sleep','rest','flight'})} = 'sleep'
    params_opt = 11;
    event_type {mustBeMember(event_type,{'PE','posterior','ripples','MUA'})} = 'posterior'
    event_num = 49
    margin = 0.2
    clim_prc = 1
end

%% IN/OUT folders
dir_OUT = 'F:\sequences\events_quantification\events_seq_radon';
mkdir(dir_OUT);

%% load data
[events, params] = decoding_load_events_quantification(exp_ID, epoch_type, params_opt, event_type);
decode = decoding_load_data(exp_ID, epoch_type, params_opt, 'load_likelihood',True);

%% arrange data
event = events(event_num);
radius_bins = round(params.radius / params.decode.pos_bin_size);
margin_bins = round(decode.Fs * margin);
IX1 = [event.start_IX:event.end_IX];
IX2 = [(event.start_IX-margin_bins):(event.end_IX+margin_bins)];
seq_model = event.seq_model; seq_model.decoder_type_str = "model"; seq_model.decoder_type_field = 'posterior';
seq_bayes = event.seq_bayes; seq_bayes.decoder_type_str = "bayes"; seq_bayes.decoder_type_field = 'likelihood';
seqs = [seq_model seq_bayes];

%% plot
fig=figure;
fig.WindowState = 'maximized';
pnl = panel();
pnl.pack('h',2,'v',[0.06 0.06 0.35 0.35 0.08]);
pnl.de.margin=5;
pnl.de.marginbottom=0;
pnl.de.margintop=3;
pnl.de.marginleft=20;
pnl.margin=[20 0 20 15];
pnl(1,4).marginbottom = 10;
pnl(2,4).marginbottom = 10;
pnl(1,3).marginbottom = 15;
pnl(2,3).marginbottom = 15;
h=pnl.title(sprintf('%s - %s - opt %d - %s - event %d',...
                    exp_ID, epoch_type, params_opt, event_type,event_num));
h.Interpreter='none';
h.FontSize = 14;
for ii_decoder = 1:length(seqs)
    % arrange data for plotting
    seq = seqs(ii_decoder);
    states_prob = squeeze(sum(decode.posterior(:,:,IX2),1));
    fn = seq.decoder_type_field;
    prob1 = squeeze(decode.(fn)(:,event.state_num,IX1));
    prob2 = squeeze(decode.(fn)(:,event.state_num,IX2));
%     prob2 = prob2 - median(prob1,'all');
%     u = ones(2*radius_bins+1,1);
%     v=1;
%     prob2 = conv2(u,v, prob2,'same');
    time1 = decode.time(IX1);
    time2 = decode.time(IX2);
    
    % plot state prob
    pnl(ii_decoder,1).select();
    hold on
    plot(time2,states_prob','LineWidth',1.5);
    xline(event.start_ts);
    xline(event.end_ts);
    rescale_plot_data('x',[1e-6 event.start_ts]);
%     xlabel('Time (s)');
    ylabel('Prob.');
    xticklabels([])
    h=legend(decode.state,'NumColumns',round(length(decode.state)/3),'Location','southoutside');
    h.Position([1 2]) = [0.75 0.94];
    h.FontSize=8;
    h.Interpreter = 'none';
    h.Box = 'on';
    
    % plot confidence
    pnl(ii_decoder,2).select();
    hold on
    [KS, HPD, sparsity] = decoding_calc_prob_confidence(prob2);
    plot(time2, sparsity,'LineWidth',2);
    plot(time2, KS,'LineWidth',2);
%     plot(time2, HPD);
    legend("sparsity","KS")
    ylabel('Conf.');
    rescale_plot_data('x',[1e-6 event.start_ts]);
    ylim([0 1.1])
    xticklabels([])
    
    % prob image
    pnl(ii_decoder,3).select();
    hold on
    imagesc(prob2,'XData',time2,'YData',decode.pos);
    xline(event.start_ts);
    xline(event.end_ts);
    xx = interp1(1:length(time1),time1,seq.xx,'linear','extrap');
    yy = interp1(1:length(decode.pos),decode.pos,seq.yy,'linear','extrap');
    plot(xx,yy,'r');
    rescale_plot_data('x',[1e-6 event.start_ts]);
    axis tight
    xlabel('Time (s)');
    ylabel('Position (m)');
    cmap = bone;
    cmap = flipud(cmap);
    colormap(gca,cmap);
    climits = prctile(prob1(:), [0 100]+[1 -1].*clim_prc);
    set(gca,'CLim', climits);
    set(gca,'TickDir','out')
    
    % Radon
    pnl(ii_decoder,4).select();
    hold on
    u = ones(2*radius_bins+1,1);
    v=1;
%     prob1 = conv2(u,v, prob1,'same');
%     prob1 = prob1 - median(prob1,'all');
    theta = linspace(0,180,1000);
    [R, xp] = radon(prob1, theta);
    imagesc('CData',R,'XData',theta,'YData',xp);
%     xline(seq.theta,'--r',seq.theta)
%     yline(seq.offset,'--r',seq.offset)
    axis tight
    colormap(gca,parula);
    xlabel('\theta (degrees)')
    ylabel('x''')
    set(gca,'TickDir','out')
    
    % Details
    pnl(ii_decoder,5).select();
    axis off
    seq_details_str = {"\theta = "+seq.theta;
                       "offset = "+seq.offset;
                       "speed = "+seq.speed;
                       "compression = "+seq.compression;
                       "confidence (sparsity) = "+seq.confidence_sparsity;
                       "confidence (KS) = "+seq.confidence_KS;
                       };
    text(0,1,seq_details_str,'HorizontalAlignment','left','VerticalAlignment','top');
end
linkaxes([pnl(1,1).axis pnl(1,2).axis pnl(1,3).axis
          pnl(2,1).axis pnl(2,2).axis pnl(2,3).axis],'x');
% save figure
figname = sprintf('%s_%s_opt_%d_%s_event_%d',...
                    exp_ID, epoch_type, params_opt, event_type,event_num);
filename = fullfile(dir_OUT, figname);
% saveas(fig, filename, 'jpg');
% close(fig);

%%






