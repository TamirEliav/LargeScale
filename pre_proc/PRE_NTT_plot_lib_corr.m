function PRE_NTT_plot_lib_corr(main_dir)

files = dir(fullfile(main_dir, '**', '*.NTT'));
for ii_file = 1:length(files)
    file_IN = fullfile(  files(ii_file).folder, files(ii_file).name );
    file_OUT = fullfile( files(ii_file).folder, strrep(files(ii_file).name, '.NTT','_lib_corr_by_units') );
    gen_lib_corr_figure_all_units(file_IN, file_OUT);
end
end


%%
function gen_lib_corr_figure_all_units(file_IN, file_OUT)

%% load NTT data
[Timestamps, CellNumbers, wvfrms, Header] = ...
     Nlx2MatSpike(file_IN, [1 0 1 0 1], 1, 1, [] );
ADBitVolts = sscanf(Header{contains(Header,'ADBitVolts')},'-ADBitVolts %f %f %f %f'); % parse header
wvfrms = wvfrms .* ADBitVolts' .* 1e6; % convert bits to uVolts

%% load library
% read params struct (should be in the same folder as the NTT file)
[FILEPATH,NAME,EXT] = fileparts(file_IN);
params_file = fullfile(FILEPATH, 'params.mat');
if exist(params_file,'file')
    load(params_file);
else
    % set defaults
    params.lib_spike_shapes = 'library_of_acceptable_spike_shapes_new.mat';
    params.lib_corr_thr = 0.9;
    params.AlignSample = 8;
end
load(params.lib_spike_shapes)

%% calculate lib corr
% For each event take the channel with the largest peak (8th point)
[~,max_ch_IX] = max(squeeze(abs(wvfrms(params.AlignSample,:,:))),[],1);
largest_waveforms = zeros(size(wvfrms,1),size(wvfrms,3));
for ch=1:4
    IX = find(max_ch_IX == ch);
    largest_waveforms(:,IX) = squeeze(wvfrms(:,ch,IX));
end
% now, let's calc the max corr
xxx_lags_shifts = [1:30; 2:31; 3:32];
ccc = [];
rrr = [];
for ii_shift = 1:size(xxx_lags_shifts,1)
    xxx_lags = xxx_lags_shifts(ii_shift, :);
    ccc = corr(largest_waveforms(xxx_lags,:), library_of_acceptable_spike_shapes(:,2:end-1)');
    rrr(ii_shift,:) = max(ccc,[],2);
end
rrr = max(rrr,[],1);

%% generate figure
CellNumbers_list = unique(CellNumbers);
nUnits = length(CellNumbers_list);
figure('Units','normalized','Position',[0 0 1 1]);
max_col = 4;
ncol = min(nUnits,max_col);
nrow = ceil(nUnits/max_col);
for ii_unit = 1:nUnits
    subplot(nrow,ncol,ii_unit);
    hold on
    unitID = CellNumbers_list(ii_unit);
    h=histogram(rrr(CellNumbers==unitID));
    plot(repelem(params.lib_corr_thr,2), get(gca,'ylim'),'m--');
    xlabel('max lib corr for all waveforms')
    ylabel('Counts')
    title(['unit ' num2str(unitID)])
end
h=suptitle(file_IN);
h.FontSize=14;
saveas(gcf, file_OUT, 'tif');
close(gcf)

end






