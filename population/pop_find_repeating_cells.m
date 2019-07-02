%%
clear
clc

%%
dir_out = 'L:\Analysis\Results\cells\repeating_cells';
mkdir(dir_out);

%% load cells
cells_t = DS_get_cells_summary();
cells_t(~ismember(cells_t.bat, [79,148,34,9861,2289] ),:) = [];
prm = PARAMS_GetAll();

%% filter cells - brain region / sorting quality
cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells = [cells.details];
cells(~contains({cells.brain_area}, 'CA1')) = [];
cells(~ismember([cells.ClusterQuality], [2])) = [];
% cells(cellfun(@isempty, {cells.stable_ts})) = [];
cells_t = cells_t({cells.cell_ID},:);

%% filter cells - mean FR
cells = cellfun(@(c)(cell_load_data(c,'stats')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells = [cells.stats];
cells = [cells.all];
cells_t([cells.meanFR_all]>prm.inclusion.interneuron_FR_thr,:) = [];

%% disp final cells table
% cells_t
whos cells_t


%% load cells FR maps
cells = cellfun(@(c)(cell_load_data(c,'FR_map')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells = arrayfun(@(c)([c.FR_map(1).all.PSTH c.FR_map(2).all.PSTH]), cells, 'UniformOutput',0);
FR_maps_all = cat(1,cells{:});
% calc map correlations
FR_maps_all(isnan(FR_maps_all)) = 0;
maps_corr = corr(FR_maps_all');

%% load cells details (mask only same bat+TT)
cells = cellfun(@(c)(cell_load_data(c,'details')), cells_t.cell_ID, 'UniformOutput',0);
cells = [cells{:}];
cells = [cells.details];
[X,Y] = meshgrid([cells.bat]);
same_bat_mask = (X==Y);
[X,Y] = meshgrid([cells.TT]);
same_TT_mask = (X==Y);
same_bat_TT_mask = same_TT_mask & same_bat_mask;

%% create null distribution for setting a threshold
mask = tril(~same_bat_TT_mask,-1);
maps_corr_diff_TT = maps_corr(mask);
mask = tril(same_bat_TT_mask,-1);
maps_corr_same_TT = maps_corr(mask);
figure('units','normalized','outerposition',[0 0 1 1])
hold on
clc
histogram(maps_corr_diff_TT,'Normalization','pdf')
histogram(maps_corr_same_TT,'Normalization','pdf')
p = [0.05 0.01 0.001];
p_quantiles = quantile(maps_corr_diff_TT, 1-p);
ylimits = get(gca,'ylim');
plot(repelem(p_quantiles',1,2), ylimits, '--m')
text(p_quantiles, repelem(ylimits(1)+0.95*range(ylimits),length(p_quantiles)), "" + p_quantiles)
legend({'different TT';'same TT'},'Location','eastoutside');
xlabel('Map Correlation')
ylabel('pdf')
title('cell pairs map correlations distribution')
figname = fullfile(dir_out, 'map_corr_hist');
saveas(gcf,figname, 'tif')
saveas(gcf,figname, 'fig')
ylim([0 2])
figname = fullfile(dir_out, 'map_corr_hist_ylim_cut');
saveas(gcf,figname, 'tif')
saveas(gcf,figname, 'fig')

%% plot corr matrix
figure
corr_thr = 0.7;
imagesc( same_bat_TT_mask.*tril(maps_corr,-1) )
colorbar
set(gca,'CLim',[corr_thr 1])

%% create figures of putative same cell pair
corr_thr_plot = 0.5;
mask = tril(same_bat_TT_mask,-1);
maps_corr_same_TT = maps_corr .* mask;
[I J] = find( maps_corr_same_TT > corr_thr_plot  );
whos I J
cell_IDs = {cells.cell_ID};
pairs_IX = [I J];
pairs = cell_IDs(pairs_IX);
for ii_pair = 1:size(pairs,1)
    %% re-order pair by date
    cell1 = cell_load_data(pairs{ii_pair,1},'details');
    cell2 = cell_load_data(pairs{ii_pair,2},'details');
    [~,sort_IX] = sort([cell1.details.date cell2.details.date]','ascend');
    pairs(ii_pair,:) = pairs(ii_pair,sort_IX);
    cell_ID1 = pairs{ii_pair,1};
    cell_ID2 = pairs{ii_pair,2};
    pair_map_corr = maps_corr_same_TT(I(ii_pair), J(ii_pair));
    fprintf('%3d:\tcorr=%.2f\t%s vs. %s\n', ii_pair , pair_map_corr , cell_ID1 , cell_ID2);
%     continue
    %% create figure with panel
    figure('units','normalized','outerposition',[.1 .1 .8 .7])
    pnl = panel();
    pnl.pack(2,2,2);
    pnl.margin = 20;
    pnl(1).de.margin = 5;
    pnl(2).de.margin = 5;
%     pnl.select('all')
%     pnl.identify()
    
    %% 
    depths = [];
    for ii_cell = 1:2
        cell_ID = pairs{ii_pair,ii_cell};
        cell = cell_load_data(cell_ID,'details','FR_map','FE','fields');
        depths(ii_cell) = cell.details.depth;
        ti = exp_get_sessions_ti(cell.details.exp_ID, 'Behave');
        for ii_dir = 1:2
            dir_color = prm.graphics.colors.flight_directions{ii_dir};
            % FR map
            pnl(ii_cell,ii_dir,1).select(); hold on
            plot(cell.FR_map(ii_dir).all.bin_centers, cell.FR_map(ii_dir).all.PSTH,...
                'Color', dir_color);
            % fields
            fields = cell.fields{ii_dir};
            for ii_field = 1:length(fields)
                field = fields(ii_field);
                if field.in_low_speed_area
                    field_color = 'r';
                else
                    field_color = 'k';
                end
                plot(field.loc, field.peak, '*', 'MarkerSize', 4, 'Color', field_color )
                plot(field.edges_href, repelem(prm.fields.width_href * field.peak,2), 'Color', field_color )
                text(field.loc, field.peak,...
                    sprintf('{\\color{magenta}%d}\n{\\color{red}%2.0f%%}\n%2.2f\n%2.2f',....
                    length(field.spikes_ts),...
                    field.num_flights_with_spikes_prc*100,...
                    field.width_href,...
                    field.width_prc),...
                    'FontSize',5,'HorizontalAlignment','center','VerticalAlignment','bottom');
            end
            xlim([0 200])
            % trajectory + spikes
            pnl(ii_cell,ii_dir,2).select(); hold on
            plot([cell.FE{ii_dir}.pos], [cell.FE{ii_dir}.ts], '.', 'color', 0.9*[1 1 1], 'MarkerSize',0.5);
            plot([cell.FE{ii_dir}.spikes_pos], [cell.FE{ii_dir}.spikes_ts], '.', 'Color', dir_color);
            rescale_plot_data('y',[1e-6/60 ti(1)]);
            n_flights = sum([cell.FE{ii_dir}.distance] > prm.flight.full_min_distance);
            text(0.99,0.9,num2str(n_flights), 'Units','normalized','FontWeight','bold','HorizontalAlignment','right');
            xlim([0 200])
        end
        h=pnl(ii_cell).title(sprintf('%s(%d)',cell_ID,cell.details.cell_num));
        h.Position = [0.5 1];
        h.Interpreter='none';
        h.FontSize = 14;
    end
    h=pnl.title(sprintf('corr=%.2f, depth_diff=%d', pair_map_corr, diff(depths)));
    h.Position = [0.95 1.1];
    h.Interpreter='none';
    h.FontSize = 12;
    figname = sprintf('corr_%.2f_diff_depth_%dum_%s(%d)_vs_%s(%d)', ...
        pair_map_corr, diff(depths),...
        cell_ID1, cell1.details.cell_num,...
        cell_ID2, cell2.details.cell_num);
    figname = strrep(figname, '.', '_');
    figname = fullfile(dir_out, figname);
    saveas(gcf, figname, 'tif')
    close(gcf)
end


%% copy important pair figures
list_excitability_changes = {
{'corr_0_81_diff_depth_20um_b0148_d170625_TT4_SS02(694)_vs_b0148_d170626_TT4_SS01(679).tif' }
{'corr_0_84_diff_depth_-40um_b0079_d160927_TT2_SS01(441)_vs_b0079_d160929_TT2_SS01(429).tif'}
{'corr_0_86_diff_depth_-20um_b0079_d160929_TT2_SS01(446)_vs_b0079_d160930_TT2_SS01(441).tif'}
{'corr_0_89_diff_depth_-40um_b0079_d160928_TT2_SS01(446)_vs_b0079_d160930_TT2_SS01(433).tif'}
{'corr_0_89_diff_depth_-60um_b0079_d160927_TT2_SS01(447)_vs_b0079_d160930_TT2_SS02(429).tif'}
{'corr_0_90_diff_depth_-60um_b0079_d160927_TT2_SS01(446)_vs_b0079_d160930_TT2_SS01(429).tif'}
{'corr_0_90_diff_depth_20um_b0034_d180313_TT4_SS02(88)_vs_b0034_d180314_TT4_SS03(73).tif'   }
{'corr_0_90_diff_depth_40um_b0034_d180312_TT4_SS01(88)_vs_b0034_d180314_TT4_SS03(56).tif'   }
{'corr_0_91_diff_depth_-20um_b0079_d160928_TT2_SS01(441)_vs_b0079_d160929_TT2_SS01(433).tif'}
{'corr_0_94_diff_depth_-20um_b0079_d160929_TT2_SS01(447)_vs_b0079_d160930_TT2_SS02(441).tif'}
{'corr_0_94_diff_depth_-40um_b0079_d160928_TT2_SS01(447)_vs_b0079_d160930_TT2_SS02(433).tif'}
};

list_not_sure = {
{'corr_0_55_diff_depth_0um_b9861_d180523_TT1_SS01(162)_vs_b9861_d180524_TT1_SS01(141).tif'   }
{'corr_0_56_diff_depth_60um_b0148_d170626_TT4_SS01(702)_vs_b0148_d170627_TT4_SS01(694).tif'  }
{'corr_0_66_diff_depth_20um_b0034_d180313_TT4_SS03(86)_vs_b0034_d180314_TT4_SS01(74).tif'    }
{'corr_0_69_diff_depth_20um_b0034_d180312_TT4_SS02(72)_vs_b0034_d180313_TT4_SS01(57).tif'    }
{'corr_0_70_diff_depth_40um_b0148_d170619_TT4_SS01(658)_vs_b0148_d170620_TT4_SS02(651).tif'  }
{'corr_0_84_diff_depth_80um_b9861_d180524_TT3_SS05(180)_vs_b9861_d180525_TT3_SS01(176).tif'  }
{'corr_0_85_diff_depth_0um_b9861_d180523_TT1_SS04(164)_vs_b9861_d180524_TT1_SS03(144).tif'   }
{'corr_0_85_diff_depth_40um_b0148_d170614_TT1_SS01(635)_vs_b0148_d170615_TT1_SS01(628).tif'  }
{'corr_0_86_diff_depth_120um_b0148_d170801_TT4_SS08(514)_vs_b0148_d170806_TT4_SS02(478).tif' }
{'corr_0_86_diff_depth_40um_b0148_d170802_TT4_SS05(499)_vs_b0148_d170803_TT4_SS04(485).tif'  }
{'corr_0_86_diff_depth_40um_b0148_d170803_TT4_SS08(489)_vs_b0148_d170806_TT4_SS02(478).tif'  }
{'corr_0_87_diff_depth_40um_b0034_d180312_TT3_SS03(85)_vs_b0034_d180314_TT3_SS02(52).tif'    }
{'corr_0_88_diff_depth_40um_b0148_d170802_TT4_SS04(498)_vs_b0148_d170803_TT4_SS08(489).tif'  }
{'corr_0_88_diff_depth_40um_b0148_d170806_TT4_SS01(477)_vs_b0148_d170807_TT4_SS04(474).tif'  }
{'corr_0_89_diff_depth_80um_b0148_d170801_TT4_SS08(514)_vs_b0148_d170803_TT4_SS08(489).tif'  }
{'corr_0_91_diff_depth_-90um_b0148_d170720_TT4_SS04(527)_vs_b0148_d170801_TT4_SS02(508).tif' }
{'corr_0_92_diff_depth_-180um_b0079_d160920_TT2_SS03(419)_vs_b0079_d160925_TT2_SS03(401).tif'}
{'corr_0_92_diff_depth_40um_b0148_d170606_TT4_SS10(592)_vs_b0148_d170607_TT4_SS05(583).tif'  }
{'corr_0_92_diff_depth_80um_b9861_d180529_TT1_SS01(266)_vs_b9861_d180601_TT1_SS01(247).tif'  }
{'corr_0_94_diff_depth_40um_b0148_d170801_TT4_SS08(514)_vs_b0148_d170802_TT4_SS04(498).tif'  }
{'corr_0_94_diff_depth_80um_b0148_d170802_TT4_SS04(498)_vs_b0148_d170806_TT4_SS02(478).tif'  }
{'corr_0_96_diff_depth_-200um_b0079_d160921_TT2_SS02(419)_vs_b0079_d160925_TT2_SS03(405).tif'}
};

list_same = {
% {'SAME_DAY!!!__need_to_merge!!__corr_0_91_diff_depth_0um_b0079_d160930_TT2_SS02(447)_vs_b0079_d160930_TT2_SS01(446).tif'}
% {'SAME_DAY!!!_consider_merging__corr_0_55_diff_depth_0um_b9861_d180526_TT3_SS04(212)_vs_b9861_d180526_TT3_SS02(210).tif'}
{'corr_0_91_diff_depth_0um_b0079_d160930_TT2_SS02(447)_vs_b0079_d160930_TT2_SS01(446).tif'          }
{'corr_0_55_diff_depth_0um_b9861_d180526_TT3_SS04(212)_vs_b9861_d180526_TT3_SS02(210).tif'          }
{'corr_0_81_diff_depth_20um_b0148_d170625_TT4_SS02(694)_vs_b0148_d170626_TT4_SS01(679).tif'                             }
{'corr_0_84_diff_depth_-40um_b0079_d160927_TT2_SS01(441)_vs_b0079_d160929_TT2_SS01(429).tif'                            }
{'corr_0_86_diff_depth_-20um_b0079_d160929_TT2_SS01(446)_vs_b0079_d160930_TT2_SS01(441).tif'                            }
{'corr_0_89_diff_depth_-40um_b0079_d160928_TT2_SS01(446)_vs_b0079_d160930_TT2_SS01(433).tif'                            }
{'corr_0_89_diff_depth_-60um_b0079_d160927_TT2_SS01(447)_vs_b0079_d160930_TT2_SS02(429).tif'                            }
{'corr_0_90_diff_depth_-60um_b0079_d160927_TT2_SS01(446)_vs_b0079_d160930_TT2_SS01(429).tif'                            }
{'corr_0_90_diff_depth_20um_b0034_d180313_TT4_SS02(88)_vs_b0034_d180314_TT4_SS03(73).tif'                               }
{'corr_0_90_diff_depth_40um_b0034_d180312_TT4_SS01(88)_vs_b0034_d180314_TT4_SS03(56).tif'                               }
{'corr_0_91_diff_depth_-20um_b0079_d160928_TT2_SS01(441)_vs_b0079_d160929_TT2_SS01(433).tif'                            }
{'corr_0_94_diff_depth_-40um_b0079_d160928_TT2_SS01(447)_vs_b0079_d160930_TT2_SS02(433).tif'                            }
{'corr_0_95_diff_depth_-20um_b0079_d160927_TT2_SS01(433)_vs_b0079_d160928_TT2_SS01(429).tif'                            }
{'corr_0_95_diff_depth_0um_b9861_d180523_TT1_SS03(163)_vs_b9861_d180524_TT1_SS02(143).tif'                              }
{'corr_0_95_diff_depth_20um_b0034_d180312_TT3_SS03(66)_vs_b0034_d180313_TT3_SS02(52).tif'                               }
{'corr_0_95_diff_depth_20um_b0034_d180313_TT3_SS02(85)_vs_b0034_d180314_TT3_SS02(66).tif'                               }
{'corr_0_95_diff_depth_20um_b0079_d160920_TT2_SS03(405)_vs_b0079_d160921_TT2_SS02(401).tif'                             }
{'corr_0_95_diff_depth_40um_b0148_d170718_TT4_SS03(533)_vs_b0148_d170720_TT4_SS05(528).tif'                             }
{'corr_0_96_diff_depth_20um_b0034_d180312_TT4_SS01(73)_vs_b0034_d180313_TT4_SS02(56).tif'                               }
};

list_excitability_changes = [list_excitability_changes{:}]';
list_not_sure = [list_not_sure{:}]';
list_same = [list_same{:}]';

util_copy_files_by_template(dir_out, fullfile(dir_out,'excitability_changes'),  list_excitability_changes);
util_copy_files_by_template(dir_out, fullfile(dir_out,'not_sure'),              list_not_sure);
util_copy_files_by_template(dir_out, fullfile(dir_out,'same'),                  list_same);




%%









%%









