function exp_calc_pos_std_y(exp_ID)

%% load exp data
exp = exp_load_data(exp_ID);
prm = PARAMS_GetAll();

%% arrange relevant data for analysis (take only data from detected flights)
pos = exp.pos.proj.pos;
ts = exp.pos.proj.ts;
pos(exp.pos.proj.outliers_IX,:)=[];
ts(exp.pos.proj.outliers_IX)=[];
ti = [[exp.flight.FE.start_ts];[exp.flight.FE.end_ts]]';
IX = get_data_in_ti(ts',ti);
ts = ts(IX);
pos = pos(IX,:);
directions = interp1(exp.pos.proc_1D.ts, sign(exp.pos.proc_1D.vel_csaps), ts);

%%
range_prc = [25 75];
trim_prc = 5;
bin_edges = 0:1:200;
bin_centers = edges2centers(bin_edges);
dir_sign = [1 -1];
pos_y_std = repelem(struct(),2);
for ii_dir = 1:2
    IX = find(directions == dir_sign(ii_dir));
    [~,~,BIN] = histcounts(pos(IX,1),bin_edges);
    
    ymean = accumarray(BIN,pos(IX,2),[],@(x)(trimmean(x,trim_prc*2)));
    ystd = accumarray(BIN,pos(IX,2),[],@(x)(trimstd(x,trim_prc*2)));
    yrange(1,:) = accumarray(BIN,pos(IX,2),[],@(x)(prctile(x,range_prc(1))));
    yrange(2,:) = accumarray(BIN,pos(IX,2),[],@(x)(prctile(x,range_prc(2))));
    
    ymean_ = nan(size(bin_centers));
    ystd_ = nan(size(bin_centers));
    yrange_ = nan(2,length(bin_centers));
    ymean_(unique(BIN)) = ymean(unique(BIN));
    ystd_(unique(BIN)) = ystd(unique(BIN));
    yrange_(:,unique(BIN)) = yrange(:,unique(BIN));
    
    pos_y_std(ii_dir).xy = pos(IX,:);
    pos_y_std(ii_dir).bin_centers = bin_centers;
    pos_y_std(ii_dir).ymean = ymean_;
    pos_y_std(ii_dir).ystd = ystd_;
    pos_y_std(ii_dir).yrange = yrange_;
    pos_y_std(ii_dir).range_prc = range_prc;
    pos_y_std(ii_dir).trim_prc = trim_prc ;
end

%% save strcut
filename_exp_pos_y_std = ['L:\Analysis\Results\exp\pos_y_std\' exp_ID '_exp_pos_y_std' ];
save(filename_exp_pos_y_std, 'pos_y_std');







end

