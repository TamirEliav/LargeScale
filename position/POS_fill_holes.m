function [bsp_full_data_pos,bsp_full_data_ts] = POS_fill_holes(raw_data_pos,raw_data_ts_usec)
% raw_data - raw bsp data from one flight: with holes that we need to fill:
%           raw_data_pos (1xN) - linearized bsp position
%           raw_data_ts_usec (1xN) - timestamps in usec
% bspFreq  - Bsp sampling frquency for this data in Hz.

interpolation_always_TH = (1/3)*10^6; %usec, max hole size for which always do linear interpolation
max_extrapolation_length = (1/3)*10^6; %usec
TH_velocity_index = 0.15;
large_hole_TH = (1.5)*10^6; %usec, max hole size for interpolation. For
% holes larger than this do extrapolation instead.
time_for_fitting_inter = (1/3)*10^6; %~5-6 samples to use from each side for the linear fitting
time_for_fitting_extra = (1/2)*10^6; %~9 samples to use from each side

diff_ts = diff(raw_data_ts_usec);
min_ts_diff = median(diff_ts) + 1; %the regular dt between consecutive samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off', 'curvefit:fit:equationBadlyConditioned')
%% remove repetitions in data: (consecutive samples with the exact same position)

vel_TH = 2; %m/s , so we will only remove repeting smaples in the middle of
% the flight, and not from the long time on the balls.
vel = [0 diff(raw_data_pos)./(diff(raw_data_ts_usec)*10^-6)];
vel_smooth_win = 2/(min_ts_diff*10^-6); %the number of samples in 2 sec
vel = smooth(vel,vel_smooth_win)';
repetitions_ind = find(and(vel>=vel_TH,diff([0 raw_data_pos])==0));
raw_data_pos(repetitions_ind+1) = [];
raw_data_ts_usec(repetitions_ind+1) = [];

%% Find holes:
diff_ts = diff(raw_data_ts_usec);
holes_ind_start = find(diff_ts > min_ts_diff);
holes_ind_end = holes_ind_start + 1;
holes_size_t = raw_data_ts_usec(holes_ind_end) - raw_data_ts_usec(holes_ind_start); %usec
%% Merge holes in consecutive indices:
% decided to solve this issues with single sample within a hole in a different way, without merging holes.
% % % % tmp1 = diff([0, holes_ind_start])>1; % 0 for holes that comes immediately after a
% % % % % hole (only one sample in between the two holes), and 1 otherwise.
% % % % holes_ind_start = holes_ind_start(tmp1);
% % % % tmp2 = find(diff(tmp1)~=0); %the transition in tmp1 from '1's to '0's and back
% % % % tmp3 = diff(tmp2); %the indices diff, so how many indices are added to the merged hole
% % % % holes_ind_end(tmp2(1:2:end)) = holes_ind_end(tmp2(1:2:end)) + tmp3(1:2:end);
% % % % holes_ind_end = holes_ind_end(tmp1);
% % % % holes_size_t = raw_data_ts_usec(holes_ind_end) - raw_data_ts_usec(holes_ind_start); %usec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill holes:


if isempty(holes_ind_start)
    bsp_full_data_pos = raw_data_pos;
    bsp_full_data_ts = raw_data_ts_usec;
    warning('off', 'curvefit:fit:equationBadlyConditioned')
    return;
end

samples_to_fill_eachside_extrapolation = round(max_extrapolation_length/min_ts_diff);

bsp_full_data_pos(1:holes_ind_start(1)) = raw_data_pos(1:holes_ind_start(1));
bsp_full_data_ts(1:holes_ind_start(1)) = raw_data_ts_usec(1:holes_ind_start(1));

hole_i = 1;
while hole_i <= length(holes_ind_start)
    total_samples_inHole = round(holes_size_t(hole_i)/min_ts_diff) - 1; %total samples in hole
    fill_in_ts = raw_data_ts_usec(holes_ind_start(hole_i)) + ((1:total_samples_inHole)*min_ts_diff);

    side1_ind = find(abs(raw_data_ts_usec(holes_ind_start(hole_i))-raw_data_ts_usec(1:holes_ind_end(hole_i)-1))<=time_for_fitting_inter);
    side2_ind = holes_ind_start(hole_i) + find(abs(raw_data_ts_usec(holes_ind_end(hole_i))-raw_data_ts_usec(holes_ind_start(hole_i)+1:end))<=time_for_fitting_inter);

    if holes_size_t(hole_i)<=(interpolation_always_TH) %small hole - ALWAYS LINEAR INTERPOLATION

        ts_for_fitting = raw_data_ts_usec([side1_ind side2_ind]);
        pos_for_fitting = raw_data_pos([side1_ind side2_ind]);

        coeffs = fit(ts_for_fitting',pos_for_fitting','poly1');
        fill_in_pos = polyval([coeffs.p1 coeffs.p2],fill_in_ts);

        bsp_full_data_pos = [bsp_full_data_pos fill_in_pos];
        bsp_full_data_ts = [bsp_full_data_ts fill_in_ts];

        if hole_i<length(holes_ind_start)
            bsp_full_data_pos = [bsp_full_data_pos raw_data_pos((holes_ind_end(hole_i)):holes_ind_start(hole_i+1))];
            bsp_full_data_ts = [bsp_full_data_ts raw_data_ts_usec((holes_ind_end(hole_i)):holes_ind_start(hole_i+1))];
        end
        hole_i = hole_i + 1;

        continue;
    end

    %% Velocity conditions:
    flag_change_hole_i = 0;
    while length(side2_ind)<2 %single sample inside a bigger hole:
        side2_ind_tmp = [side2_ind (side2_ind+1)];
        ts_side2 = raw_data_ts_usec(side2_ind_tmp);
        if diff(ts_side2) < large_hole_TH
            side2_ind = side2_ind_tmp;
        else
            flag_change_hole_i = flag_change_hole_i + 1;
            side2_ind = holes_ind_end(hole_i+flag_change_hole_i) + ...
                find(abs(raw_data_ts_usec(holes_ind_end(hole_i+flag_change_hole_i))-...
                raw_data_ts_usec(holes_ind_start(hole_i)+1:end))<=time_for_fitting_inter);
        end
    end

    if flag_change_hole_i %then holes were merged
        holes_ind_end(hole_i) = holes_ind_end(hole_i+flag_change_hole_i);
        ind_toRemove = hole_i+(1:flag_change_hole_i);
        holes_ind_start(ind_toRemove) = [];
        holes_ind_end(ind_toRemove) = [];
        holes_size_t = raw_data_ts_usec(holes_ind_end) - raw_data_ts_usec(holes_ind_start); %usec
        total_samples_inHole = round(holes_size_t(hole_i)/min_ts_diff) - 1; %total samples in hole
        fill_in_ts = raw_data_ts_usec(holes_ind_start(hole_i)) + ((1:total_samples_inHole)*min_ts_diff);
    end

    if length(side1_ind)<2 %single sample inside a bigger hole: side1.
        % the treatment of the two sides is not symmetric in the code
        % because of the order of computations, i.e. the order of the while
        % loop (it is symmetric in the actual data filling).
        side1_ind_tmp = [side1_ind (side1_ind+1)];
        side1_ind = side1_ind_tmp;
    end

    ts_side1 = raw_data_ts_usec(side1_ind);
    pos_side1 = raw_data_pos(side1_ind);
    coeffs_side1 = fit(ts_side1',pos_side1','poly1');
    vel_side1 = coeffs_side1.p1; %m/usec

    ts_side2 = raw_data_ts_usec(side2_ind);
    pos_side2 = raw_data_pos(side2_ind);
    coeffs_side2 = fit(ts_side2',pos_side2','poly1');
    vel_side2 = coeffs_side2.p1; %m/usec

    vel_hole = (raw_data_pos(holes_ind_end(hole_i)) - raw_data_pos(holes_ind_start(hole_i)))/holes_size_t(hole_i);

    side1_vs_side2 = abs((vel_side1 - vel_side2)/(vel_side1 + vel_side2 + eps)); %maybe use different TH
    hole_vs_side1 = abs((vel_hole - vel_side1)/(vel_hole + vel_side1 + eps));
    hole_vs_side2 = abs((vel_hole - vel_side2)/(vel_hole + vel_side2 + eps));

    if max([side1_vs_side2, hole_vs_side1, hole_vs_side2])>TH_velocity_index
        % no interpolation or extrapolation
        middle_nans = fill_in_ts;
        middle_nans(:) = nan;

        if flag_change_hole_i %then keep original raw samples
            [~,ind_raw_samples] = min(abs(raw_data_ts_usec(holes_ind_start(hole_i)+1:...
                holes_ind_start(hole_i)+flag_change_hole_i)' - fill_in_ts),[],2);
            fill_in_ts(ind_raw_samples) = raw_data_ts_usec(holes_ind_start(hole_i)+1:...
                holes_ind_start(hole_i)+flag_change_hole_i);
            middle_nans(ind_raw_samples) = raw_data_pos(holes_ind_start(hole_i)+1:...
                holes_ind_start(hole_i)+flag_change_hole_i);

        end
        bsp_full_data_pos = [bsp_full_data_pos middle_nans];
        bsp_full_data_ts = [bsp_full_data_ts fill_in_ts];

        if hole_i<length(holes_ind_start)
            bsp_full_data_pos = [bsp_full_data_pos raw_data_pos((holes_ind_end(hole_i)):holes_ind_start(hole_i+1))];
            bsp_full_data_ts = [bsp_full_data_ts raw_data_ts_usec((holes_ind_end(hole_i)):holes_ind_start(hole_i+1))];
        end
        hole_i = hole_i + 1;
        continue
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if holes_size_t(hole_i) < large_hole_TH  %interpolation
        ts_for_fitting = raw_data_ts_usec([side1_ind side2_ind]);
        pos_for_fitting = raw_data_pos([side1_ind side2_ind]);

        coeffs = fit(ts_for_fitting',pos_for_fitting','poly1');
        fill_in_pos = polyval([coeffs.p1 coeffs.p2],fill_in_ts);

        if flag_change_hole_i
            [~,ind_raw_samples] = min(abs(raw_data_ts_usec(holes_ind_start(hole_i)+1:...
                holes_ind_start(hole_i)+flag_change_hole_i)' - fill_in_ts),[],2);
            fill_in_ts(ind_raw_samples) = raw_data_ts_usec(holes_ind_start(hole_i)+1:...
                holes_ind_start(hole_i)+flag_change_hole_i);
            fill_in_pos(ind_raw_samples) = raw_data_pos(holes_ind_start(hole_i)+1:...
                holes_ind_start(hole_i)+flag_change_hole_i);
        end

        bsp_full_data_pos = [bsp_full_data_pos fill_in_pos];
        bsp_full_data_ts = [bsp_full_data_ts fill_in_ts];

    else %extrapolation
        side1_ind_extra = find(abs(raw_data_ts_usec(holes_ind_start(hole_i))-raw_data_ts_usec(1:holes_ind_end(hole_i)-1))<=time_for_fitting_extra);
        side2_ind_extra = holes_ind_start(hole_i) + find(abs(raw_data_ts_usec(holes_ind_end(hole_i))-raw_data_ts_usec(holes_ind_start(hole_i)+1:end))<=time_for_fitting_extra);

        if length(side1_ind_extra)<2
            side1_ind_extra = side1_ind;
        end
        if length(side2_ind_extra)<2
            side2_ind_extra = side2_ind;
        end

        ts_side1 = raw_data_ts_usec(side1_ind_extra);
        pos_side1 = raw_data_pos(side1_ind_extra);
        coeffs_side1 = fit(ts_side1',pos_side1','poly1');

        ts_side2 = raw_data_ts_usec(side2_ind_extra);
        pos_side2 = raw_data_pos(side2_ind_extra);
        coeffs_side2 = fit(ts_side2',pos_side2','poly1');

        %number of samples to extrapolate from each side of the hole
        fill_in_ts_side1 = fill_in_ts(1:samples_to_fill_eachside_extrapolation);
        fill_in_pos_side1 = polyval([coeffs_side1.p1 coeffs_side1.p2],fill_in_ts_side1);

        fill_in_ts_side2 = fill_in_ts((end-samples_to_fill_eachside_extrapolation+1):end);
        fill_in_pos_side2 = polyval([coeffs_side2.p1 coeffs_side2.p2],fill_in_ts_side2);

        middle_nans = zeros([1,total_samples_inHole-(2*samples_to_fill_eachside_extrapolation)]);
        middle_nans(:) = nan;

        fill_in_pos = [fill_in_pos_side1 middle_nans fill_in_pos_side2];
        if flag_change_hole_i
            [~,ind_raw_samples] = min(abs(raw_data_ts_usec(holes_ind_start(hole_i)+1:...
                holes_ind_start(hole_i)+flag_change_hole_i)' - fill_in_ts),[],2);
            fill_in_ts(ind_raw_samples) = raw_data_ts_usec(holes_ind_start(hole_i)+1:...
                holes_ind_start(hole_i)+flag_change_hole_i);
            fill_in_pos(ind_raw_samples) = raw_data_pos(holes_ind_start(hole_i)+1:...
                holes_ind_start(hole_i)+flag_change_hole_i);
        end

        bsp_full_data_pos = [bsp_full_data_pos fill_in_pos];
        bsp_full_data_ts = [bsp_full_data_ts fill_in_ts];%fill_in_ts_side1 fill_in_ts_side2];


    end


    if hole_i<length(holes_ind_start)
        bsp_full_data_pos = [bsp_full_data_pos raw_data_pos((holes_ind_end(hole_i)):holes_ind_start(hole_i+1))];
        bsp_full_data_ts = [bsp_full_data_ts raw_data_ts_usec((holes_ind_end(hole_i)):holes_ind_start(hole_i+1))];
    end
    hole_i = hole_i + 1;
end

bsp_full_data_pos = [bsp_full_data_pos raw_data_pos((holes_ind_end(end)):end)];
bsp_full_data_ts = [bsp_full_data_ts raw_data_ts_usec((holes_ind_end(end)):end)];
warning('off', 'curvefit:fit:equationBadlyConditioned')

end





