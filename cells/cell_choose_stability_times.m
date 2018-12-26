function cell_choose_stability_times(cell_ID)

%% load cell data
cell=cell_load_data(cell_ID,'details');
% exp=exp_load_data(cell.details.exp_ID, 'details');
ti = exp_get_sessions_ti(cell.details.exp_ID, 'Behave');

%% let the user choose the ts from the stability plot
cell_plot_map_fields(cell_ID); % consider changing this to simply plot the stability panel (so it will be large), but then we don't have the other panels...
x=[];
while length(x) < 2
    disp('choose at least two points (the last two points will be taken,sorted')
    [x,~] = ginput();
end
x(1:end-2,:)=[];
sort(x,'ascend');

%% add the ts the user chose to the plot
plot(repelem(x(1),2),get(gca,'ylim'),'-','LineWidth',2,'Color','b')
plot(repelem(x(2),2),get(gca,'ylim'),'-','LineWidth',2,'Color','b')

%% display the ts values so the user can copy it to excel (in nlg ts conversion)
stability_ts = x .* 60*1e6; % convert to usec
stability_ts = x + ti(1);   % add behave session offset
stability_ts = round(stability_ts);

fprintf('stability ts: [%d %d]\n',stability_ts);


end

















