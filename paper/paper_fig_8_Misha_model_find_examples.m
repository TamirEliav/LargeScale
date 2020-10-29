%%
clear
clc

%% load data
load('L:\Misha_attractor\20200902__new_simulations\sim.mat');
load('L:\Misha_attractor\20200902__new_simulations\sim_res.mat');

%%
dlm = 5;
seg = [ 0   80;
        80  160;
        160 240;
        240 320;
        320 400;
        0   200;
        200 400;
        0   400;
        ];
seg = seg + [
    0 0;
    0 -dlm;
    dlm 0;
    0 -dlm;
    dlm -dlm;
    dlm -dlm;
    dlm -dlm;
    dlm 0];
seg = flipud(seg)
levels = [1 2 2 3 3 3 3 3];

%%
run=10;
figure
for ii_cell = 3793:size(m,2)
    [r,c]=find(ind==ii_cell);
    clf
    subplot(211)
    plot(m(:,ii_cell,run));
    title("Cell "+ii_cell);
    subplot(212); hold on
    for ii_net = 1:size(th,2)
        levels = [3 3 3 3 3 2 2 1];
        x = th(:,ii_net);
        x(1:25)=[];
        x(end-25:end)=[];
        y = ones(size(x));
        y = y.*levels(ii_net);
        plot(x,y,'k','linewidth',2,'Color','k');
    end
%     for ii=1:size(seg,1)
%         plot(seg(ii,:),[levels([ii ii])],'k','linewidth',2,'Color','k');
%     end
    for ii_net = 1:length(r)
        net = c(ii_net);
        bn=th(r(ii_net),net);
        plot(bn,levels([net net]),'or');
    end
    hax=gca;
    hax.YLim = [0 3.5];
    hax.Visible = 'off';
%     pause
    dir_out = 'C:\Tamir\work\Projects\LargeScale\Misha_attractor\20200902__new_simulations\Cells';
    file_out = fullfile(dir_out, "Cell "+ii_cell);
    saveas(gcf,file_out,'tif'); 
end


%%
error('stop here');

%% find cells participating only in attractors 6 7 8
is_good_example_connectivity = zeros(1,size(m,2));
for ii_cell = 1:size(m,2)
    cell_num = ii_cell;
    [r,c]=find(ind==cell_num);
    if isequal(c,[6;8]) | isequal(c,[7;8]) | isequal(c,[6;7;8])
        is_good_example_connectivity(ii_cell) = 1;
        fprintf('cell %d, net %d %d\n', cell_num, c);
    end
end
sum(is_good_example_connectivity)
plot(cumsum(is_good_example_connectivity))
%% copy files for corresponding cells
tmpl_list = "Cell "+find(is_good_example_connectivity)+".tif";
tmpl_list = {tmpl_list{:}};
ext = 'tif';
dir_IN = 'L:\Misha_attractor\20200902__new_simulations\Cells\';
dir_OUT = 'L:\Misha_attractor\20200902__new_simulations\Cells\6_7_8';
util_copy_files_by_template(dir_IN, dir_OUT, tmpl_list)

%% find cells participating in as many attractors as possible
num_attractors = zeros(1,size(m,2));
for ii_cell = 1:size(m,2)
    cell_num = ii_cell;
    [r,c]=find(ind==cell_num);
    num_attractors(ii_cell) = length(c);
end
%% plot num fields vs num attractors
num_fields = reshape(res.fnumber,4000,10);
num_fields = num_fields(:,end)';
x = num_attractors;
y = num_fields;
figure
hold on
plot(x + 0.1*randn(size(x)),...
     y + 0.1*randn(size(y)),'.');
refline(1,0);
axis equal
xlabel('number of attractor the cells is participating in');
ylabel('number of fields');
saveas(gcf,'L:\zoom_meetings\20200923\num_fields_vs_num_attractors','tif');
%% plot max field size vs num fields
max_field = nan(size(num_attractors));
for ii_cell = 1:length(max_field)
    if ~isempty(res.fsizes{ii_cell})
        max_field(ii_cell) = max(res.fsizes{ii_cell});
    end
end
x = num_attractors;
% x = num_fields;
y = max_field;
x = x + 0.1*randn(size(x));
y = y + 0.1*randn(size(y));
figure
hold on
plot(x,y,'.');
[n,c] = hist3([x', y'],'Nbins',[7 8]);
contour(c{1},c{2},n',10)
xlabel('number of attractor the cells is participating in');
ylabel('largest field');
saveas(gcf,'L:\zoom_meetings\20200923\largest_field_vs_num_attractors','tif');

%% plot field ratio vs num fields
fields_ratio = reshape(res.sratio,4000,10);
fields_ratio = fields_ratio(:,end)';
% x = num_attractors;
x = num_fields;
y = fields_ratio;
x = x + 0.1*randn(size(x));
y = y + 0.1*randn(size(y));
figure
hold on
plot(x,y,'.');
[n,c] = hist3([x;y]','Nbins',[6 7]);
contour(c{1},c{2},n',logspace(1,3,20));
xlabel('num fields');
ylabel('fields ratio');
saveas(gcf,'L:\zoom_meetings\20200923\ratio_vs_num_fields','tif');


%% copy files for top 100 cells - attractor participation
[num_attractors_sorted, sort_IX]=sort(num_attractors,'descend');
tmpl_list = "Cell "+sort_IX(1:100)+".tif";
tmpl_list = {tmpl_list{:}};
ext = 'tif';
dir_IN = 'L:\Misha_attractor\20200902__new_simulations\Cells\';
dir_OUT = 'L:\Misha_attractor\20200902__new_simulations\Cells\in_many_attractors';
util_copy_files_by_template(dir_IN, dir_OUT, tmpl_list)

%% copy files for top 100 cells - number of fields
[num_attractors_sorted, sort_IX]=sort(num_fields,'descend');
tmpl_list = "Cell "+sort_IX(1:100)+".tif";
tmpl_list = {tmpl_list{:}};
ext = 'tif';
dir_IN = 'L:\Misha_attractor\20200902__new_simulations\Cells\';
dir_OUT = 'L:\Misha_attractor\20200902__new_simulations\Cells\with_many_fields';
util_copy_files_by_template(dir_IN, dir_OUT, tmpl_list)

%% copy files for top 100 cells - ratio
[ratio_sorted, sort_IX]=sort(fields_ratio ,'descend');
tmpl_list = "Cell "+sort_IX(1:100)+".tif";
tmpl_list = {tmpl_list{:}};
ext = 'tif';
dir_IN = 'L:\Misha_attractor\20200902__new_simulations\Cells\';
dir_OUT = 'L:\Misha_attractor\20200902__new_simulations\Cells\with_large_ratio';
util_copy_files_by_template(dir_IN, dir_OUT, tmpl_list)

%% copy files - large ratio + large largest field
[ratio_sorted, sort_IX]=sort(fields_ratio ,'descend');
cells_list = sort_IX(1:200);
cells_list = intersect(cells_list,  find(max_field > 40 & max_field < 80));
tmpl_list = "Cell "+cells_list+".tif";
tmpl_list = {tmpl_list{:}};
ext = 'tif';
dir_IN = 'L:\Misha_attractor\20200902__new_simulations\Cells\';
dir_OUT = 'L:\Misha_attractor\20200902__new_simulations\Cells\with_large_ratio_and_large_field';
util_copy_files_by_template(dir_IN, dir_OUT, tmpl_list)

%% copy files - large ratio + low baseline
baseline = median(m(:,:,end));
[ratio_sorted, sort_IX]=sort(fields_ratio ,'descend');
cells_list = sort_IX(1:1000);
cells_list = intersect(cells_list,  find(baseline==0));
tmpl_list = "Cell "+cells_list+".tif";
tmpl_list = {tmpl_list{:}};
ext = 'tif';
dir_IN = 'L:\Misha_attractor\20200902__new_simulations\Cells\';
dir_OUT = 'L:\Misha_attractor\20200902__new_simulations\Cells\with_large_ratio_and_zero_baseline';
util_copy_files_by_template(dir_IN, dir_OUT, tmpl_list)


%% copy files - manually selected cells
cells_list = [54 121 391 734 1016 1768 1847 1948 1969 2000 2024 2030 3988];
tmpl_list = "Cell "+cells_list+".tif";
tmpl_list = {tmpl_list{:}};
ext = 'tif';
dir_IN = 'L:\Misha_attractor\20200902__new_simulations\Cells\';
dir_OUT = 'L:\Misha_attractor\20200902__new_simulations\Cells\examples_many_field_or_large_ratio';
util_copy_files_by_template(dir_IN, dir_OUT, tmpl_list)

%%








%%
