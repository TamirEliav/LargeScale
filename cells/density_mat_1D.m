function [density_mat, data_bin] = density_mat_1D(data,bin_centers)
% this function returns the density of the input data across a space
% defined by bin_centers. it basically counts the number of data points
% that fall within each of the bins specified in bin_centers. 

% the case for which the size of data is not a vector is usually for
% shuffled spikes case. for this case thr caluclation is made separately 
% for each row of data, hence  density_mat is a matrix. 
if isempty(data)
    data_bin = nan(size(data));
    density_mat = nan(size(bin_centers));
    return;
end
size_data = size(data);
if min(size_data)==1
    
    bin_size = bin_centers(2) - bin_centers(1);
    data(isnan(data)) = [];
    data(data<bin_centers(1)-bin_size/2) = [];data(data>bin_centers(end)+bin_size/2) = [];
    
    dd = interp1(bin_centers, 1:length(bin_centers), data,'nearest');
    % dd holds the bin number for every data point
    tmp = find(isnan(dd));
    if ~isempty(tmp)
        for i = 1:length(tmp)
            if data(tmp(i))>=bin_centers(end)
                dd(tmp(i)) = length(bin_centers);
            elseif data(tmp(i))<=bin_centers(1)
                dd(tmp(i)) = 1;
            end
        end
    end
    
    data_bin = bin_centers(dd);
    % data_bin holds the bin center for every data point, e.g. the spike's bin location
    
    density_mat = accumarray(dd(:),1,[length(bin_centers),1]);
    
else
    bin_size = bin_centers(2) - bin_centers(1);
    data(data<bin_centers(1)-bin_size/2) = [];data(data>bin_centers(end)+bin_size/2) = [];
    
    dd = interp1(bin_centers, 1:length(bin_centers), data(:),'nearest');
    % dd holds the bin number for every data point
    tmp = find(isnan(dd));
    if ~isempty(tmp)
        for i = 1:length(tmp)
            if data(tmp(i))>=bin_centers(end)
                dd(tmp(i)) = length(bin_centers);
            elseif data(tmp(i))<=bin_centers(1)
                dd(tmp(i)) = 1;
            end
        end
    end
    
    data_bin = bin_centers(dd); data_bin = reshape(data_bin,size_data);
    % data_bin holds the bin center for every data point, e.g. the spike's bin location
    
    dd = reshape(dd,size_data);
    for i=1:size_data(1)
        density_mat(i,:) = accumarray(dd(i,:)',1,[length(bin_centers),1])';
    end
end

end


