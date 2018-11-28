function pos_proj = POS_calc_linearized(pos,calib)

%% project to the midline calibration
[xy,distance,t] = distance2curve(...
    calib.curvexy(1:50:end,:),...        % todo: move the sabsampling to the creation of the calibration...
    [pos(:,1) pos(:,2)],...
    'linear'); 

%% calc side of the tunnel (curve)
IX = knnsearch(calib.curvexy,xy,'K',2,'IncludeTies',false,'Distance','euclidean');
IX = sort(IX,2,'ascend'); % ascend will create the points north to the tunnel positive!
side = zeros(size(xy,1),1);
for ii_point = 1:size(xy,1)
    B1 = calib.curvexy(IX(ii_point,1),:);
    B2 = calib.curvexy(IX(ii_point,2),:);
    side(ii_point) = POS_side_of_line(B1,B2,pos(ii_point,:));
end

%% 
pos_proj = [t.*calib.tunnel_length...
            side.*distance];

end





