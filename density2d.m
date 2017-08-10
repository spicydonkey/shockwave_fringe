function ndensity=density2d(xy,edges)
% Calculates 2D density profile given 2D numcountsx2 array
% numdensity = density2d(xy, edges)

% nshot=size(xy,1);

% num_array=histcn(vertcat(xy{:}),edges{1},edges{2},edges{3});
num_array=histcn(xy,edges{1},edges{2});

dedge=cellfun(@(x) x(2)-x(1),edges);    % voxel sides from first diff elem
apixel=prod(dedge);     % pixel area

% ndensity=num_array./(apixel*nshot);     % density in voxel
ndensity=num_array./apixel;     % density
end