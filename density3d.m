function ndensity=density3d(zxy,edges)
% Calculates 3D density profile given ZXY cell-array
% numdensity = density3d(zxy, edges)

nshot=size(zxy,1);

num_array=histcn(vertcat(zxy{:}),edges{1},edges{2},edges{3});

% [EDGE_Z,EDGE_X,EDGE_Y]=meshgrid(diff(edges{1}),diff(edges{2}),diff(edges{3}));
% % BUG - for variable edge spacing

dedge=cellfun(@(x) x(2)-x(1),edges);    % voxel sides from first diff elem
vvoxel=prod(dedge);

ndensity=num_array./(vvoxel*nshot);     % density in voxel
end