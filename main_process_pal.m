%% Process PAL to analyse 1D density fringes through the jet
% configure
pal_id=1;
% jet_yz_shift=[0,0];
% jet_theta=0;
% 
jet_yz_shift=[0,-0.004];
jet_theta=0.7;

% variable 'pal' is Nx3 array
pal=vertcat(pal_zxy0{pal_id}{:});    % collate all shots in this PAL

%%% transform PAL zxy to align jet to Z-axis
% cull X --> unused
pal=pal(:,[3,1]);   % YZ

% translate in TY plane - centre to jet crossing
pal=pal+repmat(jet_yz_shift,[size(pal,1),1]);

% rotate about X-axis
pal_xrot=pal;   % preallocate
pal_xrot(:,1)=cos(jet_theta)*pal(:,1)-sin(jet_theta)*pal(:,2);
pal_xrot(:,2)=sin(jet_theta)*pal(:,1)+cos(jet_theta)*pal(:,2);

%%% evaluate 2D density
yrot_ed=edges{3};
zrot_ed=edges{1};
yrot_c=cents{3};
zrot_c=cents{1};

nn=density2d(pal_xrot,{yrot_ed,zrot_ed})';      % needs to be transposed

% visualise
hfig_ndenrot=figure();
imagesc(yrot_c,zrot_c,nn);
set(gca,'YDir','normal');
axis equal;
xlabel('Y');
ylabel('Z');

% [X,Y]=meshgrid(yrot_c,zrot_c);
% s=surf(X,Y,nn,'EdgeColor','none');
% view(2);
% shading interp;

%%% take line profile of the density thru Y-axis at X=0
[~,id_x0]=min(abs(yrot_c));

id_y=find(zrot_c>=0&zrot_c<0.02);

% nn1d=nn(id_y,id_x0);   % 1d profile

% plot
hfig_nden1d=figure();
% plot(zrot_c(id_y),nn(id_y,id_x0));
hold on;
for ii=0
    plot(zrot_c(id_y),nn(id_y,id_x0+ii));           % raw 1D
    plot(zrot_c(id_y),smooth(nn(id_y,id_x0+ii)));   % smoothed
end
box on;
xlabel('distance [m]');
ylabel('density [arb]');