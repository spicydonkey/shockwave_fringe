%% Process PAL to analyse 1D density fringes through the jet
%%% configure
jet_yz_shift=[0,-0.005];
jet_theta=0.7;

profile_width=0.5e-3;   % 1d profile perpendicularly [m]

% preallocate
nn1d=cell(pal_nseq,1);      % 1D density profile

hfig_ndenrot=figure();
hfig_ndenraw=figure();  
plot_ncol=ceil(sqrt(pal_nseq));         % configure subplots layout
plot_nrow=ceil(pal_nseq/plot_ncol);
hfig_nden1d=figure();
cc=distinguishable_colors(pal_nseq);        % colors
for pal_id=1:pal_nseq
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
    
    % plot 2D projected density profile
    figure(hfig_ndenrot);
    subplot(plot_nrow,plot_ncol,pal_id);
    imagesc(yrot_c,zrot_c,nn);
    set(gca,'YDir','normal');
    axis equal; axis tight; 
    title(sprintf('%d',pal_id));
    
    %%% take line profile of the density thru Y-axis at X=0
    [~,id_x0]=min(abs(yrot_c));
    
    id_x=find(abs(yrot_c)<profile_width);       % perpendicular indices to integrate over
    id_y=find(zrot_c>=0&zrot_c<0.02);           % y-indices ROI
    
    %%% prepare 1D density profile through shockwave
    nn_raw=nn(id_y,id_x);            % raw density in region of interest
    nn1d_temp=mean(nn_raw,2);        % integrate thru perpendicular dir
    nn1d{pal_id}=smooth(nn1d_temp,5);      % simple smoothing - moving average
    
    %%% plot result
    % raw density in region of interest
    figure(hfig_ndenraw);
    subplot(plot_nrow,plot_ncol,pal_id);
    imagesc(yrot_c(id_x),zrot_c(id_y),nn_raw);
    set(gca,'YDir','normal');
    title(sprintf('%d',pal_id));
    
    % draw rectangle for ROI on 2D density plot
    figure(hfig_ndenrot);
    subplot(plot_nrow,plot_ncol,pal_id);
    hold on;
    pos_roi=[min(yrot_c(id_x)),min(zrot_c(id_y)),max(yrot_c(id_x)),max(zrot_c(id_y))];     % build rectangle: 2 corners
    rectangle('Position',pos_roi,'EdgeColor','w');
    
    % 1D density profile
    figure(hfig_nden1d);
    hold on;
    plot(zrot_c(id_y),nn1d{pal_id},...
        'DisplayName',sprintf('%d: %0.2g, %0.2g',pal_id,fitval.y(pal_id),bec_n(pal_id)),...
        'LineWidth',1.5,'color',cc(pal_id,:));
end
%%% annotate figures
% raw density in region of interest
figure(hfig_ndenraw);

% 1D density profile
figure(hfig_nden1d);
lgd=legend('show');
title(lgd,'PAL: $n_{AL}$, $N_0$');
box on;
xlabel('distance [m]');
ylabel('density [arb]');
