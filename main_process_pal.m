%% Process PAL to analyse 1D density fringes through the jet
% transforms 3D translation (YZ) + rotation operation (X) to align fringes
% along vertical (Z') direction. take 2D density profile in Y'Z'
% projection (integrates X) - take finite width line profile along Z' to
% obtain 1d density fringes.

%% configure
% configure the 1d density profile
fringe_cfg.offset=[0,-0.006];       % yz translation for centering transformation
fringe_cfg.theta=1.1;               % rotation angle around x-axis (rad) - to align fringes in vertical Y-axis
fringe_cfg.width=2e-3;              % 1d profile perpendicularly [m]
fringe_cfg.dlim=[0,20e-3];          % line profile limits

% density - histogramming (in rotated axis)
yy_ed=edges{3};
zz_ed=edges{1};
yy_c=cents{3};
zz_c=cents{1};

nsmooth_1d_raw=5;   % moving average sample size - smoothing raw 1D profile

% % plotting
% plot_ncol=ceil(sqrt(pal_nseq));         % subplot num cols
% plot_nrow=ceil(pal_nseq/plot_ncol);     % subplot num rows

% cc=distinguishable_colors(pal_nseq);	% colors for plotting against PAL properties

%% main
% get indices for data in fringe ROI
id_yy=find(abs(yy_c)<fringe_cfg.width);       % perpendicular indices to integrate over
id_zz=find(zz_c>=fringe_cfg.dlim(1)&zz_c<fringe_cfg.dlim(2));    % y-indices ROI

% get pos vectors for ROI
yy=yy_c(id_yy);
zz=zz_c(id_zz);

% build ROI position for 'rectangle' graphical annotation
roi_x=min(yy);
roi_y=min(zz);
roi_w=diff(minmax(yy));
roi_h=diff(minmax(zz));
pos_roi=[roi_x,roi_y,roi_w,roi_h];

% preallocate
nn1d=cell(pal_nseq,1);      % 1D density profile

if vgraph>0
    hfig_ndenrot=figure();
    hfig_ndenraw=figure();
    hfig_nden1d=figure();
end
for pal_id=1:pal_nseq
    % variable 'pal' is Nx3 array
    pal=vertcat(pal_zxy0{pal_id}{:});    % collate all shots in this PAL
    
    %%% transform PAL zxy to align jet to Z-axis
    % cull X --> unused
    pal=pal(:,[3,1]);   % Y'Z'
    
    % translate in ZY plane - centre to jet crossing
    pal=pal+repmat(fringe_cfg.offset,[size(pal,1),1]);
    
    % rotate about X-axis
    pal_rot=pal;   % preallocate
    pal_rot(:,1)=cos(fringe_cfg.theta)*pal(:,1)-sin(fringe_cfg.theta)*pal(:,2);
    pal_rot(:,2)=sin(fringe_cfg.theta)*pal(:,1)+cos(fringe_cfg.theta)*pal(:,2);
    
    %%% evaluate 2D density
    nn2d=density2d(pal_rot,{yy_ed,zz_ed})';
    
    
    %%% prepare 1D density profile through shockwave
    nn_raw=nn2d(id_zz,id_yy);                 % raw 2d density in region of interest
    nn1d_temp=mean(nn_raw,2);               % integrate thru perpendicular dir
    nn1d{pal_id}=smooth(nn1d_temp,nsmooth_1d_raw);       % simple smoothing - moving average
    
    %%% plot result
    if vgraph>0
        % plot full 2D projected density profile
        figure(hfig_ndenrot);
        subplot(plot_nrow,plot_ncol,pal_id);
        imagesc(yy_c,zz_c,nn2d);
        set(gca,'YDir','normal');
        axis equal; axis tight;
        title(sprintf('%d',pal_id));
        
        % 2D density in 2D region of interest
        figure(hfig_ndenraw);
        subplot(plot_nrow,plot_ncol,pal_id);
        imagesc(yy,zz,nn_raw);
        set(gca,'YDir','normal');
        title(sprintf('%d',pal_id));
        
        % draw rectangle for ROI on 2D density plot
        figure(hfig_ndenrot);
        subplot(plot_nrow,plot_ncol,pal_id);
        hold on;
        rectangle('Position',pos_roi,'EdgeColor','w');
        
        % 1D density profile
        figure(hfig_nden1d);
        hold on;
        plot(1e3 *zz,nn1d{pal_id},...
            'DisplayName',sprintf('%d: %0.2g, %0.2g',pal_id,Nal(pal_id),N0(pal_id)),...
            'LineWidth',1.5,'color',cc(pal_id,:));
    end
end
%%% annotate figures
if vgraph>0
    % raw density in region of interest
    figure(hfig_ndenraw);
    
    % 1D density profile
    figure(hfig_nden1d);
    lgd=legend('show');
    title(lgd,'PAL: $N_{AL}$, $N_0$');
    box on;
    xlabel('distance [mm]');
    ylabel('density [arb]');
end