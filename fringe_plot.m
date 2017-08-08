%plot the fringe

figure_number=15;

rot_angle=0.55;

plot2d_tmin = 0.005;
plot2d_tmax = 0.015;
plot2d_ymin = -0.015;
plot2d_ymax = 0.015;


%main fringe
% rot_angle=0.55;
% xoffset=0*10^-3;
% line_tmin = 0.005;
% line_tmax = 0.016;
% line_ymin = -0.0005;
% line_ymax = 0.0005;

%second fringe
rot_angle=0.73;
xoffset=-3.5*10^-3;
line_tmin = 0.005;
line_tmax = 0.016;
line_ymin = -0.0005;
line_ymax = 0.0005;

sin_theta = sin(rot_angle);
cos_theta = cos(rot_angle);
atoms_centered_rot=[];
atoms_centered_rot(:,1) = atoms_centered(:,1)*cos_theta- atoms_centered(:,3)*sin_theta;
atoms_centered_rot(:,3) =xoffset+ atoms_centered(:,1)*sin_theta+ atoms_centered(:,3)*cos_theta;

tspread=files.velocity*(windows.all.tmax-windows.all.tmin)/2;

%main fringe
rot_angle=0.55;




plot2d_xmin = windows.all.xmin;
plot2d_xmax = windows.all.xmax;




line_mask=line_ymax>atoms_centered_rot(:,3)& atoms_centered_rot(:,3)>line_ymin & ...
    line_tmax>atoms_centered_rot(:,1) & atoms_centered_rot(:,1)>line_tmin;

figure(18)
%scatter(atoms_centered_rot(line_mask,3),atoms_centered_rot(line_mask,1)) 
[n,edges]=hist(atoms_centered_rot(line_mask,1),300);
n=n/(plot2d_tmax-plot2d_ymin);
plot(edges*10^3,n)
set(gcf,'Color',[1 1 1]);
title('Second Fringe')
xlabel('x(mm)')
ylabel('Linear Count Density m^{-1}')


Logscale_bool=0;
spatial_blur=0.5;
binsx=30;
binsy=350;%3508;
binsz=496;%4961;
if spatial_blur<0
    spatial_blur=0;
end
%if spatial_blur>bins/2
%    spatial_blur=0;
%end

if mod(binsx,2) == 0 %if the number of bins are even make them odd
    binsx=binsx+1;
end
if mod(binsy,2) == 0 %if the number of bins are even make them odd
    binsy=binsy+1;
end
if mod(binsz,2) == 0 %if the number of bins are even make them odd
    binsz=binsz+1;
end

if ~length(atoms_centered)==0
    XEdges=linspace(plot2d_xmin,plot2d_xmax,binsx);
    YEdges=linspace(plot2d_ymin,plot2d_ymax,binsy);
    TEdges=linspace(plot2d_tmin,plot2d_tmax,binsz);

    panes=sum([0,0,1]);
    pane_counter=1;

    figure(figure_number);
    clf
    set(gcf,'Units','normal')
    set(gca,'Position',[0 0 1 1])
    
    
    subplot(1,panes,pane_counter);
    pane_counter=pane_counter+1;
    bin_area=((plot2d_ymax-plot2d_ymin)/binsy)*((plot2d_tmax-plot2d_tmin)/binsz);
    [counts,centers]=hist3(atoms_centered_rot(:,[1,3]),'edges',{TEdges,YEdges});
    %counts=counts/bin_area;
    if  ~spatial_blur==0
        filter=fspecial('gaussian',round(10*spatial_blur),spatial_blur);
        counts=imfilter(counts, filter, 'replicate');
    end
    if Logscale_bool
        counts=log10(counts);
    end
    %imagesc seems to plot the wrong way round so we transpose here
    imagesc(10^3*centers{2},10^3*centers{1},counts)
    hold on
    xlim_val=xlim;
    center_pt=mean(xlim_val);
    x=[center_pt center_pt];
    y=[-50 50];
    line(x,y)
    ylim_val=ylim;
    center_pt=mean(ylim_val);
    y=[center_pt center_pt];
    x=[-50 50];
    line(x,y)
    hold off



    set(gca,'Ydir','normal')
    set(gcf,'Color',[1 1 1]);
    title('Spatial Dist. YT')
    xlabel('Y(mm)')
    ylabel('Z(mm)')
    axis square 
    h=colorbar;
    if Logscale_bool
       xlabel(h,'Log Count Density (m^{-1}s^{-1})')
    else
       xlabel(h,'Count Density (m^{-1}s^{-1})')
    end



%             unit16max=2^16-1;
%            scaled_counts=log10(double(counts)+10)-1;
%            scaled_counts=double(counts).^0.4;
%             
%             norm_counts=uint16(scaled_counts/max(max(scaled_counts))*unit16max);
%             rgbcounts=repmat(norm_counts,1,1,3);
%              rgbcounts= ind2rgb(norm_counts, parula(unit16max));
%              imwrite(flipud(rgbcounts),'outparula.tif')
%              
%              rgbcounts= ind2rgb(norm_counts, gray(unit16max));
%              imwrite(flipud(rgbcounts),'outgrey.tif')
%              
%              rgbcounts= ind2rgb(norm_counts, flipud(gray(unit16max)));
%              imwrite(flipud(rgbcounts),'outstar.tif')
%             
%             rgbcounts= ind2rgb(norm_counts, hot(unit16max));
%             imwrite(flipud(rgbcounts),'outhot.tif')
%             
%             rgbcounts= ind2rgb(norm_counts, copper(unit16max));
%             imwrite(flipud(rgbcounts),'outcopper.tif')
%             
%             rgbcounts= ind2rgb(norm_counts, flipud(hot(unit16max)));
%             imwrite(flipud(rgbcounts),'outhotinverse.tif')
%             
%             rgbcounts= ind2rgb(norm_counts, flipud(gray(unit16max)));
%             imwrite(flipud(rgbcounts),'outstar.tif')
% 
%             rgbcounts= ind2rgb(norm_counts, winter(unit16max));
%             imwrite(flipud(rgbcounts),'outwinter.tif')
%             
%             rgbcounts= ind2rgb(norm_counts, bone(unit16max));
%             imwrite(flipud(rgbcounts),'outbone.tif')

end
    


