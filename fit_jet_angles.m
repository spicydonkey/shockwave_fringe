%======================================100char======================================================



ycenter=mean([plot2d_ymin,plot2d_ymax]);

fit_tune=0.3;

left_mask_tmin=0.016;
left_mask_tmax=0.025;
y_inner_lim=0.007;
y_outer_lim=0.018;

left_mask_ymin=ycenter-y_outer_lim;
left_mask_ymax=ycenter-y_inner_lim;

right_mask_tmin=left_mask_tmin;
right_mask_tmax=left_mask_tmax;
right_mask_ymin=ycenter+y_inner_lim;
right_mask_ymax=ycenter+y_outer_lim;

mask_right=atoms_centered(:,3)>right_mask_ymin &...
           atoms_centered(:,3)<right_mask_ymax &...
           atoms_centered(:,1)>right_mask_tmin & ...
           atoms_centered(:,1)<right_mask_tmax;
atoms_right=atoms_centered(mask_right,:);
xdata_right=atoms_right(:,3);
figure(13)
hold on
bls_right=regress(atoms_right(:,1),[ones(size(xdata_right,1),1) xdata_right]);
robust_right=robustfit(xdata_right,atoms_right(:,1),'welsch',fit_tune);
%


x=[right_mask_ymin right_mask_ymax].*10^3;
y=robust_right(1)+robust_right(2).*x;
plot(x,y,'g')
rectangle('Position',[right_mask_ymin right_mask_tmin right_mask_ymax-right_mask_ymin right_mask_tmax-right_mask_tmin]*10^3)

pause(0.1)

right_angle=atan(robust_right(2));

mask_left=atoms_centered(:,3)>left_mask_ymin &...
           atoms_centered(:,3)<left_mask_ymax &...
           atoms_centered(:,1)>left_mask_tmin &...
           atoms_centered(:,1)<left_mask_tmax;
atoms_left=atoms_centered(mask_left,:);
xdata_left=atoms_left(:,3);

bls_left=regress(atoms_left(:,1),[ones(size(xdata_left,1),1) xdata_left]);
robust_left=robustfit(xdata_left,atoms_left(:,1),'welsch',fit_tune);


x=[left_mask_ymin left_mask_ymax].*10^3;
y=robust_left(1)+robust_left(2).*x;
plot(x,y,'g')
rectangle('Position',[left_mask_ymin left_mask_tmin left_mask_ymax-left_mask_ymin left_mask_tmax-left_mask_tmin]*10^3)

%-4.02+
hold off

figure(41)
clf
hold on
scatter(atoms_right(:,3),atoms_right(:,1),'.')
scatter(atoms_left(:,3),atoms_left(:,1),'.')
%plot(xdata_right,bls_right(1)+bls_right(2)*xdata_right,'r','LineWidth',2);
plot(xdata_right,robust_right(1)+robust_right(2)*xdata_right,'g','LineWidth',2)
%plot(xdata_left,bls_left(1)+bls_left(2)*xdata_left,'r','LineWidth',2);
plot(xdata_left,robust_left(1)+robust_left(2)*xdata_left,'g','LineWidth',2)
hold off

set(gca,'Ydir','normal')
set(gcf,'Color',[1 1 1]);
title('Spatial Dist. YT')
xlabel('Y(mm)')
ylabel('Z(mm)')
axis square 


left_angle=atan(robust_left(2));

interesect=pi-abs(left_angle)-abs(right_angle);
mach_num=1/cos(interesect/2);
fprintf('Jet angle %f deg giving mach number of %f \n',interesect*(180/pi),mach_num)





%robustdemo(xdata,atoms_right(:,1))