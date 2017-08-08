%======================================100char======================================================
%This program will fit the jet angles for the pulse sequence in Run 2 background no raman\d
%%300dpi A3 = 3508*4961 pixels
%-----------------------------------START user var----------------------------------------
 set(0, 'DefaultFigureRenderer', 'zbuffer');
 
use_txy=1;                  %if false will remake the txy_forc files
files.count_min=1000;       %min mumber of counts to read in 
files.rot_angle=0.64;       %rotation anle from the txy data 
use_saved_data=0; 
files.do_pos_correction=0;	%find the halo pos from the ceneter of bragg orders (position correction)


files.velocity=9.8*0.413;

files.path='V:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20170716_atomlaser\d';
files.save_all_points=0;    %create array that off all the points usefull for intial investigation and plot TOF and 3d
files.numstart=1;           %start file num
files.numtoimport=1000;%99;       %number of files to import

n=0;
pulse_time=0.021;

windows.all.tmin=0.462+pulse_time*n;
windows.all.tmax=0.476+pulse_time*n;
windows.all.xmin=-0.03;
windows.all.xmax=0.030;
windows.all.ymin=-0.03;
windows.all.ymax=0.030;



% files.path='C:\Tmp\Normal Trap Raman Sweep\d';
% files.save_all_points=0;    %create array that off all the points usefull for intial investigation and plot TOF and 3d
% files.numstart=1;           %start file num
% files.numtoimport=18;%99;       %number of files to import
% 
% windows.all.tmin=0.41;
% windows.all.tmax=0.8;%0.43;
% windows.all.xmin=-0.03;
% windows.all.xmax=0.030;
% windows.all.ymin=-0.03;
% windows.all.ymax=0.030;
%     tspread=(windows.all.tmax-windows.all.tmin)/2;
%     plot2d_tmin = -tspread;%-0.0227;
%     plot2d_tmax = tspread;%0.045;
%     plot2d_ymin = -0.0225;
%     plot2d_ymax = 0.0225;
%     plot2d_xmin = windows.all.xmin;
%     plot2d_xmax = windows.all.xmax;
      

% files.path='C:\Tmp\20150626-pulsed_atom_laser-biased_Btrap\data\dld_data1\pulsed_atom_laser_biased_Btrap_';
% files.save_all_points=0;    %create array that off all the points usefull for intial investigation and plot TOF and 3d
% files.numstart=1;           %start file num
% files.numtoimport=1996;       %number of files to import
% 
% windows.all.tmin=1.26;
% windows.all.tmax=1.28;
% windows.all.xmin=-0.03;
% windows.all.xmax=0.030;
% windows.all.ymin=-0.03;
% windows.all.ymax=0.030;



progress_scaling=0.1;

plot2d_hist=1;              %plot a 2d histogram(image) of the halo colapsed in time
    plot2d_hist_binw=0.0001; %bin width for hist in meters
    plot2d_hist_gauss=0.000;%apply gaussian blur in meters, for off set to zero otherwise specify blur radius



%-----------------------------------END user var----------------------------------------


[atoms_centered_cells,atoms_centered,com,all_points]=ImportShockData(files,windows,use_txy,use_saved_data,1,progress_scaling);


if ~files.do_pos_correction
    mean_atoms=mean(atoms_centered);
    atoms_centered=atoms_centered-repmat(mean_atoms,size(atoms_centered,1),1);
end


% figure(12)
% set(gcf,'Color',[1 1 1]);
% histogram(atoms_centered(:,1),1000)




















