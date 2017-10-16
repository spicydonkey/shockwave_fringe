%% configs
%% main_shockfringe
configs.flags.verbose=2;
configs.flags.savedata=0;
configs.flags.archive_txy=1;
configs.flags.force_all_stages=0;
configs.flags.graphics=1;
configs.flags.build_txy=1;

configs.misc.hbar=1.055e-34;    % reduced Planck constant [Js]
configs.misc.m=6.647e-27;       % mass of helium [kg]
configs.misc.tof=0.416;         % free-fall TOF from trap to detector [s]
configs.misc.g=9.807;           % acceleration due to gravity [ms^-2]
configs.misc.vel_z=configs.misc.g*configs.misc.tof; % free-fall velocity at detection [m/s]
configs.misc.det_qe=0.1;        % detector quantum efficiency

configs.files.path='\\AMPLPC29\He BEC Archive\EXPERIMENT-DATA\AL_Shockwaves\medium_number_al\20170717_atomlaser\d';

configs.files.dir_data=fileparts(configs.files.path);    % fullpath to data directory
configs.files.archive=fullfile(configs.files.dir_data,'archive');   % dir to archive folder
configs.files.dirout=fullfile(configs.files.dir_data,'output');      % output directory (will be time-stamped)

configs.load.version=1.1;         % TXY load stage version number

configs.load.id=1:1862;         % file id numbers to use for analysis
configs.load.mincount=1000;         % min counts in window - 0 for no min
configs.load.maxcount=Inf;          % max counts in window - Inf for no max

configs.load.rot_angle=0.61;    % detector/trap alignment

% construct window for ROI
configs.load.window{1}=[0.45,0.7];      % T [s] include all PALs
configs.load.window{2}=[-5e-3,10e-3];    % X [m]
configs.load.window{3}=[-20e-3,30e-3];    % Y [m]

configs.image.voxel_res=1e-4*[1,1,1];   % ZXY voxel resolution [m]
configs.image.size=[-30e-3,40e-3;-5e-3,5e-3; -25e-3,25e-3];   % image size/lims [T;X;Y] [m]

% PAL
configs.pal.t1=0.4658;      % estimated from the first pal flux peak on DLD front panel
configs.pal.dt=0.0212;      % estimated from peak diff of first and second
configs.pal.n=10;

configs.pal.nfitstart=2;    % AL pulse ID to begin fit from

%% main_process_pal
% configure the 1d density profile
fringe_cfg.offset=[0,-0.006];       % yz translation for centering transformation
fringe_cfg.theta=1.1;               % rotation angle around x-axis (rad) - to align fringes in vertical Y-axis
fringe_cfg.width=2e-3;              % 1d profile perpendicularly [m]
fringe_cfg.dlim=[0,20e-3];          % line profile limit

nsmooth_1d_raw=5;   % moving average sample size - smoothing raw 1D profile

%% DemoFringePeak
% background smoothing
n_sm_bgd=15;
n_sm_post=10;

% peak detection algorithm
slopethreshold=1e-5;
ampthreshold=0.05;
smoothwidth=0;          % 0 to no smoothing in peak finding
peakgroup=5;
smoothtype=1;   % smoothing is off

%% bootstrapping
% rng
rng('shuffle');

% data subset size
bootstrap_ndata=0.2;    % ratio subset size to all
% number of sampling
bootstrap_Nsamp=50;