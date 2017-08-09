clear all;

%% configs
configs.flags.verbose=1;
configs.flags.savedata=1;
configs.flags.archive_txy=1;
configs.flags.force_all_stages=0;
configs.flags.graphics=1;
configs.flags.build_txy=1;

configs.misc.vel_z=9.8*0.416;    % atom free-fall vert v at detector hit for T-to-Z conversion;
configs.misc.det_qe=0.1;

configs.files.path='C:\Users\HE BEC\Documents\lab\shockwave\20170716_atomlaser\d';

configs.files.dir_data=fileparts(configs.files.path);    % fullpath to data directory
configs.files.archive=fullfile(configs.files.dir_data,'archive');   % dir to archive folder
configs.files.dirout=fullfile(configs.files.dir_data,'output');      % output directory (will be time-stamped)

configs.load.version=1.1;         % TXY load stage version number

configs.load.id=1:3615;         % file id numbers to use for analysis
configs.load.mincount=1000;         % min counts in window - 0 for no min
configs.load.maxcount=Inf;          % max counts in window - Inf for no max

configs.load.rot_angle=0.61;    % detector/trap alignment

% construct window for ROI
configs.load.window{1}=[0.45,0.7];      % T [s] include all PALs
configs.load.window{2}=[-5e-3,10e-3];    % X [m]
configs.load.window{3}=[-20e-3,30e-3];    % Y [m]

configs.image.voxel_res=2e-4*[1,1,1];   % ZXY voxel resolution [m]
configs.image.size=[-30e-3,40e-3;-5e-3,5e-3; -25e-3,25e-3];   % image size/lims [T;X;Y] [m]

% PAL
configs.pal.t1=0.4658;      % estimated from the first pal flux peak on DLD front panel
configs.pal.dt=0.0212;      % estimated from peak diff of first and second
configs.pal.n=10;


%% initialise
do_next=configs.flags.force_all_stages;
verbose=configs.flags.verbose;
t_main_start=tic;   % for reporting process duration
datetimestr=datestr(datetime,'yyyymmdd_HHMMSS');    % timestamp when function called
configs.files.dirout=[configs.files.dirout,'_',datetimestr];
vz=configs.misc.vel_z;
det_qe=configs.misc.det_qe;
HFIG={};        % cell array for all figures generated

if configs.flags.savedata
    % output directory
    if ~isdir(configs.files.dirout)
        warning(['output directory "',configs.files.dirout,'" does not exist. Creating directory...']);
        mkdir(configs.files.dirout);
    end
end

% check for archive directory
if ~isdir(configs.files.archive)
    warning(['archive directory "',configs.files.archive,'"does not exist. Creating directory...']);
    mkdir(configs.files.archive);
end

%% load txy data
loadconfigs=configs.load;   % abbreviated configs storing only txy load part

% check for existing saved file for preloaded data
arch_mat_list=dir([configs.files.archive,'/txy_*.mat']);      % list of saved data

if ~do_next
    % check for existing saved file for preloaded data
    arch_mat_list=dir([configs.files.archive,'/txy_*.mat']);      % list of saved data
    
    do_next=1;      % to do all task until saved data found
    if size(arch_mat_list,1)~=0
        % look for file with same load configs
        for ii=1:size(arch_mat_list,1)
            this_file=arch_mat_list(ii).name;
            S_temp=load([configs.files.archive,'/',this_file],'loadconfigs');   % load configs from prev data
            
            if isequal(S_temp.loadconfigs,loadconfigs)
                % success! this is the data to be loaded
                warning('Loading TXY: Match found in archive. Loading %s.',[configs.files.archive,'/',this_file]);
                load([configs.files.archive,'/',this_file],'txy','fout');
                do_next=0;      % skip full stage
                break
            end
        end
        
        % no saved txy found
        if do_next
            warning('Loading TXY: No relevant TXY data found in archive. Setting do_next=1.');
        end
        
    end
end
clear S_temp this_file arch_mat_list;       % clean workspace

if do_next
    [txy,fout,HFIG{length(HFIG)+1}]=load_txy(configs.files.path,configs.load.id,...
        configs.load.window,configs.load.mincount,configs.load.maxcount,...
        configs.load.rot_angle,configs.flags.build_txy,verbose,configs.flags.graphics);
    
    % save loaded TXY data to archive for fast loading
    if configs.flags.archive_txy&&(~configs.flags.force_all_stages)     % don't save when debugging - creates multiplicity
        % save the loaded txy, output log, and configs used
        fpath_txy_archive=[configs.files.archive,'/txy_',datetimestr,'.mat'];
        warning('Archiving TXY data as .mat file: %s',fpath_txy_archive);
        save(fpath_txy_archive,'txy','fout','loadconfigs');
    end
end

clear loadconfigs fpath_txy_archive;  % clean workspace

%% pre-process
nshot=size(txy,1);

% txy-->zxy
zxy=cellfun(@(x) double(x.*repmat([vz,1,1],[size(x,1),1])),txy,'UniformOutput',false);

% check loaded data
if verbose>0
    hfig_all=figure();
    plot_zxy(zxy,3e4,1,'k');
    axis equal;
    view(90,0);
end

% get PAL
% build config
pal_z1=configs.pal.t1*vz;
pal_dz=configs.pal.dt*vz;
pal_nseq=configs.pal.n;

pal_zxy=capture_pal(zxy,pal_z1,pal_dz,pal_nseq);

%%% centre PAL to a common mean position
%   NOTE: centering each PAL in shot to mean position reduces fringe
%   visibility, most likely from asymmetry of PAL and shot-to-shot N
%   fluctuation

% get centre
pal_cent_avg=zeros(pal_nseq,3);
pal_cent_std=zeros(pal_nseq,3);
for ii=1:pal_nseq
    % use mean of each shot's centre as an approximation to the shot collated
    % mean position
    this_pal_cent_array=cellfun(@(x) mean(x,1),pal_zxy{ii},'UniformOutput',false);
    this_pal_cent_array=vertcat(this_pal_cent_array{:});    % form into nshotx3 array
    
    pal_cent_avg(ii,:)=mean(this_pal_cent_array,1);     % average zxy center in this pulse
    pal_cent_std(ii,:)=std(this_pal_cent_array,1);      % uncertainty in the mean centre position
end

% centre PAL
pal_zxy0=cell(pal_nseq,1);     % preallocate centered PAL array
for ii=1:pal_nseq
    % shift to evaluated centre
    pal_zxy0{ii}=cellfun(@(x) x-repmat(pal_cent_avg(ii,:),[size(x,1),1]),pal_zxy{ii},'UniformOutput',false);
end

% summarise captured PALs
if verbose>0
    hfig_pal=figure();
    plot_ncol=ceil(sqrt(pal_nseq));
    plot_nrow=ceil(pal_nseq/plot_ncol);
    for ii=1:pal_nseq
        subplot(plot_nrow,plot_ncol,ii);
        plot_zxy(pal_zxy0{ii},1e4,1,'k');
        view(90,0);
        axis equal;
        axis tight;
        box on;
        ht=sprintf('PAL: %d',ii);
        title(ht);
    end
end


%% Process PAL
%%% Evaluate atom numbers
% number in PAL
pal_n=zeros(pal_nseq,2);    % preallocate; format: [avg, std]
for ii=1:pal_nseq
    this_pal_n=cellfun(@(x) size(x,1),pal_zxy0{ii});    % number in nth PAL 'detected'
    pal_n(ii,:)=(1/det_qe)*[mean(this_pal_n),std(this_pal_n)];      % estimate actual number in PAL
end

%%%% fit to estimate number in condensate (N0) and outcoupling efficiency
% (eff_oc)
n_i=2:pal_nseq;     % TODO pulse numbers to use - don't use sat'd
param0=[1e5,0.1];   % estimate for initial param

modelfun='y~N0*(1-r)^(x1-1)*r';     % x1: PAL index; N0: BEC number; r: outcoupling frac; y: number in x1^th PAL; 
coeffnames={'N0','r'};

% define optimiser
fo = statset('TolFun',10^-6,...
    'TolX',10^-6,...
    'MaxIter',10^6,...
    'UseParallel',0);

% do the fit
fitobject=fitnlm(n_i,pal_n(n_i,1),modelfun,param0,...
    'CoefficientNames',coeffnames,'Options',fo);

% get fit results
paramfit=[fitobject.Coefficients.Estimate,fitobject.Coefficients.SE];
fitval.x=1:pal_nseq;
fitval.y=feval(fitobject,fitval.x);

% get number in BEC
bec_n=(paramfit(1,1)*(1-paramfit(2,1)).^(1:pal_nseq))';

% summarise
linewidth=1.5;
namearray={'LineWidth','MarkerFaceColor','Color'};      % error bar graphics properties
valarray={linewidth,'w','k'};                 % 90 deg (normal) data
if verbose>0
    hfig_atom_number=figure();
    % plot data
    hdata_pal_n=ploterr(1:pal_nseq,pal_n(:,1),[],pal_n(:,2),'o','hhxy',0);
    set(hdata_pal_n(1),namearray,valarray,'DisplayName','Data');
    set(hdata_pal_n(2),namearray,valarray,'DisplayName','');
    
    hold on;
    % plot fit
    hfit=plot(fitval.x,fitval.y,'r*','DisplayName','Fit');
    
    % annotate
    legend([hdata_pal_n(1),hfit]);
    xlabel('Pulse number');
    ylabel('Number in PAL');
end

%% PAL density image
%%% configure voxels
% construct edge/center vectors for each dim
edges=cell(3,1);
cents=cell(3,1);
for ii=1:3
    edges{ii}=configs.image.size(ii,1):configs.image.voxel_res(ii):configs.image.size(ii,2);
    cents{ii}=0.5*(edges{ii}(1:end-1)+edges{ii}(2:end));
end

%%% 3D density profile (full)
nden3=cellfun(@(x) density3d(x,edges),pal_zxy0,'UniformOutput',false);

%%% 2D projection - integrated
% nden2: 2D projected density (unit like m^-2), Npulse X 3(dims Z,X,Y
% int'd) cell-array
nden2=cell(pal_nseq,3);
for ii=1:pal_nseq
    for jj=1:3
        nden2{ii,jj}=squeeze(mean(nden3{ii},jj))*(edges{jj}(end)-edges{jj}(1));    % 2D density
    end
end

% summarise
if verbose>0
    % each dim plotted in separate figure
    plot_ncol=ceil(sqrt(pal_nseq));
    plot_nrow=ceil(pal_nseq/plot_ncol);
    
    for ii=1:3
        hfig_pal_nden2(ii)=figure();
        
        % get image X,Y axis - NOT CYCLIC
        ord=[1,2,3];
        ord_xy=ord([1:ii-1,ii+1:3]);    % pop the integrated dim out
        
        % plot each PAL as subfigure
        for jj=1:pal_nseq
            subplot(plot_nrow,plot_ncol,jj);
            
            % quick and dirty way to do density plot - imagesc
            imagesc(cents{ord_xy(2)},cents{ord_xy(1)},nden2{jj,ii});
            set(gca,'YDir','normal');   % orient it the correct way
            
            % annotate
            axis equal;
            box on;
            ht=sprintf('(%d) %0.2g',jj,fitval.y(jj));
            title(ht);
        end
    end
end

%% 2D slices
% animation - sweep each dim and take 2D density slice
% each dim plotted in separate figure
plot_ncol=ceil(sqrt(pal_nseq));
plot_nrow=ceil(pal_nseq/plot_ncol);

mov=cell(3,1);
for ii=1:3
    hfig_pal_nden2_sl(ii)=figure();
   
    % get figure sizes
    hpos=get(hfig_pal_nden2_sl(ii),'Position');
    hwidth=hpos(3);
    hheight=hpos(4);
    % reference frame (axis zero-intercept)
    X=cents{ii}*1e3;            % X-point array for YT slice [mm]
    [~,id_ref]=min(abs(X));     % reference at X~=0
    % preallocate movie
    mov{ii}=zeros(hheight,hwidth,1,length(X));
    
    %%% set colormap limits for PAL's at reference - !overwritten for different axis
    cmaplim=cell(pal_nseq,1);
    for jj=1:pal_nseq
        % take the slice through given axis
        this_nden3_dshift=shiftdim(nden3{jj},(ii-1));   % shift 3D matrix dim to slice each direction programmatically
        this_nden2_slice=squeeze(shiftdim(this_nden3_dshift(id_ref,:,:),mod(-(ii-1),3)));   % get 2D slice through this axis, then reverse the dim shift, squeeze to 2D
        % draw the image for setting reference colormap
        imagesc(this_nden2_slice);
        % set colormap limits to this frame (auto)
        cmaplim{jj}=minmax(minmax(this_nden2_slice)');
        cmaplim{jj}=[cmaplim{jj}(1,1),cmaplim{jj}(2,2)];
    end
    clf;    % clear graphics objects used for referencing colormap, etc.
    
    %%% plot PAL 2d sliced density profile    
    % get image X,Y axis - NOT CYCLIC
    ord=[1,2,3];
    ord_xy=ord([1:ii-1,ii+1:3]);    % pop the integrated dim out
    
    % loop thru SLICING AXIS + add frame to the movie structure
    for id = 1:length(X)
        % draw this frame - 2D sliced density for each PAL
        for jj=1:pal_nseq
            subplot(plot_nrow,plot_ncol,jj);
            
            % quick and dirty way to do density plot - imagesc
            this_nden3_dshift=shiftdim(nden3{jj},(ii-1));
            this_nden2_slice=squeeze(shiftdim(this_nden3_dshift(id,:,:),mod(-(ii-1),3)));
            
            imagesc(cents{ord_xy(2)},cents{ord_xy(1)},this_nden2_slice,cmaplim{jj});
            set(gca,'YDir','normal');   % orient it the correct way
        end
        
        % get frame as an image
        f=getframe(gcf);
        mov{ii}(:,:,1,id)=rgb2ind(f.cdata,256);     % greyscale
    end
    
    % create animated GIF
    t_gif=3;    % duration of gif [s]
    
    % save gif
    if configs.flags.savedata
        imwrite(mov{ii},[configs.files.dirout,'/density_2d_slice',int2str(ii),'.gif'], 'DelayTime', t_gif/length(X), 'LoopCount', inf);
    end
end
clear mov;

%% save results
if configs.flags.savedata
%     %%% fig
%     % TODO - need to be able to save all graphics from each subroutine
%     for ii=1:length(HFIG)
%         for jj=1:length(HFIG{ii})
%             saveas(HFIG{ii}{jj},[configs.files.dirout,'/',sprintf('fig_%d_%d',ii,jj),'.png']);
%             saveas(HFIG{ii}{jj},[configs.files.dirout,'/',sprintf('fig_%d_%d',ii,jj),'.fig']);
%         end
%     end
%     clear HFIG;     % clear HFIG graphics handle cell array from workspace
    
    %%% data
    % get all vars in workspace except graphics handles
    allvars=whos;
    tosave=cellfun(@isempty,regexp({allvars.class},'^matlab\.(ui|graphics)\.'));
    
    % doesn't save the whole workspace this way
    save([configs.files.dirout,'/',mfilename,'_data','.mat'],allvars(tosave).name);
end

%% end of code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_main_end=toc(t_main_start);   % end of code
disp('-----------------------------------------------');
fprintf('Total elapsed time (s): %7.1f\n',t_main_end);
disp('===================ALL TASKS COMPLETED===================');