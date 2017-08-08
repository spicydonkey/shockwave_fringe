clear all;

%% configs
configs.flags.verbose=1;
configs.flags.savedata=0;
configs.flags.archive_txy=1;
configs.flags.force_all_stages=0;
configs.flags.graphics=1;
configs.flags.build_txy=1;

configs.misc.vel_z=9.8*0.416;    % atom free-fall vert v at detector hit for T-to-Z conversion;

configs.files.path='C:\Users\HE BEC\Documents\lab\shockwave\20170716_atomlaser\d';

configs.files.dir_data=fileparts(configs.files.path);    % fullpath to data directory
configs.files.archive=fullfile(configs.files.dir_data,'archive');   % dir to archive folder
configs.files.dirout=fullfile(configs.files.dir_data,'output');      % output directory (will be time-stamped)

configs.load.version=1;         % TXY load stage version number

configs.load.id=1:3615;         % file id numbers to use for analysis
configs.load.mincount=0;         % min counts in window - 0 for no min
configs.load.maxcount=Inf;          % max counts in window - Inf for no max

configs.load.rot_angle=0.61;    % detector/trap alignment

% construct window for ROI
configs.load.window{1}=[0.46,0.481];      % T [s]
configs.load.window{2}=[-5e-3,10e-3];    % X [m]
configs.load.window{3}=[-20e-3,30e-3];    % Y [m]

configs.image.voxel_res=1e-4*[0.33,1,1];   % ZXY voxel resolution [m]
configs.image.size=[-30e-3,40e-3;-5e-3,5e-3; -25e-3,25e-3];   % image size/lims [T;X;Y] [m]

%% initialise
do_next=configs.flags.force_all_stages;
verbose=configs.flags.verbose;
t_main_start=tic;   % for reporting process duration
datetimestr=datestr(datetime,'yyyymmdd_HHMMSS');    % timestamp when function called
configs.files.dirout=[configs.files.dirout,'_',datetimestr];
vz=configs.misc.vel_z;
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
    hfig_pal=figure();
    plot_zxy(zxy,3e4,1,'k');
    axis equal;
    view(90,0);
end

% centre PAL to mean position (in fact this isn't used since it reduces
% fringe visibility!)
zxy_cent=cellfun(@(x) mean(x,1), zxy, 'UniformOutput', false);
zxy_cent=vertcat(zxy_cent{:});
zxy_cent_std=std(zxy_cent,1);

zxy0=cell(nshot,1);
for ii=1:nshot
    zxy0{ii}=zxy{ii}-repmat(zxy_cent(ii,:),[size(zxy{ii},1),1]);
end

% check centered
if verbose>0
    hfig_pal_cent=figure();
    plot_zxy(zxy0,3e4,1,'k');
    axis equal;
    view(90,0);
end

% study in oscillation cancellation with fringe visibility
% does oscillation cancellation affect fringe visibility? YES!
zxy_cent_avg=mean(zxy_cent,1);      % average zxy center in PAL
zxy_avg_shifted=cellfun(@(x) x-repmat(zxy_cent_avg,[size(x,1),1]),zxy,'UniformOutput',false);

%% PAL density image
% construct edge/center vectors for each dim
edges=cell(3,1);
cents=cell(3,1);
for ii=1:3
    edges{ii}=configs.image.size(ii,1):configs.image.voxel_res(ii):configs.image.size(ii,2);
    cents{ii}=0.5*(edges{ii}(1:end-1)+edges{ii}(2:end));
end

% density profile from all shots collated 
[ncounts_pal,~]=histcn(vertcat(zxy_avg_shifted{:}),edges{1},edges{2},edges{3});
vvoxel=prod(configs.image.voxel_res);       % volume of a voxel [m^3]
density_pal=ncounts_pal/(vvoxel*nshot);     % density in voxel [m^-3]

%% summarise density profiles
% density in 2D projection
density_2d=cell(3,1);   % 2D integrated XY, TY, TX - plane density profiles
for ii=1:3
    density_2d{ii}=squeeze(sum(density_pal,ii));
end

if verbose>0
    for ii=1:3
        figure();
        imagesc(density_2d{ii});
        colorbar();
    end
end

%% 2D slices animation - through X-axis
if verbose>0
    hfig_density_2d_slice=figure();
    
    % define objects at the reference frame
    X=cents{2}*1e3;     % X-point array for YT slice [mm]
    [~,id_ref]=min(abs(X));     % reference at X~=0
    
    im=imagesc(squeeze(density_pal(:,id_ref,:)));
    ht=title(sprintf('X = %0.2f mm', X(id_ref)));
    
    % set colormap to the reference frame (auto)
    f = getframe(gcf);
    [~,cmap]=rgb2ind(f.cdata, 256, 'nodither');
    
    % get figure sizes
    pos=get(hfig_density_2d_slice,'Position');
    width=pos(3);
    height=pos(4);
    
    % preallocate
    mov = zeros(height, width, 1, length(X), 'uint8');
    
    % Animate and add animation frame to the movie structure
    for id = 1:length(X)
        % Update XData and YData
        set(im,'CData',squeeze(density_pal(:,id,:)));
        set(ht, 'String', sprintf('X = %0.2f mm', X(id)));
        
        % Get frame as an image
        f = getframe(gcf);
        mov(:,:,1,id) = rgb2ind(f.cdata, cmap, 'nodither');
    end
    
    % Create animated GIF
    t_gif=3;    % duration of gif
    imwrite(mov, cmap, 'density_2d_slice.gif', 'DelayTime', t_gif/length(X), 'LoopCount', inf);
end


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