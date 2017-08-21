%% clean workspace
clear all; close all; clc;

%% configs
% User path to config
% path_config='C:\Users\HE BEC\Documents\MATLAB\shockwave_fringe\configs\config_20170716_atomlaser.m';
path_config='C:\Users\HE BEC\Documents\MATLAB\shockwave_fringe\configs\config_20170717_atomlaser.m';
% path_config='C:\Users\HE BEC\Documents\MATLAB\shockwave_fringe\configs\config_run1.m';
% path_config='C:\Users\HE BEC\Documents\MATLAB\shockwave_fringe\configs\config_run2.m';

% load config
run(path_config);

%% initialise
% flags
do_next=configs.flags.force_all_stages;
verbose=configs.flags.verbose;
vgraph=configs.flags.graphics;

% misc
t_main_start=tic;   % for reporting process duration
datetimestr=datestr(datetime,'yyyymmdd_HHMMSS');    % timestamp when function called
configs.files.dirout=[configs.files.dirout,'_',datetimestr];

HFIG={};        % cell array for all figures generated

% experimental consts/params
vz=configs.misc.vel_z;
det_qe=configs.misc.det_qe;
tof=configs.misc.tof;
hbar=configs.misc.hbar;
m=configs.misc.m;
g=configs.misc.g;

%%% PAL
pal_z1=configs.pal.t1*vz;
pal_dz=configs.pal.dt*vz;
pal_nseq=configs.pal.n;

plot_ncol=ceil(sqrt(pal_nseq));         % subplot num cols
plot_nrow=ceil(pal_nseq/plot_ncol);     % subplot num rows

%%% voxels
% construct edge/center vectors for each dim
edges=cell(3,1);
cents=cell(3,1);
for ii=1:3
    edges{ii}=configs.image.size(ii,1):configs.image.voxel_res(ii):configs.image.size(ii,2);
    cents{ii}=0.5*(edges{ii}(1:end-1)+edges{ii}(2:end));
end

%%% directories
% output
if configs.flags.savedata
    % output directory
    if ~isdir(configs.files.dirout)
        warning(['output directory "',configs.files.dirout,'" does not exist. Creating directory...']);
        mkdir(configs.files.dirout);
    end
end

% archive
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
clearvars S_temp this_file arch_mat_list;       % clean workspace

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

clearvars loadconfigs fpath_txy_archive;  % clean workspace

%% pre-process
nshot=size(txy,1);

% txy-->zxy
zxy=cellfun(@(x) double(x.*repmat([vz,1,1],[size(x,1),1])),txy,'UniformOutput',false);
clearvars txy;      % clear the biggest data from memory

% check loaded data
if vgraph>0
    if verbose>2
        hfig_all=figure();
        plot_zxy(zxy,3e4,1,'k');
        axis equal;
        view(90,0);
        
        if configs.flags.savedata
            figname=sprintf('raw_data');
            saveas(hfig_all,[configs.files.dirout,'/',figname,'.png']);
        end
    end
end

% get PAL
pal_zxy=capture_pal(zxy,pal_z1,pal_dz,pal_nseq);
clearvars zxy;      % delete zxy - not used from here

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
    pal_cent_std(ii,:)=std(this_pal_cent_array);        % fluctuations in the mean centre position
end

if verbose>0
    % summary
    disp('------------------------PAL centering------------------------');
    for ii=1:pal_nseq
    fprintf('(%d) STDEV(CENTER): [%0.2e %0.2e %0.2e]\n',ii,pal_cent_std(ii,:));
    end
    disp('-------------------------------------------------------------');
end

% centre PAL
pal_zxy0=cell(pal_nseq,1);     % preallocate centered PAL array
for ii=1:pal_nseq
    % shift to evaluated centre
    pal_zxy0{ii}=cellfun(@(x) x-repmat(pal_cent_avg(ii,:),[size(x,1),1]),pal_zxy{ii},'UniformOutput',false);
end
clearvars pal_zxy;

% summarise captured PALs
if vgraph>0
    if verbose>2
        hfig_pal=figure();
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
        
        if configs.flags.savedata
            figname=sprintf('al_zxy');
            saveas(hfig_pal,[configs.files.dirout,'/',figname,'.png']);
        end
    end
end

%% Process PAL
%%% Evaluate atom numbers
% number in PAL
num_in_pal=zeros(pal_nseq,2);    % preallocate; format: [avg, std]
for ii=1:pal_nseq
    this_num_in_pal=cellfun(@(x) size(x,1),pal_zxy0{ii});    % number in nth PAL 'detected'
    % estimate actual number in PAL - error in SD
    num_in_pal(ii,:)=(1/det_qe)*[mean(this_num_in_pal),std(this_num_in_pal)];
end
% PAL number uncertainty in SE
num_in_pal(:,2)=num_in_pal(:,2)/sqrt(nshot);

%%%% fit to estimate number in condensate and outcoupling efficiency
% MODEL:
%   x1: PAL index
%   N0: number in BEC before x1^th PAL
%   r: AL outcoupling frac
%   y: number in x1^th PAL; 
al_modelfun='y~N0*(1-r)^(x1-1)*r';     
coeffnames={'N0','r'};

idx_al_fit=configs.pal.nfitstart:pal_nseq;     % pulses to do fit
param0=[1e5,0.1];   % estimate for initial param

% define optimiser
fo = statset('TolFun',10^-6,...
    'TolX',10^-6,...
    'MaxIter',10^6,...
    'UseParallel',0);

% do the fit
N_al_fit=fitnlm(idx_al_fit,num_in_pal(idx_al_fit,1),al_modelfun,param0,...
    'CoefficientNames',coeffnames,'Options',fo);

if verbose>0
    disp(N_al_fit);
end

% get fit results
al_paramfit=[N_al_fit.Coefficients.Estimate,N_al_fit.Coefficients.SE];
Npal_fit.x=[1:pal_nseq]';
Npal_fit.y=feval(N_al_fit,Npal_fit.x);
Nal=Npal_fit.y;

N0_0=al_paramfit(1,1);      % original number in BEC
eff_al=al_paramfit(2,1);    % AL outcoupling ratio

%%% atoms numbers error analysis
% BEC
N0=(N0_0*(1-eff_al).^(1:pal_nseq))';
fit_err_rel=N_al_fit.Coefficients.SE./N_al_fit.Coefficients.Estimate;        % N0, r fit err (rel)
N0_err_fit_rel=fit_err_rel(2)*(sqrt(1:length(N0))');    % rel error for N0 - contribution from outcoupling eff
N0_err_fit_rel=sqrt(N0_err_fit_rel.^2+fit_err_rel(1)^2);    % add contribution from unc in init N0
N0_err_fit=N0.*N0_err_fit_rel;        % evaluate fit SE error (abs)

% AL
Nal_err_fit_rel=N0_err_fit_rel';  	% formula for Nal scales identically with N0
Nal_err_fit=Nal.*Nal_err_fit_rel;	% evaluate SE from fit uncertainties
% later summed in quad with detected num SE for total err

clearvars Nfitrelerr N0_err_fit_rel Nal_err_fit_rel;

if vgraph>0
    % plot AL model fit
    linewidth=1.5;
    namearray={'LineWidth','MarkerFaceColor','Color'};      % error bar graphics properties
    valarray={linewidth,'w','k'};                 % 90 deg (normal) data
    if verbose>1
        hfig_atom_number=figure();
        % plot data
        hdata_pal_n=ploterr(1:pal_nseq,num_in_pal(:,1),[],num_in_pal(:,2),'o','hhxy',0);
        set(hdata_pal_n(1),namearray,valarray,'DisplayName','Data');
        set(hdata_pal_n(2),namearray,valarray,'DisplayName','');
        
        hold on;
        % plot fit
        hfit=plot(Npal_fit.x,Nal,'r*','DisplayName','Fit');
        
        % annotate
        legend([hdata_pal_n(1),hfit]);
        xlabel('Pulse number');
        ylabel('Number in PAL');
        
        if configs.flags.savedata
            figname=sprintf('N_al_fit');
            saveas(hfig_atom_number,[configs.files.dirout,'/',figname,'.png']);
        end
    end
end

%% PAL density image
if verbose>0
    %%% 3D density profile (full)
    nden3=cellfun(@(x) density3d(x,edges),pal_zxy0,'UniformOutput',false);
    
    %%% 2D projection - projected
    % nden2: 2D projected density (unit like m^-2), Npulse X 3(dims Z,X,Y
    % int'd) cell-array
    nden2=cell(pal_nseq,3);
    for ii=1:pal_nseq
        for jj=1:3
            nden2{ii,jj}=squeeze(mean(nden3{ii},jj))*(edges{jj}(end)-edges{jj}(1));    % 2D density
        end
    end
    
    % summarise
    if vgraph>0
        if verbose>1
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
                    ht=sprintf('(%d) %0.2g / %0.2g',jj,Nal(jj),N0(jj));
                    title(ht);
                end
                
                if configs.flags.savedata
                    figname=sprintf('al_2D_%d',ii);
                    saveas(hfig_pal_nden2(ii),[configs.files.dirout,'/',figname,'.png']);
                    saveas(hfig_pal_nden2(ii),[configs.files.dirout,'/',figname,'.fig']);
                end
            end
        end
    end
end

%% 2D slices
if vgraph>1
    % animation - sweep each dim and take 2D density slice    
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
    % free up memory
    clearvars mov;
end
clearvars nden3;


%% Characterise AL
% characterise AL by a circular disc assumption: radius, width
cc=distinguishable_colors(pal_nseq);	% color set for plotting each AL

% R density distribution (YZ plane)
r_edge=linspace(0,20e-3,100);
r_cent=r_edge(1:end-1)+0.5*diff(r_edge);

% transverse dist (X-axis projection)
x_edge=linspace(0,5e-3,100);   % single sided
x_cent=x_edge(1:end-1)+0.5*diff(x_edge);

% preallocate
nden_r=zeros(pal_nseq,length(r_cent));     % radial density [dimension TODO]
pal_R=zeros(pal_nseq,1);        % charactersitic radius of AL [m]
nden_x=zeros(pal_nseq,length(x_cent));	% atom laser density along thickness [dimension TODO]
pal_Rx=zeros(pal_nseq,1);	% atom laser thickness (half-width) [m]

if vgraph>0
    hfig_al_rden=figure();
    hold on;
    hfig_al_xden=figure();
    hold on;
end
for ii=1:pal_nseq
    %%% get radial density profile
    this_nden_r=al_rad_density(vertcat(pal_zxy0{ii}{:}),r_edge);
    this_nden_r=this_nden_r/nshot;    % normalise by number of shots collated

    %%% get perp density profile
    x=vertcat(pal_zxy0{ii}{:});    % collate all shots in this PAL
    x=abs(x(:,2));     % cull Y,Z, sign on X (X-symm)
    
    this_nden_x=histcounts(x,x_edge);     % get histogram counts
    this_nden_x=this_nden_x/(nshot*2);    % norm - shots + X symmetry
    this_nden_x=this_nden_x./diff(x_edge);  % hist-->density
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get characteristic size (half maximum)
    %%% R
    [nr_max,i_max]=max(this_nden_r);
    [~,i_halfmax]=min(abs(this_nden_r(i_max:end)-nr_max/2));
    pal_R(ii)=r_cent(i_max+i_halfmax-1);      % get half maximum radius and store
    
    %%% Perp
    [nx_max,i_max]=max(this_nden_x);
    [~,i_halfmax]=min(abs(this_nden_x(i_max:end)-nx_max/2));
    pal_Rx(ii)=x_cent(i_max+i_halfmax-1);  % half maximum width
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % normalise density by the integration volume (other dims)
    this_nden_r=this_nden_r/(2*pal_Rx(ii));
    this_nden_x=this_nden_x/(pi*pal_R(ii)^2);
    
    nden_r(ii,:)=this_nden_r;
    nden_x(ii,:)=this_nden_x;
    
    if vgraph>0
        % plot
        figure(hfig_al_rden);
        plot(r_cent,nden_r(ii,:),'Color',cc(ii,:));     % profile
        scatter(pal_R(ii),max(this_nden_r)/2,'MarkerEdgeColor',cc(ii,:));   % half maximum point
        
        figure(hfig_al_xden);
        plot(x_cent,nden_x(ii,:),'Color',cc(ii,:));     % profile
        scatter(pal_Rx(ii),max(this_nden_x)/2,'MarkerEdgeColor',cc(ii,:));   % half maximum point
    end
end
if vgraph>0
    % annotate plot
    figure(hfig_al_rden);
    title('AL radial density profile');
    box on;
    xlabel('Radius [m]');
    ylabel('density [m$^{-3}$]');
    if configs.flags.savedata
        figname=sprintf('al_R_size');
        saveas(hfig_al_rden,[configs.files.dirout,'/',figname,'.png']);
        saveas(hfig_al_rden,[configs.files.dirout,'/',figname,'.fig']);
    end
    
    figure(hfig_al_xden);
    title('AL X density profile');
    box on;
    xlabel('X [m]');
    ylabel('density [m$^{-3}$]');
    if configs.flags.savedata
        figname=sprintf('al_Rx_size');
        saveas(hfig_al_xden,[configs.files.dirout,'/',figname,'.png']);
        saveas(hfig_al_xden,[configs.files.dirout,'/',figname,'.fig']);
    end
end

% % evaluate shot-to-shot variability
% % coarse edges for single shot
% r_edge_shot=linspace(0,20e-3,20);
% r_cent_shot=r_edge_shot(1:end-1)+0.5*diff(r_edge_shot);
% 
% % figure(); hold on;
% shot_pal_R=cell(pal_nseq,1);
% for ii=1:pal_nseq
%     shot_pal_R{ii}=zeros(nshot,1);
%     for jj=1:nshot
%         this_nden_r=al_rad_density(pal_zxy0{ii}{jj},r_edge_shot);
%         this_nden_r=smooth(this_nden_r,5);    % smooth out noise
%         
%         % get AL radius - HWHM
%         % get peak density and location
%         [nr_max,i_max]=max(this_nden_r);
%         % get half-maximum point
%         [~,i_halfmax]=min(abs(this_nden_r(i_max:end)-nr_max/2));
%         shot_pal_R{ii}(jj)=r_cent_shot(i_max+i_halfmax-1);      % get half maximum radius and store
%         
% %         if rand()>0.999
% %             plot(r_cent_shot,this_n_r);     % profile
% %             scatter(r_cent_shot(i_max+i_halfmax-1),this_n_r(i_max+i_halfmax-1));   % half maximum point
% %         end
%     end
% end

%%% evaluate AL volume
pal_vol=pi*pal_R.^2*2.*pal_Rx;
pal_nden_exp=Nal./pal_vol;

% figure();
% plot(1:pal_nseq,pal_nden_exp,'o');
% xlabel('AL number');
% ylabel('AL density from experiment (arb. unit)');

%% Check AL radius - N0, N_AL dependency
% R: very good!
hfig_Ral_N=figure();

rerr=(N0_err_fit./N0);
scaleexp=1/5;
ploterr(N0.^(scaleexp),pal_R,rerr*(scaleexp),[],'hhxy',0);

xlabel('$N_{0}^{1/5}$');
ylabel('R HWHM$_{AL} [m]$');

% X
hfig_Xal_N=figure();

scaleexp=1/5;
ploterr(N0.^(scaleexp),pal_Rx,rerr*(scaleexp),[],'hhxy',0);

xlabel('$N_{0}^{1/5}$');
ylabel('X HWHM$_{AL} [m]$');

%% Characterise shockwaves
pal_data=pal_zxy0;      % a copy to safely pass data to analysis scripts (not functions!) 

%%% process PAL
main_process_pal;

%%% fringe characterisation
DemoFringePeak;

% store results - overwritten in current implementation of bootstrap
Rff_peak=1e-3*peak_pos;     % peak positions (radial) [m]
lambda_ff=1e-3*peak_diff;        % peak spacing [m]
N_peak_max=max_peak_n;
cc_all=distinguishable_colors(N_peak_max-1);

%% evaluate uncertainty by bootstrapping
% configure for bootstrapping
vgraph=0;   % skip graphics

% bootstrapping
bootstrap_Nsubset=ceil(nshot*bootstrap_ndata);      % number of shots in a subset
I_subset=cell(bootstrap_Nsamp,1);

% analysis on subsets
lambda_ff_subset_cell=cell(bootstrap_Nsamp,1);
Nal_subset=zeros(pal_nseq,bootstrap_Nsamp);
for ii_bootstrap=1:bootstrap_Nsamp
    % random selection of subset
    I_subset{ii_bootstrap}=randperm(nshot,bootstrap_Nsubset);
    pal_data=cellfun(@(x) x(I_subset{ii_bootstrap}),pal_zxy0,'UniformOutput',false);
    
    % run analysis
    main_process_pal;
    DemoFringePeak;
    
    % get avg number in PAL from this subset
    Nal_subset(:,ii_bootstrap)=cellfun(@(x) mean(cellfun(@(y) size(y,1),x)),pal_data);
    
    % save result separately
    lambda_ff_subset_cell{ii_bootstrap}=peak_diff;
end
% evaluate SD of analysis outupt
% restructure subset outputs into array
lambda_ff_sub=NaN([size(lambda_ff),bootstrap_Nsamp]);   % preallocate NaN array
for ii_bootstrap=1:bootstrap_Nsamp
    this_peak_diff=lambda_ff_subset_cell{ii_bootstrap};
    lambda_ff_sub(:,1:size(this_peak_diff,2),ii_bootstrap)=this_peak_diff;
end

lambda_ff_err=std(lambda_ff_sub,0,3,'omitnan');    % std from bootstrapping
lambda_ff_err=lambda_ff_err(:,1:(N_peak_max-1));   % resize to main result

% PAL number uncertainty from bootstrapping data subsets
Nal_avg_sub=mean(Nal_subset,2);         % avg PAL number
Nal_err_sub=std(Nal_subset,0,2);        % stdev PAL number from bootstrap
Nal_err_tot=sqrt(Nal_err_sub'.^2+Nal_err_fit.^2);    % add in quadrature with fit uncertainty

clearvars pal_zxy0;     % clean workspace

%% Summary - plot with errors
%%% lambda vs number in atom laser
vgraph=configs.flags.graphics;   % reset graphics flag
if vgraph>0
    linewidth=1.5;
    namearray={'LineWidth','MarkerFaceColor'};      % error bar graphics properties
    valarray={linewidth,'w'};                 % 90 deg (normal) data
    
    hfig_dpeak_vs_Npal=figure();
    p=zeros(1,(N_peak_max-1));      % array to store figure objects for selective legend
    for ii=1:(N_peak_max-1)
        hold on;
        hdata_pal_n=ploterr(Nal,1e3*lambda_ff(:,ii),Nal_err_tot,lambda_ff_err(:,ii),'o','hhxy',0);
        set(hdata_pal_n(1),namearray,valarray,'Color',cc_all(ii,:),'DisplayName',sprintf('%d',ii));
        set(hdata_pal_n(2),namearray,valarray,'Color',cc_all(ii,:),'DisplayName','');
        set(hdata_pal_n(3),namearray,valarray,'Color',cc_all(ii,:),'DisplayName','');
        p(ii)=hdata_pal_n(1);
    end
    box on;
    lgd=legend(p);
    xlabel('$N_{AL}$');
    ylabel('Fringe spacing [mm]');
    
    if configs.flags.savedata
        figname=sprintf('dpeak_vs_Npal');
        saveas(hfig_dpeak_vs_Npal,[configs.files.dirout,'/',figname,'.png']);
        saveas(hfig_dpeak_vs_Npal,[configs.files.dirout,'/',figname,'.fig']);
    end
end

%%% lambda vs number in condensate
if vgraph>0
    linewidth=1.5;
    namearray={'LineWidth','MarkerFaceColor'};      % error bar graphics properties
    valarray={linewidth,'w'};                 % 90 deg (normal) data
    
    hfig_dpeak_vs_N0=figure();
    p=zeros(1,(N_peak_max-1));      % array to store figure objects for selective legend
    for ii=1:(N_peak_max-1)
        hold on;
        hdata_pal_n=ploterr(N0,1e3*lambda_ff(:,ii),N0_err_fit,lambda_ff_err(:,ii),'o','hhxy',0);
        set(hdata_pal_n(1),namearray,valarray,'Color',cc_all(ii,:),'DisplayName',sprintf('%d',ii));
        set(hdata_pal_n(2),namearray,valarray,'Color',cc_all(ii,:),'DisplayName','');
        set(hdata_pal_n(3),namearray,valarray,'Color',cc_all(ii,:),'DisplayName','');
        p(ii)=hdata_pal_n(1);
    end
    box on;
    lgd=legend(p);
    xlabel('$N_{0}$');
    ylabel('Fringe spacing [mm]');
    
    if configs.flags.savedata
        figname=sprintf('dpeak_vs_N0');
        saveas(hfig_dpeak_vs_N0,[configs.files.dirout,'/',figname,'.png']);
        saveas(hfig_dpeak_vs_N0,[configs.files.dirout,'/',figname,'.fig']);
    end
end

%% plot theory
t0=mean(pal_R/(tof*g));       % AL expansion time [s] TODO - explain avg
% t0=pal_R/(tof*g);

lambda_nf=lambda_ff.*t0/tof;	% near-field density modulation wavelength [m]
npeak=size(lambda_nf,2);

%%% approximate velocity
v_al_max=pal_R/tof;    % velocity at AL radius
% r_fringe=mean(Rff_peak./pal_R,1,'omitnan');   % ratio of distance from ith fringe to 
r_fringe=Rff_peak./pal_R;   % ratio of distance from ith fringe to 
v=v_al_max.*r_fringe;   % approx velocity at shock wave (row - AL; col - fringe)

%%% approximate speed of sound
% c=4.2e-12*sqrt(g*tof^6*Nal./(pi*pal_R.^5.*pal_Rx));     % approx speed of sound - uniform density approx

% approximate from radial density profile
% get nden_Rff from nden_r
nden_Rff=NaN(pal_nseq,npeak);     % AL far-field density at peak
for ii=1:pal_nseq
    for jj=1:npeak
        [dRff_min,iRff]=min(abs(Rff_peak(ii,jj)-r_cent));
        nden_Rff(ii,jj)=nden_r(ii,iRff);    % get density
    end
end
nden_Rnf=nden_Rff.*(tof./t0).^3;    % simple extrapolation to approx NF density
c=4.2e-12*sqrt(nden_Rnf);           % evaluate speed of sound

%%% plot
hfig_theory=figure();
for ii=1:npeak
%     plot((2*m/hbar*sqrt(v(:,ii).^2-c.^2)),(2*pi)./lambda_nf(:,ii),'o');
%     % for uniform speed of sound approx
    plot((2*m/hbar*sqrt(v(:,ii).^2-c(:,ii).^2)),(2*pi)./lambda_nf(:,ii),'o');
    hold on;
end
xlabel('$2 m / hbar \cdot (v^2 - c^2)^{1/2}$');
ylabel('$2 \pi / \lambda $ [m$^{-1}$]');

if configs.flags.savedata
    figname=sprintf('theory_curve');
    saveas(hfig_theory,[configs.files.dirout,'/',figname,'.png']);
    saveas(hfig_theory,[configs.files.dirout,'/',figname,'.fig']);
end

%% end of code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose>0
    t_main_end=toc(t_main_start);   % end of code
    disp('-----------------------------------------------');
    fprintf('Total elapsed time (s): %7.1f\n',t_main_end);
    disp('===================ALL TASKS COMPLETED===================');
end

%% save results
if configs.flags.savedata
    if verbose>0
        fprintf('Saving data. This may take a minute.\n');
    end
    
    %%% data
    % get all vars in workspace except graphics handles
    allvars=whos;
    tosave=cellfun(@isempty,regexp({allvars.class},'^matlab\.(ui|graphics)\.'));

    save([configs.files.dirout,'/',mfilename,'_data','.mat'],allvars(tosave).name,'-v7.3','-nocompression');
end