%% Combine fringe spacing data from different experiments
% load mat files generated by main_shockfringe.m from various experiments
%

%% CONFIG
Ndpeakplot=5;    % number of fringe spacings to plot

savefigs=0;
path_save='C:\Users\HE BEC\Desktop\shock_summary';

path_data='C:\Users\HE BEC\Documents\lab\shockwave\summary\data\ver5';
data_regexp='run_*.mat';

datetimestr=datestr(datetime,'yyyymmdd_HHMMSS');    % timestamp when function called


%%% plot
linewidth=1.5;


%% load data
% get directory
mlist=dir(fullfile(path_data,data_regexp));
mlist={mlist.name};
Nexp=length(mlist);     % number of datasets

S=cell(Nexp,1);     % S is cell array of structs
for ii=1:Nexp
    varsummary={'N_peak_max', 'lambda_ff', 'lambda_ff_err', 'Nal', 'Nal_err_tot', 'N0', 'N0_err_fit',...
        'v','c','lambda_nf'};
    S{ii}=load(fullfile(path_data,mlist{ii}),varsummary{:});
end
% tidy structure
S=cell2mat(S);      % S is now struct array

%% initialise
Npeakmax=max([S.N_peak_max]);

% plot
cc=distinguishable_colors(Npeakmax);
mm={'o','s','^','x','.','*','+','d','v','>','<','p','h'};


%% plot - fringe spacing vs Nal
hfig_dpeak_vs_Npal=figure();

namearray={'LineWidth','MarkerFaceColor'};      % error bar graphics properties
valarray={linewidth,'w'};                 % 90 deg (normal) data

for ii=1:Nexp
    thisS=S(ii);
    ndpeak=(thisS.N_peak_max-1);
    if Ndpeakplot<ndpeak
        ndpeak=Ndpeakplot;
    end
%     p=zeros(1,ndpeak);
    for jj=1:ndpeak
        hold on;
%         hdata_pal_n=ploterr(thisS.Nal,thisS.lambda_ff(:,jj),thisS.AL_N_SD,thisS.lambda_ff_err(:,jj),mm{ii},'hhxy',0);
hdata_pal_n=ploterr(thisS.Nal,1e3*thisS.lambda_ff(:,jj),thisS.Nal_err_tot,1e3*thisS.lambda_ff_err(:,jj),mm{ii},'hhxy',0);
        set(hdata_pal_n(1),namearray,valarray,'Color',cc(jj,:),'MarkerSize',6,'DisplayName',sprintf('%d',jj));
        set(hdata_pal_n(2),namearray,valarray,'Color',cc(jj,:),'DisplayName','');
        set(hdata_pal_n(3),namearray,valarray,'Color',cc(jj,:),'DisplayName','');
%         p(jj)=hdata_pal_n(1);
    end
end
box on;
% lgd=legend(p);
% title(lgd,'Fringe spacing');
xlabel('$N_{AL}$');
ylabel('Fringe spacing [mm]');

%%% save fig
if savefigs>0
    figname=sprintf('dpeak_vs_Npal_%s',datetimestr);
    saveas(hfig_dpeak_vs_Npal,fullfile(path_save,[figname,'.png']));
    saveas(hfig_dpeak_vs_Npal,fullfile(path_save,[figname,'.fig']));
end

%% plot - fringe spacing vs N0
hfig_dpeak_vs_N0=figure();

namearray={'LineWidth','MarkerFaceColor'};      % error bar graphics properties
valarray={linewidth,'w'};                 % 90 deg (normal) data

for ii=1:Nexp
    thisS=S(ii);
    ndpeak=(thisS.N_peak_max-1);
    if Ndpeakplot<ndpeak
        ndpeak=Ndpeakplot;
    end
%     p=zeros(1,ndpeak);
    for jj=1:ndpeak
        hold on;
        hdata_pal_n=ploterr(thisS.N0,1e3*thisS.lambda_ff(:,jj),thisS.N0_err_fit,1e3*thisS.lambda_ff_err(:,jj),mm{ii},'hhxy',0);
        set(hdata_pal_n(1),namearray,valarray,'Color',cc(jj,:),'MarkerSize',6,'DisplayName',sprintf('%d',jj));
        set(hdata_pal_n(2),namearray,valarray,'Color',cc(jj,:),'DisplayName','');
        set(hdata_pal_n(3),namearray,valarray,'Color',cc(jj,:),'DisplayName','');
%         p(jj)=hdata_pal_n(1);
    end
end
box on;
% lgd=legend(p);
% title(lgd,'Fringe spacing');
xlabel('$N_{0}$');
ylabel('Fringe spacing [mm]');

%%% save fig
if savefigs>0
    figname=sprintf('dpeak_vs_N0_%s',datetimestr);
    saveas(hfig_dpeak_vs_N0,fullfile(path_save,[figname,'.png']));
    saveas(hfig_dpeak_vs_N0,fullfile(path_save,[figname,'.fig']));
end

%% scaling
% % try:
% %%% density (avg)
% % TODO - this changes in AL - currently assumes uniform density
% % nal_sc ~ N_al ^ alpha * N_BEC ^ beta
% alpha=1;
% beta=-3/5;
% 
% %%% speed of sound (scaled) (avg)
% % TODO - this changes in AL
% % c_sc ~ nal_sc ^ gamma
% gamma=1/2;
% 
% %%% collision speed (avg)
% % TODO - calculated for fringes!
% % v_sc ~ CONST (FOR NOW)
% 
% %%% mach number (avg)
% % mach_sc ~ v_sc / c_sc
% 
% %%% v (flow velocity) exp
% % v_sc ~ radius_of_AL / TOF
% 
% %%% c (speed of sound scaled) exp
% % c_exp ~ 4.2 e-12 * sqrt( Nal / vol_al(t: g*t=vmax))
% 
% for ii=1:Nexp
%     % density
%     S(ii).nal_sc=(S(ii).Nal'.^alpha).*(S(ii).N0.^beta);    
%     S(ii).c_sc=S(ii).nal_sc.^gamma;        % scaled C
%     S(ii).v_sc=1;
% %     S(ii).mach_sc=(S(ii).N0.^(3/10)).*(S(ii).Nal'.^(-1/2));
%     S(ii).mach_sc=S(ii).v_sc./S(ii).c_sc;
% end
% 
% hfig_dpeak_scaled=figure();
% 
% namearray={'LineWidth','MarkerFaceColor'};      % error bar graphics properties
% valarray={linewidth,'w'};                 % 90 deg (normal) data
% 
% for ii=1:Nexp
%     thisS=S(ii);
%     ndpeak=(thisS.N_peak_max-1);
%     if Ndpeakplot<ndpeak
%         ndpeak=Ndpeakplot;
%     end
%     %     p=zeros(1,ndpeak);
%     for jj=1:ndpeak
%         hold on;
%         hdata_pal_n=ploterr(thisS.mach_sc,thisS.lambda_ff(:,jj).*S(ii).c_sc,[],[],mm{ii},'hhxy',0);
% %                 hdata_pal_n=ploterr(thisS.nal_sc,thisS.lambda_ff(:,jj).*S(ii).c_sc,[],[],mm{ii},'hhxy',0);
% %         hdata_pal_n=ploterr(thisS.nal_sc,thisS.lambda_ff(:,jj),[],[],mm{ii},'hhxy',0);
%         set(hdata_pal_n(1),namearray,valarray,'Color',cc(jj,:),'MarkerSize',6,'DisplayName',sprintf('%d',jj));
% %         set(hdata_pal_n(2),namearray,valarray,'Color',cc(jj,:),'DisplayName','');
% %         set(hdata_pal_n(3),namearray,valarray,'Color',cc(jj,:),'DisplayName','');
% %         p(jj)=hdata_pal_n(1);
%     end
% end
% box on;
% % lgd=legend(p);
% % title(lgd,'Fringe spacing');
% % xlabel('$\tilde{\rho}_{AL}$');
% % ylabel('Reduced fringe spacing');
% xlabel('inverse Mach number');
% ylabel('$\lambda_nf \cdot c_sc$');
% % 
% % %%% save fig
% % if savefigs>0
% %     figname=sprintf('dpeak_vs_N0_%s',datetimestr);
% %     saveas(hfig_dpeak_vs_N0,fullfile(path_save,[figname,'.png']));
% %     saveas(hfig_dpeak_vs_N0,fullfile(path_save,[figname,'.fig']));
% % end


%% theory - Bogoliubov-Cherenkov condition
cc=distinguishable_colors(Npeakmax);

% define constants
m=6.647e-27;
hbar=1.055e-34;

hfig_dpeak_theory=figure();

namearray={'LineWidth','MarkerFaceColor'};      % error bar graphics properties
valarray={linewidth,'w'};                 % 90 deg (normal) data

% BUILD DATA
clearvars S_new;
for ii=1:Nexp
    thisS=S(ii);
    ndpeak=(thisS.N_peak_max-1);
    if Ndpeakplot<ndpeak
        ndpeak=Ndpeakplot;
    end
%         p=zeros(1,ndpeak);
    for jj=1:ndpeak
        hold on;
        thisS.X(:,jj)=(2*m/hbar*sqrt((thisS.v(:,jj)).^2-(thisS.c(:,jj)).^2));
        thisS.Y(:,jj)=(2*pi)./thisS.lambda_nf(:,jj);
        
        % summary - collate peaks across different ALs in same run
        thisS.X_collate(jj)=mean(thisS.X(:,jj),'omitnan');
        thisS.X_collate_err(jj)=std(thisS.X(:,jj),'omitnan');
        thisS.Y_collate(jj)=mean(thisS.Y(:,jj),'omitnan');
        thisS.Y_collate_err(jj)=std(thisS.Y(:,jj),'omitnan');
    end
    
    % Errors
    % X-error = 1/2 * relerr(density @ r) + 3/2 * relerr(R AL)
    % N.B. here I assume 3.3% (same as Nal,N0) for each of those relative uncertainties above
    relerr=3.3e-2;
    thisS.Xerr=norm([1/2*relerr,3/2*relerr]).*thisS.X;
    % Y-error = relerr fringe spacing
    thisS.Yerr=(thisS.lambda_ff_err(:,1:ndpeak)./thisS.lambda_ff(:,1:ndpeak)).*thisS.Y;
    
    % update collate errors
    % error as max from stdev and sample unc
    thisS.X_collate_err_tot=max([thisS.X_collate_err;mean(thisS.Xerr,1,'omitnan')],[],1);
    thisS.Y_collate_err_tot=max([thisS.Y_collate_err;mean(thisS.Yerr,1,'omitnan')],[],1);
%     % error as mean from stdev and sample unc
%     thisS.X_collate_err_tot=mean([thisS.X_collate_err;mean(thisS.Xerr,1,'omitnan')],1,'omitnan');
%     thisS.Y_collate_err_tot=mean([thisS.Y_collate_err;mean(thisS.Yerr,1,'omitnan')],1,'omitnan');
    
    S_new(ii)=thisS;
end

% PLOTTING
for ii=1:Nexp
    thisS=S_new(ii);
    ndpeak=(thisS.N_peak_max-1);
    if Ndpeakplot<ndpeak
        ndpeak=Ndpeakplot;
    end
    
    for jj=1:ndpeak
        hdata_theory=ploterr(thisS.X(:,jj),thisS.Y(:,jj),thisS.Xerr(:,jj),thisS.Yerr(:,jj),mm{ii},'hhxy',0);
        set(hdata_theory(1),namearray,valarray,'Color',cc(jj,:),'MarkerSize',6,'DisplayName',sprintf('%d',jj));
        set(hdata_theory(2),namearray,valarray,'Color',cc(jj,:),'DisplayName','');
        set(hdata_theory(3),namearray,valarray,'Color',cc(jj,:),'DisplayName','');
%                 p(jj)=hdata_theory(1);
    end
end
box on;
axis equal;
axis square;
% lgd=legend(p);
% title(lgd,'Fringe spacing');
xlabel('$2 m / \hbar \cdot (v^2 - c^2)^{1/2}$ [m$^{-1}$]');
ylabel('$2 \pi / \lambda_{NF} $ [m$^{-1}$]');


%% collate spread from each run
% spread from each "fringe"
hfig_dpeak_theory_summ=figure();
hold on;

namearray={'LineWidth','MarkerFaceColor'};      % error bar graphics properties
valarray={linewidth,'w'};                 % 90 deg (normal) data

cc=distinguishable_colors(Nexp);
p=[];
for ii=1:Nexp
    thisS=S(ii);
    hdata_theory_summ=ploterr(S_new(ii).X_collate,S_new(ii).Y_collate,...
        S_new(ii).X_collate_err_tot,S_new(ii).Y_collate_err_tot,...
        mm{ii},'hhxy',0);
    %     set(hdata_theory_summ(1),namearray,valarray,'Color',cc(ii,:),'MarkerSize',6,'DisplayName',mlist{ii});
    set(hdata_theory_summ(1),namearray,valarray,'Color',cc(ii,:),'MarkerSize',6,'DisplayName',sprintf('%d',ii));
    set(hdata_theory_summ(2),namearray,valarray,'Color',cc(ii,:),'DisplayName','');
    set(hdata_theory_summ(3),namearray,valarray,'Color',cc(ii,:),'DisplayName','');
    p(ii)=hdata_theory_summ(1);
end
box on;
axis equal;
axis square;
lgd=legend(p,'Location','northwest');
title(lgd,'Exp');
xlabel('$2 m / \hbar \cdot (v^2 - c^2)^{1/2}$ [m$^{-1}$]');
ylabel('$2 \pi / \lambda_{NF} $ [m$^{-1}$]');