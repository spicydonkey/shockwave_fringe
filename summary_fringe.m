%% Combine fringe spacing data from different experiments
% load mat files generated by main_shockfringe.m from various experiments
%

%% CONFIG
Ndpeakplot=5;    % number of fringe spacings to plot

savefigs=0;
path_save='C:\Users\HE BEC\Desktop\shock_summary';

path_data='C:\Users\HE BEC\Documents\lab\shockwave\summary\data\ver8';
data_regexp='run_*.mat';

datetimestr=datestr(datetime,'yyyymmdd_HHMMSS');    % timestamp when function called

%%% plot
linewidth=1.5;
cgray=0.5*[1,1,1];


%% load data
% get directory
mlist=dir(fullfile(path_data,data_regexp));
mlist={mlist.name};
Nexp=length(mlist);     % number of datasets

S=cell(Nexp,1);     % S is cell array of structs
for ii=1:Nexp
    varsummary={'N_peak_max', 'lambda_ff', 'lambda_ff_err', 'Nal', 'Nal_err_tot', 'N0', 'N0_err_fit',...
        'v','c','lambda_nf', 'eff_al'};
    S{ii}=load(fullfile(path_data,mlist{ii}),varsummary{:});
end
% tidy structure
S=cell2mat(S);      % S is now struct array

%% initialise
Npeakmax=max([S.N_peak_max]);

% plot
cc=distinguishable_colors(Npeakmax);
mm={'o','s','^','d','x','.','*','+','v','>','<','p','h'};


%% plot - fringe spacing vs Nal
% hfig_dpeak_vs_Npal=figure();
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
% %     p=zeros(1,ndpeak);
%     for jj=1:ndpeak
%         hold on;
% %         hdata_pal_n=ploterr(thisS.Nal,thisS.lambda_ff(:,jj),thisS.AL_N_SD,thisS.lambda_ff_err(:,jj),mm{ii},'hhxy',0);
% hdata_pal_n=ploterr(thisS.Nal,1e3*thisS.lambda_ff(:,jj),thisS.Nal_err_tot,1e3*thisS.lambda_ff_err(:,jj),mm{ii},'hhxy',0);
%         set(hdata_pal_n(1),namearray,valarray,'Color',cc(jj,:),'MarkerSize',6,'DisplayName',sprintf('%d',jj));
%         set(hdata_pal_n(2),namearray,valarray,'Color',cc(jj,:),'DisplayName','');
%         set(hdata_pal_n(3),namearray,valarray,'Color',cc(jj,:),'DisplayName','');
% %         p(jj)=hdata_pal_n(1);
%     end
% end
% box on;
% % lgd=legend(p);
% % title(lgd,'Fringe spacing');
% xlabel('$N_{AL}$');
% ylabel('Fringe spacing [mm]');
% 
% %%% save fig
% if savefigs>0
%     figname=sprintf('dpeak_vs_Npal_%s',datetimestr);
%     saveas(hfig_dpeak_vs_Npal,fullfile(path_save,[figname,'.png']));
%     saveas(hfig_dpeak_vs_Npal,fullfile(path_save,[figname,'.fig']));
% end

%% plot - fringe spacing vs N0
% hfig_dpeak_vs_N0=figure();
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
% %     p=zeros(1,ndpeak);
%     for jj=1:ndpeak
%         hold on;
%         hdata_pal_n=ploterr(thisS.N0,1e3*thisS.lambda_ff(:,jj),thisS.N0_err_fit,1e3*thisS.lambda_ff_err(:,jj),mm{ii},'hhxy',0);
%         set(hdata_pal_n(1),namearray,valarray,'Color',cc(jj,:),'MarkerSize',6,'DisplayName',sprintf('%d',jj));
%         set(hdata_pal_n(2),namearray,valarray,'Color',cc(jj,:),'DisplayName','');
%         set(hdata_pal_n(3),namearray,valarray,'Color',cc(jj,:),'DisplayName','');
% %         p(jj)=hdata_pal_n(1);
%     end
% end
% box on;
% % lgd=legend(p);
% % title(lgd,'Fringe spacing');
% xlabel('$N_{0}$');
% ylabel('Fringe spacing [mm]');
% 
% %%% save fig
% if savefigs>0
%     figname=sprintf('dpeak_vs_N0_%s',datetimestr);
%     saveas(hfig_dpeak_vs_N0,fullfile(path_save,[figname,'.png']));
%     saveas(hfig_dpeak_vs_N0,fullfile(path_save,[figname,'.fig']));
% end

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
%     % error as max from stdev and sample unc
%     thisS.X_collate_err_tot=max([thisS.X_collate_err;mean(thisS.Xerr,1,'omitnan')],[],1);
%     thisS.Y_collate_err_tot=max([thisS.Y_collate_err;mean(thisS.Yerr,1,'omitnan')],[],1);
    % error as mean from stdev and sample unc
    thisS.X_collate_err_tot=mean([thisS.X_collate_err;mean(thisS.Xerr,1,'omitnan')],1,'omitnan');
    thisS.Y_collate_err_tot=mean([thisS.Y_collate_err;mean(thisS.Yerr,1,'omitnan')],1,'omitnan');
%     % error as from spread in samples
%     thisS.X_collate_err_tot=thisS.X_collate_err;
%     thisS.Y_collate_err_tot=thisS.Y_collate_err;
    
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


%% BCR theory - TYPE 1
% figure params
mrk_size=5;
fontsize=11;
linewidth=1.5;

papersize=[10,7];
paperposition=[0,0,papersize];


%%% Data
hfig_bcr_theory1=figure('Units','centimeters',...
    'PaperUnits','centimeters',...
    'PaperPositionMode','manual',...
    'PaperSize',papersize,...
    'PaperPosition',paperposition);
hold on;

namearray={'LineWidth'};
valarray={linewidth};

%%% approx theory - at an angle
jet_theta=0.61;      % jet angle used
xth_jet=linspace(0.5,3);
yth_jet=polyval([cos(jet_theta),0],xth_jet);

% plot
figure(hfig_bcr_theory1);
hold on;
plot(xth_jet,yth_jet,'k--','LineWidth',2);

% cc=distinguishable_colors(Nexp);
% cc_fringe=distinguishable_colors(Ndpeakplot);
% cc_fringe=gray(Ndpeakplot+2);
cc_fringe=viridis(Ndpeakplot);
p=[];
hexpconfig=[];
for ii=1:Nexp
    thisS=S_new(ii);
    
    ndpeak=(thisS.N_peak_max-1);
    if Ndpeakplot<ndpeak
        ndpeak=Ndpeakplot;
    end
    for jj=1:ndpeak
        hdata_theory_summ=ploterr(1e-6*thisS.X_collate(jj),1e-6*thisS.Y_collate(jj),...
            1e-6*thisS.X_collate_err_tot(jj),1e-6*thisS.Y_collate_err_tot(jj),...
            mm{ii},'hhxy',0);
        set(hdata_theory_summ(1),namearray,valarray,'Color',cc_fringe(jj,:),'MarkerFaceColor',cc_fringe(jj,:),'MarkerSize',mrk_size,'DisplayName',sprintf('%d',ii));
        set(hdata_theory_summ(2),namearray,valarray,'Color',cc_fringe(jj,:),'MarkerFaceColor',cc_fringe(jj,:),'DisplayName','');
        set(hdata_theory_summ(3),namearray,valarray,'Color',cc_fringe(jj,:),'MarkerFaceColor',cc_fringe(jj,:),'DisplayName','');
        %     p(ii)=hdata_theory_summ(1);
    end
    
    % legend for experimental configs
    displayname=sprintf('%0.2g, %0.2g',thisS.N0(1),thisS.eff_al);
    displayname=strrep(displayname,'+0','');
   	hexpconfig(ii)=plot(NaN,NaN,mm{ii},...
        'Color',cgray,'MarkerFaceColor',cgray,...
        'DisplayName',displayname);
    p(ii)=hexpconfig(ii);
end
set(gca,'Units','normalized',...
    'YTick',[0:0.5:3],...
    'XTick',[0:0.5:3],...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',fontsize,...
    'PlotBoxAspectRatio',[1,0.6,1]);
box on;
% axis square;
xlim([0.75,2.75]);
ylim([0.3,1.5]);

% legend
lgd=legend(p,'Location','northwest');
set(lgd,'Units','normalized',...
    'Position',[0.25 0.64 0.1 0.1],...
    'FontSize',8,...
    'Box','on');
title(lgd,'$N_{0},\eta_{\textrm{RF}}$');

% axis labels ~31/08/2017 - confusing
xlabel('$2 m / \hbar \cdot (v^2 - c^2)^{1/2}$ [$\mu$m$^{-1}$]');
ylabel('$2 \pi / \lambda_{\theta} $ [$\mu$m$^{-1}$]');

% save
if savefigs
    figname=sprintf('fig_BCR1_%s',datetimestr);
    print(hfig_bcr_theory1,fullfile(path_save,figname),'-dpdf');
end

%% BCR theory - TYPE 2
%%% Data
hfig_bcr_theory2=figure('Units','centimeters',...
    'PaperUnits','centimeters',...
    'PaperPositionMode','manual',...
    'PaperSize',papersize,...
    'PaperPosition',paperposition);
hold on;

%%% theory curve
xth_jet=linspace(0.5,Ndpeakplot+0.5);
yth_jet=polyval([0,1],xth_jet);

% plot
figure(hfig_bcr_theory2);
hold on;
plot(xth_jet,yth_jet,'k--','LineWidth',2);

% cc=distinguishable_colors(Nexp);
% cc_fringe=distinguishable_colors(Ndpeakplot);
% cc_fringe=gray(Ndpeakplot+2);
cc_fringe=viridis(Ndpeakplot);
p=[];
hexpconfig=[];
for ii=1:Nexp
    thisS=S_new(ii);
    
    ndpeak=(thisS.N_peak_max-1);
    if Ndpeakplot<ndpeak
        ndpeak=Ndpeakplot;
    end
    for jj=1:ndpeak
        hdata_theory_summ=ploterr(jj,thisS.Y_collate(jj)./(cos(jet_theta)*thisS.X_collate(jj)),...
            [],thisS.Y_collate_err_tot(jj)./(cos(jet_theta)*thisS.X_collate(jj)),...
            mm{ii},'hhxy',0);
        set(hdata_theory_summ(1),namearray,valarray,'Color',cc_fringe(jj,:),'MarkerFaceColor',cc_fringe(jj,:),'MarkerSize',mrk_size,'DisplayName',sprintf('%d',ii));
        set(hdata_theory_summ(2),namearray,valarray,'Color',cc_fringe(jj,:),'MarkerFaceColor',cc_fringe(jj,:),'DisplayName','');
%         set(hdata_theory_summ(3),namearray,valarray,'Color',cc_fringe(jj,:),'MarkerFaceColor',cc_fringe(jj,:),'DisplayName','');
        %     p(ii)=hdata_theory_summ(1);
    end
    
    % legend for experimental configs
    displayname=sprintf('%0.2g, %0.2g',thisS.N0(1),thisS.eff_al);
    displayname=strrep(displayname,'+0','');
   	hexpconfig(ii)=plot(NaN,NaN,mm{ii},...
        'Color',cgray,'MarkerFaceColor',cgray,...
        'DisplayName',displayname);
    p(ii)=hexpconfig(ii);
end
set(gca,'Units','normalized',...
    'XTick',[0:1:Ndpeakplot],...
    'YTick',[0:0.5:2],...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',fontsize,...
    'PlotBoxAspectRatio',[1,0.6,1]);
box on;
% axis square;
xlim([0.5,Ndpeakplot+0.5]);
ylim([0.5,2]);

% legend
lgd=legend(p,'Location','northeast');
set(lgd,'Units','normalized',...
    'Position',[0.7 0.63 0.1 0.1],...
    'FontSize',8,...
    'Box','on');

title(lgd,'$N_{0},\eta_{\textrm{RF}}$');

xlabel('Fringe number');
ylabel('$k_{\textrm{meas}}/k_{\textrm{theory}}$');

% save
if savefigs
    figname=sprintf('fig_BCR2_%s',datetimestr);
    print(hfig_bcr_theory2,fullfile(path_save,figname),'-dpdf');
end

%% a representative interference pattern - 1D density profile along jet
% Run 2 - the second dataset in mlist - for plotting since already used in
% paper, most fringes too

pal_id_plot=6;      % the PAL number to plot

xlim_cfg=[0,0.7];
ytick_cfg=-1:1:1;
xtick_cfg=xlim_cfg(1):0.2:xlim_cfg(2);

% figure params
mrk_size=7;
fontsize=11;
linewidth=2;

papersize=[10,3];
paperposition=[0,0,papersize];

%%% Data
hfig_nmod=figure('Units','centimeters',...
    'PaperUnits','centimeters',...
    'PaperPositionMode','manual',...
    'PaperSize',papersize,...
    'PaperPosition',paperposition);
hold on;

loadvars={'dn1d','d_1d','ppeak','pal_R'};
ss=load(fullfile(path_data,mlist{2}),loadvars{:});

% normalise the density profiles:
xx=ss.d_1d./ss.pal_R;      % normalised distance vector along jet
ppeak=ss.ppeak;     % hasn't been normed yet
dn1d=cellfun(@(x) x/max(x),ss.dn1d,'UniformOutput',false);  % diff fringe density

hold on;
for ii=1:min(length(ppeak{pal_id_plot}),Ndpeakplot+1)
%     rectangle
    line(1e-3*ppeak{pal_id_plot}(ii)/ss.pal_R(pal_id_plot)*[1,1],1.1*[-1,1],...
        'Color',cgray,'LineStyle','--','LineWidth',1.3);
end
plot(xx(pal_id_plot,:),dn1d{pal_id_plot},'Color','k','LineWidth',linewidth);
box on;
hold on;
% axis tight;
xlim(xlim_cfg);
ylim([-1.1,1.1]);

xlabel('$r/R_{\textrm{AL}}$');
ylabel('$\widetilde{n}$ (arb. u.)');

set(gca,'Units','normalized',...
    'XTick',xtick_cfg,...
    'YTick',ytick_cfg,...
    'XAxisLocation','bottom',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',fontsize);

% save
if savefigs
    figname=sprintf('nmod_%s',datetimestr);
    print(hfig_nmod,fullfile(path_save,figname),'-dpdf');
end