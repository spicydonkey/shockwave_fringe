%% Combine fringe spacing data from different experiments
% load mat files generated by main_shockfringe.m from various experiments
%

%% CONFIG
Ndpeakplot=2;    % number of fringe spacings to plot

savefigs=0;
path_save='C:\Users\HE BEC\Desktop\shock_summary';

path_data='C:\Users\HE BEC\Documents\lab\shockwave\summary\data\ver2';
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
    % load vars:
        % MAX_PEAK_N, Nal, AL_N_SD PEAK_DIFF_ALL, PEAK_DIFF_SD, N0
    varsummary={'MAX_PEAK_N', 'PEAK_DIFF_ALL', 'PEAK_DIFF_SD', 'Nal', 'Nal_SE', 'N0', 'N0_SE', 'Nfit'};
    S{ii}=load(fullfile(path_data,mlist{ii}),varsummary{:});
    
    % TODO below data analysis needs to be in main and will need to update
    % summary
%     %%% number in BEC
%     % estimate error from "fitobject" - so poorly named - Nfit
%     Nfit=S{ii}.fitobject;
%     Nfitrelerr=Nfit.Coefficients.SE./Nfit.Coefficients.Estimate;   % N0, eta_RF fit err (rel)
%     N0_SE_rel=sqrt(sum(Nfitrelerr.^2)).*(sqrt(1:length(S{ii}.N0))');     % rel error for N0
%     N0_SE=S{ii}.N0.*N0_SE_rel;       % evaluate SE
%     S{ii}.N0_SE=N0_SE;          % store into data structure
    
%     %%% number in AL
%     Nal_SE_rel=N0_SE_rel';          % scales identically with r=eta_RF + N0(0) from formula
%     Nal_SE=S{ii}.Nal.*Nal_SE_rel;	% evaluate SE from fit uncertainties
%     % add in quadrature with detected number fluctuation SE
%     Nal_SE=sqrt(S{ii}.AL_N_SD'.^2+Nal_SE.^2);
%     S{ii}.Nal_SE=Nal_SE;        % store into data structure
end
% tidy structure
S=cell2mat(S);      % S is now struct array

%% initialise
Npeakmax=max([S.MAX_PEAK_N]);

% plot
cc=distinguishable_colors(Npeakmax);
mm={'o','s','^','x','.','*','+','d','v','>','<','p','h'};


%% plot - fringe spacing vs Nal
hfig_dpeak_vs_Npal=figure();

namearray={'LineWidth','MarkerFaceColor'};      % error bar graphics properties
valarray={linewidth,'w'};                 % 90 deg (normal) data

for ii=1:Nexp
    thisS=S(ii);
    ndpeak=(thisS.MAX_PEAK_N-1);
    if Ndpeakplot<ndpeak
        ndpeak=Ndpeakplot;
    end
%     p=zeros(1,ndpeak);
    for jj=1:ndpeak
        hold on;
%         hdata_pal_n=ploterr(thisS.Nal,thisS.PEAK_DIFF_ALL(:,jj),thisS.AL_N_SD,thisS.PEAK_DIFF_SD(:,jj),mm{ii},'hhxy',0);
hdata_pal_n=ploterr(thisS.Nal,thisS.PEAK_DIFF_ALL(:,jj),thisS.Nal_SE,thisS.PEAK_DIFF_SD(:,jj),mm{ii},'hhxy',0);
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
    ndpeak=(thisS.MAX_PEAK_N-1);
    if Ndpeakplot<ndpeak
        ndpeak=Ndpeakplot;
    end
%     p=zeros(1,ndpeak);
    for jj=1:ndpeak
        hold on;
        hdata_pal_n=ploterr(thisS.N0,thisS.PEAK_DIFF_ALL(:,jj),thisS.N0_SE,thisS.PEAK_DIFF_SD(:,jj),mm{ii},'hhxy',0);
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
% try:
%%% density (avg)
% TODO - this changes in AL - currently assumes uniform density
% nal ~ N_al ^ alpha * N_BEC ^ beta
alpha=1;
beta=-3/5;

%%% speed of sound (avg)
% TODO - this changes in AL
% c ~ nal ^ gamma
gamma=1/2;

%%% collision speed (avg)
% TODO - calculated for fringes!
% v ~ CONST (FOR NOW)

%%% mach number (avg)
% mach ~ v / c

for ii=1:Nexp
    % density
    S(ii).nal=(S(ii).Nal'.^alpha).*(S(ii).N0.^beta);
    S(ii).c=S(ii).nal.^gamma;
    S(ii).v=1;
%     S(ii).mach=(S(ii).N0.^(3/10)).*(S(ii).Nal'.^(-1/2));
    S(ii).mach=S(ii).v./S(ii).c;
end

hfig_dpeak_scaled=figure();

namearray={'LineWidth','MarkerFaceColor'};      % error bar graphics properties
valarray={linewidth,'w'};                 % 90 deg (normal) data

for ii=1:Nexp
    thisS=S(ii);
    ndpeak=(thisS.MAX_PEAK_N-1);
    if Ndpeakplot<ndpeak
        ndpeak=Ndpeakplot;
    end
    %     p=zeros(1,ndpeak);
    for jj=1:ndpeak
        hold on;
        hdata_pal_n=ploterr(1./thisS.mach,thisS.PEAK_DIFF_ALL(:,jj).*S(ii).c,[],[],mm{ii},'hhxy',0);
%                 hdata_pal_n=ploterr(thisS.nal,thisS.PEAK_DIFF_ALL(:,jj).*S(ii).c,[],[],mm{ii},'hhxy',0);
%         hdata_pal_n=ploterr(thisS.nal,thisS.PEAK_DIFF_ALL(:,jj),[],[],mm{ii},'hhxy',0);
        set(hdata_pal_n(1),namearray,valarray,'Color',cc(jj,:),'MarkerSize',6,'DisplayName',sprintf('%d',jj));
%         set(hdata_pal_n(2),namearray,valarray,'Color',cc(jj,:),'DisplayName','');
%         set(hdata_pal_n(3),namearray,valarray,'Color',cc(jj,:),'DisplayName','');
%         p(jj)=hdata_pal_n(1);
    end
end
box on;
% lgd=legend(p);
% title(lgd,'Fringe spacing');
% xlabel('$\tilde{\rho}_{AL}$');
% ylabel('Reduced fringe spacing');
xlabel('inverse Mach number');
ylabel('$\lambda \cdot c$');
% 
% %%% save fig
% if savefigs>0
%     figname=sprintf('dpeak_vs_N0_%s',datetimestr);
%     saveas(hfig_dpeak_vs_N0,fullfile(path_save,[figname,'.png']));
%     saveas(hfig_dpeak_vs_N0,fullfile(path_save,[figname,'.fig']));
% end