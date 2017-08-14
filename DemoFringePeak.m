%% Fringe characterisation by peak spacing through line profile
% process captured 1D density profile with fringes to optimise or peak
% detection.
% do peak detection and output results

% %%% configure
% % background smoothing
% n_sm_bgd=15;
% n_sm_post=10;
% 
% % peak detection algorithm
% slopethreshold=1e-4;
% ampthreshold=0.1;
% smoothwidth=0;          % 0 to no smoothing in peak finding
% peakgroup=5;
% smoothtype=1;   % smoothing is off

%% main
% load line profile
d_1d=zz;     % distance along line

% preallocate
peak_list=cell(pal_nseq,1);     % full peak summary output from findpeaksG
ppeak=cell(pal_nseq,1);         % peak locations for each PAL [mm]
dppeak=cell(pal_nseq,1);        % peak spacing [mm]

if vgraph>0
    hfig_peaks=figure();
end

%% preprocess fringes for peak detection
% background subtraction - moving average
n1d_bgd=cellfun(@(x) smooth(x,n_sm_bgd),nn1d,'UniformOutput',false);

% background subtracted density profile
dn1d=cellfun(@(x,y) x-y,nn1d,n1d_bgd,'UniformOutput',false);        

% additional smoothing on diff density fringes
dn1d=cellfun(@(x) smooth(x,n_sm_post),dn1d,'UniformOutput',false);

%% Find peaks
p=zeros(10,1);  % array to store figure objects for selective legend
for ii=1:pal_nseq
    % normalise density fringe
    n=dn1d{ii};
    n=n/max(n);         % normalise to max=1
    
    % find peaks
    this_peak_list=findpeaksG(1e3*d_1d,n,slopethreshold,ampthreshold,smoothwidth,peakgroup,smoothtype);
    peak_list{ii}=this_peak_list;   % store the peaks in cell array
    
    %%% plot result
    if vgraph>0
        figure(hfig_peaks);
        hold on;
        % processed density fluctuations
        p(ii)=plot(1e3*d_1d,n,'color',cc(ii,:),'LineWidth',1.5,'DisplayName',sprintf('%d',ii));
        
        % detected peaks on profile
        scatter(this_peak_list(:,2),this_peak_list(:,3),'o','MarkerEdgeColor',cc(ii,:));
        text(this_peak_list(:,2),this_peak_list(:,3),num2str(this_peak_list(:,1)),...
            'Color',cc(ii,:),'FontSize',15);  % annotate peaks
    end
    
    % store peak positions summary
    ppeak{ii}=this_peak_list(:,2)';
    dppeak{ii}=diff(ppeak{ii});         % peak spacing [m]
end
if vgraph>0
    % annotate figures
    figure(hfig_peaks);
    box on;
    xlabel('distance [mm]');
    ylabel('density fluctuation from moving average [arb]');
    lgd=legend(p);
    title(lgd,'PAL');
end

% collate peak data in array
max_peak_n=max(cellfun(@(x)size(x,2),ppeak));       % max number of peaks found in PALs
peak_pos=NaN(pal_nseq,max_peak_n);                  % preallocate NaN array
for ii=1:pal_nseq
    this_peak_pos=ppeak{ii};
    peak_pos(ii,1:size(this_peak_pos,2))=this_peak_pos;
end
peak_diff=diff(peak_pos,1,2);   % evaluate peak spacing

cc2=distinguishable_colors(max_peak_n-1);   % for plotting per FRINGE #


% %%% plot peak spacings
% if vgraph>0
%     % peak narrowing
%     hfig_dppeaks_vs_n=figure();
%     hold on;
%     for ii=1:pal_nseq
%         plot(dppeak{ii},'o','color',cc(ii,:));
%     end
%     box on;
%     xlabel('$n$');
%     ylabel('$n$-th peak spacing [mm]');
%     
%     % peak spacing dependence on PAL
%     hfig_dppeaks_vs_pal=figure();
%     max_peak_n=max(cellfun(@(x)size(x,2),ppeak));       % max number of peaks found in PALs
%     for ipeak=1:max_peak_n-1
%         subplot(max_peak_n-1,1,ipeak);
%         xlim([1,pal_nseq]);
%         hold on;
%         box on;
%         title(sprintf('peak spacing: %d',ipeak));
%         for ii=1:pal_nseq
%             if length(dppeak{ii})>=ipeak
%                 scatter(ii,dppeak{ii}(ipeak),'MarkerEdgeColor',cc(ii,:));
%             end
%         end
%     end
% end

%% peak spacing dependance on PAL population
if vgraph>0
    linewidth=1.5;
    namearray={'LineWidth','MarkerFaceColor'};      % error bar graphics properties
    valarray={linewidth,'w'};                 % 90 deg (normal) data
    
    hfig_dpeak_vs_Npal=figure();
    p=zeros(1,(max_peak_n-1));      % array to store figure objects for selective legend
    for ii=1:(max_peak_n-1)
        hold on;
        hdata_pal_n=ploterr(Nal,peak_diff(:,ii),pal_n(:,2),[],'o','hhxy',0);
        set(hdata_pal_n(1),namearray,valarray,'Color',cc2(ii,:),'DisplayName',sprintf('%d',ii));
        set(hdata_pal_n(2),namearray,valarray,'Color',cc2(ii,:),'DisplayName','');
        p(ii)=hdata_pal_n(1);
    end
    box on;
    lgd=legend(p);
    title(lgd,'Fringe spacing');
    xlabel('$N_{AL}$');
    ylabel('Fringe spacing [mm]');
end

% %% peak spacing dependance on BEC population
% if vgraph>0
%     hfig_dpeak_vs_Nbec=figure();
%     for ii=1:(max_peak_n-1)
%         hold on;
%         plot(N0,peak_diff(:,ii),...
%             'o','Color',cc(ii,:),'DisplayName',sprintf('%d',ii));
%     end
%     box on;
%     lgd=legend('show');
%     title(lgd,'Fringe spacing');
%     xlabel('$N_{0}$');
%     ylabel('Fringe spacing [mm]');
% end