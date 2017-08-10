%% Fringe characterisation by peak spacing through line profile

%%% configure
% load line profile
d=zrot_c(id_y);     % distance along line

% peak detection algorithm
SlopeThreshold=1e-4;
AmpThreshold=0.15;
smoothwidth=0;          % 0 to no smoothing in peak finding
peakgroup=5;
smoothtype=1;   % smoothing is off

% preallocate
peak_list=cell(pal_nseq,1);         % full peak summary output from findpeaksG
ppeak=cell(pal_nseq,1);         % peak locations for each PAL
dppeak=cell(pal_nseq,1);        % peak spacing

hfig_peaks=figure();
p=zeros(10,1);
cc=distinguishable_colors(10);

% find peaks!
for ii=1:pal_nseq
    n=nn1d{ii};         % 1D density profile along line
    n=n/max(n);         % normalise density; max=1
    
    this_peak_list=findpeaksG(d,n,SlopeThreshold,AmpThreshold,smoothwidth,peakgroup,smoothtype);
    peak_list{ii}=this_peak_list;
    
    figure(hfig_peaks);
    hold on;
    p(ii)=plot(d,n,'color',cc(ii,:),'LineWidth',1.5,'DisplayName',sprintf('%d',ii));
    scatter(this_peak_list(:,2),this_peak_list(:,3),'o','MarkerEdgeColor',cc(ii,:));
    text(this_peak_list(:,2),this_peak_list(:,3),num2str(this_peak_list(:,1)),...
        'Color',cc(ii,:),'FontSize',15);  % annotate peaks
    
    ppeak{ii}=this_peak_list(:,2);
    dppeak{ii}=diff(ppeak{ii});
end
% annotate figures
figure(hfig_peaks);
box on;
xlabel('distance [m]');
ylabel('density [arb]');
legend(p);

%%% plot peak spacings
% peak narrowing
hfig_dppeaks_vs_n=figure();
hold on;
for ii=1:pal_nseq
    plot(dppeak{ii},'o','color',cc(ii,:));
end
box on;
xlabel('$n$');
ylabel('$n$-th peak spacing [m]');

% peak spacing dependence on PAL
hfig_dppeaks_vs_pal=figure();
max_peak_n=max(cellfun(@(x)size(x,1),ppeak));       % max number of peaks found in PALs
for ipeak=1:max_peak_n-1
    subplot(max_peak_n-1,1,ipeak);
    xlim([1,pal_nseq]);
    hold on;
    box on;
    title(sprintf('peak spacing: %d',ipeak));
    % improve by padding with NaN, etc.
    for ii=1:pal_nseq
        if length(dppeak{ii})>=ipeak
            scatter(ii,dppeak{ii}(ipeak),'MarkerEdgeColor',cc(ii,:));
        end
    end
end