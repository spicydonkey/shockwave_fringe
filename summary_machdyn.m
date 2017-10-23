%% Summary plot for the hypersonic transition through PAL development

%% configs
path_data='C:\Users\HE BEC\Documents\lab\shockwave\summary\data\ver7';
data_regexp='run_*.mat';

datetimestr=datestr(datetime,'yyyymmdd_HHMMSS');    % timestamp when function called

varsummary={'r_cent','tof','g','nden_r','t0','c_const','pal_R','pal_nseq'};
% varsummary={'t_dyn','t_exit','mach_dyn','pal_nseq'};

% output control
idRunToPlot=[2,3];
plotNumPal=3;

% plot
alinewidth=2;
afontsize=15;


%% main
% parse directory
mlist=dir(fullfile(path_data,data_regexp));
mlist={mlist.name};
Nexp=length(mlist);     % number of datasets

S=cell(Nexp,1);     % S is cell array of structs
for ii=1:Nexp
    load(fullfile(path_data,mlist{ii}),varsummary{:});
    
    tS=[];
    tS.pal_nseq=pal_nseq;
    v_dyn=r_cent/tof;
    tS.t_dyn=2*v_dyn/g;
    n_dyn=nden_r.*(tof./t0).^3;
    c_dyn=c_const*sqrt(n_dyn);
    tS.mach_dyn=v_dyn./c_dyn;
    tS.t_exit=2*(pal_R./tof)/g;
    
    S{ii}=tS;
end
% tidy structure
S=cell2mat(S);      % S is now struct array

%% Plot
% Mach number dynamics
hfig_dyn_mach=figure();
pLineStyle={':','-','--','-.'};

hold on;
for ii=idRunToPlot
    tS=S(ii);

    palsToPlot=1;
    if plotNumPal~=1
        palsToPlot=round(linspace(1,tS.pal_nseq,plotNumPal));
        palsToPlot=flip(palsToPlot);    % plot later pulses first for better layering
    end
    cc=gray(plotNumPal+1);	% color set for plotting each AL
    cc=flip(cc);
    
    for jj=1:plotNumPal
        idx_pal=palsToPlot(jj);
        plot(1e3*tS.t_dyn(tS.t_dyn<tS.t_exit(idx_pal)),tS.mach_dyn(idx_pal,tS.t_dyn<tS.t_exit(idx_pal)),...
            'LineStyle',pLineStyle{ii},'Color',cc(jj+1,:),'LineWidth',alinewidth);
    end
end
hold off;

% annotate
ax=gca;
ax.FontSize=afontsize;
ax.YTick=0:40:200;
box on;
axis square;
axis tight;
xlabel('t (ms)');
ylabel('Mach number');