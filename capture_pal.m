function [pal]=capture_pal(txy,t1,dt,n)
% capture pulsed atom lasers
% [pal] = capture_pal(txy, dt, n, verbose)
%
% txy:  N x 1 cell array of txy counts
% t1:   time for first pulse [s]
% dt:   pulse period [s]
% n:    number of pulses to capture from first pulse
%
% pal:  n x 1 cell array of N x 1 cell array of txy counts in each pulse

% preallocate
pal=cell(n,1);

% capture pulses
for ii=1:n
    % build temporal window to capture PAL
    t_i=t1+(ii-1)*dt-0.5*dt;
    t_f=t_i+dt;     % window size is same as pulse period
    
    % capture in 1D t-window
    pal{ii}=cellfun(@(x) boxcull(x,{[t_i,t_f],[],[]}),txy,'UniformOutput',false);
end

end