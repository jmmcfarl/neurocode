function [min_amp,min_loc,max_amp,max_loc,sm_acorr] = get_acorr_peaks(acorr,lags,smth,disp)

if nargin < 4
    disp = 0;
end
zero_lag = find(lags == 0);
sm_acorr = smooth(acorr,smth);
smth_acorr = sm_acorr(zero_lag:end);

% [min_amp,min_loc] = min(smth_acorr);
[min_amp,min_loc] = findpeaks(-smth_acorr,'npeaks',1);
[max_amp,max_loc] = findpeaks(smth_acorr(min_loc:end),'npeaks',1);

bp = zero_lag;
ep = zero_lag + min_loc + max_loc-1;
[min_amp,min_loc] = min(acorr(bp:ep));
min_loc = min_loc + bp-1;
bp = min_loc;
[max_amp,max_loc] = max(acorr(bp:end));
max_loc = max_loc + bp-1;
min_amp = acorr(min_loc);
max_amp = acorr(max_loc);

if disp==1
    figure
    plot(lags,acorr);
    hold on
    plot(lags,sm_acorr,'r','linewidth',2)
    plot(lags(min_loc),acorr(min_loc),'ko')
    plot(lags(max_loc),acorr(max_loc),'go')
    pause
    close 
end