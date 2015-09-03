function [dc_offset,maxcorr] = align_dc_ac_sigs_initial_v2(dc_sig,dc_times,ac_sig)

Fs = 1/(dc_times(2)-dc_times(1));
niqf = Fs/2;
% [b,a] = butter(2,[2/niqf 80/niqf]);
[b,a] = butter(2,2/niqf,'high');
dc_sig = filtfilt(b,a,dc_sig);

%%
Fs = 2016;
% dsf = 8;
dsf = 1;
Fsd = Fs/dsf;
niqf = Fs/2;
% [b,a] = butter(2,[2/niqf 80/niqf]);
[b,a] = butter(2,2/niqf,'high');
ac_f = filtfilt(b,a,ac_sig);
ac_fd = zscore(downsample(ac_f,dsf));
ac_t = (1:length(ac_fd))/Fsd;

%%
dc_sig_r = interp1(dc_times,dc_sig,ac_t);
dc_sig_r(isnan(dc_sig_r)) = 0;

maxlag = round(Fsd*700);
[x,l] = xcov(ac_fd,dc_sig_r,maxlag,'coef');
[maxcorr,maxloc] = max(x);
maxloc = l(maxloc);
dc_offset = maxloc/Fsd;

%%



