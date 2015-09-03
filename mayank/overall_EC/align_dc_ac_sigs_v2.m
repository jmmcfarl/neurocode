function [sweep_offsets,sweep_maxcorr,sweep_starts,sweep_stops] = align_dc_ac_sigs_v2(dc_sig,dc_times,ac_sig)

dc_Fs = 1/(dc_times(2)-dc_times(1));
niqf = dc_Fs/2;
% [b,a] = butter(2,[2/niqf 20/niqf]);
[b_dc,a_dc] = butter(2,2/niqf,'high');

%if sweeps are not relative time, correct this
if max(dc_times) > 50
    jump_inds = find(diff(dc_times) > 1e-2) + 1;
    for i = 1:length(jump_inds)-1
        dc_times(jump_inds(i):jump_inds(i+1)-1) = dc_times(jump_inds(i):jump_inds(i+1)-1) - dc_times(jump_inds(i));
    end
    dc_times(jump_inds(length(jump_inds)):end) = dc_times(jump_inds(length(jump_inds)):end) - dc_times(jump_inds(length(jump_inds)));
    sweep_starts = [1 jump_inds];
    sweep_stops = [jump_inds-1 length(dc_times)];
else    
    sweep_starts = find(dc_times==min(dc_times));
    sweep_stops = find(dc_times==max(dc_times));
    if length(sweep_starts) > length(sweep_stops)
        sweep_starts(end) = [];
    end
end
%%
Fs = 2016;
niqf = Fs/2;
% [b,a] = butter(2,[2/niqf 20/niqf]);
[b,a] = butter(2,2/niqf,'high');
ac_f = filtfilt(b,a,ac_sig);
% ac_r = resample(ac_f,Fsd,2016);
ac_r = zscore(ac_f);

%%
cur_loc = 1;
maxlag = round(30*Fs);
sweep_offsets = nan(1,length(sweep_starts));
sweep_maxcorr = nan(1,length(sweep_starts));
for t = 1:length(sweep_starts)
    fprintf('Sweep %d of %d\n',t,length(sweep_starts));
    cur_dc_samples = sweep_starts(t):sweep_stops(t);
    cur_dc_data = dc_sig(cur_dc_samples);
    cur_dc_data = filtfilt(b_dc,a_dc,cur_dc_data);
    interp_dc_taxis = dc_times(cur_dc_samples(1)):1/Fs:dc_times(cur_dc_samples(end));
    cur_dc_interp = interp1(dc_times(cur_dc_samples),cur_dc_data,interp_dc_taxis);
    dc_xcov_sig = zeros(size(ac_r));
    dc_xcov_sig(cur_loc:cur_loc+length(cur_dc_interp)-1) = cur_dc_interp;
    if length(dc_xcov_sig) > length(ac_r)
        dc_xcov_sig(length(ac_r)+1:end) = [];
    end
    [x,l] = xcov(ac_r,dc_xcov_sig,round(30*Fs),'coeff');
    [sweep_maxcorr(t),maxlag] = max(x);
    sweep_offsets(t) = (cur_loc+l(maxlag))/Fs;
    cur_loc = round(sweep_offsets(t) + range(interp_dc_taxis))*Fs;
end

