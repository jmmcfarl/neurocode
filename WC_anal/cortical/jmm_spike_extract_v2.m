
% obtain spike times from WCV

load all_eeg_data CSC1_Samples CSC8_SampleFrequencies CSC1_SampleFrequencies

[len,widwcv] = size(CSC1_Samples);

% Get the data.
wcv = reshape(CSC1_Samples,len*widwcv,1);
f = CSC1_SampleFrequencies(1);
f2 = CSC8_SampleFrequencies(1);
clear CSC1* CSC8*
% clear CSC1* CSC2* CSC3* CSC4* CSC5* CSC6* CSC7* CSC8_Ch* CSC8_Num*
wcv = -wcv; % this sign flip ensures that spikes go upwards.

% wcv = detrend(wcv); % detrend wcv.
load sync_times synct synct1id

synct = synct*10^-6;
lwcv = wcv(synct1id);
lof = 150;
lof2 = 50;

hif = 7999;

niqf = f/2;
[b,a] = butter(2,[lof/niqf hif/niqf]);
[b2,a2] = butter(2,[lof2/niqf hif/niqf]);
wcv_f = filtfilt(b,a,wcv);
wcv_f2 = filtfilt(b2,a2,wcv);

wcv_diff = ([0;diff(wcv_f)]);
wcvlen = length(wcv);
spkid = [];

spkthresh = 15; % %default 15

% %what does this step do??\
% aspk(1:alen-3) = 0.25*( aspk(1:alen-3) + aspk(2:alen-2) + aspk(3:alen-1)+ aspk(4:alen) );

wcv_diff = wcv_diff - mean(wcv_diff);
wcv_diff = wcv_diff/std(wcv_diff);
aspkid = find( (wcv_diff(1:wcvlen-1) < spkthresh) & (wcv_diff(2:wcvlen) >= spkthresh) );

spkid = [spkid aspkid'];

width_window = round(f*0.003);
window_size = round(f*0.001);
% for i = 1:length(spkid)
%     [spkamp(i), rel_peak] = max(wcv(spkid(i)-window_size:spkid(i)+window_size));
%     spkid(i) = spkid(i) + rel_peak - window_size;
%     if spkid(i) > width_window & spkid(i) < wcvlen-width_window
%            curseg = wcv_f2(spkid(i)-window_size:spkid(i)+window_size);
%            [spkpeak,pkloc] = max(curseg)
%            spkbeg = find(curseg(1:pkloc) < spkpeak*0.5,1,'last');
%            spkend = pkloc+find(curseg(pkloc:end) > spkpeak*0.5,1,'last');
%            if ~isempty(spkbeg) & ~isempty(spkend)
%                 spkdur(i) = spkend - spkbeg;
%            else
%                spkdur(i) = nan;
%            end
%     else
%            spkdur(i) = nan;
%     end
% end

% zwcv = zscore(wcv);

% spkdur = spkdur/f;
% 
% spkamp(zwcv(spkid) < 2) = [];
% spkdur(zwcv(spkid) < 2) = [];
% spkid(zwcv(spkid) < 1) = [];

dspkid = [1e6 diff(spkid)];
spkid(dspkid < f*.0003) = [];
% spkamp(dspkid < f*.001) = [];
% spkdur(dspkid < f*.001) = [];

wcv = wcv(synct1id);

spktim = zeros([wcvlen,1]);
for k = 1:ceil(f/f2)
    spktim(spkid+k-1,1) = 1;
end
spktim = spktim(synct1id);
a = find(diff(spktim) == 0);
spktim(a) = 0;
spkid = find(spktim);
% plot(synct,lwcv/max(lwcv),'r',synct,wcv/max(wcv),'c',synct(spkid),0.4*spktim(spkid),'b*');shg; 
% plot(synct,wcv/max(wcv),'r',synct(spkid),.4*spktim(spkid),'b*');

% save spike_time_jmm spkid spkamp spkdur
save spike_time_jmm spkid 

% length(spkid)
% pause;
