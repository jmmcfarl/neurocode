function [] = wcv_to_spktim()

% obtain spike times from WCV

if exist('all_eeg_data2.mat','file')
    load ./all_eeg_data2
% elseif exist('part2_eeg_data.mat')
%     load part2_eeg_data
else
    load ./all_eeg_data
end

[len,widwcv] = size(CSC1_Samples);

% Get the data.
wcv = reshape(CSC1_Samples,len*widwcv,1);
f2 = CSC8_SampleFrequencies(1);
f = CSC1_SampleFrequencies(1);

clear CSC1* CSC2* CSC3* CSC4* CSC5* CSC6* CSC7* CSC8_Ch* CSC8_Num*
wcv = -wcv; % this sign flip ensures that spikes go upwards.

load sync_times synct synct1id

synct = synct*10^-6; %convert to sec
lwcv = wcv(synct1id);

lof = 100;  %default 300
if f > 16000
    hif = 8000;
else
    hif = f/2.1;
end
nyqf = f/2;
lof = lof/nyqf; hif = hif/nyqf;
[b,a] = butter(2, [lof hif]);

wcvlen = length(wcv);

spkthresh = 10; % %default 15

wcvf = filtfilt(b,a,wcv);
dwcv = ([diff(wcvf') 0])';
alen = length(dwcv);

%********************************
% if you want to do a 4-point smoothing of the derivative
dwcv(1:alen-3) = 0.25*( dwcv(1:alen-3) + dwcv(2:alen-2) + dwcv(3:alen-1)+ dwcv(4:alen) );
%********************************

dwcv = dwcv - mean(dwcv);
dwcv = dwcv/std(dwcv);
maxspk = max(dwcv);
%     spkthresh = max([5 maxspk/2]);
spkid = find( (dwcv(1:alen-1) < spkthresh) & (dwcv(2:alen) >= spkthresh) );
spkid = spkid+1;
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

spkamp = lwcv(spkid);
save spike_time_jmm spkid spkamp
% length(spkid)
% pause;
