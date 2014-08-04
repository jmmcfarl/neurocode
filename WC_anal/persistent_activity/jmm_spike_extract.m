% obtain spike times from WCV

load all_eeg_data

[len,widwcv] = size(CSC1_Samples);

% Get the data.
wcv = reshape(CSC1_Samples,len*widwcv,1);
f = CSC1_SampleFrequencies(1);
f2 = CSC8_SampleFrequencies(1);
clear CSC1* CSC2* CSC3* CSC4* CSC5* CSC6* CSC7* CSC8_Ch* CSC8_Num*
wcv = -wcv; % this sign flip ensures that spikes go upwards.
nseg = 1; %number of segments to divide the data into

% wcv = detrend(wcv); % detrend wcv.
load sync_times synct synct1id

synct = synct*10^-6;
lwcv = wcv(synct1id);
lof = 150;

hif = 7999;

niqf = f/2;
[b,a] = butter(2,[lof/niqf hif/niqf]);

wcv = filtfilt(b,a,wcv);
wcv = ([diff(wcv') 0])';
wcvlen = length(wcv);
seglen = floor(wcvlen/nseg);
spkid = [];

spkthresh = 15; % %default 15

for k = 1:nseg
    tmpseg = 1+(k-1)*seglen;
    tmpseglen = tmpseg:tmpseg+seglen;
    if k == nseg
        tmpseglen=tmpseg:wcvlen;
    end
    aspk = wcv(tmpseglen);
    alen = length(aspk);
    
    %what does this step do??\
    aspk(1:alen-3) = 0.25*( aspk(1:alen-3) + aspk(2:alen-2) + aspk(3:alen-1)+ aspk(4:alen) );
    
    aspk = aspk - mean(aspk);
    aspk = aspk/std(aspk);
    aspkid = find( (aspk(1:alen-1) < spkthresh) & (aspk(2:alen) >= spkthresh) );
    if length(aspkid)
        aspkid = aspkid + (k-1)*seglen;
        spkid = [spkid aspkid'];
    end
end

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
