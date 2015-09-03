function [spkAmp] = getSpkAmp(spkid,wcv,Fs)

%high pass filter out low frequency signal

nyqf = Fs/2;
lif = 60/nyqf;

[b,a] = butter(4,lif,'high');
hiwcv = filtfilt(b,a,wcv);


for i=  1:length(spkid)
    
    
    begPt = spkid(i)-1;
    endPt = begPt + 100;
    spkAmp(i) = max(hiwcv(begPt:endPt));
    
end