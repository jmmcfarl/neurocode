% synchronize the times for two time vectors recorded at different sampling
% frequencies to extract the IDs (id1) of the higher sampling rate data that
% correspond to the lower sampling rate data (id2). In the following
% length(t1) > length(t2). This is an easy algorithm that assumes that no
% data ponts are dropped during the recording and time runs continuously
% and smoothly.

dt1 = median(diff(CSC1_TimeStamps))/512;
dt2 = median(diff(CSC2_TimeStamps))/512;

[len1,wid1] = size(CSC1_Samples);
clear CSC1_Samples;
[len2,wid2] = size(CSC2_Samples);
clear CSC2_Samples;

t1 = zeros([len1,wid1]);
% t1 = uint64(t1);
t2 = zeros([len2,wid2]);
% t2 = uint64(t2);

[len,t1len] = size(t1);
[len,t2len]= size(t2);

for k = 1:len
    t1(k,:) = CSC1_TimeStamps + round((k-1)*dt1);
    t2(k,:) = CSC2_TimeStamps + round((k-1)*dt2);
end

t1 = reshape(t1,512*t1len,1); 
t1endid = length(t1);
t2 = reshape(t2,512*t2len,1); 
t2endid = length(t2);


[tbegdif, m] = min(abs(t2 - t1(1)));
[tbegdif, k] = min(abs(t1 - t2(m)));

while tbegdif > 0
    k = k +1;
    [tbegdif, m] = min(abs(t2 - t1(k)));
end
synctbegin = t1(k);
synct1begid = k;
synct2begid = m;

[tenddif, m] = min(abs(t2 - t1(t1endid)));
[tenddif, k] = min(abs(t1 - t2(m)));
while tenddif > 0 
    k = k - 1;
    [tenddif, m] = min(abs(t2 - t1(k)));
end
synctend = t1(k);
synct1endid = k;
synct2endid = m;

synct1len = synct1endid - synct1begid;
synct2len = synct2endid - synct2begid;

synct1byt2 = synct1len/synct2len;

if abs(synct1byt2 - round(synct1byt2)) > 0
    'unmatched sequences'
    synct1byt2
end

synct2id = synct2begid:synct2endid;
synct1id = synct1begid:synct1endid;
synct1id = downsample(synct1id,round(synct1byt2));
synct = t2(synct2id);

save sync_times synct2id synct1id synct
