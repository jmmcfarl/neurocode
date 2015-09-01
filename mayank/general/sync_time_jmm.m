function [] = sync_time_jmm()

% synchronize the times for two time vectors recorded at different sampling
% frequencies to extract the IDs (id1) of the higher sampling rate data that
% correspond to the lower sampling rate data (id2). In the following
% length(t1) > length(t2). This is an easy algorithm that assumes that no
% data ponts are dropped during the recording and time runs continuously
% and smoothly.

if exist('all_eeg_data2.mat')
    load ./all_eeg_data2
% elseif exist('part2_eeg_data.mat')
%     load part2_eeg_data
else
    load ./all_eeg_data
end

clear CSC2* CSC3* CSC4* CSC5* CSC6* CSC7* *Channel*



%get sample period in micros for each channel
dt1 = median(diff(CSC1_TimeStamps))/512;
dt2 = median(diff(CSC8_TimeStamps))/512;

if ~isempty(find(diff(CSC1_TimeStamps) > dt1*512))
    disp('error dropped packet')
end
%     endpt = find(diff(CSC1_TimeStamps) > dt1*512)-1;
%     endtime = CSC1_TimeStamps(endpt);
%     endpt_8 = find(CSC8_TimeStamps > endtime,1,'first');
%     endtime_8 = CSC8_TimeStamps(endpt_8);
% else
    endpt = length(CSC1_TimeStamps);
    endpt_8 = length(CSC8_TimeStamps);
% end
if ~isempty(find(diff(CSC8_TimeStamps) > dt2*512))
    disp('error dropped packet csc8')
end

%make sure there is no mismatch between Timestamp and Sample vectors
CSC1_Samples(:,endpt+1:end) = [];
CSC8_Samples(:,endpt_8+1:end) = [];

size_CSC1 = size(CSC1_Samples);
size_CSC8 = size(CSC8_Samples);

%initialize time vectors
t1 = zeros(size_CSC1);
t2 = zeros(size_CSC8);

%interpolate time stamps
for k = 1:512
    t1(k,:) = CSC1_TimeStamps(1:endpt) + round((k-1)*dt1);
    t2(k,:) = CSC8_TimeStamps(1:endpt_8) + round((k-1)*dt2);
end

%reshape time stamps and find last sample number
t1 = reshape(t1,512*size_CSC1(2),1); 
t1endid = length(t1);
t2 = reshape(t2,512*size_CSC8(2),1); 
t2endid = length(t2);

%find closest CSC2 timestamp to the first CSC1 timestamp
[tbegdif, m] = min(abs(t2 - t1(1)));
%find the closest CSC1 timestamp to that CSC2 timestamp 
[tbegdif, k] = min(abs(t1 - t2(m)));
%if the match isn't perfect try subsequent CSC1 timestamps until you find a
%match
while tbegdif > 0
    k = k +1;
    [tbegdif, m] = min(abs(t2 - t1(k)));
end

%set initial sample IDs for alignment point
synctbegin = t1(k);
synct1begid = k;
synct2begid = m;

%repeat process for finding synchronous ending points
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

