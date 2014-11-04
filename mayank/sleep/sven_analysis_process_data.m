clear all

% addpath('~/James_scripts/Nlx2Mat_rel041811/')

% close all
%
% cd /Users/James/Data/Sven/2012_06_19/2012-6-19_12-21-9
%
% cd /Users/James/Data/Sven/2012_06-22/2012-6-22_14-31-1
% cd /Users/James/Data/Sven/2012_06-22/2012-6-22_18-52-27
%
% cd /Users/James/Data/Sven/2012-6-23_13-45-41
% cd /Users/James/Data/Sven/2012-6-23_14-44-12
%
% cd /Users/James/Data/Sven/2012_06_27/2012-6-27_14-35-25
% cd /Users/James/Data/Sven/2012_06_27/2012-6-27_15-36-56
% cd /Users/James/Data/Sven/2012_06_27/2012-6-27_16-42-48
% cd /Users/James/Data/Sven/2012_06_27/2012-6-27_17-42-24

% cd ~/Data/Sven/2012-6-30_16-6-48/
% cd ~/Data/Sven/2012-6-30_17-29-16/
% cd ~/Data/Sven/2012-6-30_18-44-18/

% cd ~/Data/Sven/2012-7-2_15-5-21/
% cd ~/Data/Sven/2012-7-2_16-19-30/
% cd ~/Data/Sven/2012-7-2_17-17-14/

% cd C:\WC_Germany\sleep\2012_07_03\2012-7-3_13-17-54
% cd C:\WC_Germany\sleep\2012_07_03\2012-7-3_14-17-46
% cd C:\WC_Germany\sleep\2012_07_03\2012-7-3_15-25-44

% cd C:\WC_Germany\sleep\Sleep_07-06\2012-7-6_13-47-20
% cd C:\WC_Germany\sleep\Sleep_07-06\2012-7-6_15-1-5
% cd C:\WC_Germany\sleep\Sleep_07-06\2012-7-6_16-12-14

% cd C:\WC_Germany\sleep\2012-7-5_2MiPs\2012-7-5_14-7-37
% cd C:\WC_Germany\sleep\2012-7-5_2MiPs\2012-7-5_13-16-58

% cd C:\WC_Germany\sleep\2012-7-8_2\2012-7-8_13-9-45
% cd C:\WC_Germany\sleep\2012-7-8_2\2012-7-8_14-16-26

cd C:\WC_Germany\sleep\2012-7-5_Patch\2012-7-5_17-11-29
%% LOAD RAW DATA
FieldSelectionFlags = [1 0 1 0 1];
ExtractMode = 1;
HeaderExtractionFlag = 1;

for i = 1:8
    fprintf('Loading CSC%d\n',i)
    Filename = sprintf('CSC%d.Ncs',i);
    [CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
        Nlx2MatCSC(Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
    CSC_Timestamps = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
    Fs(i) = median(CSC_SampleFrequencies);
    csc(i,:) = CSC_Samples(:)*str2num(Header{15}(13:end));
    clear CSC_Samples
end
if length(unique(Fs)) > 1
    error('Unequal sample rates')
end
Fs = unique(Fs);

%% DOWN-SAMPLE
dsf = 8;
Fsd = Fs/dsf;
for i = 1:8
    cscd(:,i) = decimate(csc(i,:),dsf);
end
t = (1:size(cscd,1))/Fsd;

%% LOAD MUA
FieldSelectionFlags = [1 0 0 1 0];
HeaderExtractionFlag = 1;
ExtractMode = 1;
Filename = 'Sc1.ntt';
[Timestamps, Samples, Header] = ...
    Nlx2MatSpike( Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
conv_factor = str2num(Header{15}(13:end));
peakvals = bsxfun(@times,Samples(1:4,:),conv_factor');
clear Samples
[ov_peak,peak_ch] = max(peakvals);
for i = 1:4
    cur_mua = find(peak_ch == i);
    mua{i} = Timestamps(cur_mua);
end
Filename = 'Sc2.ntt';
[Timestamps, Samples, Header] = ...
    Nlx2MatSpike( Filename, FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
conv_factor = str2num(Header{15}(13:end));
peakvals = bsxfun(@times,Samples(1:4,:),conv_factor');
clear Samples
[ov_peak,peak_ch] = max(peakvals);
for i = 1:4
    cur_mua = find(peak_ch == i);
    mua{i+4} = Timestamps(cur_mua);
end

%% BIN MUA
down_ts = downsample(CSC_Timestamps,dsf);
for ch = 1:8;
    cur_mua = mua{ch};
    cur_mua(cur_mua > down_ts(end) | cur_mua < down_ts(1)) = [];
    binned_spks(:,ch) = hist(cur_mua,down_ts)';
end
binned_spks(end,:) = mean(binned_spks);

%%
temp = pwd;
sl = find(temp == '\',1,'last');
dname = temp(sl+1:end);
save(dname,'binned_spks','dsf', 'Fsd', 't', 'cscd');