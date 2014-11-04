clear all
close all
cd C:\WC_Germany\sven_thomas_combined\
addpath('C:\WC_Germany\persistent_2010\')
%%
true_heka_fs = 1.99984e+004;

% cd C:\WC_Germany\sleep\2012-7-13#1\2012-7-13_13-31-49 %a number of good examples with good LFP, and OK MUA
% load ../2012_07_13_sleep_3.mat
% heka_data =  V2012_07_13_sleep_3__Read_only__Ch1.values;
% clear V2*

% cd C:\WC_Germany\sleep\2012-7-9_Sleep_Patch\2012-7-9_12-59-28
% load ../2012_07_09_sleep_1.mat
% heka_data = V2012_07_09_sleep_1_Ch1.values;
% clear V2*

cd C:\WC_Germany\sleep\2012-7-11_sleep_WC\2012-7-11_13-42-46 
load ../2012_07_11_sleep_2.mat  
heka_data = V2012_07_11_sleep_2_Ch1.values;
clear V2*

FieldSelectionFlags = [1 0 1 0 1];
ExtractMode = 1;
HeaderExtractionFlag = 1;
[CSC_Timestamps, CSC_SampleFrequencies,CSC_Samples, Header] = ...
    Nlx2MatCSC('CSC1.Ncs', FieldSelectionFlags,HeaderExtractionFlag, ExtractMode);
mp_t = interp1(1:512:numel(CSC_Samples),CSC_Timestamps,1:numel(CSC_Samples));
Fs = median(CSC_SampleFrequencies);
mp = CSC_Samples(:)*str2num(Header{15}(13:end));
mp_d = -downsample(mp,16);

dc_times = (1:length(heka_data))/true_heka_fs;

[dc_offset,dc_maxcorr] = align_dc_ac_sigs_initial_v2(heka_data,dc_times,mp_d);
dc_time = dc_times + dc_offset;
ac_time = (1:length(mp_d))/2016;
if min(diff(dc_time)) <= 0
    error('Alignment Problem!')
end
save aligned_heka heka_data dc_time ac_time dc_offset

