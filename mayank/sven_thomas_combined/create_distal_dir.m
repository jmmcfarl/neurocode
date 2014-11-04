clear all
cd C:\WC_Germany\sven_thomas_combined

distal_dir{1} = 'C:\wc_data\2012_06_12\2012-6-12_16-19-16'; %good UDS, not very clear hippocampal LFP (MAYBE). MUA peak on ch. 4.
distal_dir{2} = 'C:\wc_data\2012_06_14\2012-6-14_13-40-31'; %good UDS, somewhat clear LFP peak channels 2 and 3. No MUA peak. 12s CIS (maybe too big amp)? 
distal_dir{3} = 'C:\wc_data\2012_06_15\2012-6-15_17-47-59'; %not really any usable UDS. 12s CIs (2 different amplitudes)
distal_dir{4} = 'C:\wc_data\2012_06_17\2012-6-17_14-34-10'; %sort of questionable ctx UDS. clear hpc LFP peak channels 2-3. MUA peak 2 (and 3)  12s CIs (big amp) 
distal_dir{5} = 'C:\wc_data\2012_06_17\2012-6-17_16-0-0'; %OK UDS, not great. clear hpc LFP peak channels 2-3. MUA peak 2 
distal_dir{6} = 'C:\wc_data\2012_07_04\2012-7-4_14-2-29'; %not usable UDS
distal_dir{7} = 'C:\wc_data\2012_07_05\2012-7-5_22-3-31'; %not usable UDS
distal_dir{8} = 'C:\wc_data\2012_07_06\2012-7-6_19-3-6'; %OK UDS
distal_dir{9} = 'C:\wc_data\2012_07_08\2012-7-8_19-3-43'; %OK recording, good hpc MUA peak on ch 3, questionable UDS, probably can't use
distal_dir{10} = 'C:\wc_data\2012_07_09\2012-7-9_18-35-16'; %not clear cortical UDS, good hpc MUA peak on ch 3, can't use
distal_dir{11} = 'C:\wc_data\2012_07_10\2012-7-10_21-16-55'; %not clear cortical UDS, huge delta though
distal_dir{12} = 'C:\wc_data\2012_07_11\2012-7-11_17-8-28';
distal_dir{13} = 'C:\wc_data\2012_07_11\2012-7-11_18-14-46';
distal_dir{14} = 'C:\wc_data\2012_07_12\2012-7-12_19-59-10';

distal_heka_dir{1} = 'C:\wc_data\2012_06_12\2012_06_12_1.mat';
distal_heka_dir{2} = 'C:\wc_data\2012_06_14\2012_06_14_2.mat';
distal_heka_dir{3} = '';
distal_heka_dir{4} = 'C:\wc_data\2012_06_17\2012_06_17_1.mat';
distal_heka_dir{5} = 'C:\wc_data\2012_06_17\2012_06_17_3.mat';
distal_heka_dir{6} = '';
distal_heka_dir{7} = '';
distal_heka_dir{8} = 'C:\wc_data\2012_07_06\2012_07_06_1_v2.mat';
distal_heka_dir{9} = '';
distal_heka_dir{10} = '';
distal_heka_dir{11} = '';
distal_heka_dir{12} = '';
distal_heka_dir{13} = '';
distal_heka_dir{14} = '';

for i = 1:length(distal_dir)
%     ctx_lfp(i) = 8;
    distal_heka_type{i} = 'cont';
end
ctx_lfp = [8 5 nan 8 8 nan nan 5 5 nan 5 nan 8 8];
hpc_lfp = [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];
hpc_mua = [2 nan nan 2 2 nan nan nan 2 nan 3 nan nan 2];
% distal_usable = [1 2 4 5 8];
distal_usable = [1 2 4 5 8 9 11 13 14];
save distal_dir 