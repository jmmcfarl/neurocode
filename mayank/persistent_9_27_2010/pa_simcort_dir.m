clear all
close all

%% simultaneous cortical LFP recordings
dir_array{1} = 'F:\wc_data\2007-04-05_CWC_LFP_A\2007-4-5_14-49-1';
dir_array{2} = 'F:\wc_data\2007-04-05_CWC_LFP_B\2007-4-5_17-55-51';
dir_array{3} = 'F:\wc_data\2007-04-05_CWC_LFP_C\2007-4-5_19-31-0';
dir_array{4} = 'F:\wc_data\2007-04-10_CWC_LFP_A\2007-4-10_18-43-7';
dir_array{5} = 'F:\wc_data\2007-04-10_CWC_LFP_B\2007-4-10_19-35-5';
dir_array{6} = 'F:\wc_data\2007-04-12_CWC_LFP_A\2007-4-12_13-10-50';
dir_array{7} = 'F:\wc_data\2007-04-13_CWC_LFP_A\2007-4-13_12-19-50';
dir_array{8} = 'F:\wc_data\2007-04-13_CWC_LFP_B\2007-4-13_13-45-34';
dir_array{9} = 'F:\wc_data\2007-04-16_CWC_LFP_A\2007-4-16_13-13-21';
dir_array{10} = 'F:\wc_data\2007-04-16_CWC_LFP_B\2007-4-16_14-8-43';
dir_array{11} = 'F:\wc_data\2007-04-16_CWC_LFP_C\2007-4-16_15-18-23';
dir_array{12} = 'F:\wc_data\2007-04-17_CWC_LFP_A\2007-4-17_13-36-49';
dir_array{13} = 'F:\wc_data\2007-04-17_CWC_LFP_B\2007-4-17_14-49-45';
dir_array{14} = 'F:\wc_data\2007-04-17_CWC_LFP_C\2007-4-17_15-41-51';
dir_array{15} = 'F:\wc_data\2007-04-18_CWC_LFP_A\2007-4-18_13-40-23';
dir_array{16} = 'F:\wc_data\2007-04-18_CWC_LFP_B\2007-4-18_14-28-16';

dir_array{17} = 'F:\wc_data\2007-04-11_CWC_LFP_A\2007-4-11_16-43-19';
dir_array{18} = 'F:\wc_data\2007-04-11_CWC_LFP_B\2007-4-11_17-59-41';
dir_array{19} = 'F:\wc_data\2007-04-12_CWC_LFP_B\2007-4-12_16-55-57';
dir_array{20} = 'F:\wc_data\2007-04-12_CWC_LFP_C\2007-4-12_18-6-33';
dir_array{21} = 'F:\wc_data\2007-04-13_CWC_LFP_C\2007-4-13_16-2-15';
dir_array{22} = 'F:\wc_data\2007-04-13_CWC_LFP_D\2007-4-13_17-41-15';
dir_array{23} = 'F:\wc_data\2007-04-16_CWC_LFP_D\2007-4-16_17-34-6';
dir_array{24} = 'F:\wc_data\2007-04-16_CWC_LFP_E\2007-4-16_18-34-42';

f_names{1} = 'pre_2007-4-5_A';
f_names{2} = 'pre_2007-4-5_B';
f_names{3} = 'pre_2007-4-5_C';
f_names{4} = 'pre_2007-4-10_A';
f_names{5} = 'pre_2007-4-10_B';
f_names{6} = 'pre_2007-4-12_A';
f_names{7} = 'pre_2007-4-13_A';
f_names{8} = 'pre_2007-4-13_B';
f_names{9} = 'pre_2007-4-16_A';
f_names{10} = 'pre_2007-4-16_B';
f_names{11} = 'pre_2007-4-16_C';
f_names{12} = 'pre_2007-4-17_A';
f_names{13} = 'pre_2007-4-17_B';
f_names{14} = 'pre_2007-4-17_C';
f_names{15} = 'pre_2007-4-18_A';
f_names{16} = 'pre_2007-4-18_B';

f_names{17} = 'fro_2007-4-11_A';
f_names{18} = 'fro_2007-4-11_B';
f_names{19} = 'fro_2007-4-12_B';
f_names{20} = 'fro_2007-4-12_C';
f_names{21} = 'fro_2007-4-13_C';
f_names{22} = 'fro_2007-4-13_D';
f_names{23} = 'fro_2007-4-16_D';
f_names{24} = 'fro_2007-4-16_E';

pre = [1:16];
fro = [17:24];
%%
cd F:\WC_Germany\persistent_9_27_2010\
save pa_simcort_dir