clear all
load C:\WC_Germany\Persistent_activity\pyr_heka_dir

%% add new MEC
f_loc{18} = 'C:\WC_Germany\april_09_data\2009-04-07\2009-4-7-19\2009-04-07_CWC_LFP_spontaneous';
f_loc{19} = 'C:\WC_Germany\april_09_data\2009-04-13_A\2009-4-13-18\2009-04-13_CWC_LFP_A_spontaneous';
f_loc{20} = 'C:\WC_Germany\april_09_data\2009-04-13_B\2009-04-13_CWC_LFP_B_spontaneous';
% f_loc{21} = 'C:\WC_Germany\april_09_data\2009-04-05_CWC_LFP\2009-4-13-20-41-27_spontaneous';

%% add new
f_loc{21} = 'C:\WC_Germany\EC_MPasci\2005-01-31_HC_LFP';
f_loc{22} = 'C:\WC_Germany\EC_MPasci\2005-11-25_CWC_LFP_A';
f_loc{23} = 'C:\WC_Germany\EC_MPasci\2005-12-12_CWC_LFP_C';
f_loc{24} = 'C:\WC_Germany\EC_MPasci\2005-12-09_CWC_LFP_A';

f_loc{25} = 'C:\WC_Germany\april_09_data\2009-04-21_CWC_LFP\2009-4-21_12s_spontaneous';
f_loc{26} = 'C:\WC_Germany\april_09_data\2009-05-10_CWC_LFP_B\2009-5-10_B_spontaneous';
f_loc{27} = 'C:\WC_Germany\april_09_data\2009-05-16_CWC_LFP\2009-5-16_spontaneous';

%% add stellate 



f_names{18} = '2009-4-7';
f_names{19} = '2009-4-13_A';
f_names{20} = '2009-4-13_B';
% f_names{21} = '2009-4-5';

f_names{21} = '2005-1-31';
f_names{22} = '2005-11-25_A';
f_names{23} = '2005-12-12_C';
f_names{24} = '2005-12-09_A';

f_names{25} = '2009-4-21';
f_names{26} = '2009-5-10';
f_names{27} = '2009-5-16';


cd C:\WC_Germany\persistent_revised
save pers_revised_heka_dir

