clear all

distal_lec{1} = 'C:\wc_data\2012_06_21_A\2012-6-21_14-33-23';
distal_lec{2} = 'C:\wc_data\2012_06_21_B\2012-6-21_16-1-44'; %LEC neuron, good LFP, ?MUA, but cell wierd (depolarized MP)
distal_lec{3} = 'C:\wc_data\2012_06_23_A\2012-6-23_18-39-42'; %not usable LFP (MAYBE with LF5)
distal_lec{4} = 'C:\wc_data\2012_06_23_C\2012-6-23_19-39-37'; %usable LF5
distal_lec{5}  = 'C:\wc_data\2012_06_24_A\2012-6-24_14-49-46'; %questionable LFP
distal_lec{6} = 'C:\wc_data\2012_06_24_B\2012-6-24_15-46-53'; %good MUA, questionable LFP (at times), NEED to use LF4
distal_lec{7} = 'C:\wc_data\2012_06_24_C\2012-6-24_16-24-8'; %maybe OK MUA, questionabel LFP, (again, need LF4)
distal_lec{8} = 'C:\wc_data\2012_06_25_A\2012-6-25_12-58-45'; %not usable LFP


cd C:\WC_Germany\sven_thomas_combined
usable_distal_lec = [1 5 6 7 8];
ctx_lfp = [5 nan nan nan 5 5 5 5];
hpc_lfp = [nan nan nan nan nan nan nan nan];
hpc_mua = [nan nan nan nan nan nan nan nan];

save distal_lec_dir
