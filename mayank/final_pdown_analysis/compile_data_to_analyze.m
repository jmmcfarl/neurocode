clear all

% cd C:\WC_Germany\sven_thomas_combined\
cd ~/Analysis/Mayank/sven_thomas_combined/
load ./combined_dir_nd.mat
data_cnt = 0;

%add pyramidal L3MEC
for ii = 1:length(l3mec)
    data_cnt = data_cnt + 1;
    data(data_cnt).id = data_cnt;
    data(data_cnt).dir = combined_dir{l3mec(ii)};
    data(data_cnt).loc = 'MEC';
    data(data_cnt).layer = 3;
    data(data_cnt).ctype = 'pyr';
    if ismember(l3mec(ii),old_data_inds)
        data(data_cnt).exp = 'T'; %thomas recs
    else
        data(data_cnt).exp = 'S'; %svens recs
    end
    data(data_cnt).ep = Inf; 
    data(data_cnt).dp = Inf;
    data(data_cnt).hpc_lfp = hpc_lfp(l3mec(ii));
    data(data_cnt).heka_data = combined_heka{l3mec(ii)};
    data(data_cnt).heka_type = combined_heka_type{l3mec(ii)};
    if ismember(l3mec(ii),old_data_inds)
        data(data_cnt).is_old_type = true;
    else
        data(data_cnt).is_old_type = false;
    end
%     data(data_cnt).hpc_mua = hpc_mua(l3mec(ii));
end

%add non-pyramidal L3MEC
l3mec_np(ismember(l3mec_np,[68 69 71])) = [];
for ii = 1:length(l3mec_np)
    data_cnt = data_cnt + 1;
    data(data_cnt).id = data_cnt;
    data(data_cnt).dir = combined_dir{l3mec_np(ii)};
    data(data_cnt).loc = 'MEC';
    data(data_cnt).layer = 3;
    data(data_cnt).ctype = 'nonpyr';
    if ismember(l3mec_np(ii),old_data_inds)
        data(data_cnt).exp = 'T';
    else
        data(data_cnt).exp = 'S';
    end
    data(data_cnt).ep = Inf;
    data(data_cnt).dp = Inf;
    data(data_cnt).hpc_lfp = hpc_lfp(l3mec_np(ii));
    data(data_cnt).heka_data = combined_heka{l3mec_np(ii)};
    data(data_cnt).heka_type = combined_heka_type{l3mec_np(ii)};
    if ismember(l3mec_np(ii),old_data_inds)
        data(data_cnt).is_old_type = true;
    else
        data(data_cnt).is_old_type = false;
    end
%     data(data_cnt).hpc_mua = hpc_mua(l3mec_np(ii));
end

%add L3LEC pyramidal
for ii = 1:length(l3lec)
    data_cnt = data_cnt + 1;
     data(data_cnt).id = data_cnt;
    data(data_cnt).dir = combined_dir{l3lec(ii)};
    data(data_cnt).loc = 'LEC';
    data(data_cnt).layer = 3;
    data(data_cnt).ctype = 'pyr';
    if ismember(l3lec(ii),old_data_inds)
        data(data_cnt).exp = 'T';
    else
        data(data_cnt).exp = 'S';
    end
    data(data_cnt).ep = Inf;
    data(data_cnt).dp = Inf;
    data(data_cnt).hpc_lfp = hpc_lfp(l3lec(ii));
    data(data_cnt).heka_data = combined_heka{l3lec(ii)};
    data(data_cnt).heka_type = combined_heka_type{l3lec(ii)};
    if ismember(l3lec(ii),old_data_inds)
        data(data_cnt).is_old_type = true;
    else
        data(data_cnt).is_old_type = false;
    end
%     data(data_cnt).hpc_mua = hpc_mua(l3lec(ii));
end
%add L3LEC non-pyramidal
for ii = 1:length(l3lec_np)
    data_cnt = data_cnt + 1;
     data(data_cnt).id = data_cnt;
   data(data_cnt).dir = combined_dir{l3lec_np(ii)};
    data(data_cnt).loc = 'LEC';
    data(data_cnt).layer = 3;
    data(data_cnt).ctype = 'nonpyr';
    if ismember(l3lec_np(ii),old_data_inds)
        data(data_cnt).exp = 'T';
    else
        data(data_cnt).exp = 'S';
    end
    data(data_cnt).ep = Inf;
    data(data_cnt).dp = Inf;
    data(data_cnt).hpc_lfp = hpc_lfp(l3lec_np(ii));
    data(data_cnt).heka_data = combined_heka{l3lec_np(ii)};
    data(data_cnt).heka_type = combined_heka_type{l3lec_np(ii)};
    if ismember(l3lec_np(ii),old_data_inds)
        data(data_cnt).is_old_type = true;
    else
        data(data_cnt).is_old_type = false;
    end
%     data(data_cnt).hpc_mua = hpc_mua(l3lec_np(ii));
end

%Load new data from Sven
% cd C:\WC_Germany\persistent_downs\
cd ~/Analysis/Mayank/persistent_downs/
load ./new_pdown_dir
% cur_uset = find(new_pdown_use == 1); %data marked for use
cur_uset = 1:length(new_pdown_use);
for ii = 1:length(cur_uset)
    data_cnt = data_cnt + 1;
    data(data_cnt).id = data_cnt;
    data(data_cnt).dir = new_pdown_dir{cur_uset(ii)};
    data(data_cnt).loc = 'MEC';
    data(data_cnt).layer = nan;
    data(data_cnt).ctype = nan;
    data(data_cnt).exp = 'S';
    data(data_cnt).ep = new_pdown_ep(cur_uset(ii)); %end of reliable part of MP recording
    data(data_cnt).dp = new_pdown_dp(cur_uset(ii)); %drug injection time if using any drugs
    data(data_cnt).hpc_lfp = 2;
    data(data_cnt).heka_data = [data(data_cnt).dir '\heka_data.mat'];
    data(data_cnt).heka_type = 'cont';    
    data(data_cnt).is_old_type = false;
%     data(data_cnt).hpc_mua = new_pdown_hpcmua(cur_uset(ii))-1;
end

%% some classifications on the compiled data
clear_l3pyr = [1:41]; 
clear_l3 = [1:46];
clear_l3pyr = [clear_l3pyr 66 68 73 74 77 79 82 86 89 90 94 97 100 102 103];
clear_l3 = [clear_l3 65 66 68 69 70 73 74 77 79 80 82 86 88 89 90 94 96 97 100 102 103];
no_cell = [71 75 78 85 104 105];
unclear_uds = ...
[67  %Unstable, rec too short
76 %unclear uds during pre-atropine period (not enough LFP UDS)
83 %bad LFP artifacts
84 %bad LFP (short set before atropine, and saturation in LFP signal)
93 %unclear MP UDS (also unclear layer)
91 %bad LFP (saturation in LFP)
92 %unclear LFP and MP UDS
101 %LFP and MP UDS unclear/unstable. not the worst, but also barely long enough to begin with
104 %MP unclear/unstable [drug app anyways]
105 %MP unclear/unstable [drug app anyways]
106]; %bad LFP [drug app at 390]

%% store layer/type classifications
for ii = 1:length(data)
   if ismember(ii,clear_l3)
      data(ii).layer = 3; 
   end
   if ismember(ii,clear_l3pyr)
       data(ii).ctype = 'pyr';
   end
end

%% MEC Layer 2 [stellate] 
base_dir = 'C';
cnt = 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data\2006-04-05_CWC_LFP_A\2006-4-5_19-10-27');
l2data(cnt).heka_data = strcat(base_dir,':\wc_data\MPascii\A2006_04_05_CWC_LFP_A');
cnt = cnt + 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data\2007-05-07_CWC_LFP\2007-5-7_17-52-33');
l2data(cnt).heka_data = strcat(base_dir,':\wc_data\MPascii\A2007_05_07_CWC_LFP');
cnt = cnt + 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data\2007-05-08_CWC_LFP_A\2007-5-8_13-21-1');
l2data(cnt).heka_data = strcat(base_dir,':\wc_data\MPascii\A2007_05_08_CWC_LFP_A');
cnt = cnt + 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data\2007-05-24_CWC_LFP_A\2007-5-24_12-48-28');
l2data(cnt).heka_data = strcat(base_dir,':\wc_data\MPascii\A2007_05_24_CWC_LFP_A');
cnt = cnt + 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data\2007-05-28_CWC_LFP_A\2007-5-28_18-27-1');
l2data(cnt).heka_data = strcat(base_dir,':\wc_data\MPascii\A2007_05_28_CWC_LFP_A');
cnt = cnt + 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data\2007-05-30_CWC_LFP_A\2007-5-30_13-48-19');
l2data(cnt).heka_data = strcat(base_dir,':\wc_data\MPascii\A2007_05_30_CWC_LFP_A');
cnt = cnt + 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data\2007-06-01_CWC_LFP_B\2007-6-1_15-56-41');
l2data(cnt).heka_data = strcat(base_dir,':\wc_data\MPascii\A2007_06_01_CWC_LFP_B');
cnt = cnt + 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data\2007-06-01_CWC_LFP_C\2007-6-1_19-7-5');
l2data(cnt).heka_data = strcat(base_dir,':\wc_data\MPascii\A2007_06_01_CWC_LFP_C');
cnt = cnt + 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data\2007-06-03_CWC_LFP_A\2007-6-3_19-1-5');
l2data(cnt).heka_data = strcat(base_dir,':\wc_data\MPascii\A2007_06_03_CWC_LFP_A');
cnt = cnt + 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data\2007-06-04_CWC_LFP_A\2007-6-4_14-46-17');
l2data(cnt).heka_data = strcat(base_dir,':\wc_data\MPascii\A2007_06_04_CWC_LFP_A');
cnt = cnt + 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data\2007-06-26_CWC_LFP_B\2007-06-26_CWC_LFP_B\2007-6-26_18-47-13');
l2data(cnt).heka_data = strcat(base_dir,':\wc_data\MPascii\A2007_06_26_CWC_LFP_B');
cnt = cnt + 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data\2007-07-04_CWC_LFP_A\2007-07-04_CWC_LFP_A\2007-7-4_13-21-34');
l2data(cnt).heka_data = strcat(base_dir,':\wc_data\MPascii\A2007_07_04_CWC_LFP_A');
cnt = cnt + 1;

for ii = 1:length(l2data)
    l2data(ii).class_cert = 1; %all these recs have certain cell classification from Urliste
    l2data(ii).ctype = 'stellate'; %the above are all clear L2 stellates
end

%now add unclear cell types, or non-stellate l2
l2data(cnt).dir = strcat(base_dir,':\wc_data\2006-03-26_CWC_LFP_B\2006-3-26_21-13-56');
l2data(cnt).heka_data = nan; %havent added yet
l2data(cnt).ctype = 'pyr';
l2data(cnt).class_cert = 0;
cnt = cnt + 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data/2006-03-29_CWC_LFP/2006-3-29_23-1-1');
l2data(cnt).heka_data = nan; %havent added yet
l2data(cnt).ctype = 'unclear';
l2data(cnt).class_cert = 0;
cnt = cnt + 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data/2006-04-05_CWC_LFP_B/2006-4-5_21-23-32');
l2data(cnt).heka_data = nan; %havent added yet
l2data(cnt).ctype = 'unclear';
l2data(cnt).class_cert = 0;
cnt = cnt + 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data/2006-04-07_CWC_LFP_A/2006-4-7_14-39-42');
l2data(cnt).heka_data = nan; %havent added yet
l2data(cnt).ctype = 'interneuron';
l2data(cnt).class_cert = 0;
cnt = cnt + 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data/2006-04-07_CWC_LFP_B/2006-4-7_15-24-5');
l2data(cnt).heka_data = nan; %havent added yet
l2data(cnt).ctype = 'stellate';
l2data(cnt).class_cert = 0;
cnt = cnt + 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data/2006-08-29_CWC_LFP/2006-8-29_20-46-14');
l2data(cnt).heka_data = nan; %havent added yet
l2data(cnt).ctype = 'stellate';
l2data(cnt).class_cert = 0;
cnt = cnt + 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data/2006-09-01_CWC_LFP/2006-9-1_18-48-41');
l2data(cnt).heka_data = nan; %havent added yet
l2data(cnt).ctype = 'stellate';
l2data(cnt).class_cert = 0;
cnt = cnt + 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data/2007-05-23_CWC_LFP_A/2007-5-23_18-10-12');
l2data(cnt).heka_data = nan; %havent added yet
l2data(cnt).ctype = 'stellate';
l2data(cnt).class_cert = 0;
cnt = cnt + 1;

l2data(cnt).dir = strcat(base_dir,':\wc_data/2007-05-23_CWC_LFP_A/2007-5-23_18-10-12');
l2data(cnt).heka_data = nan; %havent added yet
l2data(cnt).ctype = 'stellate';
l2data(cnt).class_cert = 0;
cnt = cnt + 1;

% cell locations
for ii = 1:length(l2data)
    data(data_cnt + ii).id = data_cnt + ii;
    data(data_cnt + ii).dir = l2data(ii).dir;
    data(data_cnt + ii).loc = 'MEC';
    data(data_cnt + ii).layer = 2;
    data(data_cnt + ii).ctype = l2data(ii).ctype;
    data(data_cnt + ii).exp = 'T'; %thomas' recs'
    data(data_cnt + ii).ep = Inf; %recordings good til end
    data(data_cnt + ii).dp = Inf; %no drugs
    data(data_cnt + ii).heka_data = l2data(ii).heka_data;
    data(data_cnt + ii).heka_type = 'sweep';
    data(data_cnt + ii).is_old_type = true;
end
data_cnt = length(data);


%% L2 LEC
cnt = 1;
l2lec(cnt).dir = strcat(base_dir,':\wc_data\2005-02-03_HC_LFP\2005-2-3_17-11-9');
l2lec(cnt).heka_data = strcat(base_dir,':\wc_data\MPascii\A2005-02-03_HC_LFP');
l2lec(cnt).heka_type = 'sweep';
l2lec(cnt).ctype = 'fan';
l2lec(cnt).class_certainty = 0;
cnt = cnt + 1;

l2lec(cnt).dir = strcat(base_dir,':\wc_data\2005-12-09_CWC_LFP_B\2005-12-9_14-8-51');
l2lec(cnt).heka_dir = strcat(base_dir,':\wc_data\MPascii\A2005-12-09_CWC_LFP_B');
l2lec(cnt).heka_type = 'sweep';
l2lec(cnt).cell_type = 'multipolar';
l2lec(cnt).class_certainty = 0;
cnt = cnt + 1;

l2lec(cnt).dir = strcat(base_dir,':\wc_data\2005-11-28_CWC_LFP\2005-11-28_17-11-19');
l2lec(cnt).heka_dir = strcat(base_dir,':\wc_data\MPascii\A2005-11-28_CWC_LFP');
l2lec(cnt).heka_type = 'sweep';
l2lec(cnt).cell_type = 'fan';
l2lec(cnt).class_certainty = 0;
cnt = cnt + 1;

l2lec(cnt).dir = strcat(base_dir,':\wc_data\2005-11-29_CWC_LFP_A\2005-11-29_14-27-37');
l2lec(cnt).heka_dir = strcat(base_dir,':\wc_data\MPascii\A2005-11-29_CWC_LFP_A');
l2lec(cnt).heka_type = 'sweep';
l2lec(cnt).cell_type = 'multipolar';
l2lec(cnt).class_certainty = 1;
cnt = cnt + 1;

l2lec(cnt).dir = strcat(base_dir,':\wc_data\2005-11-29_CWC_LFP_B\2005-11-29_17-47-36');
l2lec(cnt).heka_dir = strcat(base_dir,':\wc_data\MPascii\A2005-11-29_CWC_LFP_B');
l2lec(cnt).heka_type = 'sweep';
l2lec(cnt).cell_type = 'multipolar';
l2lec(cnt).class_certainty = 0;
cnt = cnt + 1;

l2lec(cnt).dir = strcat(base_dir,':\wc_data\2007-01-17_CWC_LFP\2007-1-17_17-42-58');
l2lec(cnt).heka_dir = strcat(base_dir,':\wc_data\MPascii\A2007-01-17_CWC_LFP');
l2lec(cnt).heka_type = 'sweep';
l2lec(cnt).cell_type = 'fan';
l2lec(cnt).class_certainty = 1;
cnt = cnt + 1;

l2lec(cnt).dir = strcat(base_dir,':\wc_data\2009-04-24_CWC_LFP\2009-4-24_22-16-59');
l2lec(cnt).heka_dir = strcat(base_dir,':\wc_data\MPascii\2009-04-24_CWC_LFP');
l2lec(cnt).heka_type = 'cont';
l2lec(cnt).cell_type = 'multipolar';
l2lec(cnt).class_certainty = 0;
cnt = cnt + 1;

l2lec(cnt).dir = strcat(base_dir,':\wc_data\2009-04-29_CWC_LFP_A\2009-4-29_15-57-12');
l2lec(cnt).heka_dir = strcat(base_dir,':\wc_data\MPascii\2009-04-29_CWC_LFP_A');
l2lec(cnt).heka_type = 'cont';
l2lec(cnt).cell_type = 'multipolar';
l2lec(cnt).class_certainty = 1;
cnt = cnt + 1;

l2lec(cnt).dir = strcat(base_dir,':\wc_data\2009-04-29_CWC_LFP_B\2009-4-29_16-56-26');
l2lec(cnt).heka_dir = strcat(base_dir,':\wc_data\MPascii\2009-04-29_CWC_LFP_B');
l2lec(cnt).heka_type = 'cont';
l2lec(cnt).cell_type = 'fan';
l2lec(cnt).class_certainty = 1;
cnt = cnt + 1;

l2lec(cnt).dir = strcat(base_dir,':\wc_data\2010-05-26_CWC_LFP');
l2lec(cnt).heka_dir = strcat(base_dir,':\wc_data\MPascii\2010-05-26_CWC_LFP');
l2lec(cnt).heka_type = 'cont';
l2lec(cnt).cell_type = 'fan';
l2lec(cnt).class_certainty = 1;
cnt = cnt + 1;

l2lec(cnt).dir = strcat(base_dir,':\wc_data\2010-07-31_CWC_LFP_A\2010-7-31_16-49-36');
l2lec(cnt).heka_dir = strcat(base_dir,':\wc_data\MPascii\2010-07-31_CWC_LFP_A');
l2lec(cnt).heka_type = 'cont';
l2lec(cnt).cell_type = 'unclear';
l2lec(cnt).class_certainty = 0;
cnt = cnt + 1;

l2lec(cnt).dir = strcat(base_dir,':\wc_data\2010-07-31_CWC_LFP_B\2010-7-31_17-36-25');
l2lec(cnt).heka_dir = strcat(base_dir,':\wc_data\MPascii\2010-07-31_CWC_LFP_B');
l2lec(cnt).heka_type = 'cont';
l2lec(cnt).cell_type = 'fan';
l2lec(cnt).class_certainty = 1;
cnt = cnt + 1;

l2lec(cnt).dir = strcat(base_dir,':\wc_data\2010-08-01_CWC_LFP_B\2010-8-1_16-28-28');
l2lec(cnt).heka_dir = strcat(base_dir,':\wc_data\MPascii\2010-08-01_CWC_LFP_B');
l2lec(cnt).heka_type = 'cont';
l2lec(cnt).cell_type = 'multipolar';
l2lec(cnt).class_certainty = 1;
cnt = cnt + 1;

l2lec(cnt).dir = strcat(base_dir,':\wc_data\2010-08-05_CWC_LFP_A\2010-8-5_17-21-0');
l2lec(cnt).heka_dir = strcat(base_dir,':\wc_data\MPascii\2010-08-05_CWC_LFP_A');
l2lec(cnt).heka_type = 'cont';
l2lec(cnt).cell_type = 'pyramidal';
l2lec(cnt).class_certainty = 0;
cnt = cnt + 1;

l2lec(cnt).dir = strcat(base_dir,':\wc_data\2010-08-17_CWC_LFP_A\2010-8-17_17-58-5');
l2lec(cnt).heka_dir = strcat(base_dir,':\wc_data\MPascii\2010-08-17_CWC_LFP_A');
l2lec(cnt).heka_type = 'cont';
l2lec(cnt).cell_type = 'fan';
l2lec(cnt).class_certainty = 1;
cnt = cnt + 1;

% cell locations
for ii = 1:length(l2lec)
    data(data_cnt + ii).id = data_cnt + ii;
    data(data_cnt + ii).dir = l2lec(ii).dir;
    data(data_cnt + ii).loc = 'LEC';
    data(data_cnt + ii).layer = 2;
    data(data_cnt + ii).ctype = l2lec(ii).cell_type;
    data(data_cnt + ii).exp = 'T'; %thomas' recs'
    data(data_cnt + ii).ep = Inf; %recordings good til end
    data(data_cnt + ii).dp = Inf; %no drugs
    data(data_cnt + ii).heka_data = l2lec(ii).heka_data;
    data(data_cnt + ii).heka_type = l2lec(ii).heka_type;
    data(data_cnt + ii).is_old_type = true;
end
data_cnt = length(data);

%%
% cd C:\WC_Germany\final_pdown_analysis\
cd ~/Analysis/Mayank/final_pdown_analysis
save compiled_data data clear* no_cell unclear_uds
