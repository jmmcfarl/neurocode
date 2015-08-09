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
76 %unclear uds during pre-atropine period
83 %bad LFP
84 %bad LFP
93 %unclear MP UDS
91 %bad LFP
92 %unclear LFP and MP UDS
101 %MP unclear/unstable
104 %MP unclear/unstable
105 %MP unclear/unstable
106]; %bad LFP

%% store layer/type classifications
for ii = 1:length(data)
   if ismember(ii,clear_l3)
      data(ii).layer = 3; 
   end
   if ismember(ii,clear_l3pyr)
       data(ii).ctype = 'pyr';
   end
end
%%
% cd C:\WC_Germany\final_pdown_analysis\
cd ~/Analysis/Mayank/final_pdown_analysis
save compiled_data data clear* no_cell unclear_uds
