clear all
% load C:\WC_Germany\JMM_Analysis_pyr\UDS_dur_run_hist_v2\up_per_data_10_28

load C:\WC_Germany\JMM_Analysis_pyr\UDS_dur_run_hist_v2\UDS_dur_data_over_smooth
load C:\WC_Germany\JMM_Analysis_pyr\dir_tree_update
% load C:\WC_Germany\JMM_Analysis_pyr\time_domain_30\correl_data
% load C:\WC_Germany\JMM_Analysis_pyr\time_domain_30\lag_data

dsf = 8;
Fsd = 2016/dsf;
winsize = 30;
niqf = 2016/2;
hcf = 40/niqf;
lcf = 0.005/niqf;
[b,a] = butter(2,[lcf hcf]);
cnt = 0;
maxlag = round(10*Fsd);
tlags = [-maxlag:maxlag]/Fsd;

pers_thresh = 5;
back_time = 10;
forward_time = 30;
min_dur = 30;

for d = 1:length(dir_array)

    cd(dir_array{d})
    pwd

    load used_data lf8 wcv_minus_spike

    lf8 = filtfilt(b,a,lf8);
    wcv = filtfilt(b,a,wcv_minus_spike);

    lf8_d = downsample(lf8,dsf);
    wcv_d = downsample(wcv,dsf);

    lf8_d = zscore(lf8_d);
    wcv_d = zscore(wcv_d);

    total_dur = floor(length(lf8_d)/Fsd);

    %find all wcv up states lasting more than pers_dur seconds
    pers_ups = round(up_trans{d}(find(up_state_dur{d} > pers_thresh))/Fsd);
    pers_ups8 = round(up_trans8{d}(find(up_state_dur8{d} > pers_thresh))/Fsd);

    pers_begs = pers_ups - back_time;
    pers_ends = pers_ups + forward_time;
    pers_begs8= pers_ups8-back_time;
    pers_ends8 = pers_ups8+forward_time;

    pers_secs = [];
    pers_secs8 = [];
    for be = 1:length(pers_begs)
        pers_secs = [pers_secs pers_begs(be):pers_ends(be)];
    end
    for be = 1:length(pers_begs8)
        pers_secs8 = [pers_secs8 pers_begs8(be):pers_ends8(be)];
    end

    pers_secs = unique(pers_secs);
    pers_secs8 = unique(pers_secs8);
pers_secs(pers_secs <= 0) = [];
pers_secs(pers_secs > total_dur) = [];
pers_secs8(pers_secs8 <= 0) = [];
pers_secs8(pers_secs8 > total_dur) = [];

    pers_secs = setdiff(pers_secs,pers_secs8);
    pers_measure = zeros(1,total_dur);
    pers_measure(pers_secs) = 1;
    pers_diff = [0 diff(pers_measure)];
    pers_start = find(pers_diff==1);
    pers_stop = find(pers_diff == -1);

    npers_start = pers_stop;
    npers_stop = pers_start;

    %if there is some persistent activity
    if length(pers_start)>0

        if pers_measure(1) ==1
            pers_start = [1 pers_start];
        elseif pers_start(1) > 30
            npers_start = [1 npers_start];
        end
        if pers_measure(end) == 1
            pers_stop = [pers_stop total_dur];
        elseif total_dur-pers_stop(end) > 30
            npers_stop = [npers_stop total_dur];
        else
            npers_start(end) = [];
        end

        pers_dur = pers_stop-pers_start;
        npers_dur = npers_stop-npers_start;

        bad_npers = find(npers_dur < min_dur);
        bad_nstart = npers_start(bad_npers);
        bad_nstop = npers_stop(bad_npers);
        npers_start(bad_npers) = [];
        npers_stop(bad_npers) = [];
        pers_start(ismember(pers_start,bad_nstop)) = [];
        pers_stop(ismember(pers_stop,bad_nstart)) = [];

        tot_pers_dur = sum(pers_dur);
        tot_npers_dur = sum(npers_dur);
        
        perc_pers(d) = tot_pers_dur/(tot_pers_dur+tot_npers_dur);

        %convert to sample index
        pers_start = pers_start*Fsd;
        pers_stop = pers_stop*Fsd;
        npers_start = npers_start*Fsd;
        npers_stop = npers_stop*Fsd;

        pers_wcv = [];
        npers_wcv = [];
        pers_lf8 = [];
        npers_lf8 = [];

        %cycle through pers states
        %if there is some persistent activity


        num_pers_states = length(pers_start);
        num_npers_states = length(npers_start);

        %initialize
        pers_uta_mat = zeros(num_pers_states,2*maxlag+1);
        pers_dta_mat = pers_uta_mat;

        %% Cycle through all persistent states and calculate the UTA and
        %% DTA

%         for t = 1:num_pers_states
% 
%             state_beg = pers_start(t);
%             state_end = pers_stop(t);
% 
%             %get relevent lfp up transtition indices
%             cur_lf8_ups = up_trans8{d}(find(up_trans8{d} >= state_beg & up_trans8{d} < state_end));
%             cur_lf8_downs = down_trans8{d}(find(down_trans8{d} >= state_beg & down_trans8{d} < state_end));
% 
%             cur_lf8_ups(cur_lf8_ups<maxlag) = [];
%             cur_lf8_ups(total_dur*Fsd-cur_lf8_ups<maxlag) = [];
%             cur_lf8_downs(cur_lf8_downs<maxlag) = [];
%             cur_lf8_downs(total_dur*Fsd-cur_lf8_downs<maxlag) = [];
% 
%             num_8_ups = length(cur_lf8_ups);
%             num_8_downs = length(cur_lf8_downs);
% 
%             if num_8_ups > 0
% 
%                 uta_mat = zeros(num_8_ups,2*maxlag+1);
% 
%                 for u = 1:num_8_ups
% 
%                     uta_mat(u,:) = wcv_d(cur_lf8_ups(u)-maxlag:cur_lf8_ups(u)+maxlag);
% 
%                 end
% 
%                 pers_uta_mat(t,:) = mean(uta_mat);
% 
%             else
% 
%                 pers_uta_mat(t,:) = nan;
% 
%             end
% 
% 
%             if num_8_downs > 0
% 
%                 dta_mat = zeros(num_8_downs,2*maxlag+1);
% 
%                 for w = 1:num_8_downs
% 
%                     dta_mat(w,:) = wcv_d(cur_lf8_downs(w)-maxlag:cur_lf8_downs(w)+maxlag);
% 
%                 end
% 
%                 pers_dta_mat(t,:) = mean(dta_mat);
% 
%             else
% 
%                 pers_dta_mat(t,:) = nan;
% 
%             end
% 
%             clear uta_mat dta_mat
% 
%         end
% 
%         if num_pers_states > 1
%         cell_uta_pers(d,:) = nanmean(pers_uta_mat);
%         cell_dta_pers(d,:) = nanmean(pers_dta_mat);
%         else
%             cell_uta_pers(d,:) = pers_uta_mat;
%             cell_dta_pers(d,:) = pers_dta_mat;
%         end
% 


        %% Cycle through non-pers states

        %initialize
        npers_uta_mat = zeros(num_npers_states,2*maxlag+1);
        npers_dta_mat = npers_uta_mat;

%         for t = 1:num_npers_states
% 
%             state_beg = npers_start(t);
%             state_end = npers_stop(t);
% 
%             %get relevent lfp up transtition indices
%             cur_lf8_ups = up_trans8{d}(find(up_trans8{d} >= state_beg & up_trans8{d} < state_end));
%             cur_lf8_downs = down_trans8{d}(find(down_trans8{d} >= state_beg & down_trans8{d} < state_end));
% 
%             cur_lf8_ups(cur_lf8_ups<maxlag) = [];
%             cur_lf8_ups(total_dur*Fsd-cur_lf8_ups<maxlag) = [];
%             cur_lf8_downs(cur_lf8_downs<maxlag) = [];
%             cur_lf8_downs(total_dur*Fsd-cur_lf8_downs<maxlag) = [];
% 
%             num_8_ups = length(cur_lf8_ups);
% 
%             if num_8_ups > 0
% 
%                 for u = 1:num_8_ups
% 
%                     uta_mat(u,:) = wcv_d(cur_lf8_ups(u)-maxlag:cur_lf8_ups(u)+maxlag);
% 
%                 end
% 
%                 npers_uta_mat(t,:) = mean(uta_mat);
% 
%             else
% 
%                 npers_uta_mat(t,:) = nan;
% 
%             end
% 
%             num_8_downs = length(cur_lf8_downs);
% 
%             if num_8_downs > 0
% 
%                 for w = 1:num_8_downs
% 
%                     dta_mat(w,:) = wcv_d(cur_lf8_downs(w)-maxlag:cur_lf8_downs(w)+maxlag);
% 
%                 end
% 
%                 npers_dta_mat(t,:) = mean(dta_mat);
% 
%             else
% 
%                 npers_dta_mat(t,:) = nan;
% 
%             end
% 
%         end
% 
%         if num_npers_states > 1
%         cell_uta_npers(d,:) = nanmean(npers_uta_mat);
%         cell_dta_npers(d,:) = nanmean(npers_dta_mat);
%         else
%             cell_uta_npers(d,:) = npers_uta_mat;
%             cell_dta_npers(d,:) = npers_dta_mat;
%         end

    else
        disp('no persistent activity')
    end



end

% %get rid of cells without PA
% no_pa = [3 6 12];
% cell_uta_pers(no_pa,:) = [];
% cell_dta_pers(no_pa,:) = [];
% cell_uta_npers(no_pa,:) = [];
% cell_dta_npers(no_pa,:) = [];
% 
% %plot averages
% m_uta_pers = mean(cell_uta_pers);
% u_uta_pers = m_uta_pers+2*std(cell_uta_pers)/sqrt(14);
% l_uta_pers = m_uta_pers-2*std(cell_uta_pers)/sqrt(14);
% 
% m_dta_pers = mean(cell_dta_pers);
% u_dta_pers = m_dta_pers+2*std(cell_dta_pers)/sqrt(14);
% l_dta_pers = m_dta_pers-2*std(cell_dta_pers)/sqrt(14);
% 
% m_uta_npers = mean(cell_uta_npers);
% u_uta_npers = m_uta_npers+2*std(cell_uta_npers)/sqrt(14);
% l_uta_npers = m_uta_npers-2*std(cell_uta_npers)/sqrt(14);
% 
% m_dta_npers = mean(cell_dta_npers);
% u_dta_npers = m_dta_npers+2*std(cell_dta_npers)/sqrt(14);
% l_dta_npers = m_dta_npers-2*std(cell_dta_npers)/sqrt(14);
% 
% plot(tlags,m_uta_pers,'linewidth',2)
% hold on
% plot(tlags,u_uta_pers,'--')
% plot(tlags,l_uta_pers,'--')
% plot(tlags,m_uta_npers,'r','linewidth',2)
% plot(tlags,u_uta_npers,'r--')
% plot(tlags,l_uta_npers,'r--')
% xlim([-5 5])
% 
% plot(tlags,m_dta_pers,'linewidth',2)
% hold on
% plot(tlags,u_dta_pers,'--')
% plot(tlags,l_dta_pers,'--')
% plot(tlags,m_dta_npers,'r','linewidth',2)
% plot(tlags,u_dta_npers,'r--')
% plot(tlags,l_dta_npers,'r--')
% xlim([-5 5])
