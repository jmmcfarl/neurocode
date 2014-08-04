clear all

load C:\WC_Germany\JMM_Analysis_pyr\UDS_dur_run_hist_v2\UDS_dur_data_over_smooth
load C:\WC_Germany\JMM_Analysis_pyr\dir_tree_update
% load C:\WC_Germany\JMM_Analysis_pyr\time_domain_30\correl_data
% load C:\WC_Germany\JMM_Analysis_pyr\time_domain_30\lag_data

dsf = 8;
Fsd = 2016/dsf;
winsize = 30;
niqf = 2016/2;
hcf = 40/niqf;
lcf = 0.05/niqf;
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
    npers_secs = setdiff([1:total_dur],pers_secs);
    npers_secs = setdiff(npers_secs,pers_secs8);
    
    pers_measure = zeros(1,total_dur);
    pers_measure(pers_secs) = 1;
    pers_diff = [0 diff(pers_measure)];
    pers_start = find(pers_diff==1);
    pers_stop = find(pers_diff == -1);
    
    npers_measure = zeros(1,total_dur);
    npers_measure(npers_secs) = 1;
    npers_diff = [0 diff(npers_measure)];
    npers_start = find(npers_diff ==1);
    npers_stop = find(npers_diff == -1);
    if npers_measure(1) == 1
        npers_start = [1 npers_start];
    elseif pers_measure(1) == 1
        pers_start = [1 pers_start];
    end
    if npers_measure(end) == 1
        npers_stop = [npers_stop total_dur];
    elseif pers_measure(end) == 1
        pers_stop = [pers_stop total_dur];
    end
    
    pers_dur = pers_stop - pers_start;
    npers_dur = npers_stop - npers_start;
    
    bad_pers = find(pers_dur < min_dur);
    bad_npers = find(npers_dur < min_dur);
    
    pers_start(bad_pers) = [];
    pers_stop(bad_pers) = [];
    npers_start(bad_npers) = [];
    npers_stop(bad_npers) = [];
    
 pers_start_times{d} = pers_start;
 pers_stop_times{d} = pers_stop;
 npers_start_times{d} = npers_start;
 npers_stop_times{d} = npers_stop;
   %         %convert to sample index
        pers_start = pers_start*Fsd;
        pers_stop = pers_stop*Fsd;
        npers_start = npers_start*Fsd;
        npers_stop = npers_stop*Fsd;

        % find persistent ups and persistent downs
        pers_up_inds{d} = [];
        pers_down_inds{d} = [];
        if ~isempty(pers_start)
            for i = 1:length(pers_start);
                pers_up_inds{d} = [pers_up_inds{d} find(up_trans{d} > pers_start(i) & up_trans{d} < pers_stop(i))];
                pers_down_inds{d} = [pers_down_inds{d} find(down_trans{d} > pers_start(i) & down_trans{d} < pers_stop(i))];
            end
        end
        npers_up_inds{d} = find(~ismember(1:length(up_trans{d}),pers_up_inds{d}));
        npers_down_inds{d} = find(~ismember(1:length(down_trans{d}),pers_down_inds{d}));
       
        
end
    
% 
%     pers_measure8 = zeros(1,total_dur);
%     pers_measure8(pers_secs8) = 1;
%     pers_diff8 = [0 diff(pers_measure8)];
%     pers_start8 = find(pers_diff8==1);
%     pers_stop8 = find(pers_diff8 == -1);
%         
%     npers_start = pers_stop;
%     npers_stop = pers_start;
% 
%     if ~isempty(pers_start8)
%         npers_start(ismember(npers_start,pers_start8)) = [];
%         npers_stop(ismember(npers_stop,pers_stop8)) = [];
%         npers_start = [npers_start pers_stop8];
%         npers_stop = [npers_stop pers_start8];
%         npers_start = sort(npers_start);
%         npers_stop = sort(npers_stop);
%     end
    
    
    %if there is some persistent activity
%     if length(pers_start)>0
% 
% %         if pers_measure(1) ==1
% %             pers_start = [1 pers_start];
% %         elseif pers_start(1) > 30
% %             npers_start = [1 npers_start];
% %         end
% %         if pers_measure(end) == 1
% %             pers_stop = [pers_stop total_dur];
% %         elseif total_dur-pers_stop(end) > 30
% %             npers_stop = [npers_stop total_dur];
% %         else
% %             npers_start(end) = [];
% %         end
% % 
% %         pers_dur = pers_stop-pers_start;
% %         npers_dur = npers_stop-npers_start;
% % 
% %         bad_npers = find(npers_dur < min_dur);
% % %         bad_nstart = npers_start(bad_npers);
% % %         bad_nstop = npers_stop(bad_npers);
% %         npers_start(bad_npers) = [];
% %         npers_stop(bad_npers) = [];
% % %         pers_start(ismember(pers_start,bad_nstop)) = [];
% % %         pers_stop(ismember(pers_stop,bad_nstart)) = [];
% % 
% %         bad_pers = find(pers_dur < min_dur);
% % %         bad_start = pers_start(bad_pers);
% % %         bad_stop = pers_stop(bad_pers);
% %         pers_start(bad_pers) = [];
% %         pers_stop(bad_pers) = [];
% %         
% 
% 
%         %convert to sample index
%         pers_start = pers_start*Fsd;
%         pers_stop = pers_stop*Fsd;
%         npers_start = npers_start*Fsd;
%         npers_stop = npers_stop*Fsd;
% 
% 
%         %cycle through pers states
%         %if there is some persistent activity
% 
% 
%         num_pers_states = length(pers_start);
%         num_npers_states = length(npers_start);
% 
%         %initialize
%         pers_w_acorr = zeros(num_pers_states,2*maxlag+1);
%         pers_8_acorr = pers_w_acorr;
%         pers_xcorr = pers_w_acorr;
%         
%         %% Cycle through all persistent states and calculate the UTA and
%         %% DTA
% 
%         for t = 1:num_pers_states
% 
%             state_beg = pers_start(t);
%             state_end = pers_stop(t);
% 
%             wcv_seg = wcv_d(state_beg:state_end);
%             lf8_seg = lf8_d(state_beg:state_end);
%             
%             pers_w_acorr(t,:) = xcov(wcv_seg,maxlag,'coeff');
%             pers_8_acorr(t,:) = xcov(lf8_seg,maxlag,'coeff');
%             pers_xcorr(t,:) = xcov(wcv_seg,lf8_seg,maxlag,'coeff');
%             
%         end
% 
%         if num_pers_states > 1
%         ov_p_w_acorr(d,:) = mean(pers_w_acorr);
%         ov_p_8_acorr(d,:) = mean(pers_8_acorr);
%         ov_p_xcorr(d,:) = mean(pers_xcorr);
%         else
%             ov_p_w_acorr(d,:) = pers_w_acorr;
%             ov_p_8_acorr(d,:) = pers_8_acorr;
%             ov_p_xcorr(d,:) = pers_xcorr;
%         end
%         
%         %% Cycle through non-pers states
% 
%         %initialize
%         npers_w_acorr = zeros(num_npers_states,2*maxlag+1);
%         npers_8_acorr = npers_w_acorr;
%         npers_xcorr = npers_w_acorr;
%         for t = 1:num_npers_states
% 
%             state_beg = npers_start(t);
%             state_end = npers_stop(t);
%             
%             wcv_seg = wcv_d(state_beg:state_end);
%             lf8_seg = lf8_d(state_beg:state_end);
%             
%             npers_w_acorr(t,:) = xcov(wcv_seg,maxlag,'coeff');
%             npers_8_acorr(t,:) = xcov(lf8_seg,maxlag,'coeff');
%             npers_xcorr(t,:) = xcov(wcv_seg,lf8_seg,maxlag,'coeff');
%             
% 
%             
%         end
% 
%         if num_npers_states > 1
%         ov_n_w_acorr(d,:) = mean(npers_w_acorr);
%         ov_n_8_acorr(d,:) = mean(npers_8_acorr);
%         ov_n_xcorr(d,:) = mean(npers_xcorr);
%         else
%            ov_n_w_acorr(d,:) = npers_w_acorr;
%            ov_n_8_acorr(d,:) = npers_8_acorr;
%            ov_n_xcorr(d,:) = npers_xcorr;
%         end
%     else
%         disp('no persistent activity')
%     end
% 
% 
% end
% 
% %get rid of cells without PA
% no_pa = [3 6 12];
% ov_p_w_acorr(no_pa,:) = [];
% ov_p_8_acorr(no_pa,:) = [];
% ov_p_xcorr(no_pa,:) = [];
% ov_n_w_acorr(no_pa,:) = [];
% ov_n_8_acorr(no_pa,:) = [];
% ov_n_xcorr(no_pa,:) = [];
% 
% %averages and plot
% m_p_w = mean(ov_p_w_acorr);
% u_p_w = m_p_w+2*std(ov_p_w_acorr)/sqrt(14);
% l_p_w = m_p_w-2*std(ov_p_w_acorr)/sqrt(14);
% 
% m_n_w = mean(ov_n_w_acorr);
% u_n_w = m_n_w+2*std(ov_n_w_acorr)/sqrt(14);
% l_n_w = m_n_w-2*std(ov_n_w_acorr)/sqrt(14);
% 
% m_p_8 = mean(ov_p_8_acorr);
% u_p_8 = m_p_8+2*std(ov_p_8_acorr)/sqrt(14);
% l_p_8 = m_p_8-2*std(ov_p_8_acorr)/sqrt(14);
% 
% m_n_8 = mean(ov_n_8_acorr);
% u_n_8 = m_n_8+2*std(ov_n_8_acorr)/sqrt(14);
% l_n_8 = m_n_8-2*std(ov_n_8_acorr)/sqrt(14);
% 
% m_p_x = mean(ov_p_xcorr);
% u_p_x = m_p_x+2*std(ov_p_xcorr)/sqrt(14);
% l_p_x = m_p_x-2*std(ov_p_xcorr)/sqrt(14);
% 
% m_n_x = mean(ov_n_xcorr);
% u_n_x = m_n_x+2*std(ov_n_xcorr)/sqrt(14);
% l_n_x = m_n_x-2*std(ov_n_xcorr)/sqrt(14);
% 
% plot(tlags,m_p_w,'linewidth',2)
% hold on
% plot(tlags,m_n_w,'r','linewidth',2)
% legend('Persistent Activity','Non-Persistent Activity')
% plot(tlags,u_p_w,'--')
% plot(tlags,l_p_w,'--')
% plot(tlags,u_n_w,'r--')
% plot(tlags,l_n_w,'r--')
% xlim([0 10])
% xlabel('Time (s)','FontSize',14)
% ylabel('Correlation Coefficient','FontSize',14)
% 
% 
% plot(tlags,m_p_8,'linewidth',2)
% hold on
% plot(tlags,m_n_8,'r','linewidth',2)
% legend('Persistent Activity','Non-Persistent Activity')
% plot(tlags,u_p_8,'--')
% plot(tlags,l_p_8,'--')
% plot(tlags,u_n_8,'r--')
% plot(tlags,l_n_8,'r--')
% xlim([0 10])
% xlabel('Time (s)','FontSize',14)
% ylabel('Correlation Coefficient','FontSize',14)
% 
% 
% plot(tlags,m_p_x,'linewidth',2)
% hold on
% plot(tlags,m_n_x,'r','linewidth',2)
% legend('Persistent Activity','Non-Persistent Activity')
% plot(tlags,u_p_x,'--')
% plot(tlags,l_p_x,'--')
% plot(tlags,u_n_x,'r--')
% plot(tlags,l_n_x,'r--')
% xlim([-5 5])
% xlabel('Time (s)','FontSize',14)
% ylabel('Correlation Coefficient','FontSize',14)
