clear all
close all

drive_letter = 'G';
addpath(strcat(drive_letter,':\Code\smoothing\software'));
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\KPMstats\'));
addpath(strcat(drive_letter,':\Code\FullBNT-1.0.4\netlab3.3'));
addpath(strcat(drive_letter,':\WC_Germany\new_stellate_analysis\'));
addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'));
addpath(strcat(drive_letter,':\WC_Germany\parietal_cortical_2010\'));

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
% desynch_start_times(interneurons) = [];
% desynch_stop_times(interneurons) = [];
desynch_times(interneurons) = [];

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

% sess_data = sess_data(parietal);
% desynch_start_times = desynch_start_times(parietal);
% desynch_stop_times = desynch_stop_times(parietal);
% %get rid of interneurons
% interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
% sess_data(interneurons) = [];
% desynch_start_times(interneurons) = [];
% desynch_stop_times(interneurons) = [];

thomas_el = find_struct_field_vals(sess_data,'thom_elec',1);
% sess_data = sess_data(thomas_el);


for d = 1:length(sess_data)
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    
    cd(cdir)
    
    %load in MP state sequence and compute desynch times
    load hsmm_state_seq_seg_lf_4_10_10
%     desynch_times = [desynch_start_times{d}; desynch_stop_times{d}];
    time = (1:length(hsmm_bbstate_seq))/Fs_bb;
    is_desynch = zeros(size(time));
    if size(desynch_times{d},1) > 0
        for i = 1:size(desynch_times,1)
            is_desynch(time > desynch_times{d}(i,1) & time < desynch_times{d}(i,2)) = 1;
        end
    end
    is_desynch = logical(is_desynch);
%     time2 = (1:length(hsmm_state_seq))/Fs;
%     is_desynch2 = zeros(size(time2));
%     for i = 1:size(desynch_times,1)
%         is_desynch2(time2 > desynch_times(1,i) & time2 < desynch_times(2,i)) = 1;
%     end
%     is_desynch2 = logical(is_desynch2);
%     time2(1) = time(1);
%     time2(end) = time(end);

    
    load hsmm_state_seq8_seg_lf_4_10_10
%     [hsmm_bbstate_seq8_shift] = shift_state_trans(hsmm_bbstate_seq8,round(0.1*Fs_bb),round(0.04*Fs_bb));
    state_seq1 = hsmm_bbstate_seq;
    state_seq2 = hsmm_bbstate_seq8;
    [up_hist_w8(d,:,:),down_hist_w8(d,:,:),lag_hist_w8(d,:,:),rob_stats_w8(d)] = ...
        compare_synch_mp_lfp_states(state_seq1,state_seq2,desynch_times{d},Fs_bb,sess_data(d).name);
    hamming_w8(d) = pdist([state_seq1(~is_desynch)';state_seq2(~is_desynch)'],'hamming');

    if sess_data(d).thom_elec
        load hsmm_state_seq4_seg_lf_4_10_10
        %     [hsmm_bbstate_seq8_shift] = shift_state_trans(hsmm_bbstate_seq8,round(0.1*Fs_bb),round(0.04*Fs_bb));
        state_seq1 = hsmm_bbstate_seq8;
        state_seq2 = hsmm_bbstate_seq4;
        [up_hist_84(d,:,:),down_hist_84(d,:,:),lag_hist_84(d,:,:),rob_stats_84(d)] = ...
            compare_synch_mp_lfp_states(state_seq1,state_seq2,desynch_times{d},Fs_bb,sess_data(d).name);
        hamming_84(d) = pdist([state_seq1(~is_desynch)';state_seq2(~is_desynch)'],'hamming');
        
        state_seq1 = hsmm_bbstate_seq;
        state_seq2 = hsmm_bbstate_seq4;
        [up_hist_w4(d,:,:),down_hist_w4(d,:,:),lag_hist_w4(d,:,:),rob_stats_w4(d)] = ...
            compare_synch_mp_lfp_states(state_seq1,state_seq2,desynch_times{d},Fs_bb,sess_data(d).name);
        hamming_w4(d) = pdist([state_seq1(~is_desynch)';state_seq2(~is_desynch)'],'hamming');
    end
    
end

lag_axis = [-0.5 0.5];

dist_range = 0:10/Fs_bb:3;
lag_range = linspace(lag_axis(1),lag_axis(2),400);

cd G:\WC_Germany\parietal_cortical_2010
cd figures\


% cd G:\WC_Germany\parietal_cortical_2010
% save lfp_state_compare_synch *hist* hamming* *range

% figure
% plot(hamming_lf(parietal),hamming_hf(parietal),'.')
% line([0.05 0.3],[0.05 0.3],'Color','k')
% shg
% hold on
% plot(hamming_lf(prefrontal),hamming_hf(prefrontal),'r.')
% plot(hamming_lf(frontal),hamming_hf(frontal),'k.')


% figure
% errorbar(lag_range,mean(lag_hist_lf(parietal,1,:)),std(lag_hist_lf(parietal,1,:))/sqrt(length(parietal)))
% hold on
% errorbar(lag_range,mean(lag_hist_hf(parietal,1,:)),std(lag_hist_hf(parietal,1,:))/sqrt(length(parietal)),'r')
% % errorbar(lag_range,mean(lag_hist_lf(prefrontal,1,:)),std(lag_hist_lf(prefrontal,1,:))/sqrt(length(prefrontal)),'r')
% % errorbar(lag_range,mean(lag_hist_lf(frontal,1,:)),std(lag_hist_lf(frontal,1,:))/sqrt(length(frontal)),'k')
% yl = ylim;
% line([0 0],yl,'color','k')
% xlim(lag_axis)
% % ylim([0 0.25])
% 
% figure
% errorbar(lag_range,mean(lag_hist_lf(parietal,2,:)),std(lag_hist_lf(parietal,2,:))/sqrt(length(parietal)))
% hold on
% errorbar(lag_range,mean(lag_hist_hf(parietal,2,:)),std(lag_hist_lf(parietal,2,:))/sqrt(length(parietal)),'r')
% % errorbar(lag_range,mean(lag_hist_lf(prefrontal,2,:)),std(lag_hist_lf(prefrontal,2,:))/sqrt(length(prefrontal)),'r')
% % errorbar(lag_range,mean(lag_hist_lf(frontal,2,:)),std(lag_hist_lf(frontal,2,:))/sqrt(length(frontal)),'k')
% yl = ylim;
% line([0 0],yl,'color','k')
% xlim(lag_axis)
% % ylim([0 0.25])
% 
% figure
% errorbar(dist_range,mean(up_hist_lf(parietal,1,:)),std(up_hist_lf(parietal,1,:))/sqrt(length(parietal)),'g')
% hold on
% errorbar(dist_range,mean(up_hist_lf(parietal,2,:)),std(up_hist_lf(parietal,2,:))/sqrt(length(parietal)))
% errorbar(dist_range,mean(up_hist_hf(parietal,2,:)),std(up_hist_hf(parietal,2,:))/sqrt(length(parietal)),'r')
% % errorbar(dist_range,mean(up_hist_lf(prefrontal,2,:)),std(up_hist_lf(prefrontal,2,:))/sqrt(length(prefrontal)),'r')
% % errorbar(dist_range,mean(up_hist_lf(frontal,2,:)),std(up_hist_lf(frontal,2,:))/sqrt(length(frontal)),'k')
% % errorbar(dist_range,mean(up_hist_lfhf(:,2,:)),std(up_hist_lfhf(:,2,:))/sqrt(n),'k')
% yl = ylim;
% line([0 0],yl,'color','k')
% xlim([0 4])
% 
% figure
% errorbar(dist_range,mean(down_hist_lf(parietal,1,:)),std(down_hist_lf(parietal,1,:))/sqrt(length(parietal)),'g')
% hold on
% errorbar(dist_range,mean(down_hist_lf(parietal,2,:)),std(down_hist_lf(parietal,2,:))/sqrt(length(parietal)))
% errorbar(dist_range,mean(down_hist_hf(parietal,2,:)),std(down_hist_hf(parietal,2,:))/sqrt(length(parietal)),'r')
% % errorbar(dist_range,mean(down_hist_lf(prefrontal,2,:)),std(down_hist_lf(prefrontal,2,:))/sqrt(length(prefrontal)),'r')
% % errorbar(dist_range,mean(down_hist_lf(frontal,2,:)),std(down_hist_lf(frontal,2,:))/sqrt(length(frontal)),'k')
% % errorbar(dist_range,mean(down_hist_lfhf(:,2,:)),std(down_hist_lfhf(:,2,:))/sqrt(n),'k')
% yl = ylim;
% line([0 0],yl,'color','k')
% xlim([0 4])

%% 

for i = 1:length(sess_data)
    var_up_lag_w8(i) = rob_stats_w8(i).var_up_lag;
    var_down_lag_w8(i) = rob_stats_w8(i).var_down_lag;
    mean_up_lag_w8(i) = rob_stats_w8(i).mean_up_lag;
    mean_down_lag_w8(i) = rob_stats_w8(i).mean_down_lag;
    if sess_data(i).thom_elec
        var_up_lag_w4(i) = rob_stats_w4(i).var_up_lag;
        var_down_lag_w4(i) = rob_stats_w4(i).var_down_lag;
        mean_up_lag_w4(i) = rob_stats_w4(i).mean_up_lag;
        mean_down_lag_w4(i) = rob_stats_w4(i).mean_down_lag;
        var_up_lag_84(i) = rob_stats_84(i).var_up_lag;
        var_down_lag_84(i) = rob_stats_84(i).var_down_lag;
        mean_up_lag_84(i) = rob_stats_84(i).mean_up_lag;
        mean_down_lag_84(i) = rob_stats_84(i).mean_down_lag;
    end
    depth(i) = sess_data(i).depth;
    ant_post(i) = sess_data(i).ant_post;
    lateral(i) = sess_data(i).lateral;
end
thom_par = thomas_el(find(ismember(thomas_el,parietal)));
thom_pfc = setdiff(thomas_el,thom_par);
parietal = find_struct_field_vals(sess_data,'region','parietal');
frontal = setdiff(1:length(sess_data),parietal);
superficial = find_struct_field_vals(sess_data,'layer','23');
deep = setdiff(1:length(sess_data),superficial);