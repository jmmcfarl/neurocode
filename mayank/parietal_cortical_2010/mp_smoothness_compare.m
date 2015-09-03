%inputs:
% raw_data
% raw_Fs
clear all
close all

drive_letter = 'F';

addpath('F:\WC_Germany\hmm_state_detect\\')
addpath('F:\WC_Germany\hsmm_state_detection')
addpath('F:\WC_Germany\parietal_cortical_2010')
% addpath('F:\WC_Germany\new_stellate_analysis\')
addpath('F:\WC_Germany\persistent_2010\')

cd F:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load desynch_times_individual

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
Fs = 50.4;

% parietal = find_struct_field_vals(sess_data,'region','parietal');
% sess_data = sess_data(parietal);
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times_mp(interneurons) = [];

min_state_dur = 1;
max_state_dur = round(Fs*20);
dur_range = (1:max_state_dur)/Fs;

dur_hist_bincenters = (min_state_dur:max_state_dur)/Fs;

hcf_vals = [10];
% smooth_vals = [0.025 0.05 0.075 0.1 0.2 0.3];

n = length(sess_data);
for d = 1:n
    d
    direct = sess_data(d).directory;
    direct(1) = 'F';
    cd(direct)
    load used_data wcv_minus_spike
    for i = 1:length(hcf_vals)
        [hmm_bbstate_seq,hmm,state_durations{d,i}] = parietal_get_hmm_lf_test(wcv_minus_spike,raw_Fs,hcf_vals(i),desynch_times_mp{d});
        up_dur_pmf(d,i,:) = hist(state_durations{d,i}{2},dur_hist_bincenters);
        up_dur_pmf(d,i,:) = up_dur_pmf(d,i,:)/sum(up_dur_pmf(d,i,:));
        down_dur_pmf(d,i,:) = hist(state_durations{d,i}{1},dur_hist_bincenters);
        down_dur_pmf(d,i,:) = down_dur_pmf(d,i,:)/sum(down_dur_pmf(d,i,:));
        min_up_dur(d,i) = min(state_durations{d,i}{2});
        min_down_dur(d,i) = min(state_durations{d,i}{1});
        mean_up_dur(d,i) = mean(state_durations{d,i}{2});
        mean_down_dur(d,i) = mean(state_durations{d,i}{1});
        meandiff = [];
        for s = 1:hmm.Nsegs
            meandiff = [meandiff; hmm.state(2).meanfun{s}-hmm.state(1).meanfun{s}];
        end
        meandiff = mean(meandiff);
        covar1 = hmm.state(1).var;
        covar2 = hmm.state(2).var;
        mp_kl_lf(d,i) = gauss_kl_div(meandiff,covar1,covar2);
%         load hsmm_state_seq_seg_lf_pert15
%         ham_dist_mp_lf8(d,i) = compute_state_seq_seg_hamdist(hsmm_bbstate_seq,hmm_bbstate_seq8);
        clear hmm*
    end
end

% for d = 1:n
%     d
%     direct = sess_data(d).directory;
%     direct(1) = 'F';
%     cd(direct)
%     load used_data lf8
%     for i = 1:length(smooth_vals)
%         [hmm_bbstate_seq8,hmm8_hf,state_durations_hf{d,i}] = parietal_get_hmm_hf_test(lf8,raw_Fs,smooth_vals(i),desynch_times{d});
%         up_dur_hf_pmf(d,i,:) = hist(state_durations_hf{d,i}{2},dur_hist_bincenters);
%         up_dur_hf_pmf(d,i,:) = up_dur_hf_pmf(d,i,:)/sum(up_dur_hf_pmf(d,i,:));
%         down_dur_hf_pmf(d,i,:) = hist(state_durations_hf{d,i}{1},dur_hist_bincenters);
%         down_dur_hf_pmf(d,i,:) = down_dur_hf_pmf(d,i,:)/sum(down_dur_hf_pmf(d,i,:));
%         min_up_dur_hf(d,i) = min(state_durations_hf{d,i}{2});
%         min_down_dur_hf(d,i) = min(state_durations_hf{d,i}{1});
%         mean_up_dur_hf(d,i) = mean(state_durations_hf{d,i}{2});
%         mean_down_dur_hf(d,i) = mean(state_durations_hf{d,i}{1});
%         meandiff = [];
%         for s = 1:hmm8_hf.Nsegs
%             meandiff = [meandiff; hmm8_hf.state(2).meanfun{s}-hmm8_hf.state(1).meanfun{s}];
%         end
%         meandiff = mean(meandiff);
%         covar1 = hmm8_hf.state(1).var;
%         covar2 = hmm8_hf.state(2).var;
%         lf8_kl_hf(d,i) = gauss_kl_div(meandiff,covar1,covar2);
%         load hsmm_state_seq_seg_lf_pert15
%         ham_dist_mp_lf8hf(d,i) = compute_state_seq_seg_hamdist(hsmm_bbstate_seq,hmm_bbstate_seq8);
%         clear hmm*
%     end
% end
% 
cd F:\WC_Germany\parietal_cortical_2010
save smoothness_mp_results
%%
% figure
% cmap = colormap(jet(length(smooth_vals)));
% subplot(2,1,1)
% for i = 1:length(smooth_vals)
%     ebar_mat(dur_hist_bincenters,up_dur_hf_pmf(:,i,:),.001,cmap(i,:)), hold on
% end
% xlim([0 3])
% set(gca,'yscale','log')
% ylim([2e-4 5e-2])
% subplot(2,1,2)
% for i = 1:length(smooth_vals)
%     ebar_mat(dur_hist_bincenters,down_dur_hf_pmf(:,i,:),.001,cmap(i,:)), hold on
% end
% xlim([0 3])
% set(gca,'yscale','log')
% ylim([2e-4 5e-2])
% 
% figure
% cmap = colormap(jet(length(hcf_vals)));
% subplot(2,1,1)
% for i = 1:length(hcf_vals)
%     ebar_mat(dur_hist_bincenters,up_dur_pmf(:,i,:),.001,cmap(i,:)), hold on
% end
% xlim([0 3])
% set(gca,'yscale','log')
% ylim([2e-4 5e-2])
% subplot(2,1,2)
% for i = 1:length(hcf_vals)
%     ebar_mat(dur_hist_bincenters,down_dur_pmf(:,i,:),.001,cmap(i,:)), hold on
% end
% xlim([0 3])
% set(gca,'yscale','log')
% ylim([2e-4 5e-2])
% 
% figure
% errorbar(smooth_vals,mean(min_up_dur_hf),std(min_up_dur_hf)/sqrt(n)), hold on
% errorbar(smooth_vals,mean(min_down_dur_hf),std(min_down_dur_hf)/sqrt(n),'r'), hold on
% 
% figure
% errorbar(hcf_vals,mean(min_up_dur),std(min_up_dur)/sqrt(n)), hold on
% errorbar(hcf_vals,mean(min_down_dur),std(min_down_dur)/sqrt(n),'r'), hold on
% 
% figure
% errorbar(smooth_vals,mean(mean_up_dur_hf),std(mean_up_dur_hf)/sqrt(n)), hold on
% errorbar(smooth_vals,mean(mean_down_dur_hf),std(mean_down_dur_hf)/sqrt(n),'r'), hold on
% 
% figure
% errorbar(hcf_vals,mean(mean_up_dur),std(mean_up_dur)/sqrt(n)), hold on
% errorbar(hcf_vals,mean(mean_down_dur),std(mean_down_dur)/sqrt(n),'r'), hold on
% 
% figure
% errorbar(smooth_vals,mean(ham_dist_mp_lf8hf),std(ham_dist_mp_lf8hf)/sqrt(n)), hold on
% figure
% errorbar(hcf_vals,mean(ham_dist_mp_lf8),std(ham_dist_mp_lf8)/sqrt(n),'r'), hold on
% 
% figure
% errorbar(smooth_vals,mean(lf8_kl_hf),std(lf8_kl_hf)/sqrt(n)), hold on
% figure
% errorbar(hcf_vals,mean(lf8_kl_lf),std(lf8_kl_lf)/sqrt(n),'r'), hold on
% 
% 
