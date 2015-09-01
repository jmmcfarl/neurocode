%inputs:
% raw_data
% raw_Fs
clear all
close all

drive_letter = 'G';

addpath('G:\WC_Germany\hmm_state_detect\\')
addpath('G:\WC_Germany\hsmm_state_detection')
addpath('G:\WC_Germany\parietal_cortical_2010')
% addpath('G:\WC_Germany\new_stellate_analysis\')
addpath('G:\WC_Germany\persistent_2010\')
addpath('G:\WC_Germany\hsmm_uds_code\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load ./desynch_times_individual

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
Fs = 50.4;

% parietal = find_struct_field_vals(sess_data,'region','parietal');
% sess_data = sess_data(parietal);
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times_mp(interneurons) = [];
desynch_times_lf8(interneurons) = [];
desynch_times_lf4(interneurons) = [];

%restrict analysis to the datasets with a frontal LFP
thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
sess_data = sess_data(thom_el);
desynch_times_mp = desynch_times_mp(thom_el);
desynch_times_lf4 = desynch_times_lf4(thom_el);
desynch_times_lf8 = desynch_times_lf8(thom_el);

%define the range for estimating state duration distributions
min_state_dur = 1;
max_state_dur = round(Fs*30);
dur_hist_bincenters = (min_state_dur:max_state_dur)/Fs;

hcf_vals = [0.5 1 2 4 10]; %high-cutoff freqs for computing LF amplitude
smooth_vals = [0.025 0.05 0.1 0.2 0.3]; %HF power smoothing sigmas

n = length(sess_data);

%% For LF4 LF AMPLITUDE
for d = 1:n
    d
    direct = sess_data(d).directory;
    direct(1) = 'G';
    cd(direct)
    load ./used_data lf4
    %cycle over all values of the high-cutoff frequency
    for i = 1:length(hcf_vals)
        %fit a standard HMM and get the realigned ML state sequence, along
        %with the set of state durations
        [hmm_bbstate_seq4,hmm4,state_durations4{d,i}] = parietal_get_hmm_lf_test(lf4,raw_Fs,hcf_vals(i),desynch_times_lf4{d});
        %compute empirical pmfs for up and down state durations
        up_dur_pmf4(d,i,:) = hist(state_durations4{d,i}{2},dur_hist_bincenters); 
        up_dur_pmf4(d,i,:) = up_dur_pmf4(d,i,:)/sum(up_dur_pmf4(d,i,:));
        down_dur_pmf4(d,i,:) = hist(state_durations4{d,i}{1},dur_hist_bincenters);
        down_dur_pmf4(d,i,:) = down_dur_pmf4(d,i,:)/sum(down_dur_pmf4(d,i,:));
        %statistics of state durations
        min_up_dur4(d,i) = min(state_durations4{d,i}{2});
        min_down_dur4(d,i) = min(state_durations4{d,i}{1});
        mean_up_dur4(d,i) = mean(state_durations4{d,i}{2});
        mean_down_dur4(d,i) = mean(state_durations4{d,i}{1});
        
        %compute the separability (KL-Divergence) of the state-conditional
        %distributions
        meandiff = []; %vector containing the difference in state-conditional means as a function of time
        for s = 1:hmm4.Nsegs
            meandiff = [meandiff; hmm4.state(2).meanfun{s}-hmm4.state(1).meanfun{s}];
        end
%         meandiff = mean(meandiff);
        covar1 = hmm4.state(1).var; %down state variance
        covar2 = hmm4.state(2).var; %up state variance
        kl1 = gauss_kl_div(meandiff',covar1,covar2);
        kl2 = gauss_kl_div(-meandiff',covar2,covar1);
        lf4_kl_lf(d,i) = kl1+kl2;
        
% %         load hsmm_state_seq_seg_lf_4_28_10_v3
% %         ham_dist_mp_lf4(d,i) = compute_state_seq_seg_hamdist_varuds(hsmm_bbstate_seq,hmm_bbstate_seq4,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
        clear hmm*
    end
end

% %% For MP LF AMPLITUDE
% for d = 1:n
%     d
%     direct = sess_data(d).directory;
%     direct(1) = 'G';
%     cd(direct)
%     load ./used_data wcv_minus_spike
%     for i = 1:length(hcf_vals)
%         [hmm_bbstate_seq,hmm,state_durations{d,i}] = parietal_get_hmm_lf_test(wcv_minus_spike,raw_Fs,hcf_vals(i),desynch_times_mp{d});
%         up_dur_pmf(d,i,:) = hist(state_durations{d,i}{2},dur_hist_bincenters);
%         up_dur_pmf(d,i,:) = up_dur_pmf(d,i,:)/sum(up_dur_pmf(d,i,:));
%         down_dur_pmf(d,i,:) = hist(state_durations{d,i}{1},dur_hist_bincenters);
%         down_dur_pmf(d,i,:) = down_dur_pmf(d,i,:)/sum(down_dur_pmf(d,i,:));
%         min_up_dur(d,i) = min(state_durations{d,i}{2});
%         min_down_dur(d,i) = min(state_durations{d,i}{1});
%         mean_up_dur(d,i) = mean(state_durations{d,i}{2});
%         mean_down_dur(d,i) = mean(state_durations{d,i}{1});
% 
% %         meandiff = [];
% %         for s = 1:hmm.Nsegs
% %             meandiff = [meandiff; hmm.state(2).meanfun{s}-hmm.state(1).meanfun{s}];
% %         end
% % %         meandiff = mean(meandiff);
% %         covar1 = hmm.state(1).var;
% %         covar2 = hmm.state(2).var;
% %         mp_kl_lf(d,i) = gauss_kl_div(meandiff,covar1,covar2);
% % %         load hsmm_state_seq_seg_lf_4_28_10_v3
% % %         ham_dist_mp_lf4(d,i) = compute_state_seq_seg_hamdist_varuds(hsmm_bbstate_seq,hmm_bbstate_seq4,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
%         clear hmm*
%     end
% end

%% For LF4 HF POWER
for d = 1:n
    d
    direct = sess_data(d).directory;
    direct(1) = 'G';
    cd(direct)
    load used_data lf4
    for i = 1:length(smooth_vals)
        [hmm_bbstate_seq4,hmm4_hf,state_durations_hf4{d,i}] = parietal_get_hmm_hf_test(lf4,raw_Fs,smooth_vals(i),desynch_times_lf4{d});
        up_dur_hf_pmf4(d,i,:) = hist(state_durations_hf4{d,i}{2},dur_hist_bincenters);
        up_dur_hf_pmf4(d,i,:) = up_dur_hf_pmf4(d,i,:)/sum(up_dur_hf_pmf4(d,i,:));
        down_dur_hf_pmf4(d,i,:) = hist(state_durations_hf4{d,i}{1},dur_hist_bincenters);
        down_dur_hf_pmf4(d,i,:) = down_dur_hf_pmf4(d,i,:)/sum(down_dur_hf_pmf4(d,i,:));
        min_up_dur_hf4(d,i) = min(state_durations_hf4{d,i}{2});
        min_down_dur_hf4(d,i) = min(state_durations_hf4{d,i}{1});
        mean_up_dur_hf4(d,i) = mean(state_durations_hf4{d,i}{2});
        mean_down_dur_hf4(d,i) = mean(state_durations_hf4{d,i}{1});

        meandiff = [];
        for s = 1:hmm4_hf.Nsegs
            meandiff = [meandiff; hmm4_hf.state(2).meanfun{s}-hmm4_hf.state(1).meanfun{s}];
        end
%         meandiff = mean(meandiff);
        covar1 = hmm4_hf.state(1).var;
        covar2 = hmm4_hf.state(2).var;
        kl1 = gauss_kl_div(meandiff',covar1,covar2);
        kl2 = gauss_kl_div(-meandiff',covar2,covar1);
        lf4_kl_hf(d,i) = kl1+kl2;
        
% %         load hsmm_state_seq_seg_lf_4_28_10_v3
% %         ham_dist_mp_lf4hf(d,i) = compute_state_seq_seg_hamdist_varuds(hsmm_bbstate_seq,hmm_bbstate_seq4,hmm.UDS_segs,hmm4_hf.UDS_segs,Fsd);
        clear hmm*
    end
end

% %% For MP HF POWER
% for d = 1:n
%     d
%     direct = sess_data(d).directory;
%     direct(1) = 'G';
%     cd(direct)
%     load used_data wcv_minus_spike
%     for i = 1:length(smooth_vals)
%         [hmm_bbstate_seq,hmm_hf,state_durations_hf{d,i}] = parietal_get_hmm_hf_test(wcv_minus_spike,raw_Fs,smooth_vals(i),desynch_times_mp{d});
%         up_dur_hf_pmf(d,i,:) = hist(state_durations_hf{d,i}{2},dur_hist_bincenters);
%         up_dur_hf_pmf(d,i,:) = up_dur_hf_pmf(d,i,:)/sum(up_dur_hf_pmf(d,i,:));
%         down_dur_hf_pmf(d,i,:) = hist(state_durations_hf{d,i}{1},dur_hist_bincenters);
%         down_dur_hf_pmf(d,i,:) = down_dur_hf_pmf(d,i,:)/sum(down_dur_hf_pmf(d,i,:));
%         min_up_dur_hf(d,i) = min(state_durations_hf{d,i}{2});
%         min_down_dur_hf(d,i) = min(state_durations_hf{d,i}{1});
%         mean_up_dur_hf(d,i) = mean(state_durations_hf{d,i}{2});
%         mean_down_dur_hf(d,i) = mean(state_durations_hf{d,i}{1});
% 
% %         meandiff = [];
% %         for s = 1:hmm_hf.Nsegs
% %             meandiff = [meandiff; hmm_hf.state(2).meanfun{s}-hmm_hf.state(1).meanfun{s}];
% %         end
% % %         meandiff = mean(meandiff);
% %         covar1 = hmm_hf.state(1).var;
% %         covar2 = hmm_hf.state(2).var;
% %         mp_kl_hf(d,i) = gauss_kl_div(meandiff,covar1,covar2);
% % %         load hsmm_state_seq_seg_lf_4_28_10_v3
% % %         ham_dist_mp_lf4hf(d,i) = compute_state_seq_seg_hamdist_varuds(hsmm_bbstate_seq,hmm_bbstate_seq4,hmm.UDS_segs,hmm4_hf.UDS_segs,Fsd);
%         clear hmm*
%     end
% end
% 
%%
cd G:\WC_Germany\parietal_cortical_2010
save smoothness_test_results_lf4_4_6_2011

%% FOR LF4
figure
cmap = colormap(jet(length(smooth_vals)));
subplot(2,1,1)
for i = 1:length(smooth_vals)
    ebar_mat(dur_hist_bincenters,up_dur_hf_pmf4(:,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 8e-2])
subplot(2,1,2)
for i = 1:length(smooth_vals)
    ebar_mat(dur_hist_bincenters,down_dur_hf_pmf4(:,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 8e-2])

figure
cmap = colormap(jet(length(hcf_vals)));
subplot(2,1,1)
for i = 1:length(hcf_vals)
    ebar_mat(dur_hist_bincenters,up_dur_pmf4(:,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 5e-2])
subplot(2,1,2)
for i = 1:length(hcf_vals)
    ebar_mat(dur_hist_bincenters,down_dur_pmf4(:,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 5e-2])

figure
errorbar(smooth_vals,mean(min_up_dur_hf4),std(min_up_dur_hf4)/sqrt(n)), hold on
errorbar(smooth_vals,mean(min_down_dur_hf4),std(min_down_dur_hf4)/sqrt(n),'r'), hold on

figure
errorbar(hcf_vals,mean(min_up_dur4),std(min_up_dur4)/sqrt(n)), hold on
errorbar(hcf_vals,mean(min_down_dur4),std(min_down_dur4)/sqrt(n),'r'), hold on

% figure
% errorbar(smooth_vals,mean(mean_up_dur_hf4),std(mean_up_dur_hf4)/sqrt(n)), hold on
% errorbar(smooth_vals,mean(mean_down_dur_hf4),std(mean_down_dur_hf4)/sqrt(n),'r'), hold on
% 
% figure
% errorbar(hcf_vals,mean(mean_up_dur4),std(mean_up_dur4)/sqrt(n)), hold on
% errorbar(hcf_vals,mean(mean_down_dur4),std(mean_down_dur4)/sqrt(n),'r'), hold on

figure
errorbar(smooth_vals,mean(lf4_kl_hf),std(lf4_kl_hf)/sqrt(n)), hold on
figure
errorbar(hcf_vals,mean(lf4_kl_lf),std(lf4_kl_lf)/sqrt(n),'r'), hold on

%% FOR MP
figure
cmap = colormap(jet(length(smooth_vals)));
subplot(2,1,1)
for i = 1:length(smooth_vals)
    ebar_mat(dur_hist_bincenters,up_dur_hf_pmf(:,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 8e-2])
subplot(2,1,2)
for i = 1:length(smooth_vals)
    ebar_mat(dur_hist_bincenters,down_dur_hf_pmf(:,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 8e-2])

figure
cmap = colormap(jet(length(hcf_vals)));
subplot(2,1,1)
for i = 1:length(hcf_vals)
    ebar_mat(dur_hist_bincenters,up_dur_pmf(:,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 5e-2])
subplot(2,1,2)
for i = 1:length(hcf_vals)
    ebar_mat(dur_hist_bincenters,down_dur_pmf(:,i,:),.001,cmap(i,:)), hold on
end
xlim([0 3])
set(gca,'yscale','log')
ylim([2e-4 5e-2])

figure
errorbar(smooth_vals,mean(min_up_dur_hf),std(min_up_dur_hf)/sqrt(n)), hold on
errorbar(smooth_vals,mean(min_down_dur_hf),std(min_down_dur_hf)/sqrt(n),'r'), hold on

figure
errorbar(hcf_vals,mean(min_up_dur),std(min_up_dur)/sqrt(n)), hold on
errorbar(hcf_vals,mean(min_down_dur),std(min_down_dur)/sqrt(n),'r'), hold on

figure
errorbar(smooth_vals,mean(mean_up_dur_hf),std(mean_up_dur_hf)/sqrt(n)), hold on
errorbar(smooth_vals,mean(mean_down_dur_hf),std(mean_down_dur_hf)/sqrt(n),'r'), hold on

figure
errorbar(hcf_vals,mean(mean_up_dur),std(mean_up_dur)/sqrt(n)), hold on
errorbar(hcf_vals,mean(mean_down_dur),std(mean_down_dur)/sqrt(n),'r'), hold on

figure
errorbar(smooth_vals,mean(mp_kl_hf),std(mp_kl_hf)/sqrt(n)), hold on
figure
errorbar(hcf_vals,mean(mp_kl_lf),std(mp_kl_lf)/sqrt(n),'r'), hold on

