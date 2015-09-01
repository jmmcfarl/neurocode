
addpath('G:\WC_Germany\hsmm_uds_code\')

% cd G:\wc_data\2011-04-01_EEG
% load ./2011-4-1_10-11-24_all_eeg_data
cd G:\wc_data\2011-04-07-EEG
load ./2011-04-07_B_all_eeg_data
% cd G:\wc_data\2011-04-08_A_EEG
% load ./2011-04-08_A_CSC1
% cd G:\wc_data\2011-04-08_B_EEG
% load ./2011-04-08_B_CSC1

%% for 04-07-2011
% xl1 = [863.5 894.5];
xl1 = [867 894.5];
% xl4 = [560 580]; %?
% x3 = [698 740.5];
xl2 = [762 796];
% x5 = [951 986];

xl3 = [216 232];
xl4 = [450 466];
xl5 = [719 756];
xl6 = [721 790];
xl7 = [773 804];

%%
% lf8 = CSC8_Samples(:);
% lf7 = CSC7_Samples(:);
% lf6 = CSC6_Samples(:);
% lf5 = CSC5_Samples(:);
% lf4 = CSC4_Samples(:);
eeg = -CSC1_Samples(:);
Fs_orig = 2016;

%% FIRST LOCATE DESYNCHRONIZED EPOCHS
% desynch_times = [];

%% COMPUTE SIGNAL FEATURES
hmm_dsf = 40;
hmm_Fs = Fs_orig/hmm_dsf;
lf_cut_off_freqs = [0.05 2];
bb_cut_off_freqs = [0.05 20];
hf_cut_off_freqs = [20 80];
smooth_sigma = 0.05; %in seconds
% [lf6_lf,t_axis,Fs] = hsmm_uds_get_lf_features(lf6,Fs_orig,hmm_dsf,lf_cut_off_freqs); %'low-frequency amplitude'
% [lf5_lf,t_axis,Fs] = hsmm_uds_get_lf_features(lf5,Fs_orig,hmm_dsf,lf_cut_off_freqs); %'low-frequency amplitude'
% [lf7_lf,t_axis,Fs] = hsmm_uds_get_lf_features(lf7,Fs_orig,hmm_dsf,lf_cut_off_freqs); %'low-frequency amplitude'
% [lf8_lf,t_axis,Fs] = hsmm_uds_get_lf_features(lf8,Fs_orig,hmm_dsf,lf_cut_off_freqs); %'low-frequency amplitude'
% [lf5_hf,t_axis,Fs] = hsmm_uds_get_hf_features(lf5,Fs_orig,hmm_dsf,hf_cut_off_freqs,smooth_sigma); %'low-frequency amplitude'
% [lf5_bb,t_axis,Fs] = hsmm_uds_get_lf_features(lf5,Fs_orig,hmm_dsf,bb_cut_off_freqs); %'low-frequency amplitude'
% [lf4_lf,t_axis,Fs] = hsmm_uds_get_lf_features(lf4,Fs_orig,hmm_dsf,lf_cut_off_freqs); %'low-frequency amplitude'
[eeg_lf,t_axis,Fs] = hsmm_uds_get_lf_features(eeg,Fs_orig,hmm_dsf,lf_cut_off_freqs); %'low-frequency amplitude'
[eeg_bb,t_axis,Fs] = hsmm_uds_get_lf_features(eeg,Fs_orig,hmm_dsf,bb_cut_off_freqs); %'low-frequency amplitude'
[eeg_hf,t_axis,Fs] = hsmm_uds_get_hf_features(eeg,Fs_orig,hmm_dsf,hf_cut_off_freqs,smooth_sigma); %'low-frequency amplitude'

%%
used_sig = eeg_lf;

%%
[desynch_times,desynch_inds,P_mp,f,t] = locate_desynch_times_individual(eeg);

%% LOCATE THE SEGMENTS CONTAINING UDS
T = length(used_sig);
min_seg_dur = 60;
UDS_segs = hsmm_uds_get_uds_segments(desynch_times,hmm_Fs,T,min_seg_dur);

%% INITIALIZE THE HMM
params.meantype = 'variable';
params.UDS_segs = UDS_segs;
params.movingwin = [50 10];
[hmm] = hsmm_uds_initialize(used_sig,hmm_Fs,params);

%% FIT AN HMM
hmm.min_mllik = -8; %set minimum value of the maximal log likelihood for either state before switching to robust estimator.  If you don't want to use this feature, set this to -Inf
hmm.windowSize = 50;
hmm = hsmm_uds_train_hmm(hmm,used_sig);

%% DETERMINE THE VITERBI SEQUENCE
[state_seq,llik_best] = hsmm_uds_viterbi_hmm(hmm,used_sig);

%%
min_state_dur = 1;
max_state_dur = round(hmm_Fs*30);
dur_range = (1:max_state_dur)/hmm_Fs;
hmm.dur_range = dur_range;
hmm.max_state_dur = max_state_dur;
hmm.min_state_dur = min_state_dur;
smoothed_seq = hsmm_uds_smooth_stateseq(state_seq,hmm_Fs,100,100);
state_durations = hsmm_uds_compute_state_durations(smoothed_seq,hmm_Fs);

for i = 1:hmm.K
    %estimate the empirical state duration pmf
    emp_pmf(i,:) = hist(state_durations{i},dur_range);
    emp_pmf(i,:) = emp_pmf(i,:)/sum(emp_pmf(i,:));    
    p_fit = 1/nanmedian(state_durations{i}*hmm_Fs);
    geo_pmf(i,:) = geometric_pmf(1:length(dur_range),p_fit);    
    [mu,lambda] = inverse_gauss_mlfit(state_durations{i});%compute the ML parameters of the inverse gaussian dist
    ig_pmf(i,:) = inverse_gauss_pmf(dur_range,mu,lambda);%estimate the IG PMF
        [alpha,beta] = gamma_mlfit(state_durations{i}); %compute the ML parameters of the gamma
    gam_pmf(i,:) = gamma_pmf(dur_range,alpha,beta); %estimate the gamma PMF
end
% figure
% subplot(2,1,1)
% bar(dur_range,emp_pmf(1,:))
% hold on
% plot(dur_range,ig_pmf(1,:),'k',dur_range,gam_pmf(1,:),'r',dur_range,geo_pmf(i,:),'g','linewidth',2)
% xlim([0 3.5])
% title('DOWN State','fontsize',16)
% xlabel('Duration (s)','fontsize',16)
% ylabel('Relative frequency','fontsize',16)
% legend('Empirical Distribution','Inverse Gaussian','Gamma','Geometric')
% subplot(2,1,2)
% bar(dur_range,emp_pmf(2,:))
% hold on
% plot(dur_range,ig_pmf(2,:),'k',dur_range,gam_pmf(2,:),'r',dur_range,geo_pmf(i,:),'g','linewidth',2)
% xlim([0 3.5])
% title('UP State','fontsize',16)
% xlabel('Duration (s)','fontsize',16)
% ylabel('Relative frequency','fontsize',16)

%%
%select duration model.  Choices are: 'inv_gauss', and 'gamma'
hmm.state(1).dur_type = 'inv_gauss';
hmm.state(2).dur_type = 'inv_gauss';
%initialize state duration model parameters
for i = 1:hmm.K
    if strcmp(hmm.state(i).dur_type,'inv_gauss')
        hmm.state(i).dur_pars = [mu lambda];
        hmm.P(i,:) = ig_pmf(i,:);
    else
        hmm.state(i).dur_pars = [alpha beta];
        hmm.P(i,:) = gam_pmf(i,:);
    end
end
hmm.dur_range = dur_range;
hmm.max_state_dur = max_state_dur;
hmm.min_state_dur = min_state_dur;

% enforce any minimum state duration and renormalize
for i = 1:hmm.K
    if min_state_dur > 1
        hmm.P(i,1:min_state_dur-1) = zeros(1,min_state_dur-1);
        hmm.P(i,:) = hmm.P(i,:)/sum(hmm.P(i,:));
    end
end

%%
[hsmm,gamma,hmm_window_post]=hsmm_uds_train_hsmm(hmm,used_sig);
%%
[hsmm_state_seq,max_lik] = hsmm_uds_viterbi_hsmm(hmm,used_sig);

%%
dsf_bb = 8;
raw_Fs = 2016;
Fs_bb = raw_Fs/dsf_bb;
bb_lcf = 0.05;
bb_hcf = 20;
[obs_bb,t_axis,Fs_bb] = get_lf_features(eeg,raw_Fs,Fs_bb,[bb_lcf bb_hcf]);

pert_range_bb = [-0.15 0.15];
% [bb_hsmm_state_seq] = pert_optimize_hsmm_transitions_seg_dp_robust... 
%     (obs_bb,hsmm,hsmm_state_seq,gamma,Fs,Fs_bb,pert_range_bb);
[hsmm_bbstate_seq] = hsmm_uds_pert_optimize_transitions_v2(obs_bb,hsmm,hsmm_state_seq,hsmm.gamma,hsmm.Fs,Fs_bb,pert_range_bb);

%%
[fm_state_seq,fhmm] = hsmm_uds_get_fixthresh_uds_state_seq(used_sig,hmm.Fs,UDS_segs);
[fm_state_seqz,fhmmz] = hsmm_uds_get_fixthresh_uds_state_seqz(used_sig,hmm.Fs,UDS_segs);

%%
meandiff = [];
for i = 1:hsmm.Nsegs
    meandiff = [meandiff; hsmm.state(2).meanfun{i}-hsmm.state(1).meanfun{i}];
end
%         meandiff = mean(meandiff);
covar1 = hsmm.state(1).var;
covar2 = hsmm.state(2).var;
kl1 = gauss_kl_div(meandiff',covar1,covar2);
kl2 = gauss_kl_div(-meandiff',covar2,covar1);
sig_kl_lf = kl1+kl2;
sig_bd_lf = gauss_bd(meandiff',covar1,covar2);

%%
[obs_dist,obsrange,win_t] = compute_slide_dist(used_sig,50,2.5,hsmm.Fs);
figure
set(gca,'fontname','arial','fontsize',12)
imagesc(win_t,obsrange,log(obs_dist'));shading flat, caxis([-4 -0.4])
caxis([-4 0]),hold on
for i = 1:hsmm.Nsegs
    temp = (hsmm.UDS_segs(i,1):hsmm.UDS_segs(i,2))/hsmm.Fs;
    plot(temp,hsmm.state(1).meanfun{i},'w','linewidth',2)
    plot(temp,hsmm.state(2).meanfun{i},'w','linewidth',2)
end
line([0 win_t(end)],[fhmm.threshold fhmm.threshold],'color','k')
line([0 win_t(end)],[fhmmz.threshold fhmmz.threshold],'color','k')
ylim([-2.2 3.2])
xlabel('Time (s)','fontname','arial','fontsize',14)
ylabel('Amplitude (z)','fontname','arial','fontsize',14)
yl = ylim;
line([xl1(1) xl1(1)],yl,'color','k')
line([xl1(2) xl1(2)],yl,'color','k')
line([xl7(1) xl7(1)],yl,'color','k')
line([xl7(2) xl7(2)],yl,'color','k')
% line([xl2(1) xl2(1)],yl,'color','k')
% line([xl2(2) xl2(2)],yl,'color','k')
% line([xl3(1) xl3(1)],yl,'color','k')
% line([xl3(2) xl3(2)],yl,'color','k')
% line([xl4(1) xl4(1)],yl,'color','k')
% line([xl4(2) xl4(2)],yl,'color','k')
%%
[new_seg_inds] = resample_uds_seg_inds(hsmm.UDS_segs,hsmm.Fs,Fs_bb,length(obs_bb));
figure
set(gca,'fontname','arial','fontsize',12)
plot(t_axis,obs_bb)
% plot(t_axis,eeg_lf)
hold on
% plot(t_axis,lf5_bb+3,'k')
for i = 1:hsmm.Nsegs
    temp = (hsmm.UDS_segs(i,1):hsmm.UDS_segs(i,2))/hsmm.Fs;
%     plot(temp,hsmm_state_seq{i},'r','linewidth',2)
    plot(temp,fm_state_seq{i}-1.5,'k','linewidth',2)
    plot(temp,fm_state_seqz{i}-1.5,'color',[0 0.8 0],'linewidth',2)
    temp2 = (new_seg_inds(i,1):new_seg_inds(i,2))/Fs_bb;
    plot(temp2,hsmm_bbstate_seq{i}-1.5,'r','linewidth',2)
end
xlabel('Time (s)','fontname','arial','fontsize',14)
ylabel('Amplitude (z)','fontname','arial','fontsize',14)
line([0 t_axis(end)],[fhmm.threshold fhmm.threshold],'color','k')
line([0 t_axis(end)],[fhmmz.threshold fhmmz.threshold],'color','k')

%%
used_inds = [];
for i = 1:size(new_seg_inds,1)
    used_inds = [used_inds hsmm.UDS_segs(i,1):hsmm.UDS_segs(i,2)];
end
figure
x_axis = linspace(-2.2,3.2,1000);
% ksdensity(lf8_bb)
ksdensity(used_sig(used_inds),x_axis);
hold on
y1 = fhmm.priors(1)/sqrt(2*pi*fhmm.state_vars(1))*exp(-(x_axis-fhmm.state_means(1)).^2/(2*fhmm.state_vars(1)));
y2 = fhmm.priors(2)/sqrt(2*pi*fhmm.state_vars(2))*exp(-(x_axis-fhmm.state_means(2)).^2/(2*fhmm.state_vars(2)));
plot(x_axis,y1,'r')
plot(x_axis,y2,'g')
ylim([0 0.65])
yl = ylim;
% load fixmean_state_seqs4_zm_4_26_10_v2
line([fhmm.threshold fhmm.threshold],yl,'color','k')
% load fixmean_state_seqs4_4_26_10_v2
line([fhmmz.threshold fhmmz.threshold],yl,'color','k')
xlim([-2.2 3.2])

