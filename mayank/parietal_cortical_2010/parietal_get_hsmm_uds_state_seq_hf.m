function [bb_hsmm_state_seq,hsmm_state_seq,hmm_state_seq,hmm,Fs] = parietal_get_hsmm_uds_state_seq_hf...
    (raw_data,raw_Fs,desynch_times,f_names)

drive_letter = 'G';
addpath(strcat(drive_letter,':\WC_Germany\hsmm_state_detection'))

Fs_desired = 50; %(in Hz)
dsf = round(raw_Fs/Fs_desired);

Fs = raw_Fs/dsf;
niqf = raw_Fs/2;
hf_lcf = 20; %low cut-off for the low-freq filter
hf_hcf = 80; %high cut-off for the low-freq filter
hf_smooth = 0.15;
num_states = 2; %number of hidden states (must be 2 in current implementation)

meanparams.windowSize = 50; %window duration to estimate time varying state means (s)
meanparams.windowSlide = 1; %window spacing to estimate time varying state means (s)
meanparams.meantype = 'variable';

min_state_dur = 1;
max_state_dur = round(Fs*30);
min_seg_dur = 60;
dur_range = (1:max_state_dur)/Fs;
% dur_hist_range = linspace(0,10,100);

%% initializations
[obs_hf,t_axis,Fs] = get_hf_features(raw_data,raw_Fs,Fs_desired,[hf_lcf hf_hcf],hf_smooth);

%% compute UDS segments
N = length(obs_hf);
desynch_indices = round(desynch_times*Fs);
% if size(desynch_indices,1) ~= 0
%     if size(desynch_indices,1)==1
%         if desynch_indices(1,2) < N
%             UDS_segs = [desynch_indices(1,2) N];
%         else
%             UDS_segs = [];
%         end
%     elseif size(desynch_indices,1) > 1
%         UDS_segs = [];
%         for i = 1:size(desynch_indices,1)-1
%             UDS_segs = [UDS_segs; desynch_indices(i,2) desynch_indices(i+1,1)];
%         end
%     end
%     if desynch_indices(1,1) ~= 1
%         UDS_segs = [1 desynch_indices(1,1); UDS_segs];
%     end
% else
%     UDS_segs = [1 N];
% end
temp_inds = zeros(N,1);
for i = 1:size(desynch_indices,1)
    temp_inds(desynch_indices(i,1):desynch_indices(i,2)) = 1;
end
UDS_start = find(temp_inds(1:end-1)==1 & temp_inds(2:end)==0);
UDS_stop = find(temp_inds(1:end-1)==0 & temp_inds(2:end)==1);
if temp_inds(1)==0
    UDS_start = [1; UDS_start];
end
if temp_inds(end)==0
    UDS_stop = [UDS_stop; N];
end
UDS_segs = [UDS_start UDS_stop];

UDS_seg_durs = (UDS_segs(:,2)-UDS_segs(:,1))/Fs;
too_short = find(UDS_seg_durs < min_seg_dur);
UDS_segs(too_short,:) = [];

%% initialize hsmm
emiss = obs_hf;
[hmm] = initialize_uds_hsmm_seg(emiss,Fs,meanparams,UDS_segs,drive_letter);

hmm.min_mllik = -1.5; %set minimum value of the maximal log likelihood for either state before switching to robust estimator

[hmm, gamma_hmm] = hmm_uds_train_seg(hmm,emiss);

[hmm_state_seq,llik_best] = hmm_uds_viterbi_seg(hmm,emiss);

smoothed_seq = thresh_state_smooth_seg(hmm_state_seq,Fs,100,100);
[state_durations] = compute_state_durations_seg(smoothed_seq,Fs);
% for i = 1:hmm.K
%    if strcmp(hmm.state(i).dur_type,'inv_gauss')
%        [params(1),params(2)] = inverse_gauss_mlfit(state_durations{i});
%        hmm.state(i).dur_pars = params;
%        hmm.P(i,:) = inverse_gaussian_pmf(dur_range,hmm.state(i).dur_pars(1),hmm.state(i).dur_pars(2));
%    elseif strcmp(hmm.state(i).dur_type,'gamma')
%        [params(1),params(2)] = gamma_mlfit(state_durations{i});
%        hmm.state(i).dur_pars = params;
%        hmm.P(i,:) = gamma_pmf(dur_range,hmm.state(i).dur_pars(1),hmm.state(i).dur_pars(2));
%    else
%        error('invalid duration distribution')
%    end
% end
for i = 1:hmm.K
    emp_pmf = hist(state_durations{i},dur_range);
    emp_pmf = emp_pmf/sum(emp_pmf);
    [mu,lambda] = inverse_gauss_mlfit(state_durations{i});
    ig_pmf = inverse_gaussian_pmf(dur_range,mu,lambda);
    [alpha,beta] = gamma_mlfit(state_durations{i});
    gam_pmf = gamma_pmf(dur_range,alpha,beta);
    ks_ig = max(abs(cumsum(emp_pmf) - cumsum(ig_pmf)));
    ks_gam = max(abs(cumsum(emp_pmf) - cumsum(gam_pmf)));
%     if ks_ig > ks_gam
%         hmm.state(i).dur_type = 'gamma';
%         hmm.state(i).dur_pars = [alpha beta];
%         hmm.P(i,:) = gam_pmf;
%     else
        hmm.state(i).dur_type = 'inv_gauss';
        hmm.state(i).dur_pars = [mu lambda];
        hmm.P(i,:) = ig_pmf;       
%     end
end

%%

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
    
[hmm,gamma,hmm_window_post]=hsmm_uds_train_seg(hmm,emiss);

[hsmm_state_seq,max_lik] = hsmm_uds_viterbi_seg(hmm,emiss);

hmm.posteriors = hmm_window_post;
hmm.max_lik = max_lik;
hmm.gamma = gamma;

%%
dsf_bb = 8;
Fs_bb = raw_Fs/dsf_bb;
hf_smooth = 0.05;
[obs_bb,t_axis,Fs_hf] = get_hf_features(raw_data,raw_Fs,Fs_bb,[hf_lcf hf_hcf],hf_smooth);

pert_range_bb = [-0.15 0.15];
[bb_hsmm_state_seq] = pert_optimize_hsmm_transitions_seg_dp_robust... 
    (obs_bb,hmm,hsmm_state_seq,gamma,Fs,Fs_bb,pert_range_bb);

%%
% [obs_dist,obsrange,win_t] = compute_slide_dist(emiss,meanparams.windowSize,meanparams.windowSlide,Fs);
% figure
% pcolor(win_t,obsrange,log(obs_dist'));shading flat
% caxis([-4 0])
% hold on
% for i = 1:hmm.Nsegs
%     t_axis = (hmm.UDS_segs(i,1):hmm.UDS_segs(i,2))/hmm.Fs;
%     plot(t_axis,hmm.state(1).meanfun{i},'w')
%     plot(t_axis,hmm.state(2).meanfun{i},'w')
% end
% t_names = ['F:\WC_Germany\parietal_cortical_2010\state_meanfuns\hf_' f_names];
% print('-dpng',t_names), close

