clear all
close all

addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')
addpath('G:\WC_Germany\hsmm_uds_code\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_individual

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');
thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
thom_par = thom_el(find(ismember(thom_el,parietal)));
thom_pfc = setdiff(thom_el,thom_par);
parietal = find_struct_field_vals(sess_data,'region','parietal');
frontal = setdiff(1:length(sess_data),parietal);
superficial = find_struct_field_vals(sess_data,'layer','23');
deep = setdiff(1:length(sess_data),superficial);

sess_data = sess_data(parietal);

n = length(sess_data);
% d=2
for d = 1:n
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    cd(cdir)
    pwd
    load ./used_data lf8 wcv_minus_spike
    
    %% FIRST LOCATE DESYNCHRONIZED EPOCHS
    SO_thresh = -6;
    HF_thresh = -2;
    Fs_orig = 2016;
    params.dsf = 8;
    params.lcf = 0.1;
    params.hcf = 40;
    params.movingwin = [15 2.5];
    [desynch_times,desynch_ids,P,f,t] = hsmm_uds_locate_desynch_times(wcv_minus_spike,Fs_orig,params,SO_thresh,HF_thresh);
    [desynch_times4,desynch_ids4,P,f,t] = hsmm_uds_locate_desynch_times(lf8,Fs_orig,params,SO_thresh,HF_thresh);
    
    %% COMPUTE SIGNAL FEATURES
    hmm_dsf = 40;
    hmm_Fs = Fs_orig/hmm_dsf;
    lf_cut_off_freqs = [0.05 2];
    hf_cut_off_freqs = [20 80];
    smooth_sigma = 0.15; %in seconds
    [lf_features,t_axis,Fs] = hsmm_uds_get_lf_features(wcv_minus_spike,Fs_orig,hmm_dsf,lf_cut_off_freqs); %'low-frequency amplitude'
    [lf_features4,t_axis,Fs] = hsmm_uds_get_lf_features(lf8,Fs_orig,hmm_dsf,lf_cut_off_freqs); %'low-frequency amplitude'
    
    %% LOCATE THE SEGMENTS CONTAINING UDS
    T = length(lf_features);
    min_seg_dur = 60;
    UDS_segs = hsmm_uds_get_uds_segments(desynch_times,hmm_Fs,T,min_seg_dur);
    UDS_segs4 = hsmm_uds_get_uds_segments(desynch_times4,hmm_Fs,T,min_seg_dur);
    
    %% INITIALIZE THE HMM
    params.meantype = 'variable';
    params.UDS_segs = UDS_segs;
    params.movingwin = [50 10];
    params4.meantype = 'variable';
    params4.UDS_segs = UDS_segs4;
    params4.movingwin = [50 10];
    [hmm] = hsmm_uds_initialize(lf_features,hmm_Fs,params);
    [hmm4] = hsmm_uds_initialize(lf_features4,hmm_Fs,params4);
    
    %% FIT AN HMM
    hmm.min_mllik = -8; %set minimum value of the maximal log likelihood for either state before switching to robust estimator.  If you don't want to use this feature, set this to -Inf
    hmm.windowSize = 50;
%     hmm.windowSlide = 1;
    hmm4.min_mllik = -8; %set minimum value of the maximal log likelihood for either state before switching to robust estimator.  If you don't want to use this feature, set this to -Inf
    hmm4.windowSize = 50;
%     hmm4.windowSlide = 1;
    hmm = hsmm_uds_train_hmm(hmm,lf_features);
    hmm4 = hsmm_uds_train_hmm(hmm4,lf_features4);
    
    %% DETERMINE THE VITERBI SEQUENCE
    [hmm_state_seq,llik_best] = hsmm_uds_viterbi_hmm(hmm,lf_features);
    [hmm_state_seq4,llik_best] = hsmm_uds_viterbi_hmm(hmm4,lf_features4);

    dsf_bb = 8;
    Fs_bb = 2016/dsf_bb;
    bb_lcf = 0.05;
    bb_hcf = 20;
    [bb_features,t_axis,Fsd] = hsmm_uds_get_lf_features(wcv_minus_spike,Fs_orig,dsf_bb,[bb_lcf bb_hcf]); %'low-frequency amplitude'
    [bb_features4,t_axis,Fsd] = hsmm_uds_get_lf_features(lf8,Fs_orig,dsf_bb,[bb_lcf bb_hcf]); %'low-frequency amplitude'

    min_state_dur = 1;
    max_state_dur = round(hmm_Fs*30);
    dur_range = (1:max_state_dur)/hmm_Fs;
hmm_temp = hmm;
temp_p = zeros(2,length(dur_range));
temp_p(1,:) = geometric_pmf(1:length(dur_range),hmm_temp.A(2,1));
temp_p(2,:) = geometric_pmf(1:length(dur_range),hmm_temp.A(1,2));
hmm_temp.P = temp_p;
for i = 1:hmm.K
    hmm.state(i).dur_type = 'hmm';
    hmm4.state(i).dur_type = 'hmm';
end
    hmm.dur_range = dur_range;
    hmm.max_state_dur = max_state_dur;
    hmm.min_state_dur = min_state_dur;
    hmm4.dur_range = dur_range;
    hmm4.max_state_dur = max_state_dur;
    hmm4.min_state_dur = min_state_dur;

    pert_range_bb = [-0.15 0.15];
    [hmm_bb_state_seq] = hsmm_pert_optimize_transitions(bb_features,hmm,hmm_state_seq,hmm.gamma,hmm.Fs,Fs_bb,pert_range_bb);
    [hmm_bb_state_seq4] = hsmm_pert_optimize_transitions(bb_features4,hmm4,hmm_state_seq4,hmm4.gamma,hmm4.Fs,Fs_bb,pert_range_bb);
    
    %%
    smoothed_seq = hsmm_uds_smooth_stateseq(hmm_state_seq,hmm_Fs,100,100);
    state_durations = hsmm_uds_compute_state_durations(smoothed_seq,hmm_Fs);

    smoothed_seq4 = hsmm_uds_smooth_stateseq(hmm_state_seq4,hmm_Fs,100,100);
    state_durations4 = hsmm_uds_compute_state_durations(smoothed_seq4,hmm_Fs);

    %%
    
    for i = 1:hmm.K
        %estimate the empirical state duration pmf
        emp_pmf(i,:) = hist(state_durations{i},dur_range);
        emp_pmf(i,:) = emp_pmf(i,:)/sum(emp_pmf(i,:));
        emp_pmf4(i,:) = hist(state_durations4{i},dur_range);
        emp_pmf4(i,:) = emp_pmf4(i,:)/sum(emp_pmf4(i,:));
        
        p_fit = 1/nanmean(state_durations{i}*hmm_Fs);
        geo_pmf(i,:) = geometric_pmf(1:length(dur_range),p_fit);
        [mu,lambda] = inverse_gauss_mlfit(state_durations{i});%compute the ML parameters of the inverse gaussian dist
        ig_pmf(i,:) = inverse_gauss_pmf(dur_range,mu,lambda);%estimate the IG PMF
        [alpha,beta] = gamma_mlfit(state_durations{i}); %compute the ML parameters of the gamma
        gam_pmf(i,:) = gamma_pmf(dur_range,alpha,beta); %estimate the gamma PMF

        p_fit4 = 1/nanmean(state_durations4{i}*hmm_Fs);
        geo_pmf4(i,:) = geometric_pmf(1:length(dur_range),p_fit4);
        [mu4,lambda4] = inverse_gauss_mlfit(state_durations4{i});%compute the ML parameters of the inverse gaussian dist
        ig_pmf4(i,:) = inverse_gauss_pmf(dur_range,mu4,lambda4);%estimate the IG PMF
        [alpha4,beta4] = gamma_mlfit(state_durations{i}); %compute the ML parameters of the gamma
        gam_pmf4(i,:) = gamma_pmf(dur_range,alpha4,beta4); %estimate the gamma PMF
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
    hmm4.state(1).dur_type = 'inv_gauss';
    hmm4.state(2).dur_type = 'inv_gauss';
    %initialize state duration model parameters
    for i = 1:hmm.K
        if strcmp(hmm.state(i).dur_type,'inv_gauss')
            hmm.state(i).dur_pars = [mu lambda];
            hmm.P(i,:) = ig_pmf(i,:);
            hmm4.state(i).dur_pars = [mu4 lambda4];
            hmm4.P(i,:) = ig_pmf4(i,:);
        else
            hmm.state(i).dur_pars = [alpha beta];
            hmm.P(i,:) = gam_pmf(i,:);
            hmm4.state(i).dur_pars = [alpha4 beta4];
            hmm4.P(i,:) = gam_pmf4(i,:);
        end
    end
    
    % enforce any minimum state duration and renormalize
    for i = 1:hmm.K
        if min_state_dur > 1
            hmm.P(i,1:min_state_dur-1) = zeros(1,min_state_dur-1);
            hmm.P(i,:) = hmm.P(i,:)/sum(hmm.P(i,:));
            hmm4.P(i,1:min_state_dur-1) = zeros(1,min_state_dur-1);
            hmm4.P(i,:) = hmm4.P(i,:)/sum(hmm4.P(i,:));
        end
    end
    
    %%
    [hmm,gamma,hmm_window_post]=hsmm_uds_train_hsmm(hmm,lf_features);
    [hmm4,gamma4,hmm_window_post4]=hsmm_uds_train_hsmm(hmm4,lf_features4);
    % [hmm,gamma,hmm_window_post,bad_model]=hsmm_uds_train_seg(hmm,lf_features);
    %%
    [hsmm_state_seq,max_lik] = hsmm_uds_viterbi_hsmm(hmm,lf_features);
    [hsmm_state_seq4,max_lik] = hsmm_uds_viterbi_hsmm(hmm4,lf_features4);
    
    %%
    [hsmm_bb_state_seq] = hsmm_uds_pert_optimize_transitions(bb_features,hmm,hsmm_state_seq,gamma,hmm.Fs,Fs_bb,pert_range_bb);
    [hsmm_bb_state_seq4] = hsmm_uds_pert_optimize_transitions(bb_features4,hmm4,hsmm_state_seq4,gamma4,hmm4.Fs,Fs_bb,pert_range_bb);
    
    hmm_mp_lf8_ham(d) = compute_state_seq_seg_hamdist_varuds(hmm_bb_state_seq,hmm_bb_state_seq4,hmm.UDS_segs,hmm4.UDS_segs,Fs_bb);
    hsmm_mp_lf8_ham(d) = compute_state_seq_seg_hamdist_varuds(hsmm_bb_state_seq,hsmm_bb_state_seq4,hmm.UDS_segs,hmm4.UDS_segs,Fs_bb);
    
    %%
    [fixmean_sm_state_seq,fixmean_state_seq,Fs,fhmm] = parietal_get_fixthresh_uds_state_seq...
        (wcv_minus_spike,2016,desynch_times);
    [fixmean_sm_state_seq4,fixmean_state_seq4,Fs,fhmm] = parietal_get_fixthresh_uds_state_seq...
        (lf8,2016,desynch_times4);
    [fixmean_sm_state_seqz,fixmean_state_seqz,Fs,fhmm] = parietal_get_fixthresh_zm_uds_state_seq...
        (wcv_minus_spike,2016,desynch_times);
    [fixmean_sm_state_seqz4,fixmean_state_seqz4,Fs,fhmm] = parietal_get_fixthresh_zm_uds_state_seq...
        (lf8,2016,desynch_times4);
    fm_mp_lf8_ham(d) = compute_state_seq_seg_hamdist_varuds(fixmean_state_seq,fixmean_state_seq4,hmm.UDS_segs,hmm4.UDS_segs,Fs);
    fmz_mp_lf8_ham(d) = compute_state_seq_seg_hamdist_varuds(fixmean_state_seqz,fixmean_state_seqz4,hmm.UDS_segs,hmm4.UDS_segs,Fs);
    
    fprintf('Recording: %d\n',d)
    fprintf('HMM agreement: %.6f\n',hmm_mp_lf8_ham(d))
    fprintf('HSMM agreement: %.6f\n',hsmm_mp_lf8_ham(d))
    fprintf('FM agreement: %.6f\n',fm_mp_lf8_ham(d))
    fprintf('FMz agreement: %.6f\n',fmz_mp_lf8_ham(d))
end