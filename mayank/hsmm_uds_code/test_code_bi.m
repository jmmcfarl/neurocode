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

sess_data = sess_data(thom_pfc);

n = length(sess_data);
% for d = 1:n
    d=2; %d=3
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    cd(cdir)
    pwd
    load ./used_data lf4 
    sig = lf4;
    %% FIRST LOCATE DESYNCHRONIZED EPOCHS
    SO_thresh = -6;
    HF_thresh = -2;
    Fs_orig = 2016;
    params.dsf = 8;
    params.lcf = 0.1;
    params.hcf = 40;
    params.movingwin = [15 2.5];
    [desynch_times,desynch_ids,P,f,t] = hsmm_uds_locate_desynch_times(lf4,Fs_orig,params,SO_thresh,HF_thresh);
    
    %% COMPUTE SIGNAL FEATURES
    hmm_dsf = 40;
    hmm_Fs = Fs_orig/hmm_dsf;
    lf_cut_off_freqs = [0.05 2];
    hf_cut_off_freqs = [20 80];
    smooth_sigma = 0.75; %in seconds
    [lf_features,t_axis,Fs] = hsmm_uds_get_lf_features(lf4,Fs_orig,hmm_dsf,lf_cut_off_freqs); %'low-frequency amplitude'
    [hf_features,t_axis,Fs] = hsmm_uds_get_hf_features(lf4,Fs_orig,hmm_dsf,hf_cut_off_freqs,smooth_sigma); %'low-frequency amplitude'
    
    %% LOCATE THE SEGMENTS CONTAINING UDS
    T = length(lf_features);
    min_seg_dur = 60;
    UDS_segs = hsmm_uds_get_uds_segments(desynch_times,hmm_Fs,T,min_seg_dur);
    
    %% INITIALIZE THE HMM
    params.meantype = 'variable';
    params.UDS_segs = UDS_segs;
    params.movingwin = [50 10];
    [hmm] = hsmm_uds_initialize_bi([lf_features hf_features],hmm_Fs,params);
    
    %% FIT AN HMM
    hmm.min_mllik = -8; %set minimum value of the maximal log likelihood for either state before switching to robust estimator.  If you don't want to use this feature, set this to -Inf
    hmm.windowSize = 50;
    hmm = hsmm_uds_train_hmm_bi(hmm,[lf_features hf_features]);
    
    %% DETERMINE THE VITERBI SEQUENCE
    [state_seq,llik_best] = hsmm_uds_viterbi_hmm(hmm,[lf_features hf_features]);

    %%
    bb_lcf = 0.05;
    bb_hcf = 20;
    [bb_features,t_axis,Fsd] = hsmm_uds_get_lf_features(lf4,2016,hmm_dsf,[bb_lcf bb_hcf]); %'low-frequency amplitude'

    %%
    [~,temp_t,emissdist_d1,emissrange_d1] = hsmm_uds_init_state_means(lf_features,[50 5],hmm.Fs);
    [~,temp_t,emissdist_d2,emissrange_d2] = hsmm_uds_init_state_means(hf_features,[50 5],hmm.Fs);

    figure
    subplot(2,1,1)
    pcolor(temp_t,emissrange_d1,log(emissdist_d1'));shading flat
    hold on
    caxis([-4 -0.5])
    for i = 1:hmm.Nsegs
        time{i} = (UDS_segs(i,1):UDS_segs(i,2))/hmm.Fs;
        plot(time{i},hmm.state(1).meanfun{i}(:,1),'w','linewidth',2)
        plot(time{i},hmm.state(2).meanfun{i}(:,1),'w','linewidth',2)
    end
    xlabel('Time (s)','fontsize',16)
    ylabel('Amplitude (z)','fontsize',16)
    title('Sliding window density estimation','fontsize',16)
    subplot(2,1,2)
    pcolor(temp_t,emissrange_d2,log(emissdist_d2'));shading flat
    hold on
    caxis([-4 -0.5])
    for i = 1:hmm.Nsegs
        plot(time{i},hmm.state(1).meanfun{i}(:,2),'w','linewidth',2)
        plot(time{i},hmm.state(2).meanfun{i}(:,2),'w','linewidth',2)
    end
    xlabel('Time (s)','fontsize',16)
    ylabel('Amplitude (z)','fontsize',16)
    title('Sliding window density estimation','fontsize',16)
