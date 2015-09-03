clear all
close all

addpath('C:\WC_Germany\parietal_cortical_2010\')
addpath('C:\WC_Germany\hsmm_state_detection\')
addpath('C:\WC_Germany\hsmm_uds_code\')

cd C:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load C:\WC_Germany\parietal_cortical_2010\desynch_times_individual

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
    d=3; %d=3
    cdir = sess_data(d).directory;
    cdir(1) = 'C';
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
    lf_cut_off_freqs = [0.05 10];
    hf_cut_off_freqs = [20 80];
    smooth_sigma = 0.15; %in seconds
    [lf_features,t_axis,Fs] = hsmm_uds_get_lf_features(lf4,Fs_orig,hmm_dsf,lf_cut_off_freqs); %'low-frequency amplitude'
    
    %% LOCATE THE SEGMENTS CONTAINING UDS
    T = length(lf_features);
    min_seg_dur = 60;
    UDS_segs = hsmm_uds_get_uds_segments(desynch_times,hmm_Fs,T,min_seg_dur);
    
    %% INITIALIZE THE HMM
    params.meantype = 'variable';
    params.UDS_segs = UDS_segs;
    params.movingwin = [50 10];
    [hmm] = hsmm_uds_initialize_3state(lf_features,hmm_Fs,params);
    
    %% FIT AN HMM
    hmm.min_mllik = -8; %set minimum value of the maximal log likelihood for either state before switching to robust estimator.  If you don't want to use this feature, set this to -Inf
    hmm.windowSize = 50;
    hmm = hsmm_uds_train_hmm_3state(hmm,lf_features);
    
    %% DETERMINE THE VITERBI SEQUENCE
    [state_seq,llik_best] = hsmm_uds_viterbi_hmm(hmm,lf_features);

%% reassign state ids
state_vec = nan(size(lf_features));
for i = 1:hmm.Nsegs
    temp = state_seq{i};
    temp(temp==1) = 0;
    temp(temp==3) = 1;
    state_seq{i} = temp;
    state_vec(hmm.UDS_segs(i,1):hmm.UDS_segs(i,2)) = temp;
end
%%
dstate_vec = [nan; diff(state_vec)];
state_trans = find(dstate_vec ~= 0);
coming_from = nan(size(state_trans));
going_to = nan(size(state_trans));
coming_from(2:end) = state_vec(state_trans(2:end)-1);
going_to = state_vec(state_trans);

bad_trans = [];
beg_point = find(state_vec==0 | state_vec==2,1,'first');
cur_state = state_vec(beg_point);
cur_ind = beg_point;
while cur_ind < state_trans(end)
    if cur_state == 0
        leave_trans = find(going_to==2 & state_trans > cur_ind,1,'first');
        if isempty(leave_trans)
            break
        end
        prec_trans = leave_trans - 1;
        bad_trans = [bad_trans; find(state_trans > cur_ind & state_trans < state_trans(prec_trans))];
        cur_state = 2;
        cur_ind = state_trans(leave_trans);
    else
        leave_trans = find(going_to == 0 & state_trans > cur_ind,1,'first');
        if isempty(leave_trans)
            break
        end
        prec_trans = leave_trans - 1;
        bad_trans = [bad_trans; find(state_trans > cur_ind & state_trans < state_trans(prec_trans))];
        cur_state = 0;
        cur_ind = state_trans(leave_trans);
    end
end
state_trans(bad_trans) = [];
going_to(bad_trans) = [];
coming_from(bad_trans) = [];
new_state_vec = state_vec;
for i = 1:length(state_trans)-1
   new_state_vec(state_trans(i):state_trans(i+1)) = going_to(i);
end
new_state_vec(new_state_vec == 1) = 2;
new_state_vec(new_state_vec == 0) = 1;

