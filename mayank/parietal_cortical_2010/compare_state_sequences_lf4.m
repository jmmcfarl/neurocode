clear all
close all

addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_individual

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
% desynch_times_mp(interneurons) = [];
% desynch_times_lf8(interneurons) = [];
% desynch_times_lf4(interneurons) = [];

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

% thom_pfc = [thom_pfc 43 44];
sess_data = sess_data(thom_pfc);
% desynch_times_mp = desynch_times_mp(thom_pfc);
% desynch_times_lf4 = desynch_times_lf4(thom_pfc);
% desynch_times_lf8 = desynch_times_lf8(thom_pfc);

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;


n = length(sess_data);
for d = 1:9
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    cd(cdir)
    pwd
    
    load hsmm_state_seq_seg_fm_4_28_10_v3
    mp_state_seq_fm = hsmm_bbstate_seq_fm;    
    load hsmm_state_seq_seg_lf_4_28_10_v3
    load fixmean_state_seqs_zm_4_26_10_v2
    fixmean_state_seq_zm = fixmean_state_seq;
    load fixmean_state_seqs_4_26_10_v2
    mp_state_seq = hsmm_bbstate_seq;
    mp_hmm_state_seq = hmm_bbstate_seq;
    
    
        load hsmm_state_seq4_seg_fm_4_28_10_v3
    lf4_state_seq_fm = hsmm_bbstate_seq4_fm;
    load hsmm_state_seq4_seg_lf_4_28_10_v3
    load hsmm_state_seq4_seg_hf_4_28_10_v3
    load fixmean_state_seqs4_zm_4_26_10_v2
    fixmean_state_seq4_zm = fixmean_state_seq4;
    load fixmean_state_seqs4_4_26_10_v2
    load phase_state_seq4_4_26_10_v2
    load hsmm_state_seq4_seg_lfhf_4_28_10_v3_3
    lf4_state_seq = hsmm_bbstate_seq4;
    lf4hf_state_seq = hsmm_bbstate_seq4_hf;
    lf4p_state_seq = phase_state_seq4;
    lf4lfhf_state_seq = hsmm_bbstate_seq4_lfhf;
    lf4_hmm_state_seq = hmm_bbstate_seq4;
    lf4hlfhf_state_seq = hmm_state_seq4_lfhf;
    
%     fixmean_state_seq4_zm = fixmean_sm_state_seq4;
%     fixmean_state_seq4 = fixmean_sm_state_seq4;
    
    mp_lf4_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_state_seq,lf4_state_seq,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
    mp_lf4p_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_state_seq,lf4p_state_seq,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
    fmp_flf4_ham(d) = compute_state_seq_seg_hamdist_varuds(fixmean_state_seq,fixmean_state_seq4,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
    fmp_flf4_z_ham(d) = compute_state_seq_seg_hamdist_varuds(fixmean_state_seq_zm,fixmean_state_seq4_zm,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
    mp_flf4_z_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_state_seq,fixmean_state_seq4_zm,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
    mp_flf4_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_state_seq,fixmean_state_seq4,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
    mp_lf4hf_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_state_seq,lf4hf_state_seq,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
    mp_lf4lh_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_state_seq,lf4lfhf_state_seq,hmm.UDS_segs,hmm4_lfhf.UDS_segs,Fsd);
    mpfm_lf4fm_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_state_seq_fm,lf4_state_seq_fm,hmm_fm.UDS_segs,hmm4_fm.UDS_segs,Fsd);
%     [mp_lf4_fp(d),mp_lf4_fn(d)] = compute_state_seq_seg_roc_varuds(mp_state_seq,lf4_state_seq,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
%     [mp_lf4lfhf_fp(d),mp_lf4lfhf_fn(d)] = compute_state_seq_seg_roc_varuds(mp_state_seq,lf4lfhf_state_seq,hmm.UDS_segs,hmm4_lfhf.UDS_segs,Fsd);
    mph_lf4h_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_hmm_state_seq,lf4_hmm_state_seq,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
    mp_lf4h_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_state_seq,lf4_hmm_state_seq,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
    mph_lf4hlfhf_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_hmm_state_seq,lf4hlfhf_state_seq,hmm.UDS_segs,hmm4_lfhf.UDS_segs,Fsd);
    
    
    clear mp_state_seq* lf8_state_seq* lf4_state_seq* lf8p_state* lf4p_state* fixmean*
    
end

%quantify performance relative to LF HSMM
dfm = -(mp_lf4_ham - mpfm_lf4fm_ham)./mp_lf4_ham;
dlh = -(mp_lf4_ham - mp_lf4lh_ham)./mp_lf4_ham;
dh = -(mp_lf4_ham - mph_lf4h_ham)./mp_lf4_ham;
df = -(mp_lf4_ham - fmp_flf4_ham)./mp_lf4_ham;
dfz = -(mp_lf4_ham - fmp_flf4_z_ham)./mp_lf4_ham;
dhf = -(mp_lf4_ham - mp_lf4hf_ham)./mp_lf4_ham;
dp = -(mp_lf4_ham - mp_lf4p_ham)./mp_lf4_ham;

%%
figure
% boxplot([dlh' dh' dfm' df' dfz'])
boxplot([dlh' dfm' dh' df' dfz' dp' dhf'])
line([0 6],[0 0])
% ylim([-0.1 0.6])
xlim([0.5 7.5])

figure
% boxplot([mp_lf4_ham' mp_lf4lh_ham' mph_lf4h_ham' mpfm_lf4fm_ham' fmp_flf4_ham' fmp_flf4_z_ham'])
boxplot([mp_lf4_ham' mp_lf4lh_ham' mpfm_lf4fm_ham' mph_lf4h_ham' fmp_flf4_ham' fmp_flf4_z_ham' mp_lf4p_ham' mp_lf4hf_ham'])

