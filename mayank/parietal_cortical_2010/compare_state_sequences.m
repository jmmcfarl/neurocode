clear all
close all

addpath('F:\WC_Germany\parietal_cortical_2010\')
addpath('F:\WC_Germany\hsmm_state_detection\')

cd F:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load F:\WC_Germany\parietal_cortical_2010\desynch_times_individual	


% sess_data = sess_data(parietal);
% desynch_times = desynch_times(parietal);

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_times_mp(interneurons) = [];
desynch_times_lf8(interneurons) = [];
desynch_times_lf4(interneurons) = [];

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
desynch_times_mp = desynch_times_mp(thom_pfc);
desynch_times_lf4 = desynch_times_lf4(thom_pfc);
desynch_times_lf8 = desynch_times_lf8(thom_pfc);

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;


% sess_data(42:end) = [];
n = length(sess_data);
for d = 1:n
    cdir = sess_data(d).directory;
    cdir(1) = 'F';
    cd(cdir)
    pwd

    load hsmm_state_seq_seg_lf_4_10_10
    load fixmean_state_seqs_zm
    fixmean_state_seq_zm = fixmean_state_seq;
    load fixmean_state_seqs_4_10_10
%     load hsmm_state_seq_seg_lfhf
    mp_state_seq = hsmm_bbstate_seq;
    mp_hmm_state_seq = hmm_bbstate_seq;
%     mplfhf_state_seq = hsmm_bbstate_seq_lfhf;
%     for i = 1:hmm.Nsegs
%         mp_state_seq_d{i} = downsample(mp_state_seq{i},5);
%     end

    load hsmm_state_seq8_seg_lf_4_10_10
    load hsmm_state_seq8_seg_hf_4_10_10
    load fixmean_state_seqs8_4_10_10
    load phase_state_seq8_4_10_10
    load hsmm_state_seq8_seg_lfhf
    lf8_state_seq = hsmm_bbstate_seq8;
    lf8_hmm_state_seq = hmm_bbstate_seq8;
    lf8hf_state_seq = hsmm_bbstate_seq8_hf;
    lf8p_state_seq = phase_state_seq8;
    lf8lfhf_state_seq = hsmm_bbstate_seq8_lfhf;
%     for i = 1:hmm8.Nsegs
%         lf8_state_seq_d{i} = downsample(lf8_state_seq{i},5);
%     end
%     lf8_hmm_state_seq = hmm_bbstate_seq8;

    if sess_data(d).thom_elec
       load hsmm_state_seq4_seg_lf_4_10_10
       load hsmm_state_seq4_seg_hf_4_10_10
       load fixmean_state_seqs4_zm
       fixmean_state_seq4_zm = fixmean_state_seq4;
       load fixmean_state_seqs4_4_10_10
       load phase_state_seq4_4_10_10
       load hsmm_state_seq4_seg_lfhf
%        load hsmm_state_seq84_seg_lf_duallfp
       lf4_state_seq = hsmm_bbstate_seq4;
       lf4hf_state_seq = hsmm_bbstate_seq4_hf;
       lf4p_state_seq = phase_state_seq4;
       lf4lfhf_state_seq = hsmm_bbstate_seq4_lfhf;
       lf4_hmm_state_seq = hmm_bbstate_seq4;
%        lf84dual_state_seq = hsmm_bbstate_seq84;
%        for i = 1:hmm4.Nsegs
%            lf4_state_seq_d{i} = downsample(lf4_state_seq{i},5);
%        end
    end
        
    mp_lf8_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_state_seq,lf8_state_seq,hmm.UDS_segs,hmm8.UDS_segs,Fsd);
    mp_lf8hf_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_state_seq,lf8hf_state_seq,hmm.UDS_segs,hmm8.UDS_segs,Fsd);
    mp_lf8p_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_state_seq,lf8p_state_seq,hmm.UDS_segs,hmm8.UDS_segs,Fsd);
    fmp_flf8_ham(d) = compute_state_seq_seg_hamdist_varuds(fixmean_state_seq,fixmean_state_seq8,hmm.UDS_segs,hmm8.UDS_segs,Fsd);
    mp_flf8_ham(d) = compute_state_seq_seg_hamdist_varuds(fixmean_state_seq,lf8_state_seq,hmm.UDS_segs,hmm8.UDS_segs,Fsd);
    fmp_flf8hf_ham(d) = compute_state_seq_seg_hamdist_varuds(fixmean_state_seq,fixmean_state_seq8_hf,hmm.UDS_segs,hmm8.UDS_segs,Fsd);
    fmp_mp_ham(d) = compute_state_seq_seg_hamdist_varuds(fixmean_state_seq,mp_state_seq,hmm.UDS_segs,hmm.UDS_segs,Fsd);
    flf8_lf8_ham(d) = compute_state_seq_seg_hamdist_varuds(fixmean_state_seq8,lf8_state_seq,hmm8.UDS_segs,hmm8.UDS_segs,Fsd);
    lf8p_lf8_ham(d) = compute_state_seq_seg_hamdist_varuds(lf8p_state_seq,lf8_state_seq,hmm8.UDS_segs,hmm8.UDS_segs,Fsd);
    mp_lf8lh_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_state_seq,lf8lfhf_state_seq,hmm.UDS_segs,hmm8_lfhf.UDS_segs,Fsd);
%     mplh_lf8lh_ham(d) = compute_state_seq_seg_hamdist_varuds(mplfhf_state_seq,lf8lfhf_state_seq,hmm.UDS_segs,hmm8.UDS_segs,Fsd);
    mph_lf8h_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_hmm_state_seq,lf8_hmm_state_seq,hmm.UDS_segs,hmm8.UDS_segs,Fsd);

    [mp_lf8_fp(d),mp_lf8_fn(d)] = compute_state_seq_seg_roc_varuds(mp_state_seq,lf8_state_seq,hmm.UDS_segs,hmm8.UDS_segs,Fsd);
    [mp_lf8lfhf_fp(d),mp_lf8lfhf_fn(d)] = compute_state_seq_seg_roc_varuds(mp_state_seq,lf8lfhf_state_seq,hmm.UDS_segs,hmm8.UDS_segs,Fsd);
    
    if sess_data(d).thom_elec
       mp_lf4_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_state_seq,lf4_state_seq,hmm.UDS_segs,hmm4.UDS_segs,Fsd); 
       mp_lf4p_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_state_seq,lf4p_state_seq,hmm.UDS_segs,hmm4.UDS_segs,Fsd); 
       fmp_flf4_ham(d) = compute_state_seq_seg_hamdist_varuds(fixmean_state_seq,fixmean_state_seq4,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
       fmp_flf4_z_ham(d) = compute_state_seq_seg_hamdist_varuds(fixmean_state_seq_zm,fixmean_state_seq4_zm,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
        mp_flf4_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_state_seq,fixmean_state_seq4,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
       mp_lf4hf_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_state_seq,lf4hf_state_seq,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
       lf8_lf4_ham(d) = compute_state_seq_seg_hamdist_varuds(lf8_state_seq,lf4_state_seq,hmm8.UDS_segs,hmm4.UDS_segs,Fsd);
       lf8hf_lf4hf_ham(d) = compute_state_seq_seg_hamdist_varuds(lf8hf_state_seq,lf4hf_state_seq,hmm8.UDS_segs,hmm4.UDS_segs,Fsd);
       flf8_flf4_ham(d) = compute_state_seq_seg_hamdist_varuds(fixmean_state_seq8,fixmean_state_seq4,hmm8.UDS_segs,hmm4.UDS_segs,Fsd);
       lf8p_lf4p_ham(d) = compute_state_seq_seg_hamdist_varuds(lf8p_state_seq,lf4p_state_seq,hmm8.UDS_segs,hmm4.UDS_segs,Fsd);
       mp_lf4lh_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_state_seq,lf4lfhf_state_seq,hmm.UDS_segs,hmm4_lfhf.UDS_segs,Fsd);
%        mplh_lf4lh_ham(d) = compute_state_seq_seg_hamdist_varuds(mplfhf_state_seq,lf4lfhf_state_seq,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
%        lf8lh_lf4lh_ham(d) = compute_state_seq_seg_hamdist_varuds(lf8lfhf_state_seq,lf4lfhf_state_seq,hmm8.UDS_segs,hmm4.UDS_segs,Fsd);
%         mp_lf84_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_state_seq,lf84dual_state_seq,hmm.UDS_segs,hmm84.UDS_segs,Fsd);
    
        [mp_lf4_fp(d),mp_lf4_fn(d)] = compute_state_seq_seg_roc_varuds(mp_state_seq,lf4_state_seq,hmm.UDS_segs,hmm4.UDS_segs,Fsd);
        [mp_lf4lfhf_fp(d),mp_lf4lfhf_fn(d)] = compute_state_seq_seg_roc_varuds(mp_state_seq,lf4lfhf_state_seq,hmm.UDS_segs,hmm4_lfhf.UDS_segs,Fsd);
        mph_lf4h_ham(d) = compute_state_seq_seg_hamdist_varuds(mp_hmm_state_seq,lf4_hmm_state_seq,hmm.UDS_segs,hmm4.UDS_segs,Fsd);

    else
        mp_lf4_ham(d) = nan;
        mp_lf4p_ham(d) = nan;
        fmp_flf4_ham(d) = nan;
        mp_flf4_ham(d) = nan;
        mp_lf4hf_ham(d) = nan;
        lf8_lf4_ham(d) = nan;
        lf8hf_lf4hf_ham(d) = nan;
        flf8_flf4_ham(d) = nan;
        lf8p_lf4p_ham(d) = nan;
        mp_lf4lh_ham(d) = nan;
%         mplh_lf4lh_ham(d) = nan;
        lf8lh_lf4lh_ham(d) = nan;
        mp_lf84_ham(d) = nan;
        
        mp_lf4_fp(d) = nan; mp_lf4_fn(d) = nan;
        mp_lf4lfhf_fp(d) = nan; mp_lf4lfhf_fn(d) = nan;
        mph_lf4h_ham(d) = nan;
    end
    
clear mp_state_seq* lf8_state_seq* lf4_state_seq* lf8p_state* lf4p_state* fixmean*

end


thom_el = find_struct_field_vals(sess_data,'thom_elec',1);
thom_par = thom_el(find(ismember(thom_el,parietal)));
thom_pfc = setdiff(thom_el,thom_par);
parietal = find_struct_field_vals(sess_data,'region','parietal');
frontal = setdiff(1:length(sess_data),parietal);
superficial = find_struct_field_vals(sess_data,'layer','23');
deep = setdiff(1:length(sess_data),superficial);

mp_lfp_ham = [mp_lf8_ham(parietal) mp_lf4_ham(thom_pfc)];
mp_lfphf_ham = [mp_lf8hf_ham(parietal) mp_lf4hf_ham(thom_pfc)];
mp_lfplfhf_ham = [mp_lf8lh_ham(parietal) mp_lf4lh_ham(thom_pfc)];
% mplfhf_lfplfhf_ham = [mplh_lf8lh_ham(parietal) mplh_lf4lh_ham(thom_pfc)];
fmp_flfp_ham = [fmp_flf8_ham(parietal) fmp_flf4_ham(thom_pfc)];
mp_flfp_ham = [mp_flf8_ham(parietal) mp_flf4_ham(thom_pfc)];
mp_plfp_ham = [mp_lf8p_ham(parietal) mp_lf4p_ham(thom_pfc)];


