clear all

load G:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('G:\WC_Germany\parietal_cortical_2010')

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

dsf = 8;
Fsd = 2016/dsf;

niqf = 2016/2;
hcf1 = 100/niqf;
[b,a] = butter(2,[0.1/niqf 20/niqf]);

params.Fs = Fsd;
params.err = 0;
params.tapers = [2 3];
params.fpass = [0 20];

winlength = 20;
winslide = 2.5;
movingwin = [winlength winslide];


for d = 1:length(sess_data)

    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);

    load ./used_data lf8 

    %bandlimit signals
    down_8 = filtfilt(b,a,lf8);
    down_8 = downsample(down_8,dsf);
    down_8 = zscore(down_8);    %zscore
    t_axis = (1:length(down_8))/Fsd;
    old_taxis = (1:length(downsample(lf8,40)))/(2016/40);
    
    [P8{d},t{d},f{d}]=mtspecgramc(down_8,movingwin,params);   
    dt = mean(diff(t{d}));

    [peak_pow_lf8{d},peak_freq_loc] = max(10*log10(P8{d}),[],2);
    peak_freq_lf8{d} = f{d}(peak_freq_loc);
    
    load ./pa_hsmm_state_seq_new2
    mp_state_seq = hsmm_bbstate_seq;
    load ./pa_hsmm_state_seq8_new2
    lfp_state_seq = hsmm_bbstate_seq8;
    
    lfp_meanfuns = nan(2,length(old_taxis));
    mp_meanfuns = nan(2,length(old_taxis));
    for i = 1:size(hmm.UDS_segs,1)
       mp_meanfuns(1,hmm.UDS_segs(i,1):hmm.UDS_segs(i,2)) = hmm.state(1).meanfun{i}; 
       mp_meanfuns(2,hmm.UDS_segs(i,1):hmm.UDS_segs(i,2)) = hmm.state(2).meanfun{i}; 
    end
    for i = 1:size(hmm8.UDS_segs,1)
       lfp_meanfuns(1,hmm8.UDS_segs(i,1):hmm8.UDS_segs(i,2)) = hmm8.state(1).meanfun{i}; 
       lfp_meanfuns(2,hmm8.UDS_segs(i,1):hmm8.UDS_segs(i,2)) = hmm8.state(2).meanfun{i}; 
    end
    
    mp_upvar = hmm.state(2).var;
    mp_downvar = hmm.state(1).var;
    kl1 = gauss_kl_div_fun(diff(mp_meanfuns),mp_downvar,mp_upvar);
    kl2 = gauss_kl_div_fun(-diff(mp_meanfuns),mp_upvar,mp_downvar);
    mp_kl_div = kl1+kl2;

        lfp_upvar = hmm8.state(2).var;
    lfp_downvar = hmm8.state(1).var;
    kl1 = gauss_kl_div_fun(diff(lfp_meanfuns),lfp_downvar,lfp_upvar);
    kl2 = gauss_kl_div_fun(-diff(lfp_meanfuns),lfp_upvar,lfp_downvar);
    lfp_kl_div = kl1+kl2;

%     save kl_div_funs mp_kl_div lfp_kl_div old_taxis
    
    [new_seg_inds] = resample_uds_seg_inds(hmm.UDS_segs,hmm.Fs,252,length(down_8));
    [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
    [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);
    [state_durations] = compute_state_durations_seg(mp_state_seq,Fsd);
    [state_durations8] = compute_state_durations_seg(lfp_state_seq,Fsd);
    [state_dutycycles8] = compute_state_dutycycles_seg(lfp_state_seq);
    
    lf8_udsfreq{d} = nan(size(up_trans_inds));
    lf8_udspow{d} = nan(size(up_trans_inds));
    lf8_dutycyc{d} = nan(size(up_trans_inds));
    lf8_medupdur{d} = nan(size(up_trans_inds));
    lf8_meddowndur{d} = nan(size(up_trans_inds));
    for i = 1:length(t{d})
       cur_up8 = find(t_axis(up_trans_inds8) >= t{d}(i)-winlength/2 & ...
           t_axis(up_trans_inds8) < t{d}(i)+winlength/2);
       cur_down8 = find(t_axis(down_trans_inds8) >= t{d}(i)-winlength/2 & ...
           t_axis(down_trans_inds8) < t{d}(i)+winlength/2);
       cur_up = find(t_axis(up_trans_inds) >= t{d}(i)-winlength/2 & ...
           t_axis(up_trans_inds) < t{d}(i)+winlength/2);
       lf8_udsfreq{d}(cur_up) = peak_freq_lf8{d}(i);
       lf8_udspow{d}(cur_up) = peak_pow_lf8{d}(i);
       lf8_dutycyc{d}(cur_up) = nanmedian(state_dutycycles8(cur_up8));
       lf8_medupdur{d}(cur_up) = nanmedian(state_durations8{2}(cur_up8));
       lf8_meddowndur{d}(cur_up) = nanmedian(state_durations8{1}(cur_up8));
       lf8_medcycdur{d}(cur_up) = nanmedian(state_durations8{1}(cur_up8)+state_durations8{2}(cur_up8));       
    end
    
    lf8_sep{d} = nan(size(up_trans_inds));
    mp_sep{d} = nan(size(up_trans_inds));
    for i = 1:length(up_trans_inds)
        cur_point = find(old_taxis >= t_axis(up_trans_inds(i)),1,'first');
        lf8_sep{d}(i) = lfp_kl_div(cur_point);
        mp_sep{d}(i) = mp_kl_div(cur_point);
    end
    
end

cd G:\WC_Germany\persistent_9_27_2010\
save depthofanesth_data_new2 lf8_* mp_sep