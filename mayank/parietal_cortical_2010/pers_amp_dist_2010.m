cd G:\WC_Germany\persistent_2010

clear all
close all

load pers_revised_dir_2010
drive_letter = 'G';

Fs=  2016;
dsf = 8;
Fsd = Fs/dsf;
niqf = Fs/2;
lcf = 0.05/niqf;
hcf = 40/niqf;
[b,a] = butter(2,[lcf hcf]);

amp_range = linspace(-4, 4, 400);

for d = 1:length(dir_array)

    cdir = dir_array{d};
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);

    load used_data wcv_minus_spike lf8
    load hsmm_state_seq8_seg_15
    lfp_state_seq = hsmm_bbstate_seq8;
    load hsmm_state_seq_seg_15
    mp_state_seq = hsmm_bbstate_seq;
    
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);
    wcv_f = zscore(downsample(wcv_f,dsf));
    lf8_f = zscore(downsample(lf8_f,dsf));

    [new_seg_inds] = resample_uds_seg_inds(hmm.UDS_segs,hmm.Fs,Fsd,length(wcv_f));

    mp_up_marg = [];
    mp_down_marg = [];
    lfp_up_marg = [];
    lfp_down_marg = [];
    mp_up_lup = [];
    mp_up_ldown = [];
    mp_down_lup = [];
    mp_down_ldown = [];
    
    for ns = 1:hmm.Nsegs 
        cur_mp_seg = wcv_f(new_seg_inds(ns,1):new_seg_inds(ns,2));
        cur_lfp_seg = lf8_f(new_seg_inds(ns,1):new_seg_inds(ns,2));
        cur_mp_state = logical(mp_state_seq{ns}-1);
        cur_lfp_state = logical(lfp_state_seq{ns}-1);
        mp_up_marg = [mp_up_marg; cur_mp_seg(cur_mp_state)];
        mp_down_marg = [mp_down_marg; cur_mp_seg(~cur_mp_state)];
        lfp_up_marg = [lfp_up_marg; cur_lfp_seg(cur_lfp_state)];
        lfp_down_marg = [lfp_down_marg; cur_lfp_seg(~cur_lfp_state)];
        mp_up_lup = [mp_up_lup; cur_mp_seg(cur_mp_state & cur_lfp_state)];
        mp_up_ldown = [mp_up_ldown; cur_mp_seg(cur_mp_state & ~cur_lfp_state)];
        mp_down_lup = [mp_down_lup; cur_mp_seg(~cur_mp_state & cur_lfp_state)];
        mp_down_ldown = [mp_down_ldown; cur_mp_seg(~cur_mp_state & ~cur_lfp_state)];
    end
    
    d_mp_up_marg(d,:) = gpkde(mp_up_marg, -3, [-4; 4; 400]);
    d_mp_down_marg(d,:) = gpkde(mp_down_marg, -3, [-4; 4; 400]);    
    d_lfp_up_marg(d,:) = gpkde(lfp_up_marg, -3, [-4; 4; 400]);
    d_lfp_down_marg(d,:) = gpkde(lfp_down_marg, -3, [-4; 4; 400]);    
    d_mp_up_lup(d,:) = gpkde(mp_up_lup, -3, [-4; 4; 400]);
    d_mp_down_lup(d,:) = gpkde(mp_down_lup, -3, [-4; 4; 400]);    
    d_mp_up_ldown(d,:) = gpkde(mp_up_ldown, -3, [-4; 4; 400]);
    d_mp_down_ldown(d,:) = gpkde(mp_down_ldown, -3, [-4; 4; 400]); 
    
    plot(amp_range,d_mp_up_marg(d,:),'r'),hold on
    plot(amp_range,d_mp_down_marg(d,:))
    plot(amp_range,d_mp_up_lup(d,:),'g')
    plot(amp_range,d_mp_up_ldown(d,:),'k')
    plot(amp_range,d_mp_down_lup(d,:),'g--')
    plot(amp_range,d_mp_down_ldown(d,:),'k--')
   t_names = ['G:\WC_Germany\persistent_2010\state_dep_neur_dist\' f_names{d}];
    print(t_names,'-dpng')
    close
    
   clear mp* lfp*

end

cd G:\WC_Germany\persistent_2010\
save amp_dist amp_range* d_*