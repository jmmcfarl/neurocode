clear all
close all

load F:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('F:\Code\WC_anal\general\')
addpath('F:\WC_Germany\persistent_2010\')
addpath('F:\WC_Germany\down_state_theta\')
addpath('F:\Code\splinefit')
addpath('F:\WC_Germany\parietal_cortical_2010\')
addpath('F:\WC_Germany\hsmm_state_detection\')
lec = find_struct_field_vals(sess_data,'region','LEC');

Fs = 2016;
niqf = Fs/2;
lcf = 1.5/niqf;
hcf = 4.5/niqf;
% lcf = 0.05/niqf;
% hcf = 10/niqf;
[b,a] = butter(2,[lcf hcf]);
[b2,a2] = butter(2,[0.05/niqf 20/niqf]);
[b3,a3] = butter(2,[1.5/niqf 20/niqf]);
dsf = 8;
Fsd = Fs/dsf;

backlag = 5*Fsd;
forwardlag = 10*Fsd;
lags = (-backlag:forwardlag)/Fsd;

min_rsquared = 0.75;

for d = lec
    
    cdir = sess_data(d).directory;
    cdir(1) = 'F';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    
    load ./used_data lf8 wcv_minus_spike
    load ./theta_times
    load ./spike_time_jmm.mat
    spkid = round(spkid/dsf);
    
    wcv_f = filter(b,a,wcv_minus_spike);
    wcv_lf = filtfilt(b2,a2,wcv_minus_spike);
    wcv_hf = filtfilt(b3,a3,wcv_minus_spike);
    lf8_f = filtfilt(b2,a2,lf8);
    
    wcv_f = downsample(wcv_f,dsf);
    wcv_lf = downsample(wcv_lf,dsf);
    wcv_hf = downsample(wcv_hf,dsf);
    lf8_f = downsample(lf8_f,dsf);
    wcv_f = zscore(wcv_f);
    wcv_lf = zscore(wcv_lf);
    wcv_hf = zscore(wcv_hf);
    lf8_f = zscore(lf8_f);
    
    t_axis = (1:length(wcv_f))/Fsd;
    wcv_phase = angle(hilbert(wcv_f));
    binned_spikes = hist(spkid,1:length(wcv_f));
    
    
    load ./ec_hmm_state_seq
    mp_state_seq = hmm_bbstate_seq;
    [new_seg_inds] = resample_uds_seg_inds(hmm.UDS_segs,hmm.Fs,Fsd,length(wcv_f));
    mp_state_vec = nan(size(wcv_f));
    mp_utrans = [];
    mp_dtrans = [];
    for n = 1:hmm.Nsegs
        mp_state_vec(new_seg_inds(n,1):new_seg_inds(n,2)) = mp_state_seq{n};
        cur_mp_utrans = new_seg_inds(n,1) + find(mp_state_seq{n}(1:end-1) == 1 & mp_state_seq{n}(2:end) == 2);
        cur_mp_dtrans = new_seg_inds(n,1) + find(mp_state_seq{n}(1:end-1) == 2 & mp_state_seq{n}(2:end) == 1);
        cur_mp_dtrans(cur_mp_dtrans < cur_mp_utrans(1)) = [];
        cur_mp_utrans(cur_mp_utrans > cur_mp_dtrans(end)) = [];
        mp_utrans = [mp_utrans; cur_mp_utrans];
        mp_dtrans = [mp_dtrans; cur_mp_dtrans];
    end
    
    load ./ec_hmm_state_seq8
    lf8_state_seq = hmm_bbstate_seq8;
    [new_seg_inds] = resample_uds_seg_inds(hmm8.UDS_segs,hmm8.Fs,Fsd,length(wcv_f));
    lf8_state_vec = nan(size(wcv_f));
    lf8_utrans = [];
    lf8_dtrans = [];
    for n = 1:hmm8.Nsegs
        lf8_state_vec(new_seg_inds(n,1):new_seg_inds(n,2)) = lf8_state_seq{n};
        cur_lf8_utrans = new_seg_inds(n,1) + find(lf8_state_seq{n}(1:end-1) == 1 & lf8_state_seq{n}(2:end) == 2);
        cur_lf8_dtrans = new_seg_inds(n,1) + find(lf8_state_seq{n}(1:end-1) == 2 & lf8_state_seq{n}(2:end) == 1);
        cur_lf8_dtrans(cur_lf8_dtrans < cur_lf8_utrans(1)) = [];
        cur_lf8_utrans(cur_lf8_utrans > cur_lf8_dtrans(end)) = [];
        lf8_utrans = [lf8_utrans; cur_lf8_utrans];
        lf8_dtrans = [lf8_dtrans; cur_lf8_dtrans];
    end
    
    mp_uptimes = t_axis(mp_utrans);
    lf8_uptimes = t_axis(lf8_utrans);
    mp_theta_ups = [];
    lf8_theta_ups = [];
    for i = 1:length(theta_start_times)
        cur_theta_ups = find(mp_uptimes > theta_start_times(i) & mp_uptimes < theta_stop_times(i));
        mp_theta_ups = [mp_theta_ups cur_theta_ups];
        cur_theta_ups = find(lf8_uptimes > theta_start_times(i) & lf8_uptimes < theta_stop_times(i));
        lf8_theta_ups = [lf8_theta_ups cur_theta_ups];
    end
    
    %         [rlamp,rlshift,rltau,rsquared,t_50,t_90,t_10,fit_data,resid_data,used_fit] = ...
    %             get_dualsigmoid_fit_dst(mp_utrans,mp_dtrans,wcv_lf,t_axis,Fsd,[]);
    %         ord = 2;
    %                 acorr_seq2 = xcov(wcv_lf,ord);
    %         acorr_seq2 = acorr_seq2(ord+1:end);
    %         A = levinson(acorr_seq2(1:ord),ord);
    %         est_x2 = filter(-A,1,wcv_lf);
    
    mp_utrans = mp_utrans(mp_theta_ups);
    mp_dtrans = mp_dtrans(mp_theta_ups);
    lf8_utrans = lf8_utrans(lf8_theta_ups);
    lf8_dtrans = lf8_dtrans(lf8_theta_ups);
    
    
    n_mp_ups = length(mp_utrans);
    %initialize
    mp_utrig_mp_mat = nan(n_mp_ups,length(lags));
    mp_utrig_mpl_mat = nan(n_mp_ups,length(lags));
    mp_utrig_mph_mat = nan(n_mp_ups,length(lags));
    mp_utrig_lf8_mat = nan(n_mp_ups,length(lags));
    mp_utrig_spk_mat = nan(n_mp_ups,length(lags));
    %calculate mp utrigs
    for i = 1:n_mp_ups
        if mp_utrans(i) > backlag && length(wcv_f) - mp_utrans(i) > forwardlag
            mp_utrig_mp_mat(i,:) = wcv_phase(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);
            mp_utrig_mpl_mat(i,:) = wcv_lf(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);
            mp_utrig_mph_mat(i,:) = wcv_hf(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);
            mp_utrig_lf8_mat(i,:) = lf8_f(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);
            mp_utrig_spk_mat(i,:) = binned_spikes(mp_utrans(i)-backlag:mp_utrans(i)+forwardlag);
        end
    end
    
    mp_updur = (mp_dtrans - mp_utrans)/Fsd;
    [dummy,up_order] = sort(mp_updur);
    
    %% plot mp up trig matrices
    f1 = figure('visible','off');
    set(f1,'paperunits','centimeters','papersize',[20 30]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))]);
    subplot(4,1,1)
    pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_mpl_mat(up_order,:));
    shading flat; colorbar; hold on
    plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
    caxis([-2 2]);colorbar
    xlim([-2 2])
    line([0 0],[0 1],'Color','k')
    hold on
    for i = 1:n_mp_ups
        cur = find(mp_utrig_spk_mat(up_order(i),:) > 0);
        plot(lags(cur),ones(size(cur))*i/n_mp_ups,'w.','markersize',4)
    end
    subplot(4,1,2)
    pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_mph_mat(up_order,:));
    shading flat; colorbar; hold on
    plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
    caxis([-2 2]);colorbar
    xlim([-2 2])
    line([0 0],[0 1],'Color','k')
    hold on
    for i = 1:n_mp_ups
        cur = find(mp_utrig_spk_mat(up_order(i),:) > 0);
        plot(lags(cur),ones(size(cur))*i/n_mp_ups,'w.','markersize',4)
    end
    subplot(4,1,3)
    pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_mp_mat(up_order,:));
    shading flat; colorbar; hold on
    plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
    caxis([-pi pi]);colorbar
    xlim([-2 2])
    line([0 0],[0 1],'Color','k')
    subplot(4,1,4)
    pcolor(lags,(1:n_mp_ups)/n_mp_ups,mp_utrig_lf8_mat(up_order,:));
    shading flat; colorbar; hold on
    plot(mp_updur(up_order),(1:n_mp_ups)/n_mp_ups,'w','linewidth',2)
    caxis([-2 2]);colorbar
    xlim([-2 2])
    line([0 0],[0 1],'Color','k')
    f_name = ['G:\WC_Germany\overall_EC\theta_ulock\' s_name];
    print('-dpng',f_name);
    close
    
    
    n_lf8_ups = length(lf8_utrans);
    %initialize
    lf8_utrig_mp_mat = nan(n_lf8_ups,length(lags));
    lf8_utrig_mpl_mat = nan(n_lf8_ups,length(lags));
    lf8_utrig_mph_mat = nan(n_lf8_ups,length(lags));
    lf8_utrig_lf8_mat = nan(n_lf8_ups,length(lags));
    
    %calculate lf8 utrigs
    for i = 1:n_lf8_ups
        if lf8_utrans(i) > backlag && length(wcv_f) - lf8_utrans(i) > forwardlag
            lf8_utrig_mp_mat(i,:) = wcv_phase(lf8_utrans(i)-backlag:lf8_utrans(i)+forwardlag);
            lf8_utrig_mpl_mat(i,:) = wcv_lf(lf8_utrans(i)-backlag:lf8_utrans(i)+forwardlag);
            lf8_utrig_mph_mat(i,:) = wcv_hf(lf8_utrans(i)-backlag:lf8_utrans(i)+forwardlag);
            lf8_utrig_lf8_mat(i,:) = lf8_f(lf8_utrans(i)-backlag:lf8_utrans(i)+forwardlag);
        end
    end
    
    lf8_updur = (lf8_dtrans - lf8_utrans)/Fsd;
    [dummy,up_order] = sort(lf8_updur);
    
    %% plot lf8 up trig matrices
    f1 = figure('visible','off');
    set(f1,'paperunits','centimeters','papersize',[20 30]);
    set(f1,'PaperPosition',[0,0,(get(f1,'PaperSize'))]);
    subplot(4,1,1)
    pcolor(lags,(1:n_lf8_ups)/n_lf8_ups,lf8_utrig_mpl_mat(up_order,:));
    shading flat; colorbar; hold on
    plot(lf8_updur(up_order),(1:n_lf8_ups)/n_lf8_ups,'w','linewidth',2)
    caxis([-2 2]);colorbar
    xlim([-2 2])
    line([0 0],[0 1],'Color','k')
    subplot(4,1,2)
    pcolor(lags,(1:n_lf8_ups)/n_lf8_ups,lf8_utrig_mph_mat(up_order,:));
    shading flat; colorbar; hold on
    plot(lf8_updur(up_order),(1:n_lf8_ups)/n_lf8_ups,'w','linewidth',2)
    caxis([-2 2]);colorbar
    xlim([-2 2])
    line([0 0],[0 1],'Color','k')
    subplot(4,1,3)
    pcolor(lags,(1:n_lf8_ups)/n_lf8_ups,lf8_utrig_mp_mat(up_order,:));
    shading flat; colorbar; hold on
    plot(lf8_updur(up_order),(1:n_lf8_ups)/n_lf8_ups,'w','linewidth',2)
    caxis([-pi pi]);colorbar
    xlim([-2 2])
    line([0 0],[0 1],'Color','k')
    subplot(4,1,4)
    pcolor(lags,(1:n_lf8_ups)/n_lf8_ups,lf8_utrig_lf8_mat(up_order,:));
    shading flat; colorbar; hold on
    plot(lf8_updur(up_order),(1:n_lf8_ups)/n_lf8_ups,'w','linewidth',2)
    caxis([-2 2]);colorbar
    xlim([-2 2])
    line([0 0],[0 1],'Color','k')
    f_name = ['G:\WC_Germany\overall_EC\theta_ulock\lf8_' s_name];
    print('-dpng',f_name);
    close
    
end

clear *mat
cd G:\WC_Germany\lec_theta\
