clear all
close all

dsf = 8;
raw_Fs = 2016;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
[b,a] = butter(2,[0.1/niqf 45/niqf]);
[b2,a2] = butter(2,[1/niqf 10/niqf]);
[b3,a3] = butter(2,[2/niqf 4/niqf]);

load F:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('F:\Code\WC_anal\general\')
addpath('F:\WC_Germany\persistent_2010\')
addpath('F:\WC_Germany\overall_EC')
addpath('F:\WC_Germany\parietal_cortical_2010\')
addpath('F:\WC_Germany\hmm_state_detect\')
addpath('F:\WC_Germany\hmm_state_detect\maphmmbox\')
addpath('F:\WC_Germany\hsmm_state_detection\')
addpath('F:\Code\fullBNT-1.0.4\KPMstats\')
addpath('F:\Code\fullBNT-1.0.4\netlab3.3\')

load F:\WC_Germany\lec_theta\theta_classifications
pkfreqs = [];
pkamps = [];
pkpows = [];
ds_vars = [];
ds_theta = [];
ds_ltheta = [];
all_theta = [];

gmm_options(3) = 1e-15; %tolerance
gmm_options(5) = 1;%reset cov if singular values
gmm_options(14) = 100; %max iterations

params.Fs = Fsd;
params.tapers = [3 5];
params.err = 0;
movingwin = [20 1];
params.fpass = [0.01 45];

%%
for d = sup_lec
    
    cdir = sess_data(d).directory;
    cdir(1) = 'F';
    disp(sprintf('session %d',d))
    cd(cdir);
    load ./used_data wcv_minus_spike
    
    wcv = filtfilt(b,a,wcv_minus_spike);
    wcvf2 = filtfilt(b2,a2,wcv_minus_spike);
    wcvf3 = filtfilt(b3,a3,wcv_minus_spike);
    wcv = downsample(wcv,dsf);
    wcvf2 = downsample(wcvf2,dsf);
    wcvf3 = downsample(wcvf3,dsf);
    %         wcv = wcv/sess_data(d).gains(1);
    wcv = zscore(wcv);
    wcvf2 = zscore(wcvf2);
    wcvf3 = zscore(wcvf3);
    theta_amp = abs(hilbert(wcvf3));
    
    df = Fsd/length(wcv);
    t_axis = (1:length(wcv))/Fsd;
    
    load ./ec_hmm_state_seq
    mp_state_seq = hmm_bbstate_seq;
    [new_seg_inds] = resample_uds_seg_inds_v2(hmm.UDS_segs,hmm.Fs,Fsd,mp_state_seq);
    
    mp_state_vec = nan(size(wcv));
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
    
    [Sg,tg,fg]=mtspecgramc(wcv(:),movingwin,params);
    
    cur_ds_var = nan(size(tg));
    cur_ds_theta = nan(size(tg));
    cur_ds_ltheta = nan(size(tg));
    cur_all_theta = nan(size(tg));
    for i = 1:length(tg)
        cur_id = find(t_axis >= tg(i)-movingwin(1)/2 & t_axis <= tg(i)+movingwin(1)/2);
        down_id = cur_id(mp_state_vec(cur_id) == 1);
        %         up_id = cur_id(mp_state_vec(cur_id) == 2);
        if length(down_id) > 5
            cur_ds_var(i) = var(wcvf2(down_id));
            cur_ds_theta(i) = mean(theta_amp(down_id));
            cur_ds_ltheta(i) = mean(log(theta_amp(down_id)));
        end
        cur_all_theta(i) = mean(theta_amp(cur_id));
    end
    
    frange = find(fg >= 1.5 & fg <= 4.5);
    
    [cur_pkpows,pklocs] = max(10*log(Sg(:,frange)),[],2);
    [cur_pkamps,pklocs] = max(sqrt(Sg(:,frange)),[],2);
    cur_pkfreqs = fg(frange(pklocs));
    
    if sig_theta_peak(d) == 1
        pkpows = [pkpows; cur_pkpows(:)];
        pkamps = [pkamps; cur_pkamps(:)];
        pkfreqs = [pkfreqs; cur_pkfreqs(:)];
        ds_vars = [ds_vars; cur_ds_var(:)];
        ds_theta = [ds_theta; cur_ds_theta(:)];
        ds_ltheta = [ds_ltheta; cur_ds_ltheta(:)];
        all_theta = [all_theta; cur_all_theta(:)];
        
        cur_samps = cur_ds_theta(:);
        cur_samps(isnan(cur_samps)) = [];
        mix=gmm(1,2,'diag');
        [mix, netlaboptions] = gmmem(mix, cur_samps(:), gmm_options);
        meandiff = mix.centres(1)-mix.centres(2);
        kl1 = gauss_kl_div(meandiff,mix.covars(1),mix.covars(2));
        kl2 = gauss_kl_div(-meandiff,mix.covars(2),mix.covars(1));
        ds_theta_sep(d) = kl1+kl2;
        
        cur_samps = cur_ds_ltheta(:);
        cur_samps(isnan(cur_samps)) = [];
        mix=gmm(1,2,'diag');
        [mix, netlaboptions] = gmmem(mix, cur_samps(:), gmm_options);
        meandiff = mix.centres(1)-mix.centres(2);
        kl1 = gauss_kl_div(meandiff,mix.covars(1),mix.covars(2));
        kl2 = gauss_kl_div(-meandiff,mix.covars(2),mix.covars(1));
        ds_ltheta_sep(d) = kl1+kl2;
        
        cur_samps = cur_ds_var(:);
        cur_samps(isnan(cur_samps)) = [];
        mix=gmm(1,2,'diag');
        [mix, netlaboptions] = gmmem(mix, cur_samps(:), gmm_options);
        meandiff = mix.centres(1)-mix.centres(2);
        kl1 = gauss_kl_div(meandiff,mix.covars(1),mix.covars(2));
        kl2 = gauss_kl_div(-meandiff,mix.covars(2),mix.covars(1));
        ds_var_sep(d) = kl1+kl2;
        
        cur_samps = cur_pkamps(:);
        cur_samps(isnan(cur_samps)) = [];
        mix=gmm(1,2,'diag');
        [mix, netlaboptions] = gmmem(mix, cur_samps(:), gmm_options);
        meandiff = mix.centres(1)-mix.centres(2);
        kl1 = gauss_kl_div(meandiff,mix.covars(1),mix.covars(2));
        kl2 = gauss_kl_div(-meandiff,mix.covars(2),mix.covars(1));
        pkamp_sep(d) = kl1+kl2;
        
        cur_samps = cur_pkpows(:);
        cur_samps(isnan(cur_samps)) = [];
        mix=gmm(1,2,'diag');
        [mix, netlaboptions] = gmmem(mix, cur_samps(:), gmm_options);
        meandiff = mix.centres(1)-mix.centres(2);
        kl1 = gauss_kl_div(meandiff,mix.covars(1),mix.covars(2));
        kl2 = gauss_kl_div(-meandiff,mix.covars(2),mix.covars(1));
        pkpow_sep(d) = kl1+kl2;
    end
    
    %now run HMM
    cur_ds_var(isnan(cur_ds_var)) = nanmean(cur_ds_var);
    cur_ds_theta(isnan(cur_ds_theta)) = nanmean(cur_ds_theta);
    
    %     features = zscore([cur_ds_var(:) cur_pkpows(:)]);
    features = zscore([cur_all_theta(:)]);
    [hmm] = jmm_initialize_hmm(features,Fsd,2);
    [hmm] = jmm_train_hmm(hmm,features,1);
    [state_seq,lik_best]=jmm_hmm_viterbi(features,hmm);
    
    ord = 2;
    acorr_seq = xcov(wcv,ord);
    acorr_seq = acorr_seq(ord+1:end);
    A = levinson(acorr_seq(1:ord),ord);
    est_x2 = filter(-A,1,wcv);
    [Sn,tn,fn]=mtspecgramc(est_x2(:),movingwin,params);
    
    lSg1 = mean(log10(Sn(state_seq==1,:)));
    lSg2 = mean(log10(Sn(state_seq==2,:)));
    psg1 = max(lSg1(frange));
    psg2 = max(lSg2(frange));
    
    theta_state = zeros(size(state_seq));
    if psg1 > psg2
        theta_state(state_seq == 1) = 1;
    else
        theta_state(state_seq == 2) = 1;
    end
    mean_dsvar_thetastate(d) = nanmean(cur_ds_var(theta_state==1));
    mean_ldsvar_thetastate(d) = nanmean(log(cur_ds_var(theta_state==1)));
    mean_thetapow_thetastate(d) = nanmean(cur_pkpows(theta_state==1));
    mean_dsvar_nonthetastate(d) = nanmean(cur_ds_var(theta_state==0));
    mean_ldsvar_nonthetastate(d) = nanmean(log(cur_ds_var(theta_state==0)));
    mean_thetapow_nonthetastate(d) = nanmean(cur_pkpows(theta_state==0));
    
    theta_up = 1+find(theta_state(2:end) == 1 & theta_state(1:end-1) == 0);
    theta_down = 1+find(theta_state(2:end) == 0 & theta_state(1:end-1) == 1);
    theta_start_times = tg(theta_up);
    theta_stop_times = tg(theta_down);
    if theta_state(1) == 1
        theta_start_times = [0 theta_start_times];
    end
    if theta_state(end) == 1
        theta_stop_times = [theta_stop_times (tg(end)+movingwin(1))];
    end
    theta_durs = theta_stop_times-theta_start_times;
    total_theta_dur(d) = sum(theta_durs);
    theta_dur_fract(d) = total_theta_dur(d)/(length(wcvf2)/Fsd);
    
%         figure('visible','off')
%         subplot(2,1,1)
    
    if sig_theta_peak(d) == 1
        figure
        pcolor(tg,fg,log10(Sg)');shading flat
        caxis([-3.5 -0])
        ylim([0 10])
        hold on
        plot(tg,cur_ds_theta*2,'w','linewidth',1)
%         plot(tg,cur_pkpows+5,'k','linewidth',1)
        plot(tg,state_seq*2,'r','linewidth',2)
        s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer(1),'_',sess_data(d).name);
        f_name = ['F:\WC_Germany\overall_EC\normalized_mp_spectrogram\sig_all_theta_amp_2' s_name];
        print(f_name,'-dpng'), close
    end


%     % figure
% %     hist3_plot(cur_pkspows,log(cur_ds_var),[40 40],[],[],[],0,2)
% %     caxis([0 3])
% %     pause
% %     close all
%     
%     imagesc(tg,fg,10*log10(Sg)');
%     caxis([-35 0])
%     ylim([0 5])
%     set(gca,'FontName','Arial','Fontsize',12)
%     xlabel('Time (s)','Fontsize',12,'fontname','Arial')
%     ylabel('Frequency (Hz)','Fontsize',12,'fontname','Arial')
   
    
%     figure
% %     plot(tg,zscore(cur_ds_theta),'k')
% %     hold on
% %     plot(tg,zscore(cur_pkpows),'r')
%     plot(tg,cur_all_theta,'k')
%    hold on
%    plot(tg,state_seq,'r')
%    set(gca,'FontName','Arial','Fontsize',12)
%    xlabel('Time (s)','Fontsize',12,'fontname','Arial')
%    ylabel('Delta-band Power (AU)','Fontsize',12,'fontname','Arial')
% ylim([0 3])

    save theta_times_all_2 theta_start_times theta_stop_times
    
end

cd F:\WC_Germany\lec_theta\
save theta_times_class_data total_theta_dur theta_dur_fract
