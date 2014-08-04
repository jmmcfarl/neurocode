clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\overall_calcs\UDS_synch_state_dur\UDS_synch_state_dur_data

Fs = 2016;
niqf = Fs/2;
lcf = 2.5/niqf;
hcf = 10/niqf;
[b,a] = butter(1,[lcf hcf]);
dsf = 8;
Fsd = Fs/dsf;
nbins = 300;
le = -3;
re = 4;

backlag = 1*Fsd;
forwardlag = 1*Fsd;
lags = (-backlag:forwardlag)/Fsd;

for d = 1:length(over_dir)


    cd(over_dir{d})
    disp(num2str(d))

    load used_data lf8 lf2 lf3 wcv_minus_spike

    wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);
    lf2_f = filtfilt(b,a,lf2);
    lf3_f = filtfilt(b,a,lf3);

    down_w = downsample(wcv_f,dsf);
    down_8 = downsample(lf8_f,dsf);
    down_2 = downsample(lf2_f,dsf);
    down_3 = downsample(lf3_f,dsf);
    
    temp_cor = corrcoef(down_3,down_8);
    down_3s = down_3 - temp_cor(2,1)*down_8;

%     down_w = sqrt(jmm_smooth_1d(down_w.^2,10));
%         down_8 = sqrt(jmm_smooth_1d(down_8.^2,10));
%     down_3 = sqrt(jmm_smooth_1d(down_3.^2,10));
%     down_2 = sqrt(jmm_smooth_1d(down_2.^2,10));
%     down_3s = sqrt(jmm_smooth_1d(down_3s.^2,10));
    
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    down_3 = zscore(down_3);
    down_2 = zscore(down_2);
    down_3s = zscore(down_3s);

    %initialize
    mp_utrig_mp_mat = zeros(length(synch_ups{d}),length(lags));
    mp_utrig_lf8_mat = zeros(length(synch_ups{d}),length(lags));
    mp_utrig_lf3_mat = zeros(length(synch_ups{d}),length(lags));
    mp_utrig_lf3s_mat = zeros(length(synch_ups{d}),length(lags));
    mp_utrig_lf2_mat = zeros(length(synch_ups{d}),length(lags));

    lf8_utrig_mp_mat = zeros(length(synch_ups8{d}),length(lags));
    lf8_utrig_lf8_mat = zeros(length(synch_ups8{d}),length(lags));
    lf8_utrig_lf3_mat = zeros(length(synch_ups8{d}),length(lags));
    lf8_utrig_lf3s_mat = zeros(length(synch_ups8{d}),length(lags));
    lf8_utrig_lf2_mat = zeros(length(synch_ups8{d}),length(lags));

    mp_dtrig_mp_mat = zeros(length(synch_downs{d}),length(lags));
    mp_dtrig_lf8_mat = zeros(length(synch_downs{d}),length(lags));
    mp_dtrig_lf3_mat = zeros(length(synch_downs{d}),length(lags));
    mp_dtrig_lf3s_mat = zeros(length(synch_downs{d}),length(lags));
    mp_dtrig_lf2_mat = zeros(length(synch_downs{d}),length(lags));

    lf8_dtrig_mp_mat = zeros(length(synch_downs8{d}),length(lags));
    lf8_dtrig_lf8_mat = zeros(length(synch_downs8{d}),length(lags));
    lf8_dtrig_lf3_mat =zeros(length(synch_downs8{d}),length(lags));
    lf8_dtrig_lf3s_mat =zeros(length(synch_downs8{d}),length(lags));
    lf8_dtrig_lf2_mat = zeros(length(synch_downs8{d}),length(lags));


    %calculate mp utrigs
    for i = 1:length(synch_ups{d})

        if up_trans{d}(synch_ups{d}(i)) > backlag && ...
                length(down_w) - up_trans{d}(synch_ups{d}(i)) > forwardlag

            mp_utrig_mp_mat(i,:) = down_w(up_trans{d}(synch_ups{d}(i))-backlag:...
                up_trans{d}(synch_ups{d}(i))+forwardlag);
            mp_utrig_lf8_mat(i,:) = down_8(up_trans{d}(synch_ups{d}(i))-backlag:...
                up_trans{d}(synch_ups{d}(i))+forwardlag);
            mp_utrig_lf3_mat(i,:) = down_3(up_trans{d}(synch_ups{d}(i))-backlag:...
                up_trans{d}(synch_ups{d}(i))+forwardlag);
            mp_utrig_lf3s_mat(i,:) = down_3s(up_trans{d}(synch_ups{d}(i))-backlag:...
                up_trans{d}(synch_ups{d}(i))+forwardlag);
            mp_utrig_lf2_mat(i,:) = down_2(up_trans{d}(synch_ups{d}(i))-backlag:...
                up_trans{d}(synch_ups{d}(i))+forwardlag);

        else

            mp_utrig_mp_mat(i,:) = nan;
            mp_utrig_lf8_mat(i,:) = nan;
            mp_utrig_lf3_mat(i,:) = nan;
            mp_utrig_lf3s_mat(i,:) = nan;
            mp_utrig_lf2_mat(i,:) = nan;

        end

    end

    %calculate mp dtrigs
    for i = 1:length(synch_downs{d})

        if down_trans{d}(synch_downs{d}(i)) > backlag && ...
                length(down_w) - down_trans{d}(synch_downs{d}(i)) > forwardlag

            mp_dtrig_mp_mat(i,:) = down_w(down_trans{d}(synch_downs{d}(i))-backlag:...
                down_trans{d}(synch_downs{d}(i))+forwardlag);
            mp_dtrig_lf8_mat(i,:) = down_8(down_trans{d}(synch_downs{d}(i))-backlag:...
                down_trans{d}(synch_downs{d}(i))+forwardlag);
            mp_dtrig_lf3_mat(i,:) = down_3(down_trans{d}(synch_downs{d}(i))-backlag:...
                down_trans{d}(synch_downs{d}(i))+forwardlag);
            mp_dtrig_lf3s_mat(i,:) = down_3(down_trans{d}(synch_downs{d}(i))-backlag:...
                down_trans{d}(synch_downs{d}(i))+forwardlag);
            mp_dtrig_lf2_mat(i,:) = down_2(down_trans{d}(synch_downs{d}(i))-backlag:...
                down_trans{d}(synch_downs{d}(i))+forwardlag);

        else

            mp_dtrig_mp_mat(i,:) = nan;
            mp_dtrig_lf8_mat(i,:) = nan;
            mp_dtrig_lf3_mat(i,:) = nan;
            mp_dtrig_lf3s_mat(i,:) = nan;
            mp_dtrig_lf2_mat(i,:) = nan;

        end

    end

    %calculate lfp utrigs
    for i = 1:length(synch_ups8{d})

        if up_trans8{d}(synch_ups8{d}(i)) > backlag && ...
                length(down_w) - up_trans8{d}(synch_ups8{d}(i)) > forwardlag
            %
            lf8_utrig_mp_mat(i,:) = down_w(up_trans8{d}(synch_ups8{d}(i))-backlag:...
                up_trans8{d}(synch_ups8{d}(i))+forwardlag);
            lf8_utrig_lf8_mat(i,:) = down_8(up_trans8{d}(synch_ups8{d}(i))-backlag:...
                up_trans8{d}(synch_ups8{d}(i))+forwardlag);
            lf8_utrig_lf3_mat(i,:) = down_3(up_trans8{d}(synch_ups8{d}(i))-backlag:...
                up_trans8{d}(synch_ups8{d}(i))+forwardlag);
            lf8_utrig_lf3s_mat(i,:) = down_3(up_trans8{d}(synch_ups8{d}(i))-backlag:...
                up_trans8{d}(synch_ups8{d}(i))+forwardlag);
            lf8_utrig_lf2_mat(i,:) = down_2(up_trans8{d}(synch_ups8{d}(i))-backlag:...
                up_trans8{d}(synch_ups8{d}(i))+forwardlag);

        else

            lf8_utrig_mp_mat(i,:) = nan;
            lf8_utrig_lf8_mat(i,:) = nan;
            lf8_utrig_lf3_mat(i,:) = nan;
            lf8_utrig_lf3s_mat(i,:) = nan;
            lf8_utrig_lf2_mat(i,:) = nan;

        end

    end

    %calculate mp dtrigs
    for i = 1:length(synch_downs8{d})

        if down_trans8{d}(synch_downs8{d}(i)) > backlag && ...
                length(down_w) - down_trans8{d}(synch_downs8{d}(i)) > forwardlag

            lf8_dtrig_mp_mat(i,:) = down_w(down_trans8{d}(synch_downs8{d}(i))-backlag:...
                down_trans8{d}(synch_downs8{d}(i))+forwardlag);
            lf8_dtrig_lf8_mat(i,:) = down_8(down_trans8{d}(synch_downs8{d}(i))-backlag:...
                down_trans8{d}(synch_downs8{d}(i))+forwardlag);
            lf8_dtrig_lf3_mat(i,:) = down_3(down_trans8{d}(synch_downs8{d}(i))-backlag:...
                down_trans8{d}(synch_downs8{d}(i))+forwardlag);
            lf8_dtrig_lf3s_mat(i,:) = down_3(down_trans8{d}(synch_downs8{d}(i))-backlag:...
                down_trans8{d}(synch_downs8{d}(i))+forwardlag);
            lf8_dtrig_lf2_mat(i,:) = down_2(down_trans8{d}(synch_downs8{d}(i))-backlag:...
                down_trans8{d}(synch_downs8{d}(i))+forwardlag);

        else

            lf8_dtrig_mp_mat(i,:) = nan;
            lf8_dtrig_lf8_mat(i,:) = nan;
            lf8_dtrig_lf3_mat(i,:) = nan;
            lf8_dtrig_lf3s_mat(i,:) = nan;
            lf8_dtrig_lf2_mat(i,:) = nan;

        end

    end

    
  
%% calculate conditional variances
mp_utrig_mp_var(d,:) = nanvar(mp_utrig_mp_mat);
mp_dtrig_mp_var(d,:) = nanvar(mp_dtrig_mp_mat);
mp_utrig_lf8_var(d,:) = nanvar(mp_utrig_lf8_mat);
mp_dtrig_lf8_var(d,:) = nanvar(mp_dtrig_lf8_mat);
mp_utrig_lf3_var(d,:) = nanvar(mp_utrig_lf3_mat);
mp_dtrig_lf3_var(d,:) = nanvar(mp_dtrig_lf3_mat);
mp_utrig_lf2_var(d,:) = nanvar(mp_utrig_lf2_mat);
mp_dtrig_lf2_var(d,:) = nanvar(mp_dtrig_lf2_mat);
mp_utrig_lf3s_var(d,:) = nanvar(mp_utrig_lf3s_mat);
mp_dtrig_lf3s_var(d,:) = nanvar(mp_dtrig_lf3s_mat);

lf8_utrig_mp_var(d,:) = nanvar(lf8_utrig_mp_mat);
lf8_dtrig_mp_var(d,:) = nanvar(lf8_dtrig_mp_mat);
lf8_utrig_lf8_var(d,:) = nanvar(lf8_utrig_lf8_mat);
lf8_dtrig_lf8_var(d,:) = nanvar(lf8_dtrig_lf8_mat);
lf8_utrig_lf3_var(d,:) = nanvar(lf8_utrig_lf3_mat);
lf8_dtrig_lf3_var(d,:) = nanvar(lf8_dtrig_lf3_mat);
lf8_utrig_lf2_var(d,:) = nanvar(lf8_utrig_lf2_mat);
lf8_dtrig_lf2_var(d,:) = nanvar(lf8_dtrig_lf2_mat);
lf8_utrig_lf3s_var(d,:) = nanvar(lf8_utrig_lf3s_mat);
lf8_dtrig_lf3s_var(d,:) = nanvar(lf8_dtrig_lf3s_mat);

%% calculate conditional means
mp_utrig_mp_mean(d,:) = nanmean(mp_utrig_mp_mat);
mp_dtrig_mp_mean(d,:) = nanmean(mp_dtrig_mp_mat);
mp_utrig_lf8_mean(d,:) = nanmean(mp_utrig_lf8_mat);
mp_dtrig_lf8_mean(d,:) = nanmean(mp_dtrig_lf8_mat);
mp_utrig_lf3_mean(d,:) = nanmean(mp_utrig_lf3_mat);
mp_dtrig_lf3_mean(d,:) = nanmean(mp_dtrig_lf3_mat);
mp_utrig_lf2_mean(d,:) = nanmean(mp_utrig_lf2_mat);
mp_dtrig_lf2_mean(d,:) = nanmean(mp_dtrig_lf2_mat);
mp_utrig_lf3s_mean(d,:) = nanmean(mp_utrig_lf3s_mat);
mp_dtrig_lf3s_mean(d,:) = nanmean(mp_dtrig_lf3s_mat);

lf8_utrig_mp_mean(d,:) = nanmean(lf8_utrig_mp_mat);
lf8_dtrig_mp_mean(d,:) = nanmean(lf8_dtrig_mp_mat);
lf8_utrig_lf8_mean(d,:) = nanmean(lf8_utrig_lf8_mat);
lf8_dtrig_lf8_mean(d,:) = nanmean(lf8_dtrig_lf8_mat);
lf8_utrig_lf3_mean(d,:) = nanmean(lf8_utrig_lf3_mat);
lf8_dtrig_lf3_mean(d,:) = nanmean(lf8_dtrig_lf3_mat);
lf8_utrig_lf2_mean(d,:) = nanmean(lf8_utrig_lf2_mat);
lf8_dtrig_lf2_mean(d,:) = nanmean(lf8_dtrig_lf2_mat);
lf8_utrig_lf3s_mean(d,:) = nanmean(lf8_utrig_lf3s_mat);
lf8_dtrig_lf3s_mean(d,:) = nanmean(lf8_dtrig_lf3s_mat);

%% calculate marginal distributinos
    mp_marg = gpkde(down_w,-1,[le;re;nbins])';
    lf8_marg = gpkde(down_8,-1,[le;re;nbins])';
    lf3_marg = gpkde(down_3,-1,[le;re;nbins])';
    lf3s_marg = gpkde(down_3s,-1,[le;re;nbins])';
    lf2_marg = gpkde(down_2,-1,[le;re;nbins])';


%% calculate conditional density matrix estimates
    mp_utrig_mp_dens = zeros(length(lags),nbins);
    mp_dtrig_mp_dens = zeros(length(lags),nbins);
    mp_utrig_lf8_dens = zeros(length(lags),nbins);
    mp_dtrig_lf8_dens = zeros(length(lags),nbins);
    mp_utrig_lf3_dens = zeros(length(lags),nbins);
    mp_dtrig_lf3_dens = zeros(length(lags),nbins);
    mp_utrig_lf2_dens = zeros(length(lags),nbins);
    mp_dtrig_lf2_dens = zeros(length(lags),nbins);
    mp_utrig_lf3s_dens = zeros(length(lags),nbins);
    mp_dtrig_lf3s_dens = zeros(length(lags),nbins);

    lf8_utrig_mp_dens = zeros(length(lags),nbins);
    lf8_dtrig_mp_dens = zeros(length(lags),nbins);
    lf8_utrig_lf8_dens = zeros(length(lags),nbins);
    lf8_dtrig_lf8_dens = zeros(length(lags),nbins);
    lf8_utrig_lf3_dens = zeros(length(lags),nbins);
    lf8_dtrig_lf3_dens = zeros(length(lags),nbins);
    lf8_utrig_lf2_dens = zeros(length(lags),nbins);
    lf8_dtrig_lf2_dens = zeros(length(lags),nbins);
    lf8_utrig_lf3s_dens = zeros(length(lags),nbins);
    lf8_dtrig_lf3s_dens = zeros(length(lags),nbins);

    mp_utrig_mp_kl(d,:) = zeros(length(lags),1);
    mp_dtrig_mp_kl(d,:) = zeros(length(lags),1);
    mp_utrig_lf8_kl(d,:) = zeros(length(lags),1);
    mp_utrig_lf8_kl(d,:) = zeros(length(lags),1);
    mp_utrig_lf3_kl(d,:) = zeros(length(lags),1);
    mp_utrig_lf3_kl(d,:) = zeros(length(lags),1);
    mp_utrig_lf3s_kl(d,:) = zeros(length(lags),1);
    mp_utrig_lf3s_kl(d,:) = zeros(length(lags),1);
    mp_utrig_lf2_kl(d,:) = zeros(length(lags),1);
    mp_utrig_lf2_kl(d,:) = zeros(length(lags),1);

    lf8_utrig_mp_kl(d,:) = zeros(length(lags),1);
    lf8_dtrig_mp_kl(d,:) = zeros(length(lags),1);
    lf8_utrig_lf8_kl(d,:) = zeros(length(lags),1);
    lf8_utrig_lf8_kl(d,:) = zeros(length(lags),1);
    lf8_utrig_lf3_kl(d,:) = zeros(length(lags),1);
    lf8_utrig_lf3_kl(d,:) = zeros(length(lags),1);
    lf8_utrig_lf3s_kl(d,:) = zeros(length(lags),1);
    lf8_utrig_lf3s_kl(d,:) = zeros(length(lags),1);
    lf8_utrig_lf2_kl(d,:) = zeros(length(lags),1);
    lf8_utrig_lf2_kl(d,:) = zeros(length(lags),1);


    for i = 1:length(lags)

        mp_utrig_mp_dens(i,:) = gpkde(mp_utrig_mp_mat(~isnan(mp_utrig_mp_mat(:,i)),i),-1,[le; re; nbins]);
        mp_utrig_mp_kl(d,i) = k_l_divergence(mp_utrig_mp_dens(i,:),mp_marg);

        mp_dtrig_mp_dens(i,:) = gpkde(mp_dtrig_mp_mat(~isnan(mp_dtrig_mp_mat(:,i)),i),-1,[le; re; nbins]);
        mp_dtrig_mp_kl(d,i) = k_l_divergence(mp_dtrig_mp_dens(i,:),mp_marg);

        mp_utrig_lf8_dens(i,:) = gpkde(mp_utrig_lf8_mat(~isnan(mp_utrig_lf8_mat(:,i)),i),-1,[le; re; nbins]);
        mp_utrig_lf8_kl(d,i) = k_l_divergence(mp_utrig_lf8_dens(i,:),lf8_marg);

        mp_dtrig_lf8_dens(i,:) = gpkde(mp_dtrig_lf8_mat(~isnan(mp_dtrig_lf8_mat(:,i)),i),-1,[le; re; nbins]);
        mp_dtrig_lf8_kl(d,i) = k_l_divergence(mp_dtrig_lf8_dens(i,:),lf8_marg);

        mp_utrig_lf3_dens(i,:) = gpkde(mp_utrig_lf3_mat(~isnan(mp_utrig_lf3_mat(:,i)),i),-1,[le; re; nbins]);
        mp_utrig_lf3_kl(d,i) = k_l_divergence(mp_utrig_lf3_dens(i,:),lf3_marg);

        mp_dtrig_lf3_dens(i,:) = gpkde(mp_dtrig_lf3_mat(~isnan(mp_dtrig_lf3_mat(:,i)),i),-1,[le; re; nbins]);
        mp_dtrig_lf3_kl(d,i) = k_l_divergence(mp_dtrig_lf3_dens(i,:),lf3_marg);

        mp_utrig_lf2_dens(i,:) = gpkde(mp_utrig_lf2_mat(~isnan(mp_utrig_lf2_mat(:,i)),i),-1,[le; re; nbins]);
        mp_utrig_lf2_kl(d,i) = k_l_divergence(mp_utrig_lf2_dens(i,:),lf2_marg);

        mp_dtrig_lf2_dens(i,:) = gpkde(mp_dtrig_lf2_mat(~isnan(mp_dtrig_lf2_mat(:,i)),i),-1,[le; re; nbins]);
        mp_dtrig_lf2_kl(d,i) = k_l_divergence(mp_dtrig_lf2_dens(i,:),lf2_marg);

        mp_utrig_lf3s_dens(i,:) = gpkde(mp_utrig_lf3s_mat(~isnan(mp_utrig_lf3s_mat(:,i)),i),-1,[le; re; nbins]);
        mp_utrig_lf3s_kl(d,i) = k_l_divergence(mp_utrig_lf3s_dens(i,:),lf3s_marg);

        mp_dtrig_lf3s_dens(i,:) = gpkde(mp_dtrig_lf3s_mat(~isnan(mp_dtrig_lf3s_mat(:,i)),i),-1,[le; re; nbins]);
        mp_dtrig_lf3s_kl(d,i) = k_l_divergence(mp_dtrig_lf3s_dens(i,:),lf3s_marg);


        lf8_utrig_mp_dens(i,:) = gpkde(lf8_utrig_mp_mat(~isnan(lf8_utrig_mp_mat(:,i)),i),-1,[le; re; nbins]);
        lf8_utrig_mp_kl(d,i) = k_l_divergence(lf8_utrig_mp_dens(i,:),mp_marg);

        lf8_dtrig_mp_dens(i,:) = gpkde(lf8_dtrig_mp_mat(~isnan(lf8_dtrig_mp_mat(:,i)),i),-1,[le; re; nbins]);
        lf8_dtrig_mp_kl(d,i) = k_l_divergence(lf8_dtrig_mp_dens(i,:),mp_marg);

        lf8_utrig_lf8_dens(i,:) = gpkde(lf8_utrig_lf8_mat(~isnan(lf8_utrig_lf8_mat(:,i)),i),-1,[le; re; nbins]);
        lf8_utrig_lf8_kl(d,i) = k_l_divergence(lf8_utrig_lf8_dens(i,:),lf8_marg);

        lf8_dtrig_lf8_dens(i,:) = gpkde(lf8_dtrig_lf8_mat(~isnan(lf8_dtrig_lf8_mat(:,i)),i),-1,[le; re; nbins]);
        lf8_dtrig_lf8_kl(d,i) = k_l_divergence(lf8_dtrig_lf8_dens(i,:),lf8_marg);

        lf8_utrig_lf3_dens(i,:) = gpkde(lf8_utrig_lf3_mat(~isnan(lf8_utrig_lf3_mat(:,i)),i),-1,[le; re; nbins]);
        lf8_utrig_lf3_kl(d,i) = k_l_divergence(lf8_utrig_lf3_dens(i,:),lf3_marg);

        lf8_dtrig_lf3_dens(i,:) = gpkde(lf8_dtrig_lf3_mat(~isnan(lf8_dtrig_lf3_mat(:,i)),i),-1,[le; re; nbins]);
        lf8_dtrig_lf3_kl(d,i) = k_l_divergence(lf8_dtrig_lf3_dens(i,:),lf3_marg);

        lf8_utrig_lf2_dens(i,:) = gpkde(lf8_utrig_lf2_mat(~isnan(lf8_utrig_lf2_mat(:,i)),i),-1,[le; re; nbins]);
        lf8_utrig_lf2_kl(d,i) = k_l_divergence(lf8_utrig_lf2_dens(i,:),lf2_marg);

        lf8_dtrig_lf2_dens(i,:) = gpkde(lf8_dtrig_lf2_mat(~isnan(lf8_dtrig_lf2_mat(:,i)),i),-1,[le; re; nbins]);
        lf8_dtrig_lf2_kl(d,i) = k_l_divergence(lf8_dtrig_lf2_dens(i,:),lf2_marg);

        lf8_utrig_lf3s_dens(i,:) = gpkde(lf8_utrig_lf3s_mat(~isnan(lf8_utrig_lf3s_mat(:,i)),i),-1,[le; re; nbins]);
        lf8_utrig_lf3s_kl(d,i) = k_l_divergence(lf8_utrig_lf3s_dens(i,:),lf3s_marg);

        lf8_dtrig_lf3s_dens(i,:) = gpkde(lf8_dtrig_lf3s_mat(~isnan(lf8_dtrig_lf3s_mat(:,i)),i),-1,[le; re; nbins]);
        lf8_dtrig_lf3s_kl(d,i) = k_l_divergence(lf8_dtrig_lf3s_dens(i,:),lf3s_marg);

    end

%% plot results

    subplot(2,1,1)
    pcolor(lags,linspace(le,re,nbins),log(mp_utrig_mp_dens)');shading flat
    hold on
    plot(lags,mp_utrig_mp_mean(d,:),'k')
    title('MP Up Triggered')
    caxis([-5 1])
    subplot(2,1,2)
    pcolor(lags,linspace(le,re,nbins),log(mp_dtrig_mp_dens)');shading flat
    hold on
        plot(lags,mp_dtrig_mp_mean(d,:),'k')
    title('MP Down Triggered')
    caxis([-5 1])
    t_names = ['C:\WC_Germany\overall_calcs\trig_cond_dens_theta\mp_trig_mp_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    subplot(2,1,1)
    pcolor(lags,linspace(le,re,nbins),log(mp_utrig_lf8_dens)');shading flat
    hold on
    plot(lags,mp_utrig_lf8_mean(d,:),'k')
    caxis([-5 1])
    title('MP Up Triggered')
    subplot(2,1,2)
    pcolor(lags,linspace(le,re,nbins),log(mp_dtrig_lf8_dens)');shading flat
    hold on
    plot(lags,mp_dtrig_lf8_mean(d,:),'k')
    caxis([-5 1])
    title('MP Down Triggered')
    t_names = ['C:\WC_Germany\overall_calcs\trig_cond_dens_theta\mp_trig_lf8_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    subplot(2,1,1)
    pcolor(lags,linspace(le,re,nbins),log(mp_utrig_lf3_dens)');shading flat
    hold on
    plot(lags,mp_utrig_lf3_mean(d,:),'k')
    caxis([-5 1])
    title('MP Up Triggered')
    subplot(2,1,2)
    pcolor(lags,linspace(le,re,nbins),log(mp_dtrig_lf3_dens)');shading flat
    caxis([-5 1])
    hold on
    plot(lags,mp_dtrig_lf3_mean(d,:),'k')
    title('MP Down Triggered')
    t_names = ['C:\WC_Germany\overall_calcs\trig_cond_dens_theta\mp_trig_lf3_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    subplot(2,1,1)
    pcolor(lags,linspace(le,re,nbins),log(mp_utrig_lf2_dens)');shading flat
    hold on
    plot(lags,mp_utrig_lf2_mean(d,:),'k')
    caxis([-5 1])
    title('MP Up Triggered')
    subplot(2,1,2)
    pcolor(lags,linspace(le,re,nbins),log(mp_dtrig_lf2_dens)');shading flat
    hold on
    plot(lags,mp_dtrig_lf2_mean(d,:),'k')
    caxis([-5 1])
    title('MP Down Triggered')
    t_names = ['C:\WC_Germany\overall_calcs\trig_cond_dens_theta\mp_trig_lf2_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    subplot(2,1,1)
    pcolor(lags,linspace(le,re,nbins),log(mp_utrig_lf3s_dens)');shading flat
    hold on
    plot(lags,mp_utrig_lf3s_mean(d,:),'k')
    caxis([-5 1])
    title('MP Up Triggered')
    subplot(2,1,2)
    pcolor(lags,linspace(le,re,nbins),log(mp_dtrig_lf3s_dens)');shading flat
    hold on
    plot(lags,mp_dtrig_lf3s_mean(d,:),'k')
    caxis([-5 1])
    title('MP Down Triggered')
    t_names = ['C:\WC_Germany\overall_calcs\trig_cond_dens_theta\mp_trig_lf3s_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    subplot(2,1,1)
    pcolor(lags,linspace(le,re,nbins),log(lf8_utrig_mp_dens)');shading flat
    hold on
    plot(lags,lf8_utrig_mp_mean(d,:),'k')
    caxis([-5 1])
    title('LFP Up Triggered')
    subplot(2,1,2)
    pcolor(lags,linspace(le,re,nbins),log(lf8_dtrig_mp_dens)');shading flat
    caxis([-5 1])
    hold on
    plot(lags,lf8_dtrig_mp_mean(d,:),'k')
    title('LFP Down Triggered')
    t_names = ['C:\WC_Germany\overall_calcs\trig_cond_dens_theta\lf8_trig_mp_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    subplot(2,1,1)
    pcolor(lags,linspace(le,re,nbins),log(lf8_utrig_lf8_dens)');shading flat
    hold on
    plot(lags,lf8_utrig_lf8_mean(d,:),'k')
    caxis([-5 1])
    title('LFP Up Triggered')
    subplot(2,1,2)
    pcolor(lags,linspace(le,re,nbins),log(lf8_dtrig_lf8_dens)');shading flat
    caxis([-5 1])
    hold on
    plot(lags,lf8_dtrig_lf8_mean(d,:),'k')
    title('LFP Down Triggered')
    t_names = ['C:\WC_Germany\overall_calcs\trig_cond_dens_theta\lf8_trig_lf8_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    subplot(2,1,1)
    pcolor(lags,linspace(le,re,nbins),log(lf8_utrig_lf3_dens)');shading flat
    hold on
    plot(lags,lf8_utrig_lf3_mean(d,:),'k')
    caxis([-5 1])
    title('LFP Up Triggered')
    subplot(2,1,2)
    pcolor(lags,linspace(le,re,nbins),log(lf8_dtrig_lf3_dens)');shading flat
    caxis([-5 1])
    hold on
    plot(lags,lf8_dtrig_lf3_mean(d,:),'k')
    title('LFP Down Triggered')
    t_names = ['C:\WC_Germany\overall_calcs\trig_cond_dens_theta\lf8_trig_lf3_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    subplot(2,1,1)
    pcolor(lags,linspace(le,re,nbins),log(lf8_utrig_lf2_dens)');shading flat
    caxis([-5 1])
    hold on
    plot(lags,lf8_utrig_lf2_mean(d,:),'k')
    title('LFP Up Triggered')
    subplot(2,1,2)
    pcolor(lags,linspace(le,re,nbins),log(lf8_dtrig_lf2_dens)');shading flat
    caxis([-5 1])
    hold on
    plot(lags,lf8_dtrig_lf2_mean(d,:),'k')
    title('LFP Down Triggered')
    t_names = ['C:\WC_Germany\overall_calcs\trig_cond_dens_theta\lf8_trig_lf2_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    subplot(2,1,1)
    pcolor(lags,linspace(le,re,nbins),log(lf8_utrig_lf3s_dens)');shading flat
    caxis([-5 1])
    hold on
    plot(lags,lf8_utrig_lf3s_mean(d,:),'k')
    title('LFP Up Triggered')
    subplot(2,1,2)
    pcolor(lags,linspace(le,re,nbins),log(lf8_dtrig_lf3s_dens)');shading flat
    caxis([-5 1])
    hold on
    plot(lags,lf8_dtrig_lf3s_mean(d,:),'k')
    title('LFP Down Triggered')
    t_names = ['C:\WC_Germany\overall_calcs\trig_cond_dens_theta\lf8_trig_lf3s_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    plot(lags,mp_utrig_mp_kl(d,:))
    hold on
    plot(lags,mp_utrig_lf8_kl(d,:),'r')
    plot(lags,mp_utrig_lf3_kl(d,:),'g')
    plot(lags,mp_utrig_lf3s_kl(d,:),'m')
    plot(lags,mp_utrig_lf2_kl(d,:),'k')
    legend('MP','LF8','LF3','LF3s','LF2')
    t_names = ['C:\WC_Germany\overall_calcs\trig_cond_dens_theta\mp_utrig_kl_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    plot(lags,mp_dtrig_mp_kl(d,:))
    hold on
    plot(lags,mp_dtrig_lf8_kl(d,:),'r')
    plot(lags,mp_dtrig_lf3_kl(d,:),'g')
    plot(lags,mp_dtrig_lf3s_kl(d,:),'m')
    plot(lags,mp_dtrig_lf2_kl(d,:),'k')
    legend('MP','LF8','LF3','LF3s','LF2')

    t_names = ['C:\WC_Germany\overall_calcs\trig_cond_dens_theta\mp_dtrig_kl_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    plot(lags,lf8_utrig_mp_kl(d,:))
    hold on
    plot(lags,lf8_utrig_lf8_kl(d,:),'r')
    plot(lags,lf8_utrig_lf3_kl(d,:),'g')
    plot(lags,lf8_utrig_lf3s_kl(d,:),'m')
    plot(lags,lf8_utrig_lf2_kl(d,:),'k')
    legend('MP','LF8','LF3','LF3s','LF2')

    t_names = ['C:\WC_Germany\overall_calcs\trig_cond_dens_theta\lf8_utrig_kl_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    plot(lags,lf8_dtrig_mp_kl(d,:))
    hold on
    plot(lags,lf8_dtrig_lf8_kl(d,:),'r')
    plot(lags,lf8_dtrig_lf3_kl(d,:),'g')
    plot(lags,lf8_dtrig_lf3s_kl(d,:),'m')
    plot(lags,lf8_dtrig_lf2_kl(d,:),'k')
    legend('MP','LF8','LF3','LF3s','LF2')

    t_names = ['C:\WC_Germany\overall_calcs\trig_cond_dens_theta\lf8_dtrig_kl_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    
    
        plot(lags,mp_utrig_mp_var(d,:))
    hold on
    plot(lags,mp_utrig_lf8_var(d,:),'r')
    plot(lags,mp_utrig_lf3_var(d,:),'g')
    plot(lags,mp_utrig_lf3s_var(d,:),'m')
    plot(lags,mp_utrig_lf2_var(d,:),'k')
    legend('MP','LF8','LF3','LF3s','LF2')
    t_names = ['C:\WC_Germany\overall_calcs\trig_cond_dens_theta\mp_utrig_var_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    plot(lags,mp_dtrig_mp_var(d,:))
    hold on
    plot(lags,mp_dtrig_lf8_var(d,:),'r')
    plot(lags,mp_dtrig_lf3_var(d,:),'g')
    plot(lags,mp_dtrig_lf3s_var(d,:),'m')
    plot(lags,mp_dtrig_lf2_var(d,:),'k')
    legend('MP','LF8','LF3','LF3s','LF2')

    t_names = ['C:\WC_Germany\overall_calcs\trig_cond_dens_theta\mp_dtrig_var_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    plot(lags,lf8_utrig_mp_var(d,:))
    hold on
    plot(lags,lf8_utrig_lf8_var(d,:),'r')
    plot(lags,lf8_utrig_lf3_var(d,:),'g')
    plot(lags,lf8_utrig_lf3s_var(d,:),'m')
    plot(lags,lf8_utrig_lf2_var(d,:),'k')
    legend('MP','LF8','LF3','LF3s','LF2')

    t_names = ['C:\WC_Germany\overall_calcs\trig_cond_dens_theta\lf8_utrig_var_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    plot(lags,lf8_dtrig_mp_var(d,:))
    hold on
    plot(lags,lf8_dtrig_lf8_var(d,:),'r')
    plot(lags,lf8_dtrig_lf3_var(d,:),'g')
    plot(lags,lf8_dtrig_lf3s_var(d,:),'m')
    plot(lags,lf8_dtrig_lf2_var(d,:),'k')
    legend('MP','LF8','LF3','LF3s','LF2')

    t_names = ['C:\WC_Germany\overall_calcs\trig_cond_dens_theta\lf8_dtrig_var_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close

    
end

save C:\WC_Germany\overall_calcs\trig_cond_dens_theta\cond_dens_data *dens lags *kl *var *mean