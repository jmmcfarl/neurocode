
clear all
close all
%%
load G:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('G:\Code\WC_anal\general\')
addpath('G:\WC_Germany\Overall_EC\')
addpath('G:\Code\Chronux\spectral_analysis\continuous\')
addpath('G:\WC_Germany\hsmm_state_detection\\')

drive_letter = 'G';
cd G:\WC_Germany\overall_EC
load corresponding_lfp_state_data
%%
dsf = 8;
Fsd = 2016/dsf;
niqf = 2016/2;
[b,a] = butter(2,[0.05/niqf 90/niqf]);
[b2,a2] = butter(2,[0.05/niqf 2/niqf]);

params.Fs = Fsd;
params.err = [2 0.05];
params.fpass = [4 80];
movingwin = [1 1];
params.tapers = [2 3];

amp_grid = linspace(-4,4,200);

for d = 1:109
    
    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    load used_data lf8 lf3 wcv_minus_spike
    
    wcv_d = filtfilt(b,a,wcv_minus_spike);
    wcv_d = zscore(downsample(wcv_d,dsf));
    lf8_d = filtfilt(b,a,lf8);
    lf8_d = zscore(downsample(lf8_d,dsf));
    lf3_d = filtfilt(b,a,lf3);
    lf3_d = zscore(downsample(lf3_d,dsf));
    
    lf3_uds = filtfilt(b2,a2,lf3);
    lf3_uds = zscore(downsample(lf3_uds,dsf));
    
    t_axis = (1:length(wcv_d))/Fsd;
    
    %% extract up and down transition times for MP and LF8
    load ec_hmm_state_seq
    load ec_hmm_state_seq8
    mp_state_seq_c =  hmm_bbstate_seq;
    lf8_state_seq_c = hmm_bbstate_seq8;
    [new_mp_seg_inds] = round(resample_uds_seg_inds(hmm.UDS_segs,50.4,Fsd,length(t_axis)));
    [new_lf8_seg_inds] = round(resample_uds_seg_inds(hmm8.UDS_segs,50.4,Fsd,length(t_axis)));
    
    mp_state_seq = nan(size(t_axis));
    lf8_state_seq = nan(size(t_axis));
    mp_utrans = [];
    mp_dtrans = [];
    for n = 1:hmm.Nsegs
        mp_state_seq(new_mp_seg_inds(n,1):new_mp_seg_inds(n,2)) = mp_state_seq_c{n};
        cur_mp_utrans = new_mp_seg_inds(n,1) + find(mp_state_seq_c{n}(1:end-1) == 1 & mp_state_seq_c{n}(2:end) == 2);
        cur_mp_dtrans = new_mp_seg_inds(n,1) + find(mp_state_seq_c{n}(1:end-1) == 2 & mp_state_seq_c{n}(2:end) == 1);
        cur_mp_dtrans(cur_mp_dtrans < cur_mp_utrans(1)) = [];
        cur_mp_utrans(cur_mp_utrans > cur_mp_dtrans(end)) = [];
        mp_utrans = [mp_utrans; cur_mp_utrans];
        mp_dtrans = [mp_dtrans; cur_mp_dtrans];
    end
    
    lf8_utrans = [];
    lf8_dtrans = [];
    for n = 1:hmm8.Nsegs
        lf8_state_seq(new_lf8_seg_inds(n,1):new_lf8_seg_inds(n,2)) = lf8_state_seq_c{n};
        cur_lf8_utrans = new_lf8_seg_inds(n,1) + find(lf8_state_seq_c{n}(1:end-1) == 1 & lf8_state_seq_c{n}(2:end) == 2);
        cur_lf8_dtrans = new_lf8_seg_inds(n,1) + find(lf8_state_seq_c{n}(1:end-1) == 2 & lf8_state_seq_c{n}(2:end) == 1);
        cur_lf8_dtrans(cur_lf8_dtrans < cur_lf8_utrans(1)) = [];
        cur_lf8_utrans(cur_lf8_utrans > cur_lf8_dtrans(end)) = [];
        lf8_utrans = [lf8_utrans; cur_lf8_utrans];
        lf8_dtrans = [lf8_dtrans; cur_lf8_dtrans];
    end
    
    mp_up_lf8_down = nan(size(lf8_state_seq));
    mp_up_lf8_down(mp_state_seq==2 & lf8_state_seq==1) = 2;
    mp_up_lf8_down(mp_state_seq==1 | lf8_state_seq == 2) = 1;
    wupldown_utrans = find(mp_up_lf8_down(1:end-1) == 1 & mp_up_lf8_down(2:end) == 2);
    wupldown_dtrans = nan(size(wupldown_utrans));
    for ii = 1:length(wupldown_utrans)
        cur = find(mp_up_lf8_down(wupldown_utrans(ii)+1:end) == 1 | isnan(mp_up_lf8_down(wupldown_utrans(ii)+1:end)),1,'first');
        if ~isempty(cur)
            wupldown_dtrans(ii) = wupldown_utrans(ii) + cur+1;
        else
            wupldown_dtrans(ii) = length(mp_up_lf8_down);
        end
    end
    long_enough = find(wupldown_dtrans-wupldown_utrans > round(Fsd*0.05));
    wupldown_dtrans = wupldown_dtrans(long_enough)';
    wupldown_utrans = wupldown_utrans(long_enough)';
    
    mp_up_lf8_up = nan(size(lf8_state_seq));
    mp_up_lf8_up(mp_state_seq==2 & lf8_state_seq==2) = 2;
    mp_up_lf8_up(mp_state_seq==1 | lf8_state_seq == 1) = 1;
    wuplup_utrans = find(mp_up_lf8_up(1:end-1) == 1 & mp_up_lf8_up(2:end) == 2);
    wuplup_dtrans = nan(size(wuplup_utrans));
    for ii = 1:length(wuplup_utrans)
        cur = find(mp_up_lf8_up(wuplup_utrans(ii)+1:end) == 1 | isnan(mp_up_lf8_up(wuplup_utrans(ii)+1:end)),1,'first');
        if ~isempty(cur)
            wuplup_dtrans(ii) = wuplup_utrans(ii) + cur+1;
        else
            wuplup_dtrans(ii) = length(mp_up_lf8_down);
        end
    end
    long_enough = find(wuplup_dtrans-wuplup_utrans > round(Fsd*0.05));
    wuplup_dtrans = wuplup_dtrans(long_enough)';
    wuplup_utrans = wuplup_utrans(long_enough)';
    
    mp_down_lf8_down = nan(size(lf8_state_seq));
    mp_down_lf8_down(mp_state_seq==1 & lf8_state_seq==1) = 2;
    mp_down_lf8_down(mp_state_seq==2 | lf8_state_seq == 2) = 1;
    wdownldown_utrans = find(mp_down_lf8_down(1:end-1) == 1 & mp_down_lf8_down(2:end) == 2);
    wdownldown_dtrans = nan(size(wdownldown_utrans));
    for ii = 1:length(wdownldown_utrans)
        cur = find(mp_down_lf8_down(wdownldown_utrans(ii)+1:end) == 1 | isnan(mp_down_lf8_down(wdownldown_utrans(ii)+1:end)),1,'first');
        if ~isempty(cur)
            wdownldown_dtrans(ii) = wdownldown_utrans(ii) + cur+1;
        else
            wdownldown_dtrans(ii) = length(mp_up_lf8_down);
        end
    end
    long_enough = find(wdownldown_dtrans-wdownldown_utrans > round(Fsd*0.05));
    wdownldown_utrans = wdownldown_utrans(long_enough)';
    wdownldown_dtrans = wdownldown_dtrans(long_enough)';
    
    mp_down_lf8_up = nan(size(lf8_state_seq));
    mp_down_lf8_up(mp_state_seq==1 & lf8_state_seq==2) = 2;
    mp_down_lf8_up(mp_state_seq==2 | lf8_state_seq == 1) = 1;
    wdownlup_utrans = find(mp_down_lf8_up(1:end-1) == 1 & mp_down_lf8_up(2:end) == 2);
    wdownlup_dtrans = nan(size(wdownlup_utrans));
    for ii = 1:length(wdownlup_utrans)
        cur = find(mp_down_lf8_up(wdownlup_utrans(ii)+1:end) == 1 | isnan(mp_down_lf8_up(wdownlup_utrans(ii)+1:end)),1,'first');
        if ~isempty(cur)
            wdownlup_dtrans(ii) = wdownlup_utrans(ii) + cur+1;
        else
            wdownlup_dtrans(ii) = length(mp_up_lf8_down);
        end
    end
    long_enough = find(wdownlup_dtrans-wdownlup_utrans > round(Fsd*0.05));
    wdownlup_utrans = wdownlup_utrans(long_enough)';
    wdownlup_dtrans = wdownlup_dtrans(long_enough)';
    
    %%
    used_data = find(~isnan(lf8_state_seq) & ~isnan(mp_state_seq));
    
    p_lf8up(d) = sum(lf8_state_seq(used_data)==2)/length(used_data);
    p_lf8down(d) = sum(lf8_state_seq(used_data)==1)/length(used_data);
    p_mpup(d) = sum(mp_state_seq(used_data)==2)/length(used_data);
    p_mpdown(d) = sum(mp_state_seq(used_data)==1)/length(used_data);
    p_lf8up_cond_mpdown(d) = sum(lf8_state_seq(used_data)==2 & mp_state_seq(used_data)==1)/length(used_data)/p_mpdown(d);
    p_lf8down_cond_mpdown(d) = sum(lf8_state_seq(used_data)==1 & mp_state_seq(used_data)==1)/length(used_data)/p_mpdown(d);
    p_lf8up_cond_mpup(d) = sum(lf8_state_seq(used_data)==2 & mp_state_seq(used_data)==2)/length(used_data)/p_mpup(d);
    p_lf8down_cond_mpup(d) = sum(lf8_state_seq(used_data)==1 & mp_state_seq(used_data)==2)/length(used_data)/p_mpup(d);
    p_mpup_cond_lf8down(d) = sum(mp_state_seq(used_data)==2 & lf8_state_seq(used_data)==1)/length(used_data)/p_lf8down(d);
    p_mpdown_cond_lf8down(d) = sum(mp_state_seq(used_data)==1 & lf8_state_seq(used_data)==1)/length(used_data)/p_lf8down(d);
    p_mpup_cond_lf8up(d) = sum(mp_state_seq(used_data)==2 & lf8_state_seq(used_data)==2)/length(used_data)/p_lf8up(d);
    p_mpdown_cond_lf8up(d) = sum(mp_state_seq(used_data)==1 & lf8_state_seq(used_data)==2)/length(used_data)/p_lf8up(d);
    
    
    [S3_mpup(d,:), f] = mtspectrumc_unequal_length_trials(lf3_d,movingwin,params,[mp_utrans(:) mp_dtrans(:)]);
    [S3_mpdown(d,:), f] = mtspectrumc_unequal_length_trials(lf3_d,movingwin,params,[mp_dtrans(1:end-1) mp_utrans(2:end)]);
    [S3_lf8up(d,:), f] = mtspectrumc_unequal_length_trials(lf3_d,movingwin,params,[lf8_utrans(:) lf8_dtrans(:)]);
    [S3_lf8down(d,:), f] = mtspectrumc_unequal_length_trials(lf3_d,movingwin,params,[lf8_dtrans(1:end-1) lf8_utrans(2:end)]);
    [S3_mpuplf8down(d,:), f] = mtspectrumc_unequal_length_trials(lf3_d,movingwin,params,[wupldown_utrans(:) wupldown_dtrans(:)]);
    [S3_mpuplf8up(d,:), f] = mtspectrumc_unequal_length_trials(lf3_d,movingwin,params,[wuplup_utrans(:) wuplup_dtrans(:)]);
    [S3_mpdownlf8down(d,:), f] = mtspectrumc_unequal_length_trials(lf3_d,movingwin,params,[wdownldown_utrans(:) wdownldown_dtrans(:)]);
    [S3_mpdownlf8up(d,:), f] = mtspectrumc_unequal_length_trials(lf3_d,movingwin,params,[wdownlup_utrans(:) wdownlup_dtrans(:)]);
    [Sw_mpuplf8down(d,:), f] = mtspectrumc_unequal_length_trials(wcv_d,movingwin,params,[wupldown_utrans(:) wupldown_dtrans(:)]);
    [Sw_mpuplf8up(d,:), f] = mtspectrumc_unequal_length_trials(wcv_d,movingwin,params,[wuplup_utrans(:) wuplup_dtrans(:)]);
    [Sw_mpdownlf8down(d,:), f]= mtspectrumc_unequal_length_trials(wcv_d,movingwin,params,[wdownldown_utrans(:) wdownldown_dtrans(:)]);
    [Sw_mpdownlf8up(d,:), f]= mtspectrumc_unequal_length_trials(wcv_d,movingwin,params,[wdownlup_utrans(:) wdownlup_dtrans(:)]);
    
    S3_mpdiff(d,:) = log10(S3_mpup(d,:)) - log10(S3_mpdown(d,:));
    S3_lf8diff(d,:) = log10(S3_lf8up(d,:)) - log10(S3_lf8down(d,:));
    
    lf3dist_mpup(d,:) = ksdensity(lf3_uds(mp_state_seq == 2),amp_grid);
    lf3dist_mpdown(d,:) = ksdensity(lf3_uds(mp_state_seq == 1),amp_grid);
    lf3dist_lf8up(d,:) = ksdensity(lf3_uds(lf8_state_seq == 2),amp_grid);
    lf3dist_lf8down(d,:) = ksdensity(lf3_uds(lf8_state_seq == 1),amp_grid);
    lf3dist_mpup_lf8up(d,:) = ksdensity(lf3_uds(mp_state_seq == 2 & lf8_state_seq == 2),amp_grid);
    lf3dist_mpup_lf8down(d,:) = ksdensity(lf3_uds(mp_state_seq == 2 & lf8_state_seq == 1),amp_grid);
    lf3dist_mpdown_lf8up(d,:) = ksdensity(lf3_uds(mp_state_seq == 1 & lf8_state_seq == 2),amp_grid);
    lf3dist_mpdown_lf8down(d,:) = ksdensity(lf3_uds(mp_state_seq == 1 & lf8_state_seq == 1),amp_grid);
    
    mp_lf3_sepinfo(d) = compute_kl_div(lf3_uds,mp_state_seq);
    lf8_lf3_sepinfo(d) = compute_kl_div(lf3_uds,lf8_state_seq);
    
    mp_lf3_info_cond_8up = compute_kl_div(lf3_uds(lf8_state_seq == 2),mp_state_seq(lf8_state_seq == 2));
    mp_lf3_info_cond_8down = compute_kl_div(lf3_uds(lf8_state_seq == 1),mp_state_seq(lf8_state_seq == 1));
    mp_lf3_sepinfo_cond8(d) = mp_lf3_info_cond_8up*p_lf8up(d) + mp_lf3_info_cond_8down*p_lf8down(d);
    lf8_lf3_info_cond_mup = compute_kl_div(lf3_uds(mp_state_seq == 2),lf8_state_seq(mp_state_seq == 2));
    lf8_lf3_info_cond_mdown = compute_kl_div(lf3_uds(mp_state_seq == 1),lf8_state_seq(mp_state_seq == 1));
    lf8_lf3_sepinfo_condm(d) = lf8_lf3_info_cond_mup*p_mpup(d) + lf8_lf3_info_cond_mdown*p_mpdown(d);
    
    
    %     figure
    %     plot(f,S3_mpdiff(d,:)), hold on
    %     plot(f,S3_lf8diff(d,:),'k')
    %
    %     figure
    %     plot(amp_grid,lf3dist_mpup(d,:)), hold on
    %     plot(amp_grid,lf3dist_mpdown(d,:),'g')
    %     plot(amp_grid,lf3dist_lf8up(d,:),'k')
    %     plot(amp_grid,lf3dist_lf8down(d,:),'r')
    %
    %     figure
    %     plot(amp_grid,lf3dist_mpup_lf8up(d,:)), hold on
    %     plot(amp_grid,lf3dist_mpup_lf8down(d,:),'g')
    %     plot(amp_grid,lf3dist_mpdown_lf8up(d,:),'k')
    %     plot(amp_grid,lf3dist_mpdown_lf8down(d,:),'r')
    
end

%%
cd G:\WC_Germany\overall_EC\
save overall_EC_state_dep_data mp_lf3_* lf8_lf3_* lf3dist_* S3_* Sw_* f amp_grid p_*

%%

load overall_EC_coherence_data
mec = find_struct_field_vals(sess_data,'region','MEC');
layer3 = find_struct_field_vals(sess_data,'layer','3');
layer2 = find_struct_field_vals(sess_data,'layer','2');
layer23 = find_struct_field_vals(sess_data,'layer','23');
lec = find_struct_field_vals(sess_data,'region','LEC');
l3mec = intersect(mec,layer3);
l2mec = intersect(mec,layer2);
l3lec = intersect(lec,layer3);
l3mec(24:end) = [];
l23mec = intersect(mec,layer23);
l23mec = unique([l2mec l3mec l23mec]);

frange = find(f_i > 0.2 & f_i < 0.6);
avg_P83 = mean(P83(:,frange),2);
correct_lf3phase = find(avg_P83 > -0.4);
l3mec_c = l3mec(avg_P83(l3mec) > -0.4);
l3mec_w = l3mec(avg_P83(l3mec) < -0.4);
l3lec_c = l3lec(avg_P83(l3lec) > -0.4);
l3lec_w = l3lec(avg_P83(l3lec) < -0.4);

%%
figure
h = errorbar(amp_grid,nanmean(lf3dist_mpdown_lf8down(l3mec_c,:)),nanstd(lf3dist_mpdown_lf8down(l3mec_c,:))/sqrt(length(l3mec_c)));
errorbar_tick(h,.01,'units');
hold on
h = errorbar(amp_grid,nanmean(lf3dist_mpup_lf8down(l3mec_c,:)),nanstd(lf3dist_mpup_lf8down(l3mec_c,:))/sqrt(length(l3mec_c)),'r');
errorbar_tick(h,.01,'units');
h = errorbar(amp_grid,nanmean(lf3dist_mpdown_lf8up(l3mec_c,:)),nanstd(lf3dist_mpdown_lf8up(l3mec_c,:))/sqrt(length(l3mec_c)),'k');
errorbar_tick(h,.01,'units');
h = errorbar(amp_grid,nanmean(lf3dist_mpup_lf8up(l3mec_c,:)),nanstd(lf3dist_mpup_lf8up(l3mec_c,:))/sqrt(length(l3mec_c)),'c');
errorbar_tick(h,.01,'units');
xlim([-3 3])
title('L3MEC')
legend('MP UP LF8 Down','MP DOWN LF8 Down','MP DOWN LF8 UP','MP UP LF8 UP')

figure
h = errorbar(amp_grid,nanmean(lf3dist_mpdown_lf8down(l3lec_c,:)),nanstd(lf3dist_mpdown_lf8down(l3lec_c,:))/sqrt(length(l3lec_c)));
errorbar_tick(h,.01,'units');
hold on
h = errorbar(amp_grid,nanmean(lf3dist_mpup_lf8down(l3lec_c,:)),nanstd(lf3dist_mpup_lf8down(l3lec_c,:))/sqrt(length(l3lec_c)),'r');
errorbar_tick(h,.01,'units');
h = errorbar(amp_grid,nanmean(lf3dist_mpdown_lf8up(l3lec_c,:)),nanstd(lf3dist_mpdown_lf8up(l3lec_c,:))/sqrt(length(l3lec_c)),'k');
errorbar_tick(h,.01,'units');
h = errorbar(amp_grid,nanmean(lf3dist_mpup_lf8up(l3lec_c,:)),nanstd(lf3dist_mpup_lf8up(l3lec_c,:))/sqrt(length(l3lec_c)),'c');
errorbar_tick(h,.01,'units');
xlim([-3 3])
title('L3LEC')
legend('MP UP LF8 Down','MP DOWN LF8 Down','MP DOWN LF8 UP','MP UP LF8 UP')

