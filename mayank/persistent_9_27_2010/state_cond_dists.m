cd G:\WC_Germany\persistent_2010

clear all
close all

load G:\WC_Germany\overall_EC\overall_EC_dir.mat
drive_letter = 'G';
used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

Fs = 2016;
dsf = 8;
Fsd = Fs/dsf;
niqf = Fs/2;
lcf = 0.05/niqf;
hcf = 4/niqf;
[b,a] = butter(2,[lcf hcf]);
Fsd = Fs/dsf;

amp_range = linspace(-4, 4, 400);
da = amp_range(2)-amp_range(1);
for d = 1:length(sess_data)

    cdir = sess_data(d).directory;
    cdir(1) = 'G';
    disp(sprintf('session %d',d))
    cd(cdir);
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);

    load ./used_data wcv_minus_spike lf8 lf3
    load ./pa_hsmm_state_seq8_new2
    lfp_state_seq = hsmm_bbstate_seq8;
    load ./pa_hsmm_state_seq_new2
    mp_state_seq = hsmm_bbstate_seq;
    
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);
    lf3_f = filtfilt(b,a,lf3);
    wcv_f = zscore(downsample(wcv_f,dsf));
    lf8_f = zscore(downsample(lf8_f,dsf));
    lf3_f = zscore(downsample(lf3_f,dsf));

    t_axis = (1:length(wcv_f))/Fsd;
    
    [new_seg_inds] = resample_uds_seg_inds_v2(hmm.UDS_segs,hmm.Fs,Fsd,mp_state_seq);
    
    mp_state_vec = nan(size(wcv_f));
    lf8_state_vec = nan(size(lf8_f));
    for ns = 1:hmm.Nsegs
        mp_state_vec(new_seg_inds(ns,1):new_seg_inds(ns,2)) = mp_state_seq{ns};
        lf8_state_vec(new_seg_inds(ns,1):new_seg_inds(ns,2)) = lfp_state_seq{ns};
    end
    mp_state_vec = mp_state_vec-1;
    lf8_state_vec = lf8_state_vec-1;
    
    mp_up_dm(d,:) = ksdensity(wcv_f(mp_state_vec==1),amp_range);
    mp_down_dm(d,:) = ksdensity(wcv_f(mp_state_vec==0),amp_range);
    lf8_up_d8(d,:) = ksdensity(lf8_f(lf8_state_vec==1),amp_range);
    lf8_down_d8(d,:) = ksdensity(lf8_f(lf8_state_vec==0),amp_range);
    mp_up_d3(d,:) = ksdensity(lf3_f(mp_state_vec==1),amp_range);
    mp_down_d3(d,:) = ksdensity(lf3_f(mp_state_vec==0),amp_range);
    lf8_up_d3(d,:) = ksdensity(lf3_f(lf8_state_vec==1),amp_range);
    lf8_down_d3(d,:) = ksdensity(lf3_f(lf8_state_vec==0),amp_range);
    mp_up_lf8_up_d3(d,:) = ksdensity(lf3_f(mp_state_vec==1 & lf8_state_vec==1),amp_range);
    mp_up_lf8_down_d3(d,:) = ksdensity(lf3_f(mp_state_vec==1 & lf8_state_vec==0),amp_range);
    mp_down_lf8_up_d3(d,:) = ksdensity(lf3_f(mp_state_vec==0 & lf8_state_vec==1),amp_range);
    mp_down_lf8_down_d3(d,:) = ksdensity(lf3_f(mp_state_vec==0 & lf8_state_vec==0),amp_range);
    d3(d,:) = ksdensity(lf3_f,amp_range);
    p_lf8(d) = sum(lf8_state_vec==1)/sum(~isnan(lf8_state_vec));
    p_mp(d) = sum(mp_state_vec==1)/sum(~isnan(mp_state_vec));
    p_mp_cond_lf8up(d) = sum(mp_state_vec==1 & lf8_state_vec == 1)/sum(lf8_state_vec==1);
    p_mp_cond_lf8down(d) = sum(mp_state_vec==1 & lf8_state_vec==0)/sum(lf8_state_vec==0);
    p_lf8_cond_mpup(d) = sum(lf8_state_vec==1 & mp_state_vec == 1)/sum(mp_state_vec==1);
    p_lf8_cond_mpdown(d) = sum(lf8_state_vec==1 & mp_state_vec==0)/sum(mp_state_vec==0);
    
    A_d8_dm = (1-p_lf8(d))*mp_down_lf8_down_d3(d,:)*(1-p_mp_cond_lf8down(d)).*log2(mp_down_lf8_down_d3(d,:)./...
        (lf8_down_d3(d,:)));
    A_u8_dm = p_lf8(d)*mp_down_lf8_up_d3(d,:)*(1-p_mp_cond_lf8up(d)).*log2(mp_down_lf8_up_d3(d,:)./...
        (lf8_up_d3(d,:)));
    A_d8_um = (1-p_lf8(d))*mp_up_lf8_down_d3(d,:)*p_mp_cond_lf8down(d).*log2(mp_up_lf8_down_d3(d,:)./...
        (lf8_down_d3(d,:)));
    A_u8_um = p_lf8(d)*mp_up_lf8_up_d3(d,:)*p_mp_cond_lf8up(d).*log2(mp_up_lf8_up_d3(d,:)./...
        (lf8_up_d3(d,:)));
    A_d8_dm(isnan(A_d8_dm) | isinf(A_d8_dm)) = 0;
    A_u8_dm(isnan(A_u8_dm) | isinf(A_u8_dm)) = 0;
    A_d8_um(isnan(A_d8_um) | isinf(A_d8_um)) = 0;
    A_u8_um(isnan(A_u8_um) | isinf(A_u8_um)) = 0;
    
    cond_mp_info(d) = da*(trapz(A_d8_dm)+trapz(A_u8_dm)+trapz(A_d8_um)+trapz(A_u8_um));

    A_d8_dm = (1-p_mp(d))*mp_down_lf8_down_d3(d,:)*(1-p_lf8_cond_mpdown(d)).*log2(mp_down_lf8_down_d3(d,:)./...
        (mp_down_d3(d,:)));
    A_u8_dm = (1-p_mp(d))*mp_down_lf8_up_d3(d,:)*p_lf8_cond_mpdown(d).*log2(mp_down_lf8_up_d3(d,:)./...
        (mp_down_d3(d,:)));
    A_d8_um = p_mp(d)*mp_up_lf8_down_d3(d,:)*(1-p_lf8_cond_mpup(d)).*log2(mp_up_lf8_down_d3(d,:)./...
        (mp_up_d3(d,:)));
    A_u8_um = p_mp(d)*mp_up_lf8_up_d3(d,:)*p_lf8_cond_mpup(d).*log2(mp_up_lf8_up_d3(d,:)./...
        (mp_up_d3(d,:)));
    A_d8_dm(isnan(A_d8_dm) | isinf(A_d8_dm)) = 0;
    A_u8_dm(isnan(A_u8_dm) | isinf(A_u8_dm)) = 0;
    A_d8_um(isnan(A_d8_um) | isinf(A_d8_um)) = 0;
    A_u8_um(isnan(A_u8_um) | isinf(A_u8_um)) = 0;
 
    cond_lf8_info(d) = da*(trapz(A_d8_dm)+trapz(A_u8_dm)+trapz(A_d8_um)+trapz(A_u8_um));

    %%
    A1 = mp_up_d3(d,:)*p_mp(d).*log2((mp_up_d3(d,:)*p_mp(d))./(p_mp(d)*d3(d,:)));
    A2 = mp_down_d3(d,:)*(1-p_mp(d)).*log2((mp_down_d3(d,:)*(1-p_mp(d)))./((1-p_mp(d))*d3(d,:)));
    A1(isnan(A1) | isinf(A1)) = 0;
    A2(isnan(A2) | isinf(A2)) = 0;
    mp_info(d) = da*(trapz(A1)+trapz(A2));

    A1 = lf8_up_d3(d,:)*p_lf8(d).*log2((lf8_up_d3(d,:)*p_lf8(d))./(p_lf8(d)*d3(d,:)));
    A2 = lf8_down_d3(d,:)*(1-p_lf8(d)).*log2((lf8_down_d3(d,:)*(1-p_lf8(d)))./((1-p_lf8(d))*d3(d,:)));
    A1(isnan(A1) | isinf(A1)) = 0;
    A2(isnan(A2) | isinf(A2)) = 0;
    lf8_info(d) = da*(trapz(A1)+trapz(A2));

%     plot(amp_range,d3(d,:),'k'), hold on
%     plot(amp_range,lf8_up_d3(d,:),'b')
%     plot(amp_range,lf8_down_d3(d,:),'r')
%     plot(amp_range,mp_up_lf8_up_d3(d,:),'g--','linewidth',2), hold on
%     plot(amp_range,mp_down_lf8_up_d3(d,:),'c--','linewidth',2)
%     plot(amp_range,mp_up_lf8_down_d3(d,:),'g','linewidth',2)
%     plot(amp_range,mp_down_lf8_down_d3(d,:),'c','linewidth',2)
%     fprintf('mp info: %.4f\nlf8 info: %.4f\n\n',cond_mp_info(d),cond_lf8_info(d))
%     pause
%     clf
    
%    clear mp* lfp*

end

%%
syn = cond_lf8_info - lf8_info;
tot = mp_info + cond_lf8_info;
rsyn = syn./tot;
rcmp = cond_mp_info./tot;
rclfp = cond_lf8_info./tot;
rmp = mp_info./tot;
rlfp = lf8_info./tot;

mec = 1:22;
lec = 23:36;
bad_lf3 = [7 9 11 12 13 16 17 25 31 32 35];
good_lf3 = setdiff(1:36,bad_lf3);
good_lf3_mec = good_lf3(ismember(good_lf3,mec));
good_lf3_lec = good_lf3(ismember(good_lf3,lec));

