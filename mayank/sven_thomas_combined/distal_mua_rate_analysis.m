clear all
close all

%%
cd C:\WC_Germany\sven_thomas_combined\
% load ./combined_dir.mat

load .\distal_dir
distal_dir = distal_dir(distal_usable);
mctx_lfp = ctx_lfp(distal_usable);
mhpc_lfp = hpc_lfp(distal_usable);
mhpc_mua = hpc_mua(distal_usable);

load ./distal_lec_dir.mat
distal_dir = [distal_dir distal_lec(usable_distal_lec)];
ctx_lfp = [mctx_lfp ctx_lfp(usable_distal_lec)];
hpc_lfp = [mhpc_lfp hpc_lfp(usable_distal_lec)];
hpc_mua = [mhpc_mua hpc_mua(usable_distal_lec)];
distal_mec = (1:length(distal_usable));
distal_lec = (length(distal_usable)+1):(length(distal_usable)+length(usable_distal_lec));

raw_Fs = 2016;
int_width_range = [2 8]; %range of spike widths (in 32kHz samples) for interneuron spikes
pyr_width_range = [10 18]; %range of spike widths for pyramidal spikes
min_corr = 0.4; %minimum correlation coefficient with the corresponding spike template

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
lcf = 0.05;
hcf = 4;
maxlag = round(Fsd*1);
%%
for d = 1:length(distal_dir)
    cd(distal_dir{d})
    pwd
    load ./used_data lf8 lf7 lf6 lf5
    if ctx_lfp(d) == 7
        lf8 = lf7;
    elseif ctx_lfp(d) == 6
        lf8 = lf6;
    elseif ctx_lfp(d) == 5
        lf8 = lf5;
    end
    sess_dur(d) = length(lf8)/raw_Fs;
    
    if exist('./mua_data3.mat','file')
        load ./mua_data3
        
        lf8_lf = get_lf_features(lf8,raw_Fs,Fsd,[lcf hcf]);
        load ./sync_times.mat
        synct_d = downsample(synct,dsf);
        
        load ./pa_hsmm_state_seq_combined_fin_nd.mat
        load ./pa_hsmm_state_seq7_combined_fin_nd.mat
        [new_seg_inds] = resample_uds_seg_inds_v2(hsmm7.UDS_segs,hsmm7.Fs,Fsd,hsmm_bbstate_seq7);
        dur_uds = sum(diff(new_seg_inds,[],2))/Fsd;
        [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq);
        [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,hsmm_bbstate_seq7);
        
        is_sync = zeros(size(lf8_lf));
        for j = 1:size(new_seg_inds,1)
            is_sync(new_seg_inds(j,1):new_seg_inds(j,2)) = 1;
        end
        is_sync = logical(is_sync);
        for j = 1:8
            mua_rate(d,j) = length(mua_times{j})/sess_dur(d);
            mua_avg_widths(d,j) = mean(mua_widths{j});
            mua_avg_amps(d,j) = mean(mua_amps{j});
            
            mua_binned = hist(mua_times{j},synct_d);
            mua_binned([1 end]) = 0;
            cur_mua_rate = jmm_smooth_1d_cor(mua_binned,round(Fsd*0.05));
            cur_mua_rate = (cur_mua_rate - mean(cur_mua_rate(is_sync)))/std(cur_mua_rate(is_sync));
            [mua_lfp_xc(d,j,:),lags] = xcov(lf8_lf,cur_mua_rate,maxlag,'coeff');
            
            mua_u8_trig_avg(d,j,:) = get_event_trig_avg(cur_mua_rate,up_trans_inds8,maxlag,maxlag);
            mua_uw_trig_avg(d,j,:) = get_event_trig_avg(cur_mua_rate,up_trans_inds,maxlag,maxlag);
        end
        mua_ov_avg(d,:,:) = avg_waveform;
    else
        mua_rate(d,:) = nan(1,8);
    end
end

cd C:\WC_Germany\sven_thomas_combined\
save distal_mua_rate_data_nd

%%
used_recs = find(~isnan(hpc_mua));
figure
imagesc(lags/Fsd,1:7,squeeze(nanmean(mua_lfp_xc(used_recs,2:end,:))));
set(gca,'ydir','normal')
xlabel('Lag (s)','fontsize',14)
ylabel('Channel Number','fontsize',14)

figure
imagesc(lags/Fsd,1:7,squeeze(nanmean(mua_u8_trig_avg(used_recs,2:end,:))));
set(gca,'ydir','normal')
xlabel('Lag (s)','fontsize',14)
ylabel('Channel Number','fontsize',14)
figure
imagesc(lags/Fsd,1:7,squeeze(nanmean(mua_uw_trig_avg(used_recs,2:end,:))));
set(gca,'ydir','normal')
xlabel('Lag (s)','fontsize',14)
ylabel('Channel Number','fontsize',14)

figure
plot(1:7,mean(mua_rate(used_recs,2:end)),'o-')
% used_recs = find(isnan(hpc_mua));
% figure
% imagesc(lags/Fsd,1:7,squeeze(nanmean(mua_lfp_xc(used_recs,:,:))));
% set(gca,'ydir','normal')
xlabel('Chanel Number','fontsize',14)
ylabel('Average rate (Hz)','fontsize',14)
%%
used_recs = find(~isnan(hpc_mua));
for i = 1:length(used_recs)
    figure
    plot(1:7,mua_rate(used_recs(i),2:end),'o-')
    figure
    imagesc(lags/Fsd,1:7,squeeze(mua_lfp_xc(used_recs(i),2:end,:)))
    set(gca,'ydir','normal')
    hpc_mua(used_recs(i))
    pause
    close all
end
