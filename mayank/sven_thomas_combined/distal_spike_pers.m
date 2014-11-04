clear all
cd C:\WC_Germany\sven_thomas_combined\
load ./distal_core_analysis_fin_nd_np.mat

load ./distal_dir
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

dsf = 8;
raw_Fs = 2016;
Fsd = raw_Fs/dsf;
rate_sm = round(Fsd*0.05);

%%
for d = 1:length(distal_dir)
    
    cdir = distal_dir{d};
    disp(sprintf('session %d',d))
    cd(cdir);
    
    load ./spike_time_jmm
    load ./used_data wcv
    load ./sync_times.mat
    rec_dur(d) = length(wcv)/raw_Fs;
    synct_d = downsample(synct,dsf);
    if ~isnan(hpc_mua(d))
        load ./mua_data3
        hpc_mua_times = mua_times{hpc_mua(d)};
        hpc_mua_times(hpc_mua_times > synct_d(end) | hpc_mua_times < synct_d(1)) = [];
        hpc_mua_rate =hist(hpc_mua_times,synct_d)*Fsd;
        hpc_mua_rate = jmm_smooth_1d_cor(hpc_mua_rate,rate_sm);
        hpc_mua_rate = zscore(hpc_mua_rate');
        
        if ctx_lfp(d) == 8
            ctx_mua_times = mua_times{8};
        else
            ctx_mua_times = mua_times{7};
        end
        ctx_mua_rate =hist(ctx_mua_times,synct_d)*Fsd;
        if length(hpc_mua_rate) > length(synct_d)
            hpc_mua_rate = hpc_mua_rate(1:length(synct_d));
            ctx_mua_rate = ctx_mua_rate(1:length(synct_d));
        end
    else
        hpc_mua_rate = nan(size(synct_d));
        ctx_mua_rate = nan(size(synct_d));
    end
    
    mp_spike_rate = hist(synct(spkid),synct_d)*Fsd;
    
    backgnd_rate(d) = length(spkid)/length(synct_d)*Fsd;
    
%     load ./pa_hsmm_state_seq_combined_fin.mat
%     load ./pa_hsmm_state_seq7_combined_fin.mat
    load ./pa_hsmm_state_seq_combined_fin_nd.mat
    load ./pa_hsmm_state_seq7_combined_fin_nd.mat
    mp_state_seq = hsmm_bbstate_seq;
    lfp_state_seq = hsmm_bbstate_seq7;
    
    [new_seg_inds] = resample_uds_seg_inds_v2(hmm.UDS_segs,hmm.Fs,Fsd,mp_state_seq);
    seg_durs = new_seg_inds(:,2)-new_seg_inds(:,1)+1;
    
    [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
    up_state_durs = (down_trans_inds-up_trans_inds)/Fsd;
    
    wcv_up_log = nan(size(synct_d));
    lf8_up_log = nan(size(synct_d));
    
    for ns = 1:hmm.Nsegs
        cur_seg = new_seg_inds(ns,1):new_seg_inds(ns,2);
        wcv_up_log(cur_seg) = logical(mp_state_seq{ns}-1);
        lf8_up_log(cur_seg) = logical(lfp_state_seq{ns}-1);
    end
    
    mp_upstate_rate{d} = nan(size(up_trans_inds));
    for j = 1:length(up_trans_inds)
        mp_upstate_rate{d}(j) = mean(mp_spike_rate(up_trans_inds(j):down_trans_inds(j)));
    end
    
    temp_pers_mp_rate = nan(size(rt2_ups{d}));
    temp_npers_mp_rate = nan(size(nrt2_ups{d}));
    temp_pers_hpc_rate = nan(size(rt2_ups{d}));
    temp_npers_hpc_rate = nan(size(nrt2_ups{d}));
    temp_pers_during_mp = nan(size(rt2_ups{d}));
    temp_pers_around_mp = nan(size(rt2_ups{d}));
    temp_pers_during_hpc = nan(size(rt2_ups{d}));
    temp_pers_around_hpc = nan(size(rt2_ups{d}));
    for j = 1:length(rt2_ups{d})
        mp_up_inds = up_trans_inds(rt2_ups{d}(j)):down_trans_inds(rt2_ups{d}(j));
        temp_pers_mp_rate(j) = mean(mp_spike_rate(mp_up_inds));
        temp_pers_hpc_rate(j) = mean(hpc_mua_rate(mp_up_inds)); 
        
        during_inds = mp_up_inds(lf8_up_log(mp_up_inds) == 0);
        around_inds = mp_up_inds(lf8_up_log(mp_up_inds) == 1);
        if ~isempty(during_inds)
            temp_pers_during_mp(j) = mean(mp_spike_rate(during_inds));
            temp_pers_during_hpc(j) = mean(hpc_mua_rate(during_inds));
        end
        if ~isempty(around_inds)
            temp_pers_around_mp(j) = mean(mp_spike_rate(around_inds));
            temp_pers_around_hpc(j) = mean(hpc_mua_rate(around_inds));
        end
    end
    for j = 1:length(nrt2_ups{d})
        temp_npers_mp_rate(j) = mean(mp_spike_rate(up_trans_inds(nrt2_ups{d}(j)):down_trans_inds(nrt2_ups{d}(j))));
        temp_npers_hpc_rate(j) = mean(hpc_mua_rate(up_trans_inds(nrt2_ups{d}(j)):down_trans_inds(nrt2_ups{d}(j))));        
    end
    if ~isempty(rt2_ups{d})
        mp_pers_rate(d) = mean(temp_pers_mp_rate);
        hpc_pers_rate(d) = mean(temp_pers_hpc_rate);
    else
        mp_pers_rate(d) = nan;
        hpc_pers_rate(d) = nan;
    end
    hpc_npers_rate(d) = mean(temp_npers_hpc_rate);
     mp_npers_rate(d) = mean(temp_npers_mp_rate);
       
     mp_before_during_diff(d) = nanmean(temp_pers_during_mp-temp_pers_around_mp);
     hpc_before_during_diff(d) = nanmean(temp_pers_during_hpc-temp_pers_around_hpc);
     
    mp_up_rate(d) = mean(mp_spike_rate(wcv_up_log==1));
    mp_down_rate(d) = mean(mp_spike_rate(wcv_up_log==0));
    mp_up_rate_lfup(d) = mean(mp_spike_rate(wcv_up_log==1 & lf8_up_log==1));
    mp_up_rate_lfdown(d) = mean(mp_spike_rate(wcv_up_log==1 & lf8_up_log==0));
    mp_rate(d) = mean(mp_spike_rate(wcv_up_log==1 | wcv_up_log==0));
    
    if ~isnan(hpc_mua(d))
        ctx_rate_lfup(d) = mean(ctx_mua_rate(lf8_up_log==1));
        ctx_rate_lfdown(d) = mean(ctx_mua_rate(lf8_up_log==0));
        ctx_up_rate(d) = mean(ctx_mua_rate(wcv_up_log==1));
        ctx_up_rate_lfup(d) = mean(ctx_mua_rate(wcv_up_log==1 & lf8_up_log==1));
        ctx_up_rate_lfdown(d) = mean(ctx_mua_rate(wcv_up_log==1 & lf8_up_log==0));
        ctx_down_rate(d) = mean(ctx_mua_rate(wcv_up_log==0));
        ctx_down_rate_lfup(d) = mean(ctx_mua_rate(wcv_up_log==0 & lf8_up_log==1));
        ctx_down_rate_lfdown(d) = mean(ctx_mua_rate(wcv_up_log==0 & lf8_up_log==0));
        ctx_rate(d) = mean(ctx_mua_rate(wcv_up_log==1 | wcv_up_log==0));
        hpc_rate_lfup(d) = mean(hpc_mua_rate(lf8_up_log==1));
        hpc_rate_lfdown(d) = mean(hpc_mua_rate(lf8_up_log==0));
        hpc_up_rate(d) = mean(hpc_mua_rate(wcv_up_log==1));
        hpc_up_rate_lfup(d) = mean(hpc_mua_rate(wcv_up_log==1 & lf8_up_log==1));
        hpc_up_rate_lfdown(d) = mean(hpc_mua_rate(wcv_up_log==1 & lf8_up_log==0));
        hpc_down_rate(d) = mean(hpc_mua_rate(wcv_up_log==0));
        hpc_down_rate_lfup(d) = mean(hpc_mua_rate(wcv_up_log==0 & lf8_up_log==1));
        hpc_down_rate_lfdown(d) = mean(hpc_mua_rate(wcv_up_log==0 & lf8_up_log==0));
        hpc_rate(d) = mean(hpc_mua_rate(wcv_up_log==1 | wcv_up_log==0));
        
        [P,T,STATS] = anovan(hpc_mua_rate,{wcv_up_log(:) lf8_up_log(:)},'display','off');
        
    else
        ctx_rate_lfup(d) = nan;
        ctx_rate_lfdown(d) = nan;
       ctx_up_rate(d) = nan;
       ctx_up_rate_lfup(d) = nan;
       ctx_up_rate_lfdown(d) = nan;
       ctx_down_rate(d) = nan;
       ctx_down_rate_lfup(d) = nan;
       ctx_down_rate_lfdown(d) = nan;
       ctx_rate(d) = nan;
        hpc_rate_lfup(d) = nan;
        hpc_rate_lfdown(d) = nan;
       hpc_up_rate(d) = nan;
       hpc_up_rate_lfup(d) = nan;
       hpc_up_rate_lfdown(d) = nan;
       hpc_down_rate(d) = nan;
       hpc_down_rate_lfup(d) = nan;
       hpc_down_rate_lfdown(d) = nan;
       hpc_rate(d) = nan;
    end
       
end

mp_spike_mod_index = (mp_up_rate_lfup-mp_up_rate_lfdown)./(mp_up_rate_lfup+mp_up_rate_lfdown);
ctx_spike_mod_index = (ctx_up_rate_lfup-ctx_up_rate_lfdown)./(ctx_up_rate_lfup+ctx_up_rate_lfdown);
hpc_spike_mod_index = (hpc_up_rate_lfup-hpc_up_rate_lfdown)./(hpc_up_rate_lfup+hpc_up_rate_lfdown);

cd C:\WC_Germany\sven_thomas_combined\
save distal_spike_rate_data_fin_nd_np mp_* hpc_* ctx_* rec_dur backgnd_rate

%%
uset = find(~isnan(hpc_mua));
Y = [hpc_up_rate_lfup(uset)'; hpc_up_rate_lfdown(uset)'; hpc_down_rate_lfdown(uset)'; hpc_down_rate_lfup(uset)'];
G = [ones(length(uset),1); 2*ones(length(uset),1); 3*ones(length(uset),1); 4*ones(length(uset),1)];
boxplot(Y,G)
set(gca,'xtick',1:4, 'xticklabel',{'MecUP/CtxUp','MecUP/CtxDOWN','MecDOWN/CtxDOWN','MecDOWN/CtxUP'},'fontsize',12)
ylabel('Average rate (z)','fontsize',16)
set(gca,'fontsize',14,'fontname','arial')
ylim([-0.7936 1.1288])
