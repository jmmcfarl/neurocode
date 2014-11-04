%%
clear all
load C:/WC_Germany/final_pdown_analysis/compiled_data.mat

min_rec_dur = 500; %in sec
data_ids = [data(:).id];
used_dirs = find([data(:).ep] > min_rec_dur & [data(:).dp] > min_rec_dur);

% used_dirs(ismember(data_ids(used_dirs),no_cell) | ismember(data_ids(used_dirs),unclear_uds)) = [];
used_dirs(ismember(data_ids(used_dirs),unclear_uds)) = [];

data = data(used_dirs);

load C:/WC_Germany/final_pdown_analysis/fin_pdown_core_analysis.mat
if length(core_data) ~= length(data)
    error('Data mismatch');
end
cd C:\WC_Germany\final_pdown_analysis\
% load('final_cortical_state_data.mat');
load('final_cortical_state_data_fin_nobuff.mat');
%% parameters
raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 10;
lcf = 0.05;
lcf_hf = 40;
hcf_hf = 100;
hcf_sm = 0.025;
rate_sm = round(Fsd*0.05);

state_bin_size = 40;
maxlag = round(400/state_bin_size);

%%
for d = 1:length(used_dirs)
    cd(data(d).dir)
    pwd
    
    load ./used_data lf7 wcv_minus_spike
    [lfp_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf]);
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
    lfp_hf = get_lf_features(lf7,raw_Fs,Fsd,[lcf_hf hcf_hf]);

    load ./sync_times.mat
    synct_d = downsample(synct,dsf);
    
    end_time = min(data(d).ep,data(d).dp);
    ep = find(t_axis >= end_time,1);
    if ~isempty(ep)
        synct_d(ep+1:end) = []; lfp_lf(ep+1:end) = []; wcv_lf(ep+1:end) = []; t_axis(ep+1:end) = [];
        lfp_hf(ep+1:end) = [];
    else
        ep = length(t_axis);
    end
    
    %%
    load ./pa_hsmm_state_seq_combined_fin_nd.mat
    load ./pa_hsmm_state_seq7_combined_fin_nd.mat
    hsmm8 = hsmm7;
    lfp_state_seq = hsmm_bbstate_seq7;
    mp_state_seq = hsmm_bbstate_seq;
    
    [new_seg_inds] = resample_uds_seg_inds_v2(hsmm.UDS_segs,hsmm.Fs,Fsd,mp_state_seq);
    dur_uds = sum(diff(new_seg_inds,[],2))/Fsd;
    [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
    [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);
    
    bad_mp_states = find(up_trans_inds > ep | down_trans_inds > ep);
    up_trans_inds(bad_mp_states) = []; down_trans_inds(bad_mp_states) = [];
    bad_lfp_states = find(up_trans_inds8 > ep | down_trans_inds8 > ep);
    up_trans_inds8(bad_lfp_states) = []; down_trans_inds8(bad_lfp_states) = [];
    
    %%
    load ./allEC_ctx_period_data_hsmm.mat
    lfp_period_vec = nan(size(wcv_lf));
    for i = 1:size(new_seg_inds,1)
        cur_inds = new_seg_inds(i,1):new_seg_inds(i,2);
        cur_inds_used = find(cur_inds <= ep);
        lfp_period_vec(cur_inds(cur_inds_used)) = lf8_period_f{i}(cur_inds_used);
    end
    
    uds_dur = sum(~isnan(lfp_period_vec))/Fsd;
    
    %%
    mp_state_number = nan(length(wcv_lf),1);
    mp_state_vec = zeros(length(wcv_lf),1);
    lfp_state_number = nan(length(wcv_lf),1);
    lfp_state_vec = zeros(length(wcv_lf),1);
    for i = 1:length(up_trans_inds)-1
        mp_state_vec(up_trans_inds(i):down_trans_inds(i)) = 1;
        mp_state_number(up_trans_inds(i):down_trans_inds(i)) = 2*(i-1)+1;
        mp_state_number(down_trans_inds(i):up_trans_inds(i+1)) = 2*(i-1) + 2;
    end
    for i = 1:length(up_trans_inds8)-1
        lfp_state_vec(up_trans_inds8(i):down_trans_inds8(i)) = 1;
        lfp_state_number(up_trans_inds8(i):down_trans_inds8(i)) = 2*(i-1)+1;
        lfp_state_number(down_trans_inds8(i):up_trans_inds8(i+1)) = 2*(i-1) + 2;
    end
    mp_state_number(isnan(lfp_period_vec)) = nan;
    lfp_state_number(isnan(lfp_period_vec)) = nan;
    mp_state_vec(isnan(lfp_period_vec)) = nan;
    lfp_state_number(isnan(lfp_period_vec)) = nan;
    
    %% compute corresponding state transitions and transition lags
    [corresp_lf8_upinds,corresp_lf8_downinds] = find_corresponding_state_transitions_lookback(...
        up_trans_inds,down_trans_inds,up_trans_inds8,down_trans_inds8);
    [corresp_mp_upinds,corresp_mp_downinds] = find_corresponding_state_transitions_lookback(...
        up_trans_inds8,down_trans_inds8,up_trans_inds,down_trans_inds);
    
    %%
    state_edges = [0:state_bin_size:t_axis(end)];
    state_axis = 0.5*state_edges(1:end-1)+0.5*state_edges(2:end);
       
    binned_ctx_lfpow = zeros(length(state_axis),1);
    binned_ctx_hfpow = zeros(length(state_axis),1);
    binned_ctx_up_amps = zeros(length(state_axis),1);
    binned_ctx_up_ampshf = zeros(length(state_axis),1);
    binned_ctx_up_durs = zeros(length(state_axis),1);
    binned_ctx_down_durs = zeros(length(state_axis),1);
    binned_ctx_dc = zeros(length(state_axis),1);
    for ii = 1:length(state_axis)
        cur_set = find(up_trans_inds8/Fsd >= state_edges(ii) & up_trans_inds8/Fsd <= state_edges(ii+1));
        cur_inds = find(t_axis >= state_edges(ii) & t_axis <= state_edges(ii+1));
        binned_ctx_lfpow(ii) = std(lfp_lf(cur_inds));
        binned_ctx_hfpow(ii) = std(lfp_hf(cur_inds));
        if ~isempty(cur_set)
            binned_ctx_up_amps(ii) = nanmean(cfun_data(d).ctx_up_amps(cur_set));
            binned_ctx_up_ampshf(ii) = nanmean(cfun_data(d).ctx_up_amps_hf(cur_set));
            binned_ctx_up_durs(ii) = nanmean(core_data(d).lfp_up_durs(cur_set));
            binned_ctx_down_durs(ii) = nanmean(core_data(d).lfp_down_durs(cur_set));
            binned_ctx_dc(ii) = nanmean(core_data(d).lfp_dutycycs(cur_set));
        elseif ii > 1
            binned_ctx_up_amps(ii) = binned_ctx_up_amps(ii-1);
            binned_ctx_up_ampshf(ii) = binned_ctx_up_ampshf(ii-1);
            binned_ctx_up_durs(ii) = binned_ctx_up_durs(ii-1);
            binned_ctx_down_durs(ii) = binned_ctx_down_durs(ii-1);
            binned_ctx_dc(ii) = binned_ctx_dc(ii-1);
        end
    end
    
    %%
    skipped_lfp_ups = core_data(d).skipped_lfp_ups;
    skipped_lfp_downs = core_data(d).skipped_lfp_downs;
    
    rt2_down_times = t_axis(down_trans_inds(core_data(d).rt2_downs));
     rt2_up_times = t_axis(up_trans_inds(core_data(d).rt2_ups));
%     rt2_down_times = t_axis(up_trans_inds8(skipped_lfp_ups));
%      rt2_up_times = t_axis(down_trans_inds8(skipped_lfp_downs));

%      all_up_times = t_axis(up_trans_inds8);
%    all_down_times = t_axis(down_trans_inds8);
     all_up_times = t_axis(up_trans_inds(core_data(d).nrt2_ups));
   all_down_times = t_axis(down_trans_inds(core_data(d).nrt2_downs));
    
    binned_rt2_ups = histc(rt2_up_times,state_edges);
    binned_rt2_downs = histc(rt2_down_times,state_edges);
    binned_rt2_ups(end) = [];
    binned_rt2_downs(end) = [];
    binned_tot_ups = histc(all_up_times,state_edges);
    binned_tot_downs = histc(all_down_times,state_edges);
    binned_tot_ups(end) = [];
    binned_tot_downs(end) = [];
    binned_pup_frac = binned_rt2_ups(:)./binned_tot_ups(:);
    binned_pdown_frac = binned_rt2_downs(:)./binned_tot_downs(:);
    
    uset = find(~isnan(binned_pup_frac) & ~isnan(binned_pdown_frac) & ~isinf(binned_pup_frac) & ~isinf(binned_pdown_frac));
%     binned_pup_frac(cut) = [];
%     binned_pdown_frac(cut) = [];
%     binned_ctx_dc(cut) = [];
%     binned_ctx_down_durs(cut) = [];
    %%
    [xcorr_data(d).pup_acorr,lags] = xcov(binned_pup_frac(uset),maxlag,'coeff');
    [xcorr_data(d).pdown_acorr,lags] = xcov(binned_pdown_frac(uset),maxlag,'coeff');
    [xcorr_data(d).pup_pdown_xcorr,lags] = xcov(binned_pup_frac(uset),binned_pdown_frac(uset),maxlag,'coeff');
    
    xcorr_data(d).sp_corr_pdown_pup = corr(binned_pup_frac(uset),binned_pdown_frac(uset),'type','spearman');
    
    xcorr_data(d).sp_corr_pdown_dc = corr(binned_pdown_frac(uset),binned_ctx_dc(uset),'type','spearman');
    xcorr_data(d).sp_corr_pup_dc = corr(binned_pup_frac(uset),binned_ctx_dc(uset),'type','spearman');
    xcorr_data(d).sp_corr_pdown_lfpow = corr(binned_pdown_frac(uset),binned_ctx_lfpow(uset),'type','spearman');
    xcorr_data(d).sp_corr_pup_lfpow = corr(binned_pup_frac(uset),binned_ctx_lfpow(uset),'type','spearman');
    xcorr_data(d).sp_corr_pdown_downdur = corr(binned_pdown_frac(uset),binned_ctx_down_durs(uset),'type','spearman');
    xcorr_data(d).sp_corr_pup_downdur = corr(binned_pup_frac(uset),binned_ctx_down_durs(uset),'type','spearman');
    
    [xcorr_data(d).ctx_upamp_acorr] = xcov(binned_ctx_up_amps,maxlag,'coeff');
    [xcorr_data(d).ctx_upamphf_acorr] = xcov(binned_ctx_up_ampshf,maxlag,'coeff');
    [xcorr_data(d).ctx_updur_acorr] = xcov(binned_ctx_up_durs,maxlag,'coeff');
    [xcorr_data(d).ctx_downdur_acorr] = xcov(binned_ctx_down_durs,maxlag,'coeff');
    [xcorr_data(d).ctx_dc_acorr] = xcov(binned_ctx_dc,maxlag,'coeff');
    [xcorr_data(d).ctx_lfpow] = xcov(binned_ctx_lfpow,maxlag,'coeff');
    [xcorr_data(d).ctx_hfpow] = xcov(binned_ctx_hfpow,maxlag,'coeff');
   
    [xcorr_data(d).pdown_dc_xcorr] = xcov(binned_ctx_dc,binned_pdown_frac,maxlag,'coeff');
    [xcorr_data(d).pdown_ctx_downdur_xcorr] = xcov(binned_ctx_down_durs,binned_pdown_frac,maxlag,'coeff');
    
    xcorr_data(d).uds_dur = uds_dur;
end

%%
cd C:\WC_Germany\final_pdown_analysis\
save final_xcorr_data4 lags Fsd xcorr_data
%%
data_ids = [data(:).id];
l3mec = find(strcmp({data.loc},'MEC'));
l3lec = find(strcmp({data.loc},'LEC'));
l3mec(~ismember(data_ids(l3mec),clear_l3)) = [];

n_rt2_downs = cellfun(@(x) length(x),{core_data(:).rt2_downs});
n_rt2_ups = cellfun(@(x) length(x),{core_data(:).rt2_ups});

min_npers = 5;
l3mec_pups = l3mec(n_rt2_ups(l3mec) >= min_npers);
l3mec_pdowns = l3mec(n_rt2_downs(l3mec) >= min_npers);
l3mec_both = intersect(l3mec_pups,l3mec_pdowns);

%%
sp_pup_pdown_xcorr = [xcorr_data(:).sp_corr_pdown_pup];
sp_pup_dc_xcorr = [xcorr_data(:).sp_corr_pup_dc];
sp_pdown_dc_xcorr = [xcorr_data(:).sp_corr_pdown_dc];
%%
pup_acorr = [xcorr_data(:).pup_acorr]';
pdown_acorr = [xcorr_data(:).pdown_acorr]';
pupdown_xcorr = [xcorr_data(:).pup_pdown_xcorr]';
zpt = find(lags == 0);
pup_acorr(:,zpt) = nan;
pdown_acorr(:,zpt) = nan;

zlag_xcorrs = pupdown_xcorr(:,zpt);

%%
pdown_dc_xcorr = [xcorr_data(:).pdown_dc_xcorr]';
pdown_downdur_xcorr = [xcorr_data(:).pdown_ctx_downdur_xcorr]';
zpt = find(lags == 0);


%%
% h1 = figure; hold on
% errorbar(lags*state_bin_size,nanmean(pup_acorr(l3mec_pups,:)),nanstd(pup_acorr(l3mec_pups,:))/sqrt(length(l3mec_pups)));
% errorbar(lags*state_bin_size,nanmean(pdown_acorr(l3mec_pdowns,:)),nanstd(pdown_acorr(l3mec_pdowns,:))/sqrt(length(l3mec_pdowns)),'r');

h2 = figure; hold on;
% errorbar(lags*state_bin_size,nanmean(pupdown_xcorr(l3mec_both,:)),nanstd(pupdown_xcorr(l3mec_both,:))/sqrt(length(l3mec_both)));
shadedErrorBar(lags*state_bin_size,nanmean(pupdown_xcorr(l3mec_both,:)),nanstd(pupdown_xcorr(l3mec_both,:))/sqrt(length(l3mec_both)));

