
clear all
data_dir = cell(0);
data_type = cell(0);
data_exp = cell(0);
data_ep = [];
data_dp = [];
data_hpc_lfp = [];

cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir_nd.mat
for ii = 1:length(l3mec)
    data_dir = cat(2,data_dir,combined_dir{l3mec(ii)});
    data_type = cat(2,data_type,{'L3MEC'});
    if ismember(l3mec(ii),old_data_inds)
        data_exp = cat(2,data_exp,{'T'});
    else
        data_exp = cat(2,data_exp,{'S'});
    end
    data_ep = cat(2,data_ep,Inf);
    data_dp = cat(2,data_dp,nan);
    data_hpc_lfp = cat(2,data_hpc_lfp,hpc_lfp(l3mec(ii)));
end
for ii = 1:length(l3mec_np)
    data_dir = cat(2,data_dir,combined_dir{l3mec_np(ii)});
    data_type = cat(2,data_type,{'L3MEC'});
    if ismember(l3mec_np(ii),old_data_inds)
        data_exp = cat(2,data_exp,{'T'});
    else
        data_exp = cat(2,data_exp,{'S'});
    end
    data_ep = cat(2,data_ep,Inf);
    data_dp = cat(2,data_dp,nan);
    data_hpc_lfp = cat(2,data_hpc_lfp,hpc_lfp(l3mec_np(ii)));
end
for ii = 1:length(l3lec)
    data_dir = cat(2,data_dir,combined_dir{l3lec(ii)});
    data_type = cat(2,data_type,{'L3LEC'});
    if ismember(l3lec(ii),old_data_inds)
        data_exp = cat(2,data_exp,{'T'});
    else
        data_exp = cat(2,data_exp,{'S'});
    end
    data_ep = cat(2,data_ep,Inf);
    data_dp = cat(2,data_dp,nan);
    data_hpc_lfp = cat(2,data_hpc_lfp,hpc_lfp(l3lec(ii)));
end
for ii = 1:length(l3lec_np)
    data_dir = cat(2,data_dir,combined_dir{l3lec_np(ii)});
    data_type = cat(2,data_type,{'L3LEC'});
    if ismember(l3lec_np(ii),old_data_inds)
        data_exp = cat(2,data_exp,{'T'});
    else
        data_exp = cat(2,data_exp,{'S'});
    end
    data_ep = cat(2,data_ep,Inf);
    data_dp = cat(2,data_dp,nan);
    data_hpc_lfp = cat(2,data_hpc_lfp,hpc_lfp(l3lec_np(ii)));
end

cd C:\WC_Germany\persistent_downs\
load ./new_pdown_dir
cur_uset = find(new_pdown_use == 1);
for ii = 1:length(cur_uset)
    data_dir = cat(2,data_dir,new_pdown_dir{cur_uset(ii)});
    data_type = cat(2,data_type,{'L3MEC'});
    data_exp = cat(2,data_exp,{'S'});
    data_ep = cat(2,data_ep,new_pdown_ep(cur_uset(ii)));
    data_dp = cat(2,data_ep,new_pdown_dp(cur_uset(ii)));
    data_hpc_lfp = cat(2,data_hpc_lfp,2);
end

addpath('C:\WC_Germany\parietal_cortical_2010\');
addpath('C:\WC_Germany\persistent_9_27_2010\');
addpath('C:\WC_Germany\new_mec\');
addpath('C:\WC_Germany\overall_EC');
addpath('C:\WC_Germany\hsmm_state_detection');
addpath('C:\WC_Germany\sven_thomas_combined');
addpath('C:\Code\general_functions\');
%%

raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 10;
lcf = 0.05;
lcf_hf = 15;
hcf_hf = 80;
hcf_sm = 0.025;
rate_sm = round(Fsd*0.05);

maxlag = round(4*Fsd);
backlag = 4*Fsd;
forwardlag = 4*Fsd;

params.Fs = raw_Fs;
params.fpass = [0 100];
params.tapers = [3 5];
win = 25;

%robust persistence requires MP skips states longer than this
thresh_lf8_updur = 0.5;
thresh_lf8_downdur = 0.5;

%for state duration distribution calculations
up_range = [0.1 15];
down_range = [0.1 15];
numBins = 60;
log_dur_grid = logspace(log10(up_range(1)),log10(up_range(2)),numBins+1);

min_rec_dur = 500; %in sec
used_dirs = find(data_ep > min_rec_dur);

state_bin_size = 40;
maxlag = round(400/state_bin_size);

    %%
for d = 1:length(used_dirs)
    cd(data_dir{used_dirs(d)})
    pwd
    
    load ./used_data lf7 wcv_minus_spike
    [lf8_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf]);
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
    lf8_hf = get_hf_features(lf7,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
    lf8_hfr = get_lf_features(lf7,raw_Fs,Fsd,[lcf_hf hcf_hf]);
    
        load ./sync_times.mat
    synct_d = downsample(synct,dsf);
    
    end_time = min(data_ep(used_dirs(d)),data_dp(used_dirs(d)));
    ep = find(t_axis >= end_time,1);
    if ~isempty(ep)
        synct_d(ep+1:end) = []; %ac_time(ep+1:end) = [];
        lf8_lf(ep+1:end) = []; lf8_hf(ep+1:end) = []; lf8_hfr(ep+1:end) = [];
        wcv_lf(ep+1:end) = []; t_axis(ep+1:end) = []; 
    else
        ep = length(t_axis);
    end
    if exist('./all_combined_mp_uds.mat','file')
        load ./all_combined_mp_uds.mat
        load ./all_combined_lf7_uds.mat
    else
        load ./pa_hsmm_state_seq7_combined_fin_nd.mat
        load ./pa_hsmm_state_seq_combined_fin_nd.mat
    end
    hmm8 = hmm7;
    lfp_state_seq = hmm_bbstate_seq7;
    mp_state_seq = hmm_bbstate_seq;
    
    %         if you want to flip the roles of the MP and LFP
    %         temp = hsmm_bbstate_seq8;
    %         hsmm_bbstate_seq8 = hsmm_bbstate_seq;
    %         hsmm_bbstate_seq = temp;
    
    [new_seg_inds] = resample_uds_seg_inds_v2(hmm.UDS_segs,hmm.Fs,Fsd,mp_state_seq);
    dur_uds = sum(diff(new_seg_inds,[],2))/Fsd;
    [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
    [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);
    
    bad_mp_states = find(up_trans_inds > ep | down_trans_inds > ep);
    up_trans_inds(bad_mp_states) = []; down_trans_inds(bad_mp_states) = [];
    bad_lfp_states = find(up_trans_inds8 > ep | down_trans_inds8 > ep);
    up_trans_inds8(bad_lfp_states) = []; down_trans_inds8(bad_lfp_states) = [];
    
    fract_uds_dur(d) = dur_uds/t_axis(end);
    
    %%
    if exist('./allEC_ctx_period_data.mat','file')
        load ./allEC_ctx_period_data.mat
    else
        load ./combined_lf7_period_data_fin_nd
    end
    lf8_period_vec = nan(size(wcv_lf));
    for i = 1:size(new_seg_inds,1)
        cur_inds = new_seg_inds(i,1):new_seg_inds(i,2);
        cur_inds_used = find(cur_inds <= ep);
        lf8_period_vec(cur_inds(cur_inds_used)) = lf8_period_f{i}(cur_inds_used);
    end
    
    %% compute state duration distributions
    [mp_state_durations{d}] = compute_state_durations_seg(mp_state_seq,Fsd);
    [lfp_state_durations{d}] = compute_state_durations_seg(lfp_state_seq,Fsd);
    mp_state_durations{d}{1}(bad_mp_states) = []; mp_state_durations{d}{2}(bad_mp_states) = [];
    lfp_state_durations{d}{1}(bad_lfp_states) = []; lfp_state_durations{d}{2}(bad_lfp_states) = [];

    mp_updurs = mp_state_durations{d}{2};
    mp_downdurs = mp_state_durations{d}{1};
    lf8_updurs= lfp_state_durations{d}{2};
    lf8_downdurs = lfp_state_durations{d}{1};
               mp_dutycycs = mp_updurs./(mp_updurs+mp_downdurs);
    lf8_dutycycs = lf8_updurs./(lf8_updurs+lf8_downdurs);
 
    %% compute corresponding state transitions and transition lags
    [corresp_lf8_upinds{d},corresp_lf8_downinds{d}] = greedy_find_corresponding_ncx_state_transitions_simp(...
        up_trans_inds,down_trans_inds,up_trans_inds8,down_trans_inds8);
    
    %find non-skipped states
    non_skipped_mp_up_states{d} = find(~isnan(corresp_lf8_upinds{d}));
    non_skipped_mp_down_states{d} = find(~isnan(corresp_lf8_downinds{d}));
    
    %compute transition lags for non-skipped states
    mp_uplags{d} = nan(size(up_trans_inds));
    mp_downlags{d} = nan(size(up_trans_inds));
    mp_uplags{d}(non_skipped_mp_up_states{d}) = up_trans_inds(non_skipped_mp_up_states{d}) - up_trans_inds8(corresp_lf8_upinds{d}(non_skipped_mp_up_states{d}));
    mp_downlags{d}(non_skipped_mp_down_states{d}) = down_trans_inds(non_skipped_mp_down_states{d}) - down_trans_inds8(corresp_lf8_downinds{d}(non_skipped_mp_down_states{d}));
    
    %compute transition lags relative to ctx state durations
    mp_reldownlags{d} = nan(size(up_trans_inds));
    mp_reluplags{d} = nan(size(up_trans_inds));
    for i = 1:length(up_trans_inds)
        if ~isnan(corresp_lf8_downinds{d}(i))
            mp_reldownlags{d}(i) = mp_downlags{d}(i)/lf8_downdurs(corresp_lf8_downinds{d}(i))/Fsd;
        end
        if ~isnan(corresp_lf8_upinds{d}(i))
            mp_reluplags{d}(i) = mp_uplags{d}(i)/lf8_updurs(corresp_lf8_upinds{d}(i))/Fsd;
        end
    end
    
    %% compute persistence
    [mp_upskipped{d},mp_downskipped{d}] = greedy_find_skipped_ncx_states(...
        corresp_lf8_upinds{d},corresp_lf8_downinds{d},lfp_state_durations{d}{2},lfp_state_durations{d}{1},thresh_lf8_downdur,thresh_lf8_updur);
    
    rt2_ups{d} = find(mp_upskipped{d}.rnum_skipped > 0);
    t2_ups{d} = find(mp_upskipped{d}.num_skipped > 0);
    nrt2_ups{d} = find(mp_upskipped{d}.num_skipped == 0);
    fract_rt2_ups(d) = length(rt2_ups{d})/(length(rt2_ups{d}) + length(nrt2_ups{d}));
    fract_t2_ups(d) = length(t2_ups{d})/(length(t2_ups{d}) + length(nrt2_ups{d}));
    
    rt2_downs{d} = find(mp_downskipped{d}.rnum_skipped > 0);
    nrt2_downs{d} = find(mp_downskipped{d}.num_skipped==0); %number of down states that didn't skip any LFP up states
    fract_rt2_downs(d) = length(rt2_downs{d})/(length(rt2_downs{d}) + length(nrt2_downs{d}));

    state_edges = [0:state_bin_size:t_axis(end)];
    state_axis = 0.5*state_edges(1:end-1)+0.5*state_edges(2:end);
    
    %%
    ctx_up_isskipped{d} = nan(length(up_trans_inds8),1);
    for ii = 1:length(rt2_downs{d})
        cur_skipped = mp_downskipped{d}.inds{rt2_downs{d}(ii)};
        ctx_up_isskipped{d}(cur_skipped) = 1;
    end
    ctx_down_isskipped{d} = nan(length(down_trans_inds8),1);
    for ii = 1:length(rt2_ups{d})
        cur_skipped = mp_upskipped{d}.inds{rt2_ups{d}(ii)};
        ctx_down_isskipped{d}(cur_skipped) = 1;
    end
    
    %%
    [ctx_up_amps{d},down_amps] = get_state_amplitudes(lf8_lf,up_trans_inds8,down_trans_inds8);
    [ctx_up_amps_hf{d},down_amps_hf] = get_state_amplitudes(lf8_hf,up_trans_inds8,down_trans_inds8);
    ctx_up_durs{d} = lfp_state_durations{d}{2};

    binned_ctx_lfpow = zeros(length(state_axis),1);
    binned_ctx_hfpow = zeros(length(state_axis),1);
    binned_ctx_hfpowr = zeros(length(state_axis),1);
    binned_ctx_up_amps = zeros(length(state_axis),1);
    binned_ctx_up_ampshf = zeros(length(state_axis),1);
    binned_ctx_up_durs = zeros(length(state_axis),1);
    binned_ctx_down_durs = zeros(length(state_axis),1);
    binned_ctx_dc = zeros(length(state_axis),1);
    for ii = 1:length(state_axis)
       cur_set = find(up_trans_inds8/Fsd >= state_edges(ii) & up_trans_inds8/Fsd <= state_edges(ii+1));
       cur_inds = find(t_axis >= state_edges(ii) & t_axis <= state_edges(ii+1));
       binned_ctx_lfpow(ii) = std(lf8_lf(cur_inds));
       binned_ctx_hfpow(ii) = std(lf8_hf(cur_inds));
       binned_ctx_hfpowr(ii) = std(lf8_hfr(cur_inds));
       if ~isempty(cur_set)
          binned_ctx_up_amps(ii) = nanmean(ctx_up_amps{d}(cur_set));
          binned_ctx_up_ampshf(ii) = nanmean(ctx_up_amps_hf{d}(cur_set));
          binned_ctx_up_durs(ii) = nanmean(lf8_updurs(cur_set));
          binned_ctx_down_durs(ii) = nanmean(lf8_downdurs(cur_set));
          binned_ctx_dc(ii) = nanmean(lf8_dutycycs(cur_set));
       elseif ii > 1
          binned_ctx_up_amps(ii) = binned_ctx_up_amps(ii-1);
          binned_ctx_up_ampshf(ii) = binned_ctx_up_ampshf(ii-1);
          binned_ctx_up_durs(ii) = binned_ctx_up_durs(ii-1);
          binned_ctx_down_durs(ii) = binned_ctx_down_durs(ii-1);
          binned_ctx_dc(ii) = binned_ctx_dc(ii-1);
       end
    end
    %%
    all_up_times = up_trans_inds8/Fsd;
    rt2_down_times = up_trans_inds8(ctx_up_isskipped{d}==1)/Fsd;
    all_down_times = down_trans_inds8/Fsd;
    rt2_up_times = down_trans_inds8(ctx_down_isskipped{d}==1)/Fsd;
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
    
    [rt2_ups_acorr(d,:),lags] = xcov(binned_pup_frac,maxlag,'coeff');
    [rt2_downs_acorr(d,:),lags] = xcov(binned_pdown_frac,maxlag,'coeff');
    [ctx_upamp_acorr(d,:),lags] = xcov(binned_ctx_up_amps,maxlag,'coeff');
    [ctx_upamphf_acorr(d,:),lags] = xcov(binned_ctx_up_ampshf,maxlag,'coeff');
    [ctx_updur_acorr(d,:),lags] = xcov(binned_ctx_up_durs,maxlag,'coeff');
    [ctx_downdur_acorr(d,:),lags] = xcov(binned_ctx_down_durs,maxlag,'coeff');
    [ctx_lfpow_acorr(d,:),lags] = xcov(binned_ctx_lfpow,maxlag,'coeff');
    [ctx_hfpow_acorr(d,:),lags] = xcov(binned_ctx_hfpow,maxlag,'coeff');
    [ctx_hfpowr_acorr(d,:),lags] = xcov(binned_ctx_hfpowr,maxlag,'coeff');
    
    [rt2_updowns_xcorr(d,:),lags] = xcov(binned_pup_frac,binned_pdown_frac,maxlag,'coeff');
    [rt2down_upamp_xcorr(d,:),lags] = xcov(binned_ctx_up_amps,binned_pdown_frac,maxlag,'coeff');
    [rt2down_upamphf_xcorr(d,:),lags] = xcov(binned_ctx_up_ampshf,binned_pdown_frac,maxlag,'coeff');
    [rt2down_updur_xcorr(d,:),lags] = xcov(binned_ctx_up_durs,binned_rt2_downs,maxlag,'coeff');
    [rt2down_downdur_xcorr(d,:),lags] = xcov(binned_ctx_down_durs,binned_pdown_frac,maxlag,'coeff');
    [rt2down_dc_xcorr(d,:),lags] = xcov(binned_ctx_dc,binned_pdown_frac,maxlag,'coeff');
    [rt2down_lfpow_xcorr(d,:),lags] = xcov(binned_ctx_lfpow,binned_pdown_frac,maxlag,'coeff');
    [rt2down_hfpow_xcorr(d,:),lags] = xcov(binned_ctx_hfpow,binned_pdown_frac,maxlag,'coeff');
    [rt2down_hfpowr_xcorr(d,:),lags] = xcov(binned_ctx_hfpowr,binned_pdown_frac,maxlag,'coeff');
 
    [rt2up_upamp_xcorr(d,:),lags] = xcov(binned_ctx_up_amps,binned_pup_frac,maxlag,'coeff');
    [rt2up_upamphf_xcorr(d,:),lags] = xcov(binned_ctx_up_ampshf,binned_pup_frac,maxlag,'coeff');
    [upamp_upamphf_xcorr(d,:),lags] = xcov(binned_ctx_up_ampshf,binned_ctx_up_amps,maxlag,'coeff');
    [rt2up_updur_xcorr(d,:),lags] = xcov(binned_ctx_up_durs,binned_pup_frac,maxlag,'coeff');
    [rt2up_downdur_xcorr(d,:),lags] = xcov(binned_ctx_down_durs,binned_pup_frac,maxlag,'coeff');
    [rt2up_dc_xcorr(d,:),lags] = xcov(binned_ctx_dc,binned_pup_frac,maxlag,'coeff');
    [rt2up_lfpow_xcorr(d,:),lags] = xcov(binned_ctx_lfpow,binned_pup_frac,maxlag,'coeff');
    [rt2up_hfpow_xcorr(d,:),lags] = xcov(binned_ctx_hfpow,binned_pup_frac,maxlag,'coeff');
    [rt2up_hfpowr_xcorr(d,:),lags] = xcov(binned_ctx_hfpowr,binned_pup_frac,maxlag,'coeff');
    [lfpow_hfpowr_xcorr(d,:),lags] = xcov(binned_ctx_hfpowr,binned_ctx_lfpow,maxlag,'coeff');

end

%%
cd C:\WC_Germany\persistent_downs\
save new_pdown_acorr_analysis lags rt2* upamp* ctx_* fract_* data_type used_dirs state_bin_size

%%
cd C:\persDowns_paper\Figs\

 min_nstates = 15;
l3mec = find(strcmp(data_type(used_dirs),'L3MEC'));
l3lec = find(strcmp(data_type(used_dirs),'L3LEC'));

min_pdown = 0.05;
min_pup = 0.05;
used_pdowns = l3mec(fract_rt2_downs(l3mec) >= min_pdown);
used_pups = l3mec(fract_rt2_ups(l3mec) >= min_pup);
used_both = intersect(used_pdowns,used_pups);

zlag = find(lags == 0);
rt2_ups_acorr(:,zlag) = nan;
rt2_downs_acorr(:,zlag) = nan;
ctx_upamp_acorr(:,zlag) = nan;
ctx_upamphf_acorr(:,zlag) = nan;
ctx_updur_acorr(:,zlag) = nan;
ctx_downdur_acorr(:,zlag) = nan;
ctx_lfpow_acorr(:,zlag) = nan;
ctx_hfpow_acorr(:,zlag) = nan;

%%
figure
h1=shadedErrorBar(lags*state_bin_size,nanmean(rt2_ups_acorr(used_pups,:)),nanstd(rt2_ups_acorr(used_pups,:))/sqrt(length(used_pups)))
hold on
h2=shadedErrorBar(lags*state_bin_size,nanmean(rt2_downs_acorr(used_pdowns,:)),nanstd(rt2_downs_acorr(used_pdowns,:))/sqrt(length(used_pdowns)),{'color','r'})
legend([h1.mainLine h2.mainLine],{'Pers ups','Pers downs'});
xlabel('Time (s)','fontsize',14)
ylabel('Auto-Correlation','fontsize',14);
set(gca,'fontsize',12,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);

figure
shadedErrorBar(lags*state_bin_size,nanmean(rt2_updowns_xcorr(used_both,:)),nanstd(rt2_updowns_xcorr(used_both,:))/sqrt(length(used_both)));
xlabel('Time (s)','fontsize',14)
ylabel('Cross-Correlation','fontsize',14);
set(gca,'fontsize',12,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);

%%
figure
shadedErrorBar(lags*state_bin_size,nanmean(ctx_lfpow_acorr(l3mec,:)),nanstd(ctx_lfpow_acorr(l3mec,:))/sqrt(length(l3mec)))
hold on
shadedErrorBar(lags*state_bin_size,nanmean(ctx_hfpow_acorr(l3mec,:)),nanstd(ctx_hfpow_acorr(l3mec,:))/sqrt(length(l3mec)),{'color','r'})

figure
h1=shadedErrorBar(lags*state_bin_size,nanmean(ctx_upamp_acorr(l3mec,:)),nanstd(ctx_upamp_acorr(l3mec,:))/sqrt(length(l3mec)))
hold on
h2=shadedErrorBar(lags*state_bin_size,nanmean(ctx_upamphf_acorr(l3mec,:)),nanstd(ctx_upamphf_acorr(l3mec,:))/sqrt(length(l3mec)),{'color','r'})
legend([h1.mainLine h2.mainLine],{'LF amplitude','HF amplitude'});
xlabel('Time (s)','fontsize',14)
ylabel('Auto-Correlation','fontsize',14);
set(gca,'fontsize',12,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);

figure
shadedErrorBar(lags*state_bin_size,nanmean(ctx_updur_acorr(l3mec,:)),nanstd(ctx_updur_acorr(l3mec,:))/sqrt(length(l3mec)))
hold on
shadedErrorBar(lags*state_bin_size,nanmean(ctx_downdur_acorr(l3mec,:)),nanstd(ctx_downdur_acorr(l3mec,:))/sqrt(length(l3mec)),{'color','r'})

%%
figure
shadedErrorBar(lags*state_bin_size,nanmean(rt2down_upamp_xcorr(used_pdowns,:)),nanstd(rt2down_upamp_xcorr(used_pdowns,:))/sqrt(length(used_pdowns)))
hold on
shadedErrorBar(lags*state_bin_size,nanmean(rt2down_upamphf_xcorr(used_pdowns,:)),nanstd(rt2down_upamphf_xcorr(used_pdowns,:))/sqrt(length(used_pdowns)),{'color','r'})

figure
h1=shadedErrorBar(lags*state_bin_size,nanmean(rt2down_lfpow_xcorr(used_pdowns,:)),nanstd(rt2down_lfpow_xcorr(used_pdowns,:))/sqrt(length(used_pdowns)))
hold on
shadedErrorBar(lags*state_bin_size,nanmean(rt2down_hfpow_xcorr(used_pdowns,:)),nanstd(rt2down_hfpow_xcorr(used_pdowns,:))/sqrt(length(used_pdowns)),{'color','r'})
h2=shadedErrorBar(lags*state_bin_size,nanmean(rt2down_hfpowr_xcorr(used_pdowns,:)),nanstd(rt2down_hfpowr_xcorr(used_pdowns,:))/sqrt(length(used_pdowns)),{'color','b'})
legend([h1.mainLine h2.mainLine],{'LF power','HF power'});
xlabel('Time (s)','fontsize',14)
ylabel('Cross-Correlation','fontsize',14);
set(gca,'fontsize',12,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);

figure
h1=shadedErrorBar(lags*state_bin_size,nanmean(rt2down_updur_xcorr(used_pdowns,:)),nanstd(rt2down_updur_xcorr(used_pdowns,:))/sqrt(length(used_pdowns)))
hold on
shadedErrorBar(lags*state_bin_size,nanmean(rt2down_downdur_xcorr(used_pdowns,:)),nanstd(rt2down_downdur_xcorr(used_pdowns,:))/sqrt(length(used_pdowns)),{'color','r'})
h2=shadedErrorBar(lags*state_bin_size,nanmean(rt2down_dc_xcorr(used_pdowns,:)),nanstd(rt2down_dc_xcorr(used_pdowns,:))/sqrt(length(used_pdowns)),{'color','b'})
legend([h1.mainLine h2.mainLine],{'LF power','HF power'});
xlabel('Time (s)','fontsize',14)
ylabel('Cross-Correlation','fontsize',14);
set(gca,'fontsize',12,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);

%%
figure
h1=shadedErrorBar(lags*state_bin_size,nanmean(rt2up_upamp_xcorr(used_pups,:)),nanstd(rt2up_upamp_xcorr(used_pups,:))/sqrt(length(used_pups)))
hold on
h2=shadedErrorBar(lags*state_bin_size,nanmean(rt2up_upamphf_xcorr(used_pups,:)),nanstd(rt2up_upamphf_xcorr(used_pups,:))/sqrt(length(used_pups)),{'color','r'})
legend([h1.mainLine h2.mainLine],{'LF amplitude','HF amplitude'});
xlabel('Time (s)','fontsize',14)
ylabel('Cross-Correlation','fontsize',14);
set(gca,'fontsize',12,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);

figure
h1=shadedErrorBar(lags*state_bin_size,nanmean(rt2up_lfpow_xcorr(used_pups,:)),nanstd(rt2up_lfpow_xcorr(used_pups,:))/sqrt(length(used_pups)))
hold on
% h2=shadedErrorBar(lags*state_bin_size,nanmean(rt2up_hfpow_xcorr(used_pups,:)),nanstd(rt2up_hfpow_xcorr(used_pups,:))/sqrt(length(used_pups)),{'color','r'})
h2=shadedErrorBar(lags*state_bin_size,nanmean(rt2up_hfpowr_xcorr(used_pups,:)),nanstd(rt2up_hfpowr_xcorr(used_pups,:))/sqrt(length(used_pups)),{'color','b'})
legend([h1.mainLine h2.mainLine],{'LF power','HF power'});
xlabel('Time (s)','fontsize',14)
ylabel('Cross-Correlation','fontsize',14);
set(gca,'fontsize',12,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);


figure
h1=shadedErrorBar(lags*state_bin_size,nanmean(rt2up_updur_xcorr(used_pups,:)),nanstd(rt2up_updur_xcorr(used_pups,:))/sqrt(length(used_pups)))
hold on
shadedErrorBar(lags*state_bin_size,nanmean(rt2up_downdur_xcorr(used_pups,:)),nanstd(rt2up_downdur_xcorr(used_pups,:))/sqrt(length(used_pups)),{'color','r'})
h2=shadedErrorBar(lags*state_bin_size,nanmean(rt2up_dc_xcorr(used_pups,:)),nanstd(rt2up_dc_xcorr(used_pups,:))/sqrt(length(used_pups)),{'color','b'})
legend([h1.mainLine h2.mainLine],{'LF power','HF power'});
xlabel('Time (s)','fontsize',14)
ylabel('Cross-Correlation','fontsize',14);
set(gca,'fontsize',12,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);

%%
figure
h1=shadedErrorBar(lags*state_bin_size,nanmean(lfpow_hfpowr_xcorr(l3mec,:)),nanstd(lfpow_hfpowr_xcorr(l3mec,:))/sqrt(length(l3mec)))
hold on
h2=shadedErrorBar(lags*state_bin_size,nanmean(upamp_upamphf_xcorr(l3mec,:)),nanstd(upamp_upamphf_xcorr(l3mec,:))/sqrt(length(l3mec)),{'color','r'})
