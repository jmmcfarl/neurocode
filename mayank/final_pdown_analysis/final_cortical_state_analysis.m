%%
clear all
addpath('C:/Code/general/');
addpath('C:/Code/general_functions/');
addpath('C:/WC_Germany/hsmm_state_detection/');
addpath('C:/WC_Germany/persistent_downs/');
addpath('C:/WC_Germany/sven_thomas_combined/');
addpath('C:/WC_Germany/parietal_cortical_2010/');
addpath('C:/WC_Germany/final_pdown_analysis/');
addpath(genpath('C:/Code/figuremaker/'));

fig_dir = 'C:\WC_Germany\final_pdown_analysis\figures\';

load C:/WC_Germany/final_pdown_analysis/compiled_data.mat
load C:/WC_Germany/final_pdown_analysis/mua_classification.mat

min_rec_dur = 500; %in sec
data_ids = [data(:).id];
used_dirs = find([data(:).ep] > min_rec_dur & [data(:).dp] > min_rec_dur);

% used_dirs(ismember(data_ids(used_dirs),no_cell) | ismember(data_ids(used_dirs),unclear_uds)) = [];
used_dirs(ismember(data_ids(used_dirs),unclear_uds)) = [];

data = data(used_dirs);

peak_hpcmua_loc = peak_hpcmua_loc(used_dirs);
peak_hpcmua_rate = peak_hpcmua_rate(used_dirs);
usable_mua = usable_mua(used_dirs,:);

load C:/WC_Germany/final_pdown_analysis/fin_pdown_core_analysis.mat
if length(core_data) ~= length(data)
    error('Data mismatch');
end

min_hpcrate = 1;
usable_hpc_mua = ~isnan(peak_hpcmua_loc) & peak_hpcmua_rate >= min_hpcrate;

%% parameters
raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
hcf = 10;
lcf = 0.05;
lcf_hf = 30;
hcf_hf = 100;
hcf_sm = 0.025;
rate_sm = round(Fsd*0.05);

% trans_buffer = round(Fsd*0.1);
trans_buffer = round(Fsd*0);

state_amp_xx = linspace(-3,3,100);
% state_amp_mua_xx = linspace(0,5,100);
state_amp_mua_xx = linspace(-4,4,100);
%%
for d = 1:length(used_dirs)
    cd(data(d).dir)
    pwd
    
    load ./used_data lf7 wcv_minus_spike
    [lfp_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf]);
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
%     lfp_hf = get_hf_features(lf7,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
    lfp_hf = get_hf_lpower(lf7,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
    lfp_hf = zscore(lfp_hf);
    
    load ./sync_times.mat
    synct_d = downsample(synct,dsf);
        
    load ./spike_time_jmm

    end_time = min(data(d).ep,data(d).dp);
    ep = find(t_axis >= end_time,1);
    if ~isempty(ep)
        synct_d(ep+1:end) = []; lfp_lf(ep+1:end) = []; wcv_lf(ep+1:end) = []; t_axis(ep+1:end) = []; lfp_hf(ep+1:end) = []; %hpc_hf(ep+1:end) = [];
    else
        ep = length(t_axis);
    end
    
        mp_spike_rate = hist(synct(spkid),synct_d)*Fsd;
    mp_spike_rate(end) = 0;

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
           
    ctx_up_cents = (up_trans_inds8 + down_trans_inds8)/2;
    mp_up_cents = (up_trans_inds + down_trans_inds)/2;
    nearest_mp_up = nan(size(ctx_up_cents));
    for ii = 1:length(ctx_up_cents)
       [~,bb] = min(abs(ctx_up_cents(ii)-mp_up_cents));
       nearest_mp_up(ii) = bb;
    end
    %%
    if exist('./mua_data3.mat','file')
        load ./mua_data3.mat
        
        usable_ctx_mua = find(usable_mua(d,6:end)); %set of cortical MUA that can be used

        ctx_mua_times = [];
        for ii = 1:length(usable_ctx_mua)
            ctx_mua_times = [ctx_mua_times mua_times{usable_ctx_mua(ii)+5}];
        end
        ctx_mua_times = sort(ctx_mua_times);
        ctx_mua_rate = hist(ctx_mua_times,synct_d)*Fsd;
        ctx_mua_rate([1 end]) = 0;
%         ctx_mua_rate = jmm_smooth_1d_cor(ctx_mua_rate,rate_sm);
        
        if isempty(usable_ctx_mua)
            ctx_mua_rate = nan(size(synct_d));
        end        
        if length(ctx_mua_rate) > length(synct_d)
            ctx_mua_rate = ctx_mua_rate(1:length(synct_d));
        end
    else
        ctx_mua_rate = nan(size(synct_d));
    end
    avg_ctx_mua_rate(d) = nanmean(ctx_mua_rate);
    std_ctx_mua_rate(d) = nanstd(ctx_mua_rate);
    ctx_mua_rate = zscore(ctx_mua_rate);
%         ctx_mua_rate = ctx_mua_rate/avg_ctx_mua_rate(d);

    mp_spike_rate = jmm_smooth_1d_cor(mp_spike_rate,rate_sm);
    mp_spike_rate = zscore(mp_spike_rate);
    %%    
%     [ctx_up_amps,ctx_down_amps] = get_state_amplitudes(lfp_lf,up_trans_inds8,down_trans_inds8);
%     [ctx_up_amps_hf,ctx_down_amps_hf] = get_state_amplitudes(lfp_hf,up_trans_inds8,down_trans_inds8);
%     [ctx_up_amps_mua,ctx_down_amps_mua] = get_state_amplitudes(ctx_mua_rate,up_trans_inds8,down_trans_inds8);
    [ctx_up_amps,ctx_down_amps] = get_state_amplitudes_mean(lfp_lf,up_trans_inds8,down_trans_inds8,trans_buffer);
    [ctx_up_amps_hf,ctx_down_amps_hf] = get_state_amplitudes_mean(lfp_hf,up_trans_inds8,down_trans_inds8,trans_buffer);
    [ctx_up_amps_mua,ctx_down_amps_mua] = get_state_amplitudes_mean(ctx_mua_rate,up_trans_inds8,down_trans_inds8,trans_buffer);
    [mp_up_amps_spkrate,mp_down_amps_spkrate] = get_state_amplitudes_mean(mp_spike_rate,up_trans_inds,down_trans_inds,trans_buffer);
    up_avg_spkrate = mean(mp_up_amps_spkrate);
    mp_up_amps_spkrate = detrend(mp_up_amps_spkrate) + up_avg_spkrate;
    
    cfun_data(d).ctx_down_amps = ctx_down_amps;
    cfun_data(d).ctx_up_amps = ctx_up_amps;
    cfun_data(d).ctx_down_amps_hf = ctx_down_amps_hf;
    cfun_data(d).ctx_up_amps_hf = ctx_up_amps_hf;
    cfun_data(d).ctx_down_amps_mua = ctx_down_amps_mua;
    cfun_data(d).ctx_up_amps_mua = ctx_up_amps_mua;
    
    cfun_data(d).ctx_down_amp_dist = ksdensity(ctx_down_amps,state_amp_xx);
    cfun_data(d).ctx_up_amp_dist = ksdensity(ctx_up_amps,state_amp_xx);
    cfun_data(d).ctx_down_amp_hf_dist = ksdensity(ctx_down_amps_hf,state_amp_xx);
    cfun_data(d).ctx_up_amp_hf_dist = ksdensity(ctx_up_amps_hf,state_amp_xx);

%     cur_down_amp_dist = histc(ctx_down_amps,state_amp_xx);
%     cfun_data(d).ctx_down_amp_cdist = cumsum(cur_down_amp_dist)/sum(cur_down_amp_dist);
%     cur_up_amp_dist = histc(ctx_up_amps,state_amp_xx);
%     cfun_data(d).ctx_up_amp_cdist = cumsum(cur_up_amp_dist)/sum(cur_up_amp_dist);
%     cur_down_amp_hf_dist = histc(ctx_down_amps_hf,state_amp_xx);
%     cfun_data(d).ctx_down_amp_hf_cdist = cumsum(cur_down_amp_hf_dist)/sum(cur_down_amp_hf_dist);
%     cur_up_amp_hf_dist = histc(ctx_up_amps_hf,state_amp_xx);
%     cfun_data(d).ctx_up_amp_hf_cdist = cumsum(cur_up_amp_hf_dist)/sum(cur_up_amp_hf_dist);
% 
%     cur_down_amp_mua_dist = histc(ctx_down_amps_mua,state_amp_mua_xx);
%     cfun_data(d).ctx_down_amp_mua_cdist = cumsum(cur_down_amp_mua_dist)/sum(cur_down_amp_mua_dist);
%     cur_up_amp_mua_dist = histc(ctx_up_amps_mua,state_amp_mua_xx);
%     cfun_data(d).ctx_up_amp_mua_cdist = cumsum(cur_up_amp_mua_dist)/sum(cur_up_amp_mua_dist);
    
    
%     cfun_data(d).ctx_down_amp_mua_dist = ksdensity(ctx_down_amps_mua,state_amp_mua_xx);
%     cfun_data(d).ctx_up_amp_mua_dist = ksdensity(ctx_up_amps_mua,state_amp_mua_xx);
    
    %up-transition lag of MP state corresponding to each cortical up-trans
    uset = find(~isnan(corresp_mp_upinds));
    corresp_mp_uplag = nan(length(up_trans_inds8),1);
    corresp_mp_uplag(uset) = core_data(d).mp_uplags(corresp_mp_upinds(uset))/Fsd;
    cfun_data(d).corresp_mp_uplag = corresp_mp_uplag;
    
    corresp_mp_spkrate = nan(length(up_trans_inds8),1);
%     corresp_mp_spkrate(uset) = mp_up_amps_spkrate(corresp_mp_upinds(uset));
    corresp_mp_spkrate = mp_up_amps_spkrate(nearest_mp_up);
    cfun_data(d).corresp_mp_spkrate = corresp_mp_spkrate;
    
    %down-transition lag of MP state corresponding to each cortical
    %up-trans
    uset = find(~isnan(corresp_mp_upinds));
    corresp_mp_downlag = nan(length(up_trans_inds8),1);
    corresp_mp_downlag(uset) = core_data(d).mp_downlags(corresp_mp_upinds(uset))/Fsd;
    cfun_data(d).corresp_mp_downlag = corresp_mp_downlag;

    ctx_up_isskipped = nan(length(up_trans_inds8),1);
    ctx_up_isskipped(core_data(d).skipped_lfp_ups) = 1;
    ctx_up_isskipped(core_data(d).non_skipped_lfp_ups) = 0;    
    %don't use any cortical up trans that happen during an MP up state
    ctx_ups_mp_up = find(mp_state_vec(up_trans_inds8) == 1);
    ctx_up_isskipped(ctx_ups_mp_up) = nan;
    cfun_data(d).ctx_up_isskipped = ctx_up_isskipped;
    
    ctx_down_isskipped = nan(length(down_trans_inds8),1);
    ctx_down_isskipped(core_data(d).skipped_lfp_downs) = 1;
    ctx_down_isskipped(core_data(d).non_skipped_lfp_downs) = 0;    
    %don't use any cortical down-trans that happen when the MP is in a down
    %state
    ctx_downs_mp_down = find(mp_state_vec(down_trans_inds8) == 0);
    ctx_down_isskipped(ctx_ups_mp_up) = nan;
    cfun_data(d).ctx_down_isskipped = ctx_down_isskipped;
    
    %duration of preceding cortical Down state
    prev_ctx_downdur = nan(length(up_trans_inds8),1);
    prev_ctx_downdur(2:end) = core_data(d).lfp_down_durs(1:end-1);
    cfun_data(d).prev_ctx_downdur = prev_ctx_downdur;

    %time since the last MP down-transition (from current cortical up)
    TSL_mp_down = nan(size(up_trans_inds8));
    for i = 1:length(TSL_mp_down)
        prev_mp_down = find(down_trans_inds < up_trans_inds8(i),1,'last');
        if ~isempty(prev_mp_down)
            TSL_mp_down(i) = (up_trans_inds8(i) - down_trans_inds(prev_mp_down))/Fsd;
        end
    end
    cfun_data(d).TSL_mp_down = TSL_mp_down;
    
    %time since the last MP up-transition (from current cortical down)
    TSL_mp_up = nan(size(up_trans_inds8));
    for i = 1:length(TSL_mp_up)
        prev_mp_up = find(up_trans_inds < down_trans_inds8(i),1,'last');
        if ~isempty(prev_mp_up)
            TSL_mp_up(i) = (down_trans_inds8(i) - up_trans_inds(prev_mp_up))/Fsd;
        end
    end
    cfun_data(d).TSL_mp_up = TSL_mp_up;   
   
end

%%
cd C:\WC_Germany\final_pdown_analysis\
save final_cortical_state_data_fin_nobuff_wc5 cfun_data avg_ctx_mua_rate state_amp_*

%%
cd C:\WC_Germany\final_pdown_analysis\
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

min_mua_rate = 1;
used_ctx_mua = avg_ctx_mua_rate >= min_mua_rate;

%%
clear glm_pvals
for ii = 1:length(cfun_data)
   X = [cfun_data(ii).TSL_mp_down(:) cfun_data(ii).ctx_up_amps_hf(:)];
  Y = [cfun_data(ii).ctx_up_isskipped(:)];
   [b,dev,stats] = glmfit(X,Y,'binomial');
   glm_pvals_pdowns(ii,:) = stats.p(2:end);
   glm_beta_pdowns(ii,:) = b(2:end);
   X = [cfun_data(ii).TSL_mp_up(:) cfun_data(ii).ctx_down_amps_hf(:)];
  Y = [cfun_data(ii).ctx_down_isskipped(:)];
   [b,dev,stats] = glmfit(X,Y,'binomial');
   glm_pvals_pups(ii,:) = stats.p(2:end);
   glm_beta_pups(ii,:) = b(2:end);
%    yhat = glmval(b,X,'logit');

   X = [cfun_data(ii).ctx_up_amps_hf(:)];
  Y = [cfun_data(ii).corresp_mp_spkrate(:)];
  uset = find(~isnan(X) & ~isnan(Y));
  if isempty(uset)
      corr_uphf_spkrate(ii) = nan;
      corr_uphf_spkrate_p(ii) = nan;
  else
  [corr_uphf_spkrate(ii),corr_uphf_spkrate_p(ii)] = corr(X(uset),Y(uset),'type','spearman');
  end
end

%% Persistence vs Time since prev MP trans
close all
n_bins = 15;
state_dur_xx = linspace(0.1,5,100);

used_data = l3mec_pdowns;

X_fac = {cfun_data(:).TSL_mp_down};
Y_fac = {cfun_data(:).ctx_up_isskipped};

clear avg_fun std_fun N_fun sem_fun avg_xvals US_corr US_pval avg_tsd_dist
for ii = 1:length(used_data)
    if length(X_fac{used_data(ii)}) ~= length(Y_fac{used_data(ii)})
        error('!');
    end
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [US_corr(ii),US_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
    avg_tsd_dist(ii,:) = ksdensity(X_fac{used_data(ii)},state_dur_xx,'support','positive');
end

n_nonnan = sum(~isnan(avg_fun));
h = figure(); 
subplot(2,1,1);
hold on
xax = nanmean(avg_xvals);
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan));
plot(xax,nanmean(avg_fun),'ko')

used_data = l3mec_pups;
X_fac = {cfun_data(:).TSL_mp_up};
Y_fac = {cfun_data(:).ctx_down_isskipped};

clear avg_fun std_fun N_fun sem_fun avg_xvals DS_corr DS_pval avg_tsu_dist
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [DS_corr(ii),DS_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
    avg_tsu_dist(ii,:) = ksdensity(X_fac{used_data(ii)},state_dur_xx,'support','positive');
end

n_nonnan = sum(~isnan(avg_fun));
xax = nanmean(avg_xvals);
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan),{'color','r'});
plot(xax,nanmean(avg_fun),'ko')
xlabel('Time since MP transition (s)');
ylabel('Prob. skipped');

xlim([0.1 5]); 
set(gca,'xscale','log');

% f2 = figure(); hold on
subplot(2,1,2);
hold on
shadedErrorBar(state_dur_xx,nanmean(avg_tsd_dist),nanstd(avg_tsd_dist)./sqrt(sum(~isnan(avg_tsd_dist))));
shadedErrorBar(state_dur_xx,nanmean(avg_tsu_dist),nanstd(avg_tsu_dist)./sqrt(sum(~isnan(avg_tsu_dist))),{'color','r'});
xlim([0.1 5]); 
set(gca,'xscale','log');
xlabel('Time since MP transition (s)');
ylabel('Probability');

%%
fig_width = 4; rel_height = 1.6;
figufy(h);
fname = [fig_dir 'prob_skip_vs_tsltrans_fin_5s.pdf'];
exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%% HFpow BOtH
close all
n_bins = 15;

used_data = l3mec_pdowns;

X_fac = {cfun_data(:).ctx_up_amps_hf};
Y_fac = {cfun_data(:).ctx_up_isskipped};
clear avg_fun std_fun N_fun sem_fun avg_xvals XY_corr XY_pval
avg_upamp_dist = nan(length(used_data),length(state_amp_xx));
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr(ii),XY_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
    avg_upamp_dist(ii,:) = cfun_data(used_data(ii)).ctx_up_amp_hf_dist;
end
n_nonnan = sum(~isnan(avg_fun));

h = figure(); 
subplot(2,1,1);hold on
xax = nanmean(avg_xvals);
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan));
plot(xax,nanmean(avg_fun),'ko')

subplot(2,1,2); hold on
n_nonnan = sum(~isnan(avg_upamp_dist));
shadedErrorBar(state_amp_xx,nanmean(avg_upamp_dist),nanstd(avg_upamp_dist)./sqrt(n_nonnan));


used_data = l3mec_pups;

X_fac = {cfun_data(:).ctx_down_amps_hf};
Y_fac = {cfun_data(:).ctx_down_isskipped};
clear avg_fun std_fun N_fun sem_fun avg_xvals XY_corr XY_pval
avg_downamp_dist = nan(length(used_data),length(state_amp_xx));
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr(ii),XY_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
    avg_downamp_dist(ii,:) = cfun_data(used_data(ii)).ctx_down_amp_hf_dist;
end
n_nonnan = sum(~isnan(avg_fun));

% figure(h);
subplot(2,1,1);
xax = nanmean(avg_xvals);
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan),{'color','r'});
plot(xax,nanmean(avg_fun),'ko')
xlabel('Ctx state amplitude (z)');
ylabel('Prob. skipped');

% figure(h2); hold on
subplot(2,1,2);
n_nonnan = sum(~isnan(avg_downamp_dist));
shadedErrorBar(state_amp_xx,nanmean(avg_downamp_dist),nanstd(avg_downamp_dist)./sqrt(n_nonnan),{'color','r'});
xlabel('Ctx state amplitude (z)');
ylabel('Probability');

subplot(2,1,1);
xlim([-1.5 2]);
subplot(2,1,2);
xlim([-1.5 2]);

%%
fig_width = 4; rel_height = 1.6;
figufy(h);
% fname = [fig_dir 'prob_skip_vs_lfamp.pdf'];
fname = [fig_dir 'prob_skip_vs_hfamp_fin_5s.pdf'];
exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%% LF amp BOtH
close all
n_bins = 15;

used_data = l3mec_pdowns;

X_fac = {cfun_data(:).ctx_up_amps};
Y_fac = {cfun_data(:).ctx_up_isskipped};
clear avg_fun std_fun N_fun sem_fun avg_xvals XY_corr XY_pval
avg_upamp_dist = nan(length(used_data),length(state_amp_xx));
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr(ii),XY_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
    avg_upamp_dist(ii,:) = cfun_data(used_data(ii)).ctx_up_amp_hf_dist;
end
n_nonnan = sum(~isnan(avg_fun));

h = figure(); 
subplot(2,1,1);hold on
xax = nanmean(avg_xvals);
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan));
plot(xax,nanmean(avg_fun),'ko')

subplot(2,1,2); hold on
n_nonnan = sum(~isnan(avg_upamp_dist));
shadedErrorBar(state_amp_xx,nanmean(avg_upamp_dist),nanstd(avg_upamp_dist)./sqrt(n_nonnan));


used_data = l3mec_pups;

X_fac = {cfun_data(:).ctx_down_amps};
Y_fac = {cfun_data(:).ctx_down_isskipped};
clear avg_fun std_fun N_fun sem_fun avg_xvals XY_corr XY_pval
avg_downamp_dist = nan(length(used_data),length(state_amp_xx));
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr(ii),XY_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
    avg_downamp_dist(ii,:) = cfun_data(used_data(ii)).ctx_down_amp_hf_dist;
end
n_nonnan = sum(~isnan(avg_fun));

% figure(h);
subplot(2,1,1);
xax = nanmean(avg_xvals);
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan),{'color','r'});
plot(xax,nanmean(avg_fun),'ko')
xlabel('Ctx state amplitude (z)');
ylabel('Prob. skipped');

% figure(h2); hold on
subplot(2,1,2);
n_nonnan = sum(~isnan(avg_downamp_dist));
shadedErrorBar(state_amp_xx,nanmean(avg_downamp_dist),nanstd(avg_downamp_dist)./sqrt(n_nonnan),{'color','r'});
xlabel('Ctx state amplitude (z)');
ylabel('Probability');

subplot(2,1,1);
xlim([-1.5 2.5]);
subplot(2,1,2);
xlim([-1.5 2.5]);
%%
fig_width = 4; rel_height = 1.6;
figufy(h);
fname = [fig_dir 'prob_skip_vs_lfamp_5s.pdf'];
exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
used_data = l3lec;

n_bins = 20;
X_fac = {cfun_data(:).ctx_up_amps_hf};
Y_fac = {cfun_data(:).corresp_mp_spkrate};
% X_fac = {cfun_data(:).ctx_up_amps_hf};
% Y_fac = {cfun_data(:).ctx_up_amps_mua};
xname = 'Ctx up durs (s)';
yname = 'Prob up skipped';

clear avg_fun std_fun N_fun sem_fun avg_xvals XY_corr XY_pval
for ii = 1:length(used_data)
    if sum(~isnan(Y_fac{used_data(ii)})) > 0
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr(ii),XY_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
    else
       avg_fun(ii,:) = nan;
       std_fun(ii,:) = nan;
    end
end

fprintf('Avg correlation: %.4f, p=%.4g\n',mean(XY_corr),signrank(XY_corr));
nsig = sum(XY_pval < 0.05);
fprintf('%d of %d cells sig\n',nsig,length(XY_pval));

n_nonnan = sum(~isnan(avg_fun));
h = figure(); hold on
xax = nanmean(avg_xvals);
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan));
plot(xax,nanmean(avg_fun),'ko')
xlabel(xname);  ylabel(yname);

%%

% used_ctx_mua = ~isnan(avg_ctx_mua_rate);
% used_data = l3mec([core_data(l3mec).fract_rt2_downs] > min_persfrac & used_ctx_mua(l3mec));
% used_data = l3mec([core_data(l3mec).fract_rt2_ups] > min_persfrac & used_ctx_mua(l3mec));


n_bins = 20;
% X_fac = {core_data(:).lfp_up_durs};
X_fac = {cfun_data(:).ctx_up_amps_hf};
Y_fac = {cfun_data(:).ctx_up_isskipped};
xname = 'Ctx up durs (s)';
yname = 'Prob up skipped';

clear avg_fun std_fun N_fun sem_fun avg_xvals XY_corr XY_pval
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr(ii),XY_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
end

fprintf('Avg correlation: %.4f, p=%.4g\n',mean(XY_corr),signrank(XY_corr));
nsig = sum(XY_pval < 0.05);
fprintf('%d of %d cells sig\n',nsig,length(XY_pval));

n_nonnan = sum(~isnan(avg_fun));
h = figure(); hold on
xax = nanmean(avg_xvals);
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan));
plot(xax,nanmean(avg_fun),'ko')
xlabel(xname);  ylabel(yname);

%% LFamp A VS PDOWN
used_data = l3mec_pdowns;

n_bins = 20;
X_fac = {cfun_data(:).ctx_up_amps};
Y_fac = {cfun_data(:).ctx_up_isskipped};
xname = 'Ctx up-state LF';
yname = 'Prob up skipped';

clear avg_fun std_fun N_fun sem_fun avg_xvals XY_corr XY_pval
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr(ii),XY_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
end

fprintf('Avg correlation: %.4f, p=%.4g\n',mean(XY_corr),signrank(XY_corr));
nsig = sum(XY_pval < 0.05);
fprintf('%d of %d cells sig\n',nsig,length(XY_pval));

n_nonnan = sum(~isnan(avg_fun));
h = figure(); hold on
xax = nanmean(avg_xvals);
% xax = pin_cents;
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan));
plot(xax,nanmean(avg_fun),'ko')
xlabel(xname);  ylabel(yname);
xlim(xax([1 end]));

figufy(h);
% fname = [fig_dir 'pers_ups_vs_pers_downs.pdf'];
% exportfig(fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);

%% HFpow A VS PDOWN
used_data = l3mec_pdowns;

n_bins = 20;
X_fac = {cfun_data(:).ctx_up_amps_hf};
Y_fac = {cfun_data(:).ctx_up_isskipped};
xname = 'Ctx up-state HF';
yname = 'Prob up skipped';

clear avg_fun std_fun N_fun sem_fun avg_xvals XY_corr XY_pval
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr(ii),XY_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
end

fprintf('Avg correlation: %.4f, p=%.4g\n',mean(XY_corr),signrank(XY_corr));
nsig = sum(XY_pval < 0.05);
fprintf('%d of %d cells sig\n',nsig,length(XY_pval));

n_nonnan = sum(~isnan(avg_fun));
% h = figure(); hold on
xax = nanmean(avg_xvals);
% xax = pin_cents;
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan),{'color','r'});
plot(xax,nanmean(avg_fun),'ko')
xlabel(xname);  ylabel(yname);
xlim(xax([1 end]));

figufy(h);
% fname = [fig_dir 'pers_ups_vs_pers_downs.pdf'];
% exportfig(fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);

%% MUA VS PDOWN
used_data = l3mec_pdowns(used_ctx_mua(l3mec_pdowns));

n_bins = 20;
X_fac = {cfun_data(:).ctx_up_amps_mua};
Y_fac = {cfun_data(:).ctx_up_isskipped};
xname = 'Ctx up-state MUA';
yname = 'Prob up skipped';

clear avg_fun std_fun N_fun sem_fun avg_xvals XY_corr XY_pval
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr(ii),XY_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
end

fprintf('Avg correlation: %.4f, p=%.4g\n',mean(XY_corr),signrank(XY_corr));
nsig = sum(XY_pval < 0.05);
fprintf('%d of %d cells sig\n',nsig,length(XY_pval));

n_nonnan = sum(~isnan(avg_fun));
% h = figure(); hold on
xax = nanmean(avg_xvals);
% xax = pin_cents;
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan));
plot(xax,nanmean(avg_fun),'ko')
xlabel(xname);  ylabel(yname);
xlim(xax([1 end]));

figufy(h);
% fname = [fig_dir 'pers_ups_vs_pers_downs.pdf'];
% exportfig(fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);

%% MP DOWN DUR VS PDOWN
close all
used_data = l3mec_pdowns;

n_bins = 20;
% X_fac = {cfun_data(:).prev_ctx_downdur};
X_fac = {cfun_data(:).TSL_mp_down};
Y_fac = {cfun_data(:).ctx_up_isskipped};
xname = 'Previous down duration';
yname = 'Prob up skipped';

clear avg_fun std_fun N_fun sem_fun avg_xvals XY_corr XY_pval
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr(ii),XY_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
end

fprintf('Avg correlation: %.4f, p=%.4g\n',mean(XY_corr),signrank(XY_corr));
nsig = sum(XY_pval < 0.05);
fprintf('%d of %d cells sig\n',nsig,length(XY_pval));

n_nonnan = sum(~isnan(avg_fun));
h = figure(); hold on
xax = nanmean(avg_xvals);
% xax = pin_cents;
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan));
plot(xax,nanmean(avg_fun),'ko')
xlabel(xname);  ylabel(yname);
figufy(h);
set(gca,'xscale','log');
xlim([xax(1) xax(end)]);
% fname = [fig_dir 'pers_ups_vs_pers_downs.pdf'];
% exportfig(fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);

%% MP UP Dur VS PUP
close all
used_data = l3mec_pups;

n_bins = 20;
% X_fac = {cfun_data(:).prev_ctx_downdur};
X_fac = {cfun_data(:).TSL_mp_up};
Y_fac = {cfun_data(:).ctx_down_isskipped};
xname = 'Previous down duration';
yname = 'Prob up skipped';

clear avg_fun std_fun N_fun sem_fun avg_xvals XY_corr XY_pval
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr(ii),XY_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
end

fprintf('Avg correlation: %.4f, p=%.4g\n',mean(XY_corr),signrank(XY_corr));
nsig = sum(XY_pval < 0.05);
fprintf('%d of %d cells sig\n',nsig,length(XY_pval));

n_nonnan = sum(~isnan(avg_fun));
h = figure(); hold on
xax = nanmean(avg_xvals);
% xax = pin_cents;
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan));
plot(xax,nanmean(avg_fun),'ko')
xlabel(xname);  ylabel(yname);
figufy(h);
set(gca,'xscale','log');
xlim([xax(1) xax(end)]);
% fname = [fig_dir 'pers_ups_vs_pers_downs.pdf'];
% exportfig(fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);


%% LFamp A VS PUP
used_data = l3mec_pups;

n_bins = 20;
X_fac = {cfun_data(:).ctx_down_amps};
Y_fac = {cfun_data(:).ctx_down_isskipped};
xname = 'Ctx down-state';
yname = 'Prob down skipped';

clear avg_fun std_fun N_fun sem_fun avg_xvals XY_corr XY_pval
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr(ii),XY_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
end

fprintf('Avg correlation: %.4f, p=%.4g\n',mean(XY_corr),signrank(XY_corr));
nsig = sum(XY_pval < 0.05);
fprintf('%d of %d cells sig\n',nsig,length(XY_pval));

n_nonnan = sum(~isnan(avg_fun));
h = figure(); hold on
xax = nanmean(avg_xvals);
% xax = pin_cents;
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan));
plot(xax,nanmean(avg_fun),'ko')
xlabel(xname);  ylabel(yname);
xlim(xax([1 end]));

figufy(h);
% fname = [fig_dir 'pers_ups_vs_pers_downs.pdf'];
% exportfig(fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);

%% HFpow A VS PUP
used_data = l3mec_pups;

n_bins = 20;
X_fac = {cfun_data(:).ctx_down_amps_hf};
Y_fac = {cfun_data(:).ctx_down_isskipped};
xname = 'Ctx down-state HF';
yname = 'Prob down skipped';

clear avg_fun std_fun N_fun sem_fun avg_xvals XY_corr XY_pval
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr(ii),XY_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
end

fprintf('Avg correlation: %.4f, p=%.4g\n',mean(XY_corr),signrank(XY_corr));
nsig = sum(XY_pval < 0.05);
fprintf('%d of %d cells sig\n',nsig,length(XY_pval));

n_nonnan = sum(~isnan(avg_fun));
% h = figure(); hold on
xax = nanmean(avg_xvals);
% xax = pin_cents;
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan),{'color','r'});
plot(xax,nanmean(avg_fun),'ko')
xlabel(xname);  ylabel(yname);
xlim(xax([1 end]));

figufy(h);
% fname = [fig_dir 'pers_ups_vs_pers_downs.pdf'];
% exportfig(fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);

%% MUA VS PUP
used_data = l3mec_pups(used_ctx_mua(l3mec_pups));

n_bins = 20;
X_fac = {cfun_data(:).ctx_down_amps_mua};
Y_fac = {cfun_data(:).ctx_down_isskipped};
xname = 'Ctx down-state MUA';
yname = 'Prob down skipped';

clear avg_fun std_fun N_fun sem_fun avg_xvals XY_corr XY_pval
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr(ii),XY_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
end

fprintf('Avg correlation: %.4f, p=%.4g\n',mean(XY_corr),signrank(XY_corr));
nsig = sum(XY_pval < 0.05);
fprintf('%d of %d cells sig\n',nsig,length(XY_pval));

n_nonnan = sum(~isnan(avg_fun));
h = figure(); hold on
xax = nanmean(avg_xvals);
% xax = pin_cents;
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan));
plot(xax,nanmean(avg_fun),'ko')
xlabel(xname);  ylabel(yname);
xlim(xax([1 end]));

figufy(h);
% fname = [fig_dir 'pers_ups_vs_pers_downs.pdf'];
% exportfig(fname,'width',fig_width,'height',rel_heigh*fig_width,'fontmode','scaled','fontsize',1);

%% MUA BOtH
close all
n_bins = 10;

used_data = l3mec_pdowns(used_ctx_mua(l3mec_pdowns));

X_fac = {cfun_data(:).ctx_up_amps_mua};
Y_fac = {cfun_data(:).ctx_up_isskipped};
clear avg_fun std_fun N_fun sem_fun avg_xvals XY_corr XY_pval
avg_upamp_dist = nan(length(used_data),length(state_amp_xx));
for ii = 1:length(used_data)
    if sum(X_fac{used_data(ii)} == 0) >= length(X_fac{used_data(ii)})/20
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship_withzeros(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    else
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    end
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr(ii),XY_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
    avg_upamp_dist(ii,:) = cfun_data(used_data(ii)).ctx_up_amp_mua_cdist;
end
n_nonnan = sum(~isnan(avg_fun));
h = figure(); 
subplot(2,1,1);hold on
xax = nanmean(avg_xvals);
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan));
plot(xax,nanmean(avg_fun),'ko')

subplot(2,1,2);hold on
n_nonnan = sum(~isnan(avg_upamp_dist));
shadedErrorBar(state_amp_mua_xx,nanmean(avg_upamp_dist),nanstd(avg_upamp_dist)./sqrt(n_nonnan));

used_data = l3mec_pups(used_ctx_mua(l3mec_pups));

X_fac = {cfun_data(:).ctx_down_amps_mua};
Y_fac = {cfun_data(:).ctx_down_isskipped};
clear avg_fun std_fun N_fun sem_fun avg_xvals XY_corr XY_pval
avg_downamp_dist = nan(length(used_data),length(state_amp_xx));
for ii = 1:length(used_data)
    if sum(X_fac{used_data(ii)} == 0) >= length(X_fac{used_data(ii)})/20
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship_withzeros(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    else
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    end
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr(ii),XY_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
    avg_downamp_dist(ii,:) = cfun_data(used_data(ii)).ctx_down_amp_mua_cdist;
end
n_nonnan = sum(~isnan(avg_fun));

subplot(2,1,1);
xax = nanmean(avg_xvals);
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan),{'color','r'});
plot(xax,nanmean(avg_fun),'ko')
xlabel('Ctx state MUA relative rate');
ylabel('Prob. skipped');

subplot(2,1,2);
n_nonnan = sum(~isnan(avg_downamp_dist));
shadedErrorBar(state_amp_mua_xx,nanmean(avg_downamp_dist),nanstd(avg_downamp_dist)./sqrt(n_nonnan),{'color','r'});
xlabel('Ctx state MUA relative rate');
ylabel('Cumulative Probability');

% subplot(2,1,1);
% xlim([-0.01 4]);
% subplot(2,1,2);
% xlim([-0.01 4]);

subplot(2,1,1);
xlim([0.03 4]); set(gca,'xscale','log');
subplot(2,1,2);
xlim([0.03 4]); set(gca,'xscale','log');

%%
fig_width = 4; rel_height = 1.6;
figufy(h);
fname = [fig_dir 'prob_skip_vs_muaamp_log.pdf'];
exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%% State duration BOTH
close all
n_bins = 15;
state_dur_xx = linspace(0.5,5,200);

used_data = l3mec_pdowns;

X_fac = {core_data(:).lfp_up_durs};
Y_fac = {cfun_data(:).ctx_up_isskipped};

clear avg_fun std_fun N_fun sem_fun avg_xvals US_corr US_pval avg_updur_dist
for ii = 1:length(used_data)
    if length(X_fac{used_data(ii)}) ~= length(Y_fac{used_data(ii)})
        error('!');
    end
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [US_corr(ii),US_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
    avg_updur_dist(ii,:) = ksdensity(core_data(used_data(ii)).lfp_up_durs,state_dur_xx,'support','positive');
end

n_nonnan = sum(~isnan(avg_fun));
h = figure(); 
subplot(2,1,1);
hold on
xax = nanmean(avg_xvals);
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan));
plot(xax,nanmean(avg_fun),'ko')

subplot(2,1,2); hold on
n_nonnan = sum(~isnan(avg_updur_dist));
shadedErrorBar(state_dur_xx,nanmean(avg_updur_dist),nanstd(avg_updur_dist)./sqrt(n_nonnan));

used_data = l3mec_pups;
X_fac = {core_data(:).lfp_down_durs};
Y_fac = {cfun_data(:).ctx_down_isskipped};

clear avg_fun std_fun N_fun sem_fun avg_xvals DS_corr DS_pval avg_downdur_dist
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [DS_corr(ii),DS_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
    avg_downdur_dist(ii,:) = ksdensity(core_data(used_data(ii)).lfp_down_durs,state_dur_xx,'support','positive');
end

subplot(2,1,1);
n_nonnan = sum(~isnan(avg_fun));
xax = nanmean(avg_xvals);
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan),{'color','r'});
plot(xax,nanmean(avg_fun),'ko')
xlabel('State duration (s)');
ylabel('Prob. skipped');

subplot(2,1,2);
n_nonnan = sum(~isnan(avg_downamp_dist));
shadedErrorBar(state_dur_xx,nanmean(avg_downdur_dist),nanstd(avg_downdur_dist)./sqrt(n_nonnan),{'color','r'});
xlabel('State duration (s)');
ylabel('Probability');

subplot(2,1,1);
xlim([0.5 4]); 
set(gca,'xscale','log');

subplot(2,1,2);
xlim([0.5 4]); 
set(gca,'xscale','log');

%%
fig_width = 4; rel_height = 1.6;
figufy(h);
fname = [fig_dir 'prob_skip_vs_state_dur.pdf'];
exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);



%% HF pow vs up-trans lag
used_data = l3mec;

n_bins = 20;
% X_fac = {cfun_data(:).ctx_up_amps_hf};
% Y_fac = {cfun_data(:).corresp_mp_uplag};
% X_fac = {cfun_data(:).ctx_up_amps_hf};
X_fac = {cfun_data(:).ctx_down_amps_hf};
Y_fac = {cfun_data(:).corresp_mp_downlag};
xname = 'Ctx up-state HF amp';
yname = 'UP-transition lag (s)';

clear avg_fun std_fun N_fun sem_fun avg_xvals XY_corr XY_pval
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [XY_corr(ii),XY_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
end

fprintf('Avg correlation: %.4f, p=%.4g\n',mean(XY_corr),signrank(XY_corr));
nsig = sum(XY_pval < 0.05);
fprintf('%d of %d cells sig\n',nsig,length(XY_pval));

n_nonnan = sum(~isnan(avg_fun));
h = figure(); hold on
xax = nanmean(avg_xvals);
% xax = pin_cents;
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan));
plot(xax,nanmean(avg_fun),'ko')
xlabel(xname);  ylabel(yname);
xlim(xax([1 end]));

%%
fig_width = 4; rel_height = 0.8;
figufy(h);
fname = [fig_dir 'up_lag_vs_ctx_upamp_hf.pdf'];
exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
%% State duration BOTH
close all
n_bins = 15;
state_dur_xx = linspace(0.5,5,200);

used_data = l3lec;

X_fac = {cfun_data(:).ctx_up_amps_hf};
Y_fac = {cfun_data(:).corresp_mp_uplag};
Z_fac = {cfun_data(:).ctx_up_isskipped};
for ii = 1:length(used_data)
    cur = Z_fac{used_data(ii)};
    uu = find(cur ~= 0);
    Y_fac{used_data(ii)}(uu) = nan;
    X_fac{used_data(ii)}(uu) = nan;
end
clear avg_fun std_fun N_fun sem_fun avg_xvals US_corr US_pval avg_updur_dist
for ii = 1:length(used_data)
    if length(X_fac{used_data(ii)}) ~= length(Y_fac{used_data(ii)})
        error('!');
    end
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [US_corr(ii),US_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
    avg_upamp_dist(ii,:) = cfun_data(used_data(ii)).ctx_up_amp_hf_dist;
end

n_nonnan = sum(~isnan(avg_fun));
h = figure(); 
subplot(2,1,1);
hold on
xax = nanmean(avg_xvals);
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan));
plot(xax,nanmean(avg_fun),'ko')

subplot(2,1,2); hold on
n_nonnan = sum(~isnan(avg_upamp_dist));
shadedErrorBar(state_amp_xx,nanmean(avg_upamp_dist),nanstd(avg_upamp_dist)./sqrt(n_nonnan));

used_data = l3lec;
X_fac = {cfun_data(:).ctx_down_amps_hf};
Y_fac = {cfun_data(:).corresp_mp_downlag};
Z_fac = {cfun_data(:).ctx_down_isskipped};
for ii = 1:length(used_data)
    cur = Z_fac{used_data(ii)};
    uu = find(cur ~= 0);
    Y_fac{used_data(ii)}(uu) = nan;
    X_fac{used_data(ii)}(uu) = nan;
end

clear avg_fun std_fun N_fun sem_fun avg_xvals DS_corr DS_pval avg_downdur_dist
for ii = 1:length(used_data)
    [avg_fun(ii,:),std_fun(ii,:),N_fun(ii,:),pin_cents,avg_xvals(ii,:)] = ...
        get_nonparametric_relationship(Y_fac{used_data(ii)},X_fac{used_data(ii)},n_bins);
    sem_fun(ii,:) = std_fun(ii,:)./sqrt(N_fun(ii,:));
    [DS_corr(ii),DS_pval(ii)] = nancorr(X_fac{used_data(ii)},Y_fac{used_data(ii)},'spearman');
    avg_downamp_dist(ii,:) = cfun_data(used_data(ii)).ctx_down_amp_hf_dist;
end

subplot(2,1,1);
n_nonnan = sum(~isnan(avg_fun));
xax = nanmean(avg_xvals);
shadedErrorBar(xax,nanmean(avg_fun),nanstd(avg_fun)./sqrt(n_nonnan),{'color','r'});
plot(xax,nanmean(avg_fun),'ko')
xlabel('State amplitude (z)');
ylabel('Transition lag (s)');

subplot(2,1,2);
n_nonnan = sum(~isnan(avg_downamp_dist));
shadedErrorBar(state_amp_xx,nanmean(avg_downamp_dist),nanstd(avg_downamp_dist)./sqrt(n_nonnan),{'color','r'});
xlabel('State amplitude (z)');
ylabel('Probability');

subplot(2,1,1);
xlim([-1 2]); 
% % set(gca,'xscale','log');
% 
subplot(2,1,2);
xlim([-1 2]); 
% % set(gca,'xscale','log');

%%
fig_width = 4; rel_height = 1.6;
figufy(h);
fname = [fig_dir 'trans_lag_vs_hfamp_lec.pdf'];
exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
