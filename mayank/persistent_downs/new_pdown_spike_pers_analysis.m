clear all
data_dir = cell(0);
data_type = cell(0);
data_exp = cell(0);
data_ep = [];
data_dp = [];
data_hpc_lfp = [];
data_hpc_mua = [];

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
%     data_hpc_mua = cat(2,data_hpc_mua,hpc_mua(l3mec(ii)));
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
%     data_hpc_mua = cat(2,data_hpc_mua,hpc_mua(l3mec_np(ii)));
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
%     data_hpc_mua = cat(2,data_hpc_mua,hpc_mua(l3lec(ii)));
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
%     data_hpc_mua = cat(2,data_hpc_mua,hpc_mua(l3lec_np(ii)));
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
%     data_hpc_mua = cat(2,data_hpc_mua,new_pdown_hpcmua(cur_uset(ii))-1);
end

addpath('C:\WC_Germany\parietal_cortical_2010\');
addpath('C:\WC_Germany\persistent_9_27_2010\');
addpath('C:\WC_Germany\new_mec\');
addpath('C:\WC_Germany\overall_EC');
addpath('C:\WC_Germany\hsmm_state_detection');
addpath('C:\WC_Germany\sven_thomas_combined');
addpath('C:\Code\general_functions\');

min_rec_dur = 500; %in sec
used_dirs = find(data_ep > min_rec_dur);

%%
load ./mua_classification
min_hpcrate = 5;
usable_hpc_mua = ~isnan(peak_hpcmua_loc) & peak_hpcmua_rate >= min_hpcrate;

%%
dsf = 8;
raw_Fs = 2016;
Fsd = raw_Fs/dsf;
rate_sm = round(Fsd*0.05);

cd C:/WC_Germany/persistent_downs/
load ./new_down_core_analysis
%%
for d = 1:length(used_dirs)
    cd(data_dir{used_dirs(d)})
    disp(sprintf('session %d',d))
    
    load ./spike_time_jmm
    load ./used_data wcv
    load ./sync_times.mat
    [wcv_lf,t_axis] = get_lf_features(wcv,raw_Fs,Fsd,[lcf hcf]);
    rec_dur(d) = length(wcv)/raw_Fs;
    synct_d = downsample(synct,dsf);
    end_time = min(data_ep(used_dirs(d)),data_dp(used_dirs(d)));
    ep = find(t_axis >= end_time,1);
    if ~isempty(ep)
        lf8_lf(ep+1:end) = []; wcv_lf(ep+1:end) = []; t_axis(ep+1:end) = []; synct_d(ep+1:end) = [];
    else
        ep = length(t_axis);
    end
    
    %% LOAD IN MUA
    if exist('./mua_data3.mat','file')
        load ./mua_data3
        %         if ~isnan(data_hpc_mua(used_dirs(d)))
        if usable_hpc_mua(used_dirs(d))
%             hpc_mua_times = mua_times{data_hpc_mua(used_dirs(d))};
             hpc_mua_times = mua_times{peak_hpcmua_loc(used_dirs(d))};
           hpc_mua_times(hpc_mua_times > synct_d(end) | hpc_mua_times < synct_d(1)) = [];
            hpc_mua_rate =hist(hpc_mua_times,synct_d)*Fsd;
            hpc_mua_rate = jmm_smooth_1d_cor(hpc_mua_rate,rate_sm);
        else
            hpc_mua_rate = nan(size(synct_d));
        end
        usable_ctx_mua = find(usable_mua(used_dirs(d),6:end));
        used_ctx_mua_chs{d} = usable_ctx_mua;
        ctx_mua_times = [];
        for ii = 1:length(usable_ctx_mua)
            ctx_mua_times = [ctx_mua_times mua_times{usable_ctx_mua(ii)+5}];
        end
        ctx_mua_times = sort(ctx_mua_times);
%         ctx_mua_times = sort([mua_times{7} mua_times{8}]);
        ctx_mua_rate = hist(ctx_mua_times,synct_d)*Fsd;
        ctx_mua_rate = jmm_smooth_1d_cor(ctx_mua_rate,rate_sm);
        
        if isempty(usable_ctx_mua)
            ctx_mua_rate = nan(size(synct_d));
        end
        
        if length(hpc_mua_rate) > length(synct_d)
            hpc_mua_rate = hpc_mua_rate(1:length(synct_d));
            ctx_mua_rate = ctx_mua_rate(1:length(synct_d));
        end
    else
        hpc_mua_rate = nan(size(synct_d));
        ctx_mua_rate = nan(size(synct_d));
        usable_ctx_mua = [];
    end
    avg_hpc_mua_rate(d) = nanmean(hpc_mua_rate);
    avg_ctx_mua_rate(d) = nanmean(ctx_mua_rate);
    std_hpc_mua_rate(d) = nanstd(hpc_mua_rate);
    std_ctx_mua_rate(d) = nanstd(ctx_mua_rate);
%     hpc_mua_rate = zscore(hpc_mua_rate);
%     ctx_mua_rate = zscore(ctx_mua_rate);

    mp_spike_rate = hist(synct(spkid),synct_d)*Fsd;
    
    backgnd_rate(d) = length(spkid)/length(synct_d)*Fsd;

    %%
    
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
    
    [new_seg_inds] = resample_uds_seg_inds_v2(hmm.UDS_segs,hmm.Fs,Fsd,mp_state_seq);
    dur_uds = sum(diff(new_seg_inds,[],2))/Fsd;
    [up_trans_inds,down_trans_inds] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
    [up_trans_inds8,down_trans_inds8] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);
    
    bad_mp_states = find(up_trans_inds > ep | down_trans_inds > ep);
    up_trans_inds(bad_mp_states) = []; down_trans_inds(bad_mp_states) = [];
    bad_lfp_states = find(up_trans_inds8 > ep | down_trans_inds8 > ep);
    up_trans_inds8(bad_lfp_states) = []; down_trans_inds8(bad_lfp_states) = [];
    
    %%
    wcv_up_log = nan(size(synct_d));
    lf8_up_log = nan(size(synct_d));
    
    for ns = 1:hmm.Nsegs
        cur_seg = new_seg_inds(ns,1):new_seg_inds(ns,2);
        wcv_up_log(cur_seg) = logical(mp_state_seq{ns}-1);
        lf8_up_log(cur_seg) = logical(lfp_state_seq{ns}-1);
    end
    wcv_up_log(length(mp_spike_rate)+1:end) = [];
    lf8_up_log(length(mp_spike_rate)+1:end) = [];
    
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
       
     pers_lfp_ups = [];
     for i = 1:length(rt2_downs{d})
         pers_lfp_ups = [pers_lfp_ups mp_downskipped{d}.inds{rt2_downs{d}(i)}];
     end
     npers_lfp_ups = setdiff(1:length(up_trans_inds8),pers_lfp_ups);
     
    temp_dpers_hpc_rate = nan(size(pers_lfp_ups));
    temp_ndpers_hpc_rate = nan(size(npers_lfp_ups));
    for j = 1:length(pers_lfp_ups)
        lfp_up_inds = up_trans_inds8(pers_lfp_ups(j)):down_trans_inds8(pers_lfp_ups(j));
        temp_dpers_hpc_rate(j) = mean(hpc_mua_rate(lfp_up_inds)); 
    end
    for j = 1:length(npers_lfp_ups)
        lfp_up_inds = up_trans_inds8(npers_lfp_ups(j)):down_trans_inds8(npers_lfp_ups(j));
        temp_ndpers_hpc_rate(j) = mean(hpc_mua_rate(lfp_up_inds));        
    end
    if ~isempty(pers_lfp_ups)
        hpc_dpers_rate(d) = mean(temp_dpers_hpc_rate);
    else
        hpc_dpers_rate(d) = nan;
    end
    hpc_ndpers_rate(d) = mean(temp_ndpers_hpc_rate);

     
     mp_before_during_diff(d) = nanmean(temp_pers_during_mp-temp_pers_around_mp);
     hpc_before_during_diff(d) = nanmean(temp_pers_during_hpc-temp_pers_around_hpc);
     
    mp_up_rate(d) = mean(mp_spike_rate(wcv_up_log==1));
    mp_down_rate(d) = mean(mp_spike_rate(wcv_up_log==0));
    mp_up_rate_lfup(d) = mean(mp_spike_rate(wcv_up_log==1 & lf8_up_log==1));
    mp_up_rate_lfdown(d) = mean(mp_spike_rate(wcv_up_log==1 & lf8_up_log==0));
    mp_rate(d) = mean(mp_spike_rate(wcv_up_log==1 | wcv_up_log==0));
    
    dskipped_lfp_trans = [];
    for j = 1:length(up_trans_inds)
        dskipped_lfp_trans = [dskipped_lfp_trans mp_downskipped{d}.inds{j}];
    end
    non_dskipped_lfp_trans = setdiff(1:length(up_trans_inds8),dskipped_lfp_trans);

    cur_delay = round(median_uplag(d));
    hpc_ctx_uprates_lag = nan(length(up_trans_inds8),1);
    for i = 1:length(up_trans_inds8)
        cur_inds = (up_trans_inds8(i):down_trans_inds8(i)) + cur_delay;
        cur_inds(cur_inds > length(hpc_mua_rate)) = [];
        hpc_ctx_uprates_lag(i) = nanmean(hpc_mua_rate(cur_inds));
    end
    hpc_rates_pdown(d) = nanmean(hpc_ctx_uprates_lag(dskipped_lfp_trans));
    hpc_rates_npdown(d) = nanmean(hpc_ctx_uprates_lag(non_dskipped_lfp_trans));
    
    if ~isempty(usable_ctx_mua)
         ctx_rate_lfup(d) = mean(ctx_mua_rate(lf8_up_log==1));
        ctx_rate_lfdown(d) = mean(ctx_mua_rate(lf8_up_log==0));
        ctx_up_rate(d) = mean(ctx_mua_rate(wcv_up_log==1));
        ctx_up_rate_lfup(d) = mean(ctx_mua_rate(wcv_up_log==1 & lf8_up_log==1));
        ctx_up_rate_lfdown(d) = mean(ctx_mua_rate(wcv_up_log==1 & lf8_up_log==0));
        ctx_down_rate(d) = mean(ctx_mua_rate(wcv_up_log==0));
        ctx_down_rate_lfup(d) = mean(ctx_mua_rate(wcv_up_log==0 & lf8_up_log==1));
        ctx_down_rate_lfdown(d) = mean(ctx_mua_rate(wcv_up_log==0 & lf8_up_log==0));
        ctx_rate(d) = mean(ctx_mua_rate(wcv_up_log==1 | wcv_up_log==0));       
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
    end
    
%     if ~isnan(data_hpc_mua(d))
        if usable_hpc_mua(used_dirs(d))
    
        hpc_rate_lfup(d) = mean(hpc_mua_rate(lf8_up_log==1));
        hpc_rate_lfdown(d) = mean(hpc_mua_rate(lf8_up_log==0));
        hpc_up_rate(d) = mean(hpc_mua_rate(wcv_up_log==1));
        hpc_up_rate_lfup(d) = mean(hpc_mua_rate(wcv_up_log==1 & lf8_up_log==1));
        hpc_up_rate_lfdown(d) = mean(hpc_mua_rate(wcv_up_log==1 & lf8_up_log==0));
        hpc_down_rate(d) = mean(hpc_mua_rate(wcv_up_log==0));
        hpc_down_rate_lfup(d) = mean(hpc_mua_rate(wcv_up_log==0 & lf8_up_log==1));
        hpc_down_rate_lfdown(d) = mean(hpc_mua_rate(wcv_up_log==0 & lf8_up_log==0));
        hpc_rate(d) = mean(hpc_mua_rate(wcv_up_log==1 | wcv_up_log==0));
        
    else
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

cd C:\WC_Germany\persistent_downs\
save spike_rate_data_fin_nd_np2 mp_* hpc_* ctx_* rec_dur backgnd_rate avg_* std_*
% l3mec_m = l3mec(~isnan(hpc_mua(l3mec)));

%%
cd C:\persDowns_paper\Figs\

hpc_rates_npdown_z = (hpc_rates_npdown - avg_hpc_mua_rate)./std_hpc_mua_rate;
hpc_rates_pdown_z = (hpc_rates_pdown - avg_hpc_mua_rate)./std_hpc_mua_rate;

figure
plot(hpc_rates_pdown_z,hpc_rates_npdown_z,'.','markersize',8);
box off;
xlabel('Skipped ctx Up states','fontsize',12);
ylabel('Non-skipped ctx Up states','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
fillPage(gcf,'papersize',[5 5]);
line([-1 0.8],[-1 0.8],'color','k','linestyle','--');

%%
hpc_up_zrate_lfup = bsxfun(@minus,hpc_up_rate_lfup,avg_hpc_mua_rate);
hpc_up_zrate_lfup = bsxfun(@rdivide,hpc_up_zrate_lfup,std_hpc_mua_rate);
hpc_up_zrate_lfdown = bsxfun(@minus,hpc_up_rate_lfdown,avg_hpc_mua_rate);
hpc_up_zrate_lfdown = bsxfun(@rdivide,hpc_up_zrate_lfdown,std_hpc_mua_rate);
hpc_down_zrate_lfup = bsxfun(@minus,hpc_down_rate_lfup,avg_hpc_mua_rate);
hpc_down_zrate_lfup = bsxfun(@rdivide,hpc_down_zrate_lfup,std_hpc_mua_rate);
hpc_down_zrate_lfdown = bsxfun(@minus,hpc_down_rate_lfdown,avg_hpc_mua_rate);
hpc_down_zrate_lfdown = bsxfun(@rdivide,hpc_down_zrate_lfdown,std_hpc_mua_rate);

figure
Y = [hpc_up_zrate_lfup'; hpc_up_zrate_lfdown'; hpc_down_zrate_lfdown'; hpc_down_zrate_lfup'];
G = [ones(length(used_dirs),1); 2*ones(length(used_dirs),1); 3*ones(length(used_dirs),1); 4*ones(length(used_dirs),1)];
boxplot(Y,G)
set(gca,'xtick',1:4, 'xticklabel',{'MecUP/CtxUp','MecUP/CtxDOWN','MecDOWN/CtxDOWN','MecDOWN/CtxUP'})
ylabel('Average hpc MUA rate (z)','fontsize',16)
set(gca,'fontsize',8,'fontname','arial')
box off
fillPage(gcf,'papersize',[6 5]);

%%
figure
set(gca,'fontsize',14,'fontname','arial')
plot(mp_up_rate_lfup,mp_up_rate_lfdown,'r.','markersize',14)
% hold on
% plot(mp_up_rate_lfup(l3lec),mp_up_rate_lfdown(l3lec),'b.','markersize',14)
line([0 12],[0 12],'Color','k')
xlabel('MEC firing rate during neocortical UP state (Hz)','fontsize',16,'fontname','arial')
ylabel('MEC firing rate during neocortical DOWN state (Hz)','fontsize',16,'fontname','arial')


%%
% figure
% plot(ctx_up_rate_lfup(l3mec),ctx_up_rate_lfdown(l3mec),'c.','markersize',14); 
% hold on
% plot(ctx_down_rate_lfup(l3mec),ctx_down_rate_lfdown(l3mec),'g.','markersize',14); 
% xl = xlim();yl = ylim();
% line(xl,yl,'color','k')
% l3mec_m = 1:length(used_dirs);
% l3mec_m = 1:length(used_dirs);
% l3mec_m = 1:67;
% l3mec_new = 68:length(used_dirs);

mp_hpc_rate_mod = (hpc_up_rate - hpc_down_rate)./(hpc_up_rate + hpc_down_rate);
mp_hpc_rate_mod_lfup = (hpc_up_rate_lfup - hpc_down_rate_lfup)./(hpc_up_rate_lfup + hpc_down_rate_lfup);
mp_hpc_rate_mod_lfdown = (hpc_up_rate_lfdown - hpc_down_rate_lfdown)./(hpc_up_rate_lfdown + hpc_down_rate_lfdown);
ctx_hpc_rate_mod = (hpc_rate_lfup - hpc_rate_lfdown)./(hpc_rate_lfup + hpc_rate_lfdown);
ctx_hpc_rate_mod_mpup = (hpc_up_rate_lfup - hpc_up_rate_lfdown)./(hpc_up_rate_lfup + hpc_up_rate_lfdown);
ctx_hpc_rate_mod_mpdown = (hpc_down_rate_lfup - hpc_down_rate_lfdown)./(hpc_down_rate_lfup + hpc_down_rate_lfdown);

figure; hold on
plot(mp_hpc_rate_mod,ctx_hpc_rate_mod,'k.','markersize',14)
xl = xlim();yl = ylim();
line([-0.3 0.7],[-0.3 0.7],'color','k')
xlabel('MP State Depth of modulation','fontsize',12)
ylabel('Ctx State Depth of modulation','fontsize',12)
xlim([-0.3 0.7]); ylim([-0.3 0.7]);
box off;
set(gca,'fontsize',10,'fontname','arial');

% figure; hold on
% plot(mp_hpc_rate_mod_lfup,mp_hpc_rate_mod_lfdown,'k.','markersize',14)
% plot(ctx_hpc_rate_mod_mpup,ctx_hpc_rate_mod_mpdown,'r.','markersize',14)
% line([-0.2 0.7],[-0.2 0.7],'color','k')
% 
%%
mp_ctx_rate_mod = (ctx_up_rate - ctx_down_rate)./(ctx_up_rate + ctx_down_rate);
ctx_ctx_rate_mod = (ctx_rate_lfup - ctx_rate_lfdown)./(ctx_rate_lfup + ctx_rate_lfdown);

figure; hold on
box off;
plot(mp_ctx_rate_mod,ctx_ctx_rate_mod,'k.','markersize',14)
xl = xlim();yl = ylim();
line([-0.2 1],[-0.2 1],'color','k')
xlabel('MP State Depth of modulation','fontsize',12)
ylabel('Ctx State Depth of modulation','fontsize',12)
xlim([0 1]); ylim([0 1]);
box off;
set(gca,'fontsize',10,'fontname','arial');

