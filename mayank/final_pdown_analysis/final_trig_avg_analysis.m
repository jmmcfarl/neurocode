%%
clear all
load C:/WC_Germany/final_pdown_analysis/compiled_data.mat
load C:/WC_Germany/final_pdown_analysis/mua_classification2.mat

%include channels where hpc MUA peak appears on ch5
peak_hpcmua_loc([27 29 67 68 79]) = 5;

fig_dir = 'C:\WC_Germany\final_pdown_analysis\figures\';

min_rec_dur = 500; %in sec
data_ids = [data(:).id];
used_dirs = find([data(:).ep] > min_rec_dur & [data(:).dp] > min_rec_dur);
used_dirs(ismember(data_ids(used_dirs),unclear_uds)) = [];

data = data(used_dirs);

peak_hpcmua_loc = peak_hpcmua_loc(used_dirs);
peak_hpcmua_rate = peak_hpcmua_rate(used_dirs);
usable_mua = usable_mua(used_dirs,:);
avg_spkwidths = avg_spkwidths(used_dirs,:);
avg_spkwidth_fwhm = avg_spkwidth_fwhm(used_dirs,:);
avg_spkwidth_fwqm = avg_spkwidth_fwqm(used_dirs,:);

load C:/WC_Germany/final_pdown_analysis/fin_pdown_core_analysis.mat
if length(core_data) ~= length(data)
    error('Data mismatch');
end

min_hpcrate = 1;
usable_hpc_mua = ~isnan(peak_hpcmua_loc) & peak_hpcmua_rate >= min_hpcrate;

hpc_spkwidths = nan(length(usable_hpc_mua),3);
hpc_avgrate = nan(length(usable_hpc_mua),1);
for ii = 1:length(usable_hpc_mua)
    if usable_hpc_mua(ii)
        hpc_avgrate(ii) = avg_rates(ii,peak_hpcmua_loc(ii));
        hpc_spkwidths(ii,:) = [avg_spkwidths(ii,peak_hpcmua_loc(ii)) ...
            avg_spkwidth_fwhm(ii,peak_hpcmua_loc(ii)) avg_spkwidth_fwqm(ii,peak_hpcmua_loc(ii))];
    end
end

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

backlag = 2*Fsd;
forwardlag = 3*Fsd;

min_samps = 5;
%%
for d = 1:length(used_dirs)
    cd(data(d).dir)
    pwd
    
    load ./used_data lf7 lf3 lf5 wcv_minus_spike
%     if data(d).hpc_lfp == 3
        load ./used_data lf3 lf4 lf5
        if data(d).is_old_type
            lf3 = lf3 + lf5; %redefine LF3 wrt gnd
        end
        hpc_3 = lf3;
%     else
        load ./used_data lf2
        hpc_2 = lf2;
%     end

    if peak_hpcmua_loc(d) == 3
        hpc_3 = lf3;
        hpc_2 = lf2;
    elseif peak_hpcmua_loc(d) == 4
        hpc_3 = lf4;
        hpc_2 = lf3;
    end
    
    [lfp_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf]);
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
    hpc_lf3 = get_lf_features(hpc_3,raw_Fs,Fsd,[lcf hcf]);
    hpc_lf2 = get_lf_features(hpc_2,raw_Fs,Fsd,[lcf hcf]);
    
%     lfp_hf = get_hf_features(lf7,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
%     wcv_hf = get_hf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
%     hpc_hf = get_hf_features(hpc_lfp,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm);
    lfp_hf = zscore(get_hf_lpower(lf7,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm));
    wcv_hf = zscore(get_hf_lpower(wcv_minus_spike,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm));
    hpc_hf3 = zscore(get_hf_lpower(hpc_3,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm));
    hpc_hf2 = zscore(get_hf_lpower(hpc_2,raw_Fs,Fsd,[lcf_hf hcf_hf],hcf_sm));
    
    
    load ./sync_times.mat
    synct_d = downsample(synct,dsf);
    
    load ./spike_time_jmm

    end_time = min(data(d).ep,data(d).dp);
    ep = find(t_axis >= end_time,1);
    if ~isempty(ep)
        synct_d(ep+1:end) = []; lfp_lf(ep+1:end) = []; wcv_lf(ep+1:end) = []; t_axis(ep+1:end) = []; lfp_hf(ep+1:end) = [];
        hpc_lf2(ep+1:end) = []; hpc_lf3(ep+1:end) = []; 
        hpc_hf2(ep+1:end) = []; hpc_hf3(ep+1:end) = []; wcv_hf(ep+1:end) = [];
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
           
    %%
            %% CREATE INTERPOLATED MP DATA
        if exist('./aligned_heka.mat','file')
            load ./spike_time_jmm.mat
            full_t_axis = (1:length(wcv_minus_spike))/raw_Fs;
            spike_times = full_t_axis(spkid);

            load ./aligned_heka.mat
            dc_dt = median(diff(dc_time));
            dc_fs = 1/dc_dt;
            dc_dsf = 20;
            dc_fsd = dc_fs/dc_dsf;
            [dc_fb,dc_fa] = butter(4,20/(dc_fsd/2),'low');
%             [dc_fb,dc_fa] = butter(2,[0.01 10]/(dc_fsd/2));
            
            %remove spikes from DC MP signal by interpolation
            spk_back = 0.001;
            spk_for = 0.004;
            cur_spk_inds = round(interp1(dc_time,1:length(dc_time),spike_times));
            spike_times(isnan(cur_spk_inds)) = [];
            cur_spk_inds(isnan(cur_spk_inds)) = [];
            cur_spk_inds = cur_spk_inds(:); spike_times = spike_times(:); dc_time = dc_time(:);
            spike_error = abs(dc_time(cur_spk_inds) - spike_times); %use only spikes that interpolate onto actual DC data
            bad = find(spike_error > 2/Fsd);
            cur_spk_inds(bad) = [];
            blk_inds = -round(spk_back/dc_dt):round(spk_for/dc_dt); %set of times around each spike to interpolate out
            blocked_inds = bsxfun(@plus,cur_spk_inds,blk_inds);
            blocked_inds = blocked_inds(:);
            used_inds = setdiff(1:length(dc_data),blocked_inds); %non-spike samples
            dc_data = interp1(used_inds,dc_data(used_inds),1:length(dc_data)); %de-spiked data
            
            dc_data = decimate(dc_data,dc_dsf);
            dc_data = filtfilt(dc_fb,dc_fa,dc_data);
            dc_time = downsample(dc_time,dc_dsf);
            dc_time = dc_time(:);
            
            dc_interp_data = interp1(dc_time,dc_data,t_axis);
            ind_interp = ceil(interp1(dc_time,1:length(dc_time),t_axis));
            used = find(~isnan(ind_interp));
            t_error = dc_time(ind_interp(used)) - t_axis(used)';
            used = used(t_error <= 1/Fsd);
            
            dc_interp_data(setdiff(1:length(dc_interp_data),used)) = nan;
                        
            if nanstd(dc_interp_data) < 1
                dc_interp_data = dc_interp_data*100;
            end
            dc_interp_data = dc_interp_data(:);
        else
            dc_interp_data = nan(size(wcv_lf));
        end
        
        %COMPUTE AVG UP- AND DOWN-STATE CONDITIONAL DC MP VALUES FOR
        %NORMALIZATION
        trig_data(d).dc_upstate_amp = nanmean(dc_interp_data(mp_state_vec == 1));
        trig_data(d).dc_downstate_amp = nanmean(dc_interp_data(mp_state_vec == 0));
        trig_data(d).dc_avg_amp = nanmean(dc_interp_data(~isnan(mp_state_vec)));
        trig_data(d).dc_uds_amp = trig_data(d).dc_upstate_amp - trig_data(d).dc_downstate_amp;
        
        dc_interp_data_reld = dc_interp_data - trig_data(d).dc_downstate_amp;
        dc_interp_data_normd = dc_interp_data_reld/trig_data(d).dc_uds_amp;

%%
    if exist('./mua_data3.mat','file')
        load ./mua_data3
        if usable_hpc_mua(d)
            hpc_mua_times = mua_times{peak_hpcmua_loc(d)};
            hpc_mua_times(hpc_mua_times > synct_d(end) | hpc_mua_times < synct_d(1)) = [];
            hpc_mua_rate =hist(hpc_mua_times,synct_d)*Fsd;
%             hpc_mua_rate = sqrt(hpc_mua_rate);
            hpc_mua_rate = jmm_smooth_1d_cor(hpc_mua_rate,rate_sm);
        else
            hpc_mua_rate = nan(size(synct_d));
        end
        usable_ctx_mua = find(usable_mua(d,6:end));
        used_ctx_mua_chs{d} = usable_ctx_mua;
        ctx_mua_times = [];
        for ii = 1:length(usable_ctx_mua)
            ctx_mua_times = [ctx_mua_times mua_times{usable_ctx_mua(ii)+5}];
        end
        ctx_mua_times = sort(ctx_mua_times);
        ctx_mua_times(ctx_mua_times > synct_d(end) | ctx_mua_times < synct_d(1)) = [];
        ctx_mua_rate = hist(ctx_mua_times,synct_d)*Fsd;
%         ctx_mua_rate = sqrt(ctx_mua_rate);
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
    end
    
%     mp_spike_rate = sqrt(mp_spike_rate);
    mp_spike_rate = jmm_smooth_1d_cor(mp_spike_rate,rate_sm);
    
    trig_data(d).avg_hpc_mua_rate = nanmean(hpc_mua_rate);
    trig_data(d).avg_ctx_mua_rate = nanmean(ctx_mua_rate);
    trig_data(d).avg_mp_rate = nanmean(mp_spike_rate);
    trig_data(d).std_hpc_mua_rate = nanstd(hpc_mua_rate);
    trig_data(d).std_ctx_mua_rate(d) = nanstd(ctx_mua_rate);
    trig_data(d).std_mp_rate = nanstd(mp_spike_rate);
    
%     hpc_mua_rate = sqrt(hpc_mua_rate);
%     ctx_mua_rate = sqrt(ctx_mua_rate);
%     mp_spike_rate = sqrt(mp_spike_rate);   
    
    hpc_mua_rate = zscore(hpc_mua_rate);
    ctx_mua_rate = zscore(ctx_mua_rate);
    mp_spike_rate = zscore(mp_spike_rate);

    %%
    ctx_ups_mp_up = find(mp_state_vec(up_trans_inds8) == 1);
    skipped_lfp_ups = [core_data(d).skipped_lfp_ups];
    non_skipped_lfp_ups = [core_data(d).non_skipped_lfp_ups];
    
    ctx_downs_mp_down = find(mp_state_vec(down_trans_inds8) == 0);
    skipped_lfp_downs = [core_data(d).skipped_lfp_downs];
    non_skipped_lfp_downs = [core_data(d).non_skipped_lfp_downs];
          
    temp1 = find(ismember(skipped_lfp_ups,ctx_ups_mp_up));
    temp2 = find(ismember(skipped_lfp_downs,ctx_downs_mp_down));
    if ~isempty(temp1) | ~isempty(temp2)
        pause
    end
    
    trig_data(d).N_non_skipped_ups_MPup = sum(ismember(non_skipped_lfp_ups,ctx_ups_mp_up));
    trig_data(d).N_skipped_ups_MPup = sum(ismember(skipped_lfp_ups,ctx_ups_mp_up));
    skipped_lfp_ups(ismember(skipped_lfp_ups,ctx_ups_mp_up)) = [];
    non_skipped_lfp_ups(ismember(non_skipped_lfp_ups,ctx_ups_mp_up)) = [];

    trig_data(d).N_non_skipped_downs_MPdown = sum(ismember(non_skipped_lfp_downs,ctx_downs_mp_down));
    trig_data(d).N_skipped_downs_MPdown = sum(ismember(skipped_lfp_downs,ctx_downs_mp_down));
    skipped_lfp_downs(ismember(skipped_lfp_downs,ctx_downs_mp_down)) = [];
    non_skipped_lfp_downs(ismember(non_skipped_lfp_downs,ctx_downs_mp_down)) = [];

    lfp_skipped_up = false(size(lfp_state_number));
    lfp_nskipped_up = lfp_skipped_up;
    lfp_skipped_down = lfp_skipped_up;
    lfp_nskipped_down = lfp_skipped_up;
    lfp_prevskipped_up = false(size(lfp_state_number));
    lfp_nextskipped_up = false(size(lfp_state_number));
    
%     buffer_win = round(0.1*Fsd);
    buffer_win = round(0*Fsd);
    for i = 1:length(skipped_lfp_ups)
%         cur_inds = (up_trans_inds8(skipped_lfp_ups(i))+buffer_win):(down_trans_inds8(skipped_lfp_ups(i))-buffer_win);
        cur_inds = (up_trans_inds8(skipped_lfp_ups(i))+buffer_win):(down_trans_inds8(skipped_lfp_ups(i)));
        lfp_skipped_up(cur_inds) = true;
        if skipped_lfp_ups(i) > 1
            cur_inds = (down_trans_inds8(skipped_lfp_ups(i)-1)+buffer_win):(up_trans_inds8(skipped_lfp_ups(i)));
            lfp_prevskipped_up(cur_inds) = true;
        end
        if skipped_lfp_ups(i) < length(up_trans_inds8)
            cur_inds = (down_trans_inds8(skipped_lfp_ups(i))+buffer_win):(up_trans_inds8(skipped_lfp_ups(i)+1));
            lfp_nextskipped_up(cur_inds) = true;
        end
    end
    for i = 1:length(non_skipped_lfp_ups)
%         cur_inds = (up_trans_inds8(non_skipped_lfp_ups(i))+buffer_win):(down_trans_inds8(non_skipped_lfp_ups(i))-buffer_win);
        cur_inds = (up_trans_inds8(non_skipped_lfp_ups(i))+buffer_win):(down_trans_inds8(non_skipped_lfp_ups(i)));
        lfp_nskipped_up(cur_inds) = true;
    end
    for i = 1:length(skipped_lfp_downs)
        if skipped_lfp_downs(i) < length(up_trans_inds8)
%             cur_inds = (down_trans_inds8(skipped_lfp_downs(i))+buffer_win):(up_trans_inds8(skipped_lfp_downs(i)+1)-buffer_win);
            cur_inds = (down_trans_inds8(skipped_lfp_downs(i))+buffer_win):(up_trans_inds8(skipped_lfp_downs(i)+1));
        else
%             cur_inds = (down_trans_inds8(skipped_lfp_downs(i))+buffer_win):(length(lfp_state_number)-buffer_win);
            cur_inds = (down_trans_inds8(skipped_lfp_downs(i))+buffer_win):(length(lfp_state_number));
        end
        lfp_skipped_down(cur_inds) = true;
    end
    for i = 1:length(non_skipped_lfp_downs)
        if non_skipped_lfp_downs(i) < length(up_trans_inds8)
%             cur_inds = (down_trans_inds8(non_skipped_lfp_downs(i))+buffer_win):(up_trans_inds8(non_skipped_lfp_downs(i)+1)-buffer_win);
            cur_inds = (down_trans_inds8(non_skipped_lfp_downs(i))+buffer_win):(up_trans_inds8(non_skipped_lfp_downs(i)+1));
        else
%             cur_inds = (down_trans_inds8(non_skipped_lfp_downs(i))+buffer_win):(length(lfp_state_number)-buffer_win);
            cur_inds = (down_trans_inds8(non_skipped_lfp_downs(i))+buffer_win):(length(lfp_state_number));
        end
        lfp_nskipped_down(cur_inds) = true;
    end
    
    trig_data(d).avg_sk_up_ctx_hf = mean(lfp_hf(lfp_skipped_up));
    trig_data(d).avg_nsk_up_ctx_hf = mean(lfp_hf(lfp_nskipped_up));
    trig_data(d).avg_sk_down_ctx_hf = mean(lfp_hf(lfp_skipped_down));
    trig_data(d).avg_nsk_down_ctx_hf = mean(lfp_hf(lfp_nskipped_down));
    
    trig_data(d).avg_sk_up_ctx_mua = mean(ctx_mua_rate(lfp_skipped_up));
    trig_data(d).avg_nsk_up_ctx_mua = mean(ctx_mua_rate(lfp_nskipped_up));
    trig_data(d).avg_sk_down_ctx_mua = mean(ctx_mua_rate(lfp_skipped_down));
    trig_data(d).avg_nsk_down_ctx_mua = mean(ctx_mua_rate(lfp_nskipped_down));

    trig_data(d).avg_sk_up_hpc_hf2 = mean(hpc_hf2(lfp_skipped_up));
    trig_data(d).avg_nsk_up_hpc_hf2 = mean(hpc_hf2(lfp_nskipped_up));
    trig_data(d).avg_sk_down_hpc_hf2 = mean(hpc_hf2(lfp_skipped_down));
    trig_data(d).avg_nsk_down_hpc_hf2 = mean(hpc_hf2(lfp_nskipped_down));
    trig_data(d).avg_sk_up_hpc_hf3 = mean(hpc_hf3(lfp_skipped_up));
    trig_data(d).avg_nsk_up_hpc_hf3 = mean(hpc_hf3(lfp_nskipped_up));
    trig_data(d).avg_sk_down_hpc_hf3 = mean(hpc_hf3(lfp_skipped_down));
    trig_data(d).avg_nsk_down_hpc_hf3 = mean(hpc_hf3(lfp_nskipped_down));

    trig_data(d).avg_sk_up_hpc_mua = mean(hpc_mua_rate(lfp_skipped_up));
    trig_data(d).avg_nsk_up_hpc_mua = mean(hpc_mua_rate(lfp_nskipped_up));
    trig_data(d).avg_prvsk_up_hpc_mua = mean(hpc_mua_rate(lfp_prevskipped_up));
    trig_data(d).avg_nxtsk_up_hpc_mua = mean(hpc_mua_rate(lfp_nextskipped_up));
    trig_data(d).avg_sk_down_hpc_mua = mean(hpc_mua_rate(lfp_skipped_down));
    trig_data(d).avg_nsk_down_hpc_mua = mean(hpc_mua_rate(lfp_nskipped_down));
    
    lfp_state_vec_rob = lfp_state_vec;
    mp_state_vec_rob = mp_state_vec;
    rob_mindur = round(0.5/Fsd);
    for ii = 1:(length(up_trans_inds8)-1)
        cur_inds = up_trans_inds8(ii):down_trans_inds8(ii);
        if length(cur_inds) < rob_mindur
            lfp_state_vec_rob(cur_inds) = nan;
        end
        cur_inds = down_trans_inds8(ii):(up_trans_inds8(ii+1));
        if length(cur_inds) < rob_mindur
            lfp_state_vec_rob(cur_inds) = nan;
        end
    end
    for ii = 1:(length(up_trans_inds)-1)
        cur_inds = up_trans_inds(ii):down_trans_inds(ii);
        if length(cur_inds) < rob_mindur
            mp_state_vec_rob(cur_inds) = nan;
        end
        cur_inds = down_trans_inds(ii):(up_trans_inds(ii+1));
        if length(cur_inds) < rob_mindur
            mp_state_vec_rob(cur_inds) = nan;
        end
    end        
    
%     trig_data(d).avg_cd_md_hpc_mua = mean(hpc_mua_rate(lfp_state_vec == 0 & mp_state_vec == 0));
%     trig_data(d).avg_cu_mu_hpc_mua = mean(hpc_mua_rate(lfp_state_vec == 1 & mp_state_vec == 1));
%     trig_data(d).avg_cd_mu_hpc_mua = mean(hpc_mua_rate(lfp_state_vec == 0 & mp_state_vec == 1));
%     trig_data(d).avg_cu_md_hpc_mua = mean(hpc_mua_rate(lfp_state_vec == 1 & mp_state_vec == 0));
    trig_data(d).avg_cd_md_hpc_mua = mean(hpc_mua_rate(lfp_state_vec_rob == 0 & mp_state_vec_rob == 0));
    trig_data(d).avg_cu_mu_hpc_mua = mean(hpc_mua_rate(lfp_state_vec_rob == 1 & mp_state_vec_rob == 1));
    trig_data(d).avg_cd_mu_hpc_mua = mean(hpc_mua_rate(lfp_state_vec_rob == 0 & mp_state_vec_rob == 1));
    trig_data(d).avg_cu_md_hpc_mua = mean(hpc_mua_rate(lfp_state_vec_rob == 1 & mp_state_vec_rob == 0));
    
    %%

    trig_data(d).n_skipped_ups = length(skipped_lfp_ups);
    trig_data(d).n_nonskipped_ups = length(non_skipped_lfp_ups);
        
    [trig_data(d).non_sk_utrig_ctx_mua,lags] = get_event_trig_avg(ctx_mua_rate,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag);
    trig_data(d).non_sk_utrig_hpc_mua = get_event_trig_avg(hpc_mua_rate,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag);
    trig_data(d).non_sk_utrig_mp_rate = get_event_trig_avg(mp_spike_rate,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag);
%     [trig_data(d).non_sk_utrig_ctx_mua,lags] = get_event_trig_avg(ctx_mua_rate,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag,[],mp_state_number,2,[],min_samps);
%     trig_data(d).non_sk_utrig_hpc_mua = get_event_trig_avg(hpc_mua_rate,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag,[],lfp_state_number,-1,[],min_samps);
%     trig_data(d).non_sk_utrig_mp_rate = get_event_trig_avg(mp_spike_rate,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag,[],mp_state_number,2,[],min_samps);
 
    trig_data(d).non_sk_utrig_ctx_lf = get_event_trig_avg(lfp_lf,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag);
    trig_data(d).non_sk_utrig_hpc_lf2 = get_event_trig_avg(hpc_lf2,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag);
    trig_data(d).non_sk_utrig_hpc_lf3 = get_event_trig_avg(hpc_lf3,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag);
    trig_data(d).non_sk_utrig_mp_lf = get_event_trig_avg(wcv_lf,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag);
    trig_data(d).non_sk_utrig_mp_dc = get_event_trig_avg(dc_interp_data_reld,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag);

%     trig_data(d).non_sk_utrig_ctx_lf = get_event_trig_avg(lfp_lf,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag,[],mp_state_number,2,[],min_samps);
%     trig_data(d).non_sk_utrig_hpc_lf = get_event_trig_avg(hpc_lf,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag,[],lfp_state_number,-1,[],min_samps);
%     trig_data(d).non_sk_utrig_mp_lf = get_event_trig_avg(wcv_lf,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag,[],mp_state_number,2,[],min_samps);

    trig_data(d).non_sk_utrig_ctx_hf = get_event_trig_avg(lfp_hf,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag);
    trig_data(d).non_sk_utrig_hpc_hf2 = get_event_trig_avg(hpc_hf2,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag);
    trig_data(d).non_sk_utrig_hpc_hf3 = get_event_trig_avg(hpc_hf3,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag);
    trig_data(d).non_sk_utrig_mp_hf = get_event_trig_avg(wcv_hf,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag);
%     trig_data(d).non_sk_utrig_ctx_hf = get_event_trig_avg(lfp_hf,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag,[],mp_state_number,2,[],min_samps);
%     trig_data(d).non_sk_utrig_hpc_hf = get_event_trig_avg(hpc_hf,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag,[],lfp_state_number,-1,[],min_samps);
%     trig_data(d).non_sk_utrig_mp_hf = get_event_trig_avg(wcv_hf,up_trans_inds8(non_skipped_lfp_ups),backlag,forwardlag,[],mp_state_number,2,[],min_samps);

    if length(skipped_lfp_ups) > min_samps
        trig_data(d).sk_utrig_ctx_mua = get_event_trig_avg(ctx_mua_rate,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag);
        trig_data(d).sk_utrig_hpc_mua = get_event_trig_avg(hpc_mua_rate,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag);
        trig_data(d).sk_utrig_mp_rate = get_event_trig_avg(mp_spike_rate,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag);
%         trig_data(d).sk_utrig_ctx_mua = get_event_trig_avg(ctx_mua_rate,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag,[],lfp_state_number,-1,[],min_samps);
%         trig_data(d).sk_utrig_hpc_mua = get_event_trig_avg(hpc_mua_rate,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag,[],mp_state_number,0,[],min_samps);
%         trig_data(d).sk_utrig_mp_rate = get_event_trig_avg(mp_spike_rate,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag,[],mp_state_number,0,[],min_samps);
       
        trig_data(d).sk_utrig_ctx_lf = get_event_trig_avg(lfp_lf,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag);
        trig_data(d).sk_utrig_hpc_lf2 = get_event_trig_avg(hpc_lf2,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag);
        trig_data(d).sk_utrig_hpc_lf3 = get_event_trig_avg(hpc_lf3,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag);
        trig_data(d).sk_utrig_mp_lf = get_event_trig_avg(wcv_lf,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag);
        trig_data(d).sk_utrig_mp_dc = get_event_trig_avg(dc_interp_data_reld,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag);
%         trig_data(d).sk_utrig_ctx_lf = get_event_trig_avg(lfp_lf,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag,[],lfp_state_number,-1,[],min_samps);
%         trig_data(d).sk_utrig_hpc_lf = get_event_trig_avg(hpc_lf,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag,[],mp_state_number,0,[],min_samps);
%         trig_data(d).sk_utrig_mp_lf = get_event_trig_avg(wcv_lf,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag,[],mp_state_number,0,[],min_samps);
        
        trig_data(d).sk_utrig_ctx_hf = get_event_trig_avg(lfp_hf,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag);
        trig_data(d).sk_utrig_hpc_hf2 = get_event_trig_avg(hpc_hf2,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag);
        trig_data(d).sk_utrig_hpc_hf3 = get_event_trig_avg(hpc_hf3,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag);
        trig_data(d).sk_utrig_mp_hf = get_event_trig_avg(wcv_hf,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag);
%         trig_data(d).sk_utrig_ctx_hf = get_event_trig_avg(lfp_hf,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag,[],lfp_state_number,-1,[],min_samps);
%         trig_data(d).sk_utrig_hpc_hf = get_event_trig_avg(hpc_hf,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag,[],mp_state_number,0,[],min_samps);
%         trig_data(d).sk_utrig_mp_hf = get_event_trig_avg(wcv_hf,up_trans_inds8(skipped_lfp_ups),backlag,forwardlag,[],mp_state_number,0,[],min_samps);
    else
        trig_data(d).sk_utrig_ctx_mua = nan(length(lags),1);
        trig_data(d).sk_utrig_hpc_mua = nan(length(lags),1);
        trig_data(d).sk_utrig_mp_rate = nan(length(lags),1);
        trig_data(d).sk_utrig_ctx_lf = nan(length(lags),1);
        trig_data(d).sk_utrig_hpc_lf2 = nan(length(lags),1);
        trig_data(d).sk_utrig_hpc_lf3 = nan(length(lags),1);
        trig_data(d).sk_utrig_mp_lf = nan(length(lags),1);
        trig_data(d).sk_utrig_mp_dc = nan(length(lags),1);
        trig_data(d).sk_utrig_ctx_hf = nan(length(lags),1);
        trig_data(d).sk_utrig_hpc_hf2 = nan(length(lags),1);
        trig_data(d).sk_utrig_hpc_hf3 = nan(length(lags),1);
        trig_data(d).sk_utrig_mp_hf = nan(length(lags),1);
    end
        
    %%
    skipped_downs_nextups = skipped_lfp_downs + 1;
    skipped_downs_nextups(skipped_downs_nextups > length(up_trans_inds8)) = [];
    non_skipped_downs_nextups = non_skipped_lfp_downs + 1;
    non_skipped_downs_nextups(non_skipped_downs_nextups > length(up_trans_inds8)) = [];
    
    trig_data(d).n_skipped_downs = length(skipped_lfp_downs);
    trig_data(d).n_nonskipped_downs = length(non_skipped_lfp_downs);
    
    [trig_data(d).non_sk_dtrig_ctx_mua,lags] = get_event_trig_avg(ctx_mua_rate,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag);
    trig_data(d).non_sk_dtrig_hpc_mua = get_event_trig_avg(hpc_mua_rate,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag);
    trig_data(d).non_sk_dtrig_mp_rate = get_event_trig_avg(mp_spike_rate,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag);
%     [trig_data(d).non_sk_dtrig_ctx_mua,lags] = get_event_trig_avg(ctx_mua_rate,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag,[],lfp_state_number,-1,[],min_samps);
%     trig_data(d).non_sk_dtrig_hpc_mua = get_event_trig_avg(hpc_mua_rate,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag,[],mp_state_number,2,[],min_samps);
%     trig_data(d).non_sk_dtrig_mp_rate = get_event_trig_avg(mp_spike_rate,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag,[],mp_state_number,2,[],min_samps);
 
%     [trig_data(d).non_sk_dtrignext_ctx_mua,lags] = get_event_trig_avg(ctx_mua_rate,up_trans_inds8(non_skipped_downs_nextups),backlag,forwardlag,[],lfp_state_number,-1,[],min_samps);
%     [trig_data(d).non_sk_dtrignext_ctx_hf,lags] = get_event_trig_avg(lfp_hf,up_trans_inds8(non_skipped_downs_nextups),backlag,forwardlag,[],lfp_state_number,-1,[],min_samps);
%     [trig_data(d).non_sk_dtrignext_ctx_lf,lags] = get_event_trig_avg(lfp_lf,up_trans_inds8(non_skipped_downs_nextups),backlag,forwardlag,[],lfp_state_number,-1,[],min_samps);

        trig_data(d).non_sk_dtrig_ctx_lf = get_event_trig_avg(lfp_lf,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag);
    trig_data(d).non_sk_dtrig_hpc_lf2 = get_event_trig_avg(hpc_lf2,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag);
    trig_data(d).non_sk_dtrig_hpc_lf3 = get_event_trig_avg(hpc_lf3,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag);
    trig_data(d).non_sk_dtrig_mp_lf = get_event_trig_avg(wcv_lf,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag);
    trig_data(d).non_sk_dtrig_mp_dc = get_event_trig_avg(dc_interp_data_reld,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag);
%     trig_data(d).non_sk_dtrig_ctx_lf = get_event_trig_avg(lfp_lf,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag,[],lfp_state_number,-1,[],min_samps);
%     trig_data(d).non_sk_dtrig_hpc_lf = get_event_trig_avg(hpc_lf,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag,[],mp_state_number,2,[],min_samps);
%     trig_data(d).non_sk_dtrig_mp_lf = get_event_trig_avg(wcv_lf,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag,[],mp_state_number,2,[],min_samps);

    trig_data(d).non_sk_dtrig_ctx_hf = get_event_trig_avg(lfp_hf,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag);
    trig_data(d).non_sk_dtrig_hpc_hf2 = get_event_trig_avg(hpc_hf2,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag);
    trig_data(d).non_sk_dtrig_hpc_hf3 = get_event_trig_avg(hpc_hf3,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag);
    trig_data(d).non_sk_dtrig_mp_hf = get_event_trig_avg(wcv_hf,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag);
%     trig_data(d).non_sk_dtrig_ctx_hf = get_event_trig_avg(lfp_hf,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag,[],lfp_state_number,-1,[],min_samps);
%     trig_data(d).non_sk_dtrig_hpc_hf = get_event_trig_avg(hpc_hf,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag,[],mp_state_number,2,[],min_samps);
%     trig_data(d).non_sk_dtrig_mp_hf = get_event_trig_avg(wcv_hf,down_trans_inds8(non_skipped_lfp_downs),backlag,forwardlag,[],mp_state_number,2,[],min_samps);

    if length(skipped_lfp_downs) > min_samps
        trig_data(d).sk_dtrig_ctx_mua = get_event_trig_avg(ctx_mua_rate,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag);
        trig_data(d).sk_dtrig_hpc_mua = get_event_trig_avg(hpc_mua_rate,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag);
        trig_data(d).sk_dtrig_mp_rate = get_event_trig_avg(mp_spike_rate,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag);
%         trig_data(d).sk_dtrig_ctx_mua = get_event_trig_avg(ctx_mua_rate,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag,[],lfp_state_number,-1,[],min_samps);
%         trig_data(d).sk_dtrig_hpc_mua = get_event_trig_avg(hpc_mua_rate,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag,[],mp_state_number,0,[],min_samps);
%         trig_data(d).sk_dtrig_mp_rate = get_event_trig_avg(mp_spike_rate,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag,[],mp_state_number,0,[],min_samps);

%         trig_data(d).sk_dtrignext_ctx_mua = get_event_trig_avg(ctx_mua_rate,up_trans_inds8(skipped_downs_nextups),backlag,forwardlag,[],lfp_state_number,-1,[],min_samps);
%         trig_data(d).sk_dtrignext_ctx_hf = get_event_trig_avg(lfp_hf,up_trans_inds8(skipped_downs_nextups),backlag,forwardlag,[],lfp_state_number,-1,[],min_samps);
%         trig_data(d).sk_dtrignext_ctx_lf = get_event_trig_avg(lfp_lf,up_trans_inds8(skipped_downs_nextups),backlag,forwardlag,[],lfp_state_number,-1,[],min_samps);

trig_data(d).sk_dtrig_ctx_lf = get_event_trig_avg(lfp_lf,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag);
trig_data(d).sk_dtrig_hpc_lf2 = get_event_trig_avg(hpc_lf2,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag);
        trig_data(d).sk_dtrig_hpc_lf3 = get_event_trig_avg(hpc_lf3,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag);
        trig_data(d).sk_dtrig_mp_lf = get_event_trig_avg(wcv_lf,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag);
        trig_data(d).sk_dtrig_mp_dc = get_event_trig_avg(dc_interp_data_reld,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag);
%         trig_data(d).sk_dtrig_ctx_lf = get_event_trig_avg(lfp_lf,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag,[],lfp_state_number,-1,[],min_samps);
%         trig_data(d).sk_dtrig_hpc_lf = get_event_trig_avg(hpc_lf,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag,[],mp_state_number,0,[],min_samps);
%         trig_data(d).sk_dtrig_mp_lf = get_event_trig_avg(wcv_lf,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag,[],mp_state_number,0,[],min_samps);
        
        trig_data(d).sk_dtrig_ctx_hf = get_event_trig_avg(lfp_hf,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag);
        trig_data(d).sk_dtrig_hpc_hf2 = get_event_trig_avg(hpc_hf2,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag);
        trig_data(d).sk_dtrig_hpc_hf3 = get_event_trig_avg(hpc_hf3,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag);
        trig_data(d).sk_dtrig_mp_hf = get_event_trig_avg(wcv_hf,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag);
%         trig_data(d).sk_dtrig_ctx_hf = get_event_trig_avg(lfp_hf,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag,[],lfp_state_number,-1,[],min_samps);
%         trig_data(d).sk_dtrig_hpc_hf = get_event_trig_avg(hpc_hf,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag,[],mp_state_number,0,[],min_samps);
%         trig_data(d).sk_dtrig_mp_hf = get_event_trig_avg(wcv_hf,down_trans_inds8(skipped_lfp_downs),backlag,forwardlag,[],mp_state_number,0,[],min_samps);
    else
        trig_data(d).sk_dtrig_ctx_mua = nan(length(lags),1);
        trig_data(d).sk_dtrig_hpc_mua = nan(length(lags),1);
        trig_data(d).sk_dtrig_mp_rate = nan(length(lags),1);
        trig_data(d).sk_dtrig_ctx_lf = nan(length(lags),1);
        trig_data(d).sk_dtrig_hpc_lf2 = nan(length(lags),1);
        trig_data(d).sk_dtrig_hpc_lf3 = nan(length(lags),1);
        trig_data(d).sk_dtrig_mp_lf = nan(length(lags),1);
        trig_data(d).sk_dtrig_mp_dc = nan(length(lags),1);
        trig_data(d).sk_dtrig_ctx_hf = nan(length(lags),1);
        trig_data(d).sk_dtrig_hpc_hf2 = nan(length(lags),1);
        trig_data(d).sk_dtrig_hpc_hf3 = nan(length(lags),1);
        trig_data(d).sk_dtrig_mp_hf = nan(length(lags),1);
    end
        
end 

%%
cd C:\WC_Germany\final_pdown_analysis\
save final_trig_avg_data4_nocon_nobuff_newpeaks_wch5_dcmp_noORDrej trig_data lags Fsd buffer_win

%%
% close all
 good_ctx_mua = ~isnan([trig_data(:).avg_ctx_mua_rate]);
good_hpc_mua = ~isnan([trig_data(:).avg_hpc_mua_rate]);

data_ids = [data(:).id];
l3mec = find(strcmp({data.loc},'MEC'));
l3lec = find(strcmp({data.loc},'LEC'));
l3mec(~ismember(data_ids(l3mec),clear_l3)) = [];

%%
% n_skipped_ups = [trig_data(:).n_skipped_ups];
% n_skipped_downs = [trig_data(:).n_skipped_downs];
n_pers_downs = cellfun(@(x) length(x),{core_data(:).rt2_downs});
n_pers_ups = cellfun(@(x) length(x),{core_data(:).rt2_ups});

min_rate = 1;
avg_hpc_rate = [trig_data(:).avg_hpc_mua_rate];
avg_ctx_rate = [trig_data(:).avg_ctx_mua_rate];
l3mec_ctxmua = l3mec(good_ctx_mua(l3mec) & avg_ctx_rate(l3mec) >= min_rate);
l3mec_hpcmua = l3mec(good_hpc_mua(l3mec) & avg_hpc_rate(l3mec) >= min_rate);

min_states = 5;
% l3mec_pups = l3mec(n_skipped_downs(l3mec) >= min_states);
% l3mec_pdowns = l3mec(n_skipped_ups(l3mec) >= min_states);
l3mec_pups = l3mec(n_pers_ups(l3mec) >= min_states);
l3mec_pdowns = l3mec(n_pers_downs(l3mec) >= min_states);
l3mec_both = intersect(l3mec_pups,l3mec_pdowns);

% l3mec_hpcmua_pdowns = l3mec_hpcmua(n_skipped_ups(l3mec_hpcmua) >= min_states);
% l3mec_hpcmua_pups = l3mec_hpcmua(n_skipped_downs(l3mec_hpcmua) >= min_states);
% l3mec_ctxmua_pdowns = l3mec_ctxmua(n_skipped_ups(l3mec_ctxmua) >= min_states);
% l3mec_ctxmua_pups = l3mec_ctxmua(n_skipped_downs(l3mec_ctxmua) >= min_states);
l3mec_hpcmua_pdowns = l3mec_hpcmua(ismember(l3mec_hpcmua,l3mec_pdowns));
l3mec_hpcmua_pups = l3mec_hpcmua(ismember(l3mec_hpcmua,l3mec_pups));
l3mec_hpcmua_both = l3mec_hpcmua(ismember(l3mec_hpcmua,l3mec_both));
l3mec_ctxmua_pdowns = l3mec_ctxmua(ismember(l3mec_ctxmua,l3mec_pdowns));
l3mec_ctxmua_pups = l3mec_ctxmua(ismember(l3mec_ctxmua,l3mec_pups));

%%
sk_up_ctx_hf = [trig_data(:).avg_sk_up_ctx_hf];
nonsk_up_ctx_hf = [trig_data(:).avg_nsk_up_ctx_hf];
sk_up_ctx_mua = [trig_data(:).avg_sk_up_ctx_mua];
nonsk_up_ctx_mua = [trig_data(:).avg_nsk_up_ctx_mua];
sk_up_hpc_mua = [trig_data(:).avg_sk_up_hpc_mua];
nonsk_up_hpc_mua = [trig_data(:).avg_nsk_up_hpc_mua];
% prevsk_up_hpc_mua = [trig_data(:).avg_prvsk_up_hpc_mua];
% nextsk_up_hpc_mua = [trig_data(:).avg_nxtsk_up_hpc_mua];

sk_down_ctx_hf = [trig_data(:).avg_sk_down_ctx_hf];
nonsk_down_ctx_hf = [trig_data(:).avg_nsk_down_ctx_hf];
sk_down_ctx_mua = [trig_data(:).avg_sk_down_ctx_mua];
nonsk_down_ctx_mua = [trig_data(:).avg_nsk_down_ctx_mua];
sk_down_hpc_mua = [trig_data(:).avg_sk_down_hpc_mua];
nonsk_down_hpc_mua = [trig_data(:).avg_nsk_down_hpc_mua];

avg_cd_md_hpc_mua = [trig_data(:).avg_cd_md_hpc_mua];
avg_cu_mu_hpc_mua = [trig_data(:).avg_cu_mu_hpc_mua];
avg_cd_mu_hpc_mua = [trig_data(:).avg_cd_mu_hpc_mua];
avg_cu_md_hpc_mua = [trig_data(:).avg_cu_md_hpc_mua];

avg_pdown_effect = sk_up_hpc_mua-nonsk_up_hpc_mua;
avg_pup_effect = sk_down_hpc_mua-nonsk_down_hpc_mua;
avg_ctxup_effect = sk_up_hpc_mua - avg_cd_md_hpc_mua;
avg_ctxdown_effect = sk_down_hpc_mua - avg_cu_mu_hpc_mua;


%%
% all_data = [nonsk_down_hpc_mua(l3mec_hpcmua_pups)'; sk_down_hpc_mua(l3mec_hpcmua_pups)'; ...
%     nonsk_up_hpc_mua(l3mec_hpcmua_pdowns)'; sk_up_hpc_mua(l3mec_hpcmua_pdowns)'];
% all_data = [avg_cd_md_hpc_mua(l3mec_hpcmua)'; avg_cd_mu_hpc_mua(l3mec_hpcmua)'; ...
%     avg_cu_mu_hpc_mua(l3mec_hpcmua)'; avg_cu_md_hpc_mua(l3mec_hpcmua)'];
% all_data = [avg_pdown_effect(l3mec_hpcmua)'; avg_pup_effect(l3mec_hpcmua)'; ...
%     avg_ctxup_effect(l3mec_hpcmua)'; avg_ctxdown_effect(l3mec_hpcmua)'];
% groups = [ones(length(l3mec_hpcmua),1); 2*ones(length(l3mec_hpcmua),1); ...
%     3*ones(length(l3mec_hpcmua),1); 4*ones(length(l3mec_hpcmua),1)];
% all_data = [nonsk_down_hpc_mua(l3mec_hpcmua_both)'; sk_down_hpc_mua(l3mec_hpcmua_both)'; nonsk_up_hpc_mua(l3mec_hpcmua_both)'; sk_up_hpc_mua(l3mec_hpcmua_both)'];
% groups = [ones(length(l3mec_hpcmua_both),1); 2*ones(length(l3mec_hpcmua_both),1); 3*ones(length(l3mec_hpcmua_both),1); 4*ones(length(l3mec_hpcmua_both),1)];

all_data = [avg_pdown_effect(l3mec_hpcmua_pdowns)'; avg_pup_effect(l3mec_hpcmua_pups)'; ...
    avg_ctxup_effect(l3mec_hpcmua_pdowns)'; avg_ctxdown_effect(l3mec_hpcmua_pups)'];
groups = [ones(length(l3mec_hpcmua_pdowns),1); 2*ones(length(l3mec_hpcmua_pups),1); ...
    3*ones(length(l3mec_hpcmua_pdowns),1); 4*ones(length(l3mec_hpcmua_pups),1)];

h = figure();
boxplot(all_data,groups);
ylabel('Hpc MUA (z)');
xl = xlim();
line(xl,[0 0],'color','k','linestyle','--');

%%
fig_width = 3.5; rel_height = 0.8;
figufy(h);
fname = [fig_dir 'hpcmua_boxplot_5s_np_wch5.pdf'];
exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h);

%% UTRIG SPIKE RATES
close all
xr = [-1.5 3];

% l3mec_non_sk_utrig_ctxmua = reshape([trig_data(l3mec_ctxmua).non_sk_utrig_ctx_mua],[length(lags) length(l3mec_ctxmua)]);
% l3mec_sk_utrig_ctxmua = reshape([trig_data(l3mec_ctxmua).sk_utrig_ctx_mua],[length(lags) length(l3mec_ctxmua)]);
l3mec_non_sk_utrig_ctxmua = reshape([trig_data(l3mec_ctxmua_pdowns).non_sk_utrig_ctx_mua],[length(lags) length(l3mec_ctxmua_pdowns)]);
l3mec_sk_utrig_ctxmua = reshape([trig_data(l3mec_ctxmua_pdowns).sk_utrig_ctx_mua],[length(lags) length(l3mec_ctxmua_pdowns)]);
f1 = figure(); 
subplot(3,1,1);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_utrig_ctxmua,2),nanstd(l3mec_non_sk_utrig_ctxmua,[],2)./sqrt(sum(~isnan(l3mec_non_sk_utrig_ctxmua),2)));
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_ctxmua,2),nanstd(l3mec_sk_utrig_ctxmua,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_ctxmua),2)),{'color','r'});
xlim(xr);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
title('Ctx MUA');
ylabel('MUA rate (z)');

l3mec_non_sk_utrig_hpcmua = reshape([trig_data(l3mec_hpcmua_pdowns).non_sk_utrig_hpc_mua],[length(lags) length(l3mec_hpcmua_pdowns)]);
l3mec_sk_utrig_hpcmua = reshape([trig_data(l3mec_hpcmua_pdowns).sk_utrig_hpc_mua],[length(lags) length(l3mec_hpcmua_pdowns)]);
subplot(3,1,2);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_utrig_hpcmua,2),nanstd(l3mec_non_sk_utrig_hpcmua,[],2)./sqrt(sum(~isnan(l3mec_non_sk_utrig_hpcmua),2)),{'color','k'});
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_hpcmua,2),nanstd(l3mec_sk_utrig_hpcmua,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_hpcmua),2)),{'color','r'});
xlim(xr);
ylim([-0.5 0.5])
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
title('Hpc MUA');
ylabel('MUA rate (z)');

l3mec_non_sk_utrig_mprate = reshape([trig_data(l3mec_pdowns).non_sk_utrig_mp_rate],[length(lags) length(l3mec_pdowns)]);
l3mec_sk_utrig_mprate = reshape([trig_data(l3mec_pdowns).sk_utrig_mp_rate],[length(lags) length(l3mec_pdowns)]);
l3mec_meanrate = [trig_data(l3mec_pdowns).avg_mp_rate];
l3mec_stdrate = [trig_data(l3mec_pdowns).std_mp_rate];
l3mec_non_sk_utrig_mprate = bsxfun(@times,l3mec_non_sk_utrig_mprate,l3mec_stdrate);
l3mec_non_sk_utrig_mprate = bsxfun(@plus,l3mec_non_sk_utrig_mprate,l3mec_meanrate);
l3mec_sk_utrig_mprate = bsxfun(@times,l3mec_sk_utrig_mprate,l3mec_stdrate);
l3mec_sk_utrig_mprate = bsxfun(@plus,l3mec_sk_utrig_mprate,l3mec_meanrate);


% l3mec_non_sk_utrig_mprate = reshape([trig_data(l3mec_hpcmua_pdowns).non_sk_utrig_mp_rate],[length(lags) length(l3mec_hpcmua_pdowns)]);
% l3mec_sk_utrig_mprate = reshape([trig_data(l3mec_hpcmua_pdowns).sk_utrig_mp_rate],[length(lags) length(l3mec_hpcmua_pdowns)]);
subplot(3,1,3);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_utrig_mprate,2),nanstd(l3mec_non_sk_utrig_mprate,[],2)./sqrt(sum(~isnan(l3mec_non_sk_utrig_mprate),2)),{'color','k'});
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_mprate,2),nanstd(l3mec_sk_utrig_mprate,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_mprate),2)),{'color','r'});
xlim(xr);
% ylim([-0.75 0.75])
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
title('MP spiking');
ylabel('spike rate (z)');
xlabel('Time (s)');

% figure
% shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_hpcmua - l3mec_non_sk_utrig_hpcmua,2),nanstd(l3mec_sk_utrig_hpcmua - l3mec_non_sk_utrig_hpcmua,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_hpcmua),2)));
% hold on
% shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_ctxmua - l3mec_non_sk_utrig_ctxmua,2),nanstd(l3mec_sk_utrig_ctxmua - l3mec_non_sk_utrig_ctxmua,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_ctxmua),2)),{'color','r'});

%%
fig_width = 3.5; rel_height = 2.4;
figufy(f1);
fname = [fig_dir 'utrig_spikerates_fin_5s_np_wch5_v2.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);
%% DOWN-triggered spike rates
xr = [-1.5 3];

l3mec_non_sk_dtrig_ctxmua = reshape([trig_data(l3mec_ctxmua_pups).non_sk_dtrig_ctx_mua],[length(lags) length(l3mec_ctxmua_pups)]);
l3mec_sk_dtrig_ctxmua = reshape([trig_data(l3mec_ctxmua_pups).sk_dtrig_ctx_mua],[length(lags) length(l3mec_ctxmua_pups)]);

f1 = figure(); 
subplot(3,1,1); hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_dtrig_ctxmua,2),nanstd(l3mec_non_sk_dtrig_ctxmua,[],2)./sqrt(sum(~isnan(l3mec_non_sk_dtrig_ctxmua),2)));
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_dtrig_ctxmua,2),nanstd(l3mec_sk_dtrig_ctxmua,[],2)./sqrt(sum(~isnan(l3mec_sk_dtrig_ctxmua),2)),{'color','r'});
xlim(xr);
ylim([-1 1.5])
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
title('Ctx MUA');
ylabel('MUA rate (z)');

l3mec_non_sk_dtrig_hpcmua = reshape([trig_data(l3mec_hpcmua_pups).non_sk_dtrig_hpc_mua],[length(lags) length(l3mec_hpcmua_pups)]);
l3mec_sk_dtrig_hpcmua = reshape([trig_data(l3mec_hpcmua_pups).sk_dtrig_hpc_mua],[length(lags) length(l3mec_hpcmua_pups)]);

subplot(3,1,2); hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_dtrig_hpcmua,2),nanstd(l3mec_non_sk_dtrig_hpcmua,[],2)./sqrt(sum(~isnan(l3mec_non_sk_dtrig_hpcmua),2)));
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_dtrig_hpcmua,2),nanstd(l3mec_sk_dtrig_hpcmua,[],2)./sqrt(sum(~isnan(l3mec_sk_dtrig_hpcmua),2)),{'color','r'});
xlim(xr);
xlim(xr);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
title('Hpc MUA');
ylabel('MUA rate (z)');

l3mec_non_sk_dtrig_mprate = reshape([trig_data(l3mec_pups).non_sk_dtrig_mp_rate],[length(lags) length(l3mec_pups)]);
l3mec_sk_dtrig_mprate = reshape([trig_data(l3mec_pups).sk_dtrig_mp_rate],[length(lags) length(l3mec_pups)]);
l3mec_meanrate = [trig_data(l3mec_pups).avg_mp_rate];
l3mec_stdrate = [trig_data(l3mec_pups).std_mp_rate];
l3mec_non_sk_dtrig_mprate = bsxfun(@times,l3mec_non_sk_dtrig_mprate,l3mec_stdrate);
l3mec_non_sk_dtrig_mprate = bsxfun(@plus,l3mec_non_sk_dtrig_mprate,l3mec_meanrate);
l3mec_sk_dtrig_mprate = bsxfun(@times,l3mec_sk_dtrig_mprate,l3mec_stdrate);
l3mec_sk_dtrig_mprate = bsxfun(@plus,l3mec_sk_dtrig_mprate,l3mec_meanrate);

subplot(3,1,3); hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_dtrig_mprate,2),nanstd(l3mec_non_sk_dtrig_mprate,[],2)./sqrt(sum(~isnan(l3mec_non_sk_dtrig_mprate),2)));
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_dtrig_mprate,2),nanstd(l3mec_sk_dtrig_mprate,[],2)./sqrt(sum(~isnan(l3mec_sk_dtrig_mprate),2)),{'color','r'});
xlim(xr);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
title('MP spiking');
ylabel('spike rate (z)');
xlabel('Time (s)');

%%
fig_width = 3.5; rel_height = 2.4;
figufy(f1);
fname = [fig_dir 'dtrig_spikerates_fin_5s_np_wch5_v2.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

%% UTRIG HF POWER
% close all
xr = [-1.5 3];

% good_hpc_lfp = ~isnan([data(:).hpc_lfp]);
% good_hpc_lfp = true(length(data),1);
good_hpc_lfp = ~isnan(peak_hpcmua_loc);
l3mec_pdown_hpclfp = l3mec_pdowns(good_hpc_lfp(l3mec_pdowns));

l3mec_non_sk_utrig_ctxhf = reshape([trig_data(l3mec_pdowns).non_sk_utrig_ctx_hf],[length(lags) length(l3mec_pdowns)]);
l3mec_sk_utrig_ctxhf = reshape([trig_data(l3mec_pdowns).sk_utrig_ctx_hf],[length(lags) length(l3mec_pdowns)]);
f1 = figure(); 
subplot(3,1,1);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_utrig_ctxhf,2),nanstd(l3mec_non_sk_utrig_ctxhf,[],2)./sqrt(sum(~isnan(l3mec_non_sk_utrig_ctxhf),2)));
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_ctxhf,2),nanstd(l3mec_sk_utrig_ctxhf,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_ctxhf),2)),{'color','r'});
xlim(xr);
ylim([-1 1.5])
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
title('Ctx HF');
ylabel('HF power (z)');


% l3mec_non_sk_utrig_hpchf = reshape([trig_data(l3mec_pdowns).non_sk_utrig_hpc_hf2],[length(lags) length(l3mec_pdowns)]);
% l3mec_sk_utrig_hpchf = reshape([trig_data(l3mec_pdowns).sk_utrig_hpc_hf2],[length(lags) length(l3mec_pdowns)]);
l3mec_non_sk_utrig_hpchf = reshape([trig_data(l3mec_pdown_hpclfp).non_sk_utrig_hpc_hf3],[length(lags) length(l3mec_pdown_hpclfp)]);
l3mec_sk_utrig_hpchf = reshape([trig_data(l3mec_pdown_hpclfp).sk_utrig_hpc_hf3],[length(lags) length(l3mec_pdown_hpclfp)]);
subplot(3,1,2);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_utrig_hpchf,2),nanstd(l3mec_non_sk_utrig_hpchf,[],2)./sqrt(sum(~isnan(l3mec_non_sk_utrig_hpchf),2)),{'color','k'});
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_hpchf,2),nanstd(l3mec_sk_utrig_hpchf,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_hpchf),2)),{'color','r'});
xlim(xr);
xlim(xr);
ylim([-0.75 1.25])
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
title('Hpf HF');
ylabel('HF power (z)');

l3mec_non_sk_utrig_mphf = reshape([trig_data(l3mec_pdowns).non_sk_utrig_mp_hf],[length(lags) length(l3mec_pdowns)]);
l3mec_sk_utrig_mphf = reshape([trig_data(l3mec_pdowns).sk_utrig_mp_hf],[length(lags) length(l3mec_pdowns)]);
subplot(3,1,3);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_utrig_mphf,2),nanstd(l3mec_non_sk_utrig_mphf,[],2)./sqrt(sum(~isnan(l3mec_non_sk_utrig_mphf),2)),{'color','k'});
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_mphf,2),nanstd(l3mec_sk_utrig_mphf,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_mphf),2)),{'color','r'});
xlim(xr);
xlim(xr);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
title('MP HF');
ylabel('HF power (z)');

%%
fig_width = 3.5; rel_height = 2.4;
figufy(f1);
fname = [fig_dir 'utrig_HFpow_fin_5s_np_wch5.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

%% DTRIG HF POW
xr = [-1.5 3];
close all

% good_hpc_lfp = ~isnan([data(:).hpc_lfp]);
good_hpc_lfp = ~isnan(peak_hpcmua_loc);
l3mec_pup_hpclfp = l3mec_pups(good_hpc_lfp(l3mec_pups));

l3mec_non_sk_dtrig_ctxhf = reshape([trig_data(l3mec_pups).non_sk_dtrig_ctx_hf],[length(lags) length(l3mec_pups)]);
l3mec_sk_dtrig_ctxhf = reshape([trig_data(l3mec_pups).sk_dtrig_ctx_hf],[length(lags) length(l3mec_pups)]);
f1 = figure(); 
subplot(3,1,1);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_dtrig_ctxhf,2),nanstd(l3mec_non_sk_dtrig_ctxhf,[],2)./sqrt(sum(~isnan(l3mec_non_sk_dtrig_ctxhf),2)));
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_dtrig_ctxhf,2),nanstd(l3mec_sk_dtrig_ctxhf,[],2)./sqrt(sum(~isnan(l3mec_sk_dtrig_ctxhf),2)),{'color','r'});
xlim(xr);
ylim([-1 1.25])
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
title('Ctx HF');
ylabel('HF power (z)');

l3mec_non_sk_dtrig_hpchf = reshape([trig_data(l3mec_pup_hpclfp).non_sk_dtrig_hpc_hf3],[length(lags) length(l3mec_pup_hpclfp)]);
l3mec_sk_dtrig_hpchf = reshape([trig_data(l3mec_pup_hpclfp).sk_dtrig_hpc_hf3],[length(lags) length(l3mec_pup_hpclfp)]);
% l3mec_non_sk_dtrig_hpchf = reshape([trig_data(l3mec_pups).non_sk_dtrig_hpc_hf3],[length(lags) length(l3mec_pups)]);
% l3mec_sk_dtrig_hpchf = reshape([trig_data(l3mec_pups).sk_dtrig_hpc_hf3],[length(lags) length(l3mec_pups)]);
subplot(3,1,2);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_dtrig_hpchf,2),nanstd(l3mec_non_sk_dtrig_hpchf,[],2)./sqrt(sum(~isnan(l3mec_non_sk_dtrig_hpchf),2)),{'color','k'});
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_dtrig_hpchf,2),nanstd(l3mec_sk_dtrig_hpchf,[],2)./sqrt(sum(~isnan(l3mec_sk_dtrig_hpchf),2)),{'color','r'});
xlim(xr);
ylim([-0.75 1])
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
title('Hpc HF');
ylabel('HF power (z)');

l3mec_non_sk_dtrig_mphf = reshape([trig_data(l3mec_pups).non_sk_dtrig_mp_hf],[length(lags) length(l3mec_pups)]);
l3mec_sk_dtrig_mphf = reshape([trig_data(l3mec_pups).sk_dtrig_mp_hf],[length(lags) length(l3mec_pups)]);
subplot(3,1,3);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_dtrig_mphf,2),nanstd(l3mec_non_sk_dtrig_mphf,[],2)./sqrt(sum(~isnan(l3mec_non_sk_dtrig_mphf),2)),{'color','k'});
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_dtrig_mphf,2),nanstd(l3mec_sk_dtrig_mphf,[],2)./sqrt(sum(~isnan(l3mec_sk_dtrig_mphf),2)),{'color','r'});
xlim(xr);
ylim([-0.75 1])
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
title('MP HF');
ylabel('HF power (z)');

%%
fig_width = 3.5; rel_height = 2.4;
figufy(f1);
fname = [fig_dir 'dtrig_HFpow_fin_5s_np_wch5.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

%% UTRIG LF AMPLITUDE
% close all
xr = [-1.5 3];

% good_hpc_lfp = ~isnan([data(:).hpc_lfp]);
good_hpc_lfp = ~isnan(peak_hpcmua_loc);
l3mec_pdown_hpclfp = l3mec_pdowns(good_hpc_lfp(l3mec_pdowns));
% l3mec_pdown_hpclfp = l3mec_pdowns;

l3mec_non_sk_utrig_ctx = reshape([trig_data(l3mec_pdowns).non_sk_utrig_ctx_lf],[length(lags) length(l3mec_pdowns)]);
l3mec_sk_utrig_ctx = reshape([trig_data(l3mec_pdowns).sk_utrig_ctx_lf],[length(lags) length(l3mec_pdowns)]);
f1 = figure(); 
subplot(3,1,1);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_utrig_ctx,2),nanstd(l3mec_non_sk_utrig_ctx,[],2)./sqrt(sum(~isnan(l3mec_non_sk_utrig_ctx),2)));
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_ctx,2),nanstd(l3mec_sk_utrig_ctx,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_ctx),2)),{'color','r'});
xlim(xr);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
title('Ctx LF');
ylabel('LF Amp (z)');


l3mec_non_sk_utrig_hpc = reshape([trig_data(l3mec_pdown_hpclfp).non_sk_utrig_hpc_lf3],[length(lags) length(l3mec_pdown_hpclfp)]);
l3mec_sk_utrig_hpc = reshape([trig_data(l3mec_pdown_hpclfp).sk_utrig_hpc_lf3],[length(lags) length(l3mec_pdown_hpclfp)]);
subplot(3,1,2);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_utrig_hpc,2),nanstd(l3mec_non_sk_utrig_hpc,[],2)./sqrt(sum(~isnan(l3mec_non_sk_utrig_hpc),2)),{'color','k'});
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_hpc,2),nanstd(l3mec_sk_utrig_hpc,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_hpc),2)),{'color','r'});
xlim(xr);
xlim(xr);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
title('Hpc LF');
ylabel('LF Amp (z)');

% l3mec_non_sk_utrig_mp = reshape([trig_data(l3mec_pdowns).non_sk_utrig_mp_lf],[length(lags) length(l3mec_pdowns)]);
% l3mec_sk_utrig_mp = reshape([trig_data(l3mec_pdowns).sk_utrig_mp_lf],[length(lags) length(l3mec_pdowns)]);
l3mec_non_sk_utrig_mp = reshape([trig_data(l3mec_pdowns).non_sk_utrig_mp_dc],[length(lags) length(l3mec_pdowns)]);
l3mec_sk_utrig_mp = reshape([trig_data(l3mec_pdowns).sk_utrig_mp_dc],[length(lags) length(l3mec_pdowns)]);
subplot(3,1,3);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_utrig_mp,2),nanstd(l3mec_non_sk_utrig_mp,[],2)./sqrt(sum(~isnan(l3mec_non_sk_utrig_mp),2)),{'color','k'});
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_utrig_mp,2),nanstd(l3mec_sk_utrig_mp,[],2)./sqrt(sum(~isnan(l3mec_sk_utrig_mp),2)),{'color','r'});
xlim(xr);
xlim(xr);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
title('MP LF');
ylabel('LF Amp (z)');

%%
fig_width = 3.5; rel_height = 2.4;
figufy(f1);
fname = [fig_dir 'utrig_LFamp_5s_wch5_dcmp.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

%% DTRIG LF Amp
xr = [-1.5 3];

% good_hpc_lfp = ~isnan([data(:).hpc_lfp]);
good_hpc_lfp = ~isnan(peak_hpcmua_loc);
l3mec_pup_hpclfp = l3mec_pups(good_hpc_lfp(l3mec_pups));

l3mec_non_sk_dtrig_ctx = reshape([trig_data(l3mec_pups).non_sk_dtrig_ctx_lf],[length(lags) length(l3mec_pups)]);
l3mec_sk_dtrig_ctx = reshape([trig_data(l3mec_pups).sk_dtrig_ctx_lf],[length(lags) length(l3mec_pups)]);
f1 = figure(); 
subplot(3,1,1);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_dtrig_ctx,2),nanstd(l3mec_non_sk_dtrig_ctx,[],2)./sqrt(sum(~isnan(l3mec_non_sk_dtrig_ctx),2)));
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_dtrig_ctx,2),nanstd(l3mec_sk_dtrig_ctx,[],2)./sqrt(sum(~isnan(l3mec_sk_dtrig_ctx),2)),{'color','r'});
xlim(xr);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
title('Ctx LF');
ylabel('LF Amp (z)');

l3mec_non_sk_dtrig_hpc = reshape([trig_data(l3mec_pup_hpclfp).non_sk_dtrig_hpc_lf3],[length(lags) length(l3mec_pup_hpclfp)]);
l3mec_sk_dtrig_hpc = reshape([trig_data(l3mec_pup_hpclfp).sk_dtrig_hpc_lf3],[length(lags) length(l3mec_pup_hpclfp)]);
subplot(3,1,2);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_dtrig_hpc,2),nanstd(l3mec_non_sk_dtrig_hpc,[],2)./sqrt(sum(~isnan(l3mec_non_sk_dtrig_hpc),2)),{'color','k'});
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_dtrig_hpc,2),nanstd(l3mec_sk_dtrig_hpc,[],2)./sqrt(sum(~isnan(l3mec_sk_dtrig_hpc),2)),{'color','r'});
xlim(xr);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
title('Hpc LF');
ylabel('LF Amp (z)');

% l3mec_non_sk_dtrig_mp = reshape([trig_data(l3mec_pups).non_sk_dtrig_mp_lf],[length(lags) length(l3mec_pups)]);
% l3mec_sk_dtrig_mp = reshape([trig_data(l3mec_pups).sk_dtrig_mp_lf],[length(lags) length(l3mec_pups)]);
l3mec_non_sk_dtrig_mp = reshape([trig_data(l3mec_pups).non_sk_dtrig_mp_dc],[length(lags) length(l3mec_pups)]);
l3mec_sk_dtrig_mp = reshape([trig_data(l3mec_pups).sk_dtrig_mp_dc],[length(lags) length(l3mec_pups)]);
subplot(3,1,3);hold on
shadedErrorBar(lags/Fsd,nanmean(l3mec_non_sk_dtrig_mp,2),nanstd(l3mec_non_sk_dtrig_mp,[],2)./sqrt(sum(~isnan(l3mec_non_sk_dtrig_mp),2)),{'color','k'});
shadedErrorBar(lags/Fsd,nanmean(l3mec_sk_dtrig_mp,2),nanstd(l3mec_sk_dtrig_mp,[],2)./sqrt(sum(~isnan(l3mec_sk_dtrig_mp),2)),{'color','r'});
xlim(xr);
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
line(xr,[0 0],'color','k','linestyle','--');
title('MP LF');
ylabel('LF Amp (z)');

%%
fig_width = 3.5; rel_height = 2.4;
figufy(f1);
fname = [fig_dir 'dtrig_LFamp_5s_wch5_dcmp.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

