%%
clear all

load ~/Analysis/Mayank/final_pdown_analysis/compiled_data.mat
load ~/Analysis/Mayank/final_pdown_analysis/mua_classification_fin.mat

min_rec_dur = 500; %minimum total duration of recording (in sec)
data_ids = [data(:).id];
used_dirs = find([data(:).ep] > min_rec_dur & [data(:).dp] > min_rec_dur);
used_dirs(used_dirs == 72) = [];
% used_dirs(ismember(data_ids(used_dirs),no_cell) | ismember(data_ids(used_dirs),unclear_uds)) = [];

data = data(used_dirs);


%% parameters
raw_Fs = 2016;
dsf = 8;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
lcf = 0.05; hcf = 10; %freq range of bandpass filter for 'LF amplitude' signals

dc_dsf = 20; %down-sample factor for DC MP signal
dc_hcf = 20; %high-cut freq for DC MP signal
%range of times around detected spike peaks to interpolate out from MP sig
spk_back = 0.001; 
spk_for = 0.004;

bimod_nBoot = 100;
seg_win = 50; %duration of time windows to compute segmented bimodality 

%%
for d = 1:length(data)
%     cd(data(d).dir)
    cd(data(d).new_dir)
    pwd
    
    dist_data(d).data_id = data(d).id;
    
    load ./used_data lf7 wcv_minus_spike

    %get LF amplitude signals
    [lfp_lf,t_axis] = get_lf_features(lf7,raw_Fs,Fsd,[lcf hcf]);
    wcv_lf = get_lf_features(wcv_minus_spike,raw_Fs,Fsd,[lcf hcf]);
           
    %get time stamps
    load ./sync_times.mat
    synct_d = downsample(synct,dsf);
    
    %get MP spike times
    load ./spike_time_jmm

    %exclude any data that was out of bounds
    end_time = min(data(d).ep,data(d).dp); %final good time of the rec
    ep = find(t_axis >= end_time,1); %final index 
    if ~isempty(ep)
        synct_d(ep+1:end) = []; lfp_lf(ep+1:end) = []; wcv_lf(ep+1:end) = []; 
        t_axis(ep+1:end) = []; 
    else
        ep = length(t_axis);
    end

    %% load in UDS classified data
    load ./pa_hsmm_state_seq_combined_fin_nd.mat
    load ./pa_hsmm_state_seq7_combined_fin_nd.mat
    hsmm_ctx = hsmm7;
    lfp_state_seq = hsmm_bbstate_seq7;
    mp_state_seq = hsmm_bbstate_seq;
        
    [new_seg_inds] = resample_uds_seg_inds_v2(hsmm.UDS_segs,hsmm.Fs,Fsd,mp_state_seq); %get indices of UDS segs
    dur_uds = sum(diff(new_seg_inds,[],2))/Fsd; %total duration of UDS segs
    
    %get state transition indices
    [up_trans_inds_mp,down_trans_inds_mp] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
    [up_trans_inds_lfp,down_trans_inds_lfp] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);
    
    bad_mp_states = find(up_trans_inds_mp > ep | down_trans_inds_mp > ep);
    up_trans_inds_mp(bad_mp_states) = []; down_trans_inds_mp(bad_mp_states) = [];
    bad_lfp_states = find(up_trans_inds_lfp > ep | down_trans_inds_lfp > ep);
    up_trans_inds_lfp(bad_lfp_states) = []; down_trans_inds_lfp(bad_lfp_states) = [];
        
    %% load cortical UDS period data and get the vector of time measured in cortical UDS cycles
    load ./allEC_ctx_period_data_hsmm.mat
    lfp_period_vec = nan(size(wcv_lf));
    for i = 1:size(new_seg_inds,1)
        cur_inds = new_seg_inds(i,1):new_seg_inds(i,2);
        cur_inds_used = find(cur_inds <= ep);
        lfp_period_vec(cur_inds(cur_inds_used)) = lf8_period_f{i}(cur_inds_used);
    end
    
    %% get MP and LFP UDS state vectors
    mp_state_number = nan(length(wcv_lf),1);
    mp_state_vec = zeros(length(wcv_lf),1);
    lfp_state_number = nan(length(wcv_lf),1);
    lfp_state_vec = zeros(length(wcv_lf),1);
    for i = 1:length(up_trans_inds_mp)-1
        mp_state_vec(up_trans_inds_mp(i):down_trans_inds_mp(i)) = 1;
        mp_state_number(up_trans_inds_mp(i):down_trans_inds_mp(i)) = 2*(i-1)+1;
        mp_state_number(down_trans_inds_mp(i):up_trans_inds_mp(i+1)) = 2*(i-1) + 2;
    end
    for i = 1:length(up_trans_inds_lfp)-1
        lfp_state_vec(up_trans_inds_lfp(i):down_trans_inds_lfp(i)) = 1;
        lfp_state_number(up_trans_inds_lfp(i):down_trans_inds_lfp(i)) = 2*(i-1)+1;
        lfp_state_number(down_trans_inds_lfp(i):up_trans_inds_lfp(i+1)) = 2*(i-1) + 2;
    end
    mp_state_number(isnan(lfp_period_vec)) = nan;
    lfp_state_number(isnan(lfp_period_vec)) = nan;
    mp_state_vec(isnan(lfp_period_vec)) = nan;
    lfp_state_number(isnan(lfp_period_vec)) = nan;
               
    %% CREATE INTERPOLATED MP DATA
    if exist('./aligned_heka.mat','file')
        load ./spike_time_jmm.mat
        full_t_axis = (1:length(wcv_minus_spike))/raw_Fs;
        spike_times = full_t_axis(spkid);
        
        load ./aligned_heka.mat
        dc_dt = median(diff(dc_time));
        dc_fs = 1/dc_dt;
        dc_fsd = dc_fs/dc_dsf;
%         [dc_fb,dc_fa] = butter(4,dc_hcf/(dc_fsd/2),'low'); %make low-pass filter
                    [dc_fb,dc_fa] = butter(2,[0.01 dc_hcf]/(dc_fsd/2));
        
        %remove spikes from DC MP signal by interpolation
        cur_spk_inds = round(interp1(dc_time,1:length(dc_time),spike_times)); %indices of spikes
        spike_times(isnan(cur_spk_inds)) = [];
        cur_spk_inds(isnan(cur_spk_inds)) = [];
        cur_spk_inds = cur_spk_inds(:); spike_times = spike_times(:); dc_time = dc_time(:);
        spike_error = abs(dc_time(cur_spk_inds) - spike_times); %use only spikes that interpolate onto actual DC data
        bad = find(spike_error > 2/Fsd); %if interpolated time is too far from the real time get rid of these
        cur_spk_inds(bad) = [];
        blk_inds = -round(spk_back/dc_dt):round(spk_for/dc_dt); %set of times around each spike to interpolate out
        blocked_inds = bsxfun(@plus,cur_spk_inds,blk_inds);
        blocked_inds = blocked_inds(:);
        used_inds = setdiff(1:length(dc_data),blocked_inds); %non-spike samples
        dc_data = interp1(used_inds,dc_data(used_inds),1:length(dc_data)); %de-spiked data
        
        dc_data = decimate(dc_data,dc_dsf); %down-sample de-spiked signal
        dc_data = filtfilt(dc_fb,dc_fa,dc_data); %low-pass filter
        dc_time = downsample(dc_time,dc_dsf); %down-sample dc time axis
        dc_time = dc_time(:);
        
        dc_interp_data = interp1(dc_time,dc_data,t_axis); %interpolate dc data onto same time axis
        ind_interp = ceil(interp1(dc_time,1:length(dc_time),t_axis));
        used = find(~isnan(ind_interp));
        t_error = dc_time(ind_interp(used)) - t_axis(used)'; %find how far off interpolated time points were from dc time samples
        used = used(t_error <= 1/Fsd); %only use data that was interpolated
        dc_interp_data(setdiff(1:length(dc_interp_data),used)) = nan; %set the rest to nan
        
        if nanstd(dc_interp_data) < 1 %if using units of 1/100 V (as is case in some recs)
            dc_interp_data = dc_interp_data*100;
        end
        dc_interp_data = dc_interp_data(:);
    else
        dc_interp_data = nan(size(wcv_lf));
    end

    %% loop over data segments and check MP bimodality in each seg
    seg_durs = diff(new_seg_inds,[],2)/Fsd;
    cnt = 1;
    clear MP_bimod_stats MP_bimod_p
    for ii = 1:size(new_seg_inds,1)
       cur_Nwins =  floor(seg_durs(ii)/seg_win);
       for jj = 1:cur_Nwins
            cur_inds = new_seg_inds(ii,1) + (1:round(seg_win*Fsd)) + (jj-1)*round(seg_win*Fsd);
            cur_inds(cur_inds > length(dc_interp_data)) = [];
            cur_dc_data = downsample(dc_interp_data(cur_inds),5);
            cur_dc_data(isnan(cur_dc_data)) = [];
            if ~isnan(cur_dc_data)
            [MP_bimod_stats(cnt),MP_bimod_p(cnt)] = hartigansdipsigniftest(cur_dc_data,bimod_nBoot);
            cnt = cnt + 1;
            end
       end
    end
    
    dist_data(d).bimod_stats = MP_bimod_stats;
    dist_data(d).bimod_p = MP_bimod_p;
    
    %check overall bimodality
    cur_dc_data = downsample(dc_interp_data(~isnan(mp_state_vec)),10);
    cur_dc_data(isnan(cur_dc_data)) = [];
    [dist_data(d).ov_bimod_stats,dist_data(d).ov_bimod_p] = hartigansdipsigniftest(cur_dc_data,bimod_nBoot);
    
    dist_data(d).tot_UDS_dur = sum(seg_durs);
    
end

%%
tot_UDS_durs = [dist_data(:).tot_UDS_dur];
temp = arrayfun(@(x) sum(x.bimod_p < 0.05),dist_data);

cd ~/Analysis/Mayank/final_pdown_analysis/

