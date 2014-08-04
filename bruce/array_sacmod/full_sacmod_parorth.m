clear all
% close all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';

Expt_nums =  [86 88 89];

%% SET ANALYSIS PARAMETERS
stim_fs = 100; %in Hz

dt = 0.005;

min_trial_dur = 0.5;
beg_buffer = 0.3;

backlag = round(0.25/dt);
forwardlag = round(0.55/dt);

sua_sm_sig = (0.01/dt);
mua_sm_sig = (0.005/dt);

%% LOAD OVERALL SU DATA
load ~/Analysis/bruce/summary_analysis/su_data.mat

%% LOOP OVER EXPTS
for ex = 1:length(Expt_nums)
    % for ex = 3
    clear expt* included_type
    
    Expt_name = sprintf('G0%d',Expt_nums(ex));
    data_dir = [dir_prefix '/Data/bruce/' Expt_name];
    addpath([dir_prefix '/James_scripts/bruce/G081/'])
    
    cd(data_dir);
    
    load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
    load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data
    
    save_dir = [dir_prefix '/Analysis/bruce/' Expt_name '/sac_mod'];
    if ~exist(save_dir,'dir')
        system(['mkdir ' save_dir]);
    end
    
    if exist('./fin_aclust_data.mat','file')
        load('./fin_aclust_data.mat');
        [n_aclust_expts,n_aclust_probes] = size(autoclust);
    else
        disp('No fin_aclust_data found.');
        autoclust = [];
        n_aclust_expts = 0; n_aclust_probes = 0;
    end
    
    trial_dur = 4;
    
    %%
    include_expts = {'rls.seXFaRC','rls.seXFaXFrRC','rls.orXFaXsl'};
    for ii = 1:length(Expts)
        expt_names{ii} = Expts{ii}.Header.expname;
        expt_dds(ii) = Expts{ii}.Stimvals.dd;
        expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
        expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
        expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
        expt_imback(ii) = isfield(Expts{ii}.Trials,'imi');
        included_type(ii) = any(strcmp(expt_names{ii},include_expts));
        expt_ntrials(ii) = length(Expts{ii}.Trials);
    end
    expt_has_ds(isnan(expt_has_ds)) = 0;
    expt_has_ds(expt_has_ds == -1) = 0;
    expt_binoc(isnan(expt_binoc)) = 0;
    cur_block_set = find(included_type & expt_ntrials' >= 10);
    
    if strcmp(Expt_name,'G088')
        cur_block_set(cur_block_set == 22) = [];
    end
    
    %% load overall su data
    cur_expt_id = find(su_data.expt_nums == Expt_nums(ex));
    su_probes = find(su_data.is_su(cur_expt_id,:));
    mua_probes = setdiff(1:96,su_probes); %probes with ONLY MU
    aclust_probenums = [autoclust(cur_block_set(1),:).probenum];
    autoclust = autoclust(:,ismember(aclust_probenums,su_probes));
    
    %% COMPUTE TRIAL DATA
    cd(data_dir);
    
    fprintf('Computing prep data\n');
    trial_cnt = 0;
    
    all_stim_times = [];
    all_t_axis = [];
    all_t_bin_edges = [];
    all_tsince_start = [];
    all_blockvec = [];
    all_trialvec = [];
    all_trial_wi = [];
    all_frame_dur = [];
    all_bar_or = [];
    all_sac_dir = [];
    all_trial_start_times = [];
    all_trial_end_times = [];
    all_bin_edge_pts = [];
    all_spk_times = cell(96,1);
    all_clust_ids = cell(length(su_probes),1);
    for ee = 1:length(cur_block_set);
        cur_block = cur_block_set(ee);
        fname = sprintf('Expt%dClusterTimesDetails.mat',cur_block);
        load(fname);
        for cc = 1:96
            all_spk_times{cc} = cat(1,all_spk_times{cc},ClusterDetails{cc}.t');
            cur_ind = find(su_probes == cc);
            if ~isempty(cur_ind)
                all_clust_ids{cur_ind} = cat(1,all_clust_ids{cur_ind},autoclust(cur_block_set(ee),cur_ind).idx(:));
            end
        end
        
        trial_start_times = [Expts{cur_block}.Trials(:).TrialStart]/1e4;
        trial_end_times = [Expts{cur_block}.Trials(:).TrueEnd]/1e4;
        trial_durs = [Expts{cur_block}.Trials(:).dur]/1e4;
        trial_ids = [Expts{cur_block}.Trials(:).id];
        
        [un_ids,id_inds] = unique(trial_ids);
        if length(un_ids) < length(trial_ids)
            fprintf('Warning, repeat trial inds detected!\n');
        end
        
        use_trials = find(trial_durs >= min_trial_dur);
        
        all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
        all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)');
                
        if isfield(Expts{cur_block}.Trials(1),'Fr')
            trial_Fr = [Expts{cur_block}.Trials(:).Fr];
            all_frame_dur = [all_frame_dur; trial_Fr(use_trials)'];
        elseif isfield(Expts{cur_block}.Trials(1),'sl')
            trial_sl = [Expts{cur_block}.Trials(:).sl];
            trial_sl(trial_sl == 0) = 1;
            all_frame_dur = [all_frame_dur; trial_sl(use_trials)'];
        else
            all_frame_dur = [all_frame_dur; ones(length(use_trials),1)*expt_Fr(cur_block)];
        end
        if isfield(Expts{cur_block}.Trials(1),'or')
            trial_or = [Expts{cur_block}.Trials(:).or];
            all_bar_or = [all_bar_or; trial_or(use_trials)'];
        else
            all_bar_or = [all_bar_or; ones(length(use_trials),1)*expt_bar_ori(cur_block)];
        end
        if isfield(Expts{cur_block}.Trials(1),'Fa')
            trial_Fa = [Expts{cur_block}.Trials(:).Fa];
            all_sac_dir = [all_sac_dir; trial_Fa(use_trials)'];
        else
            all_sac_dir = [all_sac_dir; ones(length(use_trials),1)*expt_sac_dir(cur_block)];
        end

        n_trials = length(use_trials);
        for tt = 1:n_trials
            cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start'/1e4;
            if length(cur_stim_times) == 1
                cur_stim_times = trial_start_times(use_trials(tt)):1/stim_fs:trial_end_times(use_trials(tt));
            end
            cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
            cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
            
            cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
                        
            all_stim_times = [all_stim_times; cur_stim_times'];
            all_t_axis = [all_t_axis; cur_t_axis];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
        trial_cnt = trial_cnt + n_trials;
    end
    all_sac_dir = mod(all_sac_dir,180);
    
    %%
    mahal_thresh = su_data.mah_thresh;
    all_binned_spikes = nan(length(all_t_axis),96);
    su_used_blocks = false(length(cur_block_set),length(su_probes));
    %for only-MU probes
    for cc = 1:96
        if ~ismember(cc,su_probes) %if probe doesn't have an SU
            [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc},all_t_bin_edges);
            cur_spkhist(all_bin_edge_pts) = [];
            all_binned_spikes(:,cc) = cur_spkhist;
        end
    end
    %for SU probes
    for ss = 1:length(su_probes)
        %for MUA
        cur_muahist = histc(all_spk_times{su_probes(ss)},all_t_bin_edges);
        cur_muahist(all_bin_edge_pts) = [];
        
        all_su_inds = find(all_clust_ids{ss} == 1);
        cur_suahist = histc(all_spk_times{su_probes(ss)}(all_su_inds),all_t_bin_edges);
        cur_suahist(all_bin_edge_pts) = [];
        
        spk_block_inds = round(interp1(all_t_axis,all_blockvec,all_spk_times{su_probes(ss)}));
        for ee = 1:length(cur_block_set);
            cur_block_inds = find(all_blockvec == ee);
            %if SU is isolated in this block
            if autoclust(cur_block_set(ee),ss).mahal_d > mahal_thresh || autoclust(cur_block_set(ee),ss).man_code == 4
                su_used_blocks(ee,ss) = true;
                all_binned_spikes(cur_block_inds,su_probes(ss)) = cur_suahist(cur_block_inds);
            else %otherwise use MUA
                all_binned_spikes(cur_block_inds,su_probes(ss)) = cur_muahist(cur_block_inds);
            end
        end
    end
    %% DEFINE DATA USED FOR ANALYSIS
    used_inds = find(all_tsince_start >= beg_buffer);
    %% SMOOTH AND NORMALIZE SPIKING DATA
    all_spike_rate = nan(size(all_binned_spikes));
    for cc = 1:96
        all_spike_rate(:,cc) = jmm_smooth_1d_cor(all_binned_spikes(:,cc),mua_sm_sig);
    end
    all_spike_rate_norm = nan(size(all_binned_spikes));
    block_mean_rates = nan(length(cur_block_set),96);
    block_n_spikes = nan(length(cur_block_set),96);
    for ee = 1:length(cur_block_set)
        cur_block_inds = used_inds(all_blockvec(used_inds) == ee);
        block_mean_rates(ee,:) = mean(all_spike_rate(cur_block_inds,:));
        block_n_spikes(ee,:) = sum(all_binned_spikes(cur_block_inds,:));
    end
    
    mua_avg_rates = nan(96,1);
    sua_avg_rates = nan(length(su_probes),1);
    mua_tot_nspikes = nan(96,1);
    sua_tot_nspikes = nan(96,1);
    
    mua_avg_rates(mua_probes) = mean(block_mean_rates(:,mua_probes));
    mua_tot_nspikes(mua_probes) = sum(block_n_spikes(:,mua_probes));
    for ss = 1:length(su_probes)
        sua_avg_rates(ss) = mean(block_mean_rates(su_used_blocks(:,ss),su_probes(ss)));
        sua_tot_nspikes(ss) = sum(block_n_spikes(su_used_blocks(:,ss),su_probes(ss)));
        mu_blocks = find(su_used_blocks(:,ss) == 0);
        if ~isempty(mu_blocks)
            mua_avg_rates(su_probes(ss)) = mean(block_mean_rates(mu_blocks,su_probes(ss)));
            mua_tot_nspikes(su_probes(ss)) = sum(block_n_spikes(mu_blocks,su_probes(ss)));
            for bb = 1:length(mu_blocks)
                cur_block_inds = used_inds(all_blockvec(used_inds) == mu_blocks(bb));
                all_spike_rate(cur_block_inds,su_probes(ss)) = jmm_smooth_1d_cor(all_binned_spikes(cur_block_inds,su_probes(ss)),sua_sm_sig);
            end
        end
    end
    %normalized firing rates (smoothed)
    for ee = 1:length(cur_block_set)
        cur_block_inds = used_inds(all_blockvec(used_inds) == ee);
        all_spike_rate_norm(cur_block_inds,:) = bsxfun(@rdivide,all_spike_rate(cur_block_inds,:),...
            block_mean_rates(ee,:));
    end
    
    %% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)
    
    emfile = ['jbe' Expt_name '.em.mat'];
    load(emfile);
    
    all_eye_vals = [];
    all_eye_speed = [];
    all_eye_ts = [];
    all_eye_blockvec = [];
    eye_smooth = 3;
    for ee = 1:length(cur_block_set);
        fprintf('Loading ET data for expt %d, block %d of %d\n',Expt_nums(ex),ee,length(cur_block_set));
        cur_set = find(all_blockvec==ee);
        [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data_new(Expt,all_t_axis(cur_set([1 end])));
        
        eye_dt = Expt.Header.CRrates(1);
        eye_fs = 1/eye_dt;
        lEyeXY = eye_vals_interp(:,1:2);
        rEyeXY = eye_vals_interp(:,3:4);
        
        %slight smoothing before computing speed
        sm_avg_eyepos = lEyeXY; eye_vel = lEyeXY; %initialization
        sm_avg_eyepos(:,1) = smooth(lEyeXY(:,1),eye_smooth);
        sm_avg_eyepos(:,2) = smooth(lEyeXY(:,2),eye_smooth);
        eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/eye_dt;
        eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/eye_dt;
        
        eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
        
        all_eye_vals = [all_eye_vals; lEyeXY rEyeXY];
        all_eye_speed = [all_eye_speed; eye_speed];
        all_eye_ts = [all_eye_ts; eye_ts_interp'];
        all_eye_blockvec = [all_eye_blockvec; ee*ones(size(eye_speed))];
    end
    
    back_pts = 1 + find(diff(all_eye_ts) <= 0);
    double_samples = [];
    for i = 1:length(back_pts)
        next_forward = find(all_eye_ts > all_eye_ts(back_pts(i)-1),1,'first');
        double_samples = [double_samples back_pts(i):next_forward];
    end
    all_eye_ts(double_samples) = [];
    all_eye_speed(double_samples) = [];
    all_eye_vals(double_samples,:) = [];
    
    interp_eye_vals = interp1(all_eye_ts,all_eye_vals,all_t_axis);
    interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);
    
    %%
    sac_thresh = 10;
    peak_sig = [0; diff(sign(diff(all_eye_speed))); 0];
    saccade_inds = find(peak_sig == -2 & all_eye_speed > sac_thresh);
    
    peri_thresh = 3; %threshold eye speed for defining saccade boundary inds
    thresh_cross_up = 1 + find(all_eye_speed(1:end-1) < peri_thresh & all_eye_speed(2:end) >= peri_thresh);
    thresh_cross_down = 1 + find(all_eye_speed(1:end-1) >= peri_thresh & all_eye_speed(2:end) < peri_thresh);
    sac_start_inds = nan(size(saccade_inds));
    sac_stop_inds = nan(size(saccade_inds));
    for ii = 1:length(saccade_inds)
        next_tc = find(thresh_cross_down > saccade_inds(ii),1,'first');
        if ~isempty(next_tc)
            sac_stop_inds(ii) = thresh_cross_down(next_tc);
        end
        prev_tc = find(thresh_cross_up < saccade_inds(ii),1,'last');
        if ~isempty(prev_tc)
            sac_start_inds(ii) = thresh_cross_up(prev_tc);
        end
        
    end
    
    %get rid of double-peaks
    min_isi = 0.05; max_isi = Inf;
    isis = [Inf; diff(sac_start_inds)]/eye_fs;
    bad_isis = (isis < min_isi | isis > max_isi);
    bad_sacs = find(isnan(sac_stop_inds) | isnan(sac_start_inds) | bad_isis);
    saccade_inds(bad_sacs) = []; isis(bad_sacs) = []; sac_start_inds(bad_sacs) = []; sac_stop_inds(bad_sacs) = [];
    
    saccade_times = all_eye_ts(saccade_inds);
    sac_start_times = all_eye_ts(sac_start_inds);
    sac_stop_times = all_eye_ts(sac_stop_inds);
    sac_durs = sac_stop_times - sac_start_times;
    
    sac_dbuff = round(0.005/eye_dt);
    pre_inds = saccade_inds - sac_dbuff;
    pre_inds(pre_inds < 1) = 1;
    sac_pre_pos = all_eye_vals(pre_inds,:);
    post_inds = saccade_inds + sac_dbuff;
    post_inds(post_inds > length(all_eye_ts)) = length(all_eye_ts);
    sac_post_pos = all_eye_vals(post_inds,:);
    
    %use only left-eye signal here
    sac_delta_pos = sac_post_pos(:,1:2) - sac_pre_pos(:,1:2);
    sac_amps = sqrt(sum(sac_delta_pos.^2,2));
    sac_dirs = atan2(sac_delta_pos(:,2),sac_delta_pos(:,1));
    
    temp = ones(length(saccade_times),1);
    saccades = struct('peak_time',mat2cell(saccade_times,temp),'start_time',mat2cell(sac_start_times,temp),...
        'stop_time',mat2cell(sac_stop_times,temp),'isi',mat2cell(isis,temp),...
        'duration',mat2cell(sac_durs,temp),'amplitude',mat2cell(sac_amps,temp),'direction',mat2cell(sac_dirs,temp),...
        'pre_Lx',mat2cell(sac_pre_pos(:,1),temp),'post_Lx',mat2cell(sac_post_pos(:,1),temp),...
        'pre_Ly',mat2cell(sac_pre_pos(:,2),temp),'post_Ly',mat2cell(sac_post_pos(:,2),temp),...
        'pre_Rx',mat2cell(sac_pre_pos(:,3),temp),'post_Rx',mat2cell(sac_post_pos(:,3),temp),...
        'pre_Ry',mat2cell(sac_pre_pos(:,4),temp),'post_Ry',mat2cell(sac_post_pos(:,4),temp));
    
    
    %% PICK OUT SACCADES FOR ANALYSIS
    sac_start_times = [saccades(:).start_time];
    interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
    interp_sac_start_inds(isnan(interp_sac_start_inds)) = 1;
    sac_error = abs(sac_start_times - all_t_axis(interp_sac_start_inds)');
    bad_sacs = find(isnan(interp_sac_start_inds) | sac_error > dt);
    sac_start_times(bad_sacs) = [];
    saccades(bad_sacs) = [];
    interp_sac_start_inds(bad_sacs) = [];
    
    sac_stop_times = [saccades(:).stop_time];
    interp_sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_stop_times));
    interp_sac_stop_inds(isnan(interp_sac_stop_inds)) = 1;
    
    sac_peak_times = [saccades(:).peak_time];
    interp_sac_peak_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_peak_times));
    interp_sac_peak_inds(isnan(interp_sac_peak_inds)) = 1;
    
    max_dur = 0.1;
    sac_durs = [saccades(:).duration];
    is_blink = sac_durs > max_dur;
    
    isis = [saccades(:).isi];
    is_sacburst = false(length(saccades),1);
    is_sacburst(isis(1:end-1) < 0.15 | isis(2:end) < 0.15) = true;
    
    sac_amps = [saccades(:).amplitude];
    is_micro = sac_amps < 1;
    
    sac_post_Lx = [saccades(:).post_Lx];
    sac_post_Ly = [saccades(:).post_Ly];
    sac_pre_Lx = [saccades(:).pre_Lx];
    sac_pre_Ly = [saccades(:).pre_Ly];
    sac_deltaX = sac_post_Lx - sac_pre_Lx;
    sac_deltaY = sac_post_Ly - sac_pre_Ly;
    sac_dirs = [0 90];
    sim_sac_times = [0.7 1.4 2.1 2.8 3.5];
    sac_thresh = 1;
    
    sac_start_expts = all_blockvec(interp_sac_start_inds);
    delta_sacpar = nan(size(saccades));
    for bb = 1:length(sac_dirs)
        cur_bar_expts = find(expt_sac_dir(cur_block_set) == sac_dirs(bb));
        cur_sac_set = find(ismember(sac_start_expts,cur_bar_expts));
        
        cur_delta_sac_par = abs(sac_deltaX(cur_sac_set)*cos(sac_dirs(bb)*pi/180) + sac_deltaY(cur_sac_set)*sin(sac_dirs(bb)*pi/180));
        delta_sacpar(cur_sac_set) = cur_delta_sac_par;
    end
    is_gsac = delta_sacpar' > sac_thresh;
        
    %define which saccades to use
    used_msacs = find(is_micro & ~is_blink);
    used_gsacs = find(is_gsac & ~is_blink);
    
    gsac_start_inds = interp_sac_start_inds(used_gsacs);
    msac_start_inds = interp_sac_start_inds(used_msacs);
    
    trial_start_inds = find(all_tsince_start < dt);
    
    p100_trials = find(all_frame_dur == 1 & all_sac_dir == all_bar_or);
    p30_trials = find(all_frame_dur == 3 & all_sac_dir == all_bar_or);
    o100_trials = find(all_frame_dur == 1 & all_sac_dir ~= all_bar_or);
    o30_trials = find(all_frame_dur == 3 & all_sac_dir ~= all_bar_or);
    
    p100_gsacs = gsac_start_inds(ismember(all_trialvec(gsac_start_inds),p100_trials));
    p100_msacs = msac_start_inds(ismember(all_trialvec(msac_start_inds),p100_trials));
    p30_gsacs = gsac_start_inds(ismember(all_trialvec(gsac_start_inds),p30_trials));
    p30_msacs = msac_start_inds(ismember(all_trialvec(msac_start_inds),p30_trials));
    o100_gsacs = gsac_start_inds(ismember(all_trialvec(gsac_start_inds),o100_trials));
    o100_msacs = msac_start_inds(ismember(all_trialvec(msac_start_inds),o100_trials));
    o30_gsacs = gsac_start_inds(ismember(all_trialvec(gsac_start_inds),o30_trials));
    o30_msacs = msac_start_inds(ismember(all_trialvec(msac_start_inds),o30_trials));
      
    p100_inds = find(ismember(all_trialvec,p100_trials));
    p30_inds = find(ismember(all_trialvec,p30_trials));
    o100_inds = find(ismember(all_trialvec,o100_trials));
    o30_inds = find(ismember(all_trialvec,o30_trials));
    
    %% COMPUTE TRIG AVGS FOR MUA
    nboot = [];
    clear mua_data
    for cc = 1:96
        
        fprintf('Computing trig avgs for MUA %d of %d\n',cc,96);        
        su_probe_ind = find(su_probes == cc);
        if ~isempty(su_probe_ind)
            cur_used_blocks = find(su_used_blocks(:,su_probe_ind) == 0); %blocks when NO SU
            cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
        else
            cur_poss_inds = 1:length(all_trialvec);
        end
        cur_used_inds = used_inds(ismember(used_inds,cur_poss_inds));
        
        %general averages
        [mua_data(cc).p100_gsac_avg,lags] = get_event_trig_avg(all_spike_rate_norm(:,cc),p100_gsacs,backlag,forwardlag,nboot,all_trialvec,0,cur_used_inds);
        [mua_data(cc).p100_msac_avg,lags] = get_event_trig_avg(all_spike_rate_norm(:,cc),p100_msacs,backlag,forwardlag,nboot,all_trialvec,0,cur_used_inds);
        [mua_data(cc).p30_gsac_avg,lags] = get_event_trig_avg(all_spike_rate_norm(:,cc),p30_gsacs,backlag,forwardlag,nboot,all_trialvec,0,cur_used_inds);
        [mua_data(cc).p30_msac_avg,lags] = get_event_trig_avg(all_spike_rate_norm(:,cc),p30_msacs,backlag,forwardlag,nboot,all_trialvec,0,cur_used_inds);
        [mua_data(cc).o100_gsac_avg,lags] = get_event_trig_avg(all_spike_rate_norm(:,cc),o100_gsacs,backlag,forwardlag,nboot,all_trialvec,0,cur_used_inds);
        [mua_data(cc).o100_msac_avg,lags] = get_event_trig_avg(all_spike_rate_norm(:,cc),o100_msacs,backlag,forwardlag,nboot,all_trialvec,0,cur_used_inds);
        [mua_data(cc).o30_gsac_avg,lags] = get_event_trig_avg(all_spike_rate_norm(:,cc),o30_gsacs,backlag,forwardlag,nboot,all_trialvec,0,cur_used_inds);
        [mua_data(cc).o30_msac_avg,lags] = get_event_trig_avg(all_spike_rate_norm(:,cc),o30_msacs,backlag,forwardlag,nboot,all_trialvec,0,cur_used_inds);
             
        cur_p100_inds = intersect(p100_inds,cur_used_inds);
        cur_p30_inds = intersect(p30_inds,cur_used_inds);
        cur_o100_inds = intersect(o100_inds,cur_used_inds);
        cur_o30_inds = intersect(o30_inds,cur_used_inds);
        p100_avg = mean(all_spike_rate_norm(cur_p100_inds,cc));
        p30_avg = mean(all_spike_rate_norm(cur_p30_inds,cc));
        o100_avg = mean(all_spike_rate_norm(cur_o100_inds,cc));
        o30_avg = mean(all_spike_rate_norm(cur_o30_inds,cc));
        mua_data(cc).cond_avg = [p100_avg p30_avg o100_avg o30_avg];
        mua_data(cc).avg_rate = mean(all_binned_spikes(cur_used_inds,cc));
        mua_data(cc).tot_nspikes = sum(all_binned_spikes(cur_used_inds,cc));
    end
    
    %% COMPUTE TRIG AVGS FOR SUA
    nboot = 100;
    %set trial numbers to Inf so they don't get included in trig averaging
    clear sua_data
    for ss = 1:length(su_probes)
        
        fprintf('Computing trig avgs for SU %d of %d\n',ss,length(su_probes));
        
        cur_used_blocks = find(su_used_blocks(:,ss)); %blocks when SU
        cur_used_inds = used_inds(ismember(all_blockvec(used_inds),cur_used_blocks));
        
        %general averages
        [sua_data(ss).p100_gsac_avg,lags] = get_event_trig_avg(all_spike_rate_norm(:,su_probes(ss)),p100_gsacs,backlag,forwardlag,nboot,all_trialvec,0,cur_used_inds);
        [sua_data(ss).p100_msac_avg,lags] = get_event_trig_avg(all_spike_rate_norm(:,su_probes(ss)),p100_msacs,backlag,forwardlag,nboot,all_trialvec,0,cur_used_inds);
        [sua_data(ss).p30_gsac_avg,lags] = get_event_trig_avg(all_spike_rate_norm(:,su_probes(ss)),p30_gsacs,backlag,forwardlag,nboot,all_trialvec,0,cur_used_inds);
        [sua_data(ss).p30_msac_avg,lags] = get_event_trig_avg(all_spike_rate_norm(:,su_probes(ss)),p30_msacs,backlag,forwardlag,nboot,all_trialvec,0,cur_used_inds);
        [sua_data(ss).o100_gsac_avg,lags] = get_event_trig_avg(all_spike_rate_norm(:,su_probes(ss)),o100_gsacs,backlag,forwardlag,nboot,all_trialvec,0,cur_used_inds);
        [sua_data(ss).o100_msac_avg,lags] = get_event_trig_avg(all_spike_rate_norm(:,su_probes(ss)),o100_msacs,backlag,forwardlag,nboot,all_trialvec,0,cur_used_inds);
        [sua_data(ss).o30_gsac_avg,lags] = get_event_trig_avg(all_spike_rate_norm(:,su_probes(ss)),o30_gsacs,backlag,forwardlag,nboot,all_trialvec,0,cur_used_inds);
        [sua_data(ss).o30_msac_avg,lags] = get_event_trig_avg(all_spike_rate_norm(:,su_probes(ss)),o30_msacs,backlag,forwardlag,nboot,all_trialvec,0,cur_used_inds);
        
        cur_p100_inds = intersect(p100_inds,cur_used_inds);
        cur_p30_inds = intersect(p30_inds,cur_used_inds);
        cur_o100_inds = intersect(o100_inds,cur_used_inds);
        cur_o30_inds = intersect(o30_inds,cur_used_inds);
        p100_avg = mean(all_spike_rate_norm(cur_p100_inds,cc));
        p30_avg = mean(all_spike_rate_norm(cur_p30_inds,cc));
        o100_avg = mean(all_spike_rate_norm(cur_o100_inds,cc));
        o30_avg = mean(all_spike_rate_norm(cur_o30_inds,cc));

         sua_data(ss).cond_avg = [p100_avg p30_avg o100_avg o30_avg];
       sua_data(ss).avg_rate = mean(all_binned_spikes(cur_used_inds,su_probes(ss)));
        sua_data(ss).tot_nspikes = sum(all_binned_spikes(cur_used_inds,su_probes(ss)));
        %start and stop trig avgs
        %micro sac burst and non-burst
    end
    
    %%
    
    anal_params.dt = dt;
    anal_params.mua_smooth_sig = mua_sm_sig*dt;
    anal_params.sua_smooth_sig = sua_sm_sig*dt;
    anal_params.backlag = backlag*dt;
    anal_params.forwardlag = forwardlag*dt;
    anal_params.beg_buffer = beg_buffer;
    anal_params.min_trial_dur = min_trial_dur;
    anal_params.lags = lags;
    cd(save_dir)
    sname = 'full_parorth_data';
    save(sname,'sua_data','su_probes','mua_data','anal_params');
    clear sua_data su_probes mua_data 
end