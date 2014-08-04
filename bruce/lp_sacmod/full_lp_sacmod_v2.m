clear all
% close all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';

Expt_nums = [266];

%% SET ANALYSIS PARAMETERS
stim_fs = 100; %in Hz

dt = 0.005;

min_trial_dur = 0.5;
beg_buffer = 0.3;

backlag = round(0.25/dt);
forwardlag = round(0.55/dt);

sua_sm_sig = (0.01/dt);
mua_sm_sig = (0.005/dt);

n_probes = 24;
trial_dur = 4;

%% LOOP OVER EXPTS
for ex = 1:length(Expt_nums)
    % for ex = 8
    clear expt* included_type
    
    Expt_name = sprintf('M%d',Expt_nums(ex));
    data_dir = [dir_prefix '/Data/bruce/' Expt_name];
    
    cd(data_dir);
    
    load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct
    
    save_dir = [dir_prefix '/Analysis/bruce/' Expt_name '/sac_mod'];
    if ~exist(save_dir,'dir')
        system(['mkdir ' save_dir]);
    end
    
    % LOAD REFCLUSTERS
    cluster_dir = [dir_prefix '/Analysis/bruce/' Expt_name '/clustering'];
    fname = [cluster_dir '/final_cluster.mat'];
    if exist(fname,'file')
        load(fname);
        SU_numbers = unique(SU_ID_mat(~isnan(SU_ID_mat)));
        for ii = 1:length(SU_numbers)
            SU_tot_nblocks(ii) = sum(SU_ID_mat(:) == SU_numbers(ii));
        end
        fprintf('%d SUs Clustered\n',length(SU_numbers));
        
    else
        disp('No Cluster data found.');
    end
    
    %%
    include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaXFaXFs'};
    for ii = 1:length(Expts)
        if ~isempty(Expts{ii})
            expt_names{ii} = Expts{ii}.Header.expname;
            expt_dds(ii) = Expts{ii}.Stimvals.dd;
            expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
            expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
            expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
            expt_imback(ii) = isfield(Expts{ii}.Trials,'imi');
            included_type(ii) = any(strcmp(expt_names{ii},include_expts));
        else
            included_type(ii) = 0;
        end
    end
    expt_binoc = zeros(size(expt_bar_ori));
    expt_binoc(isnan(expt_binoc)) = 0;
    cur_block_set = find(included_type & ~expt_binoc & expt_Fr == 1);
    
    if strcmp(Expt_name,'M270')
        cur_block_set(cur_block_set == 5) = [];
    end
    
    %     sim_sac_expts = find(~expt_has_ds(cur_block_set));
    %     imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
    %     grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');
    
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
    all_trial_Fs = [];
    all_trial_wi = [];
    all_trial_start_times = [];
    all_trial_end_times = [];
    all_bin_edge_pts = [];
    all_spk_times = cell(n_probes,1);
    all_clust_ids = cell(n_probes,1);
    trial_toffset = zeros(length(cur_block_set),1);
    cur_toffset = 0;
    for ee = 1:length(cur_block_set);
        cur_block = cur_block_set(ee);
        
        fname = [cluster_dir sprintf('/Block%d_Clusters.mat',cur_block)];
        load(fname,'Clusters');
        for cc = 1:n_probes
            all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times + cur_toffset);
            all_clust_ids{cc} = cat(1,all_clust_ids{cc},Clusters{cc}.spike_clusts);
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
        
        all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
        all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
        
        if isfield(Expts{cur_block}.Trials(1),'Fs')
            trial_Fs = [Expts{cur_block}.Trials(:).Fs];
        else
            trial_Fs = nan(1,length(trial_durs));
        end
        %         trial_fs = trial_fs(id_inds);
        all_trial_Fs = cat(1,all_trial_Fs,trial_Fs(use_trials)');
        
        n_trials = length(use_trials);
        for tt = 1:n_trials
            cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start'/1e4;
            if length(cur_stim_times) == 1
                cur_stim_times = trial_start_times(use_trials(tt)):1/stim_fs:trial_end_times(use_trials(tt));
            end
            cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
            cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
            
            cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
            
            all_stim_times = [all_stim_times; cur_stim_times' + cur_toffset];
            all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
        trial_cnt = trial_cnt + n_trials;
        trial_toffset(ee) = all_t_bin_edges(end);
        cur_toffset = trial_toffset(ee);
    end
    
    %%
    all_binned_mua = nan(length(all_t_axis),n_probes);
    %for only-MU probes
    for cc = 1:n_probes
        cur_mua_inds = all_clust_ids{cc} > 0;
        [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc}(cur_mua_inds),all_t_bin_edges);
        cur_spkhist(all_bin_edge_pts) = [];
        all_binned_mua(:,cc) = cur_spkhist;
    end
    
    %for SU probes
    fprintf('Using %d SUs\n',length(SU_numbers));
    all_binned_sua = nan(length(all_t_axis),length(SU_numbers));
    cur_used_blocks = find(ismember(SU_target_blocks,cur_block_set));
    SU_ID_mat = SU_ID_mat(cur_used_blocks,:);
    [CC,BB] = meshgrid(1:length(SU_clust_data),1:length(cur_used_blocks));
    
    for ss = 1:length(SU_numbers)
        used_clust_set = unique(CC(SU_ID_mat==ss)); %set of clusters used to capture this SU
        SU_block_probes(ss,:) = nan(1,length(cur_block_set));
        for cc = 1:length(used_clust_set)
            cur_clust = used_clust_set(cc);
            cur_probe = SU_clust_data(cur_clust).probe_num;
            cur_clust_label = SU_clust_data(cur_clust).cluster_label;
            cur_blocks = find(SU_ID_mat(:,cur_clust) == ss);
            SU_block_probes(ss,cur_blocks) = cur_probe;
            
            all_su_inds = all_clust_ids{cur_probe} == cur_clust_label;
            all_su_spk_times = all_spk_times{cur_probe}(all_su_inds);
            spk_block_inds = round(interp1(all_t_axis,all_blockvec,all_su_spk_times));
            all_su_spk_times = all_su_spk_times(ismember(spk_block_inds,cur_blocks));
            
            cur_suahist = histc(all_su_spk_times,all_t_bin_edges);
            cur_suahist(all_bin_edge_pts) = [];
            cur_id_set = ismember(all_blockvec,cur_blocks);
            all_binned_sua(cur_id_set,ss) = cur_suahist(cur_id_set);
        end
    end
    
    %% DEFINE DATA USED FOR ANALYSIS
    used_inds = find(all_tsince_start >= beg_buffer);
    
    %% SMOOTH AND NORMALIZE SPIKING DATA
    all_mua_rate = nan(size(all_binned_mua));
    for cc = 1:24
        all_mua_rate(:,cc) = jmm_smooth_1d_cor(all_binned_mua(:,cc),mua_sm_sig);
    end
    all_mua_rate_norm = nan(size(all_binned_mua));
    all_sua_rate = nan(size(all_binned_sua));
    for cc = 1:length(SU_numbers)
        all_sua_rate(:,cc) = jmm_smooth_1d_cor(all_binned_sua(:,cc),sua_sm_sig);
    end
    all_sua_rate_norm = nan(size(all_binned_sua));
    
    block_mean_muarates = nan(length(cur_block_set),n_probes);
    block_n_muaspikes = nan(length(cur_block_set),n_probes);
    for ee = 1:length(cur_block_set)
        cur_block_inds = used_inds(all_blockvec(used_inds) == ee);
        block_mean_muarates(ee,:) = mean(all_mua_rate(cur_block_inds,:));
        block_n_muaspikes(ee,:) = sum(all_binned_mua(cur_block_inds,:));
    end
    block_mean_suarates = nan(length(cur_block_set),length(SU_numbers));
    block_n_suaspikes = nan(length(cur_block_set),length(SU_numbers));
    for ee = 1:length(cur_block_set)
        cur_block_inds = used_inds(all_blockvec(used_inds) == ee);
        block_mean_suarates(ee,:) = mean(all_sua_rate(cur_block_inds,:));
        block_n_suaspikes(ee,:) = sum(all_binned_sua(cur_block_inds,:));
    end
    
    mua_avg_rates = nanmean(block_mean_muarates);
    mua_tot_nspikes = nansum(block_n_muaspikes);
    sua_avg_rates = nanmean(block_mean_suarates);
    sua_tot_nspikes = nansum(block_n_suaspikes);
    
    %normalized firing rates (smoothed)
    for ee = 1:length(cur_block_set)
        cur_block_inds = used_inds(all_blockvec(used_inds) == ee);
        all_mua_rate_norm(cur_block_inds,:) = bsxfun(@rdivide,all_mua_rate(cur_block_inds,:),...
            block_mean_muarates(ee,:));
        all_sua_rate_norm(cur_block_inds,:) = bsxfun(@rdivide,all_sua_rate(cur_block_inds,:),...
            block_mean_suarates(ee,:));
    end
    
    %% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)
    
    all_eye_vals = [];
    all_eye_speed = [];
    all_eye_ts = [];
    all_eye_blockvec = [];
    eye_smooth = 3;
    for ee = 1:length(cur_block_set);
        
        emfile = sprintf('lem%s.%d.em.mat',Expt_name,cur_block_set(ee));
        load(emfile);
        
        fprintf('Loading ET data for expt %d, block %d of %d\n',Expt_nums(ex),ee,length(cur_block_set));
        cur_set = find(all_blockvec==ee);
        if ee > 1
            cur_toffset = trial_toffset(ee-1);
        else
            cur_toffset = 0;
        end
        [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data_new(Expt,all_t_axis(cur_set([1 end])) - cur_toffset);
        
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
        all_eye_ts = [all_eye_ts; eye_ts_interp' + cur_toffset];
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
    
    %     if strcmp(Expt_name,'G081')
    %         sac_dirs = [0 45 90 135];
    %         sim_sac_times = [0.7 1.4];
    %         sac_thresh = 0.5;
    %     else
    %         sac_dirs = [0 90];
    sim_sac_times = [0.7 1.4 2.1 2.8 3.5];
    sac_thresh = 1;
    %     end
    sac_start_expts = all_blockvec(interp_sac_start_inds);
    %     delta_sacpar = nan(size(saccades));
    delta_sac = sqrt(sac_deltaX.^2 + sac_deltaY.^2);
    %     for bb = 1:length(sac_dirs)
    %         cur_bar_expts = find(expt_bar_ori(cur_block_set) == sac_dirs(bb));
    %         cur_sac_set = find(ismember(sac_start_expts,cur_bar_expts));
    %
    %         cur_delta_sac_par = abs(sac_deltaX(cur_sac_set)*cos(sac_dirs(bb)*pi/180) + sac_deltaY(cur_sac_set)*sin(sac_dirs(bb)*pi/180));
    %         delta_sacpar(cur_sac_set) = cur_delta_sac_par;
    %     end
    %     is_gsac = delta_sacpar' > sac_thresh;
    is_gsac = delta_sac > sac_thresh;
    
    for ii = 1:length(cur_block_set)
        cur_eyevar(ii) = var(interp_eye_vals(all_blockvec == ii,2));
    end
    sim_sac_expts = find(cur_eyevar < 1);
    
    all_sim_sacs = [];
    sim_expt_inds = find(ismember(all_blockvec,sim_sac_expts));
    sim_sacs = cell(length(sim_sac_times),1);
    for ii = 1:length(sim_sac_times)
        sim_sacs{ii} = sim_expt_inds(all_tsince_start(sim_expt_inds(1:end-1)) < sim_sac_times(ii) & ...
            all_tsince_start(sim_expt_inds(2:end)) >= sim_sac_times(ii));
        all_sim_sacs = [all_sim_sacs; sim_sacs{ii}];
    end
    all_sim_msacs = all_sim_sacs(all_trial_Fs(all_trialvec(all_sim_sacs)) == 0.2);
    all_sim_gsacs = all_sim_sacs(all_trial_Fs(all_trialvec(all_sim_sacs)) == 3);
    
    imback_gs_expts = find(~ismember(1:length(cur_block_set),sim_sac_expts) & expt_imback(cur_block_set));
    grayback_gs_expts = find(~ismember(1:length(cur_block_set),sim_sac_expts) & ~expt_imback(cur_block_set));
    
    
    %define which saccades to use
    used_msacs = find(is_micro & ~is_blink);
    used_gsacs = find(is_gsac & ~is_blink);
    %     hori_expts = find(expt_bar_ori(cur_block_set) == 0);
    %     ver_expts = find(expt_bar_ori(cur_block_set) == 90);
    
    trial_start_inds = find(all_tsince_start < dt);
    simsac_trial_starts = trial_start_inds(ismember(all_blockvec(trial_start_inds),sim_sac_expts));
    %     grayback_trial_starts = trial_start_inds(ismember(all_blockvec(trial_start_inds),grayback_gs_expts));
    %     imback_trial_starts = trial_start_inds(ismember(all_blockvec(trial_start_inds),imback_gs_expts));
    
    used_msacs_grayback = used_msacs(ismember(all_blockvec(interp_sac_start_inds(used_msacs)),grayback_gs_expts));
    used_msacs_imback = used_msacs(ismember(all_blockvec(interp_sac_start_inds(used_msacs)),imback_gs_expts));
    used_gsacs_grayback = used_gsacs(ismember(all_blockvec(interp_sac_start_inds(used_gsacs)),grayback_gs_expts));
    used_gsacs_imback = used_gsacs(ismember(all_blockvec(interp_sac_start_inds(used_gsacs)),imback_gs_expts));
    %
    %     used_msacs_hori = used_msacs(ismember(all_blockvec(interp_sac_start_inds(used_msacs)),hori_expts));
    %     used_msacs_ver = used_msacs(ismember(all_blockvec(interp_sac_start_inds(used_msacs)),ver_expts));
    %     used_gsacs_hori = used_gsacs(ismember(all_blockvec(interp_sac_start_inds(used_gsacs)),hori_expts));
    %     used_gsacs_ver = used_gsacs(ismember(all_blockvec(interp_sac_start_inds(used_gsacs)),ver_expts));
    %
    %     simsacs_hori = all_sim_sacs(ismember(all_blockvec(all_sim_sacs),hori_expts));
    %     simsacs_ver = all_sim_sacs(ismember(all_blockvec(all_sim_sacs),ver_expts));
    
    %microsacs excluding bursts
    non_burst_msacs = used_msacs(~is_sacburst(used_msacs));
    burst_msacs = used_msacs(is_sacburst(used_msacs));
    
    %% FOR MODEL-BASED SAC-MOD ANALYSIS
    %     fprintf('Computing model-based MUA analysis\n');
    %     sac_bin_width = 1;
    %     sac_binspace = sac_bin_width*dt;
    %     sac_bin_edges = -(backlag*dt-sac_binspace/2):sac_binspace:(forwardlag*dt+sac_binspace/2);
    %     sac_bin_cents = 0.5*sac_bin_edges(1:end-1) + 0.5*sac_bin_edges(2:end);
    %     n_sac_bins = length(sac_bin_cents);
    %
    %     L2_params = create_L2_params([],[1 n_sac_bins],n_sac_bins);
    %     L2_params = create_L2_params(L2_params,n_sac_bins + [1 n_sac_bins],n_sac_bins);
    %     L2_params = create_L2_params(L2_params,2*n_sac_bins + [1 n_sac_bins],n_sac_bins);
    %     L2_params = create_L2_params(L2_params,[1 3*n_sac_bins],3*n_sac_bins,0);
    %
    %     gsac_inds = interp_sac_start_inds(used_gsacs);
    %     msac_inds = interp_sac_start_inds(used_msacs);
    %     trial_simsac_mat = zeros(length(all_t_axis),n_sac_bins);
    %     trial_gsac_mat = zeros(length(all_t_axis),n_sac_bins);
    %     trial_msac_mat = zeros(length(all_t_axis),n_sac_bins);
    %     for i = 1:n_sac_bins
    %         for cc = 1:sac_bin_width
    %             cur_inds = all_sim_sacs + round(sac_bin_cents(i)/dt) - floor(sac_bin_width/2) + cc;
    %             cur_inds(cur_inds < 1 | cur_inds > length(all_t_axis)) = [];
    %             trial_simsac_mat(cur_inds,i) = 1;
    %             cur_inds = gsac_inds + round(sac_bin_cents(i)/dt) - floor(sac_bin_width/2) + cc;
    %             cur_inds(cur_inds < 1 | cur_inds > length(all_t_axis)) = [];
    %             trial_gsac_mat(cur_inds,i) = 1;
    %             cur_inds = msac_inds + round(sac_bin_cents(i)/dt) - floor(sac_bin_width/2) + cc;
    %             cur_inds(cur_inds < 1 | cur_inds > length(all_t_axis)) = [];
    %             trial_msac_mat(cur_inds,i) = 1;
    %         end
    %     end
    %
    %     Xexpt = zeros(length(all_t_axis),length(cur_block_set)-1);
    %     for i = 1:length(cur_block_set)-1
    %         cur_set = find(all_blockvec==i);
    %         Xexpt(cur_set,i) = 1;
    %     end
    %
    %     cur_Xmat = [trial_gsac_mat trial_msac_mat trial_simsac_mat Xexpt];
    % %     clear trial_gsac_mat trial_msac_mat trial_simsac_mat Xexpt
    % %     cur_Xmat = [trial_gsac_mat trial_msac_mat Xexpt];
    % %     clear trial_gsac_mat trial_msac_mat Xexpt
    %
    %     exclude_mua_probes = [su_probes];
    %     include_mua_probes = setdiff(1:24,exclude_mua_probes);
    % %     ov_mua_rate = mean(all_spike_rate_norm(:,include_mua_probes),2);
    % %     [gen_data.mua_msac_avg,lags] = get_event_trig_avg(ov_mua_rate,interp_sac_start_inds(used_msacs),backlag,forwardlag,0,all_trialvec,0,used_inds);
    % %     [gen_data.mua_gsac_avg,lags] = get_event_trig_avg(ov_mua_rate,interp_sac_start_inds(used_gsacs),backlag,forwardlag,0,all_trialvec,0,used_inds);
    % %     [gen_data.mua_simsac_avg,lags] = get_event_trig_avg(ov_mua_rate,all_sim_sacs,backlag,forwardlag,0,all_trialvec,0,used_inds);
    % %     [gen_data.mua_nb_msac,lags] = get_event_trig_avg(ov_mua_rate,interp_sac_start_inds(non_burst_msacs),backlag,forwardlag,0,all_trialvec,0,used_inds);
    % %     [gen_data.mua_b_msac,lags] = get_event_trig_avg(ov_mua_rate,interp_sac_start_inds(burst_msacs),backlag,forwardlag,0,all_trialvec,0,used_inds);
    % for ii = 1:length(include_mua_probes)
    %     fprintf('Fitting kernels for probe %d of %d\n',ii,length(include_mua_probes));
    %     Robs = sum(all_binned_spikes(used_inds,include_mua_probes(ii)),2);
    %     [fitmod] = regGLM_fit(cur_Xmat(used_inds,:),Robs,L2_params,[1000 1000 1000 50],[],[],0);
    %     kern_fits(ii).gsac_kern = fitmod.K((1:n_sac_bins));
    %     kern_fits(ii).msac_kern = fitmod.K((1:n_sac_bins) + n_sac_bins);
    %     kern_fits(ii).simsac_kern = fitmod.K((1:n_sac_bins) + 2*n_sac_bins);
    % %     gen_data.mua_gsac_kern = fitmod.K((1:n_sac_bins));
    % %     gen_data.mua_msac_kern = fitmod.K((1:n_sac_bins) + n_sac_bins);
    % %     gen_data.mua_simsac_kern = fitmod.K((1:n_sac_bins) + 2*n_sac_bins);
    % end
    %     clear cur_Xmat
    %
    %
    %
    %     L2_params = create_L2_params([],[1 n_sac_bins],n_sac_bins);
    %     L2_params = create_L2_params(L2_params,[1 n_sac_bins],n_sac_bins,0);
    %     cur_Xmat = [trial_msac_mat Xexpt];
    %     for ii = 1:length(include_mua_probes)
    %         fprintf('Fitting kernels for probe %d of %d\n',ii,length(include_mua_probes));
    %         Robs = sum(all_binned_spikes(used_inds,include_mua_probes(ii)),2);
    %         [fitmod] = regGLM_fit(cur_Xmat(used_inds,:),Robs,L2_params,[500 5],[],[],1);
    %
    %         [LL, penLL, pred_rate, G] = regGLM_eval(fitmod,Robs,cur_Xmat(used_inds,:));
    %         [simsac_predavg(ii,:),lags] = get_event_trig_avg(pred_rate,find(ismember(used_inds,all_sim_msacs)),backlag,forwardlag,nboot,used_trialvec(used_inds),0);
    %         [simsac_obsavg(ii,:),lags] = get_event_trig_avg(all_spike_rate(:,ii),all_sim_msacs,backlag,forwardlag,nboot,used_trialvec,0);
    %
    %         %     gen_data.mua_gsac_kern = fitmod.K((1:n_sac_bins));
    %         %     gen_data.mua_msac_kern = fitmod.K((1:n_sac_bins) + n_sac_bins);
    %         %     gen_data.mua_simsac_kern = fitmod.K((1:n_sac_bins) + 2*n_sac_bins);
    %     end
    
    %%     SACCADE TIMING ANALYSIS
    binned_msacs = hist(interp_sac_start_inds(used_msacs),1:length(all_t_axis));
    binned_gsacs = hist(interp_sac_start_inds(used_gsacs),1:length(all_t_axis));
    binned_msac_sm = jmm_smooth_1d_cor(binned_msacs,sua_sm_sig);
    binned_gsac_sm = jmm_smooth_1d_cor(binned_gsacs,sua_sm_sig);
    
    %trial averages (USE ALL OF TRIAL HERE)
    %     [gen_data.msac_simsac_trial_avg,trial_lags] = get_event_trig_avg(binned_msac_sm,simsac_trial_starts,0,round(trial_dur/dt),0,all_trialvec,0);
    %     [gen_data.msac_grayback_trial_avg,trial_lags] = get_event_trig_avg(binned_msac_sm,grayback_trial_starts,0,round(trial_dur/dt),0,all_trialvec,0);
    %     [gen_data.msac_imback_trial_avg,trial_lags] = get_event_trig_avg(binned_msac_sm,imback_trial_starts,0,round(trial_dur/dt),0,all_trialvec,0);
    %     [gen_data.gsac_simsac_trial_avg,trial_lags] = get_event_trig_avg(binned_gsac_sm,simsac_trial_starts,0,round(trial_dur/dt),0,all_trialvec,0);
    %     [gen_data.gsac_grayback_trial_avg,trial_lags] = get_event_trig_avg(binned_gsac_sm,grayback_trial_starts,0,round(trial_dur/dt),0,all_trialvec,0);
    %     [gen_data.gsac_imback_trial_avg,trial_lags] = get_event_trig_avg(binned_gsac_sm,imback_trial_starts,0,round(trial_dur/dt),0,all_trialvec,0);
    
    maxlag = round(0.5/dt);
    [gen_data.msac_acorr,acorr_lags] = xcov(binned_msac_sm,maxlag,'coeff');
    [gen_data.gsac_acorr,acorr_lags] = xcov(binned_gsac_sm,maxlag,'coeff');
    [gen_data.msac_gsac_xcorr,acorr_lags] = xcov(binned_gsac_sm,binned_msac_sm,maxlag,'coeff');
    
    %% COMPUTE TRIG AVGS FOR MUA
    nboot = [];
    %set trial numbers to Inf so they don't get included in trig averaging
    used_trialvec = all_trialvec; used_trialvec(used_inds == 0) = Inf;
    clear mua_data
    for cc = 1:n_probes
        
        fprintf('Computing trig avgs for MUA %d of %d\n',cc,24);
        
        cur_used_inds = used_inds;
        
        %         %trial averages (USE ALL OF TRIAL HERE)
        %         [mua_data(cc).simsac_trial_avg,trial_lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),simsac_trial_starts,0,round(trial_dur/dt),nboot,all_trialvec,0,cur_poss_inds);
        %         [mua_data(cc).grayback_trial_avg,trial_lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),grayback_trial_starts,0,round(trial_dur/dt),nboot,all_trialvec,0,cur_poss_inds);
        %         [mua_data(cc).imback_trial_avg,trial_lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),imback_trial_starts,0,round(trial_dur/dt),nboot,all_trialvec,0,cur_poss_inds);
        
        %general averages
        [mua_data(cc).msac_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),interp_sac_start_inds(used_msacs),backlag,forwardlag,nboot,used_trialvec,0,cur_used_inds);
        [mua_data(cc).gsac_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),interp_sac_start_inds(used_gsacs),backlag,forwardlag,nboot,used_trialvec,0,cur_used_inds);
        [mua_data(cc).simsac_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),all_sim_sacs,backlag,forwardlag,nboot,used_trialvec,0,cur_used_inds);
        [mua_data(cc).sim_msac_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),all_sim_msacs,backlag,forwardlag,nboot,used_trialvec,0,cur_used_inds);
        [mua_data(cc).sim_gsac_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),all_sim_gsacs,backlag,forwardlag,nboot,used_trialvec,0,cur_used_inds);
        
        %         [mua_data(cc).msac_end_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),interp_sac_stop_inds(used_msacs),backlag,forwardlag,nboot,used_trialvec,0,cur_used_inds);
        %         [mua_data(cc).gsac_end_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),interp_sac_stop_inds(used_gsacs),backlag,forwardlag,nboot,used_trialvec,0,cur_used_inds);
        %
        %background dependent
        [mua_data(cc).msac_gray_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),interp_sac_start_inds(used_msacs_grayback),backlag,forwardlag,nboot,used_trialvec,0,cur_used_inds);
        [mua_data(cc).gsac_gray_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),interp_sac_start_inds(used_gsacs_grayback),backlag,forwardlag,nboot,used_trialvec,0,cur_used_inds);
        [mua_data(cc).msac_im_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),interp_sac_start_inds(used_msacs_imback),backlag,forwardlag,nboot,used_trialvec,0,cur_used_inds);
        [mua_data(cc).gsac_im_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),interp_sac_start_inds(used_gsacs_imback),backlag,forwardlag,nboot,used_trialvec,0,cur_used_inds);
        %
        %         %saccade direction dependent
        %         [mua_data(cc).msac_ver_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),interp_sac_start_inds(used_msacs_ver),backlag,forwardlag,nboot,used_trialvec,0,cur_used_inds);
        %         [mua_data(cc).gsac_ver_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),interp_sac_start_inds(used_gsacs_ver),backlag,forwardlag,nboot,used_trialvec,0,cur_used_inds);
        %         [mua_data(cc).simsac_ver_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),simsacs_ver,backlag,forwardlag,nboot,used_trialvec,0,cur_used_inds);
        %         [mua_data(cc).msac_hor_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),interp_sac_start_inds(used_msacs_hori),backlag,forwardlag,nboot,used_trialvec,0,cur_used_inds);
        %         [mua_data(cc).gsac_hor_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),interp_sac_start_inds(used_gsacs_hori),backlag,forwardlag,nboot,used_trialvec,0,cur_used_inds);
        %         [mua_data(cc).simsac_hor_avg,lags] = get_event_trig_avg(all_mua_rate_norm(:,cc),simsacs_hori,backlag,forwardlag,nboot,used_trialvec,0,cur_used_inds);
        
        mua_data(cc).avg_rate = mean(all_binned_mua(cur_used_inds,cc));
        mua_data(cc).tot_nspikes = sum(all_binned_mua(cur_used_inds,cc));
        
        %start and stop trig avgs
        %micro sac burst and non-burst
    end
    
    %%
    close all
    sim_msac_avgs = [mua_data(:).sim_msac_avg];
    sim_gsac_avgs = [mua_data(:).sim_gsac_avg];
    msac_avgs = [mua_data(:).msac_avg];
    gsac_avgs = [mua_data(:).gsac_avg];
    msac_gray_avgs = [mua_data(:).msac_gray_avg];
    gsac_gray_avgs = [mua_data(:).gsac_gray_avg];
    msac_im_avgs = [mua_data(:).msac_im_avg];
    gsac_im_avgs = [mua_data(:).gsac_im_avg];
    
    figure; hold on
    h1=shadedErrorBar(lags*dt,mean(gsac_avgs,2),std(gsac_avgs,[],2)/sqrt(24),{'color','r'});
    h2=shadedErrorBar(lags*dt,mean(msac_avgs,2),std(msac_avgs,[],2)/sqrt(24),{'color','b'});
    h3=shadedErrorBar(lags*dt,mean(sim_gsac_avgs,2),std(sim_gsac_avgs,[],2)/sqrt(24),{'color','k'});
    h4=shadedErrorBar(lags*dt,mean(sim_msac_avgs,2),std(sim_msac_avgs,[],2)/sqrt(24),{'color','m'});
    legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'Guided','Micro','Sim-guided','Sim-micro'});
    xlim([-0.2 0.45]);
    yl = ylim();
    line([0 0],yl,'color','k','linestyle','--');
    line([-0.2 0.45],[1 1],'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',12);
    ylabel('Relative firing rate','fontsize',12);
    
    figure; hold on
    h1=shadedErrorBar(lags*dt,mean(gsac_gray_avgs,2),std(gsac_gray_avgs,[],2)/sqrt(24),{'color','r'});
    h2=shadedErrorBar(lags*dt,mean(msac_gray_avgs,2),std(msac_gray_avgs,[],2)/sqrt(24),{'color','b'});
    h3=shadedErrorBar(lags*dt,mean(gsac_im_avgs,2),std(gsac_im_avgs,[],2)/sqrt(24),{'color','k'});
    h4=shadedErrorBar(lags*dt,mean(msac_im_avgs,2),std(msac_im_avgs,[],2)/sqrt(24),{'color','m'});
    legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'Guided-gray','Micro-gray','Guided-im','Micro-im'});
    xlim([-0.2 0.45]);
    yl = ylim();
    line([0 0],yl,'color','k','linestyle','--');
    line([-0.2 0.45],[1 1],'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',12);
    ylabel('Relative firing rate','fontsize',12);
    
    figure; hold on
    subplot(2,2,1);
    imagesc(lags*dt,1:24,gsac_avgs'); caxis([0.75 1.25]);
    xlim([-0.2 0.45]);
    yl = ylim();
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',12);
    ylabel('Probe number','fontsize',12);
    title('Guided')
    subplot(2,2,2);
    imagesc(lags*dt,1:24,msac_avgs'); caxis([0.75 1.25]);
    xlim([-0.2 0.45]);
    yl = ylim();
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',12);
    ylabel('Probe number','fontsize',12);
    title('Micro')
    subplot(2,2,3);
    imagesc(lags*dt,1:24,sim_gsac_avgs'); caxis([0.75 1.25]);
    xlim([-0.2 0.45]);
    yl = ylim();
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',12);
    ylabel('Probe number','fontsize',12);
    title('Sim-Guided')
    subplot(2,2,4);
    imagesc(lags*dt,1:24,sim_msac_avgs'); caxis([0.75 1.25]);
    xlim([-0.2 0.45]);
    yl = ylim();
    line([0 0],yl,'color','k','linestyle','--');
    xlabel('Time since saccade onset (s)','fontsize',12);
    ylabel('Probe number','fontsize',12);
    title('Sim-Micro')
    
    %% COMPUTE TRIG AVGS FOR SUA
    nboot = 50;
    %set trial numbers to Inf so they don't get included in trig averaging
    used_trialvec = all_trialvec; used_trialvec(used_inds == 0) = Inf;
    clear sua_data
    for ss = 1:length(SU_numbers)
        
        fprintf('Computing trig avgs for SU %d of %d\n',ss,length(SU_numbers));
        cur_used_blocks = find(~isnan(SU_block_probes(ss,:)));
        cur_use_inds = find(ismember(all_blockvec,cur_used_blocks));
        
        if length(cur_used_blocks) >= 3
            sua_data(ss).used = true;
            %         %trial averages (USE ALL OF TRIAL HERE)
            %         [sua_data(ss).simsac_trial_avg,trial_lags,sua_data(ss).simsac_trial_sem,sua_data(ss).nused_simsac_trials] = get_event_trig_avg(all_sua_rate_norm(:,ss),simsac_trial_starts,0,round(trial_dur/dt),nboot,all_trialvec,0,cur_use_inds);
            %         [sua_data(ss).grayback_trial_avg,trial_lags,sua_data(ss).grayback_trial_sem,sua_data(ss).nused_grayback_trials] = get_event_trig_avg(all_sua_rate_norm(:,ss),grayback_trial_starts,0,round(trial_dur/dt),nboot,all_trialvec,0,cur_use_inds);
            %         [sua_data(ss).imback_trial_avg,trial_lags,sua_data(ss).imback_trial_sem,sua_data(ss).nused_imback_trials] = get_event_trig_avg(all_sua_rate_norm(:,ss),imback_trial_starts,0,round(trial_dur/dt),nboot,all_trialvec,0,cur_use_inds);
            
            %general averages
            [sua_data(ss).msac_avg,lags,sua_data(ss).msac_sem,sua_data(ss).nused_msac] = get_event_trig_avg(all_sua_rate_norm(:,ss),interp_sac_start_inds(used_msacs),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
            [sua_data(ss).gsac_avg,lags,sua_data(ss).gsac_sem,sua_data(ss).nused_gsac] = get_event_trig_avg(all_sua_rate_norm(:,ss),interp_sac_start_inds(used_gsacs),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
            [sua_data(ss).simsac_avg,lags,sua_data(ss).simsac_sem,sua_data(ss).nused_simsac] = get_event_trig_avg(all_sua_rate_norm(:,ss),all_sim_sacs,backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
            [sua_data(ss).sim_msac_avg,lags,sua_data(ss).sim_msac_sem,sua_data(ss).nused_sim_msac] = get_event_trig_avg(all_sua_rate_norm(:,ss),all_sim_msacs,backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
            [sua_data(ss).sim_gsac_avg,lags,sua_data(ss).sim_gsac_sem,sua_data(ss).nused_sim_gsac] = get_event_trig_avg(all_sua_rate_norm(:,ss),all_sim_gsacs,backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
            
            %         [sua_data(ss).msac_end_avg,lags,sua_data(ss).msac_end_sem,sua_data(ss).nused_msac] = get_event_trig_avg(all_sua_rate_norm(:,ss),interp_sac_stop_inds(used_msacs),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
            %         [sua_data(ss).gsac_end_avg,lags,sua_data(ss).gsac_end_sem,sua_data(ss).nused_gsac] = get_event_trig_avg(all_sua_rate_norm(:,ss),interp_sac_stop_inds(used_gsacs),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
            
            %background dependent
            [sua_data(ss).msac_gray_avg,lags,sua_data(ss).msac_gray_sem,sua_data(ss).nused_msac_gray] = get_event_trig_avg(all_sua_rate_norm(:,ss),interp_sac_start_inds(used_msacs_grayback),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
            [sua_data(ss).gsac_gray_avg,lags,sua_data(ss).gsac_gray_sem,sua_data(ss).nused_gsac_gray] = get_event_trig_avg(all_sua_rate_norm(:,ss),interp_sac_start_inds(used_gsacs_grayback),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
            [sua_data(ss).msac_im_avg,lags,sua_data(ss).msac_im_sem,sua_data(ss).nused_msac_im] = get_event_trig_avg(all_sua_rate_norm(:,ss),interp_sac_start_inds(used_msacs_imback),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
            [sua_data(ss).gsac_im_avg,lags,sua_data(ss).gsac_im_sem,sua_data(ss).nused_gsac_im] = get_event_trig_avg(all_sua_rate_norm(:,ss),interp_sac_start_inds(used_gsacs_imback),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
            
            %         %sassade direction dependent
            %         [sua_data(ss).msac_ver_avg,lags,sua_data(ss).msac_ver_sem,sua_data(ss).nused_msac_ver] = get_event_trig_avg(all_sua_rate_norm(:,ss),interp_sac_start_inds(used_msacs_ver),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
            %         [sua_data(ss).gsac_ver_avg,lags,sua_data(ss).gsac_ver_sem,sua_data(ss).nused_gsac_ver] = get_event_trig_avg(all_sua_rate_norm(:,ss),interp_sac_start_inds(used_gsacs_ver),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
            %         [sua_data(ss).simsac_ver_avg,lags,sua_data(ss).simsac_ver_sem,sua_data(ss).nused_simsac_ver] = get_event_trig_avg(all_sua_rate_norm(:,ss),simsacs_ver,backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
            %         [sua_data(ss).msac_hor_avg,lags,sua_data(ss).msac_hor_sem,sua_data(ss).nused_msac_hor] = get_event_trig_avg(all_sua_rate_norm(:,ss),interp_sac_start_inds(used_msacs_hori),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
            %         [sua_data(ss).gsac_hor_avg,lags,sua_data(ss).gsac_hor_sem,sua_data(ss).used_gsac_hor] = get_event_trig_avg(all_sua_rate_norm(:,ss),interp_sac_start_inds(used_gsacs_hori),backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
            %         [sua_data(ss).simsac_hor_avg,lags,sua_data(ss).simsac_hor_sem,sua_data(ss).nused_simsac_hor] = get_event_trig_avg(all_sua_rate_norm(:,ss),simsacs_hori,backlag,forwardlag,nboot,used_trialvec,0,cur_use_inds);
            
            sua_data(ss).avg_rate = mean(all_binned_sua(cur_use_inds,ss));
            sua_data(ss).tot_nspikes = sum(all_binned_sua(cur_use_inds,ss));
            %start and stop trig avgs
            %micro sac burst and non-burst
        else
            sua_data(ss).used = false;
        end
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
    %     anal_params.trial_lags = trial_lags;
    %     anal_params.acorr_lags = acorr_lags;
    %     anal_params.sac_bin_cents = sac_bin_cents;
    cd(save_dir)
    sname = 'full_sacmod_data2';
    save(sname,'sua_data','SU_block_probes','mua_data','anal_params','gen_data');
    clear sua_data su_probes mua_data gen_data
end