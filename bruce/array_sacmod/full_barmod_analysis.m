clear all
% close all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';

% Expt_nums = [85 86 87 88 89 91 92 93];
Expt_nums = 95;

%% SET ANALYSIS PARAMETERS
stim_fs = 100; %in Hz
dt = 0.01;

min_trial_dur = 0.5;
beg_buffer = 0.3;

backlag = round(0.25/dt);
forwardlag = round(0.55/dt);

sua_sm_sig = (0.01/dt);
mua_sm_sig = (0.005/dt);

flen = 12;
nPix = 24;
stim_params = NIMcreate_stim_params([flen nPix],dt);
Fr = 1;

%% LOAD OVERALL SU DATA
load ~/Analysis/bruce/summary_analysis/su_data.mat

%% LOOP OVER EXPTS
for ex = 1:length(Expt_nums)
    % for ex = 8
    clear expt* included_type
    
    Expt_name = sprintf('G0%d',Expt_nums(ex));
    data_dir = [dir_prefix '/Data/bruce/' Expt_name];
    addpath([dir_prefix '/James_scripts/bruce/G081/'])
    
    cd(data_dir);
    
    load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
    if ~strcmp(Expt_name,'G081')
        load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data
    end
    
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
    
    if strcmp(Expt_name,'G081')
        trial_dur = 2;
    else
        trial_dur = 4;
    end
    
    %%
    if strcmp(Expt_name,'G081')
        for i = 1:length(Expts)
            if strcmp(Expts{i}.Header.expname,'grating.OpXseRC') | strcmp(Expts{i}.Header.expname,'grating.OpRC')
                is_bar_expt(i) = 1;
            else
                is_bar_expt(i) = 0;
            end
            
            if strcmp(Expts{i}.Stimvals.Bs,'image')
                expt_imback(i) = 1;
            else
                expt_imback(i) = 0;
            end
            
            expt_sim_sacs(i) = Expts{i}.Stimvals.ijump;
            expt_bar_ori(i) = Expts{i}.Stimvals.or;
        end
        expt_has_ds = (expt_sim_sacs == 0);
        expt_bar_ori(expt_bar_ori == -45) = 135;
        expt_binoc = zeros(size(expt_bar_ori));
        expt_imback = expt_imback';
        cur_block_set = find(is_bar_expt);
    else
        if strcmp(Expt_name,'G093')
            include_expts = {'rls.FaXwi','rls.FaXwiXimi'};
        elseif strcmp(Expt_name,'G095')
            include_expts = {'rls.Fa','rls.FaXimi','rls.FaXFaXFs'};
        else
            include_expts = {'rls.Fa', 'rls.FaXimi'};
        end
        for ii = 1:length(Expts)
            expt_names{ii} = Expts{ii}.Header.expname;
            expt_dds(ii) = Expts{ii}.Stimvals.dd;
            expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
            expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
            expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
            expt_imback(ii) = isfield(Expts{ii}.Trials,'imi');
            included_type(ii) = any(strcmp(expt_names{ii},include_expts));
        end
        expt_has_ds(isnan(expt_has_ds)) = 0;
        expt_has_ds(expt_has_ds == -1) = 0;
        expt_binoc(isnan(expt_binoc)) = 0;
        cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1);
        
        if strcmp(Expt_name,'G087')
            cur_block_set(cur_block_set == 15) = []; %only 6 trials and causes problems
        end
    end
    
    if strcmp(Expt_name,'G087')
        cur_block_set(cur_block_set == 15) = []; %only 6 trials and causes problems
    end
    if strcmp(Expt_name,'G093')
        cur_block_set(cur_block_set ==  28) = []; %only 6 trials and causes problems
    end
    
    sim_sac_expts = find(~expt_has_ds(cur_block_set));
    imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
    grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');
    
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
%     all_Xmat = [];
all_stim_mat = [];
    all_t_axis = [];
    all_t_bin_edges = [];
    all_tsince_start = [];
    all_blockvec = [];
    all_trialvec = [];
    all_trial_wi = [];
    all_trial_start_times = [];
    all_trial_end_times = [];
    all_bin_edge_pts = [];
    all_spk_times = cell(96,1);
    all_clust_ids = cell(length(su_probes),1);
    for ee = 1:length(cur_block_set);
        if ismember(ee,grayback_gs_expts)
            fprintf('Expt %d Block %d of %d; grayback GS, ori:%d\n',Expt_nums(ex),ee,length(cur_block_set),expt_bar_ori(cur_block_set(ee)));
        elseif ismember(ee,imback_gs_expts)
            fprintf('Expt %d Block %d of %d; imback GS, ori:%d\n',Expt_nums(ex),ee,length(cur_block_set),expt_bar_ori(cur_block_set(ee)));
        elseif ismember(ee,sim_sac_expts)
            fprintf('Expt %d Block %d of %d; SimSac, ori:%d\n',Expt_nums(ex),ee,length(cur_block_set),expt_bar_ori(cur_block_set(ee)));
        else
            fprintf('Expt %d Block %d of %d;  UNMATCHED EXPT TYPE\n',Expt_nums(ex),ee,length(cur_block_set));
        end
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
        rpt_trials = false;
        if length(un_ids) < length(trial_ids)
            rpt_trials = true;
            fprintf('Warning, repeat trial inds detected!\n');
            use_trials = [];
        else
            use_trials = find(trial_durs >= min_trial_dur);
        end
        %         trial_durs = trial_durs(id_inds);
        %         trial_start_times = trial_start_times(id_inds);
        %         trial_end_times = trial_end_times(id_inds);
        
        
        all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
        all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)');
        
        if strcmp(Expt_name,'G093')
            trial_wi = [Expts{cur_block}.Trials(:).wi];
            trial_wi = trial_wi(id_inds);
            all_trial_wi = cat(1,all_trial_wi,trial_wi(use_trials)');
        end
        
        fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
        load(fname);
        buffer_pix = floor((expt_npix(cur_block) - nPix)/2);
        cur_use_pix = (1:nPix) + buffer_pix;
        
        
        n_trials = length(use_trials);
        for tt = 1:n_trials
            cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start'/1e4;
            n_frames = size(left_stim_mats{use_trials(tt)},1);
            if length(cur_stim_times) == 1
                cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
                cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
            end
            %             cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
            cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
            
            cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
            
            if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > min_trial_dur/dt
                use_frames = min(length(cur_stim_times),n_frames);
                cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
%                 bar_Xmat = create_time_embedding(cur_stim_mat,stim_params);
                
                all_stim_times = [all_stim_times; cur_stim_times];
                all_t_axis = [all_t_axis; cur_t_axis];
%                 all_Xmat = [all_Xmat; bar_Xmat];
all_stim_mat = [all_stim_mat; cur_stim_mat];
                all_t_bin_edges = [all_t_bin_edges; cur_t_edges];
                all_tsince_start = [all_tsince_start; cur_tsince_start];
                all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
                all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
                all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
            end
        end
        trial_cnt = trial_cnt + n_trials;
    end
    all_Xmat = create_time_embedding(all_stim_mat,stim_params);
    
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
    if strcmp(Expt_name,'G093')
        un_wi_vals = unique(all_trial_wi);
        use_wi_trials = find(all_trial_wi == un_wi_vals(2));
        used_inds = used_inds(ismember(all_trialvec(used_inds),use_wi_trials));
    end
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
        end
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
        if ~isempty(cur_set)
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
    
    if strcmp(Expt_name,'G081')
        sac_dirs = [0 45 90 135];
        sim_sac_times = [0.7 1.4];
        sac_thresh = 0.5;
    else
        sac_dirs = [0 90];
        sim_sac_times = [0.7 1.4 2.1 2.8 3.5];
        sac_thresh = 1;
    end
    sac_start_expts = all_blockvec(interp_sac_start_inds);
    delta_sacpar = nan(size(saccades));
    for bb = 1:length(sac_dirs)
        cur_bar_expts = find(expt_bar_ori(cur_block_set) == sac_dirs(bb));
        cur_sac_set = find(ismember(sac_start_expts,cur_bar_expts));
        
        cur_delta_sac_par = abs(sac_deltaX(cur_sac_set)*cos(sac_dirs(bb)*pi/180) + sac_deltaY(cur_sac_set)*sin(sac_dirs(bb)*pi/180));
        delta_sacpar(cur_sac_set) = cur_delta_sac_par;
    end
    is_gsac = delta_sacpar' > sac_thresh;
    
    all_sim_sacs = [];
    sim_expt_inds = find(ismember(all_blockvec,sim_sac_expts));
    sim_sacs = cell(length(sim_sac_times),1);
    for ii = 1:length(sim_sac_times)
        sim_sacs{ii} = sim_expt_inds(all_tsince_start(sim_expt_inds(1:end-1)) < sim_sac_times(ii) & ...
            all_tsince_start(sim_expt_inds(2:end)) >= sim_sac_times(ii));
        all_sim_sacs = [all_sim_sacs; sim_sacs{ii}];
    end
    
    %define which saccades to use
    used_msacs = find(is_micro & ~is_blink);
    used_gsacs = find(is_gsac & ~is_blink);
    msac_start_inds = interp_sac_start_inds(used_msacs);
    gsac_start_inds = interp_sac_start_inds(used_gsacs);
    
    %% CREATE SACCADE PREDICTOR MATS
    sac_backlag = 0.2;
    sac_forlag = 0.5;
    sac_bin_width = 2; %number of dt bins
    sac_binspace = sac_bin_width*dt;
    sac_bin_edges = -(sac_backlag-sac_binspace/2):sac_binspace:(sac_forlag+sac_binspace/2);
    sac_bin_cents = 0.5*sac_bin_edges(1:end-1) + 0.5*sac_bin_edges(2:end);
    n_sac_bins = length(sac_bin_cents);
    
    trial_gsac_mat = zeros(length(all_t_axis),n_sac_bins);
    trial_msac_mat = zeros(length(all_t_axis),n_sac_bins);
    for i = 1:n_sac_bins
        for bb = 1:sac_bin_width
            cur_inds = gsac_start_inds + round(sac_bin_cents(i)/dt) - floor(sac_bin_width/2) + bb;
            cur_inds(cur_inds < 1 | cur_inds > length(all_t_axis)) = [];
            trial_gsac_mat(cur_inds,i) = 1;
            
            cur_inds = msac_start_inds + round(sac_bin_cents(i)/dt) - floor(sac_bin_width/2) + bb;
            cur_inds(cur_inds < 1 | cur_inds > length(all_t_axis)) = [];
            trial_msac_mat(cur_inds,i) = 1;
        end
    end
    
    %% MUA STA/STC ANALYSIS
    hor_blocks = find(expt_bar_ori(cur_block_set) == 0);
    ver_blocks = find(expt_bar_ori(cur_block_set) == 90);
    nneg = 5; npos = 5;
    stc_thresh = -1.5e-3;
    
    for cc = 1:96
        
        fprintf('Computing trig avgs for MUA %d of %d\n',cc,96);
        
        %for horizontal bar blocks
        su_probe_ind = find(su_probes == cc);
        if ~isempty(su_probe_ind)
            cur_used_blocks = hor_blocks(su_used_blocks(hor_blocks,su_probe_ind) == 0); %blocks when NO SU
            cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
        else
            cur_used_blocks = hor_blocks; %blocks when NO SU
            cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
        end
        cur_used_inds = used_inds(ismember(used_inds,cur_poss_inds));
        
        if ~isempty(cur_used_inds)
            spikebins = convert_to_spikebins(all_binned_spikes(cur_used_inds,cc));
            spike_cond_stim = all_Xmat(cur_used_inds(spikebins),:);
            sta      = mean(spike_cond_stim) - mean(all_Xmat(cur_used_inds,:));
            sta = sta/norm(sta);
            proj_mat = sta'/(sta*sta')*sta;
            stim_proj = all_Xmat(cur_used_inds,:) - all_Xmat(cur_used_inds,:)*proj_mat;
            % stim_proj = stim_emb;
            stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
            [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
            stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
            
            mua_data(cc).stc_hor_use = true;
            mua_data(cc).sta_hor = sta;
            mua_data(cc).stcs_hor = stcs;
            mua_data(cc).evecs_hor = diag(evals);
            cur_evec_diff = diff(flipud(mua_data(cc).evecs_hor));
            mua_data(cc).npos_stc_hor = find(cur_evec_diff > stc_thresh,1,'first')-1;
            mua_data(cc).nneg_stc_hor = find(cur_evec_diff(end:-1:1) > stc_thresh,1,'first')-1;
        else
            mua_data(cc).stc_hor_use = false;
            mua_data(cc).sta_hor = nan;
            mua_data(cc).stcs_hor = nan;
            mua_data(cc).evecs_hor = nan;
            mua_data(cc).npos_stc_hor = nan;
            mua_data(cc).nneg_stc_hor = nan;
        end
        
        %for vertical bar blocks
        su_probe_ind = find(su_probes == cc);
        if ~isempty(su_probe_ind)
            cur_used_blocks = ver_blocks(su_used_blocks(ver_blocks,su_probe_ind) == 0); %blocks when NO SU
            cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
        else
            cur_used_blocks = ver_blocks; %blocks when NO SU
            cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
        end
        cur_used_inds = used_inds(ismember(used_inds,cur_poss_inds));
        
        if ~isempty(cur_used_inds)
            spikebins = convert_to_spikebins(all_binned_spikes(cur_used_inds,cc));
            spike_cond_stim = all_Xmat(cur_used_inds(spikebins),:);
            sta      = mean(spike_cond_stim) - mean(all_Xmat(cur_used_inds,:));
            sta = sta/norm(sta);
            proj_mat = sta'/(sta*sta')*sta;
            stim_proj = all_Xmat(cur_used_inds,:) - all_Xmat(cur_used_inds,:)*proj_mat;
            % stim_proj = stim_emb;
            stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
            [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
            stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
            
            mua_data(cc).stc_ver_use = true;
            mua_data(cc).sta_ver = sta;
            mua_data(cc).stcs_ver = stcs;
            mua_data(cc).evecs_ver = diag(evals);
            cur_evec_diff = diff(flipud(mua_data(cc).evecs_ver));
            mua_data(cc).npos_stc_ver = find(cur_evec_diff > stc_thresh,1,'first')-1;
            mua_data(cc).nneg_stc_ver = find(cur_evec_diff(end:-1:1) > stc_thresh,1,'first')-1;
        else
            mua_data(cc).stc_ver_use = false;
            mua_data(cc).sta_ver = nan;
            mua_data(cc).stcs_ver = nan;
            mua_data(cc).evecs_ver = nan;
            mua_data(cc).npos_stc_ver = nan;
            mua_data(cc).nneg_stc_ver = nan;
        end
    end
    
    %% SUA STA/STC ANALYSIS
    
    for ss = 1:length(su_probes)
        
        fprintf('Computing trig avgs for SUA %d of %d\n',ss,length(su_probes));
        
        %for horizontal bars
        cur_used_blocks = hor_blocks(su_used_blocks(hor_blocks,ss)); %blocks when SU
        cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
        cur_used_inds = used_inds(ismember(used_inds,cur_poss_inds));
        
        if ~isempty(cur_used_inds)
            spikebins = convert_to_spikebins(all_binned_spikes(cur_used_inds,su_probes(ss)));
            spike_cond_stim = all_Xmat(cur_used_inds(spikebins),:);
            sta      = mean(spike_cond_stim) - mean(all_Xmat(cur_used_inds,:));
            sta = sta/norm(sta);
            proj_mat = sta'/(sta*sta')*sta;
            stim_proj = all_Xmat(cur_used_inds,:) - all_Xmat(cur_used_inds,:)*proj_mat;
            % stim_proj = stim_emb;
            stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
            [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
            stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
            
            sua_data(ss).stc_hor_use = true;
            sua_data(ss).sta_hor = sta;
            sua_data(ss).stcs_hor = stcs;
            sua_data(ss).evecs_hor = diag(evals);
            cur_evec_diff = diff(flipud(sua_data(ss).evecs_hor));
            sua_data(ss).npos_stc_hor = find(cur_evec_diff > stc_thresh,1,'first')-1;
            sua_data(ss).nneg_stc_hor = find(cur_evec_diff(end:-1:1) > stc_thresh,1,'first')-1;
        else
            sua_data(ss).stc_hor_use = false;
            sua_data(ss).sta_hor = nan;
            sua_data(ss).stcs_hor = nan;
            sua_data(ss).evecs_hor = nan;
            sua_data(ss).npos_stc_hor = nan;
            sua_data(ss).nneg_stc_hor = nan;
        end
        
%         %for vertical bars
%         cur_used_blocks = ver_blocks(su_used_blocks(ver_blocks,ss)); %blocks when SU
%         cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
%         cur_used_inds = used_inds(ismember(used_inds,cur_poss_inds));
%         
%         if ~isempty(cur_used_inds)
%             spikebins = convert_to_spikebins(all_binned_spikes(cur_used_inds,su_probes(ss)));
%             spike_cond_stim = all_Xmat(cur_used_inds(spikebins),:);
%             sta      = mean(spike_cond_stim) - mean(all_Xmat(cur_used_inds,:));
%             sta = sta/norm(sta);
%             proj_mat = sta'/(sta*sta')*sta;
%             stim_proj = all_Xmat(cur_used_inds,:) - all_Xmat(cur_used_inds,:)*proj_mat;
%             % stim_proj = stim_emb;
%             stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
%             [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
%             stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
%             
%             sua_data(ss).stc_ver_use = true;
%             sua_data(ss).sta_ver = sta;
%             sua_data(ss).stcs_ver = stcs;
%             sua_data(ss).evecs_ver = diag(evals);
%             cur_evec_diff = diff(flipud(sua_data(ss).evecs_ver));
%             sua_data(ss).npos_stc_ver = find(cur_evec_diff > stc_thresh,1,'first')-1;
%             sua_data(ss).nneg_stc_ver = find(cur_evec_diff(end:-1:1) > stc_thresh,1,'first')-1;
%         else
%             sua_data(ss).stc_ver_use = false;
%             sua_data(ss).sta_ver = nan;
%             sua_data(ss).stcs_ver = nan;
%             sua_data(ss).evecs_ver = nan;
%             sua_data(ss).npos_stc_ver = nan;
%             sua_data(ss).nneg_stc_ver = nan;
%         end
    end
    
    %% SUA MODEL FITTING
    l_d2XT = 500;
    l_L1 = 50;
    ql_L1 = 1;
    ql_d2XT = 50;
    
    stim_params = NIMcreate_stim_params([flen nPix],dt);
    reg_params = NIMcreate_reg_params('lambda_d2XT',l_d2XT,'lambda_L1',l_L1);
    init_reg_params = NIMcreate_reg_params('lambda_d2XT',l_d2XT);
    silent = 1;
    init_optim_params.progTol = 1e-7;
    init_optim_params.optTol = 1e-4;
    init_optim_params.maxIter = 100;
    
    optim_params.progTol = 1e-7;
    optim_params.optTol = 1e-4;
    sub_samp_fac = 10;
    
    poss_spks = 0:10;
    n_symbs = length(poss_spks);
    
    for ss = 1:length(su_probes)
        
        %% for horizontal bars
        cur_used_blocks = hor_blocks(su_used_blocks(hor_blocks,ss)); %blocks when SU
        cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
        cur_used_inds = used_inds(ismember(used_inds,cur_poss_inds));
        NT = length(cur_used_inds);
        if ~isempty(cur_used_inds)
            Robs = all_binned_spikes(cur_used_inds,su_probes(ss));
            fprintf('Fitting unit %d, with %d exc and %d sup squared\n',su_probes(ss),sua_data(ss).npos_stc_hor,sua_data(ss).nneg_stc_hor);
            mod_signs = [1 ones(1,sua_data(ss).npos_stc_hor) -1*ones(1,sua_data(ss).nneg_stc_hor)];
            NL_types = cat(2,'lin', repmat({'quad'},1,sua_data(ss).npos_stc_hor + sua_data(ss).nneg_stc_hor));
            
            fit0 = NIMinitialize_model(stim_params,mod_signs,NL_types,init_reg_params); %initialize NIM
            fit0 = NIMadjust_regularization(fit0,[2:length(mod_signs)],'lambda_d2XT',ql_d2XT);
            
            sub_sample = randperm(NT);
            sub_sample = sub_sample(1:round(NT/sub_samp_fac));
            
            fit0 = NIMfit_filters(fit0,Robs(sub_sample),all_Xmat(cur_used_inds(sub_sample),:),[],[],silent,init_optim_params); %fit stimulus filters
            fit0 = NIMfit_filters(fit0,Robs,all_Xmat(cur_used_inds,:),[],[],silent,init_optim_params); %fit stimulus filters
            
            fit1 = NIMadjust_regularization(fit0,1,'lambda_L1',l_L1);
            fit1 = NIMadjust_regularization(fit1,[2:length(mod_signs)],'lambda_L1',ql_L1);
            fit1 = NIMfit_filters(fit1,Robs,all_Xmat(cur_used_inds,:),[],[],silent,optim_params); %fit stimulus filters
            
            [LL, penLL, pred_rate, G, gint, fgint] = NIMmodel_eval(fit1,Robs,all_Xmat(cur_used_inds,:));
            G = G - fit1.spk_NL_params(1);
            fit2 = NIMfit_logexp_spkNL_withoffset(fit1,Robs,all_Xmat(cur_used_inds,:),[],0,[]);
            ov_spk_nl_params = fit2.spk_NL_params;
            [g_dist,g_x] = ksdensity(G,100);
            ov_pred_resp = logexp_fun(g_x,ov_spk_nl_params);
            
            sua_data(ss).nimfit_hor = fit2;
            sua_data(ss).g_dist_hor = g_dist;
            sua_data(ss).g_x_hor = g_x;
            sua_data(ss).resp_fun_hor = ov_pred_resp;
            
            %ESTIMATE SPK NL
            avg_gsacta = zeros(n_sac_bins,1);
            gsacdep_spk_nl = nan(n_sac_bins,4);
            gsacdep_resp = nan(n_sac_bins,length(g_x));
            avg_msacta = zeros(n_sac_bins,1);
            msacdep_spk_nl = nan(n_sac_bins,4);
            msacdep_resp = nan(n_sac_bins,length(g_x));
            for tt = 1:n_sac_bins
                cur_inds = find(trial_gsac_mat(cur_used_inds,tt) == 1);
                avg_sacta(tt) = mean(Robs(cur_inds));
                cur_fit = NIMfit_logexp_spkNL_withoffset(fit2,Robs(cur_inds),all_Xmat(cur_used_inds(cur_inds),:),[],0,[1 2]);
                sacdep_spk_nl(tt,:) = cur_fit.spk_NL_params;
                sacdep_resp(tt,:) = logexp_fun(g_x,sacdep_spk_nl(tt,:));
            end
            sua_data(ss).sacdep_spk_nl_hor = sacdep_spk_nl;
            sua_data(ss).sac_trigavg_hor = avg_sacta;
            
            %INFO CALS
            avg_rate = mean(Robs);
            n_pts = 100;
            [G_hy,G_hx] = hist(G,n_pts); G_hy = G_hy/sum(G_hy);
            h_pred_resp = logexp_fun(G_hx,ov_spk_nl_params);
            poiss_prob = poisspdf(repmat(poss_spks,n_pts,1),repmat(h_pred_resp',1,n_symbs));
            marg_prob = sum(bsxfun(@times,poiss_prob,G_hy'));
            ov_info = sum(G_hy'.*sum(poiss_prob.*log2(bsxfun(@rdivide,poiss_prob,marg_prob)),2));
            sacdep_info = nan(n_sac_bins,1);
            for tt = 1:n_sac_bins
                cur_pred_resp = logexp_fun(G_hx,sacdep_spk_nl(tt,:));
                poiss_prob = poisspdf(repmat(poss_spks,n_pts,1),repmat(cur_pred_resp',1,n_symbs));
                marg_prob = sum(bsxfun(@times,poiss_prob,G_hy'));
                sacdep_info(tt) = sum(G_hy'.*sum(poiss_prob.*log2(bsxfun(@rdivide,poiss_prob,marg_prob)),2));
            end
            sua_data(ss).ov_info_hor = ov_info;
            sua_data(ss).sacdep_info_hor = sacdep_info;
            
            
            
            %% for vertical bars
            cur_used_blocks = ver_blocks(su_used_blocks(ver_blocks,ss)); %blocks when SU
            cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
            cur_used_inds = used_inds(ismember(used_inds,cur_poss_inds));
            NT = length(cur_used_inds);
            if ~isempty(cur_used_inds)
                Robs = all_binned_spikes(cur_used_inds,su_probes(ss));
                fprintf('Fitting unit %d, with %d exc and %d sup squared\n',su_probes(ss),sua_data(ss).npos_stc_ver,sua_data(ss).nneg_stc_ver);
                mod_signs = [1 ones(1,sua_data(ss).npos_stc_ver) -1*ones(1,sua_data(ss).nneg_stc_ver)];
                NL_types = cat(2,'lin', repmat({'quad'},1,sua_data(ss).npos_stc_ver + sua_data(ss).nneg_stc_ver));
                
                fit0 = NIMinitialize_model(stim_params,mod_signs,NL_types,init_reg_params); %initialize NIM
                fit0 = NIMadjust_regularization(fit0,[2:length(mod_signs)],'lambda_d2XT',ql_d2XT);
                
                sub_sample = randperm(NT);
                sub_sample = sub_sample(1:round(NT/sub_samp_fac));
                
                fit0 = NIMfit_filters(fit0,Robs(sub_sample),all_Xmat(cur_used_inds(sub_sample),:),[],[],silent,init_optim_params); %fit stimulus filters
                fit0 = NIMfit_filters(fit0,Robs,all_Xmat(cur_used_inds,:),[],[],silent,init_optim_params); %fit stimulus filters
                
                fit1 = NIMadjust_regularization(fit0,1,'lambda_L1',l_L1);
                fit1 = NIMadjust_regularization(fit1,[2:length(mod_signs)],'lambda_L1',ql_L1);
                fit1 = NIMfit_filters(fit1,Robs,all_Xmat(cur_used_inds,:),[],[],silent,optim_params); %fit stimulus filters
                
                [LL, penLL, pred_rate, G, gint, fgint] = NIMmodel_eval(fit1,Robs,all_Xmat(cur_used_inds,:));
                G = G - fit1.spk_NL_params(1);
                fit2 = NIMfit_logexp_spkNL_withoffset(fit1,Robs,all_Xmat(cur_used_inds,:),[],0,[]);
                ov_spk_nl_params = fit2.spk_NL_params;
                [g_dist,g_x] = ksdensity(G,100);
                ov_pred_resp = logexp_fun(g_x,ov_spk_nl_params);
                
                sua_data(ss).nimfit_ver = fit2;
                sua_data(ss).g_dist_ver = g_dist;
                sua_data(ss).g_x_ver = g_x;
                sua_data(ss).resp_fun_ver = ov_pred_resp;
                
                %ESTIMATE SPK NL
                avg_sacta = zeros(n_sac_bins,1);
                sacdep_spk_nl = nan(n_sac_bins,4);
                sacdep_resp = nan(n_sac_bins,length(g_x));
                for tt = 1:n_sac_bins
                    cur_inds = find(trial_gsac_mat(cur_used_inds,tt) == 1);
                    avg_sacta(tt) = mean(Robs(cur_inds));
                    cur_fit = NIMfit_logexp_spkNL_withoffset(fit2,Robs(cur_inds),all_Xmat(cur_used_inds(cur_inds),:),[],0,[1 2]);
                    sacdep_spk_nl(tt,:) = cur_fit.spk_NL_params;
                    sacdep_resp(tt,:) = logexp_fun(g_x,sacdep_spk_nl(tt,:));
                end
                sua_data(ss).sacdep_spk_nl_ver = sacdep_spk_nl;
                sua_data(ss).sac_trigavg_ver = avg_sacta;
                
                %INFO CALS
                avg_rate = mean(Robs);
                n_pts = 100;
                [G_hy,G_hx] = hist(G,n_pts); G_hy = G_hy/sum(G_hy);
                h_pred_resp = logexp_fun(G_hx,ov_spk_nl_params);
                poiss_prob = poisspdf(repmat(poss_spks,n_pts,1),repmat(h_pred_resp',1,n_symbs));
                marg_prob = sum(bsxfun(@times,poiss_prob,G_hy'));
                ov_info = sum(G_hy'.*sum(poiss_prob.*log2(bsxfun(@rdivide,poiss_prob,marg_prob)),2));
                sacdep_info = nan(n_sac_bins,1);
                for tt = 1:n_sac_bins
                    cur_pred_resp = logexp_fun(G_hx,sacdep_spk_nl(tt,:));
                    poiss_prob = poisspdf(repmat(poss_spks,n_pts,1),repmat(cur_pred_resp',1,n_symbs));
                    marg_prob = sum(bsxfun(@times,poiss_prob,G_hy'));
                    sacdep_info(tt) = sum(G_hy'.*sum(poiss_prob.*log2(bsxfun(@rdivide,poiss_prob,marg_prob)),2));
                end
                sua_data(ss).ov_info_ver = ov_info;
                sua_data(ss).sacdep_info_ver = sacdep_info;
            end
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
    anal_params.trial_lags = trial_lags;
    anal_params.acorr_lags = acorr_lags;
    anal_params.sac_bin_cents = sac_bin_cents;
    cd(save_dir)
    sname = 'full_sacmod_data';
    save(sname,'sua_data','su_probes','mua_data','anal_params','gen_data');
    clear sua_data su_probes mua_data gen_data
end