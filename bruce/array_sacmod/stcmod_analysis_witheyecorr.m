clear all
% close all

dir_prefix = '~';
% dir_prefix = '/Volumes/james';

Expt_num = 86;

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
use_nPix = 24;
full_nPix = 36;
stim_params = NIMcreate_stim_params([flen full_nPix],dt);
Fr = 1;
trial_dur = 4;
bar_ori = 0;

anal_name = 'hbar_stc_sacmod_analysis';

%% LOAD OVERALL SU DATA
load ~/Analysis/bruce/summary_analysis/su_data.mat


clear expt* included_type

Expt_name = sprintf('G0%d',Expt_num);
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
anal_dir = ['~/Analysis/bruce/' Expt_name];
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

%%
if strcmp(Expt_name,'G093')
    include_expts = {'rls.FaXwi','rls.FaXwiXimi'};
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
% cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1);
cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori);

if strcmp(Expt_name,'G087')
    cur_block_set(cur_block_set == 15) = []; %only 6 trials and causes problems
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
cur_expt_id = find(su_data.expt_nums == Expt_num);
su_probes = find(su_data.is_su(cur_expt_id,:));
mua_probes = setdiff(1:96,su_probes); %probes with ONLY MU
aclust_probenums = [autoclust(cur_block_set(1),:).probenum];
autoclust = autoclust(:,ismember(aclust_probenums,su_probes));

%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_stim_times = [];
all_Xmat = [];
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
        fprintf('Expt %d Block %d of %d; grayback GS, ori:%d\n',Expt_num,ee,length(cur_block_set),expt_bar_ori(cur_block_set(ee)));
    elseif ismember(ee,imback_gs_expts)
        fprintf('Expt %d Block %d of %d; imback GS, ori:%d\n',Expt_num,ee,length(cur_block_set),expt_bar_ori(cur_block_set(ee)));
    elseif ismember(ee,sim_sac_expts)
        fprintf('Expt %d Block %d of %d; SimSac, ori:%d\n',Expt_num,ee,length(cur_block_set),expt_bar_ori(cur_block_set(ee)));
    else
        fprintf('Expt %d Block %d of %d;  UNMATCHED EXPT TYPE\n',Expt_num,ee,length(cur_block_set));
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
    buffer_pix = floor((expt_npix(cur_block) - full_nPix)/2);
    cur_use_pix = (1:full_nPix) + buffer_pix;
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start'/1e4;
        n_frames = size(left_stim_mats{use_trials(tt)},1);
        if n_frames > 0
            if length(cur_stim_times) == 1
                cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
                cur_stim_times(cur_stim_times > trial_end_times(use_trials(tt))) = [];
                cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
            end
        end
        %             cur_t_edges = (trial_start_times(use_trials(tt)):dt:trial_end_times(use_trials(tt)))';
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > min_trial_dur/dt
            use_frames = min(length(cur_stim_times),n_frames);
            cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            bar_Xmat = create_time_embedding(cur_stim_mat,stim_params);
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times];
            all_t_axis = [all_t_axis; cur_t_axis];
            all_Xmat = [all_Xmat; bar_Xmat];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
    end
    trial_cnt = trial_cnt + n_trials;
end

%% BIN SPIKES FOR MU AND SU
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
%% SMOOTH AND NORMALIZE BINNED SPIKE DATA
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
            if ~isempty(cur_block_inds)
                all_spike_rate(cur_block_inds,su_probes(ss)) = jmm_smooth_1d_cor(all_binned_spikes(cur_block_inds,su_probes(ss)),sua_sm_sig);
            end
        end
    end
end

%normalized firing rates (smoothed)
for ee = 1:length(cur_block_set)
    cur_block_inds = used_inds(all_blockvec(used_inds) == ee);
    if ~isempty(cur_block_inds)
        all_spike_rate_norm(cur_block_inds,:) = bsxfun(@rdivide,all_spike_rate(cur_block_inds,:),...
            block_mean_rates(ee,:));
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
    fprintf('Loading ET data for expt %d, block %d of %d\n',Expt_num,ee,length(cur_block_set));
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

%% DETECT SACCADES
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

sac_dbuff = round(0.05/eye_dt);
sac_pre_pos = nan(length(saccade_inds),4);
sac_post_pos = nan(length(saccade_inds),4);
for ii = 1:length(saccade_inds)
    pre_inds = (sac_start_inds(ii) - sac_dbuff):sac_start_inds(ii);
    pre_inds(pre_inds < 1) = 1;
    sac_pre_pos(ii,:) = median(all_eye_vals(pre_inds,:),1);
    post_inds = sac_stop_inds(ii):(sac_stop_inds(ii) + sac_dbuff);
    post_inds(post_inds > length(all_eye_ts)) = length(all_eye_ts);
    sac_post_pos(ii,:) = median(all_eye_vals(post_inds,:),1);
end

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
sac_bin_width = 1; %number of dt bins
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

time_to_gsac = nan(length(all_t_axis),1);
time_to_msac = nan(length(all_t_axis),1);
max_dist = round(0.5/dt);
rel_tax = (-max_dist:max_dist)';
for ii = 1:length(used_gsacs)
    cur_inds = (gsac_start_inds(ii) - max_dist):(gsac_start_inds(ii) + max_dist);
    cur_uset = find(cur_inds > 1 & cur_inds < length(all_t_axis));
    cur_uset(abs(time_to_gsac(cur_inds(cur_uset))) > abs(rel_tax(cur_uset))) = [];
    time_to_gsac(cur_inds(cur_uset)) = rel_tax(cur_uset);
end
for ii = 1:length(used_msacs)
    cur_inds = (msac_start_inds(ii) - max_dist):(msac_start_inds(ii) + max_dist);
    cur_uset = find(cur_inds > 1 & cur_inds < length(all_t_axis));
    cur_uset(abs(time_to_gsac(cur_inds(cur_uset))) > abs(rel_tax(cur_uset))) = [];
    time_to_msac(cur_inds(cur_uset)) = rel_tax(cur_uset);
end
time_to_gsac = time_to_gsac*dt;
time_to_msac = time_to_msac*dt;

tsig = 1.5; gsig = 2;
[TT,GG] = meshgrid(-2*tsig:2*tsig,-2*gsig:2*gsig);
DD = TT.^2/(2*tsig^2) + GG.^2/(2*gsig^2);
kkern = exp(-DD.^2); kkern = kkern/sum(kkern(:));

%% BLOCK PREDICTOR
Xblock = zeros(length(all_stim_times),length(cur_block_set));
for i = 1:length(cur_block_set)
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end

%% LOAD INFERRED EYE POSITIONS
if bar_ori == 0
    eye_cor_data_name = 'monoc_eyecorr_hbar';
elseif bar_ori == 90
    eye_cor_data_name = 'monoc_eyecorr_vbar';
end
cd(anal_dir)
load(eye_cor_data_name);

inferred_eyepos_cor = it_fshift_cor{end};
interp_eyepos_cor = round(interp1(inferred_eyepos_taxis,inferred_eyepos_cor,all_t_axis(used_inds)));
used_inds(isnan(interp_eyepos_cor)) = [];
interp_eyepos_cor(isnan(interp_eyepos_cor)) = [];

NT = length(used_inds);
resh_X = reshape(all_Xmat(used_inds,:)',[flen full_nPix NT]);
for ii = 1:NT
    resh_X(:,:,ii) = shift_matrix_Nd(resh_X(:,:,ii), -interp_eyepos_cor(ii),2);
end
X_sh = reshape(resh_X,flen*full_nPix,NT)';

[Xinds,Tinds] = meshgrid(1:full_nPix,1:flen);
buffer_pix = floor((full_nPix - use_nPix)/2);
cur_use_pix = (1:use_nPix) + buffer_pix;
use_kInds = find(ismember(Xinds(:),cur_use_pix));

Xmat_raw = all_Xmat(used_inds,use_kInds);
Xmat_sh = X_sh(:,use_kInds);
clear X_sh

%% SUA STA/STC ANALYSIS
addpath('~/James_scripts/TentBasis2D/')

back_t = 0.2;
for_t = 0.5;
n_Gbins = 35;
Xtick = -(back_t-dt/2):(dt):(for_t+dt/2);
n_sbins = length(Xtick);

nneg = 5; npos = 5;
stc_thresh = -1.5e-3;
nonpar_bins = 40;
if bar_ori == 0
    cur_blocks = find(expt_bar_ori(cur_block_set) == 0);
else
    cur_blocks = find(expt_bar_ori(cur_block_set) == 90);
end

stim_params = NMMcreate_stim_params([flen use_nPix],dt);
stim_params(2) = NMMcreate_stim_params(length(cur_block_set),dt);
max_squared_filts = 3;
lambda_d2XT = 500;
for ss = 1:length(su_probes)
    
    fprintf('Computing STC for SU %d of %d\n',ss,length(su_probes));
    
    cur_used_blocks = cur_blocks(su_used_blocks(cur_blocks,ss)); %blocks when SU
    cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
    cur_used_inds = find(ismember(used_inds,cur_poss_inds));
    
    if ~isempty(cur_used_inds)
        Robs = all_binned_spikes(used_inds(cur_used_inds),su_probes(ss));
        Robs_sm = all_spike_rate_norm(used_inds(cur_used_inds),su_probes(ss));
        avg_rate = mean(Robs);
        avg_rate_sm = mean(Robs_sm);
        sua_data(ss).block_avg_rates = block_mean_rates(cur_used_blocks,su_probes(ss));
        sua_data(ss).avg_rate = mean(Robs);
        sua_data(ss).nspks = sum(Robs);
        cur_n_grayback_blocks = sum(ismember(cur_used_blocks,grayback_gs_expts));
        cur_n_imback_blocks = sum(ismember(cur_used_blocks,imback_gs_expts));
        sua_data(ss).nblocks = [cur_n_grayback_blocks cur_n_imback_blocks];
        cur_n_gsacs = sum(ismember(gsac_start_inds,cur_used_inds));
        cur_n_msacs = sum(ismember(gsac_start_inds,cur_used_inds));
        sua_data(ss).nsacs = [cur_n_gsacs cur_n_msacs];
        sua_data(ss).stc_use = true;
        
        %% UNCOR DATA
        spikebins = convert_to_spikebins(Robs);
        spike_cond_stim = Xmat_raw(cur_used_inds(spikebins),:);
        sta = mean(spike_cond_stim) - mean(Xmat_raw(cur_used_inds,:));
        sta = sta/norm(sta);
        %         proj_mat = sta'/(sta*sta')*sta;
        %         stim_proj = Xmat_raw(cur_used_inds,:) - Xmat_raw(cur_used_inds,:)*proj_mat;
        stim_proj = Xmat_raw(cur_used_inds,:);
        stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
        [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
        stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
        evals = diag(evals);
        sua_data(ss).sta = sta;
        sua_data(ss).stcs = stcs;
        sua_data(ss).evals = evals;
        cur_eval_diff = diff(flipud(sua_data(ss).evals));
        sua_data(ss).npos_stc = min(5,find(cur_eval_diff > stc_thresh,1,'first')-1);
        sua_data(ss).nneg_stc = min(5,find(cur_eval_diff(end:-1:1) > stc_thresh,1,'first')-1);
        
        %fit GLM
        g_sta = Xmat_raw(cur_used_inds,:)*sua_data(ss).sta';
        g_posstc = zeros(length(cur_used_inds),sua_data(ss).npos_stc);
        g_negstc = zeros(length(cur_used_inds),sua_data(ss).nneg_stc);
        for gg = 1:sua_data(ss).npos_stc
            g_posstc(:,gg) = Xmat_raw(cur_used_inds,:)*sua_data(ss).stcs(:,gg);
        end
        for gg = 1:sua_data(ss).nneg_stc
            g_negstc(:,gg) =  Xmat_raw(cur_used_inds,:)*sua_data(ss).stcs(:,end-nneg+gg);
        end
        cur_Xmat = [g_sta g_posstc.^2 g_negstc.^2 Xblock(cur_used_inds,:)];
        [fitmod] = regGLM_fit(cur_Xmat,Robs,[],[],[],[],1);
        [~, ~, ~, G] = regGLM_eval(fitmod,Robs,cur_Xmat);
        G = G - fitmod.theta;
        [g_dist,g_x] = ksdensity(G); g_dist = g_dist/sum(g_dist);
        sua_data(ss).stc_glm = fitmod;
        
        TB_stim = [time_to_gsac(used_inds(cur_used_inds)) G];
        %         Ytick = linspace(my_prctile(TB_stim(:,2),0.5),my_prctile(TB_stim(:,2),100-0.5),n_Gbins);
        Ytick = linspace(my_prctile(TB_stim(:,2),0.05),my_prctile(TB_stim(:,2),100-0.05),n_Gbins);
        TB = TentBasis2D(Xtick, Ytick);
        used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
            TB_stim(:,2) >= Ytick(1) & TB_stim(:,2) <= Ytick(end));
        [TB_Xmat,TB_counts] = TB.InputNL2D(TB_stim(used_data,:));
        L2_params = create_L2_params([],[1 n_sbins*n_Gbins],[n_sbins n_Gbins],2,3,[Inf Inf],[0.1 1]);
        TB_fitmod = regGLM_fit(TB_Xmat,Robs(used_data),L2_params,10,[],[],1);
        TB_K = reshape(TB_fitmod.K,n_sbins,n_Gbins)';
        bin_areas = TB.GetBinAreas();
        gsac_TB_dist = TB_counts./bin_areas;
        gsac_TB_dist = gsac_TB_dist'/sum(gsac_TB_dist(:));
        gsac_TB_rate = log(1 + exp(TB_K + TB_fitmod.theta));
        
        %INFO CALS
        cur_avg_rate = mean(Robs(used_data));
        marg_gdist = sum(gsac_TB_dist,2);
        marg_sdist = sum(gsac_TB_dist);
        marg_gsacrate = sum(gsac_TB_dist.*gsac_TB_rate)./marg_sdist;
        marg_grate = sum(gsac_TB_dist.*gsac_TB_rate,2)./marg_gdist;
        ov_info_gsac = sum(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate));
        gsacdep_info = nan(n_sac_bins,1);
        for tt = 1:n_sbins
            gsacdep_info(tt) = sum(gsac_TB_dist(:,tt).*gsac_TB_rate(:,tt).*log2(gsac_TB_rate(:,tt)/marg_gsacrate(tt)))/sum(gsac_TB_dist(:,tt));
        end
        sua_data(ss).gsac_TB_dist = gsac_TB_dist;
        sua_data(ss).gsac_TB_rate = gsac_TB_rate;
        sua_data(ss).Gtick = Ytick;
        sua_data(ss).Stick = Xtick;
        sua_data(ss).ov_info_gsac = ov_info_gsac;
        sua_data(ss).gsacdep_info = gsacdep_info;
        sua_data(ss).marg_gsacrate = marg_gsacrate;
        
        %         TB_stim = [time_to_msac(used_inds(cur_used_inds)) G];
        % %         Ytick = linspace(my_prctile(TB_stim(:,2),0.5),my_prctile(TB_stim(:,2),100-0.5),n_Gbins);
        %         Ytick = linspace(my_prctile(TB_stim(:,2),0.05),my_prctile(TB_stim(:,2),100-0.05),n_Gbins);
        %         TB = TentBasis2D(Xtick, Ytick);
        %         used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
        %             TB_stim(:,2) >= Ytick(1) & TB_stim(:,2) <= Ytick(end));
        %         [TB_Xmat,TB_counts] = TB.InputNL2D(TB_stim(used_data,:));
        %         L2_params = create_L2_params([],[1 n_sbins*n_Gbins],[n_sbins n_Gbins],2,3,[Inf Inf],[0.1 1]);
        %         TB_fitmod = regGLM_fit(TB_Xmat,Robs(used_data),L2_params,10,[],[],1);
        %         TB_K = reshape(TB_fitmod.K,n_sbins,n_Gbins)';
        %         bin_areas = TB.GetBinAreas();
        %         msac_TB_dist = TB_counts./bin_areas;
        %         msac_TB_dist = msac_TB_dist'/sum(msac_TB_dist(:));
        %         msac_TB_rate = log(1 + exp(TB_K + TB_fitmod.theta));
        %
        %         %INFO CALS
        %         cur_avg_rate = mean(Robs(used_data));
        %         marg_gdist = sum(msac_TB_dist,2);
        %         marg_sdist = sum(msac_TB_dist);
        %         marg_msacrate = sum(msac_TB_dist.*msac_TB_rate)./marg_sdist;
        %         marg_grate = sum(msac_TB_dist.*msac_TB_rate,2)./marg_gdist;
        %         ov_info_msac = sum(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate));
        %         msacdep_info = nan(n_sac_bins,1);
        %         for tt = 1:n_sbins
        %             msacdep_info(tt) = sum(msac_TB_dist(:,tt).*msac_TB_rate(:,tt).*log2(msac_TB_rate(:,tt)/marg_msacrate(tt)))/sum(msac_TB_dist(:,tt));
        %         end
        %         sua_data(ss).msac_TB_dist = msac_TB_dist;
        %         sua_data(ss).msac_TB_rate = msac_TB_rate;
        %         sua_data(ss).Gtick = Ytick;
        %         sua_data(ss).Stick = Xtick;
        %         sua_data(ss).ov_info_msac = ov_info_msac;
        %         sua_data(ss).msacdep_info = msacdep_info;
        
        %% EYE-COR DATA
        spikebins = convert_to_spikebins(Robs);
        spike_cond_stim = Xmat_sh(cur_used_inds(spikebins),:);
        sta_cor = mean(spike_cond_stim) - mean(Xmat_sh(cur_used_inds,:));
        sta_cor = sta_cor/norm(sta_cor);
        %         proj_mat = sta_cor'/(sta_cor*sta_cor')*sta_cor;
        %         stim_proj = Xmat_sh(cur_used_inds,:) - Xmat_sh(cur_used_inds,:)*proj_mat;
        stim_proj = Xmat_sh(cur_used_inds,:);
        stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
        [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
        stcs_cor = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]);
        stcs_cor = stcs_cor(:,end:-1:1);
        evals_cor = diag(evals);
        sua_data(ss).sta_cor = sta_cor;
        sua_data(ss).stcs_cor = stcs_cor;
        sua_data(ss).evals_cor = evals_cor;
        cur_eval_diff = diff(flipud(sua_data(ss).evals_cor));
        sua_data(ss).npos_stc_cor = min(5,find(cur_eval_diff > stc_thresh,1,'first')-1);
        sua_data(ss).nneg_stc_cor = min(5,find(cur_eval_diff(end:-1:1) > stc_thresh,1,'first')-1);
        
        %fit GLM
        g_sta = Xmat_sh(cur_used_inds,:)*sua_data(ss).sta_cor';
        g_posstc = zeros(length(cur_used_inds),sua_data(ss).npos_stc_cor);
        g_negstc = zeros(length(cur_used_inds),sua_data(ss).nneg_stc_cor);
        for gg = 1:sua_data(ss).npos_stc_cor
            g_posstc(:,gg) = Xmat_sh(cur_used_inds,:)*sua_data(ss).stcs_cor(:,gg);
        end
        for gg = 1:sua_data(ss).nneg_stc_cor
            g_negstc(:,gg) =  Xmat_sh(cur_used_inds,:)*sua_data(ss).stcs_cor(:,end-nneg+gg);
        end
        cur_Xmat = [g_sta g_posstc.^2 g_negstc.^2 Xblock(cur_used_inds,:)];
        [fitmod] = regGLM_fit(cur_Xmat,Robs,[],[],[],[],1);
        [~, ~, ~, G] = regGLM_eval(fitmod,Robs,cur_Xmat);
        G = G - fitmod.theta;
        [g_dist,g_x] = ksdensity(G); g_dist = g_dist/sum(g_dist);
        sua_data(ss).stc_glm_cor = fitmod;
        
%         %FIT NMM
%         cur_nsq_filts = min(max_squared_filts,sua_data(ss).npos_stc_cor);
%         reg_params = NMMcreate_reg_params('lambda_d2XT',[lambda_d2XT ones(1,cur_nsq_filts)*sqrt(lambda_d2XT) 0]);
%         mod_signs = ones(1,cur_nsq_filts+1+1);
%         NL_types = [{'lin'} repmat({'quad'},1,cur_nsq_filts) {'lin'}];
%         Xtargs = [ones(1,cur_nsq_filts+1) 2];
%         Xcell{1} = Xmat_sh(cur_used_inds,:);
%         Xcell{2} = Xblock(used_inds(cur_used_inds,:),:);
%         init_filts = cell(length(mod_signs),1);
%         init_filts{end} = zeros(length(cur_block_set),1);
%         fit0 = NMMinitialize_model(stim_params,mod_signs,NL_types,reg_params,Xtargs,init_filts);
%         fit0 = NMMfit_filters(fit0,Robs,Xcell,[],[],0);
%         [LL,~,~,nmmG] = NMMmodel_eval(fit0,Robs,Xcell);
        
        TB_stim = [time_to_gsac(used_inds(cur_used_inds)) G];
%         TB_stim = [time_to_gsac(used_inds(cur_used_inds)) nmmG];
%         TB_stim = [time_to_gsac(used_inds(cur_used_inds)) g_sta];
%         TB_stim = [time_to_gsac(used_inds(cur_used_inds)) sum(g_posstc.^2,2)];
%         TB_stim = [time_to_gsac(used_inds(cur_used_inds)) sum(g_negstc.^2,2)];
        %         Ytick = linspace(my_prctile(TB_stim(:,2),0.5),my_prctile(TB_stim(:,2),100-0.5),n_Gbins);
        Ytick = linspace(my_prctile(TB_stim(:,2),0.05),my_prctile(TB_stim(:,2),100-0.05),n_Gbins);
        TB = TentBasis2D(Xtick, Ytick);
        used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
            TB_stim(:,2) >= Ytick(1) & TB_stim(:,2) <= Ytick(end));
        [TB_Xmat,TB_counts] = TB.InputNL2D(TB_stim(used_data,:));
        L2_params = create_L2_params([],[1 n_sbins*n_Gbins],[n_sbins n_Gbins],2,3,[Inf Inf],[0.1 1]);
        TB_fitmod = regGLM_fit(TB_Xmat,Robs(used_data),L2_params,10,[],[],1);
        TB_K = reshape(TB_fitmod.K,n_sbins,n_Gbins)';
        bin_areas = TB.GetBinAreas();
        gsac_TB_dist = TB_counts./bin_areas;
        gsac_TB_dist = gsac_TB_dist'/sum(gsac_TB_dist(:));
        gsac_TB_rate = log(1 + exp(TB_K + TB_fitmod.theta));
        
        %INFO CALS
        cur_avg_rate = mean(Robs(used_data));
        marg_gdist = sum(gsac_TB_dist,2);
        marg_sdist = sum(gsac_TB_dist);
        marg_gsacrate = sum(gsac_TB_dist.*gsac_TB_rate)./marg_sdist;
        marg_grate = sum(gsac_TB_dist.*gsac_TB_rate,2)./marg_gdist;
        ov_info_gsac = sum(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate));
        gsacdep_info = nan(n_sac_bins,1);
        for tt = 1:n_sbins
            gsacdep_info(tt) = sum(gsac_TB_dist(:,tt).*gsac_TB_rate(:,tt).*log2(gsac_TB_rate(:,tt)/marg_gsacrate(tt)))/sum(gsac_TB_dist(:,tt));
        end
        sua_data(ss).gsac_TB_dist_cor = gsac_TB_dist;
        sua_data(ss).gsac_TB_rate_cor = gsac_TB_rate;
        sua_data(ss).Gtick_cor = Ytick;
        sua_data(ss).Stick_cor = Xtick;
        sua_data(ss).ov_info_gsac_cor = ov_info_gsac;
        sua_data(ss).gsacdep_info_cor = gsacdep_info;
        sua_data(ss).marg_gsacrate_cor = marg_gsacrate;
        sua_data(ss).marg_gdist_cor = marg_gdist;
        sua_data(ss).marg_gsac_dist_cor = marg_sdist;
        
        TB_stim = [time_to_msac(used_inds(cur_used_inds)) G];
%                 Ytick = linspace(my_prctile(TB_stim(:,2),0.5),my_prctile(TB_stim(:,2),100-0.5),n_Gbins);
        Ytick = linspace(my_prctile(TB_stim(:,2),0.05),my_prctile(TB_stim(:,2),100-0.05),n_Gbins);
        TB = TentBasis2D(Xtick, Ytick);
        used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
            TB_stim(:,2) >= Ytick(1) & TB_stim(:,2) <= Ytick(end));
        [TB_Xmat,TB_counts] = TB.InputNL2D(TB_stim(used_data,:));
        L2_params = create_L2_params([],[1 n_sbins*n_Gbins],[n_sbins n_Gbins],2,3,[Inf Inf],[0.1 1]);
        TB_fitmod = regGLM_fit(TB_Xmat,Robs(used_data),L2_params,10,[],[],1);
        TB_K = reshape(TB_fitmod.K,n_sbins,n_Gbins)';
        bin_areas = TB.GetBinAreas();
        msac_TB_dist = TB_counts./bin_areas;
        msac_TB_dist = msac_TB_dist'/sum(msac_TB_dist(:));
        msac_TB_rate = log(1 + exp(TB_K + TB_fitmod.theta));
        
        %INFO CALS
        cur_avg_rate = mean(Robs(used_data));
        marg_gdist = sum(msac_TB_dist,2);
        marg_sdist = sum(msac_TB_dist);
        marg_msacrate = sum(msac_TB_dist.*msac_TB_rate)./marg_sdist;
        marg_grate = sum(msac_TB_dist.*msac_TB_rate,2)./marg_gdist;
        ov_info_msac = sum(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate));
        msacdep_info = nan(n_sac_bins,1);
        for tt = 1:n_sbins
            msacdep_info(tt) = sum(msac_TB_dist(:,tt).*msac_TB_rate(:,tt).*log2(msac_TB_rate(:,tt)/marg_msacrate(tt)))/sum(msac_TB_dist(:,tt));
        end
        sua_data(ss).msac_TB_dist_cor = msac_TB_dist;
        sua_data(ss).msac_TB_rate_cor = msac_TB_rate;
        sua_data(ss).Gtick_cor = Ytick;
        sua_data(ss).Stick_cor = Xtick;
        sua_data(ss).ov_info_msac_cor = ov_info_msac;
        sua_data(ss).msacdep_info_cor = msacdep_info;
        sua_data(ss).marg_msacrate_cor = marg_msacrate;
        sua_data(ss).marg_gdist_msac_cor = marg_gdist;
        sua_data(ss).marg_msac_dist_cor = marg_sdist;
        
    else
        sua_data(ss).stc_hor_use = false;
    end
end


%% MUA STA/STC ANALYSIS
for cc = 1:96
    
    fprintf('Computing STC for MU %d of %d\n',cc,96);
    
    fprintf('Computing trig avgs for MUA %d of %d\n',cc,96);
    su_probe_ind = find(su_probes == cc);
    if ~isempty(su_probe_ind)
        cur_used_blocks = find(su_used_blocks(:,su_probe_ind) == 0); %blocks when NO SU
    else
        cur_used_blocks = 1:length(cur_block_set);
    end
    
    cur_used_blocks = cur_used_blocks(ismember(cur_used_blocks,cur_blocks)); %blocks when SU
    cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
    cur_used_inds = find(ismember(used_inds,cur_poss_inds));
    
    if ~isempty(cur_used_inds)
        Robs = all_binned_spikes(used_inds(cur_used_inds),cc);
        Robs_sm = all_spike_rate_norm(used_inds(cur_used_inds),cc);
        avg_rate = mean(Robs);
        avg_rate_sm = mean(Robs_sm);
        mua_data(cc).block_avg_rates = block_mean_rates(cur_used_blocks,cc);
        mua_data(cc).avg_rate = mean(Robs);
        mua_data(cc).nspks = sum(Robs);
        cur_n_grayback_blocks = sum(ismember(cur_used_blocks,grayback_gs_expts));
        cur_n_imback_blocks = sum(ismember(cur_used_blocks,imback_gs_expts));
        mua_data(cc).nblocks = [cur_n_grayback_blocks cur_n_imback_blocks];
        cur_n_gsacs = sum(ismember(gsac_start_inds,cur_used_inds));
        cur_n_msacs = sum(ismember(gsac_start_inds,cur_used_inds));
        mua_data(cc).nsacs = [cur_n_gsacs cur_n_msacs];
        mua_data(cc).stc_use = true;
        
        %% UNCOR DATA
        spikebins = convert_to_spikebins(Robs);
        spike_cond_stim = Xmat_raw(cur_used_inds(spikebins),:);
        sta = mean(spike_cond_stim) - mean(Xmat_raw(cur_used_inds,:));
        sta = sta/norm(sta);
        %         proj_mat = sta'/(sta*sta')*sta;
        %         stim_proj = Xmat_raw(cur_used_inds,:) - Xmat_raw(cur_used_inds,:)*proj_mat;
        stim_proj = Xmat_raw(cur_used_inds,:);
        stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
        [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
        stcs  = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]); stcs  = stcs(:,end:-1:1);
        evals = diag(evals);
        mua_data(cc).sta = sta;
        mua_data(cc).stcs = stcs;
        mua_data(cc).evals = evals;
        cur_eval_diff = diff(flipud(mua_data(cc).evals));
        mua_data(cc).npos_stc = min(5,find(cur_eval_diff > stc_thresh,1,'first')-1);
        mua_data(cc).nneg_stc = min(5,find(cur_eval_diff(end:-1:1) > stc_thresh,1,'first')-1);
        
        %fit GLM
        g_sta = Xmat_raw(cur_used_inds,:)*mua_data(cc).sta';
        g_posstc = zeros(length(cur_used_inds),mua_data(cc).npos_stc);
        g_negstc = zeros(length(cur_used_inds),mua_data(cc).nneg_stc);
        for gg = 1:mua_data(cc).npos_stc
            g_posstc(:,gg) = Xmat_raw(cur_used_inds,:)*mua_data(cc).stcs(:,gg);
        end
        for gg = 1:mua_data(cc).nneg_stc
            g_negstc(:,gg) =  Xmat_raw(cur_used_inds,:)*mua_data(cc).stcs(:,end-nneg+gg);
        end
        cur_Xmat = [g_sta g_posstc.^2 g_negstc.^2 Xblock(cur_used_inds,:)];
        [fitmod] = regGLM_fit(cur_Xmat,Robs,[],[],[],[],1);
        [~, ~, ~, G] = regGLM_eval(fitmod,Robs,cur_Xmat);
        G = G - fitmod.theta;
        [g_dist,g_x] = ksdensity(G); g_dist = g_dist/sum(g_dist);
        mua_data(cc).stc_glm = fitmod;
        
        TB_stim = [time_to_gsac(used_inds(cur_used_inds)) G];
        %         Ytick = linspace(my_prctile(TB_stim(:,2),0.5),my_prctile(TB_stim(:,2),100-0.5),n_Gbins);
        Ytick = linspace(my_prctile(TB_stim(:,2),0.05),my_prctile(TB_stim(:,2),100-0.05),n_Gbins);
        TB = TentBasis2D(Xtick, Ytick);
        used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
            TB_stim(:,2) >= Ytick(1) & TB_stim(:,2) <= Ytick(end));
        [TB_Xmat,TB_counts] = TB.InputNL2D(TB_stim(used_data,:));
        L2_params = create_L2_params([],[1 n_sbins*n_Gbins],[n_sbins n_Gbins],2,3,[Inf Inf],[0.1 1]);
        TB_fitmod = regGLM_fit(TB_Xmat,Robs(used_data),L2_params,10,[],[],1);
        TB_K = reshape(TB_fitmod.K,n_sbins,n_Gbins)';
        bin_areas = TB.GetBinAreas();
        gsac_TB_dist = TB_counts./bin_areas;
        gsac_TB_dist = gsac_TB_dist'/sum(gsac_TB_dist(:));
        gsac_TB_rate = log(1 + exp(TB_K + TB_fitmod.theta));
        
        %INFO CALS
        cur_avg_rate = mean(Robs(used_data));
        marg_gdist = sum(gsac_TB_dist,2);
        marg_sdist = sum(gsac_TB_dist);
        marg_gsacrate = sum(gsac_TB_dist.*gsac_TB_rate)./marg_sdist;
        marg_grate = sum(gsac_TB_dist.*gsac_TB_rate,2)./marg_gdist;
        ov_info_gsac = sum(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate));
        gsacdep_info = nan(n_sac_bins,1);
        for tt = 1:n_sbins
            gsacdep_info(tt) = sum(gsac_TB_dist(:,tt).*gsac_TB_rate(:,tt).*log2(gsac_TB_rate(:,tt)/marg_gsacrate(tt)))/sum(gsac_TB_dist(:,tt));
        end
        mua_data(cc).gsac_TB_dist = gsac_TB_dist;
        mua_data(cc).gsac_TB_rate = gsac_TB_rate;
        mua_data(cc).Gtick = Ytick;
        mua_data(cc).Stick = Xtick;
        mua_data(cc).ov_info_gsac = ov_info_gsac;
        mua_data(cc).gsacdep_info = gsacdep_info;
        mua_data(cc).marg_gsacrate = marg_gsacrate;
        
        %         TB_stim = [time_to_msac(used_inds(cur_used_inds)) G];
        % %         Ytick = linspace(my_prctile(TB_stim(:,2),0.5),my_prctile(TB_stim(:,2),100-0.5),n_Gbins);
        %         Ytick = linspace(my_prctile(TB_stim(:,2),0.05),my_prctile(TB_stim(:,2),100-0.05),n_Gbins);
        %         TB = TentBasis2D(Xtick, Ytick);
        %         used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
        %             TB_stim(:,2) >= Ytick(1) & TB_stim(:,2) <= Ytick(end));
        %         [TB_Xmat,TB_counts] = TB.InputNL2D(TB_stim(used_data,:));
        %         L2_params = create_L2_params([],[1 n_sbins*n_Gbins],[n_sbins n_Gbins],2,3,[Inf Inf],[0.1 1]);
        %         TB_fitmod = regGLM_fit(TB_Xmat,Robs(used_data),L2_params,10,[],[],1);
        %         TB_K = reshape(TB_fitmod.K,n_sbins,n_Gbins)';
        %         bin_areas = TB.GetBinAreas();
        %         msac_TB_dist = TB_counts./bin_areas;
        %         msac_TB_dist = msac_TB_dist'/sum(msac_TB_dist(:));
        %         msac_TB_rate = log(1 + exp(TB_K + TB_fitmod.theta));
        %
        %         %INFO CALS
        %         cur_avg_rate = mean(Robs(used_data));
        %         marg_gdist = sum(msac_TB_dist,2);
        %         marg_sdist = sum(msac_TB_dist);
        %         marg_msacrate = sum(msac_TB_dist.*msac_TB_rate)./marg_sdist;
        %         marg_grate = sum(msac_TB_dist.*msac_TB_rate,2)./marg_gdist;
        %         ov_info_msac = sum(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate));
        %         msacdep_info = nan(n_sac_bins,1);
        %         for tt = 1:n_sbins
        %             msacdep_info(tt) = sum(msac_TB_dist(:,tt).*msac_TB_rate(:,tt).*log2(msac_TB_rate(:,tt)/marg_msacrate(tt)))/sum(msac_TB_dist(:,tt));
        %         end
        %         mua_data(cc).msac_TB_dist = msac_TB_dist;
        %         mua_data(cc).msac_TB_rate = msac_TB_rate;
        %         mua_data(cc).Gtick = Ytick;
        %         mua_data(cc).Stick = Xtick;
        %         mua_data(cc).ov_info_msac = ov_info_msac;
        %         mua_data(cc).msacdep_info = msacdep_info;
        
        %% EYE-COR DATA
        spikebins = convert_to_spikebins(Robs);
        spike_cond_stim = Xmat_sh(cur_used_inds(spikebins),:);
        sta_cor = mean(spike_cond_stim) - mean(Xmat_sh(cur_used_inds,:));
        sta_cor = sta_cor/norm(sta_cor);
        %         proj_mat = sta_cor'/(sta_cor*sta_cor')*sta_cor;
        %         stim_proj = Xmat_sh(cur_used_inds,:) - Xmat_sh(cur_used_inds,:)*proj_mat;
        stim_proj = Xmat_sh(cur_used_inds,:);
        stvcv = cov(stim_proj(spikebins,:));  utvcv = cov(stim_proj);
        [evecs,evals] = eig(stvcv-utvcv); evs   = diag(evals);
        stcs_cor = evecs(:,[nneg:-1:1,length(evs)-npos+1:end]);
        stcs_cor = stcs_cor(:,end:-1:1);
        evals_cor = diag(evals);
        mua_data(cc).sta_cor = sta_cor;
        mua_data(cc).stcs_cor = stcs_cor;
        mua_data(cc).evals_cor = evals_cor;
        cur_eval_diff = diff(flipud(mua_data(cc).evals_cor));
        mua_data(cc).npos_stc_cor = min(5,find(cur_eval_diff > stc_thresh,1,'first')-1);
        mua_data(cc).nneg_stc_cor = min(5,find(cur_eval_diff(end:-1:1) > stc_thresh,1,'first')-1);
        
        %fit GLM
        g_sta = Xmat_sh(cur_used_inds,:)*mua_data(cc).sta_cor';
        g_posstc = zeros(length(cur_used_inds),mua_data(cc).npos_stc_cor);
        g_negstc = zeros(length(cur_used_inds),mua_data(cc).nneg_stc_cor);
        for gg = 1:mua_data(cc).npos_stc_cor
            g_posstc(:,gg) = Xmat_sh(cur_used_inds,:)*mua_data(cc).stcs_cor(:,gg);
        end
        for gg = 1:mua_data(cc).nneg_stc_cor
            g_negstc(:,gg) =  Xmat_sh(cur_used_inds,:)*mua_data(cc).stcs_cor(:,end-nneg+gg);
        end
        cur_Xmat = [g_sta g_posstc.^2 g_negstc.^2 Xblock(cur_used_inds,:)];
        [fitmod] = regGLM_fit(cur_Xmat,Robs,[],[],[],[],1);
        [~, ~, ~, G] = regGLM_eval(fitmod,Robs,cur_Xmat);
        G = G - fitmod.theta;
        [g_dist,g_x] = ksdensity(G); g_dist = g_dist/sum(g_dist);
        mua_data(cc).stc_glm_cor = fitmod;
        
        TB_stim = [time_to_gsac(used_inds(cur_used_inds)) G];
        %         Ytick = linspace(my_prctile(TB_stim(:,2),0.5),my_prctile(TB_stim(:,2),100-0.5),n_Gbins);
        Ytick = linspace(my_prctile(TB_stim(:,2),0.05),my_prctile(TB_stim(:,2),100-0.05),n_Gbins);
        TB = TentBasis2D(Xtick, Ytick);
        used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
            TB_stim(:,2) >= Ytick(1) & TB_stim(:,2) <= Ytick(end));
        [TB_Xmat,TB_counts] = TB.InputNL2D(TB_stim(used_data,:));
        L2_params = create_L2_params([],[1 n_sbins*n_Gbins],[n_sbins n_Gbins],2,3,[Inf Inf],[0.1 1]);
        TB_fitmod = regGLM_fit(TB_Xmat,Robs(used_data),L2_params,10,[],[],1);
        TB_K = reshape(TB_fitmod.K,n_sbins,n_Gbins)';
        bin_areas = TB.GetBinAreas();
        gsac_TB_dist = TB_counts./bin_areas;
        gsac_TB_dist = gsac_TB_dist'/sum(gsac_TB_dist(:));
        gsac_TB_rate = log(1 + exp(TB_K + TB_fitmod.theta));
        
        %INFO CALS
        cur_avg_rate = mean(Robs(used_data));
        marg_gdist = sum(gsac_TB_dist,2);
        marg_sdist = sum(gsac_TB_dist);
        marg_gsacrate = sum(gsac_TB_dist.*gsac_TB_rate)./marg_sdist;
        marg_grate = sum(gsac_TB_dist.*gsac_TB_rate,2)./marg_gdist;
        ov_info_gsac = sum(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate));
        gsacdep_info = nan(n_sac_bins,1);
        for tt = 1:n_sbins
            gsacdep_info(tt) = sum(gsac_TB_dist(:,tt).*gsac_TB_rate(:,tt).*log2(gsac_TB_rate(:,tt)/marg_gsacrate(tt)))/sum(gsac_TB_dist(:,tt));
        end
        mua_data(cc).gsac_TB_dist_cor = gsac_TB_dist;
        mua_data(cc).gsac_TB_rate_cor = gsac_TB_rate;
        mua_data(cc).Gtick_cor = Ytick;
        mua_data(cc).Stick_cor = Xtick;
        mua_data(cc).ov_info_gsac_cor = ov_info_gsac;
        mua_data(cc).gsacdep_info_cor = gsacdep_info;
        mua_data(cc).marg_gsacrate_cor = marg_gsacrate;
        
        TB_stim = [time_to_msac(used_inds(cur_used_inds)) G];
        %         Ytick = linspace(my_prctile(TB_stim(:,2),0.5),my_prctile(TB_stim(:,2),100-0.5),n_Gbins);
        Ytick = linspace(my_prctile(TB_stim(:,2),0.05),my_prctile(TB_stim(:,2),100-0.05),n_Gbins);
        TB = TentBasis2D(Xtick, Ytick);
        used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
            TB_stim(:,2) >= Ytick(1) & TB_stim(:,2) <= Ytick(end));
        [TB_Xmat,TB_counts] = TB.InputNL2D(TB_stim(used_data,:));
        L2_params = create_L2_params([],[1 n_sbins*n_Gbins],[n_sbins n_Gbins],2,3,[Inf Inf],[0.1 1]);
        TB_fitmod = regGLM_fit(TB_Xmat,Robs(used_data),L2_params,10,[],[],1);
        TB_K = reshape(TB_fitmod.K,n_sbins,n_Gbins)';
        bin_areas = TB.GetBinAreas();
        msac_TB_dist = TB_counts./bin_areas;
        msac_TB_dist = msac_TB_dist'/sum(msac_TB_dist(:));
        msac_TB_rate = log(1 + exp(TB_K + TB_fitmod.theta));
        
        %INFO CALS
        cur_avg_rate = mean(Robs(used_data));
        marg_gdist = sum(msac_TB_dist,2);
        marg_sdist = sum(msac_TB_dist);
        marg_msacrate = sum(msac_TB_dist.*msac_TB_rate)./marg_sdist;
        marg_grate = sum(msac_TB_dist.*msac_TB_rate,2)./marg_gdist;
        ov_info_msac = sum(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate));
        msacdep_info = nan(n_sac_bins,1);
        for tt = 1:n_sbins
            msacdep_info(tt) = sum(msac_TB_dist(:,tt).*msac_TB_rate(:,tt).*log2(msac_TB_rate(:,tt)/marg_msacrate(tt)))/sum(msac_TB_dist(:,tt));
        end
        mua_data(cc).msac_TB_dist_cor = msac_TB_dist;
        mua_data(cc).msac_TB_rate_cor = msac_TB_rate;
        mua_data(cc).Gtick_cor = Ytick;
        mua_data(cc).Stick_cor = Xtick;
        mua_data(cc).ov_info_msac_cor = ov_info_msac;
        mua_data(cc).msacdep_info_cor = msacdep_info;
        mua_data(cc).marg_msacrate_cor = marg_msacrate;
        
    else
        mua_data(cc).stc_ver_use = false;
    end
end

%%
cd(anal_dir)
save(anal_name,'sua_data','mua_data');
