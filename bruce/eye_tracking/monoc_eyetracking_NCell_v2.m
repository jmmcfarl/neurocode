clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

Expt_num = 86;
Expt_name = sprintf('G0%d',Expt_num);
data_dir = ['~/Data/bruce/' Expt_name];
cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

anal_name = 'monoc_eyecorr_hbar_Ncell';
mod_data_name = 'monoc_eyecorr_hbar_mods_v4';
recompute_init_mods = 0;

use_measured_pos = 0;
%%
stim_fs = 100; %in Hz
dt = 0.01;
flen = 10;
use_nPix = 24;
full_nPix = 36;
stim_params = NIMcreate_stim_params([flen full_nPix],dt);
Fr = 1;
bar_ori = 0;
min_trial_dur = 0.75;

beg_buffer = 0.2;
end_buffer = 0.05;
trial_dur = 4;

use_right_eye = false;

n_use_blocks = Inf;

spatial_usfac = 2;
use_nPix_us = use_nPix*spatial_usfac;

sac_backlag = round(0.05/dt);
sac_forlag = round(0.3/dt);
sac_bincents = -sac_backlag:sac_forlag;
n_sac_bins = length(sac_bincents);

sp_dx = 0.05/spatial_usfac;
max_shift = round(15*spatial_usfac);
dshift = 1;
shifts = -max_shift:dshift:max_shift;
n_shifts = length(shifts);
zero_frame = find(shifts==0);

%% LOAD OVERALL SU DATA
load ~/Analysis/bruce/summary_analysis/su_data.mat

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
% cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori & expt_sac_dir == bar_ori);
cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori);

if exist('./fin_aclust_data.mat','file')
    load('./fin_aclust_data.mat');
    [n_aclust_expts,n_aclust_probes] = size(autoclust);
else
    disp('No fin_aclust_data found.');
    autoclust = [];
    n_aclust_expts = 0; n_aclust_probes = 0;
end


if strcmp(Expt_name,'G087')
    cur_block_set(cur_block_set == 15) = []; %only 6 trials and causes problems
end

if strcmp(Expt_name,'G087')
    cur_block_set(cur_block_set == 15) = []; %only 6 trials and causes problems
end
if strcmp(Expt_name,'G093')
    cur_block_set(cur_block_set ==  28) = []; %only 6 trials and causes problems
end

n_blocks = length(cur_block_set);

sim_sac_expts = find(~expt_has_ds(cur_block_set));
imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');

if length(cur_block_set) > n_use_blocks
    cur_block_set = cur_block_set(1:n_use_blocks);
end
n_blocks = length(cur_block_set);

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
% all_Xmat = [];
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
for ee = 1:n_blocks;
    if ismember(ee,grayback_gs_expts)
        fprintf('Expt %d Block %d of %d; grayback GS, ori:%d\n',Expt_num,ee,n_blocks,expt_bar_ori(cur_block_set(ee)));
    elseif ismember(ee,imback_gs_expts)
        fprintf('Expt %d Block %d of %d; imback GS, ori:%d\n',Expt_num,ee,n_blocks,expt_bar_ori(cur_block_set(ee)));
    elseif ismember(ee,sim_sac_expts)
        fprintf('Expt %d Block %d of %d; SimSac, ori:%d\n',Expt_num,ee,n_blocks,expt_bar_ori(cur_block_set(ee)));
    else
        fprintf('Expt %d Block %d of %d;  UNMATCHED EXPT TYPE\n',Expt_num,ee,n_blocks);
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
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > min_trial_dur/dt
            use_frames = min(length(cur_stim_times),n_frames);
            cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            %             bar_Xmat = create_time_embedding(cur_stim_mat,stim_params);
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times];
            all_t_axis = [all_t_axis; cur_t_axis];
            %             all_Xmat = [all_Xmat; bar_Xmat];
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

%%
full_nPix_us = spatial_usfac*full_nPix;
if spatial_usfac == 2
    all_stimmat_up = zeros(size(all_stim_mat,1),full_nPix_us);
    for ii = 1:size(all_stim_mat,2)
        all_stimmat_up(:,2*(ii-1)+1) = all_stim_mat(:,ii);
        all_stimmat_up(:,2*(ii-1)+2) = all_stim_mat(:,ii);
    end
elseif spatial_usfac == 1
    all_stimmat_up = all_stim_mat;
else
    error('Unsupported spatial Up-sample factor!');
end
stim_params_us = NMMcreate_stim_params([flen full_nPix_us],dt);
all_Xmat = create_time_embedding(all_stim_mat,stim_params);
all_Xmat_us = create_time_embedding(all_stimmat_up,stim_params_us);

%% select submatrix with central pixels
[Xinds,Tinds] = meshgrid(1:full_nPix,1:flen);
buffer_pix = floor((full_nPix - use_nPix)/2);
cur_use_pix = (1:use_nPix) + buffer_pix;
use_kInds = find(ismember(Xinds(:),cur_use_pix));

[Xinds_up,Tinds] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
if spatial_usfac > 1
    use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1)-1/spatial_usfac & Xinds_up(:) <= cur_use_pix(end));
else
    use_kInds_up = use_kInds;
end
use_kInds_back = find(ismember(Xinds_up(use_kInds_up),cur_use_pix));

%% BIN SPIKES FOR MU AND SU
mahal_thresh = su_data.mah_thresh;
all_binned_spikes = nan(length(all_t_axis),96);
su_used_blocks = false(n_blocks,length(su_probes));
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
    for ee = 1:n_blocks;
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
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer);
if strcmp(Expt_name,'G093')
    un_wi_vals = unique(all_trial_wi);
    use_wi_trials = find(all_trial_wi == un_wi_vals(2));
    used_inds = used_inds(ismember(all_trialvec(used_inds),use_wi_trials));
end
NT = length(used_inds);

%% CREATE EVENT PREDICTORS FOR REAL AND SIM SACCADES (POOLING MICRO AND MACRO SACS)
Xblock = zeros(length(all_stim_times),n_blocks);
for i = 1:n_blocks
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end

%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)

emfile = ['jbe' Expt_name '.em.mat'];
load(emfile);

all_eye_vals = [];
all_eye_speed = [];
all_eye_ts = [];
all_eye_blockvec = [];
eye_smooth = 3;
for ee = 1:n_blocks;
    fprintf('Loading ET data for expt %d, block %d of %d\n',Expt_num,ee,n_blocks);
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

%% BIAS-CORRECTION ON EYE-TRACKING DATA AND REFINE USED DATA TO EXCLUDE DATA WITH EYES OUTSIDE FIXATION WINDOW

%correct for median offsets
corrected_interp_eyevals = interp_eye_vals;
corrected_interp_eyevals = bsxfun(@minus,corrected_interp_eyevals,nanmedian(corrected_interp_eyevals(used_inds,:)));
corrected_eyevals = all_eye_vals;
corrected_eyevals = bsxfun(@minus,corrected_eyevals,nanmedian(interp_eye_vals(used_inds,:)));

%now correct for dependence of orthoganol position on parallel position
par_thresh = 4;
orth_thresh = 1.5;
if bar_ori == 0
    if use_right_eye
        fit_inds = used_inds(abs(corrected_interp_eyevals(used_inds,4)) < orth_thresh & ...
            abs(corrected_interp_eyevals(used_inds,3)) < par_thresh);
        b = robustfit(corrected_interp_eyevals(fit_inds,3),corrected_interp_eyevals(fit_inds,4));
        pred_eyevals = corrected_interp_eyevals(used_inds,3)*b(2) + b(1);
        corrected_interp_eyevals(used_inds,4) = corrected_interp_eyevals(used_inds,4) - pred_eyevals;
        pred_eyevals = corrected_eyevals(:,3)*b(2) + b(1);
        corrected_eyevals(:,4) = corrected_eyevals(:,4) - pred_eyevals;
    end
    
    fit_inds = used_inds(abs(corrected_interp_eyevals(used_inds,2)) < orth_thresh & ...
        abs(corrected_interp_eyevals(used_inds,1)) < par_thresh);
    b = robustfit(corrected_interp_eyevals(fit_inds,1),corrected_interp_eyevals(fit_inds,2));
    pred_eyevals = corrected_interp_eyevals(used_inds,1)*b(2) + b(1);
    corrected_interp_eyevals(used_inds,2) = corrected_interp_eyevals(used_inds,2) - pred_eyevals;
    pred_eyevals = corrected_eyevals(:,1)*b(2) + b(1);
    corrected_eyevals(:,2) = corrected_eyevals(:,2) - pred_eyevals;
else
    if use_right_eye
        fit_inds = used_inds(abs(corrected_interp_eyevals(used_inds,3)) < orth_thresh & ...
            abs(corrected_interp_eyevals(used_inds,4)) < par_thresh);
        b = robustfit(corrected_interp_eyevals(fit_inds,4),corrected_interp_eyevals(fit_inds,3));
        pred_eyevals = corrected_interp_eyevals(used_inds,4)*b(2) + b(1);
        corrected_interp_eyevals(used_inds,3) = corrected_interp_eyevals(used_inds,3) - pred_eyevals;
        pred_eyevals = corrected_eyevals(:,4)*b(2) + b(1);
        corrected_eyevals(:,3) = corrected_eyevals(:,3) - pred_eyevals;
    end
    
    fit_inds = used_inds(abs(corrected_interp_eyevals(used_inds,1)) < orth_thresh & ...
        abs(corrected_interp_eyevals(used_inds,2)) < par_thresh);
    b = robustfit(corrected_interp_eyevals(fit_inds,2),corrected_interp_eyevals(fit_inds,1));
    pred_eyevals = corrected_interp_eyevals(used_inds,2)*b(2) + b(1);
    corrected_interp_eyevals(used_inds,1) = corrected_interp_eyevals(used_inds,1) - pred_eyevals;
    pred_eyevals = corrected_eyevals(:,2)*b(2) + b(1);
    corrected_eyevals(:,1) = corrected_eyevals(:,1) - pred_eyevals;
end
if bar_ori == 0
    outside_window_left = find(abs(corrected_interp_eyevals(used_inds,1)) > par_thresh | ...
        abs(corrected_interp_eyevals(used_inds,2)) > orth_thresh);
    outside_window_right = find(abs(corrected_interp_eyevals(used_inds,3)) > par_thresh | ...
        abs(corrected_interp_eyevals(used_inds,4)) > orth_thresh);
else
    outside_window_left = find(abs(corrected_interp_eyevals(used_inds,2)) > par_thresh | ...
        abs(corrected_interp_eyevals(used_inds,1)) > orth_thresh);
    outside_window_right = find(abs(corrected_interp_eyevals(used_inds,4)) > par_thresh | ...
        abs(corrected_interp_eyevals(used_inds,3)) > orth_thresh);
end
if use_right_eye
    outside_window = union(outside_window_left,outside_window_right);
else
    outside_window = outside_window_left;
end

%% DETECT TIMES WHEN TRIAL FAILURES OCCUR AND FORCE-END TRIALS
par_thresh = 4;
orth_thresh = 1.2;
back_samps = 1;

n_trials = length(all_trial_start_times);
out_of_trial = [];
for ii = 1:n_trials
    cur_trial_inds = used_inds(all_trialvec(used_inds) == ii);
    if ~isempty(cur_trial_inds)
        if bar_ori == 0
            cur_out = find(abs(corrected_interp_eyevals(cur_trial_inds,2)) >= orth_thresh | ...
                abs(corrected_interp_eyevals(cur_trial_inds,1)) >= par_thresh,1,'first');
        elseif bar_ori == 90
            cur_out = find(abs(corrected_interp_eyevals(cur_trial_inds,1)) >= orth_thresh | ...
                abs(corrected_interp_eyevals(cur_trial_inds,2)) >= par_thresh,1,'first');
        end
        if ~isempty(cur_out)
            if cur_out > back_samps
                cur_out = cur_out - back_samps;
            end
            cur_out = cur_out:length(cur_trial_inds);
        end
        out_of_trial = cat(1,out_of_trial,cur_trial_inds(cur_out(:)));
    end
end
fract_out = length(out_of_trial)/length(used_inds);
fprintf('Eliminating %.4f of data out of window\n',fract_out);
used_inds(ismember(used_inds,out_of_trial)) = [];
NT = length(used_inds);

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
    %     sac_pre_pos(ii,:) = median(all_eye_vals(pre_inds,:),1);
    sac_pre_pos(ii,:) = median(corrected_eyevals(pre_inds,:),1);
    post_inds = sac_stop_inds(ii):(sac_stop_inds(ii) + sac_dbuff);
    post_inds(post_inds > length(all_eye_ts)) = length(all_eye_ts);
    %     sac_post_pos(ii,:) = median(all_eye_vals(post_inds,:),1);
    sac_post_pos(ii,:) = median(corrected_eyevals(post_inds,:),1);
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

sac_start_times = [saccades(:).start_time];
interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_start_inds(isnan(interp_sac_start_inds)) = 1;
sac_error = abs(sac_start_times - all_t_axis(interp_sac_start_inds)');
bad_sacs = find(isnan(interp_sac_start_inds) | sac_error > dt);
sac_start_times(bad_sacs) = [];
saccades(bad_sacs) = [];
interp_sac_start_inds(bad_sacs) = [];

%% CREATE SACCADE PREDICTOR MATS
saccade_start_inds = find(ismember(used_inds,interp_sac_start_inds));
used_saccade_set = find(ismember(interp_sac_start_inds,used_inds));

sac_amps = [saccades(:).amplitude];
is_micro = sac_amps(used_saccade_set) < 1;
big_sacs = find(~is_micro);
micro_sacs = find(is_micro);

saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));
Xsac = zeros(NT,n_sac_bins);
Xmsac = zeros(NT,n_sac_bins);
for ii = 1:n_sac_bins
    cur_sac_target = saccade_start_inds(big_sacs) + sac_bincents(ii);
    uu = find(cur_sac_target > 1 & cur_sac_target < NT);
    cur_sac_target = cur_sac_target(uu);
    cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_trial_inds(big_sacs(uu))) = [];
    Xsac(cur_sac_target,ii) = 1;
    
    cur_sac_target = saccade_start_inds(micro_sacs) + sac_bincents(ii);
    uu = find(cur_sac_target > 1 & cur_sac_target < NT);
    cur_sac_target = cur_sac_target(uu);
    cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_trial_inds(micro_sacs(uu))) = [];
    Xmsac(cur_sac_target,ii) = 1;
end

%% DEFINE POINTS TO ALLOW RAPID CHANGES (TRIAL STARTS AFTER SACCADES)
trial_start_inds = [1; find(diff(all_trialvec(used_inds)) ~= 0) + 1];

%push the effects of saccades forward in time
sac_shift = round(0.05/dt);

%define times when to resort to independent prior
use_prior = false(size(used_inds));
use_prior(trial_start_inds) = true;
for i = 1:length(saccade_start_inds)
    next_trial = trial_start_inds(find(trial_start_inds > saccade_start_inds(i),1,'first'));
    %mark a point (forward in time) as a saccade if it occurs within the
    %same trial as the saccade
    if next_trial > saccade_start_inds(i) + sac_shift
        use_prior(saccade_start_inds(i) + sac_shift) = 1;
    end
end

jump_inds = [find(use_prior == 1); length(used_inds)+1];
n_fixs = length(jump_inds)-1;
fix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = jump_inds(ii):(jump_inds(ii+1)-1);
    fix_ids(cur_inds) = ii;
end

%%
all_Xmat_us = create_time_embedding(all_stimmat_up,stim_params_us);
all_Xmat = create_time_embedding(all_stim_mat,stim_params);

%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITSi init_stim_params = NMMcreate_stim_params([flen use_nPix],dt);
init_stim_params(2) = NMMcreate_stim_params([n_blocks 1],dt);

fin_stim_params(1) = NMMcreate_stim_params([flen use_nPix_us],dt);
fin_stim_params(2) = NMMcreate_stim_params([n_blocks 1],dt);
fin_stim_params(3) = NMMcreate_stim_params(n_sac_bins,dt);
fin_stim_params(4) = NMMcreate_stim_params(n_sac_bins,dt);

null_stim_params = fin_stim_params(2:end);

block_L2 = 1;
silent = 1;
sac_d2t = 100;

base_lambda_d2XT = 8;
base_lambda_L1 = 5;

init_optim_p.optTol = 1e-5; init_optim_p.progTol = 1e-9;

sac_reg_params = NMMcreate_reg_params('lambda_d2T',sac_d2t);
null_reg_params = NMMcreate_reg_params('lambda_d2T',[0; sac_d2t; sac_d2t],'lambda_L2',[block_L2; 0; 0]);

n_squared_filts = 2;
mod_signs = ones(1,n_squared_filts+2);
NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts) {'lin'}];
init_d2XT = [ones(n_squared_filts+1,1); 0;];
init_L2 = [zeros(n_squared_filts+1,1); block_L2];
init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_d2XT,'lambda_L2',init_L2);
init_Xtargs = [ones(n_squared_filts+1,1); 2];

init_filts = cell(length(mod_signs),1);
cd(anal_dir);

if ~exist(['./' mod_data_name '.mat'],'file') || recompute_init_mods
    tot_nUnits = length(su_probes) + 96;
    all_mod_SU = zeros(tot_nUnits,1);
    for ss = 1:96;
        fprintf('Computing base LLs for MU %d of %d\n',ss,96);
        su_probe_ind = find(su_probes == ss);
        if ~isempty(su_probe_ind)
            cur_used_blocks = find(su_used_blocks(:,su_probe_ind) == 0); %blocks when NO SU
            cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
        else
            cur_used_blocks = 1:n_blocks; %blocks when NO SU
            cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
        end
        cur_used_inds = find(ismember(used_inds,cur_poss_inds));
        cur_NT = length(cur_used_inds);
        if ~isempty(cur_used_inds)
            Robs = all_binned_spikes(used_inds(cur_used_inds),ss);
            
            cur_X{1} = all_Xmat(used_inds(cur_used_inds),use_kInds);
            cur_X{2} = Xblock(used_inds(cur_used_inds),:);
            cur_X{3} = Xsac(cur_used_inds,:);
            cur_X{4} = Xmsac(cur_used_inds,:);
            
            null_mod = NMMinitialize_model(null_stim_params,[1 1 1],repmat({'lin'},1,3),null_reg_params,[1 2 3]);
            null_mod = NMMfit_filters(null_mod,Robs,cur_X(2:end));
            all_nullmod(ss) = null_mod;
            
            gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
            gqm1 = NMMfit_filters(gqm1,Robs,cur_X,[],[],silent,init_optim_p);
            
            cur_X{1} = all_Xmat_us(used_inds(cur_used_inds),use_kInds_up);
            %spatial up-sampling of filter estimates
            base_filts = reshape([gqm1.mods(find(init_Xtargs == 1)).filtK],[flen use_nPix n_squared_filts+1]);
            if spatial_usfac == 2
                base_filts_up = zeros(flen,use_nPix_us,n_squared_filts+1);
                for ii = 1:use_nPix
                    base_filts_up(:,2*(ii-1)+1,:) = 0.5*base_filts(:,ii,:);
                    base_filts_up(:,2*(ii-1)+2,:) = 0.5*base_filts(:,ii,:);
                end
            elseif spatial_usfac == 1
                base_filts_up = base_filts;
            else
                error('unsupported')
            end
            base_filts_up = reshape(base_filts_up,use_nPix_us*flen,n_squared_filts+1);
            
            init_filts{end} = gqm1.mods(find(init_Xtargs==2)).filtK;
            for ii = 1:n_squared_filts+1
                init_filts{ii} = base_filts_up(:,ii);
            end
            gqm2 = NMMinitialize_model(fin_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs,init_filts);
            gqm2.spk_NL_params(1) = gqm1.spk_NL_params(1);
            [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm2,Robs,cur_X);
            gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
            gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
            gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,3,zeros(n_sac_bins,1),sac_reg_params); %sac term
            gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,4,zeros(n_sac_bins,1),sac_reg_params); %nsac term
            gqm2 = NMMfit_filters(gqm2,Robs,cur_X,[],[],silent);
            
            all_mod_fits(ss) = gqm2;
            LL = NMMmodel_eval(gqm2,Robs,cur_X);
            all_mod_LLimp(ss) = (LL-null_mod.LL_seq(end))/log(2);
            
            all_mod_fits_withspkNL(ss) = NMMfit_logexp_spkNL(gqm2,Robs,cur_X);
        end
    end
    
    for ss = 1:length(su_probes);
        fprintf('Computing base LLs for SU %d of %d\n',ss,length(su_probes));
        all_mod_SU(ss+96) = su_probes(ss);
        cur_used_blocks = find(su_used_blocks(:,ss)); %blocks when SU
        cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
        cur_used_inds = find(ismember(used_inds,cur_poss_inds));
        cur_NT = length(cur_used_inds);
        if ~isempty(cur_used_inds)
            Robs = all_binned_spikes(used_inds(cur_used_inds),su_probes(ss));
            
            cur_X{1} = all_Xmat(used_inds(cur_used_inds),use_kInds);
            cur_X{2} = Xblock(used_inds(cur_used_inds),:);
            cur_X{3} = Xsac(cur_used_inds,:);
            cur_X{4} = Xmsac(cur_used_inds,:);
            
            null_mod = NMMinitialize_model(null_stim_params,[1 1 1],repmat({'lin'},1,3),null_reg_params,[1 2 3]);
            null_mod = NMMfit_filters(null_mod,Robs,cur_X(2:end));
            all_nullmod(ss+96) = null_mod;
            
            gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
            gqm1 = NMMfit_filters(gqm1,Robs,cur_X,[],[],silent,init_optim_p);
            
            cur_X{1} = all_Xmat_us(used_inds(cur_used_inds),use_kInds_up);
            %spatial up-sampling of filter estimates
            base_filts = reshape([gqm1.mods(find(init_Xtargs == 1)).filtK],[flen use_nPix n_squared_filts+1]);
            if spatial_usfac == 2
                base_filts_up = zeros(flen,use_nPix_us,n_squared_filts+1);
                for ii = 1:use_nPix
                    base_filts_up(:,2*(ii-1)+1,:) = 0.5*base_filts(:,ii,:);
                    base_filts_up(:,2*(ii-1)+2,:) = 0.5*base_filts(:,ii,:);
                end
            elseif spatial_usfac == 1
                base_filts_up = base_filts;
            else
                error('unsupported')
            end
            base_filts_up = reshape(base_filts_up,use_nPix_us*flen,n_squared_filts+1);
            
            init_filts{end} = gqm1.mods(find(init_Xtargs==2)).filtK;
            for ii = 1:n_squared_filts+1
                init_filts{ii} = base_filts_up(:,ii);
            end
            gqm2 = NMMinitialize_model(fin_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs,init_filts);
            gqm2.spk_NL_params(1) = gqm1.spk_NL_params(1);
            [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm2,Robs,cur_X);
            gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
            gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
            gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,3,zeros(n_sac_bins,1),sac_reg_params); %sac term
            gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,4,zeros(n_sac_bins,1),sac_reg_params); %nsac term
            gqm2 = NMMfit_filters(gqm2,Robs,cur_X,[],[],silent);
            
            all_mod_fits(ss+96) = gqm2;
            LL = NMMmodel_eval(gqm2,Robs,cur_X);
            all_mod_LLimp(ss+96) = (LL-null_mod.LL_seq(end))/log(2);
            
            all_mod_fits_withspkNL(ss+96) = NMMfit_logexp_spkNL(gqm2,Robs,cur_X);
        end
    end
    save(mod_data_name,'all_mod*','all_nullmod','su_probes');
else
    fprintf('Loading pre-computed initial models\n');
    load(mod_data_name);
end

%% INITIALIZE TRANSITION PRIORS FOR HMM
%overall prior on shifts
eps_prior_sigma = 0.15; %0.125 start
leps_prior = -(shifts*sp_dx).^2/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize

%% SELECT USABLE UNITS AND make Robs_mat
LL_imp_thresh = 5e-4;
usable_units = find(all_mod_LLimp >= LL_imp_thresh);
max_units = length(usable_units);

n_iter = 3;

Ncells = [max_units -1 90 80 70 60 50 40 30 20 10];
for cc = 1:length(Ncells)
    if Ncells(cc) > 0
        rperm = randperm(length(usable_units));
        tr_set = usable_units(rperm(1:Ncells(cc)));
    else
        tr_set = usable_units(all_mod_SU(usable_units) == 0);
    end
    Ncell_tr_set{cc} = tr_set;
    n_used_sus = sum(all_mod_SU(usable_units) ~= 0);
    n_used_mus = sum(all_mod_SU(usable_units) == 0);
    fprintf('Using %d SUs and %d MUs for analysis\n',n_used_sus,n_used_mus);
    n_tr_chs = length(tr_set);
    
    % make Robs_mat
    poss_blocks = 1:(n_blocks-1);
    Robs_mat = nan(length(used_inds),n_tr_chs);
    for ss = 1:n_tr_chs
        if all_mod_SU(tr_set(ss)) == 0
            su_probe_ind = find(su_probes == tr_set(ss));
            if ~isempty(su_probe_ind)
                cur_used_blocks = find(su_used_blocks(:,su_probe_ind) == 0); %blocks when NO SU
            else
                cur_used_blocks = 1:n_blocks; %blocks when NO SU
            end
            cur_use = find(ismember(all_blockvec(used_inds),cur_used_blocks));
            Robs_mat(cur_use,ss) = all_binned_spikes(used_inds(cur_use),tr_set(ss));
        else
            su_ind = find(su_probes == all_mod_SU(tr_set(ss)));
            cur_used_blocks = find(su_used_blocks(:,su_ind)); %blocks when SU
            cur_use = find(ismember(all_blockvec(used_inds),cur_used_blocks));
            Robs_mat(cur_use,ss) = all_binned_spikes(used_inds(cur_use),all_mod_SU(tr_set(ss)));
        end
    end
    
    %%
    it_ref_mods{cc}{1} = all_mod_fits;
    it_ref_mods_spkNL{cc}{1} = all_mod_fits_withspkNL;
    it_ref_LL_imp{cc}(1,:) = all_mod_LLimp;
    for nn = 2:(n_iter+1)
        %% PREPROCESS MODEL COMPONENTS
        use_nPix_up = use_nPix*spatial_usfac;
        klen = flen*use_nPix_up;
        filt_bank = zeros(n_tr_chs,klen,n_squared_filts+1);
        lin_kerns = nan(n_tr_chs,n_blocks);
        sac_kerns = nan(n_tr_chs,n_sac_bins);
        msac_kerns = nan(n_tr_chs,n_sac_bins);
        mod_spkNL_params = nan(n_tr_chs,3);
        for ss = 1:n_tr_chs
            cur_Xtargs = [it_ref_mods{cc}{nn-1}(tr_set(ss)).mods(:).Xtarget];
            cur_k = [it_ref_mods{cc}{nn-1}(tr_set(ss)).mods(cur_Xtargs == 1).filtK];
            n_used_filts = size(cur_k,2);
            filt_bank(ss,:,1:n_used_filts) = cur_k;
            mod_spkNL_params(ss,:) = it_ref_mods_spkNL{cc}{nn-1}(tr_set(ss)).spk_NL_params;
            lin_kerns(ss,:) = it_ref_mods{cc}{nn-1}(tr_set(ss)).mods(cur_Xtargs == 2).filtK;
            sac_kerns(ss,:) = it_ref_mods{cc}{nn-1}(tr_set(ss)).mods(cur_Xtargs == 3).filtK;
            msac_kerns(ss,:) = it_ref_mods{cc}{nn-1}(tr_set(ss)).mods(cur_Xtargs == 4).filtK;
        end
        filt_bank = permute(filt_bank,[2 1 3]);
        filt_bank = reshape(filt_bank,[flen use_nPix_up n_tr_chs n_squared_filts+1]);
        
        %indicator predictions
        block_out = Xblock(used_inds,:)*lin_kerns';
        sac_out = Xsac*sac_kerns';
        msac_out = Xmsac*msac_kerns';
        
        %% ESTIMATE LL for each shift in each stimulus frame
        cur_Xmat = all_Xmat_us(used_inds,use_kInds_up);
        
        %precompute LL at all shifts for all units
        frame_LLs = nan(NT,n_shifts);
        for xx = 1:length(shifts)
            fprintf('Shift %d of %d\n',xx,n_shifts);
            cur_filt_shift = shift_matrix_Nd(filt_bank,shifts(xx),2);
            cur_filt_shift = reshape(cur_filt_shift,[klen n_tr_chs n_squared_filts+1]);
            
            %outputs of stimulus models at current X-matrix shift
            gfuns = ones(NT,n_tr_chs);
            gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,1)');
            gfuns = gfuns + cur_Xmat*squeeze(cur_filt_shift(:,:,1));
            for ff = 2:(n_squared_filts+1)
                gfuns = gfuns + (cur_Xmat*squeeze(cur_filt_shift(:,:,ff))).^2;
            end
            
            %add contributions from extra lin kernels
            gfuns = gfuns + block_out + sac_out + msac_out;
            
            %incorporate beta
            gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,2)');
            
            too_large = gfuns > 50;
            pred_rate = log(1+exp(gfuns));
            pred_rate(too_large) = gfuns(too_large);
            
            %incorporate alpha
            pred_rate = bsxfun(@times,pred_rate,mod_spkNL_params(:,3)');
            
            pred_rate(pred_rate < 1e-20) = 1e-20;
            
            frame_LLs(:,xx) = squeeze(nansum(Robs_mat.*log(pred_rate) - pred_rate,2));
        end
        
        %% INFER MICRO-SAC SEQUENCE
        jump_inds = [find(use_prior == 1); length(used_inds)+1];
        n_fixs = length(jump_inds)-1;
        fix_LLs = nan(n_fixs,n_shifts);
        fix_ids = nan(NT,1);
        for ii = 1:n_fixs
            cur_inds = jump_inds(ii):(jump_inds(ii+1)-1);
            fix_LLs(ii,:) = sum(frame_LLs(cur_inds,:));
            fix_ids(cur_inds) = ii;
        end
        
        lPost = bsxfun(@plus,fix_LLs,leps_prior);
        lPost = bsxfun(@minus,lPost,logsumexp(lPost,2));
        Post = exp(lPost);
        fix_post_mean{cc}(nn,:) = sum(bsxfun(@times,Post,shifts),2);
        % cur_post_std = sqrt(sum(bsxfun(@times,Post,shifts.^2),2) - cur_post_mean.^2);
        cur_diff = bsxfun(@minus,fix_post_mean{cc}(nn,:)',shifts).^2;
        fix_post_std{cc}(nn,:) = sqrt(sum(cur_diff.*Post,2));
        
        all_fix_post_mean = fix_post_mean{cc}(nn,fix_ids);
        all_fix_post_std = fix_post_std{cc}(nn,fix_ids);
        
        %back-project saccade-times
        all_fix_post_mean_cor = all_fix_post_mean;
        %align jumps to saccade times
        for i = length(saccade_start_inds):-1:1
            cur_inds = saccade_start_inds(i):(saccade_start_inds(i)+sac_shift);
            all_fix_post_mean_cor(cur_inds) = all_fix_post_mean(cur_inds(end));
            cur_back_inds = (saccade_start_inds(i)-(flen-sac_shift)):(saccade_start_inds(i)-1);
            all_fix_post_mean_cor(cur_back_inds) = all_fix_post_mean(saccade_start_inds(i));
        end
        %%
        all_shift_stimmat_up = all_stimmat_up;
        for i=1:NT
            all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-round(all_fix_post_mean_cor(i)),2);
        end
        all_Xmat_up_fixcor = create_time_embedding(all_shift_stimmat_up,stim_params_us);
        all_Xmat_up_fixcor = all_Xmat_up_fixcor(used_inds,:);
        
        %% REFIT ALL CELLS
        silent = 1;
        % cur_fit_set = setdiff(tr_set,poss_xv_set);
        cur_fit_set = tr_set;
        for ss = 1:length(cur_fit_set)
            tr_cell = cur_fit_set(ss);
            fprintf('Refitting model for tr cell %d of %d\n',ss,length(cur_fit_set));
            cur_unit_ind = find(tr_set == tr_cell);
            cur_uset = find(~isnan(Robs_mat(:,cur_unit_ind)));
            
            cur_X{1} = all_Xmat_up_fixcor(cur_uset,use_kInds_up);
            cur_X{2} = Xblock(used_inds(cur_uset),:);
            cur_X{3} = Xsac(cur_uset,:);
            cur_X{4} = Xmsac(cur_uset,:);
            
            it_ref_mods{cc}{nn}(tr_cell) = it_ref_mods{cc}{nn-1}(tr_cell);
            it_ref_mods{cc}{nn}(tr_cell) = NMMfit_filters(it_ref_mods{cc}{nn}(tr_cell),Robs_mat(cur_uset,cur_unit_ind),...
                cur_X,[],[],silent); %fit stimulus filters
            
            it_ref_LL_imp{cc}(nn,tr_cell) = (it_ref_mods{cc}{nn}(tr_cell).LL_seq(end) - all_nullmod(tr_cell).LL_seq(end))/log(2);
            
            fprintf('Original: %.4f  New: %.4f\n',all_mod_LLimp(tr_cell),it_ref_LL_imp{cc}(nn,tr_cell));
            
            %refit spk NL
            it_ref_mods_spkNL{cc}{nn}(tr_cell) = NMMfit_logexp_spkNL(it_ref_mods{cc}{nn}(tr_cell),Robs_mat(cur_uset,cur_unit_ind),cur_X);
            
        end        
    end
end

%%
cd(anal_dir)
save(anal_name,'fix_*','Ncell*','all_mod*','it_ref*');
%%
min_std = 0.05;
ww_std = fix_post_std{1}(end,:);
% min_std = max(ww_std);

ww_std(ww_std < min_std) = min_std;
ww = 1./ww_std;
% ww = ones(n_fixs,1);
ww = ww/sum(ww);

for nn = 1:length(Ncells);
    sigx = fix_post_mean{1}(end,:);
    sigy = fix_post_mean{nn}(end,:);
    mx = sum(ww.*sigx);
    my = sum(ww.*sigy);
    cov_xy = sum(ww.*(sigx - mx).*(sigy - my));
    cov_xx = sum(ww.*(sigx - mx).*(sigx - mx));
    cov_yy = sum(ww.*(sigy - my).*(sigy - my));
    
    it_weight_r(nn) = cov_xy./sqrt(cov_xx*cov_yy);
    it_r(nn) = corr(sigx',sigy');
    it_r_sp(nn) = corr(sigx',sigy','type','spearman');
end

use_Ncells = Ncells;
use_Ncells(use_Ncells < 0) = sum(all_mod_SU(Ncell_tr_set{1}) == 0);

uset = 1:length(use_Ncells);
uset(2) = [];
figure
plot(use_Ncells(uset),it_weight_r(uset),'ko-','markersize',6,'linewidth',1.5);
hold on
plot(use_Ncells(uset),it_r(uset),'ro-','markersize',6,'linewidth',1.5);
plot(use_Ncells(uset),it_r_sp(uset),'bo-','markersize',6,'linewidth',1.5)
legend('Weighted Pearson','Unweighted Pearson','Spearman','Location','Southeast');
xlabel('Number of units','fontsize',12);
ylabel('Correlation','fontsize',12);
set(gca,'fontsize',8,'fontname','arial');
box off
fillPage(gcf,'papersize',[4 4]);

%%
for nn = 1:length(Ncells)
    LL_imp_set = it_ref_LL_imp{nn}(end,Ncell_tr_set{nn})./all_mod_LLimp(Ncell_tr_set{nn});
   avg_LL_imp(nn) = mean(LL_imp_set);
   std_LL_imp(nn) = std(LL_imp_set);
end
%%
jitter_std = 0.01;
NF = size(fix_post_mean{1},2);
msize = 3;
cmap = jet(9);
figure
subplot(3,1,1);hold on
line([-0.8 0.8],[-0.8 0.8],'color','k');
line([0 0],[-0.8 0.8],'color','k','linestyle','--');
line([-0.8 0.8],[0 0],'color','k','linestyle','--');
plot(fix_post_mean{1}(end,:)*sp_dx + randn(1,NF)*jitter_std,fix_post_mean{4}(end,:)*sp_dx + randn(1,NF)*jitter_std,'.','color',[1 0 0],'markersize',msize);
xlim([-0.8 0.8]); ylim([-0.8 0.8]);
box off
xlabel('Inferred position, all cells (deg)','fontsize',12);
ylabel('Inferred position, subset (deg)','fontsize',12);
title('80 units','fontsize',12);
set(gca,'fontname','arial','fontsize',10);

subplot(3,1,2);hold on
line([-0.8 0.8],[-0.8 0.8],'color','k');
line([0 0],[-0.8 0.8],'color','k','linestyle','--');
line([-0.8 0.8],[0 0],'color','k','linestyle','--');
plot(fix_post_mean{1}(end,:)*sp_dx + randn(1,NF)*jitter_std,fix_post_mean{7}(end,:)*sp_dx + randn(1,NF)*jitter_std,'.','color',[0.2 0.8 0.2],'markersize',msize);
xlim([-0.8 0.8]); ylim([-0.8 0.8]);
box off
xlabel('Inferred position, all cells (deg)','fontsize',12);
ylabel('Inferred position, subset (deg)','fontsize',12);
title('50 units','fontsize',12);
set(gca,'fontname','arial','fontsize',10);


subplot(3,1,3);hold on
line([-0.8 0.8],[-0.8 0.8],'color','k');
line([0 0],[-0.8 0.8],'color','k','linestyle','--');
line([-0.8 0.8],[0 0],'color','k','linestyle','--');
plot(fix_post_mean{1}(end,:)*sp_dx + randn(1,NF)*jitter_std,fix_post_mean{10}(end,:)*sp_dx + randn(1,NF)*jitter_std,'.','color',[0.8 0.8 0.2],'markersize',msize);
xlim([-0.8 0.8]); ylim([-0.8 0.8]);
box off
xlabel('Inferred position, all cells (deg)','fontsize',12);
ylabel('Inferred position, subset (deg)','fontsize',12);
title('20 units','fontsize',12);
set(gca,'fontname','arial','fontsize',10);


fillPage(gcf,'papersize',[4 12]);
fname = 'Ncell_examp_scatters';


%%
fix_corr_mean1 = fix_post_mean{1}(end,fix_ids);
fix_corr_std1 = fix_post_std{1}(end,fix_ids); %will want to use corrected versions of the uncertainty too!
fix_corr_mean2 = fix_post_mean{4}(end,fix_ids);
fix_corr_std2 = fix_post_std{4}(end,fix_ids); %will want to use corrected versions of the uncertainty too!
fix_corr_mean3 = fix_post_mean{7}(end,fix_ids);
fix_corr_std3 = fix_post_std{7}(end,fix_ids); %will want to use corrected versions of the uncertainty too!

close all
cmap = jet(9);
n_trials = length(unique(all_trialvec));
for tt = 96:n_trials
% for tt = [96 137 154 179 376 409]
    uu = find(all_trialvec(used_inds) == tt);
    if ~isempty(uu)
        bt = all_t_axis(used_inds(uu(1)));
        et = all_t_axis(used_inds(uu(end)));
        dur = et-bt;
        if dur > 3.5
            hold on
            h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fix_corr_mean1(uu)*sp_dx,fix_corr_std1(uu)*sp_dx,{'color','k'});
%             h2=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fix_corr_mean2(uu)*sp_dx,fix_corr_std2(uu)*sp_dx,{'color','m'});
%             h3=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fix_corr_mean3(uu)*sp_dx,fix_corr_std3(uu)*sp_dx,{'color','g'});
            if bar_ori == 0
                h4=plot(all_t_axis(used_inds(uu))-bt,corrected_interp_eyevals(used_inds(uu),2),'r','linewidth',2);
            else
                h4=plot(all_t_axis(used_inds(uu))-bt,corrected_interp_eyevals(used_inds(uu),1),'r','linewidth',2);
            end            
            for ii = 4:3:length(Ncells)
                plot(all_t_axis(used_inds(uu))-bt,fix_post_mean{ii}(end,fix_ids(uu))*sp_dx,'color',cmap(ii-2,:),'linewidth',1.5);
            end

            
%             legend([h1.mainLine h2.mainLine h3.mainLine h4],{'All units','MU only','Left eye'})
            xlim([0 dur]);
            ylim([-0.45 0.45]);
            xlabel('Time (s)','fontsize',10);
            ylabel('Orthoganol position (deg)','fontsize',10);
            title(sprintf('Trial %d',tt));
            set(gca,'fontsize',8,'fontname','arial');
            fillPage(gcf,'papersize',[8 5]);
            pause
            clf
        end
    end
end

%%
close all
nn = 4;
use_su = 104;
mod1 = it_ref_mods{1}{nn}(use_su);
xtargs = [mod1.mods(:).Xtarget];
mod1_kmat = [mod1.mods(xtargs == 1).filtK];

mod2 = it_ref_mods{4}{nn}(use_su);
xtargs = [mod2.mods(:).Xtarget];
mod2_kmat = [mod2.mods(xtargs == 1).filtK];

mod3 = it_ref_mods{7}{nn}(use_su);
xtargs = [mod3.mods(:).Xtarget];
mod3_kmat = [mod3.mods(xtargs == 1).filtK];

mod4 = it_ref_mods{10}{nn}(use_su);
xtargs = [mod4.mods(:).Xtarget];
mod4_kmat = [mod4.mods(xtargs == 1).filtK];

orig_mod = all_mod_fits(use_su);
xtargs = [orig_mod.mods(:).Xtarget];
orig_kmat = [orig_mod.mods(xtargs == 1).filtK];

%compute all filter output distros

n_stimfilts = size(mod1_kmat,2);
for ii = 1:n_stimfilts
    subplot(n_stimfilts,4,(ii-1)*4+1)
    imagesc(reshape(orig_kmat(:,ii),flen,use_nPix_us));
    ca = max(abs(mod1_kmat(:,ii))); caxis([-ca ca]);
    set(gca,'xtick',[],'ytick',[]);
    if ii == 1
        title('Original','fontsize',10);
    end
    subplot(n_stimfilts,4,(ii-1)*4+2)
    imagesc(reshape(mod4_kmat(:,ii),flen,use_nPix_us));
    ca = max(abs(mod1_kmat(:,ii))); caxis([-ca ca]);
    set(gca,'xtick',[],'ytick',[]);
    if ii == 1
        title('20 units','fontsize',10);
    end
    
    subplot(n_stimfilts,4,(ii-1)*4+3)
    imagesc(reshape(mod3_kmat(:,ii),flen,use_nPix_us));
    ca = max(abs(mod1_kmat(:,ii))); caxis([-ca ca]);
    set(gca,'xtick',[],'ytick',[]);
    if ii == 1
        title('50 units','fontsize',10);
    end
% 
%     subplot(n_stimfilts,5,(ii-1)*5+4)
%     imagesc(reshape(mod2_kmat(:,ii),flen,use_nPix_us));
%     ca = max(abs(mod1_kmat(:,ii))); caxis([-ca ca]);
%     set(gca,'xtick',[],'ytick',[]);
%     if ii == 1
%         title('80 units','fontsize',10);
%     end
    
    subplot(n_stimfilts,4,(ii-1)*4+4)
    imagesc(reshape(mod1_kmat(:,ii),flen,use_nPix_us));
    ca = max(abs(mod1_kmat(:,ii))); caxis([-ca ca]);
    set(gca,'xtick',[],'ytick',[]);
    if ii == 1
        title('All units','fontsize',10);
    end    
    
end
colormap(gray);
fillPage(gcf,'papersize',[9 6]);
%% FOR REFITTING A UNIT USING A PARTICULAR INFERRED FIXATION SEQUENCE
cell_num = 104;
cc = 10;
nn = 4;
su_ind = find(su_probes == all_mod_SU(cell_num));
cur_used_blocks = find(su_used_blocks(:,su_ind)); %blocks when SU
cur_use = find(ismember(all_blockvec(used_inds),cur_used_blocks));
Robs_mat(cur_use,cell_num) = all_binned_spikes(used_inds(cur_use),all_mod_SU(cell_num));

all_fix_post_mean = fix_post_mean{cc}(nn,fix_ids);
%back-project saccade-times
all_fix_post_mean_cor = all_fix_post_mean;
%align jumps to saccade times
for i = length(saccade_start_inds):-1:1
    cur_inds = saccade_start_inds(i):(saccade_start_inds(i)+sac_shift);
    all_fix_post_mean_cor(cur_inds) = all_fix_post_mean(cur_inds(end));
    cur_back_inds = (saccade_start_inds(i)-(flen-sac_shift)):(saccade_start_inds(i)-1);
    all_fix_post_mean_cor(cur_back_inds) = all_fix_post_mean(saccade_start_inds(i));
end
all_shift_stimmat_up = all_stimmat_up;
for i=1:NT
    all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-round(all_fix_post_mean_cor(i)),2);
end
all_Xmat_up_fixcor = create_time_embedding(all_shift_stimmat_up,stim_params_us);
all_Xmat_up_fixcor = all_Xmat_up_fixcor(used_inds,:);

silent = 0;
fprintf('Refitting model for tr cell %d\n',cell_num);
cur_uset = find(~isnan(Robs_mat(:,cell_num)));

cur_X{1} = all_Xmat_up_fixcor(cur_uset,use_kInds_up);
cur_X{2} = Xblock(used_inds(cur_uset),:);
cur_X{3} = Xsac(cur_uset,:);
cur_X{4} = Xmsac(cur_uset,:);

it_ref_mods{cc}{nn}(cell_num) = all_mod_fits(cell_num);
it_ref_mods{cc}{nn}(cell_num) = NMMfit_filters(it_ref_mods{cc}{nn}(cell_num),Robs_mat(cur_uset,cell_num),...
    cur_X,[],[],silent); %fit stimulus filters

% it_ref_LL_imp{cc}(nn,cell_num) = (it_ref_mods{cc}{nn}(cell_num).LL_seq(end) - all_nullmod(cell_num).LL_seq(end))/log(2);



