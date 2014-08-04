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

anal_name = 'monoc_eyecorr_hbar_LOOxv';
mod_data_name = 'monoc_eyecorr_hbar_mods_fin3';
% mod_data_name = 'monoc_eyecorr_hbar_mods_v4';

recompute_init_mods = 0;
use_measured_pos = 0;
use_sac_kerns = 1;

%dont fit stim models using these blocks
ignore_blocks = [16 17 28 30]; %G086
%%
xv_frac = 0.25;

n_fix_inf_it = 5;
n_drift_inf_it = 3;
fix_prior_sigma = 0.175;
drift_sigma = 0.01; %.015
drift_jump_sigma = 0.06; %0.05 start

stim_fs = 100; %in Hz
dt = 0.01;
flen = 12;
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
klen_us = use_nPix_us*flen;

sac_backlag = round(0.05/dt);
sac_forlag = round(0.3/dt);
sac_bincents = -sac_backlag:sac_forlag;
n_sac_bins = length(sac_bincents);

sp_dx = 0.0565/spatial_usfac;
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

%LIST OF BLOCKS TO EXCLUDE
if strcmp(Expt_name,'G086')
    %     cur_block_set(cur_block_set == 17) = []; %fails to align stim rc files.
    %     cur_block_set(cur_block_set == 28) = []; %clearly problem with the stimulus but don't know what it is!
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
all_trial_blocknums = [];
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
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
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

if use_sac_kerns
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
end

%% INCORPORATE MEASURED EYE-POSITIONS
if use_measured_pos
    fprintf('Incorporating measured eye-corrections\n');
    if bar_ori == 0
        measured_orth_pos = round(corrected_interp_eyevals(used_inds,2)/sp_dx);
    elseif bar_ori == 90
        measured_orth_pos = round(corrected_interp_eyevals(used_inds,1)/sp_dx);
    end
    max_meas_shift = 30;
    measured_orth_pos(measured_orth_pos > max_meas_shift) = max_meas_shift;
    measured_orth_pos(measured_orth_pos < -max_meas_shift) = -max_meas_shift;
    measured_orth_pos_ds = round(measured_orth_pos/spatial_usfac);
    for ii = 1:NT
        all_stimmat_up(used_inds(ii),:) = shift_matrix_Nd(all_stimmat_up(used_inds(ii),:),-measured_orth_pos(ii),2);
        all_stim_mat(used_inds(ii),:) = shift_matrix_Nd(all_stim_mat(used_inds(ii),:),-measured_orth_pos_ds(ii),2);
    end
    all_Xmat_us = create_time_embedding(all_stimmat_up,stim_params_us);
    all_Xmat = create_time_embedding(all_stim_mat,stim_params);
end

%% Create set of TR and XV trials
use_trials = unique(all_trialvec(used_inds));
nuse_trials = length(use_trials);

poss_tr_trials = find(~ismember(cur_block_set(all_trial_blocknums),ignore_blocks));
nposs_tr_trials = length(poss_tr_trials);
n_xv_trials = round(xv_frac*nposs_tr_trials);
xv_trials = randperm(nposs_tr_trials);
xv_trials(n_xv_trials+1:end) = [];
xv_trials = poss_tr_trials(xv_trials);
tr_trials = setdiff(poss_tr_trials,xv_trials);

tr_inds = find(ismember(all_trialvec(used_inds),tr_trials));
xv_inds = find(ismember(all_trialvec(used_inds),xv_trials));

tr_blocks = find(~ismember(cur_block_set,ignore_blocks));

%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITS
init_stim_params = NMMcreate_stim_params([flen use_nPix],dt);
init_stim_params(2) = NMMcreate_stim_params([n_blocks 1],dt);

fin_stim_params(1) = NMMcreate_stim_params([flen use_nPix_us],dt);
fin_stim_params(2) = NMMcreate_stim_params([n_blocks 1],dt);
if use_sac_kerns
    fin_stim_params(3) = NMMcreate_stim_params(n_sac_bins,dt);
    fin_stim_params(4) = NMMcreate_stim_params(n_sac_bins,dt);
end

null_stim_params = fin_stim_params(2:end);

block_L2 = 1;
silent = 1;
sac_d2t = 100;

if expt_dds(cur_block_set(1)) == 67
    base_lambda_d2XT = 60;
    base_lambda_L1 = 20;
else
    base_lambda_d2XT = 10;
    base_lambda_L1 = 10;
end
init_optim_p.optTol = 1e-5; init_optim_p.progTol = 1e-9;

sac_reg_params = NMMcreate_reg_params('lambda_d2T',sac_d2t);
if use_sac_kerns
    null_reg_params = NMMcreate_reg_params('lambda_d2T',[0; sac_d2t; sac_d2t],'lambda_L2',[block_L2; 0; 0]);
else
    null_reg_params = NMMcreate_reg_params('lambda_L2',block_L2);
end

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
            cur_used_blocks = tr_blocks; %blocks when NO SU
            cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
        end
        cur_tr_inds = tr_inds(ismember(used_inds(tr_inds),cur_poss_inds));
        cur_xv_inds = xv_inds(ismember(used_inds(xv_inds),cur_poss_inds));
        tr_NT = length(cur_tr_inds);
        xv_NT = length(cur_xv_inds);
        if ~isempty(cur_tr_inds)
            Robs = all_binned_spikes(used_inds(cur_tr_inds),ss);
            Robsxv = all_binned_spikes(used_inds(cur_xv_inds),ss);
            
            tr_X{1} = all_Xmat(used_inds(cur_tr_inds),use_kInds);
            tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
            xv_X{1} = all_Xmat(used_inds(cur_xv_inds),use_kInds);
            xv_X{2} = Xblock(used_inds(cur_xv_inds),:);
            
            if use_sac_kerns
                tr_X{3} = Xsac(cur_tr_inds,:);
                tr_X{4} = Xmsac(cur_tr_inds,:);
                xv_X{3} = Xsac(cur_xv_inds,:);
                xv_X{4} = Xmsac(cur_xv_inds,:);
                null_mod = NMMinitialize_model(null_stim_params,[1 1 1],repmat({'lin'},1,3),null_reg_params,[1 2 3]);
            else
                null_mod = NMMinitialize_model(null_stim_params,[1],{'lin'},null_reg_params,[1]);
            end
            
            null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));
            all_nullmod(ss) = null_mod;
            
            gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
            gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent,init_optim_p);
            
            tr_X{1} = all_Xmat_us(used_inds(cur_tr_inds),use_kInds_up);
            xv_X{1} = all_Xmat_us(used_inds(cur_xv_inds),use_kInds_up);
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
            [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm2,Robs,tr_X);
            gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
            gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
            if use_sac_kerns
                gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,3,zeros(n_sac_bins,1),sac_reg_params); %sac term
                gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,4,zeros(n_sac_bins,1),sac_reg_params); %nsac term
            end
            gqm2 = NMMfit_filters(gqm2,Robs,tr_X,[],[],silent);
            
            all_mod_fits(ss) = gqm2;
            all_mod_fits_withspkNL(ss) = NMMfit_logexp_spkNL(gqm2,Robs,tr_X);
            
            LL = NMMmodel_eval(all_mod_fits_withspkNL(ss),Robs,tr_X);
            pLL = NMMmodel_eval(all_mod_fits(ss),Robs,tr_X);
            null_LL(ss) = NMMmodel_eval(null_mod,Robs,tr_X(2:end));
            all_mod_LLimp(ss) = (LL-null_LL(ss))/log(2);
            all_mod_pLLimp(ss) = (pLL-null_LL(ss))/log(2);
            
            all_xvmod_fits(ss) = NMMfit_filters(gqm2,Robsxv,xv_X,[],[],silent);
            all_xvmod_fits_withspkNL(ss) = NMMfit_logexp_spkNL(all_xvmod_fits(ss),Robsxv,xv_X);
            xvLL = NMMmodel_eval(all_xvmod_fits_withspkNL(ss),Robsxv,xv_X);
            truexvLL = NMMmodel_eval(all_mod_fits_withspkNL(ss),Robsxv,xv_X);
            pxvLL = NMMmodel_eval(all_xvmod_fits(ss),Robsxv,xv_X);
            ptruexvLL = NMMmodel_eval(all_mod_fits(ss),Robsxv,xv_X);
            null_xvLL(ss) = NMMmodel_eval(null_mod,Robsxv,xv_X(2:end));
            all_mod_xvLLimp(ss) = (xvLL-null_xvLL(ss))/log(2);
            all_mod_truexvLLimp(ss) = (truexvLL-null_xvLL(ss))/log(2);
            all_mod_pxvLLimp(ss) = (pxvLL-null_xvLL(ss))/log(2);
            all_mod_ptruexvLLimp(ss) = (ptruexvLL-null_xvLL(ss))/log(2);
        end
    end
    
    for ss = 1:length(su_probes);
        fprintf('Computing base LLs for SU %d of %d\n',ss,length(su_probes));
        all_mod_SU(ss+96) = su_probes(ss);
        cur_used_blocks = tr_blocks(su_used_blocks(tr_blocks,ss)); %blocks when SU
        cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
        cur_tr_inds = tr_inds(ismember(used_inds(tr_inds),cur_poss_inds));
        cur_xv_inds = xv_inds(ismember(used_inds(xv_inds),cur_poss_inds));
        tr_NT = length(cur_tr_inds);
        xv_NT = length(cur_xv_inds);
        if ~isempty(cur_tr_inds)
            Robs = all_binned_spikes(used_inds(cur_tr_inds),su_probes(ss));
            Robsxv = all_binned_spikes(used_inds(cur_xv_inds),su_probes(ss));
            
            tr_X{1} = all_Xmat(used_inds(cur_tr_inds),use_kInds);
            tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
            xv_X{1} = all_Xmat(used_inds(cur_xv_inds),use_kInds);
            xv_X{2} = Xblock(used_inds(cur_xv_inds),:);
            if use_sac_kerns
                tr_X{3} = Xsac(cur_tr_inds,:);
                tr_X{4} = Xmsac(cur_tr_inds,:);
                xv_X{3} = Xsac(cur_xv_inds,:);
                xv_X{4} = Xmsac(cur_xv_inds,:);
                null_mod = NMMinitialize_model(null_stim_params,[1 1 1],repmat({'lin'},1,3),null_reg_params,[1 2 3]);
            else
                null_mod = NMMinitialize_model(null_stim_params,[1],{'lin'},null_reg_params,[1]);
            end
            null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));
            all_nullmod(ss+96) = null_mod;
            
            gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
            gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent,init_optim_p);
            
            tr_X{1} = all_Xmat_us(used_inds(cur_tr_inds),use_kInds_up);
            xv_X{1} = all_Xmat_us(used_inds(cur_xv_inds),use_kInds_up);
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
            [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm2,Robs,tr_X);
            gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
            gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
            if use_sac_kerns
                gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,3,zeros(n_sac_bins,1),sac_reg_params); %sac term
                gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,4,zeros(n_sac_bins,1),sac_reg_params); %nsac term
            end
            gqm2 = NMMfit_filters(gqm2,Robs,tr_X,[],[],silent);
            
            all_mod_fits(ss+96) = gqm2;
            
            all_mod_fits_withspkNL(ss+96) = NMMfit_logexp_spkNL(gqm2,Robs,tr_X);
            
            LL = NMMmodel_eval(all_mod_fits_withspkNL(ss+96),Robs,tr_X);
            pLL = NMMmodel_eval(all_mod_fits(ss+96),Robs,tr_X);
            null_LL(ss+96) = NMMmodel_eval(null_mod,Robs,tr_X(2:end));
            all_mod_LLimp(ss+96) = (LL-null_LL(ss+96))/log(2);
            all_mod_pLLimp(ss+96) = (pLL-null_LL(ss+96))/log(2);
            
            all_xvmod_fits(ss+96) = NMMfit_filters(gqm2,Robsxv,xv_X,[],[],silent);
            all_xvmod_fits_withspkNL(ss+96) = NMMfit_logexp_spkNL(all_xvmod_fits(ss+96),Robsxv,xv_X);
            xvLL = NMMmodel_eval(all_xvmod_fits_withspkNL(ss+96),Robsxv,xv_X);
            truexvLL = NMMmodel_eval(all_mod_fits_withspkNL(ss+96),Robsxv,xv_X);
            pxvLL = NMMmodel_eval(all_xvmod_fits(ss+96),Robsxv,xv_X);
            ptruexvLL = NMMmodel_eval(all_mod_fits(ss+96),Robsxv,xv_X);
            null_xvLL(ss+96) = NMMmodel_eval(null_mod,Robsxv,xv_X(2:end));
            all_mod_xvLLimp(ss+96) = (xvLL-null_xvLL(ss+96))/log(2);
            all_mod_truexvLLimp(ss+96) = (truexvLL-null_xvLL(ss+96))/log(2);
            all_mod_pxvLLimp(ss+96) = (xvLL-null_xvLL(ss+96))/log(2);
            all_mod_ptruexvLLimp(ss+96) = (truexvLL-null_xvLL(ss+96))/log(2);
            
        end
    end
    if xv_frac > 0
        save(mod_data_name,'all_mod*','all_nullmod','all_xv*','su_probes','null_xvLL','null_LL','*_trials');
    else
        save(mod_data_name,'all_mod*','all_nullmod','su_probes','null_xvLL','null_LL','*_trials');
    end
else
    fprintf('Loading pre-computed initial models\n');
    load(mod_data_name);
    
    tr_inds = find(ismember(all_trialvec(used_inds),tr_trials));
    xv_inds = find(ismember(all_trialvec(used_inds),xv_trials));
end

%% SELECT USABLE UNITS AND make Robs_mat
if xv_frac == 0
    LL_imp_thresh = 5e-3;
    usable_units = find(all_mod_LLimp >= LL_imp_thresh);
else
    LL_imp_thresh = 1e-3;
    usable_units = find(all_mod_truexvLLimp >= LL_imp_thresh);
end
n_used_sus = sum(all_mod_SU(usable_units) ~= 0);
n_used_mus = sum(all_mod_SU(usable_units) == 0);
fprintf('Using %d SUs and %d MUs for analysis\n',n_used_sus,n_used_mus);
tr_set = usable_units;
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
%% DIAGNOSTICS
cur_X{1} = all_Xmat_us(used_inds,use_kInds_up);
cur_X{2} = Xblock(used_inds,:);
if use_sac_kerns
    cur_X{3} = Xsac;
    cur_X{4} = Xmsac;
end
full_prates = nan(NT,n_tr_chs);
% full_xvprates = nan(NT,n_tr_chs);
full_nullrates = nan(NT,n_tr_chs);
for cc = 1:n_tr_chs
    [~,~,full_prates(:,cc)] = NMMmodel_eval(all_mod_fits_withspkNL(tr_set(cc)),Robs_mat(:,cc),cur_X);
    %     [~,~,full_xvprates(:,cc)] = NMMmodel_eval(all_xvmod_fits(tr_set(cc)),Robs_mat(:,cc),cur_X);
    [~,~,full_nullrates(:,cc)] = NMMmodel_eval(all_nullmod(tr_set(cc)),Robs_mat(:,cc),cur_X(2:end));
end
full_modLL = Robs_mat.*log(full_prates) - full_prates;
full_nullLL = Robs_mat.*log(full_nullrates) - full_nullrates;
full_LLimp = full_modLL-full_nullLL;

trial_blocknums = nan(nuse_trials,1);
trial_LLimp = nan(nuse_trials,n_tr_chs);
trial_meanrate = nan(nuse_trials,n_tr_chs);
trial_nspks = nan(nuse_trials,n_tr_chs);
trial_durs = nan(nuse_trials,1);
trial_eyeerr = nan(nuse_trials,1);
for tt = 1:nuse_trials
    cur_used_inds = find(all_trialvec(used_inds) == use_trials(tt));
    trial_LLimp(tt,:) = nanmean(full_LLimp(cur_used_inds,:));
    trial_meanrate(tt,:) = nanmean(Robs_mat(cur_used_inds,:));
    trial_nspks(tt,:) = nansum(Robs_mat(cur_used_inds,:));
    trial_blocknums(tt) = unique(all_blockvec(used_inds(cur_used_inds)));
    trial_durs(tt) = length(cur_used_inds);
    trial_eyeerr(tt) = mean(corrected_interp_eyevals(used_inds(cur_used_inds),2));
end

block_LLimp = nan(n_blocks,n_tr_chs);
block_meanrate = nan(n_blocks,n_tr_chs);
block_nspks = nan(n_blocks,n_tr_chs);
for tt = 1:n_blocks
    cur_used_inds = find(all_blockvec(used_inds) == tt);
    block_LLimp(tt,:) = nanmean(full_LLimp(cur_used_inds,:));
    block_meanrate(tt,:) = nanmean(Robs_mat(cur_used_inds,:));
    block_nspks(tt,:) = nansum(Robs_mat(cur_used_inds,:));
end

tr_trialset = find(ismember(use_trials,tr_trials));
xv_trialset = find(ismember(use_trials,xv_trials));

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

leps_prior = -(shifts*sp_dx).^2/(2*fix_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize

%generate shift matrices. Must be applied to the stimulus (not the filters)
It = speye(flen);
shift_mat = cell(n_shifts,1);
for xx = 1:n_shifts
    temp = spdiags( ones(full_nPix_us,1), -shifts(xx), full_nPix_us, full_nPix_us);
    temp = kron(temp,It);
    shift_mat{xx} = temp(:,use_kInds_up);
end

%generate shift matrices. Must be applied to the stimulus (not the filters)
It = speye(flen);
filt_shift_mat = cell(n_shifts,1);
for xx = 1:n_shifts
    temp = spdiags( ones(use_nPix_us,1), shifts(xx), use_nPix_us, use_nPix_us);
    temp = kron(temp,It);
    filt_shift_mat{xx} = temp;
end

%% ITERATE FIXATION-BASED CORRECTIONS
fprintf('Inferring fixation corrs using all units\n');
%% PREPROCESS MODEL COMPONENTS
n_tr_chs = length(tr_set);

filt_bank = zeros(n_tr_chs,klen_us,n_squared_filts+1);
lin_kerns = nan(n_tr_chs,n_blocks);
if use_sac_kerns
    sac_kerns = nan(n_tr_chs,n_sac_bins);
    msac_kerns = nan(n_tr_chs,n_sac_bins);
end
mod_spkNL_params = nan(n_tr_chs,3);
for ss = 1:n_tr_chs
    cur_Xtargs = [all_mod_fits(tr_set(ss)).mods(:).Xtarget];
    cur_k = [all_mod_fits(tr_set(ss)).mods(cur_Xtargs == 1).filtK];
    n_used_filts = size(cur_k,2);
    filt_bank(ss,:,1:n_used_filts) = cur_k;
    mod_spkNL_params(ss,:) = all_mod_fits_withspkNL(tr_set(ss)).spk_NL_params;
    lin_kerns(ss,:) = all_mod_fits(tr_set(ss)).mods(cur_Xtargs == 2).filtK;
    if use_sac_kerns
        sac_kerns(ss,:) = all_mod_fits(tr_set(ss)).mods(cur_Xtargs == 3).filtK;
        msac_kerns(ss,:) = all_mod_fits(tr_set(ss)).mods(cur_Xtargs == 4).filtK;
    end
end
filt_bank = permute(filt_bank,[2 1 3]);
filt_bank = reshape(filt_bank,[flen use_nPix_us n_tr_chs n_squared_filts+1]);

%indicator predictions
block_out = Xblock(used_inds,:)*lin_kerns';
if use_sac_kerns
    sac_out = Xsac*sac_kerns';
    msac_out = Xmsac*msac_kerns';
end
%% ESTIMATE LL for each shift in each stimulus frame
cur_Xmat = all_Xmat_us(used_inds,use_kInds_up);

cur_filt_shift = nan(klen_us,n_tr_chs,n_squared_filts+1);
%precompute LL at all shifts for all units
frame_LLs = nan(NT,n_shifts);
for xx = 1:length(shifts)
    fprintf('Shift %d of %d\n',xx,n_shifts);
    cur_filt_shift = shift_matrix_Nd(filt_bank,shifts(xx),2);
    cur_filt_shift = reshape(cur_filt_shift,[klen_us n_tr_chs n_squared_filts+1]);
    
    %outputs of stimulus models at current X-matrix shift
    gfuns = ones(NT,n_tr_chs);
    gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,1)');
    gfuns = gfuns + cur_Xmat*squeeze(cur_filt_shift(:,:,1));
    for ff = 2:(n_squared_filts+1)
        gfuns = gfuns + (cur_Xmat*squeeze(cur_filt_shift(:,:,ff))).^2;
    end
    
    %add contributions from extra lin kernels
    gfuns = gfuns + block_out;
    if use_sac_kerns
        gfuns = gfuns + sac_out + msac_out;
    end
    
    %incorporate beta
    gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,2)');
    
    %handle numerical overflow with log(1+exp)
    too_large = gfuns > 50;
    pred_rate = log(1+exp(gfuns));
    pred_rate(too_large) = gfuns(too_large);
    
    %incorporate alpha
    pred_rate = bsxfun(@times,pred_rate,mod_spkNL_params(:,3)');
    
    %enforce min predicted rate
    pred_rate(pred_rate < 1e-50) = 1e-50;
    
    frame_LLs(:,xx) = squeeze(nansum(Robs_mat.*log(pred_rate) - pred_rate,2));
end

%% INFER MICRO-SAC SEQUENCE
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
fix_post_mean = sum(bsxfun(@times,Post,shifts),2);
cur_diff = bsxfun(@minus,fix_post_mean,shifts).^2;
fix_post_std = sqrt(sum(cur_diff.*Post,2));

%back-project saccade-times
all_fix_post_mean_cor = fix_post_mean(fix_ids);
%align jumps to saccade times
for i = length(saccade_start_inds):-1:1
    cur_inds = saccade_start_inds(i):(saccade_start_inds(i)+sac_shift);
    cur_inds(cur_inds > NT) = [];
    all_fix_post_mean_cor(cur_inds) = all_fix_post_mean_cor(cur_inds(end));
    cur_back_inds = (saccade_start_inds(i)-(flen-sac_shift)):(saccade_start_inds(i)-1);
    cur_back_inds(cur_back_inds < 1) = [];
    all_fix_post_mean_cor(cur_back_inds) = all_fix_post_mean_cor(saccade_start_inds(i));
end

%% RECOMPUTE XMAT
all_shift_stimmat_up = all_stimmat_up;
for i=1:NT
    all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-round(all_fix_post_mean_cor(i)),2);
end
all_Xmat_up_fixcor = create_time_embedding(all_shift_stimmat_up,stim_params_us);
all_Xmat_up_fixcor = all_Xmat_up_fixcor(used_inds,:);

%% REFIT XV CELLS
cur_X{1} = all_Xmat_up_fixcor(:,use_kInds_up);
cur_X{2} = Xblock(used_inds,:);
if use_sac_kerns
    cur_X{3} = Xsac;
    cur_X{4} = Xmsac;
end

silent = 1;
for ss = 1:length(tr_set)
    cur_cell = tr_set(ss);
    fprintf('Refitting model for tr cell %d of %d\n',ss,length(tr_set));
    cur_unit_ind = find(tr_set == cur_cell);
    cur_tr_uset = tr_inds(~isnan(Robs_mat(tr_inds,cur_unit_ind)));
    cur_xv_uset = xv_inds(~isnan(Robs_mat(xv_inds,cur_unit_ind)));
    
    tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
    xv_X = get_Xcell_tInds(cur_X,cur_xv_uset);
    
    ref_mods(cur_cell) = all_mod_fits(cur_cell);
    ref_mods(cur_cell) = NMMfit_filters(ref_mods(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),...
        tr_X,[],[],silent); %fit stimulus filters
    
    %refit spk NL
    ref_mods_spkNL(cur_cell) = NMMfit_logexp_spkNL(ref_mods(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
    
    newLL = NMMmodel_eval(ref_mods_spkNL(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
    
    ref_xvmods(cur_cell) = all_xvmod_fits(cur_cell);
    ref_xvmods(cur_cell) = NMMfit_filters(ref_xvmods(cur_cell),Robs_mat(cur_xv_uset,cur_unit_ind),...
        xv_X,[],[],silent); %fit stimulus filters
    %refit spk NL
    ref_xvmods_spkNL(cur_cell) = NMMfit_logexp_spkNL(ref_xvmods(cur_cell),Robs_mat(cur_xv_uset,cur_unit_ind),xv_X);
    
    newxvLL = NMMmodel_eval(ref_xvmods_spkNL(cur_cell),Robs_mat(cur_xv_uset,cur_unit_ind),xv_X);
    
    it_LLimp(cur_cell) = (newLL - null_LL(cur_cell))/log(2);
    it_xvLLimp(cur_cell) = (newxvLL - null_xvLL(cur_cell))/log(2);
    
    fprintf('Original: %.4f  New: %.4f\n',all_mod_LLimp(cur_cell),it_LLimp(cur_cell));
    fprintf('XV Original: %.4f  New: %.4f\n',all_mod_xvLLimp(cur_cell),it_xvLLimp(cur_cell));
end

%%
cur_Xmat = all_Xmat_us(used_inds,use_kInds_up);

for nn = 1:length(tr_set)
    fprintf('Inferring fixation corrections, XV %d of %d\n',nn,length(tr_set));
    
    %% PREPROCESS MODEL COMPONENTS
    cur_xv_cells = tr_set(nn);
    cur_tr_cells = setdiff(tr_set,cur_xv_cells);
    cur_tr_cinds = find(ismember(tr_set,cur_tr_cells));
    cur_xv_cinds = find(ismember(tr_set,cur_xv_cells));
    cur_n_tr_chs = length(cur_tr_cells);
    
    filt_bank = zeros(n_squared_filts+1,klen_us);
    cur_Xtargs = [all_mod_fits(cur_xv_cells).mods(:).Xtarget];
    cur_k = [all_mod_fits(cur_xv_cells).mods(cur_Xtargs == 1).filtK];
    n_used_filts = size(cur_k,2);
    filt_bank(1:n_used_filts,:) = cur_k';
    mod_spkNL_params = all_mod_fits_withspkNL(cur_xv_cells).spk_NL_params;
    lin_kerns = all_mod_fits(cur_xv_cells).mods(cur_Xtargs == 2).filtK;
    if use_sac_kerns
        sac_kerns = all_mod_fits(cur_xv_cells).mods(cur_Xtargs == 3).filtK;
        msac_kerns = all_mod_fits(cur_xv_cells).mods(cur_Xtargs == 4).filtK;
    end
    
    %indicator predictions
    block_out = Xblock(used_inds,:)*lin_kerns;
    if use_sac_kerns
        sac_out = Xsac*sac_kerns;
        msac_out = Xmsac*msac_kerns;
    end
    %% ESTIMATE LL for each shift in each stimulus frame
    
    %precompute LL at all shifts for all units
    xv_frame_LLs = nan(NT,n_shifts);
    for xx = 1:length(shifts)
        %         fprintf('Shift %d of %d\n',xx,n_shifts);
        %         cur_stim_shift = cur_Xmat*shift_mat{xx};
        cur_filt_shift = (filt_bank*filt_shift_mat{xx})';
        
        %outputs of stimulus models at current X-matrix shift
        gfuns = ones(NT,1)*mod_spkNL_params(1);
        gfuns = gfuns + cur_Xmat*cur_filt_shift(:,1);
        for ff = 2:(n_squared_filts+1)
            gfuns = gfuns + (cur_Xmat*cur_filt_shift(:,ff)).^2;
        end
        
        %add contributions from extra lin kernels
        gfuns = gfuns + block_out;
        if use_sac_kerns
            gfuns = gfuns + sac_out + msac_out;
        end
        
        %incorporate beta
        gfuns = gfuns*mod_spkNL_params(2);
        
        %handle numerical overflow with log(1+exp)
        too_large = gfuns > 50;
        pred_rate = log(1+exp(gfuns));
        pred_rate(too_large) = gfuns(too_large);
        
        %incorporate alpha
        pred_rate = pred_rate*mod_spkNL_params(3);
        
        %enforce min predicted rate
        pred_rate(pred_rate < 1e-50) = 1e-50;
        
        xv_frame_LLs(:,xx) = squeeze(nansum(Robs_mat(:,cur_xv_cinds).*log(pred_rate) - pred_rate,2));
    end
    
    %% INFER MICRO-SAC SEQUENCE
    cur_frame_LLs = frame_LLs - xv_frame_LLs;
    
    fix_LLs = nan(n_fixs,n_shifts);
    fix_ids = nan(NT,1);
    for ii = 1:n_fixs
        cur_inds = jump_inds(ii):(jump_inds(ii+1)-1);
        fix_LLs(ii,:) = sum(cur_frame_LLs(cur_inds,:));
        fix_ids(cur_inds) = ii;
    end
    
    lPost = bsxfun(@plus,fix_LLs,leps_prior);
    lPost = bsxfun(@minus,lPost,logsumexp(lPost,2));
    Post = exp(lPost);
    xv_fix_post_mean(nn,:) = sum(bsxfun(@times,Post,shifts),2);
    cur_diff = bsxfun(@minus,xv_fix_post_mean(nn,:)',shifts).^2;
    xv_fix_post_std(nn,:) = sqrt(sum(cur_diff.*Post,2));
    
    %back-project saccade-times
    all_fix_post_mean_cor = xv_fix_post_mean(nn,fix_ids);
    %align jumps to saccade times
    for i = length(saccade_start_inds):-1:1
        cur_inds = saccade_start_inds(i):(saccade_start_inds(i)+sac_shift);
        cur_inds(cur_inds > NT) = [];
        all_fix_post_mean_cor(cur_inds) = all_fix_post_mean_cor(cur_inds(end));
        cur_back_inds = (saccade_start_inds(i)-(flen-sac_shift)):(saccade_start_inds(i)-1);
        cur_back_inds(cur_back_inds < 1) = [];
        all_fix_post_mean_cor(cur_back_inds) = all_fix_post_mean_cor(saccade_start_inds(i));
    end
    
    %% RECOMPUTE XMAT
    all_shift_stimmat_up = all_stimmat_up;
    for i=1:NT
        all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-round(all_fix_post_mean_cor(i)),2);
    end
    all_Xmat_up_fixcor = create_time_embedding(all_shift_stimmat_up,stim_params_us);
    all_Xmat_up_fixcor = all_Xmat_up_fixcor(used_inds,:);
    
    %% REFIT XV CELLS
    cur_X{1} = all_Xmat_up_fixcor(:,use_kInds_up);
    
    silent = 1;
    cur_cell = cur_xv_cells;
    fprintf('Refitting model \n');
    cur_unit_ind = find(tr_set == cur_cell);
    cur_tr_uset = tr_inds(~isnan(Robs_mat(tr_inds,cur_unit_ind)));
    cur_xv_uset = xv_inds(~isnan(Robs_mat(xv_inds,cur_unit_ind)));
    
    tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
    xv_X = get_Xcell_tInds(cur_X,cur_xv_uset);
    
    xv_mods(cur_cell) = all_mod_fits(cur_cell);
    xv_mods(cur_cell) = NMMfit_filters(xv_mods(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),...
        tr_X,[],[],silent); %fit stimulus filters
    
    %refit spk NL
    xv_mods_spkNL(cur_cell) = NMMfit_logexp_spkNL(xv_mods(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
    
    newLL = NMMmodel_eval(xv_mods_spkNL(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
    pnewLL = NMMmodel_eval(xv_mods(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
    
    xv_xvmods(cur_cell) = all_xvmod_fits(cur_cell);
    xv_xvmods(cur_cell) = NMMfit_filters(xv_xvmods(cur_cell),Robs_mat(cur_xv_uset,cur_unit_ind),...
        xv_X,[],[],silent); %fit stimulus filters
    %refit spk NL
    xv_xvmods_spkNL(cur_cell) = NMMfit_logexp_spkNL(xv_xvmods(cur_cell),Robs_mat(cur_xv_uset,cur_unit_ind),xv_X);
    
    newxvLL = NMMmodel_eval(xv_xvmods_spkNL(cur_cell),Robs_mat(cur_xv_uset,cur_unit_ind),xv_X);
    
    xv_LLimp(cur_cell) = (newLL - null_LL(cur_cell))/log(2);
    xv_xvLLimp(cur_cell) = (newxvLL - null_xvLL(cur_cell))/log(2);
    
end

%%
leps_prior = -(shifts*sp_dx).^2/(2*drift_jump_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
lA_tflip = repmat(leps_prior,n_shifts,1);

cdist = squareform(pdist(shifts'*sp_dx));
lA = -cdist.^2/(2*drift_sigma^2);
lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize
fprintf('Inferring drift corrections\n');

%% RECOMPUTE XMAT
all_fix_post_mean_cor = fix_post_mean(fix_ids);
all_shift_stimmat_up = all_stimmat_up;
for i=1:NT
    all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-round(all_fix_post_mean_cor(i)),2);
end
all_Xmat_up_fixcor = create_time_embedding(all_shift_stimmat_up,stim_params_us);
all_Xmat_up_fixcor = all_Xmat_up_fixcor(used_inds,:);

%% PREPROCESS MODEL COMPONENTS
filt_bank = zeros(n_tr_chs,klen_us,n_squared_filts+1);
lin_kerns = nan(n_tr_chs,n_blocks);
if use_sac_kerns
    sac_kerns = nan(n_tr_chs,n_sac_bins);
    msac_kerns = nan(n_tr_chs,n_sac_bins);
end
mod_spkNL_params = nan(n_tr_chs,3);
for ss = 1:n_tr_chs
    cur_Xtargs = [ref_mods(tr_set(ss)).mods(:).Xtarget];
    cur_k = [ref_mods(tr_set(ss)).mods(cur_Xtargs == 1).filtK];
    n_used_filts = size(cur_k,2);
    filt_bank(ss,:,1:n_used_filts) = cur_k;
    mod_spkNL_params(ss,:) = ref_mods_spkNL(tr_set(ss)).spk_NL_params;
    lin_kerns(ss,:) = ref_mods(tr_set(ss)).mods(cur_Xtargs == 2).filtK;
    if use_sac_kerns
        sac_kerns(ss,:) = ref_mods(tr_set(ss)).mods(cur_Xtargs == 3).filtK;
        msac_kerns(ss,:) = ref_mods(tr_set(ss)).mods(cur_Xtargs == 4).filtK;
    end
end
filt_bank = permute(filt_bank,[2 1 3]);
filt_bank = reshape(filt_bank,[flen use_nPix_us n_tr_chs n_squared_filts+1]);

%indicator predictions
block_out = Xblock(used_inds,:)*lin_kerns';
if use_sac_kerns
    sac_out = Xsac*sac_kerns';
    msac_out = Xmsac*msac_kerns';
end
%% ESTIMATE LL for each shift in each stimulus frame

cur_Xmat = all_Xmat_up_fixcor(:,use_kInds_up);

%precompute LL at all shifts for all units
frame_LLs = nan(NT,n_shifts);
for xx = 1:length(shifts)
    fprintf('Shift %d of %d\n',xx,n_shifts);
    cur_filt_shift = shift_matrix_Nd(filt_bank,shifts(xx),2);
    cur_filt_shift = reshape(cur_filt_shift,[klen_us n_tr_chs n_squared_filts+1]);
    
    gfuns = ones(NT,1);
    gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,1)');
    gfuns = gfuns + cur_Xmat*squeeze(cur_filt_shift(:,:,1));
    for ff = 2:(n_squared_filts+1)
        gfuns = gfuns + (cur_Xmat*squeeze(cur_filt_shift(:,:,ff))).^2;
    end
    
    %add contributions from extra lin kernels
    gfuns = gfuns + block_out;
    if use_sac_kerns
        gfuns = gfuns + sac_out + msac_out;
    end
    
    %incorporate beta
    gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,2)');
    
    %handle numerical overflow with log(1+exp)
    too_large = gfuns > 50;
    pred_rate = log(1+exp(gfuns));
    pred_rate(too_large) = gfuns(too_large);
    
    %incorporate alpha
    pred_rate = bsxfun(@times,pred_rate,mod_spkNL_params(:,3)');
    
    %enforce min predicted rate
    pred_rate(pred_rate < 1e-50) = 1e-50;
    
    frame_LLs(:,xx) = squeeze(nansum(Robs_mat.*log(pred_rate) - pred_rate,2));
end


%% INFER DRIFT CORRECTIONS

lalpha=zeros(NT,n_shifts);
lbeta = zeros(NT,n_shifts);
lscale=zeros(NT,1); %initialize rescaling parameters
%compute rescaled forward messages
lalpha(1,:) = leps_prior + frame_LLs(1,:);
lscale(1)=logsumexp(lalpha(1,:));
lalpha(1,:) = lalpha(1,:) - lscale(1);
for t=2:NT
    if use_prior(t)
        cur_lA = lA_tflip;
    else
        cur_lA = lA;
    end
    lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + frame_LLs(t,:);
    lscale(t) = logsumexp(lalpha(t,:));
    lalpha(t,:) = lalpha(t,:) - lscale(t);
end

%compute rescaled backward messages
lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
for t=NT-1:-1:1
    if use_prior(t+1)
        cur_lA = lA_tflip;
    else
        cur_lA = lA;
    end
    lf1 = lbeta(t+1,:) + frame_LLs(t+1,:);
    lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
end
lgamma= lalpha + lbeta;
lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));

gamma = exp(lgamma);
drift_post_mean = sum(bsxfun(@times,gamma,shifts),2);
drift_post_std = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - drift_post_mean.^2);

%back-project saccade-times
drift_post_mean_cor = drift_post_mean;

%% RECOMPUTE XMAT

drift_Xmat = reshape(all_Xmat_up_fixcor,[NT flen full_nPix_us]);
for ii = 1:NT
    drift_Xmat(ii,:,:) = shift_matrix_Nd(drift_Xmat(ii,:,:),-round(drift_post_mean_cor(ii)),3);
end
drift_Xmat = reshape(drift_Xmat,[NT flen*full_nPix_us]);
%% REFIT ALL CELLS
cur_X{1} = drift_Xmat(:,use_kInds_up);
cur_X{2} = Xblock(used_inds,:);
if use_sac_kerns
    cur_X{3} = Xsac;
    cur_X{4} = Xmsac;
end

silent = 1;
for ss = 1:length(tr_set)
    cur_cell = tr_set(ss);
    fprintf('Refitting model for tr cell %d of %d\n',ss,length(tr_set));
    cur_unit_ind = find(tr_set == cur_cell);
    cur_tr_uset = tr_inds(~isnan(Robs_mat(tr_inds,cur_unit_ind)));
    cur_xv_uset = xv_inds(~isnan(Robs_mat(xv_inds,cur_unit_ind)));
    
    tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
    xv_X = get_Xcell_tInds(cur_X,cur_xv_uset);
    
    dref_mods(cur_cell) = ref_mods(cur_cell);
    dref_mods(cur_cell) = NMMfit_filters(dref_mods(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),...
        tr_X,[],[],silent); %fit stimulus filters
    
    dref_mods_spkNL(cur_cell) = NMMfit_logexp_spkNL(dref_mods(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
    
    newLL = NMMmodel_eval(dref_mods_spkNL(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
    
    dref_xvmods(cur_cell) = ref_xvmods(cur_cell);
    dref_xvmods(cur_cell) = NMMfit_filters(dref_xvmods(cur_cell),Robs_mat(cur_xv_uset,cur_unit_ind),...
        xv_X,[],[],silent); %fit stimulus filters
    
    dref_xvmods_spkNL(cur_cell) = NMMfit_logexp_spkNL(dref_xvmods(cur_cell),Robs_mat(cur_xv_uset,cur_unit_ind),xv_X);
    
    newxvLL = NMMmodel_eval(dref_xvmods_spkNL(cur_cell),Robs_mat(cur_xv_uset,cur_unit_ind),xv_X);
    
    dit_LLimp(cur_cell) = (newLL - null_LL(cur_cell))/log(2);
    dit_xvLLimp(cur_cell) = (newxvLL - null_xvLL(cur_cell))/log(2);
    
    fprintf('Original: %.5f  Prev: %.5f  New: %.5f\n',all_mod_LLimp(cur_cell),it_LLimp(cur_cell),dit_LLimp(cur_cell));
    fprintf('XV Original: %.4f  Prev: %.4f  New: %.4f\n',all_mod_xvLLimp(cur_cell),it_xvLLimp(cur_cell),dit_xvLLimp(cur_cell));
end

%%

for nn = 1:length(tr_set)
    fprintf('Inferring fixation corrections, XV %d of %d\n',nn,length(tr_set));
    
    %% PREPROCESS MODEL COMPONENTS
    cur_xv_cells = tr_set(nn);
    cur_tr_cells = setdiff(tr_set,cur_xv_cells);
    cur_tr_cinds = find(ismember(tr_set,cur_tr_cells));
    cur_xv_cinds = find(ismember(tr_set,cur_xv_cells));
    cur_n_tr_chs = length(cur_tr_cells);
    
    filt_bank = zeros(n_squared_filts+1,klen_us);
    cur_Xtargs = [ref_mods(cur_xv_cells).mods(:).Xtarget];
    cur_k = [ref_mods(cur_xv_cells).mods(cur_Xtargs == 1).filtK];
    n_used_filts = size(cur_k,2);
    filt_bank(1:n_used_filts,:) = cur_k';
    mod_spkNL_params = ref_mods_spkNL(cur_xv_cells).spk_NL_params;
    lin_kerns = ref_mods(cur_xv_cells).mods(cur_Xtargs == 2).filtK;
    if use_sac_kerns
        sac_kerns = ref_mods(cur_xv_cells).mods(cur_Xtargs == 3).filtK;
        msac_kerns = ref_mods(cur_xv_cells).mods(cur_Xtargs == 4).filtK;
    end
    
    %indicator predictions
    block_out = Xblock(used_inds,:)*lin_kerns;
    if use_sac_kerns
        sac_out = Xsac*sac_kerns;
        msac_out = Xmsac*msac_kerns;
    end
    %% ESTIMATE LL for each shift in each stimulus frame
    
    %precompute LL at all shifts for all units
    xv_frame_LLs = nan(NT,n_shifts);
    for xx = 1:length(shifts)
        %         fprintf('Shift %d of %d\n',xx,n_shifts);
        %         cur_stim_shift = cur_Xmat*shift_mat{xx};
        cur_filt_shift = (filt_bank*filt_shift_mat{xx})';
        
        %outputs of stimulus models at current X-matrix shift
        gfuns = ones(NT,1)*mod_spkNL_params(1);
        gfuns = gfuns + cur_Xmat*cur_filt_shift(:,1);
        for ff = 2:(n_squared_filts+1)
            gfuns = gfuns + (cur_Xmat*cur_filt_shift(:,ff)).^2;
        end
        
        %add contributions from extra lin kernels
        gfuns = gfuns + block_out;
        if use_sac_kerns
            gfuns = gfuns + sac_out + msac_out;
        end
        
        %incorporate beta
        gfuns = gfuns*mod_spkNL_params(2);
        
        %handle numerical overflow with log(1+exp)
        too_large = gfuns > 50;
        pred_rate = log(1+exp(gfuns));
        pred_rate(too_large) = gfuns(too_large);
        
        %incorporate alpha
        pred_rate = pred_rate*mod_spkNL_params(3);
        
        %enforce min predicted rate
        pred_rate(pred_rate < 1e-50) = 1e-50;
        
        xv_frame_LLs(:,xx) = squeeze(nansum(Robs_mat(:,cur_xv_cinds).*log(pred_rate) - pred_rate,2));
    end
    
    %% INFER MICRO-SAC SEQUENCE
    cur_frame_LLs = frame_LLs - xv_frame_LLs;
    
    lalpha=zeros(NT,n_shifts);
    lbeta = zeros(NT,n_shifts);
    lscale=zeros(NT,1); %initialize rescaling parameters
    %compute rescaled forward messages
    lalpha(1,:) = leps_prior + cur_frame_LLs(1,:);
    lscale(1)=logsumexp(lalpha(1,:));
    lalpha(1,:) = lalpha(1,:) - lscale(1);
    for t=2:NT
        if use_prior(t)
            cur_lA = lA_tflip;
        else
            cur_lA = lA;
        end
        lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + cur_frame_LLs(t,:);
        lscale(t) = logsumexp(lalpha(t,:));
        lalpha(t,:) = lalpha(t,:) - lscale(t);
    end
    
    %compute rescaled backward messages
    lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
    for t=NT-1:-1:1
        if use_prior(t+1)
            cur_lA = lA_tflip;
        else
            cur_lA = lA;
        end
        lf1 = lbeta(t+1,:) + cur_frame_LLs(t+1,:);
        lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
    end
    lgamma= lalpha + lbeta;
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    
    gamma = exp(lgamma);
    post_mean(nn,:) = sum(bsxfun(@times,gamma,shifts),2);
    post_std(nn,:) = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - post_mean(nn,:)'.^2);
    
    %% RECOMPUTE XMAT
    drift_Xmat = reshape(all_Xmat_up_fixcor,[NT flen full_nPix_us]);
    for ii = 1:NT
        drift_Xmat(ii,:,:) = shift_matrix_Nd(drift_Xmat(ii,:,:),-round(drift_post_mean_cor(ii)),3);
    end
    drift_Xmat = reshape(drift_Xmat,[NT flen*full_nPix_us]);
    
    %% REFIT XV CELLS
    cur_X{1} = all_Xmat_up_fixcor(:,use_kInds_up);
    
    silent = 1;
    cur_cell = cur_xv_cells;
    fprintf('Refitting model \n');
    cur_unit_ind = find(tr_set == cur_cell);
    cur_tr_uset = tr_inds(~isnan(Robs_mat(tr_inds,cur_unit_ind)));
    cur_xv_uset = xv_inds(~isnan(Robs_mat(xv_inds,cur_unit_ind)));
    
    tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
    xv_X = get_Xcell_tInds(cur_X,cur_xv_uset);
    
    drift_xv_mods(cur_cell) = xv_mods(cur_cell);
    drift_xv_mods(cur_cell) = NMMfit_filters(drift_xv_mods(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),...
        tr_X,[],[],silent); %fit stimulus filters
    
    %refit spk NL
    drift_xv_mods_spkNL(cur_cell) = NMMfit_logexp_spkNL(drift_xv_mods(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
    
    newLL = NMMmodel_eval(drift_xv_mods_spkNL(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
    
    drift_xv_xvmods(cur_cell) = all_xvmod_fits(cur_cell);
    drift_xv_xvmods(cur_cell) = NMMfit_filters(drift_xv_xvmods(cur_cell),Robs_mat(cur_xv_uset,cur_unit_ind),...
        xv_X,[],[],silent); %fit stimulus filters
    %refit spk NL
    drift_xv_xvmods_spkNL(cur_cell) = NMMfit_logexp_spkNL(drift_xv_xvmods(cur_cell),Robs_mat(cur_xv_uset,cur_unit_ind),xv_X);
    
    newxvLL = NMMmodel_eval(drift_xv_xvmods_spkNL(cur_cell),Robs_mat(cur_xv_uset,cur_unit_ind),xv_X);
    
    drift_xv_LLimp(cur_cell) = (newLL - null_LL(cur_cell))/log(2);
    drift_xv_xvLLimp(cur_cell) = (newxvLL - null_xvLL(cur_cell))/log(2);
    
end

%%
cd(anal_dir);
save(anal_name,'drift_xv*','*LLimp','it_*','drift_post*','post_*','tr_set','tr_inds','jump_inds','xv_fix*','fix_post*');
%%
figure
plot(all_mod_LLimp(tr_set),it_LLimp(tr_set),'.','markersize',8)
hold on
plot(all_mod_LLimp(tr_set),xv_LLimp(tr_set),'r.','markersize',8);
xlim([0 0.5]);
ylim([0 0.5]);
legend('All units','Leave-one-out','Location','Southeast');
line([0 0.5],[0 0.5],'color','k');
b1 = regress(it_LLimp(tr_set)',all_mod_LLimp(tr_set)');
b2 = regress(xv_LLimp(tr_set)',all_mod_LLimp(tr_set)');
xx = linspace(0,0.5,100);
plot(xx,b1*xx);
plot(xx,b2*xx,'r')
xlabel('Initial LL (bits/spk)','fontsize',12);
ylabel('Final LL (bits/spk)','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
box off
fillPage(gcf,'papersize',[4 4]);
fname = 'fix_LOOxv_LLscatter';

%%
figure
plot(it_LLimp(tr_set),dit_LLimp(tr_set),'.','markersize',8)
hold on
plot(it_LLimp(tr_set),drift_xv_LLimp(tr_set),'r.','markersize',8)
legend('All units','Leave-one-out','Location','Southeast');
xlim([0 0.7]);
ylim([0 0.7]);
line([0 0.7],[0 0.7],'color','k');
b1 = regress(dit_LLimp(tr_set)',it_LLimp(tr_set)');
b2 = regress(drift_xv_LLimp(tr_set)',it_LLimp(tr_set)');
xx = linspace(0,0.7,100);
plot(xx,b1*xx);
plot(xx,b2*xx,'r')
xlabel('Initial LL (bits/spk)','fontsize',12);
ylabel('Final LL (bits/spk)','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
box off
fillPage(gcf,'papersize',[4 4]);
fname = 'drift_LOOxv_LLscatter';

%%
use_ids = 1:length(tr_set);
use_ids = find(dit_LLimp(tr_set) > 0.1);
fix_pop_mean = mean(xv_fix_post_mean(use_ids,:));
fix_pop_std = std(xv_fix_post_mean(use_ids,:));
fix_pop_mean = fix_pop_mean(fix_ids);
fix_pop_std = fix_pop_std(fix_ids);
full_pop_mean = mean(post_mean(use_ids,:));
full_pop_std = std(post_mean(use_ids,:));

close all
n_trials = length(unique(all_trialvec));
for tt = 50:n_trials
    % for tt = [96 137 154 179 376 409]
    uu = find(all_trialvec(used_inds) == tt);
    if ~isempty(uu)
        bt = all_t_axis(used_inds(uu(1)));
        et = all_t_axis(used_inds(uu(end)));
        dur = et-bt;
        if dur > 3.5
            hold on
            h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fix_pop_mean(uu)*sp_dx,fix_pop_std(uu)*sp_dx,{'color','b'});
            h2=shadedErrorBar(all_t_axis(used_inds(uu))-bt,full_pop_mean(uu)*sp_dx,full_pop_std(uu)*sp_dx,{'color','k'});
            if bar_ori == 0
                h3=plot(all_t_axis(used_inds(uu))-bt,corrected_interp_eyevals(used_inds(uu),2),'r','linewidth',2);
%                 h4=plot(all_t_axis(used_inds(uu))-bt,corrected_interp_eyevals(used_inds(uu),4)-median(corrected_interp_eyevals(used_inds(uu),4)),'color',[0.2 0.8 0.2],'linewidth',2);
            else
                h3=plot(all_t_axis(used_inds(uu))-bt,corrected_interp_eyevals(used_inds(uu),1),'r','linewidth',2);
%                 h4=plot(all_t_axis(used_inds(uu))-bt,corrected_interp_eyevals(used_inds(uu),3)-median(corrected_interp_eyevals(used_inds(uu),3)),'color',[0.2 0.8 0.2],'linewidth',2);
            end
            
%             plot(all_t_axis(used_inds(uu))-bt,nanmean(Robs_mat(uu,:),2)/5,'k');
            
            %             legend([h1.mainLine h2.mainLine h3 h4],{'Fixation corrections','Drift corrections','Left eye','Right eye'})
            xlim([0 dur]);
            ylim([-0.5 0.5]);
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
avg_drift = mean(post_mean(:,tr_inds));
var_drift = bsxfun(@minus,post_mean(:,tr_inds),avg_drift);

tr_fixs = find(ismember(jump_inds,tr_inds));
avg_fix = mean(xv_fix_post_mean(:,tr_fixs));
var_fix = bsxfun(@minus,xv_fix_post_mean(:,tr_fixs),avg_fix);