clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

Expt_num = 93;
Expt_name = sprintf('G0%d',Expt_num);
data_dir = ['~/Data/bruce/' Expt_name];
cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

anal_name = 'monoc_eyecorr_vbar_v4';
mod_data_name = 'monoc_eyecorr_vbar_mods_v4';
outdata_name = 'monoc_eyecorr_vbar_refitmods';

use_measured_pos = 0;
%%
stim_fs = 100; %in Hz
dt = 0.01;
flen = 12;
use_nPix = 20;
full_nPix = 36;
stim_params = NIMcreate_stim_params([flen full_nPix],dt);
Fr = 1;
bar_ori = 90;
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

%% LOAD IN EYE_TRACKING DATA
cd(anal_dir);
load(anal_name);

sac_shift = round(0.05/dt);

trial_start_inds = [1; find(diff(all_trialvec(used_inds)) ~= 0) + 1];
cjump_inds = sort([saccade_start_inds; trial_start_inds; NT]);
n_fixs = length(cjump_inds)-1;
fpost_mean_cor2 = fpost_mean;
for ii = 1:n_fixs
    cur_inds = cjump_inds(ii):(cjump_inds(ii+1)-1);
%     cur_inds(cur_inds > NT) = [];
    fpost_mean_cor2(cur_inds(1:end-sac_shift + 1)) = fpost_mean(cur_inds(sac_shift:end));
    cur_inds = cjump_inds(ii):(cjump_inds(ii)+sac_shift); cur_inds(cur_inds > NT) = [];
%     cur_inds(cur_inds > NT) = [];
    fpost_mean_cor2(cur_inds) = fpost_mean(cur_inds(end));
end

fix_corrections = it_all_fix_post_mean_cor;
full_corrections = it_all_fix_post_mean_cor + fpost_mean_cor2;


%% RECOMPUTE XMAT
all_shift_stimmat_up = all_stimmat_up;
for i=1:NT
    all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-round(full_corrections(i)),2);
end
all_Xmat_up_fixcor = create_time_embedding(all_shift_stimmat_up,stim_params_us);
all_Xmat_up_fixcor = all_Xmat_up_fixcor(used_inds,:);

all_Xmat_up = create_time_embedding(all_stimmat_up,stim_params_us);
all_Xmat_up = all_Xmat_up(used_inds,:);

%%
tr_set = 1:length(all_mod_SU);
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

%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITS
% fin_stim_params(1) = NMMcreate_stim_params([flen use_nPix_us],dt);
% fin_stim_params(2) = NMMcreate_stim_params([n_blocks 1],dt);
% 
% null_stim_params = fin_stim_params(2:end);
% 
% block_L2 = 1;
% silent = 0;
% 
% % base_lambda_d2XT = 6;
% % base_lambda_L1 = 6;
% base_lambda_d2XT = 60;
% base_lambda_L1 = 15;
% 
% init_optim_p.optTol = 1e-5; init_optim_p.progTol = 1e-8; init_optim_p.maxIter = 30;
% 
% null_reg_params = NMMcreate_reg_params('lambda_L2',[block_L2]);
% 
% uu = 105;
% target_su = tr_set(uu);
% 
% n_squared_filts = 4;
% mod_signs = [1 1 1 -1 -1 1];
% NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts) {'lin'}];
% init_d2XT = [ones(n_squared_filts+1,1); 0;];
% 
% init_L2 = [zeros(n_squared_filts+1,1); block_L2];
% init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_d2XT,'lambda_L2',init_L2,...
%     'mixing_prop',repmat([1 0.75],n_squared_filts+2,1),'boundary_conds',repmat([0 0 0],n_squared_filts+2,1));
% init_Xtargs = [ones(n_squared_filts+1,1); 2];
% 
% cur_uset = find(~isnan(Robs_mat(:,uu)));
% cur_X{1} = all_Xmat_up_fixcor(cur_uset,use_kInds_up);
% cur_X{2} = Xblock(used_inds(cur_uset),:);
% gqm1 = NMMinitialize_model(fin_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
% gqm1 = NMMfit_filters(gqm1,Robs_mat(cur_uset,uu),cur_X,[],[],silent,init_optim_p);
% [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,Robs_mat(cur_uset,uu),cur_X);
% gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
% gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
% cor_gqm2 = NMMfit_filters(gqm1,Robs_mat(cur_uset,uu),cur_X,[],[],silent);
% [LL, penLL, pred_rate, G, cor_gint] = NMMmodel_eval(cor_gqm2,Robs_mat(cur_uset,uu),cur_X);
% 
% % NO COR
% cur_X{1} = all_Xmat_up(cur_uset,use_kInds_up);
% cur_X{2} = Xblock(used_inds(cur_uset),:);
% gqm2 = NMMfit_filters(cor_gqm2,Robs_mat(cur_uset,uu),cur_X,[],[],silent);
% [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm2,Robs_mat(cur_uset,uu),cur_X);
% 
% %% GENERATE FULL MODEL COMPARISON FIGURE
% lag_axis = (0:flen-1)*dt;
% pix_axis = (-floor(use_nPix_us/2):floor(use_nPix_us/2))*.05/spatial_usfac;
% 
% cor_filts = [cor_gqm2.mods(find(init_Xtargs == 1)).filtK];
% uncor_filts = [gqm2.mods(find(init_Xtargs == 1)).filtK];
% max_vals = max(abs(cor_filts));
% uc_max_vals = max(abs(uncor_filts));
% 
% cor_filts = reshape(cor_filts,[flen use_nPix_us n_squared_filts+1]);
% uncor_filts = reshape(uncor_filts,[flen use_nPix_us n_squared_filts+1]);
% 
% figure
% n_stimfilts = n_squared_filts + 1;
% for ii = 1:n_stimfilts
% subplot(n_stimfilts,3,(ii-1)*3+1);
% imagesc(pix_axis,lag_axis,squeeze(cor_filts(:,:,ii))); caxis([-max_vals(ii) max_vals(ii)]);
% xlabel('Position (deg)','fontsize',6);
% ylabel('Lag (s)','fontsize',6);
% set(gca,'xtick',-0.5:0.25:0.5,'fontsize',5,'fontname','arial');
% subplot(n_stimfilts,3,(ii-1)*3+2);
% imagesc(pix_axis,lag_axis,squeeze(uncor_filts(:,:,ii))); caxis([-max_vals(ii) max_vals(ii)]);
% xlabel('Position (deg)','fontsize',6);
% ylabel('Lag (s)','fontsize',6);
% set(gca,'xtick',-0.5:0.25:0.5,'fontsize',5,'fontname','arial');
% subplot(n_stimfilts,3,(ii-1)*3+3);
% [h_y,h_x] = hist(cor_gint(:,ii),100);
% h_yunc = hist(gint(:,ii),h_x);
% if ii == 1
%     cur_mody = h_x;
% else
%     cur_mody = h_x.^2;
% end
% cur_mody = cur_mody-min(cur_mody);
% cur_mody = cur_mody/max(cur_mody)*max(h_yunc);
% plot(h_x,h_y,h_x,h_yunc,'r','linewidth',2)
% hold on
% plot(h_x,cur_mody,'k')
% xlim([h_x(1) -h_x(1)])
% xlabel('Generating signal','fontsize',6);
% set(gca,'xtick',0,'ytick',0,'fontsize',5,'fontname','arial');
% end
% colormap(gray);
% subplot(n_stimfilts,3,1);
% title('Corrected','fontsize',8);
% subplot(n_stimfilts,3,2);
% title('Uncorrected','fontsize',8);
% fillPage(gcf,'papersize',[5 9]);
% 

%%
% max_vals = max(max(abs(cor_filts)));
% uc_max_vals = max(max(abs(uncor_filts)));
% figure
% subplot(2,1,1)
% imagesc(pix_axis,lag_axis,squeeze(cor_filts(:,:,1))); caxis([-max_vals(1) max_vals(1)]);
% xlabel('Position (deg)','fontsize',8);
% ylabel('Lag (s)','fontsize',8);
% set(gca,'xtick',-0.5:0.25:0.5,'fontsize',6,'fontname','arial');
% title('Corrected','fontsize',10);
% set(gca,'ydir','normal');
% subplot(2,1,2)
% imagesc(pix_axis,lag_axis,squeeze(uncor_filts(:,:,1))); caxis([-uc_max_vals(1) uc_max_vals(1)]);
% xlabel('Position (deg)','fontsize',8);
% ylabel('Lag (s)','fontsize',8);
% title('Uncorrected','fontsize',10);
% set(gca,'xtick',-0.5:0.25:0.5,'fontsize',6,'fontname','arial');
% fillPage(gcf,'papersize',[3 5]);
% colormap(gray);
% set(gca,'ydir','normal');
% 
% use_filt = 2;
% figure
% subplot(2,1,1)
% imagesc(pix_axis,lag_axis,squeeze(cor_filts(:,:,use_filt))); caxis([-max_vals(use_filt) max_vals(use_filt)]);
% xlabel('Position (deg)','fontsize',8);
% ylabel('Lag (s)','fontsize',8);
% set(gca,'xtick',-0.5:0.25:0.5,'fontsize',6,'fontname','arial');
% title('Corrected','fontsize',10);
% set(gca,'ydir','normal');
% subplot(2,1,2)
% imagesc(pix_axis,lag_axis,squeeze(uncor_filts(:,:,use_filt))); caxis([-uc_max_vals(use_filt) uc_max_vals(use_filt)]);
% xlabel('Position (deg)','fontsize',8);
% ylabel('Lag (s)','fontsize',8);
% title('Uncorrected','fontsize',10);
% set(gca,'xtick',-0.5:0.25:0.5,'fontsize',6,'fontname','arial');
% fillPage(gcf,'papersize',[3 5]);
% colormap(gray);
% set(gca,'ydir','normal');
% 
% figure
% subplot(2,1,1)
% imagesc(Fx,Ft,squeeze(cor_ffts(:,:,1)));
% xlabel('Spatial frequency (cyc/deg)','fontsize',10);
% ylabel('Temporal frequency (Hz)','fontsize',10);
% set(gca,'fontsize',8,'fontname','arial');
% title('Corrected','fontsize',10);
% xlim([-12 12]);
% ylim([-40 40]);
% line([-12 12],[0 0],'color','w');
% line([0 0],[-40 40],'color','w');
% set(gca,'ydir','normal');
% subplot(2,1,2)
% imagesc(Fx,Ft,squeeze(uncor_ffts(:,:,1)));
% xlabel('Spatial frequency (cyc/deg)','fontsize',10);
% ylabel('Temporal frequency (Hz)','fontsize',10);
% xlim([-12 12]);
% ylim([-40 40]);
% set(gca,'fontsize',8,'fontname','arial');
% fillPage(gcf,'papersize',[3 5]);
% colormap(gray);
% line([-12 12],[0 0],'color','w');
% line([0 0],[-40 40],'color','w');
% set(gca,'ydir','normal');
% 
% figure
% subplot(2,1,1)
% imagesc(Fx,Ft,squeeze(cor_ffts(:,:,use_filt)));
% xlabel('Spatial frequency (cyc/deg)','fontsize',10);
% ylabel('Temporal frequency (Hz)','fontsize',10);
% set(gca,'fontsize',8,'fontname','arial');
% title('Corrected','fontsize',10);
% set(gca,'ydir','normal');
% xlim([-12 12]);
% ylim([-40 40]);
% line([-12 12],[0 0],'color','w');
% line([0 0],[-40 40],'color','w');
% subplot(2,1,2)
% imagesc(Fx,Ft,squeeze(uncor_ffts(:,:,use_filt)));
% xlabel('Spatial frequency (cyc/deg)','fontsize',10);
% ylabel('Temporal frequency (Hz)','fontsize',10);
% xlim([-12 12]);
% ylim([-40 40]);
% set(gca,'fontsize',8,'fontname','arial');
% fillPage(gcf,'papersize',[3 5]);
% colormap(gray);
% line([-12 12],[0 0],'color','w');
% line([0 0],[-40 40],'color','w');
% set(gca,'ydir','normal');

%% REFIT ALL MODELS BEFORE AND AFTER
fin_stim_params(1) = NMMcreate_stim_params([flen use_nPix_us],dt);
fin_stim_params(2) = NMMcreate_stim_params([n_blocks 1],dt);

null_stim_params = fin_stim_params(2:end);

block_L2 = 1;
silent = 1;

% base_lambda_d2XT = 6;
% base_lambda_L1 = 6;
base_lambda_d2XT = 60;
base_lambda_L1 = 15;

init_optim_p.optTol = 1e-5; init_optim_p.progTol = 1e-8; init_optim_p.maxIter = 30;

null_reg_params = NMMcreate_reg_params('lambda_L2',[block_L2]);

cant_use = false(length(tr_set),1);
for uu = 1:length(tr_set)
    fprintf('Fitting corrected model unit %d of %d\n',uu,length(tr_set));
    target_su = tr_set(uu);
    
    n_squared_filts = 2;
    mod_signs = [1 1 1 1];
    NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts) {'lin'}];
    init_d2XT = [ones(n_squared_filts+1,1); 0;];
    
    init_L2 = [zeros(n_squared_filts+1,1); block_L2];
    init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_d2XT,'lambda_L2',init_L2,...
        'mixing_prop',repmat([1 0.75],n_squared_filts+2,1),'boundary_conds',repmat([0 0 0],n_squared_filts+2,1));
    init_Xtargs = [ones(n_squared_filts+1,1); 2];
    
    cur_uset = find(~isnan(Robs_mat(:,uu)));
    data_dur(uu) = length(cur_uset)*dt;
    if data_dur(uu) > 60*5
        cur_X{1} = all_Xmat_up_fixcor(cur_uset,use_kInds_up);
        cur_X{2} = Xblock(used_inds(cur_uset),:);
        gqm1 = NMMinitialize_model(fin_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
        gqm1 = NMMfit_filters(gqm1,Robs_mat(cur_uset,uu),cur_X,[],[],silent,init_optim_p);
        [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,Robs_mat(cur_uset,uu),cur_X);
        gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
        gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
        cor_mod(uu) = NMMfit_filters(gqm1,Robs_mat(cur_uset,uu),cur_X,[],[],silent);
        
        fprintf('Fitting UNcorrected model unit %d of %d\n',uu,length(tr_set));
        % NO COR
        cur_X{1} = all_Xmat_up(cur_uset,use_kInds_up);
        cur_X{2} = Xblock(used_inds(cur_uset),:);
        uncor_mod(uu) = NMMfit_filters(cor_mod(uu),Robs_mat(cur_uset,uu),cur_X,[],[],silent);
    else
        cant_use(uu) = 1;
    end
end
%% COMPUTE STATS ACROSS ALL FILTERS

zpad_factor = 5;
dx =sp_dx;
Ft = linspace(-1/dt/2,1/dt/2,zpad_factor*flen);
Fx = linspace(-1/dx/2,1/dx/2,zpad_factor*use_nPix_us);

sig_Ft = 0.1; sig_Fx = 0.1;
[FFx,FFt] = meshgrid(Fx,Ft);
gauss_kern = FFt.^2/(2*sig_Ft^2) + FFx.^2/(2*sig_Fx^2);
gauss_kern = exp(-gauss_kern);
gauss_kern = gauss_kern/sum(gauss_kern(:));


for cc = 1:length(tr_set)
    if ~isempty(cor_mod(cc).LL_seq)
    cc
    cur_cor_mod = cor_mod(cc);
    cur_uncor_mod = uncor_mod(cc);
    cor_Xtargs = [cur_cor_mod.mods(:).Xtarget];
    uncor_Xtargs = [cur_uncor_mod.mods(:).Xtarget];
    
    cor_filts = reshape([cur_cor_mod.mods(find(cor_Xtargs == 1)).filtK],[flen use_nPix_us n_squared_filts+1]);
    uncor_filts = reshape([cur_uncor_mod.mods(find(uncor_Xtargs == 1)).filtK],[flen use_nPix_us n_squared_filts+1]);
    
    cor_spatial_profiles = squeeze(std(cor_filts));
    cor_meanlocs(cc,:) = sum(bsxfun(@times,cor_spatial_profiles,(1:use_nPix_us)'))./sum(cor_spatial_profiles);
    cor_diff = bsxfun(@minus,repmat(cor_meanlocs(cc,:),use_nPix_us,1),(1:use_nPix_us)');
    cor_std(cc,:) = sqrt(sum(cor_spatial_profiles.*cor_diff.^2)./sum(cor_spatial_profiles));
    lin_spatial_profiles = cor_spatial_profiles(:,1);
    if max(lin_spatial_profiles) > 0
        [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',lin_spatial_profiles);
        cor_lin_mean(cc) = fit_params(1);
        cor_lin_std(cc) = fit_params(2);
    else
        cor_lin_mean(cc) = nan;
        cor_lin_std(cc) = nan;
    end
    avg_spatial_profiles = mean(cor_spatial_profiles(:,2:3),2);
    if max(avg_spatial_profiles) > 0
        [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',avg_spatial_profiles);
        cor_sq_mean(cc) = fit_params(1);
        cor_sq_std(cc) = fit_params(2);
    else
        cor_sq_mean(cc) = nan;
        cor_sq_std(cc) = nan;
    end
    uncor_spatial_profiles = squeeze(std(uncor_filts));
    uncor_meanlocs(cc,:) = sum(bsxfun(@times,uncor_spatial_profiles,(1:use_nPix_us)'))./sum(uncor_spatial_profiles);
    uncor_diff = bsxfun(@minus,repmat(uncor_meanlocs(cc,:),use_nPix_us,1),(1:use_nPix_us)');
    uncor_std(cc,:) = sqrt(sum(uncor_spatial_profiles.*cor_diff.^2)./sum(uncor_spatial_profiles));
    lin_spatial_profiles = uncor_spatial_profiles(:,1);
    if max(lin_spatial_profiles) > 0
        [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',lin_spatial_profiles);
        uncor_lin_mean(cc) = fit_params(1);
        uncor_lin_std(cc) = fit_params(2);
    else
        uncor_lin_mean(cc) = nan;
        uncor_lin_std(cc) = nan;
    end
    avg_spatial_profiles = mean(uncor_spatial_profiles(:,2:3),2);
    if max(avg_spatial_profiles) > 0
        [fit_params,fit_z] = fitGaussianCurve((1:use_nPix_us)',avg_spatial_profiles);
        uncor_sq_mean(cc) = fit_params(1);
        uncor_sq_std(cc) = fit_params(2);
    else
        uncor_sq_mean(cc) = nan;
        uncor_sq_std(cc) = nan;
    end
%     cor_temporal_profiles = squeeze(std(cor_filts,[],2));
%     uncor_temporal_profiles = squeeze(std(uncor_filts,[],2));
%     [~,cor_peak_lags] = max(cor_temporal_profiles);
%     [~,uncor_peak_lags] = max(uncor_temporal_profiles);
    
    cor_max_pow(cc,:) = zeros(1,n_squared_filts+1);
    cor_max_ploc(cc,:) = zeros(1,n_squared_filts+1);
    cor_tot_pow(cc,:) = zeros(1,n_squared_filts+1);
    cor_filts_zpad = zeros(zpad_factor*flen,zpad_factor*use_nPix_us,n_squared_filts+1);
    cor_filts_zpad((flen+1):2*flen,(use_nPix_us+1):2*use_nPix_us,:) = cor_filts;
    cor_ffts = zeros(size(cor_filts_zpad));
    for ii = 1:n_squared_filts+1
        cur_ffts = abs(fftshift(fft2(cor_filts_zpad(:,:,ii))));
        cur_ffts = conv2(squeeze(cur_ffts),gauss_kern,'same');
        cor_max_pow(cc,ii) = max(cur_ffts(:));
        cor_max_ploc(cc,ii) = find(cur_ffts == cor_max_pow(cc,ii),1);
        cor_tot_pow(cc,ii) = sum(cur_ffts(:));
        cor_ffts(:,:,ii) = cur_ffts;
    end
    cor_FFx(cc,:) = abs(FFx(cor_max_ploc(cc,:)));
    cor_FFt(cc,:) = abs(FFt(cor_max_ploc(cc,:)));
    
    uncor_max_pow(cc,:) = zeros(n_squared_filts+1,1);
    uncor_max_ploc(cc,:) = zeros(n_squared_filts+1,1);
    uncor_tot_pow(cc,:) = zeros(n_squared_filts+1,1);
    uncor_filts_zpad = zeros(zpad_factor*flen,zpad_factor*use_nPix_us,n_squared_filts+1);
    uncor_filts_zpad((flen+1):2*flen,(use_nPix_us+1):2*use_nPix_us,:) = uncor_filts;
    uncor_ffts = zeros(size(uncor_filts_zpad));
    for ii = 1:n_squared_filts+1
        cur_ffts = abs(fftshift(fft2(uncor_filts_zpad(:,:,ii))));
        cur_ffts = conv2(squeeze(cur_ffts),gauss_kern,'same');
        uncor_max_pow(cc,ii) = max(cur_ffts(:));
        uncor_max_ploc(cc,ii) = find(cur_ffts == uncor_max_pow(cc,ii),1);
        uncor_tot_pow(cc,ii) = sum(cur_ffts(:));
        uncor_ffts(:,:,ii) = cur_ffts;
    end
    uncor_FFx(cc,:) = abs(FFx(uncor_max_ploc(cc,:)));
    uncor_FFt(cc,:) = abs(FFt(uncor_max_ploc(cc,:)));
    else
        cor_FFx(cc,:) = nan;
        uncor_FFx(cc,:) = nan;
        cor_tot_pow(cc,:) = nan;
        uncor_tot_pow(cc,:) = nan;
    end
end
%%
cd(anal_dir)
save(outdata_name,'cor*','uncor*','flen','use_nPix','tr_set')
%%
is_su = find(all_mod_SU(tr_set) > 0);

cd(anal_dir);
xname = 'monoc_eyecorr_vbar_refitmods';
yname = 'monoc_eyecorr_hbar_refitmods';
load(xname);
clear x_cor_stats x_uncor_stats
for i = 1:length(tr_set)
    x_cor_stats(i).lin_mean = cor_lin_mean(i);
    x_cor_stats(i).lin_std = cor_lin_std(i);
    x_cor_stats(i).sq_mean = cor_sq_mean(i);
    x_cor_stats(i).sq_std = cor_sq_std(i);
    x_cor_stats(i).lin_pow = cor_tot_pow(i,1);
    x_cor_stats(i).sq_pow = mean(cor_tot_pow(i,2:end),2);
    x_cor_stats(i).lin_FFx = cor_FFx(i,1);
    x_cor_stats(i).sq_FFx = mean(cor_FFx(i,2:end),2);

    x_uncor_stats(i).lin_mean = uncor_lin_mean(i);
    x_uncor_stats(i).lin_std = uncor_lin_std(i);
    x_uncor_stats(i).sq_mean = uncor_sq_mean(i);
    x_uncor_stats(i).sq_std = uncor_sq_std(i);
    x_uncor_stats(i).lin_pow = uncor_tot_pow(i,1);
    x_uncor_stats(i).sq_pow = mean(uncor_tot_pow(i,2:end),2);
    x_uncor_stats(i).lin_FFx = uncor_FFx(i,1);
    x_uncor_stats(i).sq_FFx = mean(uncor_FFx(i,2:end),2);
    
    if ~isempty(cor_mod(i).LL_seq)
        x_LL(i) = cor_mod(i).LL_seq(end);
    else
        x_LL(i) = -Inf;
    end
end


load(yname);
clear y_cor_stats y_uncor_stats
for i = 1:length(tr_set)
    y_cor_stats(i).lin_mean = cor_lin_mean(i);
    y_cor_stats(i).lin_std = cor_lin_std(i);
    y_cor_stats(i).sq_mean = cor_sq_mean(i);
    y_cor_stats(i).sq_std = cor_sq_std(i);
    y_cor_stats(i).lin_pow = cor_tot_pow(i,1);
    y_cor_stats(i).sq_pow = mean(cor_tot_pow(i,2:end),2);
    y_cor_stats(i).lin_FFx = cor_FFx(i,1);
    y_cor_stats(i).sq_FFx = mean(cor_FFx(i,2:end),2);

    y_uncor_stats(i).lin_mean = uncor_lin_mean(i);
    y_uncor_stats(i).lin_std = uncor_lin_std(i);
    y_uncor_stats(i).sq_mean = uncor_sq_mean(i);
    y_uncor_stats(i).sq_std = uncor_sq_std(i);
    y_uncor_stats(i).lin_pow = uncor_tot_pow(i,1);
    y_uncor_stats(i).sq_pow = mean(uncor_tot_pow(i,2:end),2);
    y_uncor_stats(i).lin_FFx = uncor_FFx(i,1);
    y_uncor_stats(i).sq_FFx = mean(uncor_FFx(i,2:end),2);
    
    if ~isempty(cor_mod(i).LL_seq)
        y_LL(i) = cor_mod(i).LL_seq(end);
    else
        y_LL(i) = -Inf;
    end
end

cor_stats = x_cor_stats;
cor_stats(y_LL > x_LL) = y_cor_stats(y_LL > x_LL);
uncor_stats = x_uncor_stats;
uncor_stats(y_LL > x_LL) = y_uncor_stats(y_LL > x_LL);

%replace MU with SU on electrodes that have it
su_inds = all_mod_SU(is_su);
x_cor_stats(su_inds) = x_cor_stats(97:end);
x_uncor_stats(su_inds) = x_uncor_stats(97:end);
y_cor_stats(su_inds) = y_cor_stats(97:end);
y_uncor_stats(su_inds) = y_uncor_stats(97:end);
cor_stats(su_inds) = cor_stats(97:end);
uncor_stats(su_inds) = uncor_stats(97:end);

x_cor_stats(97:end) = [];
x_uncor_stats(97:end) = [];
y_cor_stats(97:end) = [];
y_uncor_stats(97:end) = [];
cor_stats(97:end) = [];
uncor_stats(97:end) = [];


load ~/Data/bruce/general_array_data/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;

x_sq_stds = [x_uncor_stats.sq_std]*dx;
y_sq_stds = [y_uncor_stats.sq_std]*dx;
x_sq_means = ([x_uncor_stats.sq_mean]-use_nPix_us/2)*dx + 0.34;
y_sq_means = ([y_uncor_stats.sq_mean]-use_nPix_us/2)*dx - 0.43;

for i = 1:96
    tempx(Y_pos(i),X_pos(i)) = x_sq_means(i);
    tempy(Y_pos(i),X_pos(i)) = y_sq_means(i);
    tempinds(Y_pos(i),X_pos(i)) = i;
end
weights = ones(10,10);
weights(tempx==0) = 0;
xpos_interp = smoothn(tempx,weights,'robust');
ypos_interp = smoothn(tempy,weights,'robust');

figure;
subplot(2,1,1); hold on
for ii = 1:96
    ellipse(2*x_sq_stds(ii),2*y_sq_stds(ii),0,x_sq_means(ii),y_sq_means(ii),'k','linewidth',0.5);
end
for ii = 1:length(su_inds)
    ellipse(2*x_sq_stds(su_inds(ii)),2*y_sq_stds(su_inds(ii)),0,x_sq_means(su_inds(ii)),y_sq_means(su_inds(ii)),'r','linewidth',2);
end
xlim([0 0.9]);ylim([-0.9 0]);
xlabel('Horizontal position (deg)','fontsize',10);
ylabel('Vertical position (deg)','fontsize',10);
legend('MU','SU');
set(gca,'fontsize',8,'fontname','arial');

subplot(2,1,2); hold on
plot(xpos_interp,ypos_interp,'ko','markersize',4)
plot(xpos_interp(su_inds),ypos_interp(su_inds),'ro','markersize',4);
xlim([0 0.9]);ylim([-0.9 0]);
xlabel('Horizontal position (deg)','fontsize',10);
ylabel('Vertical position (deg)','fontsize',10);
legend('MU','SU');
set(gca,'fontsize',8,'fontname','arial');
fillPage(gcf,'papersize',[4 8]);

uncor_SFs = [uncor_stats(:).sq_FFx];
uncor_stds = [uncor_stats(:).sq_std];
figure
subplot(2,1,1)
hist(uncor_SFs,100);
xlim([0 10]);
xlabel('Preferred SF (cyc/deg)','fontsize',10);
ylabel('Number of units','fontsize',10);
set(gca,'fontsize',8,'fontname','arial');
box off
subplot(2,1,2)
hist(uncor_stds*dx*2,100);
xlim([0.07 0.17]);
xlabel('RF width (deg)','fontsize',10);
ylabel('Number of units','fontsize',10);
set(gca,'fontsize',8,'fontname','arial');
box off
fillPage(gcf,'papersize',[4 7]);

% figure
% plot(x_squared_SFs,y_squared_SFs,'ko','markersize',4);
% hold on
% plot(x_squared_SFs(su_inds),y_squared_SFs(su_inds),'ro','linewidth',2,'markersize',4);
% xlim([0 9]); ylim([0 9]);
% xlabel('Preferred SF_x (cyc/deg)','fontsize',10);
% ylabel('Preferred SF_y (cyc/deg)','fontsize',10);
% set(gca,'fontsize',8,'fontname','arial');
% fillPage(gcf,'papersize',[4 4]);
% 
% figure
% plot(2*x_squared_stds,2*y_squared_stds,'ko','markersize',4);
% hold on
% plot(2*x_squared_stds(su_inds),2*y_squared_stds(su_inds),'ro','linewidth',2,'markersize',4);
% xlim([0.08 0.18]); ylim([0.08 0.18]);
% xlabel('RF x-width (deg)','fontsize',10);
% ylabel('RF y-width (deg)','fontsize',10);
% set(gca,'fontsize',8,'fontname','arial');
% fillPage(gcf,'papersize',[4 4]);

%%
FFx_jitter = 0.05;

cor_lin_tot_pow = [cor_stats(:).lin_pow];
uncor_lin_tot_pow = [uncor_stats(:).lin_pow];
used_linfilts = find(lin_tot_pow > 500);
used_linfilts_su = su_inds(ismember(su_inds,used_linfilts));
used_linfilts(ismember(used_linfilts,used_linfilts_su)) = [];

figure
subplot(2,3,1)
plot(sqrt(uncor_lin_tot_pow(used_linfilts)),sqrt(cor_lin_tot_pow(used_linfilts)),'k.')
hold on
plot(sqrt(uncor_lin_tot_pow(used_linfilts_su)),sqrt(cor_lin_tot_pow(used_linfilts_su)),'r.')
xlim([10 90]);ylim([10 90]);
line([10 90],[10 90],'color','b')
xlabel('Uncorrected power','fontsize',10);
ylabel('Corrected power','fontsize',10);
set(gca,'fontname','arial','fontsize',8,'xtick',0,'ytick',0)
box off;
title('Linear filters','fontsize',10);
% fillPage(gcf,'papersize',[4 4]);
% print('Pre_post_linpow','-dpdf');
% close

cor_lin_FFx = [cor_stats(:).lin_FFx];
uncor_lin_FFx = [uncor_stats(:).lin_FFx];
subplot(2,3,2)
plot(uncor_lin_FFx(used_linfilts)+randn(1,length(used_linfilts))*FFx_jitter,...
    cor_lin_FFx(used_linfilts)+randn(1,length(used_linfilts))*FFx_jitter,'k.')
hold on
plot(uncor_lin_FFx(used_linfilts_su)+randn(1,length(used_linfilts_su))*FFx_jitter,...
    cor_lin_FFx(used_linfilts_su)+randn(1,length(used_linfilts_su))*FFx_jitter,'r.')
xlim([0 5]);ylim([0 5]);
line([0 5],[0 5],'color','b')
xlabel('Uncorrected SF (cyc/deg)','fontsize',10);
ylabel('Corrected SF (cyc/deg)','fontsize',10);
% set(gca,'fontname','arial','fontsize',8,'xtick',0,'ytick',0)
set(gca,'fontname','arial','fontsize',8)
box off;
title('Linear filters','fontsize',10);
% fillPage(gcf,'papersize',[4 4]);
% print('Pre_post_linFFx','-dpdf');
% close

cor_lin_std = [cor_stats(:).lin_std];
uncor_lin_std = [uncor_stats(:).lin_std];
subplot(2,3,3)
plot(2*uncor_lin_std(used_linfilts)*dx,2*cor_lin_std(used_linfilts)*dx,'k.')
hold on
plot(2*uncor_lin_std(used_linfilts_su)*dx,2*cor_lin_std(used_linfilts_su)*dx,'r.')
xlim([0 0.4]);ylim([0 0.4]);
line([0 0.5],[0 0.5],'color','b')
xlabel('Uncorrected RF width (deg)','fontsize',10);
ylabel('Corrected RF width (deg)','fontsize',10);
% set(gca,'fontname','arial','fontsize',8,'xtick',0,'ytick',0)
set(gca,'fontname','arial','fontsize',8)
box off;
title('Linear filters','fontsize',10);
% fillPage(gcf,'papersize',[4 4]);
% print('Pre_post_linstd','-dpdf');
% close

used_sq_inds = setdiff(1:96,su_inds);
cor_sq_tot_pow = [cor_stats(:).sq_pow];
uncor_sq_tot_pow = [uncor_stats(:).sq_pow];
subplot(2,3,4)
plot(sqrt(uncor_sq_tot_pow(used_sq_inds)),sqrt(cor_sq_tot_pow(used_sq_inds)),'k.')
hold on
plot(sqrt(uncor_sq_tot_pow(su_inds)),sqrt(cor_sq_tot_pow(su_inds)),'r.')
xlim([30 70]);ylim([30 70]);
line([0 100],[0 100],'color','b')
xlabel('Uncorrected power','fontsize',10);
ylabel('Corrected power','fontsize',10);
set(gca,'fontname','arial','fontsize',8,'xtick',0,'ytick',0)
box off;
title('Squared filters','fontsize',10);
% fillPage(gcf,'papersize',[4 4]);
% print('Pre_post_quadpow','-dpdf');
% close

cor_sq_FFx = [cor_stats(:).sq_FFx];
uncor_sq_FFx = [uncor_stats(:).sq_FFx];
subplot(2,3,5)
plot(uncor_sq_FFx(used_sq_inds)+randn(1,length(used_sq_inds))*FFx_jitter,...
    cor_sq_FFx(used_sq_inds) +randn(1,length(used_sq_inds))*FFx_jitter,'k.')
hold on
plot(uncor_sq_FFx(su_inds) +randn(1,length(su_inds))*FFx_jitter,...
    cor_sq_FFx(su_inds) +randn(1,length(su_inds))*FFx_jitter,'r.')
xlim([0 10]);ylim([0 10]);
line([0 10],[0 10],'color','b')
xlabel('Uncorrected SF (cyc/deg)','fontsize',10);
ylabel('Corrected SF (cyc/deg)','fontsize',10);
% set(gca,'fontname','arial','fontsize',8,'xtick',0,'ytick',0)
set(gca,'fontname','arial','fontsize',8)
title('Squared filters','fontsize',10);
box off;
% fillPage(gcf,'papersize',[4 4]);
% print('Pre_post_quadFFx','-dpdf');
% close

cor_sq_std = [cor_stats(:).sq_std];
uncor_sq_std = [uncor_stats(:).sq_std];
subplot(2,3,6)
plot(2*uncor_sq_std(used_sq_inds)*dx,2*cor_sq_std(used_sq_inds)*dx,'k.')
hold on
plot(2*uncor_sq_std(su_inds)*dx,2*cor_sq_std(su_inds)*dx,'r.')
xlim([0 0.2]);ylim([0 0.2]);
line([0 20]*dx,[0 20]*dx,'color','b')
xlabel('Uncorrected RF width (deg)','fontsize',10);
ylabel('Corrected RF width (deg)','fontsize',10);
% set(gca,'fontname','arial','fontsize',8,'xtick',0,'ytick',0)
set(gca,'fontname','arial','fontsize',8)
title('Squared filters','fontsize',10);
box off;
% fillPage(gcf,'papersize',[4 4]);
% print('Pre_post_quadstd','-dpdf');
% close

fillPage(gcf,'papersize',[9 6]);

%%
su_set = find(all_mod_SU > 0);

for ii = su_set'
    if ~isempty(cor_mod(ii).LL_seq)
    NMMdisplay_model(cor_mod(ii),[],[],1);
    fprintf('FFx: %.3f\n', cor_FFx(ii,1));
    pause
    close all
    end
end