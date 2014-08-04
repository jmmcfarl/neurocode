clear all
% close all

dir_prefix = '/home/james/';
% dir_prefix = '/Volumes/james/';
fig_dir = [dir_prefix 'Analysis/bruce/ET_final/'];

Expt_num = 88;
Expt_name = sprintf('G0%d',Expt_num);
data_dir = [dir_prefix 'Data/bruce/' Expt_name];
cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = [dir_prefix 'Analysis/bruce/' Expt_name '/ET_final'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

bar_ori = 0;

if bar_ori == 90
    left_anal_name = 'left_eyecorr_vbar';
    right_anal_name = 'right_eyecorr_vbar';
else
    left_anal_name = 'left_eyecorr_hbar3';
    right_anal_name = 'right_eyecorr_hbar3';
end

recompute_init_mods = 0;
use_measured_pos = 0;
use_sac_kerns = 1;
use_LOOXV = 0;
use_coils = [0 0]; %[L R]

if any(use_coils > 0)
   anal_name = [anal_name '_Cprior'];
end
if use_measured_pos == 1
    mod_data_name = [mod_data_name '_Cinit'];
    anal_name = [anal_name '_Cinit'];
end

%dont fit stim models using these blocks
if Expt_num == 86
    ignore_blocks = [16 17 28 30]; %G086
elseif Expt_num == 87
    ignore_blocks = [15];
elseif Expt_num == 93
    ignore_blocks = [28];
else
   ignore_blocks = []; 
end

%%
xv_frac = 0.2;

flen = 12;
use_nPix = 16;

n_fix_inf_it = 3; %3
n_drift_inf_it = 1; %3

fix_prior_sigma = 0.15;
fix_noise_sigma = 0.1;
drift_noise_sigma = 0.003;
drift_prior_sigma = 0.005;
drift_jump_sigma = 0.05; %0.05 start
drift_dsf = 3;

min_trial_dur = 0.75;

spatial_usfac = 2;

%%
eps = -1e3;

stim_fs = 100; %in Hz
dt = 0.01;
full_nPix = 36;
stim_params = NIMcreate_stim_params([flen full_nPix],dt);
Fr = 1;

beg_buffer = 0.2;
end_buffer = 0.05;
trial_dur = 4;

n_probes = 96;

use_right_eye = false;

n_use_blocks = Inf;

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
temp_dx = [-2*max_shift:dshift:2*max_shift];
shift_dx_edges = [temp_dx*sp_dx-dshift*sp_dx/2 temp_dx(end)*sp_dx+dshift*sp_dx/2];
ovshift_dx_edges = [shifts*sp_dx-dshift*sp_dx/2 shifts(end)*sp_dx+dshift*sp_dx/2];

max_Dshift = round(8*spatial_usfac);
Dshifts = -max_Dshift:dshift:max_Dshift;
n_Dshifts = length(Dshifts);
zero_Dframe = find(Dshifts==0);
temp_dx = [-2*max_Dshift:dshift:2*max_Dshift];
Dshift_dx_edges = [temp_dx*sp_dx-dshift*sp_dx/2 temp_dx(end)*sp_dx+dshift*sp_dx/2];
ovDshift_dx_edges = [Dshifts*sp_dx-dshift*sp_dx/2 Dshifts(end)*sp_dx+dshift*sp_dx/2];

max_Tshift = max_shift + max_Dshift;
Tshifts = -max_Tshift:dshift:max_Tshift;

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
% cur_block_set = find(included_type & expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori);
cur_block_set = find(included_type & expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori & expt_dds == 12);

cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];

n_blocks = length(cur_block_set);

sim_sac_expts = find(~expt_has_ds(cur_block_set));
imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');

if length(cur_block_set) > n_use_blocks
    cur_block_set = cur_block_set(1:n_use_blocks);
end
n_blocks = length(cur_block_set);


%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_stim_times = [];
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
all_spk_times = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);

trial_toffset = zeros(length(cur_block_set),1);
cur_spkind_offset = 0;
cur_toffset = 0;
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
            if ~isempty(all_stim_times)
                if any(cur_stim_times < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times];
            all_t_axis = [all_t_axis; cur_t_axis];
%             all_stim_mat = [all_stim_mat; cur_stim_mat];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
    end
    trial_cnt = trial_cnt + n_trials;
end


%% DEFINE DATA USED FOR ANALYSIS
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer);
NT = length(used_inds);

%% CREATE EVENT PREDICTORS FOR REAL AND SIM SACCADES (POOLING MICRO AND MACRO SACS)
Xblock = zeros(length(all_stim_times),n_blocks);
for i = 1:n_blocks
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end

%% PROCESS EYE TRACKING DATA
trial_toffset = zeros(length(cur_block_set),1);
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

%compute corrected eye data in bar-oriented frame
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,expt_bar_ori(cur_block_set),used_inds);

[saccades,et_params] = detect_saccades(corrected_eye_vals,all_eye_speed,all_eye_ts,et_params);

par_thresh = 4;
orth_thresh = 1;
[out_of_range] = detect_bad_fixation(corrected_eye_vals_interp,all_trialvec,used_inds,par_thresh,orth_thresh);
fract_out = length(out_of_range)/length(used_inds);
fprintf('Eliminating %.4f of data out of window\n',fract_out);
used_inds(ismember(used_inds,out_of_range)) = [];
NT = length(used_inds);

sac_start_times = [saccades(:).start_time];
sac_stop_times = [saccades(:).stop_time];
interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_start_inds(isnan(interp_sac_start_inds)) = 1;
interp_sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_stop_times));
interp_sac_stop_inds(isnan(interp_sac_stop_inds)) = 1;
sac_error = abs(sac_start_times - all_t_axis(interp_sac_start_inds)');
bad_sacs = find(isnan(interp_sac_start_inds) | sac_error > dt);
sac_start_times(bad_sacs) = [];
sac_stop_times(bad_sacs) = [];
saccades(bad_sacs) = [];
interp_sac_start_inds(bad_sacs) = [];
interp_sac_stop_inds(bad_sacs) = [];

%% CREATE SACCADE PREDICTOR MATS
saccade_start_inds = find(ismember(used_inds,interp_sac_start_inds));
used_saccade_set = find(ismember(interp_sac_start_inds,used_inds));
%nearest index in the used data set of the saccade stop time
saccade_stop_inds = round(interp1(used_inds,1:length(used_inds),interp_sac_stop_inds(used_saccade_set)))';

saccade_stop_inds(isnan(saccade_stop_inds)) = length(used_inds);

sac_amps = [saccades(:).amplitude];
is_micro = sac_amps(used_saccade_set) < 1;
big_sacs = find(~is_micro);
micro_sacs = find(is_micro);

saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));

%% DEFINE FIXATION POINTS
trial_start_inds = [1; find(diff(all_trialvec(used_inds)) ~= 0) + 1];
trial_end_inds = [find(diff(all_trialvec(used_inds)) ~= 0); NT];

fix_start_inds = sort([trial_start_inds; saccade_stop_inds]);
fix_stop_inds = sort([trial_end_inds; saccade_start_inds]);
fix_durs = (fix_stop_inds-fix_start_inds)*dt;
fix_start_inds(fix_durs<=0) = [];
fix_stop_inds(fix_durs <= 0) = [];
fix_durs(fix_durs <= 0) = [];
n_fixs = length(fix_start_inds);

%push the effects of saccades forward in time
sac_shift = round(0.05/dt);
pfix_start_inds = fix_start_inds;
pfix_stop_inds = fix_stop_inds;
for i = 1:length(saccade_start_inds)
    next_trial = trial_start_inds(find(trial_start_inds >= fix_start_inds(i),1,'first'));
    if next_trial > fix_start_inds(i) + sac_shift
        pfix_start_inds(i) = fix_start_inds(i) + sac_shift;
    end
    next_trial = trial_start_inds(find(trial_start_inds >= fix_stop_inds(i),1,'first'));
    if next_trial > fix_stop_inds(i) + sac_shift
        pfix_stop_inds(i) = fix_stop_inds(i) + sac_shift;
    end
end

%for all times within the (forward-projected) saccade, use 'prior'
%state-transition model
use_prior = zeros(NT,1);
for i = 1:n_fixs-1
    use_prior((pfix_stop_inds(i)+1):pfix_start_inds(i+1)) = 1;
end
fix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
end

%% Create set of TR and XV trials
use_trials = unique(all_trialvec(used_inds));
nuse_trials = length(use_trials);

n_xv_trials = round(xv_frac*nuse_trials);
xv_trials = randperm(nuse_trials);
xv_trials(n_xv_trials+1:end) = [];
xv_trials = use_trials(xv_trials);
tr_trials = setdiff(use_trials,xv_trials);

tr_inds = find(ismember(all_trialvec(used_inds),tr_trials));
xv_inds = find(ismember(all_trialvec(used_inds),xv_trials));

full_inds = sort([tr_inds; xv_inds]);

%don't use separate xv set for eye-tracking
tr_inds = full_inds;

%% DEFINE POINTS TO ALLOW RAPID CHANGES (TRIAL STARTS AFTER SACCADES)
%push the effects of saccades forward in time
sac_shift = round(0.05/dt);
pfix_start_inds = fix_start_inds;
pfix_stop_inds = fix_stop_inds;
for i = 1:length(fix_start_inds)
    next_trial = trial_start_inds(find(trial_start_inds >= fix_start_inds(i),1,'first'));
    if next_trial > fix_start_inds(i) + sac_shift
        pfix_start_inds(i) = fix_start_inds(i) + sac_shift;
    end
    next_trial = trial_start_inds(find(trial_start_inds >= fix_stop_inds(i),1,'first'));
    if next_trial > fix_stop_inds(i) + sac_shift
        pfix_stop_inds(i) = fix_stop_inds(i) + sac_shift;
    end
end

%for all times within the (forward-projected) saccade, use 'prior'
%state-transition model
use_prior = zeros(NT,1);
for i = 1:n_fixs-1
    use_prior((pfix_stop_inds(i)+1):pfix_start_inds(i+1)) = 1;
end
fix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
end

trial_ids = nan(NT,1);
for ii = 1:n_trials
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    trial_ids(cur_inds) = ii;
end


%%
use_fix_n = 3;
use_drift_n = 1;

cd(anal_dir);
load(left_anal_name);
left_it_post_mean = it_fix_post_mean(use_fix_n,:);
left_it_post_std = it_fix_post_std(use_fix_n,:);
left_drift_post_mean = drift_post_mean(use_drift_n,:);
left_drift_post_std = drift_post_std(use_drift_n,:);

left_hr_name = 'left_eyecorr_hbar_highres';
load(left_hr_name,'drift_*');
left_hr_drift_mean = drift_post_mean(end,:);
left_hr_drift_std = drift_post_std(end,:);

cd(anal_dir);
load(right_anal_name);
right_it_post_mean = it_fix_post_mean(use_fix_n,:);
right_it_post_std = it_fix_post_std(use_fix_n,:);
right_drift_post_mean = drift_post_mean(use_drift_n,:);
right_drift_post_std = drift_post_std(use_drift_n,:);

right_hr_name = 'right_eyecorr_hbar_highres';
load(right_hr_name,'drift_*');
right_hr_drift_mean = drift_post_mean(end,:);
right_hr_drift_std = drift_post_std(end,:);

%% MAKE FINAL DRIFT AND FIXATION CORRECTIONS

finL_fix_corr = nan(NT,1);
finL_fix_std = nan(NT,1);
finR_fix_corr = nan(NT,1);
finR_fix_std = nan(NT,1);
finL_fix_corr(~isnan(fix_ids)) = left_it_post_mean(fix_ids(~isnan(fix_ids)));
finL_fix_corr = interp1(find(~isnan(fix_ids)),finL_fix_corr(~isnan(fix_ids)),1:NT);
finR_fix_corr(~isnan(fix_ids)) = right_it_post_mean(fix_ids(~isnan(fix_ids)));
finR_fix_corr = interp1(find(~isnan(fix_ids)),finR_fix_corr(~isnan(fix_ids)),1:NT);

finL_fix_std(~isnan(fix_ids)) = left_it_post_std(fix_ids(~isnan(fix_ids)));
finL_fix_std = interp1(find(~isnan(fix_ids)),finL_fix_std(~isnan(fix_ids)),1:NT);
finR_fix_std(~isnan(fix_ids)) = right_it_post_std(fix_ids(~isnan(fix_ids)));
finR_fix_std = interp1(find(~isnan(fix_ids)),finR_fix_std(~isnan(fix_ids)),1:NT);

finL_fix_corr = finL_fix_corr*sp_dx;
finL_fix_std = finL_fix_std*sp_dx;
finR_fix_corr = finR_fix_corr*sp_dx;
finR_fix_std = finR_fix_std*sp_dx;

% finL_drift_corr = left_drift_post_mean*sp_dx;
% finL_drift_std = left_drift_post_std*sp_dx;
% finR_drift_corr = right_drift_post_mean*sp_dx;
% finR_drift_std = right_drift_post_std*sp_dx;
finL_drift_corr = left_hr_drift_mean*sp_dx/2;
finL_drift_std = left_hr_drift_std*sp_dx/2;
finR_drift_corr = right_hr_drift_mean*sp_dx/2;
finR_drift_std = right_hr_drift_std*sp_dx/2;
min_fix_dur = 0.1;
fix_inds = [];
long_fix_inds = [];
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
%     if length(cur_inds) > sac_shift
%         finL_drift_corr(cur_inds(1:end-sac_shift+1)) = finL_drift_corr(cur_inds(sac_shift:end));
%         finL_drift_std(cur_inds(1:end-sac_shift+1)) = finL_drift_std(cur_inds(sac_shift:end));
%         finR_drift_corr(cur_inds(1:end-sac_shift+1)) = finR_drift_corr(cur_inds(sac_shift:end));
%         finR_drift_std(cur_inds(1:end-sac_shift+1)) = finR_drift_std(cur_inds(sac_shift:end));
%     end
    fix_inds = [fix_inds cur_inds];
    if length(cur_inds)*dt >= min_fix_dur
       long_fix_inds = [long_fix_inds cur_inds]; 
    end
end

for ii = 1:length(trial_start_inds)
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    finL_drift_corr(cur_inds(1:end-sac_shift)) = finL_drift_corr(cur_inds(sac_shift+1:end));
    finR_drift_corr(cur_inds(1:end-sac_shift)) = finR_drift_corr(cur_inds(sac_shift+1:end));
    finL_drift_std(cur_inds(1:end-sac_shift)) = finL_drift_std(cur_inds(sac_shift+1:end));
    finR_drift_std(cur_inds(1:end-sac_shift)) = finR_drift_std(cur_inds(sac_shift+1:end));
end
finL_drift_corr = interp1(find(~isnan(fix_ids)),finL_drift_corr(~isnan(fix_ids)),1:NT);
finR_drift_corr = interp1(find(~isnan(fix_ids)),finR_drift_corr(~isnan(fix_ids)),1:NT);
finL_drift_std = interp1(find(~isnan(fix_ids)),finL_drift_std(~isnan(fix_ids)),1:NT);
finR_drift_std = interp1(find(~isnan(fix_ids)),finR_drift_std(~isnan(fix_ids)),1:NT);


finL_tot_corr = finL_fix_corr + finL_drift_corr;
finR_tot_corr = finR_fix_corr + finR_drift_corr;
finL_tot_std = sqrt(finL_fix_std.^2 + finL_drift_std.^2);
finR_tot_std = sqrt(finR_fix_std.^2 + finR_drift_std.^2);

%%
smooth_eyepos = corrected_eye_vals_interp;
smooth_eyepos(isnan(smooth_eyepos)) = 0;

% %smooth out fast transients in eye signal
% eye_smooth_sig = round(0.025/dt);
% interp_inds = [];
% for ii = 1:n_fixs
%     cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
%     if length(cur_inds) > eye_smooth_sig*5;
%         smooth_eyepos(used_inds(cur_inds),1) = jmm_smooth_1d_cor(smooth_eyepos(used_inds(cur_inds),1),eye_smooth_sig,2);
%         smooth_eyepos(used_inds(cur_inds),2) = jmm_smooth_1d_cor(smooth_eyepos(used_inds(cur_inds),2),eye_smooth_sig,2);
%         smooth_eyepos(used_inds(cur_inds),3) = jmm_smooth_1d_cor(smooth_eyepos(used_inds(cur_inds),3),eye_smooth_sig,2);
%         smooth_eyepos(used_inds(cur_inds),4) = jmm_smooth_1d_cor(smooth_eyepos(used_inds(cur_inds),4),eye_smooth_sig,2);
%     end
%     interp_inds = [interp_inds; cur_inds'];
% end
% 
% interp_inds = unique(interp_inds);
% smooth_eyepos(used_inds,:) = interp1(used_inds(interp_inds),smooth_eyepos(used_inds(interp_inds),:),used_inds);
% 
% %subtract off within-trial median
% for ii = 1:length(trial_start_inds)
%     cur_inds = trial_start_inds(ii):trial_end_inds(ii);
%     smooth_eyepos(used_inds(cur_inds),:) = bsxfun(@minus,smooth_eyepos(used_inds(cur_inds),:),median(smooth_eyepos(used_inds(cur_inds),:)));
% end

for ii = 1:length(cur_block_set)
    cur_inds = used_inds(all_blockvec(used_inds) == ii);
    smooth_eyepos(cur_inds,:) = bsxfun(@minus,smooth_eyepos(cur_inds,:),nanmedian(smooth_eyepos(cur_inds,:)));
end
    
smooth_eyepos = smooth_eyepos(used_inds,:);

%% GENERATE EXAMPLE COMPARISON FIGURE
%Trial 12 is current example trial

close all
n_trials = length(unique(all_trialvec));

H = figure();
% for tt = 1:n_trials
% tt = 12;
    fprintf('Trial %d of %d\n',tt,n_trials);
    for tt = [12 15 29 170 192 242]
    uu = find(all_trialvec(used_inds) == tt);
    if ~isempty(uu)
        bt = all_t_axis(used_inds(uu(1)));
        et = all_t_axis(used_inds(uu(end)));
        dur = et-bt;
        if dur > 3
            cur_sac_inds = find(ismember(saccade_start_inds,uu));
            rel_sac_start_times = all_t_axis(used_inds(saccade_start_inds(cur_sac_inds))) - bt;
            rel_sac_end_times = all_t_axis(used_inds(saccade_stop_inds(cur_sac_inds))) - bt;
            
            hold on
            h2=shadedErrorBar(all_t_axis(used_inds(uu))-bt,finR_tot_corr(uu),finR_tot_std(uu),{'color','r'},0);
            h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,finL_tot_corr(uu),finL_tot_std(uu),{'color','b'},0);
            h3=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),2),'k');
            delete(h1.edge); delete(h2.edge);
            line([0 dur],[0 0],'color','k','linestyle','--');
            legend([h1.mainLine h2.mainLine h3],{'Left-eye inferred','Right-eye inferred','Coil measured'})
            xlim([0 dur]);
            ylim([-0.4 0.4]);
            xlabel('Time (s)');
            ylabel('Orthoganol position (deg)');
            title(sprintf('Trial %d',tt));
            
            yl = ylim();
            for ii = 1:length(rel_sac_start_times)
                line(rel_sac_start_times([ii ii]),yl,'color',[0.2 0.8 0.2],'linestyle','--');
                %                 line(rel_sac_end_times([ii ii]),yl,'color',[0.2 0.8 0.2],'linestyle','--');
            end
            
            figufy(H);
            pause
            clf(H);
        end
    end
end

%% PRINT EXAMPLE COMPARISON FIGURE
fig_width = 4.86; %3.27 4.86 6.83
rel_height = 0.6;

fname = [fig_dir 'binoc_example_trace_t29.pdf'];
exportfig(H,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);


%% Compute left-right difference and plot joint distribution
LRdiff = finL_tot_corr(fix_inds) - finR_tot_corr(fix_inds);
xr = [-0.3 0.3];
sqrt_sc = false;
H = figure;
if sqrt_sc
    [h,dens_det] = DensityPlot_jmm(finL_tot_corr(fix_inds),finR_tot_corr(fix_inds),'ynormal','xrange',xr,'yrange',xr,'sd',[3 3],'sqrtsc');
else
    [h,dens_det] = DensityPlot_jmm(finL_tot_corr(fix_inds),finR_tot_corr(fix_inds),'ynormal','xrange',xr,'yrange',xr,'sd',[3 3]);
end
line(xr,xr,'color','w');
xlabel('Left eye position (deg)','fontsize',12);
ylabel('Right eye position (deg)','fontsize',12);
set(gca,'fontname','arial','fontsize',10);
axis square
fillPage(gcf,'papersize',[5 5]);
set(gca,'xtick',xr(1):0.1:xr(2),'ytick',xr(1):0.1:xr(2));

n_bins = 45;
h2 = figure; hold on
bin_ax = linspace(xr(1),xr(2),n_bins + 1);
y = histc(finL_tot_corr(fix_inds),bin_ax);
y = y/sum(y);
stairs(bin_ax,y);

y = histc(finR_tot_corr(fix_inds),bin_ax);
y = y/sum(y);
stairs(bin_ax,y,'r');

y = histc(LRdiff,bin_ax);
y = y/sum(y);
stairs(bin_ax,y,'k');

xlim(bin_ax([1 end]));
set(gca,'xtick',[],'ytick',[]);
yl = ylim();
line([0 0],yl,'color','k');
title('Left position','fontsize',10);
fillPage(gcf,'papersize',[4 4]);

%% PRINT OUT BINOC DISTRIBUTION FIGURE
fig_width = 3.27; %3.27 4.86 6.83
rel_height = 0.8;

figufy(H);
fname = [fig_dir 'binoc_LRjoint_dist2.pdf'];
exportfig(H,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

figufy(h2);
fname = [fig_dir 'binoc_LRmarg_dist2.pdf'];
exportfig(h2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
sparse_blocks = find(expt_dds(cur_block_set) == 12);
sparse_inds = find(ismember(all_blockvec(used_inds),sparse_blocks));
sparse_fix_inds = sparse_inds(ismember(sparse_inds,fix_inds));
sparse_long_fix_inds = sparse_inds(ismember(sparse_inds,long_fix_inds));

% min_unc = 0.001;
% avg_unc = sqrt(finL_tot_std.^2 + finR_tot_std.^2);
% avg_unc(avg_unc < min_unc) = min_unc;
% weights = 1./avg_unc;
% [weighted_r,raw_r] = weighted_correlation(finR_tot_corr(sparse_inds),finL_tot_corr(sparse_inds),weights(sparse_inds));

%% ROBUST STATS FIG
close all

max_unc = 0.04;
% max_unc = Inf;
tot_unc = sqrt(finR_tot_std(sparse_inds).^2 + finL_tot_std(sparse_inds).^2);
cert_data = find(tot_unc < max_unc);

xx = linspace(-0.7,0.7,100);
inf_bdisp = finL_tot_corr(sparse_inds) - finR_tot_corr(sparse_inds);

yobs = ksdensity(inf_bdisp,xx);
ycert = ksdensity(inf_bdisp(cert_data),xx);
g1 = normpdf(xx,0,std(inf_bdisp));
g2 = normpdf(xx,0,robust_std_dev(inf_bdisp));

h1=figure;
plot(xx,yobs);
hold on
plot(xx,ycert,'g');
plot(xx,g1,'r',xx,g2,'k');
set(gca,'yscale','log');
ylim([1e-3 15])
xlim(xx([1 end]));
xlabel('Vergence error (deg)');
ylabel('Probability');

xx = linspace(-0.25,0.25,100);
yobs = ksdensity(inf_bdisp,xx);
ycert = ksdensity(inf_bdisp(cert_data),xx);
g1 = normpdf(xx,0,std(inf_bdisp));
g2 = normpdf(xx,0,robust_std_dev(inf_bdisp));
h2=figure;
plot(xx,yobs);
hold on
plot(xx,ycert,'g');
plot(xx,g1,'r',xx,g2,'k');
xlim([-0.2 0.2]);
xlabel('Vergence error (deg)');
ylabel('Probability');

h3=figure;
[n,x] = hist(tot_unc,100);
n = n/sum(n);
bar(x,n);
yl = ylim();
line([max_unc max_unc],yl,'color','r');
xlabel('Uncertainty (deg)');
ylabel('Relative frequency');

%% PRINT OUT BINOC DISTRIBUTION FIGURE
fig_width = 3.27; %3.27 4.86 6.83
rel_height = 0.8;

figufy(h1);
fname = [fig_dir 'verg_err_robust_dist_log.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h1)

figufy(h2);
fname = [fig_dir 'verg_err_robust_dist_lin.pdf'];
exportfig(h2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h2)

figufy(h3);
fname = [fig_dir 'verg_err_robust_unc_dist.pdf'];
exportfig(h3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h3)

%% GENERATE BINOC DISPARITY DISTRIBUTION FIGURE
close all
% max_unc = 0.04;
% max_unc = Inf;
% cert_data = sparse_inds(finR_tot_std(sparse_inds) < max_unc & finL_tot_std(sparse_inds) < max_unc);
% cert_data = long_fix_inds;
cert_data = fix_inds;

unc_corr_fix = corr(finL_fix_corr(cert_data)',finR_fix_corr(cert_data)','type','spearman');
unc_corr_tot = corr(finL_tot_corr(cert_data)',finR_tot_corr(cert_data)','type','spearman');

% measured_pos = smooth_eyepos(:,2);
measured_pos = corrected_eye_vals_interp(used_inds,2);

% %for corrected eye position data, use only complete trials to avoid
% %under-estimating due to trial-based median subtraction
% all_trial_durs = (all_trial_end_times - all_trial_start_times);
% usable_trials = find(all_trial_durs > 3.5);
% usable_trial_inds = find(ismember(all_trialvec(used_inds),usable_trials));

% bin_ax = linspace(-0.15,0.15,50);
bin_ax = linspace(-0.3,0.3,50);
LRdiff = finL_tot_corr(cert_data) - finR_tot_corr(cert_data);
LRdiff_unc = sqrt(finL_tot_std(cert_data).^2 + finR_tot_std(cert_data).^2);

H=figure; hold on
yD = histc(LRdiff,bin_ax);
yD = yD/sum(yD);
stairs(bin_ax,yD,'g','linewidth',1.5)
yM = histc(measured_pos(cert_data),bin_ax);
yM = yM/sum(yM);
yL = histc(finL_tot_corr(cert_data),bin_ax);
yL = yL/sum(yL);
yR = histc(finR_tot_corr(cert_data),bin_ax);
yR = yR/sum(yR);
% stairs(bin_ax,yM,'r','linewidth',1.5)
stairs(bin_ax,yL,'b','linewidth',1.5)
stairs(bin_ax,yR,'r','linewidth',1.5)
xlim(bin_ax([1 end]));
yl = ylim();
line([0 0],yl,'color','k','linestyle','--');
xlabel('Position difference (deg)');
ylabel('Relative frequency');
legend('Inferred binocular disparity','Measured fixation disparity');


figure
bin_ax = linspace(0,0.3,50);
yD = histc(abs(LRdiff),bin_ax);
yD = cumsum(yD)/length(LRdiff);
plot(bin_ax,yD);
xlim([0 0.3]);
line([1/60 1/60],[0 1])
line([2/60 2/60],[0 1])
%% PRINT OUT BINOC DISTRIBUTION FIGURE
fig_width = 3.27; %3.27 4.86 6.83
rel_height = 0.8;

figufy(H);
fname = [fig_dir 'binoc_dist_compare.pdf'];
exportfig(H,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);


%%
figure;
y = histc(finR_tot_corr(sparse_inds),bin_ax);
y = y/sum(y);
stairs(bin_ax,y);
xlim(bin_ax([1 end]));
set(gca,'xtick',[],'ytick',[]);
yl = ylim();
line([0 0],yl,'color','k');
title('Right position','fontsize',10);
fillPage(gcf,'papersize',[4 4]);
fname = sprintf('Right_dist_bar%d',bar_ori);

LRdiff = finL_tot_corr(sparse_inds) - finR_tot_corr(sparse_inds);
figure;
y = histc(LRdiff,bin_ax);
y = y/sum(y);
stairs(bin_ax,y);
xlim(bin_ax([1 end]));
set(gca,'xtick',[],'ytick',[]);
yl = ylim();
line([0 0],yl,'color','k');
title('Diff position','fontsize',10);
fillPage(gcf,'papersize',[4 4]);
fname = sprintf('LR_diff_dist_bar%d',bar_ori);

