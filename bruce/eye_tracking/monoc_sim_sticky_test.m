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

bar_ori = 0;


% noise_model = 'gauss';
noise_model = 'log_exp';
n_squared_filts = 2;

recompute_init_mods = 0;
use_measured_pos = 0;
use_sac_kerns = 1;

%dont fit stim models using these blocks
if Expt_num == 86
    ignore_blocks = [16 17 28 30]; %G086
else
    ignore_blocks = [];
end

%%
xv_frac = 0;
n_fix_inf_it = 5;
n_drift_inf_it = 3;

eps = -1e3;
fix_prior_sigma = 0.15;
drift_prior_sigma = 0.004;
drift_jump_sigma = 0.05; %0.05 start
drift_dsf = 2;

% fix_prior_sigma = 0.15*sim_gain;
% drift_sigma = 0.015*sim_gain; %.015
% drift_jump_sigma = 0.05*sim_gain; %0.05 start

stim_fs = 100; %in Hz
dt = 0.01;
flen = 12;
use_nPix = 16;
full_nPix = 36;
stim_params = NIMcreate_stim_params([flen full_nPix],dt);
Fr = 1;
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
max_shift = round(16*spatial_usfac); %up to 18
dshift = 1;
shifts = -max_shift:dshift:max_shift;
n_shifts = length(shifts);
zero_frame = find(shifts==0);

temp_dx = [-2*max_shift:dshift:2*max_shift];
shift_dx_edges = [temp_dx*sp_dx-dshift*sp_dx/2 temp_dx(end)*sp_dx+dshift*sp_dx/2];
ovshift_dx_edges = [shifts*sp_dx-dshift*sp_dx/2 shifts(end)*sp_dx+dshift*sp_dx/2];


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
all_trial_blocknums = [];

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
n_unique_stims = 500;
stim_set = ceil(rand(n_unique_stims,1)*size(all_stim_mat,1));
unique_stims = all_stim_mat(stim_set,:);

rand_stim_seq = ceil(rand(size(all_stim_mat,1),1)*n_unique_stims);
new_stim_mat = unique_stims(rand_stim_seq,:);

%%
full_nPix_us = spatial_usfac*full_nPix;
if spatial_usfac == 2
new_stimmat_up = zeros(size(all_stim_mat,1),full_nPix_us);
for ii = 1:size(all_stim_mat,2)
    new_stimmat_up(:,2*(ii-1)+1) = new_stim_mat(:,ii);
    new_stimmat_up(:,2*(ii-1)+2) = new_stim_mat(:,ii);
end
elseif spatial_usfac == 1
    new_stimmat_up = new_stim_mat;
else
    error('Unsupported spatial Up-sample factor!');
end
stim_params_us = NMMcreate_stim_params([flen full_nPix_us],dt);

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
% saccade_stop_inds = find(ismember(used_inds,interp_sac_stop_inds));
used_saccade_set = find(ismember(interp_sac_start_inds,used_inds));
saccade_stop_inds = round(interp1(used_inds,1:length(used_inds),interp_sac_stop_inds(used_saccade_set)))';

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
trial_end_inds = [find(diff(all_trialvec(used_inds)) ~= 0); NT];

fix_start_inds = sort([trial_start_inds; saccade_stop_inds]);
fix_stop_inds = sort([trial_end_inds; saccade_start_inds]);
fix_durs = fix_stop_inds-fix_start_inds;
fix_start_inds(fix_durs==0) = [];
fix_stop_inds(fix_durs == 0) = [];
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

%% INCORPORATE INFERRED EYE-POSITIONS AND MAKE SIMULATED ERRORS
cd(anal_dir)
true_eyedata_name = 'monoc_eyecorr_hbar_FinFin2.mat';
load(true_eyedata_name,'*mods*','et_tr_set','it_*');

true_mods = dit_mods{end};
tr_set = et_tr_set;
n_tr_chs = length(tr_set);

true_moddata_name = 'monoc_eyecorr_hbar_mods_FinFin.mat';
load(true_moddata_name,'all_mod_SU');
true_mod_SU = all_mod_SU;

sim_eyepos = nan(NT,1);
sim_eyepos(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
sim_eyepos = interp1(find(~isnan(fix_ids)),sim_eyepos(~isnan(fix_ids)),1:NT);
sim_eyepos(isnan(sim_eyepos)) = 0;

sim_eyepos_raw = sim_eyepos;

sim_eyepos_us_rnd = round(sim_eyepos);
sim_eyepos_rnd = round(sim_eyepos/spatial_usfac);

new_stimmat_true = new_stimmat_up;
for ii = 1:length(sim_eyepos)
    new_stimmat_true(ii,:) = shift_matrix_Nd(new_stimmat_up(ii,:),-sim_eyepos_us_rnd(ii),2);
end
all_Xmat_true = create_time_embedding(new_stimmat_true,stim_params_us);

%%
flen = 1;
stim_params_us = NMMcreate_stim_params([flen full_nPix_us],dt);
stim_params = NMMcreate_stim_params([flen full_nPix],dt);

[Xinds_up,Tinds] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
if spatial_usfac > 1
use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1)-1/spatial_usfac & Xinds_up(:) <= cur_use_pix(end));
else
    use_kInds_up = use_kInds;
end
use_kInds_back = find(ismember(Xinds_up(use_kInds_up),cur_use_pix));


[Xinds,Tinds] = meshgrid(1:full_nPix,1:flen);
mod_use_kInds = find(ismember(Xinds(:),cur_use_pix));

[Xinds_up,Tinds] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
if spatial_usfac > 1
    mod_use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1)-1/spatial_usfac & Xinds_up(:) <= cur_use_pix(end));
else
    mod_use_kInds_up = mod_use_kInds;
end
mod_use_kInds_back = find(ismember(Xinds_up(mod_use_kInds_up),cur_use_pix));

%%

all_Xmat_us = create_time_embedding(new_stimmat_up,stim_params_us);
all_Xmat = create_time_embedding(new_stim_mat,stim_params);

%% SIMULATE SPIKING RESPONSES
tr_X{1} = all_Xmat_true(used_inds,use_kInds_up);
tr_X{2} = Xblock(used_inds,:);
tr_X{3} = Xsac;
tr_X{4} = Xmsac;

% make Robs_mat
prate_mat = nan(length(used_inds),n_tr_chs);
Robs_mat = nan(length(used_inds),n_tr_chs);
for ss = 1:n_tr_chs
    [~, ~, pred_rate] = NMMmodel_eval(true_mods(tr_set(ss)),[], tr_X);
    if all_mod_SU(tr_set(ss)) == 0
        su_probe_ind = find(su_probes == tr_set(ss));
        if ~isempty(su_probe_ind)
            cur_used_blocks = find(su_used_blocks(:,su_probe_ind) == 0); %blocks when NO SU
        else
            cur_used_blocks = 1:n_blocks; %blocks when NO SU
        end
        cur_use = find(ismember(all_blockvec(used_inds),cur_used_blocks));
    else
        su_ind = find(su_probes == all_mod_SU(tr_set(ss)));
        cur_used_blocks = find(su_used_blocks(:,su_ind)); %blocks when SU
        cur_use = find(ismember(all_blockvec(used_inds),cur_used_blocks));
    end
    prate_mat(cur_use,ss) = pred_rate(cur_use);
end
Robs_mat = poissrnd(prate_mat);

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

base_lambda_d2X = 10;
base_lambda_L1 = 0;
if strcmp(noise_model,'gauss')
    base_lambda_L1 = 0;
end
init_optim_p.optTol = 1e-5; init_optim_p.progTol = 1e-9;

sac_reg_params = NMMcreate_reg_params('lambda_d2T',sac_d2t);
if use_sac_kerns
    null_reg_params = NMMcreate_reg_params('lambda_d2T',[0; sac_d2t; sac_d2t],'lambda_L2',[block_L2; 0; 0]);
else
    null_reg_params = NMMcreate_reg_params('lambda_L2',block_L2);
end

mod_signs = ones(1,n_squared_filts+2);
NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts) {'lin'}];
init_d2XT = [ones(n_squared_filts+1,1); 0;];
init_L2 = [zeros(n_squared_filts+1,1); block_L2];
init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_d2XT,'lambda_L2',init_L2);
init_Xtargs = [ones(n_squared_filts+1,1); 2];
if strcmp(noise_model,'gauss')
    spk_NL = 'linear';
else
    spk_NL = 'logexp';
end
init_filts = cell(length(mod_signs),1);
cd(anal_dir);

% if ~exist(['./' mod_data_name '.mat'],'file') || recompute_init_mods == 1
    for ss = 1:length(tr_set)
        fprintf('Computing base LLs for unit %d of %d\n',ss,n_tr_chs);
        cur_tr_inds = tr_inds(~isnan(Robs_mat(tr_inds,ss)));
        tr_NT = length(cur_tr_inds);
        if ~isempty(cur_tr_inds)
            Robs = Robs_mat(cur_tr_inds,ss);
            
            tr_X{1} = all_Xmat(used_inds(cur_tr_inds),mod_use_kInds);
            tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
            
            if use_sac_kerns
                tr_X{3} = Xsac(cur_tr_inds,:);
                tr_X{4} = Xmsac(cur_tr_inds,:);
                null_mod = NMMinitialize_model(null_stim_params,[1 1 1],repmat({'lin'},1,3),null_reg_params,[1 2 3]);
            else
                null_mod = NMMinitialize_model(null_stim_params,[1],{'lin'},null_reg_params,[1]);
            end
            
            null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));
            all_nullmod(ss) = null_mod;
            
            gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs,[],spk_NL);
            gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent,init_optim_p);
            
            tr_X{1} = all_Xmat_us(used_inds(cur_tr_inds),mod_use_kInds_up);
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
            gqm2 = NMMinitialize_model(fin_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs,init_filts,spk_NL);
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
            if ~strcmp(noise_model,'gauss')
                all_mod_fits_withspkNL(ss) = NMMfit_logexp_spkNL(gqm2,Robs,tr_X);
            else
                all_mod_fits_withspkNL(ss) = all_mod_fits(ss);
            end
            
            LL = NMMmodel_eval(all_mod_fits_withspkNL(ss),Robs,tr_X);
            null_LL(ss) = NMMmodel_eval(null_mod,Robs,tr_X(2:end));
            all_mod_LLimp(ss) = (LL-null_LL(ss))/log(2);
        end
    end
    save('temp_mods','all_mod*','all_nullmod','su_probes','null_LL','*_trials');
% else
%     fprintf('Loading pre-computed initial models\n');
%     load(mod_data_name);
%     
%     tr_inds = find(ismember(all_trialvec(used_inds),tr_trials));
% end

%%

%%
%overall prior on shifts
fix_prior = diff(erf(ovshift_dx_edges/(sqrt(2)*fix_prior_sigma)));
fix_prior = log(fix_prior/sum(fix_prior));
eps = -1e3;
fix_prior(fix_prior < eps) = eps;

%generate shift matrices. Must be applied to the stimulus (not the filters)
It = speye(flen);
shift_mat = cell(n_shifts,1);
for xx = 1:n_shifts
    temp = spdiags( ones(full_nPix_us,1), -shifts(xx), full_nPix_us, full_nPix_us);
    temp = kron(temp,It);
    shift_mat{xx} = temp(:,use_kInds_up);
end

%%
clear it_*
it_mods{1} = all_mod_fits;
it_mods_spkNL{1} = all_mod_fits_withspkNL;
it_LLimp(1,:) = all_mod_LLimp;
it_fix_sigma(1) = fix_prior_sigma;
for nn = 1:n_fix_inf_it
    fprintf('Inferring fixation corrections, iter %d of %d\n',nn,n_fix_inf_it);
    
    %% PREPROCESS MODEL COMPONENTS
    filt_bank = zeros(n_tr_chs,klen_us,n_squared_filts+1);
    lin_kerns = nan(n_tr_chs,n_blocks);
    if use_sac_kerns
        sac_kerns = nan(n_tr_chs,n_sac_bins);
        msac_kerns = nan(n_tr_chs,n_sac_bins);
    end
    mod_spkNL_params = nan(n_tr_chs,3);
    for ss = 1:n_tr_chs
        cur_Xtargs = [it_mods{nn}(ss).mods(:).Xtarget];
        cur_k = [it_mods{nn}(ss).mods(cur_Xtargs == 1).filtK];
        n_used_filts = size(cur_k,2);
        filt_bank(ss,:,1:n_used_filts) = cur_k;
        mod_spkNL_params(ss,:) = it_mods_spkNL{nn}(ss).spk_NL_params;
        lin_kerns(ss,:) = it_mods{nn}(ss).mods(cur_Xtargs == 2).filtK;
        if use_sac_kerns
            sac_kerns(ss,:) = it_mods{nn}(ss).mods(cur_Xtargs == 3).filtK;
            msac_kerns(ss,:) = it_mods{nn}(ss).mods(cur_Xtargs == 4).filtK;
        end
    end
    filt_bank = permute(filt_bank,[2 1 3]);
    
    %indicator predictions
    block_out = Xblock(used_inds,:)*lin_kerns';
    if use_sac_kerns
        sac_out = Xsac*sac_kerns';
        msac_out = Xmsac*msac_kerns';
    end

    %% ESTIMATE LL for each shift in each stimulus frame
    cur_Xmat = all_Xmat_us(used_inds,:);
    %precompute LL at all shifts for all units
    frame_LLs = nan(NT,n_shifts);
    for xx = 1:length(shifts)
        fprintf('Shift %d of %d\n',xx,n_shifts);
        cur_stim_shift = cur_Xmat*shift_mat{xx};
        
        %outputs of stimulus models at current X-matrix shift
        gfuns = ones(NT,n_tr_chs);
        gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,1)');
        gfuns = gfuns + cur_stim_shift*squeeze(filt_bank(:,:,1));
        for ff = 2:(n_squared_filts+1)
            gfuns = gfuns + (cur_stim_shift*squeeze(filt_bank(:,:,ff))).^2;
        end
        
        %add contributions from extra lin kernels
        gfuns = gfuns + block_out;
        if use_sac_kerns
            gfuns = gfuns + sac_out + msac_out;
        end
        
        %incorporate beta
        if strcmp(noise_model,'gauss')
            pred_rate = gfuns;
            frame_LLs(:,xx) = -squeeze(nansum((Robs_mat-pred_rate).^2,2));
        else
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
    end
    
    %% INFER MICRO-SAC SEQUENCE
    fix_LLs = nan(n_fixs,n_shifts);
    for ii = 1:n_fixs
        cur_inds = pfix_start_inds(ii):pfix_stop_inds(ii);
        fix_LLs(ii,:) = sum(frame_LLs(cur_inds,:));
    end
    
    lPost = bsxfun(@plus,fix_LLs,fix_prior);
    lPost = bsxfun(@minus,lPost,logsumexp(lPost,2));
    Post = exp(lPost);
    it_fix_sigma(nn) = sqrt(mean(sum(bsxfun(@times,Post,shifts*sp_dx).^2,2)));

    it_fix_post_mean(nn,:) = sum(bsxfun(@times,Post,shifts),2);
    cur_diff = bsxfun(@minus,it_fix_post_mean(nn,:)',shifts).^2;
    it_fix_post_std(nn,:) = sqrt(sum(cur_diff.*Post,2));
    
    %back-project saccade-times
    all_fix_post_mean_cor = nan(NT,1);
    all_fix_post_mean_cor(~isnan(fix_ids)) = it_fix_post_mean(nn,fix_ids(~isnan(fix_ids)));
    all_fix_post_mean_cor = interp1(find(~isnan(fix_ids)),all_fix_post_mean_cor(~isnan(fix_ids)),1:NT);

    %% RECOMPUTE XMAT
    all_shift_stimmat_up = all_stimmat_up;
    for i=1:NT
        all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-round(all_fix_post_mean_cor(i)),2);
    end
    all_Xmat_up_fixcor = create_time_embedding(all_shift_stimmat_up,stim_params_us);
    all_Xmat_up_fixcor = all_Xmat_up_fixcor(used_inds,:);
    
    %% REFIT ALL CELLS
    cur_X{1} = all_Xmat_up_fixcor(:,use_kInds_up);
    cur_X{2} = Xblock(used_inds,:);
    if use_sac_kerns
        cur_X{3} = Xsac;
        cur_X{4} = Xmsac;
    end
    uncor_X = cur_X;
    uncor_X{1} = all_Xmat_us(used_inds,use_kInds_up);
    
    silent = 1;
    for ss = 1:length(tr_set)
        fprintf('Refitting model for tr cell %d of %d\n',ss,length(tr_set));
        cur_tr_uset = tr_inds(~isnan(Robs_mat(tr_inds,ss)));
        
        tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
        
        it_mods{nn+1}(ss) = it_mods{nn}(ss);
        it_mods{nn+1}(ss) = NMMfit_filters(it_mods{nn+1}(ss),Robs_mat(cur_tr_uset,ss),tr_X,[],[],silent); %fit stimulus filters
        
        %refit spk NL
        if ~strcmp(noise_model,'gauss')
            it_mods_spkNL{nn+1}(ss) = NMMfit_logexp_spkNL(it_mods{nn+1}(ss),Robs_mat(cur_tr_uset,ss),tr_X);
        else
            it_mods_spkNL{nn+1}(ss) = it_mods{nn+1}(ss);
        end
        newLL = NMMmodel_eval(it_mods_spkNL{nn+1}(ss),Robs_mat(cur_tr_uset,ss),tr_X);
        it_LLimp(nn+1,ss) = (newLL - null_LL(ss))/log(2);
        
        if nn > 1
            fprintf('Original: %.4f  Prev: %.4f  New: %.4f\n',all_mod_LLimp(ss),it_LLimp(nn,ss),it_LLimp(nn+1,ss));
        else
            fprintf('Original: %.4f  New: %.4f\n',all_mod_LLimp(ss),it_LLimp(nn+1,ss));
        end
    end
end

%% NOW INFER DRIFT CORRECTIONS
%overall prior on shifts
drift_jump_prior = diff(erf(ovshift_dx_edges/(sqrt(2)*drift_jump_sigma)));
drift_jump_prior = log(drift_jump_prior/sum(drift_jump_prior));
eps = -1e3;
drift_jump_prior(drift_jump_prior < eps) = eps;
lA_tflip = repmat(drift_jump_prior,n_shifts,1);

% cdist = squareform(pdist(shifts'*sp_dx));
% lA = -cdist.^2/(2*drift_sigma^2);
% lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize
gprobs = diff(erf(shift_dx_edges/(drift_prior_sigma*sqrt(2)*drift_dsf)));
lA = toeplitz(gprobs(n_shifts:end),gprobs(n_shifts:-1:1));
lA = log(bsxfun(@rdivide,lA,sum(lA,2)));
lA(lA < eps) = eps;

dit_mods{1} = it_mods{n_fix_inf_it+1};
dit_mods_spkNL{1} = it_mods_spkNL{n_fix_inf_it+1};
dit_LLimp(1,:) = it_LLimp(n_fix_inf_it+1,:);
dit_drift_sigma(1) = drift_prior_sigma;
dit_jump_sigma(1) = drift_jump_sigma;
% dit_lA(1,:,:) = lA;
for nn = 1:n_drift_inf_it
    fprintf('Inferring drift corrections, iter %d of %d\n',nn,n_drift_inf_it);
    
    %% PREPROCESS MODEL COMPONENTS
    filt_bank = zeros(n_tr_chs,klen_us,n_squared_filts+1);
    lin_kerns = nan(n_tr_chs,n_blocks);
    if use_sac_kerns
        sac_kerns = nan(n_tr_chs,n_sac_bins);
        msac_kerns = nan(n_tr_chs,n_sac_bins);
    end
    mod_spkNL_params = nan(n_tr_chs,3);
    for ss = 1:n_tr_chs
        cur_Xtargs = [dit_mods{nn}(ss).mods(:).Xtarget];
        cur_k = [dit_mods{nn}(ss).mods(cur_Xtargs == 1).filtK];
        n_used_filts = size(cur_k,2);
        filt_bank(ss,:,1:n_used_filts) = cur_k;
        mod_spkNL_params(ss,:) = dit_mods_spkNL{nn}(ss).spk_NL_params;
        lin_kerns(ss,:) = dit_mods{nn}(ss).mods(cur_Xtargs == 2).filtK;
        if use_sac_kerns
            sac_kerns(ss,:) = dit_mods{nn}(ss).mods(cur_Xtargs == 3).filtK;
            msac_kerns(ss,:) = dit_mods{nn}(ss).mods(cur_Xtargs == 4).filtK;
        end
    end
    filt_bank = permute(filt_bank,[2 1 3]);
    
    %indicator predictions
    block_out = Xblock(used_inds,:)*lin_kerns';
    if use_sac_kerns
        sac_out = Xsac*sac_kerns';
        msac_out = Xmsac*msac_kerns';
    end

    %% ESTIMATE LL for each shift in each stimulus frame
    
    %precompute LL at all shifts for all units
    frame_LLs = nan(NT,n_shifts);
    for xx = 1:length(shifts)
        fprintf('Shift %d of %d\n',xx,n_shifts);
        cur_stim_shift = all_Xmat_up_fixcor*shift_mat{xx};
        
        %outputs of stimulus models at current X-matrix shift
        gfuns = ones(NT,n_tr_chs);
        gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,1)');
        gfuns = gfuns + cur_stim_shift*squeeze(filt_bank(:,:,1));
        for ff = 2:(n_squared_filts+1)
            gfuns = gfuns + (cur_stim_shift*squeeze(filt_bank(:,:,ff))).^2;
        end
        
        %add contributions from extra lin kernels
        gfuns = gfuns + block_out;
        if use_sac_kerns
            gfuns = gfuns + sac_out + msac_out;
        end
        
        %incorporate beta
        if strcmp(noise_model,'gauss')
            pred_rate = gfuns;
            frame_LLs(:,xx) = -squeeze(nansum((Robs_mat-pred_rate).^2,2));
        else
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
    end
    
    %% INFER DRIFT CORRECTIONS
        
    lgamma = nan(NT,n_shifts);
    for ff = 1:n_fixs
        if mod(ff,100)==0
        fprintf('Fixation %d of %d\n',ff,n_fixs);
        end
            
        tset = find(fix_ids==ff)';
        ntset = length(tset);
        nt_pts = ceil(ntset/drift_dsf);
        tset_inds = 1+floor((0:(ntset-1))/drift_dsf);
        talpha=zeros(nt_pts,n_shifts);
        tbeta = zeros(nt_pts,n_shifts);
        
        tpt_loc = ceil(1:drift_dsf:nt_pts*drift_dsf);
        tpt_loc(end) = ntset;
        
        cur_LL_set = frame_LLs(tset,:);
        if mod(ntset,drift_dsf) ~= 0
            dangling_pts = nt_pts*drift_dsf-ntset;
            cur_LL_set = cat(1,cur_LL_set,zeros(dangling_pts,n_shifts));
        end
        cur_LL_set = reshape(cur_LL_set,[drift_dsf nt_pts n_shifts]);
        cur_LL_set = squeeze(sum(cur_LL_set,1));
        
        talpha(1,:) = drift_jump_prior + cur_LL_set(1,:);
        for t = 2:nt_pts
            talpha(t,:) = logmulexp(talpha(t-1,:),lA) + cur_LL_set(t,:);
        end
        
        tbeta(end,:)=log(ones(1,n_shifts));
        for t = (nt_pts-1):-1:1
            lf1 = tbeta(t+1,:) + cur_LL_set(t+1,:);
            tbeta(t,:) = logmulexp(lf1,lA');
        end
        temp_gamma = talpha + tbeta;
        if drift_dsf > 1
            if nt_pts > 1
                int_gamma = interp1(tpt_loc,temp_gamma,1:ntset);
                lgamma(tset,:) = int_gamma;
            else
                lgamma(tset,:) = repmat(temp_gamma,ntset,1);
            end
        else
            lgamma(tset,:) = temp_gamma;
        end
    end
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    
    gamma = exp(lgamma);
    drift_post_mean(nn,:) = sum(bsxfun(@times,gamma,shifts),2);
    drift_post_std(nn,:) = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - drift_post_mean(nn,:)'.^2);
    
    drift_post_mean(nn,:) = interp1(find(~isnan(fix_ids)),drift_post_mean(nn,~isnan(fix_ids)),1:NT);
    drift_post_std(nn,:) = interp1(find(~isnan(fix_ids)),drift_post_std(nn,~isnan(fix_ids)),1:NT);

    drift_post_mean_cor = squeeze(drift_post_mean(nn,:));
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
        fprintf('Refitting model for tr cell %d of %d\n',ss,length(tr_set));
        cur_tr_uset = tr_inds(~isnan(Robs_mat(tr_inds,ss)));
        
        tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
        
        dit_mods{nn+1}(ss) = dit_mods{nn}(ss);
        dit_mods{nn+1}(ss) = NMMfit_filters(dit_mods{nn+1}(ss),Robs_mat(cur_tr_uset,ss),tr_X,[],[],silent); %fit stimulus filters
        if strcmp(noise_model,'gauss')
            dit_mods_spkNL{nn+1}(ss) = dit_mods{nn+1}(ss);
        else
            dit_mods_spkNL{nn+1}(ss) = NMMfit_logexp_spkNL(dit_mods{nn+1}(ss),Robs_mat(cur_tr_uset,ss),tr_X);
        end
        newLL = NMMmodel_eval(dit_mods_spkNL{nn+1}(ss),Robs_mat(cur_tr_uset,ss),tr_X);
        dit_LLimp(nn+1,ss) = (newLL - null_LL(ss))/log(2);        
        fprintf('Original: %.5f  Prev: %.5f  New: %.5f\n',all_mod_LLimp(ss),dit_LLimp(nn,ss),dit_LLimp(nn+1,ss));
    end
end
%%
cd(anal_dir)
save(anal_name,'sim_eyepos*','it_*','dit_*','sim_gain','all_mod*','tr_set','drift_post*');

%%
n_fixs = length(fix_start_inds);
fix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
end

fin_fix_corr = nan(NT,1);
fin_fix_std = nan(NT,1);
fin_fix_corr(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
fin_fix_std(~isnan(fix_ids)) = it_fix_post_std(end,fix_ids(~isnan(fix_ids)));
fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);

fin_fix_corr = fin_fix_corr*sp_dx;
fin_fix_std = fin_fix_std*sp_dx;

fin_drift_corr = drift_post_mean(end,:)*sp_dx;
fin_drift_std = drift_post_std(end,:)*sp_dx;
all_fix_inds = [];
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > sac_shift
        fin_drift_corr(cur_inds(1:end-sac_shift+1)) = fin_drift_corr(cur_inds(sac_shift:end));
        fin_drift_std(cur_inds(1:end-sac_shift+1)) = fin_drift_std(cur_inds(sac_shift:end));
    end
    all_fix_inds = [all_fix_inds cur_inds];
end

fin_tot_corr = fin_fix_corr + fin_drift_corr;
fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2);

usable_inds = find(~ismember(all_blockvec(used_inds),ignore_blocks));

fin_tot_err = fin_tot_corr - sim_eyepos(used_inds)';
usable_fix_inds = all_fix_inds(ismember(all_fix_inds,usable_inds));
fin_sure_inds = usable_fix_inds(fin_tot_std(usable_fix_inds) <= 0.05);

%%
figure; hold on
xx = linspace(0,1,500);
fin_nn = hist(abs(fin_tot_err(usable_fix_inds)),xx);
fin_nn = cumsum(fin_nn)/sum(fin_nn);
plot(xx,1-fin_nn,'r')
fin_nn = hist(abs(fin_tot_err(fin_sure_inds)),xx);
fin_nn = cumsum(fin_nn)/sum(fin_nn);
plot(xx,1-fin_nn,'k')

err_mad = median(abs(fin_tot_err(usable_fix_inds)));
%%

% measured_seqL = corrected_eye_vals_interp(used_inds,2);
measured_seqL = sim_eyepos(used_inds);
measured_seqR = corrected_eye_vals_interp(used_inds,4);

min_fix_dur = 0.15;
inferred_drift = nan(size(fin_tot_corr));
measured_driftL = nan(size(fin_tot_corr));
measured_driftR = nan(size(fin_tot_corr));
inferred_fix_avg = nan(n_fixs,1);
measured_fix_avgL = nan(n_fixs,1);
measured_fix_avgR = nan(n_fixs,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds)*dt >= min_fix_dur & ismember(used_inds(cur_inds(1)),usable_inds)
        cur_inf = fin_tot_corr(cur_inds);
        inferred_fix_avg(ii) = median(fin_tot_corr(cur_inds));
%         inferred_fix_avg(ii) = mean(fin_tot_corr(cur_inds));
        inferred_drift(cur_inds) = cur_inf - inferred_fix_avg(ii);
        
        measured_fix_avgL(ii) = median(measured_seqL(cur_inds));
        measured_fix_avgR(ii) = median(measured_seqR(cur_inds));
%         measured_fix_avgL(ii) = mean(measured_seqL(cur_inds));
%         measured_fix_avgR(ii) = mean(measured_seqR(cur_inds));
        measured_driftL(cur_inds) = measured_seqL(cur_inds) - measured_fix_avgL(ii);
        measured_driftR(cur_inds) = measured_seqR(cur_inds) - measured_fix_avgR(ii);
     end
end

u = find(~isnan(measured_driftL) & ~isnan(inferred_drift));
[drift_corrs,drif_pvals] = corr([measured_driftL(u)' measured_driftR(u)'],inferred_drift(u)','type','spearman');
[drift_corrs_pear,drif_pvals] = corr([measured_driftL(u)' measured_driftR(u)'],inferred_drift(u)');
u = find(~isnan(measured_fix_avgL) & ~isnan(inferred_fix_avg));
[fix_corrs,fix_pvals] = corr([measured_fix_avgL(u) measured_fix_avgR(u)],inferred_fix_avg(u),'type','spearman');
[tot_corrs,tot_pvals] = corr([measured_seqL(usable_inds) measured_seqR(usable_inds)],fin_tot_corr(usable_inds)','type','spearman');
[tot_corrs_pear,tot_pvals_pear] = corr([measured_seqL(usable_inds) measured_seqR(usable_inds)],fin_tot_corr(usable_inds)','type','pearson');


%%


%%
figure
xx = linspace(-1,1,500);
lin_nn = hist(lin_tot_err(used_fix_inds),xx);
lin_nn = cumsum(lin_nn)/sum(lin_nn);
fin_nn = hist(fin_tot_err(used_fix_inds),xx);
fin_nn = cumsum(fin_nn)/sum(fin_nn);
gauss_nn = hist(gauss_tot_err(used_fix_inds),xx);
gauss_nn = cumsum(gauss_nn)/sum(gauss_nn);
plot(xx,lin_nn,xx,fin_nn,'r',xx,gauss_nn,'k')

figure
xx = linspace(-1,1,500);
lin_nn = hist(lin_tot_err(lin_sure_inds),xx);
lin_nn = cumsum(lin_nn)/sum(lin_nn);
fin_nn = hist(fin_tot_err(fin_sure_inds),xx);
fin_nn = cumsum(fin_nn)/sum(fin_nn);
plot(xx,lin_nn,xx,fin_nn,'r')

%%

%%
sMarkers = [cjump_inds(used_fixs) ejump_inds(used_fixs)];
params.Fs = 1/dt;
params.tapers = [3 5];
params.err = 0;
movingwin = [0.2 0.2];
[C,~,~,S,f] = coherencyc_unequal_length_trials([fin_tot_corr(:) sim_eyepos_us(used_inds)*sp_dx],movingwin,params,sMarkers);
[Clin,~,~,Slin,f] = coherencyc_unequal_length_trials([lin_tot_corr(:) sim_eyepos_us(used_inds)*sp_dx],movingwin,params,sMarkers);
[Cgauss,~,~,Sgauss,f] = coherencyc_unequal_length_trials([gauss_tot_corr(:) sim_eyepos_us(used_inds)*sp_dx],movingwin,params,sMarkers);

figure
plot(f,Clin,f,C,'r',f,Cgauss,'k')

%%
movingwin = [0.5 0.5];
sMarkers = [trial_start_inds(use_tt) trial_end_inds(use_tt)];
cur_trial_durs = (trial_end_inds(use_tt)-trial_start_inds(use_tt))*dt;
sMarkers(cur_trial_durs <= movingwin(1),:) = [];
[C_full,~,~,S,f_full] = coherencyc_unequal_length_trials([fin_tot_corr(:) sim_eyepos_us(used_inds)*sp_dx],movingwin,params,sMarkers);
[Clin_full,~,~,Slin,f_full] = coherencyc_unequal_length_trials([lin_tot_corr(:) sim_eyepos_us(used_inds)*sp_dx],movingwin,params,sMarkers);
[Cgauss_full,~,~,Sgauss,f_full] = coherencyc_unequal_length_trials([gauss_tot_corr(:) sim_eyepos_us(used_inds)*sp_dx],movingwin,params,sMarkers);

figure
plot(f_full,Clin_full,f_full,C_full,'r',f_full,Cgauss_full,'k')

%%
cd ~/Analysis/bruce/summary_analysis/eyetrack_figs/

close all
for tt = 1:n_trials
    uu = find(all_trialvec(used_inds) == tt);
    if ~isempty(uu)
        bt = all_t_axis(used_inds(uu(1)));
        et = all_t_axis(used_inds(uu(end)));
        dur = et-bt;
        if dur > 3.5
            hold on
%             h3=shadedErrorBar(all_t_axis(used_inds(uu))-bt,lin_fin_corr(uu),lin_fin_std(uu),{'color','b'});
            h3=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_tot_corr(uu),fin_tot_std(uu),{'color','r'});
                h3=plot(all_t_axis(used_inds(uu))-bt,sim_eyepos(used_inds(uu)),'k','linewidth',2);
            xlim([0 dur]);
            ylim([-0.75 0.75]);
%             ylim([-0.3 0.3]);
            xlabel('Time (s)','fontsize',10);
            ylabel('Orthoganol position (deg)','fontsize',10);
            title(sprintf('Trial %d',tt));
            set(gca,'fontsize',8,'fontname','arial');
            
            
            fillPage(gcf,'papersize',[8 8]);
            pause
            clf
        end
    end
end

% n_bins = 100;
% b_ax = linspace(-0.6,0.6,n_bins);
% nx = histc(measured_seqL(usable_inds)',b_ax);
% nx = nx/sum(nx);
% figure; hold on
% stairs(b_ax,nx,'b');
% yl = ylim();
% line([0 0],yl,'color','k')
% xlim(b_ax([1 end])); 
% box off;
% title('Inferred position','fontsize',10);
% fillPage(gcf,'papersize',[4 4]);
% fname = sprintf('Measured_orthdist_ori%d',bar_ori);

