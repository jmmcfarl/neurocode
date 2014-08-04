clear all

addpath('~/James_scripts/bruce/eye_tracking/');
addpath('~/James_scripts/bruce/processing/');

global Expt_name bar_ori use_LOOXV

Expt_name = 'G086';
% use_LOOXV = 0; %[0 no LOOXV; 1 SU LOOXV; 2 all LOOXV]

bar_ori = 0; %bar orientation to use (only for UA recs)

%%

Expt_num = str2num(Expt_name(2:end));
if Expt_name(1) == 'M'
    rec_type = 'LP';
elseif Expt_name(1) == 'G'
    rec_type = 'UA';
end

if strcmp(rec_type,'LP')
    switch Expt_num
        case 266
            bar_ori = 80;
        case 270
            bar_ori = 60;
        case 275
            bar_ori = 135;
        case 277
            bar_ori = 70;
        case 281
            bar_ori = 140;
        case 287
            bar_ori = 90;
        case 289
            bar_ori = 160;
        case 294
            bar_ori = 40;
    end
end

if Expt_num > 280 && Expt_num < 289
    data_dir = ['/media/NTlab_data2/Data/bruce/' Expt_name];
elseif Expt_num >= 289
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
else
    data_dir = ['~/Data/bruce/' Expt_name];
end

cd(data_dir);

if strcmp(rec_type,'LP')
    if Expt_num >= 275
        rpt_seed = 1001; %M275 M277 M281
    else
        rpt_seed = 1e4; %m270 and 266
    end
    load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct
    good_coils = [1 1]; %which coils are usable
    use_coils = [1 1]; %[L R] Use info from coils?
    n_probes = 24;
    use_nPix = 32;
    use_LOOXV = 2;
elseif strcmp(rec_type,'UA')
    rpt_seed = nan;
    load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
    good_coils = [1 0]; %which coils are usable
    use_coils = [0 0]; %[L R] Use info from coils?
    n_probes = 96;
    use_nPix = 16;
    use_LOOXV = 1;
end

load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

et_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];

et_mod_data_name = 'full_eyetrack_initmods';
et_anal_name = 'full_eyetrack';

%%
%dont fit stim models using these blocks
ignore_blocks = [];
switch Expt_num
    case 270
        ignore_blocks = [5 19];
    case 289
        ignore_blocks = [27]; %27 is off somehow
    case 294
        ignore_blocks = [37 38 39]; %37-39 have slightly different dw used in these blocks
    case 86
        ignore_blocks = [16 17 28 30];
    case 87
        ignore_blocks = [15];
    case 93
        ignore_blocks = [28];
end

%problem with M270 where vd was wrong, need this correction factor to get
%correct units
if Expt_num==270
    scale_fac = 1.72;
else
    scale_fac = 1;
end

%if using coil info
if any(use_coils > 0)
    et_anal_name = [et_anal_name '_Cprior'];
end

%%
if ~isnan(rpt_seed)
    xv_type = 'rpt';
    xv_frac = nan;
else
    xv_type = 'uni';
    xv_frac = 0.2;
end

flen = 15;
spatial_usfac = 2;

%these recs have larger bar widths
if ismember(Expt_num,[287 289 294])
    use_nPix = 15;
    spatial_usfac = 4;
end

min_trial_dur = 0.75;

stim_fs = 100; %in Hz
dt = 0.01;
Fr = 1;

full_nPix=36;
switch Expt_num
    case 270
        full_nPix=32;
    case  287
        full_nPix = 22;
    case 289
        full_nPix = 22;
    case 294
        full_nPix = 20;
end

%exclude data at beginning and end of each trial
beg_buffer = 0.2;
end_buffer = 0.05;
trial_dur = 4;

use_nPix_us = use_nPix*spatial_usfac;
klen_us = use_nPix_us*flen;

%saccade kernel time axis
sac_backlag = round(0.05/dt);
sac_forlag = round(0.3/dt);
sac_bincents = -sac_backlag:sac_forlag;
n_sac_bins = length(sac_bincents);

stim_params = NIMcreate_stim_params([flen full_nPix],dt);

%% SELECT BLOCKS FOR ANALYSIS
include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaXFaXFs','rls.AllSac','rls.imiXFa','rls.FaXwi','rls.FaXwiXimi'};
expt_names = cell(1,length(Expts));
expt_dds = nan(1,length(Expts));
expt_bar_ori = nan(1,length(Expts));
expt_sac_dir = nan(1,length(Expts));
expt_Fr = nan(1,length(Expts));
expt_imback = nan(1,length(Expts));
included_type = false(1,length(Expts));
for ii = 1:length(Expts)
    if ~isempty(Expts{ii})
        expt_names{ii} = Expts{ii}.Header.expname;
        expt_dds(ii) = Expts{ii}.Stimvals.dd;
        expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
        expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
        expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
        expt_imback(ii) = isfield(Expts{ii}.Trials,'imi');
        included_type(ii) = any(strcmp(expt_names{ii},include_expts));
    end
end
expt_has_ds(isnan(expt_has_ds)) = 0;
expt_has_ds(expt_has_ds == -1) = 0;
expt_binoc(isnan(expt_binoc)) = 0;

if strcmp(rec_type,'LP')
    expt_bar_ori(expt_bar_ori > 360) = bar_ori;
end

cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori);

cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];

sim_sac_expts = find(~expt_has_ds(cur_block_set));
imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');

n_blocks = length(cur_block_set);

all_nfs = cellfun(@(x) x.Stimvals.nf,Expts(cur_block_set));
if length(unique(all_nfs)) > 1
    fprintf('Warning, multiple different nfs detected: %.4f\n',all_nfs);
end
%%
all_dws = cellfun(@(x) x.Stimvals.dw,Expts(cur_block_set));
base_sp_dx = mode(all_dws);
if length(unique(all_dws)) > 1
    fprintf('Warning, multiple different dws detected, using %.3f\n',base_sp_dx);
end
sp_dx = base_sp_dx/spatial_usfac/scale_fac; %model dx in deg

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
all_trial_Se = [];
all_trial_wi = [];
all_trial_blk = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
trial_toffset = zeros(length(cur_block_set),1);
cur_toffset = 0;
cur_spkind_offset = 0;
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
    
    fname = [cluster_dir sprintf('/Block%d_Clusters.mat',cur_block)];
    load(fname,'Clusters');
    for cc = 1:n_probes
        all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times + cur_toffset);
        all_spk_inds{cc} = cat(1,all_spk_inds{cc},Clusters{cc}.spk_inds + cur_spkind_offset);
        all_clust_ids{cc} = cat(1,all_clust_ids{cc},Clusters{cc}.spike_clusts);
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
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    trial_Se = [Expts{cur_block}.Trials(:).se];
    trial_Se = trial_Se(id_inds);
    all_trial_Se = cat(1,all_trial_Se,trial_Se(use_trials)');
    all_trial_blk = cat(1,all_trial_blk,ones(length(use_trials),1)*ee);
    
    if isfield(Expts{cur_block}.Trials,'wi')
        trial_wi = [Expts{cur_block}.Trials(:).wi];
        trial_wi = trial_wi(id_inds);
        all_trial_wi = cat(1,all_trial_wi,trial_wi(use_trials)');
    else
        all_trial_wi = cat(1,all_trial_wi,nan(length(use_trials),1));
    end
    
    fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
    load(fname);
    buffer_pix = floor((expt_npix(cur_block) - full_nPix)/2);
    if buffer_pix == -1
        for ii = 1:length(left_stim_mats)
            left_stim_mats{ii} = [zeros(size(left_stim_mats{ii},1),1) left_stim_mats{ii} zeros(size(left_stim_mats{ii},1),1)];
        end
        buffer_pix = 0;
    end
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
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times+cur_toffset < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times + cur_toffset];
            all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
            all_stim_mat = [all_stim_mat; cur_stim_mat];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
    end
    trial_cnt = trial_cnt + n_trials;
    
    %need to keep track of block time offsets for LP recordings
    if strcmp(rec_type,'LP')
        trial_toffset(ee) = all_t_bin_edges(end);
        cur_toffset = trial_toffset(ee);
        cur_spkind_offset = cur_spkind_offset+ round(32e3*(max(trial_end_times)-min(trial_start_times) + 50));
    end
end

%%
full_nPix_us = spatial_usfac*full_nPix;
if spatial_usfac > 1
    all_stimmat_up = zeros(size(all_stim_mat,1),full_nPix_us);
    for ii = 1:size(all_stim_mat,2)
        for jj = 1:spatial_usfac
            all_stimmat_up(:,spatial_usfac*(ii-1)+jj) = all_stim_mat(:,ii);
        end
    end
elseif spatial_usfac == 1
    all_stimmat_up = all_stim_mat;
end
stim_params_us = NMMcreate_stim_params([flen full_nPix_us],dt);

%% select submatrix with central pixels

%pick out the predictors (of the time-embedded stimulus-matrices)
%corresponding to the central pixels we're using in the models
[Xinds,~] = meshgrid(1:full_nPix,1:flen);
buffer_pix = floor((full_nPix - use_nPix)/2);
cur_use_pix = (1:use_nPix) + buffer_pix;
use_kInds = find(ismember(Xinds(:),cur_use_pix)); %use this set of dimensions in the Xmatrix

%repeat for up-sampled versions of the Xmatrix
[Xinds_up,~] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
cur_use_pix = (1/spatial_usfac:1/spatial_usfac:use_nPix) + buffer_pix;
use_kInds_up = find(ismember(Xinds_up(:),cur_use_pix));

%% BIN SPIKES FOR MU AND SU
clust_params.n_probes = n_probes;
if strcmp(rec_type,'LP')
    clust_params.exclude_adjacent = true;
else
    clust_params.exclude_adjacent = false;
end
[all_binned_mua,all_binned_sua,Clust_data] = ...
    get_binned_spikes(cluster_dir,all_spk_times,all_clust_ids,all_spk_inds,...
    all_t_axis,all_t_bin_edges,all_bin_edge_pts,cur_block_set,all_blockvec,clust_params);
su_probes = Clust_data.su_probes;
SU_numbers = Clust_data.SU_numbers;

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
% [all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data_v2(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset,good_coils);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

%compute corrected eye data in bar-oriented frame
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,expt_bar_ori(cur_block_set),used_inds);

[saccades,et_params] = detect_saccades_v2(corrected_eye_vals,all_eye_vals,all_eye_speed,all_eye_ts,et_params);

is_blink = detect_blinks(all_eye_ts,all_eye_vals,saccades,et_params);

[saccades,is_blink] = merge_blinks(saccades,is_blink);

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
is_blink(bad_sacs) = [];
interp_sac_start_inds(bad_sacs) = [];
interp_sac_stop_inds(bad_sacs) = [];

%% CREATE SACCADE PREDICTOR MATS
saccade_start_inds = find(ismember(used_inds,interp_sac_start_inds));
used_saccade_set = find(ismember(interp_sac_start_inds,used_inds));
%nearest index in the used data set of the saccade stop time
saccade_stop_inds = round(interp1(used_inds,1:length(used_inds),interp_sac_stop_inds(used_saccade_set)))';

saccade_stop_inds(isnan(saccade_stop_inds)) = length(used_inds);
used_is_blink = is_blink(used_saccade_set);

sac_amps = [saccades(:).amplitude];
is_micro = sac_amps(used_saccade_set) < 1;
big_sacs = find(~is_micro & ~used_is_blink');
micro_sacs = find(is_micro & ~used_is_blink');
sac_durs = [saccades(:).duration];

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

%% DEFINE FIXATION POINTS
trial_start_inds = [1; find(diff(all_trialvec(used_inds)) ~= 0) + 1];
trial_end_inds = [find(diff(all_trialvec(used_inds)) ~= 0); NT];

fix_start_inds = sort([trial_start_inds; saccade_stop_inds]);
fix_stop_inds = sort([trial_end_inds; saccade_start_inds]);
fix_durs = fix_stop_inds-fix_start_inds;
fix_start_inds(fix_durs<=0) = [];
fix_stop_inds(fix_durs <= 0) = [];
fix_post_blink = ismember(fix_start_inds,saccade_stop_inds(used_is_blink));
n_fixs = length(fix_start_inds);

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
pfix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
    cur_inds = pfix_start_inds(ii):(pfix_stop_inds(ii));
    pfix_ids(cur_inds) = ii;
end

%push the effects of saccades forward in time
trial_ids = nan(NT,1);
for ii = 1:n_trials
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    trial_ids(cur_inds) = ii;
end


%%
cd(et_dir)
load(et_mod_data_name,'all_mod*');
load(et_anal_name,'drift*','it_*','et_tr_set');
tr_set = et_tr_set;
%%
n_chs = size(all_binned_mua,2) + size(all_binned_sua,2);

% make Robs_mat
poss_blocks = 1:(n_blocks-1);
Robs_mat = nan(length(used_inds),n_chs);
for ss = 1:n_chs
    if all_mod_SU(ss) > 0
        su_probe_ind = find(Clust_data.SU_numbers == all_mod_SUnum(ss));
        Robs_mat(:,ss) = all_binned_sua(used_inds,su_probe_ind);
    else
        Robs_mat(:,ss) = all_binned_mua(used_inds,ss);
    end
end

%% Create set of TR and XV trials
use_trials = unique(all_trialvec(used_inds));
nuse_trials = length(use_trials);

if strcmp(xv_type,'rpt')
    xv_trials = find(all_trial_Se==rpt_seed);
    n_xv_trials = length(xv_trials);
else
    n_xv_trials = round(xv_frac*nuse_trials);
    xv_trials = randperm(nuse_trials);
    xv_trials(n_xv_trials+1:end) = [];
    xv_trials = use_trials(xv_trials);
end
tr_trials = setdiff(use_trials,xv_trials);
n_tr_trials = length(tr_trials);
fprintf('Initializing models with %d training trials and %d xval trials\n',n_tr_trials,n_xv_trials);

tr_inds = find(ismember(all_trialvec(used_inds),tr_trials));
xv_inds = find(ismember(all_trialvec(used_inds),xv_trials));

%%
if use_LOOXV == 2
    loo_set = 1:length(tr_set);
elseif use_LOOXV == 1
    loo_set = find(all_mod_SU(tr_set) > 0);
else
    loo_set = [];
end

%%
backlag = round(0.05/dt);
forlag = round(0.2/dt);
slags = -backlag:forlag;
n_sac_bins = length(slags);

Xsac = zeros(NT,length(slags));
Xmsac = zeros(NT,length(slags));
for ii = 1:n_sac_bins
    cur_sac_target = saccade_start_inds(big_sacs) + slags(ii);
    uu = find(cur_sac_target > 1 & cur_sac_target < NT);
    cur_sac_target = cur_sac_target(uu);
    cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_trial_inds(big_sacs(uu))) = [];
    Xsac(cur_sac_target,ii) = 1;
    
    cur_sac_target = saccade_start_inds(micro_sacs) + slags(ii);
    uu = find(cur_sac_target > 1 & cur_sac_target < NT);
    cur_sac_target = cur_sac_target(uu);
    cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_trial_inds(micro_sacs(uu))) = [];
    Xmsac(cur_sac_target,ii) = 1;
end

%% RECONSTRUCT RETINAL STIM FOR DESIRED UNIT
fprintf('Reconstructing retinal stim for unit %d\n',cc);
fin_fix_corr = nan(NT,1);
fin_fix_std = nan(NT,1);
fin_fix_corr(~isnan(fix_ids)) = squeeze(it_fix_post_mean(end,fix_ids(~isnan(fix_ids))));
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
fin_fix_std(~isnan(fix_ids)) = squeeze(it_fix_post_std(end,fix_ids(~isnan(fix_ids))));
fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);

fin_drift_corr = squeeze(drift_post_mean(end,:));
fin_drift_std = squeeze(drift_post_std(end,:));

for ii = 1:length(trial_start_inds)
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
    fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
end
fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT);
fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);

fin_tot_corr = fin_fix_corr + fin_drift_corr;
fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2);

fin_shift_cor = round(fin_tot_corr);
fin_shift_cor(isnan(fin_shift_cor)) = 0;

%RECOMPUTE XMAT
all_shift_stimmat_up = all_stimmat_up;
for i=1:NT
    %     all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_post_cor(i)};
    all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
end
all_Xmat_shift = create_time_embedding(all_shift_stimmat_up,stim_params_us);
all_Xmat_shift = all_Xmat_shift(used_inds,use_kInds_up);

%%
for cc = [97 98 99 100 101 102 103]
    % for cc = [101 102 103]
    tr_cc = find(tr_set(loo_set) == cc);
    cc_uinds = find(~isnan(Robs_mat(:,cc)));
    
    
    %%
    fprintf('Fitting stim model for unit %d\n',cc);
    
    if expt_dds(cur_block_set(1)) == 12
        base_lambda_d2XT = 50;
        base_lambda_L1 = 5;
        base_qlambda_d2XT = 20;
        base_qlambda_L1 = 2;
    else
        base_lambda_d2XT = 200;
        base_lambda_L1 = 20;
        base_qlambda_d2XT = 80;
        base_qlambda_L1 = 8;
    end
    init_reg_params = NMMcreate_reg_params('lambda_d2XT',5);
    fin_stim_params = NMMcreate_stim_params([flen use_nPix_us],dt);
    
    silent = 0;
    
    cur_Robs = Robs_mat(cc_uinds,cc);
    
    stim_mod_signs = [1 1 1 -1 -1];
    stim_NL_types = {'lin','quad','quad','quad','quad'};
    qfilts = find(strcmp('quad',stim_NL_types));
    nqfilts = find(~strcmp('quad',stim_NL_types));
    init_mod = NMMinitialize_model( fin_stim_params, stim_mod_signs, stim_NL_types, init_reg_params);
    Xtargs = [init_mod.mods(:).Xtarget];
    
    init_fitN = ceil(length(cc_uinds)/10);
    init_fit_subset = randperm(length(cc_uinds));
    init_fit_subset = init_fit_subset(1:init_fitN);
    
    init_mod = NMMfit_filters(init_mod,cur_Robs(init_fit_subset),all_Xmat_shift(cc_uinds(init_fit_subset),:),[],[],silent);
    [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(init_mod,cur_Robs,all_Xmat_shift(cc_uinds,:));
    init_mod = NMMadjust_regularization(init_mod,nqfilts,'lambda_d2XT',base_lambda_d2XT./var(gint(:,nqfilts))');
    init_mod = NMMadjust_regularization(init_mod,nqfilts,'lambda_L1',base_lambda_L1./std(gint(:,nqfilts))');
    if ~isempty(qfilts)
        init_mod = NMMadjust_regularization(init_mod,qfilts,'lambda_d2XT',base_qlambda_d2XT./var(gint(:,qfilts))');
        init_mod = NMMadjust_regularization(init_mod,qfilts,'lambda_L1',base_qlambda_L1./std(gint(:,qfilts))');
    end
    cor_gqm = NMMfit_filters(init_mod,cur_Robs,all_Xmat_shift(cc_uinds,:),[],[],silent);
    all_stim_filts = reshape([cor_gqm.mods(:).filtK],flen,use_nPix_us,length(stim_mod_signs));
    all_temp_kerns{cc} = squeeze(nanstd(all_stim_filts,[],2));
    all_gqm_mods(cc) = cor_gqm;
    
end

%%
all_Xmat_us = create_time_embedding(all_stimmat_up,stim_params_us);
all_Xmat_us = all_Xmat_us(used_inds,:);

max_shift = round(0.8475/sp_dx);
dshift = 1;
shifts = -max_shift:dshift:max_shift;
n_shifts = length(shifts);
zero_frame = find(shifts==0);

%generate shift matrices. Must be applied to the stimulus (not the filters)
It = speye(flen);
shift_mat = cell(n_shifts,1);
fshift_mat = cell(n_shifts,1);
for xx = 1:n_shifts
    temp = spdiags( ones(full_nPix_us,1), -shifts(xx), full_nPix_us, full_nPix_us);
    temp = kron(temp,It);
    shift_mat{xx} = temp(:,use_kInds_up);
    fshift_mat{xx} = temp;
end
%%
T = tabulate(all_trialvec(used_inds));
comp_trials = find(T(:,2) == 375);

test_stim_trial = comp_trials(ceil(rand*length(comp_trials)));

comp_trials = comp_trials(randperm(length(comp_trials)));
test_stim_inds = find(all_trialvec(used_inds) == test_stim_trial);
test_stim = all_Xmat_us(test_stim_inds,:);

N_rpts = 100;
rpt_stim = repmat(test_stim,[N_rpts 1]);
rpt_et_trials = comp_trials(1:N_rpts);
rpt_eye_pos = nan(375,N_rpts);
for ii = 1:N_rpts
    rpt_eye_pos(:,ii) = fin_tot_corr(all_trialvec(used_inds) == rpt_et_trials(ii));
end
vec_rpt_eye_pos = rpt_eye_pos(:);
vec_rpt_eye_pos = round(vec_rpt_eye_pos);
vec_rpt_eye_pos(vec_rpt_eye_pos > max_shift) = max_shift;
vec_rpt_eye_pos(vec_rpt_eye_pos < -max_shift) = -max_shift;

vec_rpt_eye_pos = vec_rpt_eye_pos + max_shift + 1;
ep_test_stim = nan(N_rpts*375,length(use_kInds_up));
for ii = 1:N_rpts*375
    ep_test_stim(ii,:) = rpt_stim(ii,:)*shift_mat{vec_rpt_eye_pos(ii)};
end
avg_ep = nanmean(rpt_eye_pos,2);
%%
for cc = [97 98 99 100 101 102 103]

rate_noise_frac = 0;
rate_sc = 1;
trial_noise_fac = 0;

[LL, penLL, pred_rate,G] = NMMmodel_eval(all_gqm_mods(cc),[],ep_test_stim);
pred_rate = pred_rate*rate_sc;
avg_rate = mean(pred_rate);

stim_rate_SD = std(pred_rate);

rate_noise_SD = stim_rate_SD*rate_noise_frac;

tbt_pred_rate = reshape(pred_rate,[375 N_rpts])';
tbt_obs_rate = tbt_pred_rate + normrnd(zeros(size(tbt_pred_rate)),1)*rate_noise_SD;
trial_noise_SD = trial_noise_fac*stim_rate_SD;
trial_noise_rate = randn(N_rpts,1)*trial_noise_SD;
tbt_obs_rate = bsxfun(@plus,tbt_obs_rate,trial_noise_rate);

tbt_obs_rate(tbt_obs_rate < 0) = 0;

em_var = nanmean(nanvar(tbt_pred_rate));
psth_var = nanvar(nanmean(tbt_pred_rate));
true_alpha = psth_var/(em_var+psth_var);
tot_var = nanvar(tbt_obs_rate(:));

%%
tbt_spikes = poissrnd(tbt_obs_rate);
vec_spikes = reshape(tbt_spikes',[],1);
obs_psth = nanmean(tbt_spikes);

%% estimate spline models
% use_eye_range = [-0.4 0.4];
% n_pts = 500;
% eye_xx = linspace(use_eye_range(1),use_eye_range(2),n_pts);
% n_knots = 25;
% knot_pts = linspace(-0.4,0.4,n_knots);
% all_trials = 1:N_rpts;
% xv_frac = 0.2;
% xv_trials = randperm(length(all_trials));
% xv_trials(round(xv_frac*N_rpts)+1:end) = [];
% xv_trials = all_trials(xv_trials);
% tr_trials = setdiff(all_trials,xv_trials);
% [RR,TT] = meshgrid(1:375,1:N_rpts);
% TTvec = reshape(TT',[],1);
% 
% poss_lambdas = [1 5 10 50 500 1000];
% for ll = 1:length(poss_lambdas)
%     ll
%     lambda = poss_lambdas(ll);
%     
%     smoothed_eye_obs_fun = nan(length(eye_xx),375);
%     np_prate = nan(N_rpts,375);
%     sp = NMMcreate_stim_params([1 n_knots]);
%     spr = NMMcreate_reg_params('lambda_d2X',lambda,'boundary_conds',[Inf Inf Inf]);
%     
%     for ii = 1:375
% %         ii
%         xx =rpt_eye_pos(ii,:)*sp_dx;
%         xx(xx > use_eye_range(2) | xx < use_eye_range(1)) = nan;
%         yy_obs = squeeze(tbt_spikes(:,ii));
%         
%         uu = all_trials(~isnan(xx(all_trials)) & ~isnan(yy_obs(all_trials)'));
%         %         knot_pts = prctile(xx(uu),linspace(0,100,n_knots));
%         
%         uu_tr_set = ismember(TT(uu,ii),tr_trials);
%         uu_tr = uu(ismember(TT(uu,ii),tr_trials));
%         uu_xv = uu(ismember(TT(uu,ii),xv_trials));
%         if sum(yy_obs(uu_tr)>0) < 3
%             np_prate(uu,ii) = nan(size(uu));
%             smoothed_eye_obs_fun(:,ii) = 0;
%         else
%             temp = fnval(csape(knot_pts,eye(length(knot_pts)),'var'),xx(uu));
%             mod = NMMinitialize_model(sp,[1],{'lin'},spr,[],[],'logexp');
%             mod = NMMfit_filters(mod,yy_obs(uu_tr),temp(:,uu_tr_set)');
%             mod_spk = NMMfit_logexp_spkNL(mod,yy_obs(uu),temp',[],[],[3]);
%             
%             filtK = mod_spk.mods(1).filtK;
%             pfun = csape(knot_pts,filtK,xx(uu),'var');
%             fit_spk_nl = mod_spk.spk_NL_params;
%             smoothed_eye_obs_fun(:,ii) = fit_spk_nl(3)*log(1+exp(fit_spk_nl(2)*(fit_spk_nl(1) + fnval(pfun,eye_xx))));
%             
%             beta = fit_spk_nl(2)*(fit_spk_nl(1)+fnval(pfun,xx(uu)));
%             too_big = beta > 50;
%             cur_eyepos_p = fit_spk_nl(3)*log(1+exp(beta));
%             cur_eyepos_p(too_big) = fit_spk_nl(3)*beta(too_big);
%             np_prate(uu,ii) = cur_eyepos_p;
%             
%             out_bounds = find(eye_xx > max(xx(uu)) | eye_xx < min(xx(uu)));
%             smoothed_eye_obs_fun(out_bounds,ii) = nan;
%             
%             udata = find(ismember(TTvec,xv_trials));
%             vec_np_prate = reshape(np_prate',[],1);
%             spline_xvLL(ll) = nansum(vec_spikes(udata).*log(vec_np_prate(udata)) - vec_np_prate(udata))./nansum(vec_spikes(udata));
%             
%         end
%     end
%     
% end

%%
ep_noise = 0.03;

n_tests = 100;

clear shift_em_var shift_psth_var* shift_alpha*
for tt = 1:n_tests
    tbt_spikes = poissrnd(tbt_obs_rate);
    tbt_spikes = bsxfun(@minus,tbt_spikes,mean(tbt_spikes,2));
%     tbt_spikes = bsxfun(@rdivide,tbt_spikes,mean(tbt_spikes,2));
%     tbt_spikes(isnan(tbt_spikes)) = 0;
    vec_spikes = reshape(tbt_spikes',[],1);
    obs_psth = nanmean(tbt_spikes);
    
    
    [II,JJ] = meshgrid(1:N_rpts);
    Tinds = [II(:) JJ(:)];
    
    use_neur_eye_pos = rpt_eye_pos*sp_dx;
    use_neur_eye_pos = use_neur_eye_pos + randn(size(rpt_eye_pos))*ep_noise;
    
%     eye_noise = trnd(3,[size(rpt_eye_pos,1) size(rpt_eye_pos,2)]);
%     eye_noise = eye_noise/std(eye_noise(:))*ep_noise;
%     use_neur_eye_pos = use_neur_eye_pos + eye_noise;
%     [eye_b,eye_a] = butter(4,10/(1/(2*dt)),'low');
%     eye_noise = randn(375*N_rpts,1);
%     eye_noise = filtfilt(eye_b,eye_a,eye_noise);
%     eye_noise = eye_noise/std(eye_noise)*ep_noise;
%     eye_noise = reshape(eye_noise,375,N_rpts);
%     use_neur_eye_pos = use_neur_eye_pos + eye_noise;
    
    [a,b] = sort(use_neur_eye_pos,2);
    
    rand_eyepos = randn(size(use_neur_eye_pos));
    rand_eyepos(isnan(use_neur_eye_pos)) = nan;
    [arand,brand] = sort(rand_eyepos,2);
    
    not_usable = isnan(use_neur_eye_pos');
    
    b = b';a = a';
    arand = arand'; brand = brand';
    btshuff = b;
    
    max_adiff = sp_dx;
    adiff = [zeros(1,375); diff(a)];
    
    shift_size = 1;
    
    full_psth_eyeperm1 = tbt_spikes;
    full_psth_eyeperm2 = tbt_spikes;
    full_psth_rperm1 = tbt_spikes;
    full_psth_rperm2 = tbt_spikes;
    
    for ii = 1:375
        temp = squeeze(tbt_spikes(:,ii));
        tempr = temp(b(:,ii));
        full_psth_eyeperm1(:,ii) = tempr;
        
        temp = squeeze(tbt_spikes(:,ii));
        tempr(1:end-shift_size,:) = temp(b((1+shift_size):end,ii),:);
        tempr((end-shift_size+1):end,:) = nan;
        tempr(adiff(2:end,ii) > max_adiff,:) = nan;
        full_psth_eyeperm2(:,ii) = tempr;
        
        temp = squeeze(tbt_spikes(:,ii));
        tempr = temp(brand(:,ii),:);
        full_psth_rperm1(:,ii) = tempr;
        
        temp = squeeze(tbt_spikes(:,ii));
        tempr(1:end-1,:) = temp(brand(2:end,ii),:);
        tempr(adiff(2:end,ii) > max_adiff,:) = nan;
        full_psth_rperm2(:,ii) = tempr;
    end
    
    vec_obs_eyeperm1 = reshape(full_psth_eyeperm1',[],1);
    vec_obs_eyeperm2 = reshape(full_psth_eyeperm2',[],1);
    vec_obs_rperm1 = reshape(full_psth_rperm1',[],1);
    vec_obs_rperm2 = reshape(full_psth_rperm2',[],1);
    
    shift_em_var(tt) = nanmean(vec_obs_eyeperm1.*vec_obs_eyeperm2) - nanmean(vec_obs_eyeperm1).*nanmean(vec_obs_eyeperm2);
    shift_psth_var(tt) = nanmean(vec_obs_rperm1.*vec_obs_rperm2) - nanmean(vec_obs_rperm1).*nanmean(vec_obs_rperm2);
    shift_psth_var2(tt) = var(obs_psth);
    shift_alpha(tt) = shift_psth_var(tt)/shift_em_var(tt);
    shift_alpha2(tt) = var(obs_psth)/shift_em_var(tt);
    
    psth_noise = mean(var(tbt_spikes,[],2))/N_rpts;
    shift_alpha3(tt) = (shift_psth_var2(tt)-psth_noise)/shift_em_var(tt);
       
end

fprintf('Est alpha: Mean: %.4f SD: %.4f\n',mean(shift_alpha3),std(shift_alpha3));
fprintf('True alpha: %.4f\n',true_alpha);
avg_shift_alpha(cc) = mean(shift_alpha3);
all_true_alphas(cc) = true_alpha;
std_shift_alpha(cc) = std(shift_alpha3);
all_psth_var(cc) = psth_var;
all_em_var(cc) = em_var;
avg_shift_em_var(cc) = nanmean(shift_em_var);
all_avg_rates(cc) = avg_rate;
end
%%
srange = find(shifts*sp_dx >= use_eye_range(1) & shifts*sp_dx <= use_eye_range(2));
eyefun_pred_out = nan(length(srange),375);
for ss = 1:length(srange)
    fprintf('%d of %d\n',ss,length(srange));
    shift_data = test_stim*shift_mat{srange(ss)};
    [~, ~, temp_prate] = NMMmodel_eval(all_gqm_mods(cc),[],shift_data);
    eyefun_pred_out(ss,:) = temp_prate;
end

eyefun_pred_out_int = interp1(-shifts(srange)*sp_dx,eyefun_pred_out,eye_xx');












%%
T = tabulate(all_trialvec(used_inds));
comp_trials = find(T(:,2) == 375);

test_stim_trial = comp_trials(ceil(rand*length(comp_trials)));

test_stim_inds = find(all_trialvec(used_inds) == test_stim_trial);
test_stim = all_Xmat_us(test_stim_inds,:);

cc = 98;
rate_noise_frac = 0;
rate_sc = 1;
N_rpts = 50;
poss_shifts = linspace(-0.4,0.4,20);
clear obs_psth
for pp = 1:length(poss_shifts)
    fprintf('Shift %d of %d\n',pp,length(poss_shifts));
    
    [~,cur_shift] = min(abs(shifts*sp_dx - poss_shifts(pp)));
    shift_test_stim = test_stim*fshift_mat{cur_shift};
    
    comp_trials = comp_trials(randperm(length(comp_trials)));
    rpt_stim = repmat(shift_test_stim,[N_rpts 1]);
    rpt_et_trials = comp_trials(1:N_rpts);
    rpt_eye_pos = nan(375,N_rpts);
    for ii = 1:N_rpts
        rpt_eye_pos(:,ii) = fin_tot_corr(all_trialvec(used_inds) == rpt_et_trials(ii));
    end
    vec_rpt_eye_pos = rpt_eye_pos(:);
    vec_rpt_eye_pos = round(vec_rpt_eye_pos);
    vec_rpt_eye_pos(vec_rpt_eye_pos > max_shift) = max_shift;
    vec_rpt_eye_pos(vec_rpt_eye_pos < -max_shift) = -max_shift;
    
    vec_rpt_eye_pos = vec_rpt_eye_pos + max_shift + 1;
    ep_test_stim = nan(N_rpts*375,length(use_kInds_up));
    for ii = 1:N_rpts*375
        ep_test_stim(ii,:) = rpt_stim(ii,:)*shift_mat{vec_rpt_eye_pos(ii)};
    end
    
    [LL, penLL, pred_rate,G] = NMMmodel_eval(all_gqm_mods(cc),[],ep_test_stim);
    pred_rate = pred_rate*rate_sc;
    avg_rate = mean(pred_rate);
    
    stim_rate_SD = std(pred_rate);
    rate_noise_SD = stim_rate_SD*rate_noise_frac;
    tbt_pred_rate = reshape(pred_rate,[375 N_rpts])';
    tbt_obs_rate = tbt_pred_rate + normrnd(zeros(size(tbt_pred_rate)),1)*rate_noise_SD;
    tbt_obs_rate(tbt_obs_rate < 0) = 0;
    
    tbt_spikes = poissrnd(tbt_obs_rate);
    obs_psth(pp,:) = mean(tbt_spikes);
    
end

%%
lags = -length(poss_shifts)/2:length(poss_shifts)/2;
spacing = median(diff(poss_shifts));

obs_corr = corr(obs_psth');
clear temp
for ii = 1:length(lags)
    temp(ii) = mean(diag(obs_corr,lags(ii)));
end
lags = lags*spacing;

% temp(lags == 0) = nan;
temp = temp/nansum(temp);
temp_SD = sqrt(nansum(temp.*lags.^2));
eye_SD = nanstd(fin_tot_corr*sp_dx);
eye_SD2 = robust_std_dev(fin_tot_corr*sp_dx);
eye_SD = 0.5*eye_SD + 0.5*eye_SD2;
pe = normpdf(lags,0,eye_SD);

%%
eyefun_pred_out_int(isnan(eyefun_pred_out_int)) = 0;

lags2 = -200:200;
ecorr = corr(eyefun_pred_out_int');
clear temp2
for ii = 1:length(lags2)
    temp2(ii) = nanmean(diag(ecorr,lags2(ii)));
end
spacing = median(diff(eye_xx));
lags2 = lags2*spacing;

temp2 = temp2/sum(temp2);
temp2_SD = sqrt(sum(temp2.*lags2.^2));

%%
figure
plot(lags,pe/max(pe),'k')
hold on
plot(lags,temp/max(temp),'r')
plot(lags2,temp2/max(temp2),'b')

est_SD = sqrt(temp_SD.^2-eye_SD.^2)/sqrt(2);
t1 = normpdf(lags2,0,est_SD);
plot(lags2,t1/max(t1),'c')
t2 = normpdf(lags,0,temp_SD);
plot(lags,t2/max(t2),'r--')



%%
freqs = linspace(-20,20,100);
true_pspec = sqrt(normpdf(freqs,0,1/(2*pi*temp2_SD)));
eye_pspec = normpdf(freqs,0,1/(2*pi*eye_SD));
eye_pspec = eye_pspec/trapz(freqs,eye_pspec);
