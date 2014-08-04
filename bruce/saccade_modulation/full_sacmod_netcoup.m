clear all

addpath('~/James_scripts/bruce/eye_tracking/');
addpath('~/James_scripts/bruce/processing/');
addpath('~/James_scripts/bruce/saccade_modulation/');

global Expt_name bar_ori use_LOOXV use_MUA

Expt_name = 'G095';
% Expt_name = 'M266';
% use_LOOXV = 0; %[0 no LOOXV; 1 SU LOOXV; 2 all LOOXV]

use_MUA = false;
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

if strcmp(rec_type,'UA') && bar_ori == 90
et_mod_data_name = 'full_eyetrack_initmods_vbars';
et_anal_name = 'full_eyetrack_vbars';
end

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

is_TBT_expt = false;
if Expt_num >= 275
    is_TBT_expt = true;
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
SU_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;
SU_isodist = Clust_data.SU_isodist;
SU_Lratio = Clust_data.SU_Lratio;

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
sac_direction = [saccades(:).direction];
is_micro = sac_amps(used_saccade_set) < 1;
big_sacs = find(~is_micro & ~used_is_blink');
micro_sacs = find(is_micro & ~used_is_blink');
sac_durs = [saccades(:).duration];
sac_prepos = reshape([saccades(:).pre_pos],[],length(saccades));
sac_postpos = reshape([saccades(:).post_pos],[],length(saccades));
sac_deltaX = sac_postpos(1,:) - sac_prepos(1,:);

sac_burst_isi = 0.15;
sacburst_set = find([saccades(used_saccade_set).isi] < sac_burst_isi | [saccades(used_saccade_set).next_isi] < sac_burst_isi);
micro_sacs(ismember(micro_sacs,sacburst_set)) = [];

saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));

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
    elseif ~isnan(all_mod_SU(ss))
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
backlag = round(0.1/dt);
forlag = round(0.25/dt);
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


%% 

anal_dir = ['/home/james/Analysis/bruce/' Expt_name '/sac_mod/'];
fname = 'sacStimProc';
if strcmp(rec_type,'UA') && bar_ori == 90
    fname = [fname '_vbars'];
end
cd(anal_dir)
load(fname);
%%
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
best_shift_stimmat_up = all_stimmat_up;
for i=1:NT
    best_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
end

%%
if use_MUA
    targs = tr_set;
else
    targs = tr_set(all_mod_SU(tr_set) > 0);
end

for cc = targs
    tr_cc = find(tr_set(loo_set) == cc);
    cc_uinds = find(~isnan(Robs_mat(:,cc)));
    cur_Robs = Robs_mat(cc_uinds,cc);
        
    %% COMPUTE UNIT DATA
    unit_data.avg_rate = mean(cur_Robs)/dt;
    unit_data.tot_spikes = sum(cur_Robs);
    unit_data.N_used_samps = length(cur_Robs);
    used_blocks = unique(all_blockvec(used_inds(cc_uinds)));
    unit_data.n_used_blocks = length(used_blocks);
    block_rates = nan(unit_data.n_used_blocks,1);
    for ii = 1:unit_data.n_used_blocks
        block_rates(ii) = mean(cur_Robs(all_blockvec(used_inds(cc_uinds)) == used_blocks(ii)));
    end
    unit_data.stability_cv = std(block_rates)/mean(block_rates);
    
    if cc > n_probes
        unit_data.SU_number = SU_numbers(cc-n_probes);
        unit_data.probe_number = SU_probes(cc-n_probes);
        unit_data.SU_Lratio = SU_Lratio(cc-n_probes);
        unit_data.SU_isodist = SU_isodist(cc-n_probes);
    else
        unit_data.SU_number = nan;
        unit_data.probe_number = cc;
        unit_data.SU_Lratio = nan;
        unit_data.SU_isodist = nan;
    end
        
    sacStimProc(cc).unit_data = unit_data;
    
     %% RECON RETINAL STIM
%      if cc > n_probes
%          fprintf('Reconstructing retinal stim for unit %d\n',cc);
%          fin_fix_corr = nan(NT,1);
%          fin_fix_std = nan(NT,1);
%          fin_fix_corr(~isnan(fix_ids)) = squeeze(it_fix_post_mean_LOO(tr_cc,end,fix_ids(~isnan(fix_ids))));
%          fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
%          fin_fix_std(~isnan(fix_ids)) = squeeze(it_fix_post_std_LOO(tr_cc,end,fix_ids(~isnan(fix_ids))));
%          fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);
%          
%          fin_drift_corr = squeeze(drift_post_mean_LOO(tr_cc,end,:));
%          fin_drift_std = squeeze(drift_post_std_LOO(tr_cc,end,:));
%          
%          for ii = 1:length(trial_start_inds)
%              cur_inds = trial_start_inds(ii):trial_end_inds(ii);
%              fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
%              fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
%          end
%          fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT);
%          fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);
%          
%          fin_tot_corr = fin_fix_corr + fin_drift_corr;
%          fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2);
%          
%          fin_shift_cor = round(fin_tot_corr);
%          fin_shift_cor(isnan(fin_shift_cor)) = 0;
%          
%          %RECOMPUTE XMAT
%          all_shift_stimmat_up = all_stimmat_up;
%          for i=1:NT
%              %     all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_post_cor(i)};
%              all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
%          end
%          all_Xmat_shift = create_time_embedding(all_shift_stimmat_up,stim_params_us);
%          all_Xmat_shift = all_Xmat_shift(used_inds(cc_uinds),use_kInds_up);
%          
%      else
         all_Xmat_shift = create_time_embedding(best_shift_stimmat_up,stim_params_us);
         all_Xmat_shift = all_Xmat_shift(used_inds(cc_uinds),use_kInds_up);
%      end
    %% FIT STIM-PROCESSING MODEL
%     fprintf('Fitting stim model for unit %d\n',cc);
%     
%     if expt_dds(cur_block_set(1)) == 12
%         base_lambda_d2XT = 50;
%         base_lambda_L1 = 5;
%         base_qlambda_d2XT = 25;
%         base_qlambda_L1 = 5;
%     else
%         base_lambda_d2XT = 200;
%         base_lambda_L1 = 5;
%         base_qlambda_d2XT = 100;
%         base_qlambda_L1 = 5;
%     end
%     init_reg_params = NMMcreate_reg_params('lambda_d2XT',5);
%     fin_stim_params(1) = NMMcreate_stim_params([flen use_nPix_us],dt);
%     
%     silent = 1;
%     
%     stim_mod_signs = [1 1 1 -1 -1];
%     stim_Xtargs = [ones(1,length(stim_mod_signs))];
%     stim_NL_types = {'lin','quad','quad','quad','quad'};
%     qfilts = find(strcmp('quad',stim_NL_types) & stim_Xtargs == 1);
%     nqfilts = find(~strcmp('quad',stim_NL_types) & stim_Xtargs == 1);
%     init_mod = NMMinitialize_model( fin_stim_params, stim_mod_signs, stim_NL_types, init_reg_params,stim_Xtargs);
%     Xtargs = [init_mod.mods(:).Xtarget];
%     
%     init_fitN = ceil(length(cc_uinds)/5);
%     init_fit_subset = randperm(length(cc_uinds));
%     init_fit_subset = init_fit_subset(1:init_fitN);
%     
%     clear optim_params
%     init_mod = NMMfit_filters(init_mod,cur_Robs(init_fit_subset),all_Xmat_shift(init_fit_subset,:),[],[],silent);
%     [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(init_mod,cur_Robs,all_Xmat_shift);
%     vgint = var(gint);
%     sgint = std(gint);
%     init_mod = NMMadjust_regularization(init_mod,nqfilts,'lambda_d2XT',base_lambda_d2XT./vgint(nqfilts));
%     init_mod = NMMadjust_regularization(init_mod,nqfilts,'lambda_L1',base_lambda_L1./sgint(nqfilts));
%     if ~isempty(qfilts)
%         init_mod = NMMadjust_regularization(init_mod,qfilts,'lambda_d2XT',base_qlambda_d2XT./vgint(qfilts));
%         init_mod = NMMadjust_regularization(init_mod,qfilts,'lambda_L1',base_qlambda_L1./sgint(qfilts));
%     end
%     cor_gqm = NMMfit_filters(init_mod,cur_Robs,all_Xmat_shift,[],[],silent);
% 
%     %convert linear filter into separate thresh-lin filters and refit
%         new_mod = cor_gqm;
%     cor_gqm.mods = [cor_gqm.mods cor_gqm.mods(1)];
%     cor_gqm.mods(end).sign = -1;
%     cor_gqm.mods(end).filtK = -cor_gqm.mods(end).filtK;
%     cor_gqm.mods(1).NLtype = 'threshlin'; cor_gqm.mods(end).NLtype = 'threshlin';
%     stim_NL_types{1} = 'threshlin'; stim_NL_types = cat(2,stim_NL_types,{'threshlin'});
%     stim_mod_signs = [cor_gqm.mods(:).sign];
%     cor_gqm = NMMfit_filters(cor_gqm,cur_Robs,all_Xmat_shift,[],[],silent);
% 
%     sacStimProc(cc).stimMod = cor_gqm;
%     [LL, penLL, pred_rate, G, gint,fgint,nullLL] = NMMmodel_eval(cor_gqm,cur_Robs,all_Xmat_shift);
%     sacStimProc(cc).LLimp = (LL-nullLL)/log(2);
%     sacStimProc(cc).rel_filt_contr = std(fgint);
%     
%     sacStimProc(cc).lags = slags;
%     sacStimProc(cc).dt = dt;

cor_gqm = sacStimProc(cc).stimMod;
stim_mod_signs = [cor_gqm.mods(:).sign];

    %% FIT UPSTREAM STIM-MODULATION
    any_sac_inds = find(any(Xmsac(cc_uinds,:) > 0,2));
    cur_Xsac = Xmsac(cc_uinds,:);

    other_units = Robs_mat(cc_uinds,:);    
    other_units(:,cc) = [];
    other_units(isnan(other_units)) = 0;
    if cc > n_probes
        other_units(:,unit_data.probe_number) = [];
    end
    other_units = zscore(other_units);
%     other_units = nanmean(other_units,2);
    
    X{1} = all_Xmat_shift;
    X{2} = other_units;
    cor_gqm_net = cor_gqm;
    cor_gqm_net.stim_params(2) = NMMcreate_stim_params(size(other_units,2));
    cor_gqm_net = NMMadd_NLinput(cor_gqm_net,{'lin'},1,2);
    cor_gqm_net.mods(end).filtK(:) = 0;
    cor_gqm_net.mods(end).reg_params = NMMcreate_reg_params();
    cor_gqm_net = NMMfit_filters(cor_gqm_net,cur_Robs,X,[],length(cor_gqm_net.mods),0);
    
    
    %% FIT POST-INTEGRATION GAIN
    [LL, penLL, pred_rate, G,gint,fgint] = NMMmodel_eval(cor_gqm, cur_Robs, X);
    fgint(:,stim_mod_signs == -1) = -fgint(:,stim_mod_signs == -1);
    
    lambda_d2T = 50;
    lambda_L2 = 5;
%     lambda_d2T = 5;
%     lambda_L2 = 0.1;
    sac_reg_params = NMMcreate_reg_params('lambda_d2T',lambda_d2T,'lambda_L2',lambda_L2,'boundary_conds',[0 0 0]);
    
    optim_params.optTol = 1e-6;
    optim_params.progTol = 1e-10;
    optim_params.Method = 'lbfgs';
    optim_params.verbose = 1;        
    
    G = sum(fgint(:,1:end-1),2);
    g_tot = G - cor_gqm.spk_NL_params(1);
    Xsac_tot = bsxfun(@times,cur_Xsac,g_tot);
 
    net_out = X{2}*cor_gqm_net.mods(end).filtK;
    Xsac_net = bsxfun(@times,cur_Xsac,net_out);

    silent = 0;
    
    clear tr_stim
    tr_stim{1} = [g_tot];
    tr_stim{2} = cur_Xsac;
    tr_stim{3} = Xsac_tot;
    tr_stim{4} = net_out;
    tr_stim{5} = Xsac_net;
    
    clear sac_stim_params
    sac_stim_params(1) = NMMcreate_stim_params(1);
    sac_stim_params(2) = NMMcreate_stim_params([size(Xsac_tot,2)]);
    sac_stim_params(3) = NMMcreate_stim_params([size(Xsac_tot,2)]);
    sac_stim_params(4) = NMMcreate_stim_params(1);
    sac_stim_params(5) = NMMcreate_stim_params([size(Xsac_tot,2)]);

    mod_signs = [1 1 1 1 1];
    Xtargets = [1 2 3 4 5];
    NL_types = {'lin','lin','lin','lin','lin'};
    spost_gsac_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
    spost_gsac_mod.mods(4).filtK(:) = 0;
    spost_gsac_mod.mods(5).filtK(:) = 0;
    spost_gsac_mod = NMMfit_filters(spost_gsac_mod,cur_Robs,tr_stim,[],[1:3],silent);
    spost_gsac_mod = NMMfit_filters(spost_gsac_mod,cur_Robs,tr_stim,[],[4:5],silent);

    sacStimProc(cc).gsac_post_singmod = spost_gsac_mod;
    [spost_LL,~,spost_pred_rate,~,~,~,nullLL] = NMMmodel_eval( spost_gsac_mod, cur_Robs, tr_stim);

end
%     %% COMPUTE MODEL-BASED INFORMATION
    
    %     [LL,pred_rate] = eval_sacgain_mod( sacGainMod, cur_Robs, all_Xmat_shift, cur_Xsac);
    %     cur_info = nan(length(slags),1);
    [spost_LL,~,spost_pred_rate,~,~,~,nullLL] = NMMmodel_eval( spost_gsac_mod, cur_Robs, tr_stim);
    [sac_spost_info,sac_info,sac_LL,sac_spost_LL,sac_nullLL,sac_Nspks] = deal(nan(length(slags),1));
    for ii = 1:length(slags)
        temp = find(cur_Xsac(:,ii) == 1);
        cur_avg_rate = mean(cur_Robs(temp))*ones(size(temp));
        sac_nullLL(ii) = nansum(cur_Robs(temp).*log2(cur_avg_rate) - cur_avg_rate);
        sac_Nspks(ii) = sum(cur_Robs(temp));
        
        %         sac_LL(ii) = nansum(cur_Robs(temp).*log2(pred_rate(temp)) - pred_rate(temp));
        %         sac_info(ii) = nanmean(pred_rate(temp).*log2(pred_rate(temp)/mean(pred_rate(temp))))/mean(pred_rate(temp));
        
        sac_spost_LL(ii) = nansum(cur_Robs(temp).*log2(spost_pred_rate(temp)) - spost_pred_rate(temp));
        sac_spost_info(ii) = nanmean(spost_pred_rate(temp).*log2(spost_pred_rate(temp)/mean(spost_pred_rate(temp))))/mean(spost_pred_rate(temp));
    end
    
    %    sacStimProc(cc).gsac_LL_info = (sac_LL - sac_nullLL)./sac_Nspks;
    %     sacStimProc(cc).gsac_ovLL_info = (LL-nullLL)/log(2);
    %    sacStimProc(cc).gsac_ov_info = mean(pred_rate/mean(pred_rate).*log2(pred_rate/mean(pred_rate)));
    %    sacStimProc(cc).gsac_mod_info = sac_info;
    
    sacStimProc(cc).gsac_spostLL_info = (sac_spost_LL - sac_nullLL)./sac_Nspks;
    sacStimProc(cc).gsac_ov_spost_LLinfo = (spost_LL-nullLL)/log(2);
    sacStimProc(cc).gsac_spostmod_info = sac_spost_info;
    sacStimProc(cc).gsac_spost_ov_info = mean(spost_pred_rate/mean(spost_pred_rate).*log2(spost_pred_rate/mean(spost_pred_rate)));
%     
% 
%    %% COMPUTE MODEL SEQUENCE
% %    %compute best time lag over E filts
% %    [Xinds,Tinds] = meshgrid(1:use_nPix_us,1:flen);
% %    cur_filts = reshape([cor_gqm.mods(1:3).filtK],[flen use_nPix_us 3]);
% %    cur_tfilt = squeeze(mean(std(cur_filts,[],2),3));
% %    [~,best_lag] = max(cur_tfilt);
% %    uk = find(Tinds == best_lag);
% %    
% %    newX = reshape(all_Xmat_shift(any_sac_inds,uk),[length(any_sac_inds) 1 use_nPix_us]);
% %    X{1} = cur_Xsac(any_sac_inds,:);
% %    X{2} = reshape(bsxfun(@times,newX,cur_Xsac(any_sac_inds,:)),length(any_sac_inds),[]);
% %    
% %    cur_stim_params(1) = NMMcreate_stim_params(length(slags));
% %    cur_stim_params(2) = NMMcreate_stim_params([length(slags) use_nPix_us]);
% %    
% %    optim_params.optTol = 1e-5;
% %    optim_params.progTol = 1e-9;
% %    mod_signs = [1 1 1 1];
% %    Xtargs = [1 2 2 2];
% %    NL_types = {'lin','lin','quad','quad'};
% %    reg_params = NMMcreate_reg_params('lambda_d2XT',[0 50 25 25]','lambda_d2T',[10 0 0 0]');
% %    init_mod = NMMinitialize_model(cur_stim_params,mod_signs,NL_types,reg_params,Xtargs);
% %    init_mod = NMMfit_filters(init_mod,cur_Robs(any_sac_inds),X,[],[],0,optim_params);
% %    
% %    cur_filts = reshape([init_mod.mods(2:end).filtK],[length(slags) use_nPix_us 3]);
% %    phasedep_filt = squeeze(cur_filts(:,:,1));
% %    phaseind_filt = squeeze(sqrt(sum(cur_filts(:,:,2:3).^2,3)));
% %    sacStimProc(cc).gsac_phaseDep_filt = phasedep_filt;
% %    sacStimProc(cc).gsac_phaseInd_filt = phaseind_filt;
% 
%     %% %% EVALUATE TIME-DEPENDENT INFO
% %     cur_Xsac = Xsac(cc_uinds,:);
% %     [best_LL,pred_rate] = eval_sacgain_mod( sacGainMod, cur_Robs, all_Xmat_shift,cur_Xsac);
% %     
% %     temp = Xinds_up;
% %     temp(use_kInds_up) = nan;
% %     used_pix = find(isnan(temp(1,:)));
% %     used_stimmat_up = all_shift_stimmat_up(:,used_pix);
% %     new_stim_params = NMMcreate_stim_params([flen use_nPix_us]);
% % 
% %     shuf_LL = nan(length(slags),1);
% %     for ii = 1:length(slags)
% %         fprintf('Testing pre-scramble at gsac lag %d of %d\n',ii,length(slags));
% %         to_shuffle = find(t_since_sac_start == slags(ii));
% %         shuf_stim = used_stimmat_up;
% %         shuf_stim(used_inds(to_shuffle),:) = used_stimmat_up(randi(NT,length(to_shuffle),1),:);
% %         shuf_X = create_time_embedding(shuf_stim,new_stim_params);
% %         shuf_X = shuf_X(used_inds(cc_uinds),:);
% % 
% %         [shuf_LL(ii)] = eval_sacgain_mod( sacGainMod, cur_Robs, shuf_X, cur_Xsac);
% %         n_frames(ii) = sum(t_since_sac_start == slags(ii));
% %     end
% %     sacStimProc(cc).gsac_stimshuff_info = (best_LL-shuf_LL')./n_frames*length(cc_uinds)/log(2);
%     
%     %% compute post-shuffle info    
% %     for ii = 1:length(slags)
% %         fprintf('Testing post-scramble at gsac lag %d of %d\n',ii,length(slags));
% %         to_shuffle = find(t_since_sac_start(cc_uinds)==slags(ii));
% %         shuf_X = all_Xmat_shift;
% %         shuf_X(to_shuffle,:) = shuf_X(randi(length(cc_uinds),length(to_shuffle),1),:);
% %                 
% %         [shuf_LL(ii)] = eval_sacgain_mod( sacGainMod, cur_Robs, shuf_X, cur_Xsac);
% %         n_frames(ii) = sum(t_since_sac_start == slags(ii));
% %     end
% %     sacStimProc(cc).gsac_shuff_info = (best_LL-shuf_LL')./n_frames*length(cc_uinds)/log(2);
% 
%     %% CREATE TENT_BASIS MODEL OF SACCADE-MODULATION
%     n_Gbins = 35;
%     TB_lambda = 20;
%     Xtick = -(backlag-1/2):(1):(forlag+1/2);
%     n_sbins = length(Xtick);
%     addpath('~/James_scripts/TentBasis2D/');
%     [LL, penLL, pred_rate, G,gint,fgint] = NMMmodel_eval(cor_gqm, cur_Robs, all_Xmat_shift);
%     g_tot = G - cor_gqm.spk_NL_params(1);
%         
%     cur_sac_starts = saccade_start_inds(big_sacs);
%     cur_sac_stops = saccade_stop_inds(big_sacs);
%     t_since_sac_start = nan(NT,1);
%     for ii = 1:length(cur_sac_starts)
%         prev_tstart = find(trial_start_inds <= cur_sac_starts(ii),1,'last');
%         next_tstop = find(trial_end_inds >= cur_sac_starts(ii),1,'first');
%         cur_inds = (cur_sac_starts(ii) - backlag):(cur_sac_starts(ii) + forlag);
%         cur_uset = find(cur_inds > trial_start_inds(prev_tstart) & cur_inds < trial_end_inds(next_tstop));
%         t_since_sac_start(cur_inds(cur_uset)) = slags(cur_uset);
%     end
%     
%     TB_stim = [t_since_sac_start(cc_uinds) g_tot];
%     Ytick = linspace(my_prctile(TB_stim(:,2),0.1),my_prctile(TB_stim(:,2),100-1),n_Gbins);
%     TB = TentBasis2D(Xtick, Ytick);
%     
%     used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
%         TB_stim(:,2) >= Ytick(1) & TB_stim(:,2) <= Ytick(end));
%     [TB_Xmat,TB_counts] = TB.InputNL2D(TB_stim(used_data,:));
%     L2_params = create_L2_params([],[1 n_sbins*n_Gbins],[n_sbins n_Gbins],2,3,[Inf Inf],[0.025 1]);
%     TB_fitmod = regGLM_fit(TB_Xmat,cur_Robs(used_data),L2_params,TB_lambda,[],[],1);
%     [LL, penLL, pred_rate, G] = regGLM_eval(TB_fitmod,cur_Robs(used_data),TB_Xmat);
%     TB_K = reshape(TB_fitmod.K,n_sbins,n_Gbins)';
%     bin_areas = TB.GetBinAreas();
%     gsac_TB_dist = TB_counts./bin_areas;
%     gsac_TB_dist = gsac_TB_dist'/sum(gsac_TB_dist(:));
%     gsac_TB_rate = log(1 + exp(TB_K + TB_fitmod.theta));
%     sacStimProc(cc).gsac_TB_rate = gsac_TB_rate;
%     
%     %INFO CALS
%     cur_avg_rate = mean(cur_Robs(used_data));
%     marg_gdist = sum(gsac_TB_dist,2);
%     marg_sdist = sum(gsac_TB_dist);
%     marg_gsacrate = sum(gsac_TB_dist.*gsac_TB_rate)./marg_sdist;
%     marg_grate = sum(gsac_TB_dist.*gsac_TB_rate,2)./marg_gdist;
%     gsacdep_info = nan(1,n_sac_bins);
%     for tt = 1:n_sbins
%         gsacdep_info(tt) = sum(gsac_TB_dist(:,tt).*gsac_TB_rate(:,tt).*log2(gsac_TB_rate(:,tt)/marg_gsacrate(tt)))/sum(gsac_TB_dist(:,tt));
%     end
%     gcumdist = cumsum(marg_gdist)/sum(marg_gdist);
%     
%     sacStimProc(cc).gsac_ov_TB_info = sum(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate))/cur_avg_rate;
%     
%     sacStimProc(cc).gsac_avg_rate = marg_gsacrate;
%     sacStimProc(cc).gsac_TB_info = gsacdep_info./marg_gsacrate;
%     sacStimProc(cc).gsac_TB_gdist = marg_gdist;
%     sacStimProc(cc).gsac_TB_grate = marg_grate;
%     
%     %% NOW FOR MSACS
%     %% FIT UPSTREAM STIM-MODULATION
%     any_sac_inds = find(any(Xmsac(cc_uinds,:) > 0,2));
%     cur_Xsac = Xmsac(cc_uinds,:);
%     tr_sac_inds = any_sac_inds(ismember(cc_uinds(any_sac_inds),tr_inds));
%     xv_sac_inds = any_sac_inds(ismember(cc_uinds(any_sac_inds),xv_inds));
% 
% %     lambda_d2T = 50;
% %     lambda_L2 = 5;
% %     Xsac_mat = cur_Xsac(any_sac_inds,:);
% %     sacGainMod = fit_sacgain_model(cor_gqm,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:),Xsac_mat,lambda_d2T,lambda_L2);
% %     sacStimProc(cc).msacGainMod = sacGainMod;
%       
%             %% FIT POST-INTEGRATION GAIN
%     [LL, penLL, pred_rate, G,gint,fgint] = NMMmodel_eval(cor_gqm, cur_Robs, all_Xmat_shift);
%     fgint(:,stim_mod_signs == -1) = -fgint(:,stim_mod_signs == -1);
%     
%     lambda_d2T = 50;
%     lambda_L2 = 5;
%     sac_reg_params = NMMcreate_reg_params('lambda_d2T',lambda_d2T,'lambda_L2',lambda_L2,'boundary_conds',[0 0 0]);
%     
%     g_exc = sum(fgint(:,stim_mod_signs==1),2);
%     g_inh = sum(fgint(:,stim_mod_signs==-1),2);
%     Xsac_estim = bsxfun(@times,cur_Xsac,g_exc);
%     Xsac_istim = bsxfun(@times,cur_Xsac,g_inh);
%     
%     optim_params.optTol = 1e-6;
%     optim_params.progTol = 1e-10;
%     optim_params.Method = 'lbfgs';
%     optim_params.verbose = 1;
%     
%     clear tr_stim
%     tr_stim{1} = [g_exc g_inh];
%     tr_stim{2} = cur_Xsac;
%     tr_stim{3} = Xsac_estim;
%     tr_stim{4} = Xsac_istim;
%     sac_stim_params(1) = NMMcreate_stim_params(2);
%     sac_stim_params(2) = NMMcreate_stim_params([size(Xsac_estim,2)]);
%     sac_stim_params(3) = NMMcreate_stim_params([size(Xsac_estim,2)]);
%     sac_stim_params(4) = NMMcreate_stim_params([size(Xsac_estim,2)]);
%     mod_signs = [1 1 1 1];
%     Xtargets = [1 2 3 4];
%     NL_types = {'lin','lin','lin','lin'};
%     post_msac_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
%     post_msac_mod.mods(1).reg_params = NMMcreate_reg_params();
%     post_msac_mod = NMMfit_filters(post_msac_mod,cur_Robs,tr_stim,[],[],silent,optim_params);
%     
%     sacStimProc(cc).msac_post_fullmod = post_msac_mod;
% 
%     g_tot = G - cor_gqm.spk_NL_params(1);
%     Xsac_tot = bsxfun(@times,cur_Xsac,g_tot);
%     clear tr_stim
%     tr_stim{1} = [g_tot];
%     tr_stim{2} = cur_Xsac;
%     tr_stim{3} = Xsac_tot;
%     clear sac_stim_params
%     sac_stim_params(1) = NMMcreate_stim_params(1);
%     sac_stim_params(2) = NMMcreate_stim_params([size(Xsac_tot,2)]);
%     sac_stim_params(3) = NMMcreate_stim_params([size(Xsac_tot,2)]);
% 
%     mod_signs = [1 1 1];
%     Xtargets = [1 2 3];
%     NL_types = {'lin','lin','lin'};
%     spost_msac_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
%     spost_msac_mod = NMMfit_filters(spost_msac_mod,cur_Robs,tr_stim,[],[],silent);
% 
%     sacStimProc(cc).msac_post_singmod = spost_msac_mod;
% 
%     %%
%     %     [LL,pred_rate] = eval_sacgain_mod( sacGainMod, cur_Robs, all_Xmat_shift, cur_Xsac);
%     %     cur_info = nan(length(slags),1);
%     [spost_LL,~,spost_pred_rate,~,~,~,nullLL] = NMMmodel_eval( spost_msac_mod, cur_Robs, tr_stim);
%     [sac_spost_info,sac_info,sac_LL,sac_spost_LL,sac_nullLL,sac_Nspks] = deal(nan(length(slags),1));
%     for ii = 1:length(slags)
%         temp = find(cur_Xsac(:,ii) == 1);
%         cur_avg_rate = mean(cur_Robs(temp))*ones(size(temp));
%         sac_nullLL(ii) = nansum(cur_Robs(temp).*log2(cur_avg_rate) - cur_avg_rate);
%         sac_Nspks(ii) = sum(cur_Robs(temp));
%         
%         %         sac_LL(ii) = nansum(cur_Robs(temp).*log2(pred_rate(temp)) - pred_rate(temp));
%         %         sac_info(ii) = nanmean(pred_rate(temp).*log2(pred_rate(temp)/mean(pred_rate(temp))))/mean(pred_rate(temp));
%         
%         sac_spost_LL(ii) = nansum(cur_Robs(temp).*log2(spost_pred_rate(temp)) - spost_pred_rate(temp));
%         sac_spost_info(ii) = nanmean(spost_pred_rate(temp).*log2(spost_pred_rate(temp)/mean(spost_pred_rate(temp))))/mean(spost_pred_rate(temp));
%     end
%     
%     %    sacStimProc(cc).msac_LL_info = (sac_LL - sac_nullLL)./sac_Nspks;
%     %     sacStimProc(cc).msac_ovLL_info = (LL-nullLL)/log(2);
%     %    sacStimProc(cc).msac_ov_info = mean(pred_rate/mean(pred_rate).*log2(pred_rate/mean(pred_rate)));
%     %    sacStimProc(cc).msac_mod_info = sac_info;
%     
%     sacStimProc(cc).msac_spostLL_info = (sac_spost_LL - sac_nullLL)./sac_Nspks;
%     sacStimProc(cc).msac_spost_ovLL_info = (spost_LL-nullLL)/log(2);
%     sacStimProc(cc).msac_spostmod_info = sac_spost_info;
%     sacStimProc(cc).msac_spost_ov_info = mean(spost_pred_rate/mean(spost_pred_rate).*log2(spost_pred_rate/mean(spost_pred_rate)));
%     
% 
%    %% COMPUTE MODEL SEQUENCE
% %    %compute best time lag over E filts
% %    [Xinds,Tinds] = meshgrid(1:use_nPix_us,1:flen);
% %    cur_filts = reshape([cor_gqm.mods(1:3).filtK],[flen use_nPix_us 3]);
% %    cur_tfilt = squeeze(mean(std(cur_filts,[],2),3));
% %    [~,best_lag] = max(cur_tfilt);
% %    uk = find(Tinds == best_lag);
% %    
% %    newX = reshape(all_Xmat_shift(any_sac_inds,uk),[length(any_sac_inds) 1 use_nPix_us]);
% %    X{1} = cur_Xsac(any_sac_inds,:);
% %    X{2} = reshape(bsxfun(@times,newX,cur_Xsac(any_sac_inds,:)),length(any_sac_inds),[]);
% %    
% %    cur_stim_params(1) = NMMcreate_stim_params(length(slags));
% %    cur_stim_params(2) = NMMcreate_stim_params([length(slags) use_nPix_us]);
% %    
% %    optim_params.optTol = 1e-5;
% %    optim_params.progTol = 1e-9;
% %    mod_signs = [1 1 1 1];
% %    Xtargs = [1 2 2 2];
% %    NL_types = {'lin','lin','quad','quad'};
% %    reg_params = NMMcreate_reg_params('lambda_d2XT',[0 50 25 25]','lambda_d2T',[10 0 0 0]');
% %    init_mod = NMMinitialize_model(cur_stim_params,mod_signs,NL_types,reg_params,Xtargs);
% %    init_mod = NMMfit_filters(init_mod,cur_Robs(any_sac_inds),X,[],[],0,optim_params);
% %    
% %    cur_filts = reshape([init_mod.mods(2:end).filtK],[length(slags) use_nPix_us 3]);
% %    phasedep_filt = squeeze(cur_filts(:,:,1));
% %    phaseind_filt = squeeze(sqrt(sum(cur_filts(:,:,2:3).^2,3)));
% %    sacStimProc(cc).msac_phaseDep_filt = phasedep_filt;
% %    sacStimProc(cc).msac_phaseInd_filt = phaseind_filt;
%    
%     %%
%     n_Gbins = 35;
%     TB_lambda = 20;
%     Xtick = -(backlag-1/2):(1):(forlag+1/2);
%     n_sbins = length(Xtick);
%     addpath('~/James_scripts/TentBasis2D/');
%     [LL, penLL, pred_rate, G,gint,fgint] = NMMmodel_eval(cor_gqm, cur_Robs, all_Xmat_shift);
%     g_tot = G - cor_gqm.spk_NL_params(1);
%     
%     cur_sac_starts = saccade_start_inds(micro_sacs);
%     cur_sac_stops = saccade_stop_inds(micro_sacs);
%     t_since_sac_start = nan(NT,1);
%     for ii = 1:length(cur_sac_starts)
%         prev_tstart = find(trial_start_inds <= cur_sac_starts(ii),1,'last');
%         next_tstop = find(trial_end_inds >= cur_sac_starts(ii),1,'first');
%         cur_inds = (cur_sac_starts(ii) - backlag):(cur_sac_starts(ii) + forlag);
%         cur_uset = find(cur_inds > trial_start_inds(prev_tstart) & cur_inds < trial_end_inds(next_tstop));
%         t_since_sac_start(cur_inds(cur_uset)) = slags(cur_uset);
%     end
%     
%     TB_stim = [t_since_sac_start(cc_uinds) g_tot];
%     Ytick = linspace(my_prctile(TB_stim(:,2),0.1),my_prctile(TB_stim(:,2),100-1),n_Gbins);
%     TB = TentBasis2D(Xtick, Ytick);
%     
%     used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
%         TB_stim(:,2) >= Ytick(1) & TB_stim(:,2) <= Ytick(end));
%     [TB_Xmat,TB_counts] = TB.InputNL2D(TB_stim(used_data,:));
%     L2_params = create_L2_params([],[1 n_sbins*n_Gbins],[n_sbins n_Gbins],2,3,[Inf Inf],[0.025 1]);
%     TB_fitmod = regGLM_fit(TB_Xmat,cur_Robs(used_data),L2_params,TB_lambda,[],[],1);
%     [LL, penLL, pred_rate, G] = regGLM_eval(TB_fitmod,cur_Robs(used_data),TB_Xmat);
%     TB_K = reshape(TB_fitmod.K,n_sbins,n_Gbins)';
%     bin_areas = TB.GetBinAreas();
%     msac_TB_dist = TB_counts./bin_areas;
%     msac_TB_dist = msac_TB_dist'/sum(msac_TB_dist(:));
%     msac_TB_rate = log(1 + exp(TB_K + TB_fitmod.theta));
%     sacStimProc(cc).msac_TB_rate = msac_TB_rate;
%     
%     %INFO CALS
%     cur_avg_rate = mean(cur_Robs(used_data));
%     marg_gdist = sum(msac_TB_dist,2);
%     marg_sdist = sum(msac_TB_dist);
%     marg_msacrate = sum(msac_TB_dist.*msac_TB_rate)./marg_sdist;
%     marg_grate = sum(msac_TB_dist.*msac_TB_rate,2)./marg_gdist;
%     ov_info_msac = sum(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate));
%     msacdep_info = nan(1,n_sac_bins);
%     for tt = 1:n_sbins
%         msacdep_info(tt) = sum(msac_TB_dist(:,tt).*msac_TB_rate(:,tt).*log2(msac_TB_rate(:,tt)/marg_msacrate(tt)))/sum(msac_TB_dist(:,tt));
%     end
%     gcumdist = cumsum(marg_gdist)/sum(marg_gdist);
%     
%     sacStimProc(cc).msac_avg_rate = marg_msacrate;
%     sacStimProc(cc).msac_TB_info = msacdep_info;
%     sacStimProc(cc).msac_TB_gdist = marg_gdist;
%     sacStimProc(cc).msac_TB_grate = marg_grate;
% end

%%
% anal_dir = ['/home/james/Analysis/bruce/' Expt_name '/sac_mod/'];
% fname = 'sacStimProc2';
% if strcmp(rec_type,'UA') && bar_ori == 90
%     fname = [fname '_vbars'];
% end
% cd(anal_dir)
% % save(fname,'targs','slags','dt','all_gqm_mods','flen','use_nPix*','*gsac*','*msac*');
% save(fname,'targs','slags','dt','sacStimProc');
% 
