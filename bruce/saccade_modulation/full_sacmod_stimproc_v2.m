%
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements//');
addpath('~/James_scripts/bruce/processing/');
addpath('~/James_scripts/bruce/saccade_modulation/');
addpath('~/James_scripts/TentBasis2D/');

global Expt_name bar_ori use_MUA


Expt_name = 'M296';
% Expt_name = 'G093';
% use_MUA = false;
bar_ori = 0; %bar orientation to use (only for UA recs)


fit_unCor = false;
fit_subMod = true;
fitModSeq = false;
fitUpstream = true;
fitSTA = true;
fitMsacs = false;

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
        case 296
            bar_ori = 45;
    end
end

if Expt_num >= 280
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
    use_LOOXV = 1;
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
mod_data_dir = ['~/Analysis/bruce/' Expt_name '/models'];

% et_mod_data_name = 'full_eyetrack_initmods';
% et_anal_name = 'full_eyetrack';
et_mod_data_name = 'full_eyetrack_initmods_Rinit';
et_anal_name = 'full_eyetrack_Rinit';
mod_data_name = 'corrected_models';

if strcmp(rec_type,'UA') && bar_ori == 90
    et_mod_data_name = 'full_eyetrack_initmods_vbars';
    et_anal_name = 'full_eyetrack_vbars';
    mod_data_name = 'corrected_models_vbars';
end

if fit_unCor
    mod_data_name = [mod_data_name '_unCor'];
end

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

flen = 15;
spatial_usfac = 2;

%these recs have larger bar widths
if ismember(Expt_num,[287 289 294])
    use_nPix = 15;
    spatial_usfac = 4;
elseif ismember(Expt_num,[296])
    use_nPix = 22;
    spatial_usfac = 2;
end

min_trial_dur = 0.75;

stim_fs = 100; %in Hz
dt = 0.01;
Fr = 1;

backlag = round(0.1/dt);
forlag = round(0.3/dt);
slags = -backlag:forlag;
n_sac_bins = length(slags);


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
    case 296
        full_nPix = 54;
end

%exclude data at beginning and end of each trial
beg_buffer = 0.2;
end_buffer = 0.05;
trial_dur = 4;

use_nPix_us = use_nPix*spatial_usfac;
klen_us = use_nPix_us*flen;

stim_params = NIMcreate_stim_params([flen full_nPix],dt);

%% SELECT BLOCKS FOR ANALYSIS
include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaXFaXFs','rls.AllSac','rls.imiXFa','rls.FaXwi','rls.FaXwiXimi','rls.AllSacB'};
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
buffer_pix = floor((full_nPix - use_nPix)/2);

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
fprintf('Loading ET data\n');
cd(et_dir)
load(et_mod_data_name,'all_mod*');
load(et_anal_name,'drift*','it_*','et_tr_set');
tr_set = et_tr_set;

fprintf('Loading model fits\n');
cd(mod_data_dir)
load(mod_data_name);

%%
n_chs = size(all_binned_mua,2) + size(all_binned_sua,2);

% make Robs_mat
poss_blocks = 1:(n_blocks-1);
Robs_mat = nan(length(used_inds),n_chs);
for ss = 1:n_chs
    if ss > n_probes
        su_probe_ind = ss - n_probes;
        Robs_mat(:,ss) = all_binned_sua(used_inds,su_probe_ind);
    else
        Robs_mat(:,ss) = all_binned_mua(used_inds,ss);
    end
end

%% Create set of TR and XV trials
use_trials = unique(all_trialvec(used_inds));
nuse_trials = length(use_trials);

rpt_trials = find(all_trial_Se==rpt_seed);
n_rpt_trials = length(rpt_trials);
tr_trials = setdiff(use_trials,rpt_trials);
n_tr_trials = length(tr_trials);
tr_inds = find(ismember(all_trialvec(used_inds),tr_trials));

%%

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
%set of units where we have LOOXV on eye-tracking
if use_LOOXV == 2
    loo_set = tr_set;
elseif use_LOOXV == 1
    loo_set = tr_set(all_mod_SU(tr_set) > 0);
else
    loo_set = [];
end

%for model fitting
if use_MUA
    targs = 1:n_chs; %SU and MU
else
    targs = setdiff(1:n_chs,1:n_probes); %SU only
end

%% Recon retinal stim for non LOO data
cur_fix_post_mean = squeeze(it_fix_post_mean(end,:));
cur_fix_post_std = squeeze(it_fix_post_std(end,:));
cur_drift_post_mean = squeeze(drift_post_mean(end,:));
cur_drift_post_std = squeeze(drift_post_std(end,:));
[fin_tot_corr,fin_tot_std] = construct_eye_position(cur_fix_post_mean,cur_fix_post_std,...
    cur_drift_post_mean,cur_drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);

fin_shift_cor = round(fin_tot_corr);
fin_shift_cor(isnan(fin_shift_cor)) = 0;

%RECOMPUTE XMAT
best_shift_stimmat_up = all_stimmat_up;
if ~fit_unCor
    for i=1:NT
        best_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
    end
end

%% Model-fitting params
silent = 0;
lambda_d2T = 50;
lambda_L2 = 0;
optim_params.optTol = 1e-6;
optim_params.progTol = 1e-10;

n_Gbins = 35;
TB_lambda = 20;


%%
anal_dir = ['/home/james/Analysis/bruce/' Expt_name '/sac_mod/'];
if ~exist(anal_dir,'dir')
    mkdir(anal_dir)
end
fname = 'sacStimProc_v2';
if strcmp(rec_type,'UA') && bar_ori == 90
    fname = [fname '_vbars'];
end
cd(anal_dir)
% load(fname)
% fprintf('Loading %s\n',fname);
%%

for cc = targs
    
    fprintf('Starting model fits for unit %d\n',cc);
    
    loo_cc = find(loo_set == cc);
    cc_uinds = tr_inds(~isnan(Robs_mat(tr_inds,cc)));
    cur_Robs = Robs_mat(cc_uinds,cc);
    if ~isempty(cc_uinds)
        
        %         other_spike_ids = 1:size(Robs_mat,2);
        %         other_spikes = Robs_mat(cc_uinds,:);
        %         other_spikes(:,cc) = [];
        %         other_spike_ids(cc) = [];
        %         for oo = 1:size(other_spikes,2)
        %             other_spikes(isnan(other_spikes(:,oo)),oo) = nanmean(other_spikes(:,oo));
        %         end
        %         other_spikes = zscore(other_spikes);
        %         bad_spkids = find(any(isnan(other_spikes)));
        %         other_spikes(:,bad_spkids) = [];
        %         other_spike_ids(bad_spkids) = [];
        
        cur_GQM = ModData(cc).rectGQM;
        unit_data = ModData(cc).unit_data;
        sacStimProc(cc).ModData = ModData(cc);
        
        
        sacStimProc(cc).used = true;
        
        fprintf('Reconstructing retinal stim for unit %d\n',cc);
        if ismember(cc,loo_set)
            cur_fix_post_mean = squeeze(it_fix_post_mean_LOO(loo_cc,end,:));
            cur_fix_post_std = squeeze(it_fix_post_std_LOO(loo_cc,end,:));
            cur_drift_post_mean = squeeze(drift_post_mean_LOO(loo_cc,end,:));
            cur_drift_post_std = squeeze(drift_post_std_LOO(loo_cc,end,:));
            [fin_tot_corr,fin_tot_std] = construct_eye_position(cur_fix_post_mean,cur_fix_post_std,...
                cur_drift_post_mean,cur_drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);
            
            fin_shift_cor = round(fin_tot_corr);
            fin_shift_cor(isnan(fin_shift_cor)) = 0;
            
            %RECOMPUTE XMAT
            all_shift_stimmat_up = all_stimmat_up;
            if ~fit_unCor
                for i=1:NT
                    all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
                end
            end
            all_Xmat_shift = create_time_embedding(all_shift_stimmat_up,stim_params_us);
            all_Xmat_shift = all_Xmat_shift(used_inds(cc_uinds),use_kInds_up);
            
        else
            all_Xmat_shift = create_time_embedding(best_shift_stimmat_up,stim_params_us);
            all_Xmat_shift = all_Xmat_shift(used_inds(cc_uinds),use_kInds_up);
        end
        
        %% FIT SPK NL PARAMS
        cur_GQM = NMMfit_logexp_spkNL(cur_GQM,cur_Robs,all_Xmat_shift);
%          new_mod = NMMfit_logexp_spkNL(cur_GQM,cur_Robs,all_Xmat_shift,[],1);
%         [~, ~, cprate, tot_G] = NMMmodel_eval(new_mod, cur_Robs, all_Xmat_shift);

%% Adjust to optimal filter amplitudes (offset effects of regularization)
        [~, ~, cprate, tot_G,ind_Gints,fgint] = NMMmodel_eval(cur_GQM, cur_Robs, all_Xmat_shift);
        used_filts = find(sum(abs(fgint)) > 0);
        
        stim_mod_signs = [cur_GQM.mods(:).sign];
        tempStimParams = NMMcreate_stim_params(length(used_filts),1);
        tempMod = NMMinitialize_model(tempStimParams,1,{'lin'});
        tempMod = NMMfit_filters(tempMod,cur_Robs,fgint(:,used_filts));
        
        gainAdj = zeros(length(cur_GQM.mods),1);
        gainAdj(used_filts) = tempMod.mods(1).filtK;
        fgint = bsxfun(@times,fgint,gainAdj');
        stimG = sum(fgint,2);
        sacStimProc(cc).gainAdj = gainAdj;
        
        for jj = 1:length(cur_GQM.mods)
            if strcmp(cur_GQM.mods(jj).NLtype,'quad')
                cur_adj = sign(gainAdj(jj))*sqrt(abs(gainAdj(jj)));
                %                 cur_adj = sqrt(abs(gainAdj(jj)));
            else
                cur_adj = abs(gainAdj(jj));
            end
            cur_GQM.mods(jj).filtK = cur_adj*cur_GQM.mods(jj).filtK;
        end
        
        sacStimProc(cc).ov_beta = robust_std_dev(tot_G*cur_GQM.spk_NL_params(2));
        sacStimProc(cc).ov_theta = cur_GQM.spk_NL_params(1);
                
        %% FOR GSACS
        any_sac_inds = find(any(Xsac(cc_uinds,:) > 0,2));
        cur_Xsac = Xsac(cc_uinds,:);
        
        %% FIT UPSTREAM STIM-MODULATION
        if fitUpstream
            fprintf('Fitting upstream saccade kernel\n');
            Xsac_mat = Xsac(cc_uinds(any_sac_inds),:);
            sacGainMod = fit_sacgain_model(cur_GQM,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:),Xsac_mat,lambda_d2T,lambda_L2);
            sacStimProc(cc).gsacGainMod = sacGainMod;
        end
        sacGainMod = sacStimProc(cc).gsacGainMod;
        %% FIT POST-INTEGRATION GAIN
        fprintf('Fitting post-filter models\n');
        
        sac_reg_params = NMMcreate_reg_params('lambda_d2T',lambda_d2T,'lambda_L2',lambda_L2,'boundary_conds',[0 0 0]);
        
        %         cur_nmods = length(cur_GQM.mods);
        %         clear tr_stim
        %         tr_stim{1} = stimG;
        %         tr_stim{2} = cur_Xsac;
        %         for kk = 1:cur_nmods
        %             tr_stim{kk+2} = bsxfun(@times,cur_Xsac,fgint(:,kk));
        %         end
        %         sac_stim_params(1) = NMMcreate_stim_params(1);
        %         sac_stim_params(2:cur_nmods+2) = NMMcreate_stim_params(size(cur_Xsac,2));
        %         mod_signs = [ones(1,cur_nmods+2)];
        %         Xtargets = [1:(cur_nmods+2)];
        %         NL_types = repmat({'lin'},1,cur_nmods+2);
        %         post_gsac_Fmod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
        %         post_gsac_Fmod.mods(1).reg_params = NMMcreate_reg_params();
        %         post_gsac_Fmod = NMMfit_filters(post_gsac_Fmod,cur_Robs,tr_stim,[],[],silent,optim_params);
        %         [post_gsac_Fmod_LL,~,post_Fmod_predrate,~,~,~,nullLL] = NMMmodel_eval(post_gsac_Fmod,cur_Robs,tr_stim);
        %         exc_filters = find(stim_mod_signs == 1);
        %         inh_filters = find(stim_mod_signs == -1);
        %         all_post_gain_filters = [post_gsac_Fmod.mods(3:end).filtK];
        %         if ~isempty(exc_filters)
        %             post_gsac_Egains = mean(all_post_gain_filters(:,exc_filters),2);
        %         else
        %             post_gsac_Egains = nan(length(slags),1);
        %         end
        %         if ~isempty(inh_filters)
        %             post_gsac_Igains = mean(all_post_gain_filters(:,inh_filters),2);
        %         else
        %             post_gsac_Igains = nan(length(slags),1);
        %         end
        %         sacStimProc(cc).gsac_post_fullmod = post_gsac_Fmod;
        %         sacStimProc(cc).gsac_post_Egains = post_gsac_Egains;
        %         sacStimProc(cc).gsac_post_Igains = post_gsac_Igains;
        
        if sum(stim_mod_signs == 1) > 0 && sum(stim_mod_signs == -1) > 0
            g_exc = sum(fgint(:,stim_mod_signs==1),2);
            g_inh = sum(fgint(:,stim_mod_signs==-1),2);
            Xsac_estim = bsxfun(@times,cur_Xsac,g_exc);
            Xsac_istim = bsxfun(@times,cur_Xsac,g_inh);
            clear tr_stim
            tr_stim{1} = [g_exc g_inh];
            tr_stim{2} = cur_Xsac;
            tr_stim{3} = Xsac_estim;
            tr_stim{4} = Xsac_istim;
            sac_stim_params(1) = NMMcreate_stim_params(2);
            sac_stim_params(2:4) = NMMcreate_stim_params([size(Xsac_estim,2)]);
            mod_signs = [1 1 1 1];
            Xtargets = [1 2 3 4];
            NL_types = {'lin','lin','lin','lin'};
            post_gsac_EImod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
            post_gsac_EImod.mods(1).reg_params = NMMcreate_reg_params();
            post_gsac_EImod = NMMfit_filters(post_gsac_EImod,cur_Robs(any_sac_inds),get_Xcell_tInds(tr_stim,any_sac_inds),[],[],silent,optim_params);
            post_gsac_EImod = NMMfit_logexp_spkNL(post_gsac_EImod,cur_Robs(any_sac_inds),get_Xcell_tInds(tr_stim,any_sac_inds));
            [post_gsac_EImod_LL,~,post_EImod_predrate] = NMMmodel_eval(post_gsac_EImod,cur_Robs,tr_stim);
            sacStimProc(cc).gsac_post_EImod = post_gsac_EImod;
            sacStimProc(cc).gsac_post_Egains = post_gsac_EImod.mods(3).filtK;
            sacStimProc(cc).gsac_post_Igains = post_gsac_EImod.mods(4).filtK;
        end
        
        %         g_tot = stimG - cur_GQM.spk_NL_params(1);
        g_tot = stimG;
        Xsac_tot = bsxfun(@times,cur_Xsac,g_tot);
        clear tr_stim
        tr_stim{1} = [g_tot];
        tr_stim{2} = cur_Xsac;
        tr_stim{3} = Xsac_tot;
        clear sac_stim_params
        sac_stim_params(1) = NMMcreate_stim_params(1);
        sac_stim_params(2:3) = NMMcreate_stim_params([size(Xsac_tot,2)]);
        mod_signs = [1 1 1];
        Xtargets = [1 2 3];
        NL_types = {'lin','lin','lin'};
        post_gsac_Smod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
        post_gsac_Smod.mods(1).reg_params = NMMcreate_reg_params();
        post_gsac_Smod = NMMfit_filters(post_gsac_Smod,cur_Robs,tr_stim,[],[],silent);
        post_gsac_Smod = NMMfit_logexp_spkNL(post_gsac_Smod,cur_Robs,tr_stim);
        post_gsac_Smod_LL = NMMmodel_eval(post_gsac_Smod,cur_Robs(any_sac_inds),get_Xcell_tInds(tr_stim,any_sac_inds));
        [post_Smod_LL,~,post_Smod_predrate,~,~,~,nullLL] = NMMmodel_eval( post_gsac_Smod, cur_Robs, tr_stim);
        sacStimProc(cc).gsac_post_singmod = post_gsac_Smod;
        
                
        %% COMPUTE SUBSPAC SAC MODEL
        if fit_subMod
            fprintf('Estimating sac-dep subspace model\n');
            
            X{1} = cur_Xsac;
            X{2} = reshape(bsxfun(@times,cur_Xsac,reshape(ind_Gints,length(cc_uinds),1,[])),length(cc_uinds),[]);
            
            cur_stim_params(1) = NMMcreate_stim_params(length(slags));
            cur_stim_params(2) = NMMcreate_stim_params([length(slags) length(cur_GQM.mods)]);
            
            optim_params.optTol = 1e-5;
            optim_params.progTol = 1e-9;
            modSeq_d2T = 25;
            modSeq_L2 = 0;
            mod_signs = [1 1 1 1];
            Xtargs = [1 2 2 2];
            NL_types = {'lin','lin','quad','quad'};
            reg_params = NMMcreate_reg_params('lambda_d2T',[lambda_d2T repmat(modSeq_d2T,1,3)]','lambda_L2',[lambda_L2 repmat(modSeq_L2,1,3)]');
            init_mod = NMMinitialize_model(cur_stim_params,mod_signs,NL_types,reg_params,Xtargs);
            init_mod = NMMfit_filters(init_mod,cur_Robs(any_sac_inds),get_Xcell_tInds(X,any_sac_inds),[],[],silent,optim_params);
            subspace_mod = NMMfit_filters(init_mod,cur_Robs(any_sac_inds),get_Xcell_tInds(X,any_sac_inds),[],[],silent,optim_params);
            subspace_mod = NMMfit_logexp_spkNL(subspace_mod,cur_Robs(any_sac_inds),get_Xcell_tInds(X,any_sac_inds));
            [subspace_LL,~,subspace_predrate] = NMMmodel_eval(subspace_mod,cur_Robs,X);
            
            cur_filts = reshape([subspace_mod.mods(2:end).filtK],[length(slags) length(cur_GQM.mods) 3]);
            stim_filts = reshape([cur_GQM.mods.filtK],[flen*use_nPix_us length(cur_GQM.mods)]);
            sacdep_filts = nan(3,length(slags),flen*use_nPix_us);
            for jj = 1:3
                sacdep_filts(jj,:,:) = squeeze(cur_filts(:,:,jj))*stim_filts';
            end
            sacStimProc(cc).gsac_phaseDep_subfilt = squeeze(sacdep_filts(1,:,:));
            sacStimProc(cc).gsac_phaseInd_subfilt = squeeze(sqrt(sum(sacdep_filts(2:3,:,:).^2)));
        end
        
        %% COMPUTE MODEL-BASED INFORMATION
        fprintf('Computing saccade-triggered information\n');
        if fitUpstream
            [gainLL,gain_pred_rate] = eval_sacgain_mod( sacGainMod, cur_Robs, all_Xmat_shift, cur_Xsac);
        end
        [sac_fpost_info,sac_spost_info,sac_subpost_info,sac_info,sac_LL,sac_fpost_LL,sac_spost_LL,sac_subpost_LL,sac_nullLL,sac_Nspks,sac_avgrate] = deal(nan(length(slags),1));
        for ii = 1:length(slags)
            temp = find(cur_Xsac(:,ii) == 1);
            sac_avgrate(ii) = mean(cur_Robs(temp));
            cur_avg_rate = sac_avgrate(ii)*ones(size(temp));
            sac_nullLL(ii) = nansum(cur_Robs(temp).*log2(sac_avgrate(ii)) - sac_avgrate(ii));
            sac_Nspks(ii) = sum(cur_Robs(temp));
            
            if fitUpstream
                sac_LL(ii) = nansum(cur_Robs(temp).*log2(gain_pred_rate(temp)) - gain_pred_rate(temp));
                sac_info(ii) = nanmean(gain_pred_rate(temp).*log2(gain_pred_rate(temp)/mean(gain_pred_rate(temp))))/mean(gain_pred_rate(temp));
            end
            
            %             sac_fpost_LL(ii) = nansum(cur_Robs(temp).*log2(post_Fmod_predrate(temp)) - post_Fmod_predrate(temp));
            %             sac_fpost_info(ii) = nanmean(post_Fmod_predrate(temp).*log2(post_Fmod_predrate(temp)/mean(post_Fmod_predrate(temp))))/mean(post_Fmod_predrate(temp));
            sac_spost_LL(ii) = nansum(cur_Robs(temp).*log2(post_Smod_predrate(temp)) - post_Smod_predrate(temp));
            sac_spost_info(ii) = nanmean(post_Smod_predrate(temp).*log2(post_Smod_predrate(temp)/mean(post_Smod_predrate(temp))))/mean(post_Smod_predrate(temp));
            
            if fit_subMod
                sac_subpost_LL(ii) = nansum(cur_Robs(temp).*log2(subspace_predrate(temp)) - subspace_predrate(temp));
                sac_subpost_info(ii) = nanmean(subspace_predrate(temp).*log2(subspace_predrate(temp)/mean(subspace_predrate(temp))))/mean(subspace_predrate(temp));
            end
        end
        
        if fitUpstream
            sacStimProc(cc).gsac_LLinfo = (sac_LL - sac_nullLL)./sac_Nspks;
            sacStimProc(cc).gsac_ov_LLinfo = (gainLL-nullLL)/log(2);
            sacStimProc(cc).gsac_ov_modinfo = mean(gain_pred_rate/mean(gain_pred_rate).*log2(gain_pred_rate/mean(gain_pred_rate)));
            sacStimProc(cc).gsac_modinfo = sac_info;
        end
        
        sacStimProc(cc).gsac_avg_rate = sac_avgrate;
        %         sacStimProc(cc).gsac_fpost_LLinfo = (sac_fpost_LL - sac_nullLL)./sac_Nspks;
        %         sacStimProc(cc).gsac_fpost_ov_LLinfo = (post_gsac_Fmod_LL-nullLL)/log(2);
        %         sacStimProc(cc).gsac_fpost_modinfo = sac_fpost_info;
        %         sacStimProc(cc).gsac_fpost_ov_modinfo = mean(post_Fmod_predrate/mean(post_Fmod_predrate).*log2(post_Fmod_predrate/mean(post_Fmod_predrate)));
        sacStimProc(cc).gsac_spost_LLinfo = (sac_spost_LL - sac_nullLL)./sac_Nspks;
        sacStimProc(cc).gsac_spost_ov_LLinfo = (post_gsac_Smod_LL-nullLL)/log(2);
        sacStimProc(cc).gsac_spost_modinfo = sac_spost_info;
        sacStimProc(cc).gsac_spost_ov_modinfo = mean(post_Smod_predrate/mean(post_Smod_predrate).*log2(post_Smod_predrate/mean(post_Smod_predrate)));
        
        if fit_subMod
            sacStimProc(cc).gsac_sub_LLinfo = (sac_subpost_LL - sac_nullLL)./sac_Nspks;
            sacStimProc(cc).gsac_sub_ov_LLinfo = (subspace_LL-nullLL)/log(2);
            sacStimProc(cc).gsac_sub_modinfo = sac_subpost_info;
            sacStimProc(cc).gsac_sub_ov_modinfo = mean(subspace_predrate/mean(subspace_predrate).*log2(subspace_predrate/mean(subspace_predrate)));
        end
        %% COMPUTE MODEL SEQUENCE
        if fitSTA
            fprintf('Estimating sac-dep STAs\n');
            lag_dep_sta = nan(length(slags),use_nPix_us*flen);
            lag_dep_asta = nan(length(slags),use_nPix_us*flen);
            for ii = 1:length(slags)
                temp = find(cur_Xsac(:,ii) == 1);
                cur_sta = sum(bsxfun(@times,all_Xmat_shift(temp,:),cur_Robs(temp)))/sum(cur_Robs(temp));
                base_a = mean(all_Xmat_shift(temp,:));
                lag_dep_sta(ii,:) = cur_sta-base_a;
                
                cur_sta = sum(bsxfun(@times,abs(all_Xmat_shift(temp,:)),cur_Robs(temp)))/sum(cur_Robs(temp));
                base_a = mean(abs(all_Xmat_shift(temp,:)));
                lag_dep_asta(ii,:) = cur_sta-base_a;
            end
            sacStimProc(cc).gsac_phaseDep_sta = reshape(lag_dep_sta,length(slags),flen, use_nPix_us);
            sacStimProc(cc).gsac_phaseInd_sta = reshape(lag_dep_asta,length(slags),flen,use_nPix_us);
        end
        
        if fitModSeq
            fprintf('Estimating sac-dep model sequence\n');
            
            %compute best time lag over E filts
            [Xinds,Tinds] = meshgrid(1:use_nPix_us,1:flen);
            cur_filts = reshape([cur_GQM.mods(1:end).filtK],[flen use_nPix_us length(cur_GQM.mods)]);
            cur_tfilt = squeeze(mean(std(cur_filts,[],2),3));
            [~,best_lag] = max(cur_tfilt);
            uk = find(Tinds == best_lag);
            newX = reshape(all_Xmat_shift(any_sac_inds,uk),[length(any_sac_inds) 1 use_nPix_us]);
            X{1} = cur_Xsac(any_sac_inds,:);
            X{2} = reshape(bsxfun(@times,newX,cur_Xsac(any_sac_inds,:)),length(any_sac_inds),[]);
            
            cur_stim_params(1) = NMMcreate_stim_params(length(slags));
            cur_stim_params(2) = NMMcreate_stim_params([length(slags) use_nPix_us]);
            
            optim_params.optTol = 1e-4;
            optim_params.progTol = 1e-8;
            modSeq_d2X = 100;
            modSeq_d2T = 20;
            modSeq_init_d2X = 200;
            mod_signs = [1 1 1 1];
            Xtargs = [1 2 2 2];
            NL_types = {'lin','lin','quad','quad'};
            reg_params = NMMcreate_reg_params('lambda_d2X',[0 modSeq_init_d2X modSeq_init_d2X modSeq_init_d2X]','lambda_d2T',[modSeq_d2T 0 0 0]');
            init_mod = NMMinitialize_model(cur_stim_params,mod_signs,NL_types,reg_params,Xtargs);
            init_mod = NMMfit_filters(init_mod,cur_Robs(any_sac_inds),X,[],[],silent,optim_params);
            %             [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(init_mod,cur_Robs(any_sac_inds),X);
            %             vgint = var(gint);
            %             sgint = std(gint);
            %             init_mod = NMMadjust_regularization(init_mod,2:4,'lambda_d2XT',modSeq_d2T./vgint(2:end));
            %             init_mod = NMMfit_filters(init_mod,cur_Robs(any_sac_inds),X,[],[],silent,optim_params);
            
            cur_filts = reshape([init_mod.mods(2:end).filtK],[length(slags) use_nPix_us 3]);
            phasedep_filt = squeeze(cur_filts(:,:,1));
            phaseind_filt = squeeze(sqrt(sum(cur_filts(:,:,2:3).^2,3)));
            sacStimProc(cc).gsac_phaseDep_filt = phasedep_filt;
            sacStimProc(cc).gsac_phaseInd_filt = phaseind_filt;
        end
        
        %% CREATE TENT_BASIS MODEL OF SACCADE-MODULATION
            fprintf('Estimating tent-basis model\n');
            
            Xtick = -(backlag-1/2):(1):(forlag+1/2);
            n_sbins = length(Xtick);
            [LL, penLL, pred_rate, G,gint,fgint] = NMMmodel_eval(cur_GQM, cur_Robs, all_Xmat_shift);
            %         g_tot = G - cur_GQM.spk_NL_params(1);
            
            cur_sac_starts = saccade_start_inds(big_sacs);
            cur_sac_stops = saccade_stop_inds(big_sacs);
            t_since_sac_start = nan(NT,1);
            for ii = 1:length(cur_sac_starts)
                prev_tstart = find(trial_start_inds <= cur_sac_starts(ii),1,'last');
                next_tstop = find(trial_end_inds >= cur_sac_starts(ii),1,'first');
                cur_inds = (cur_sac_starts(ii) - backlag):(cur_sac_starts(ii) + forlag);
                cur_uset = find(cur_inds > trial_start_inds(prev_tstart) & cur_inds < trial_end_inds(next_tstop));
                t_since_sac_start(cur_inds(cur_uset)) = slags(cur_uset);
            end
            
            TB_stim = [t_since_sac_start(cc_uinds) g_tot];
            Ytick = linspace(my_prctile(TB_stim(any_sac_inds,2),0.1),my_prctile(TB_stim(any_sac_inds,2),100-1),n_Gbins);
            TB = TentBasis2D(Xtick, Ytick);
            
            used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
                TB_stim(:,2) >= Ytick(1) & TB_stim(:,2) <= Ytick(end));
            [TB_Xmat,TB_counts] = TB.InputNL2D(TB_stim(used_data,:));
            L2_params = create_L2_params([],[1 n_sbins*n_Gbins],[n_sbins n_Gbins],2,3,[Inf Inf],[0.025 1]);
            TB_fitmod = regGLM_fit(TB_Xmat,cur_Robs(used_data),L2_params,TB_lambda,[],[],silent);
            [LL, penLL, pred_rate, G] = regGLM_eval(TB_fitmod,cur_Robs(used_data),TB_Xmat);
            TB_K = reshape(TB_fitmod.K,n_sbins,n_Gbins)';
            bin_areas = TB.GetBinAreas();
            gsac_TB_dist = TB_counts./bin_areas;
            gsac_TB_dist = gsac_TB_dist'/sum(gsac_TB_dist(:));
            gsac_TB_rate = log(1 + exp(TB_K + TB_fitmod.theta));
            sacStimProc(cc).gsac_TB_rate = gsac_TB_rate;
            
            %INFO CALS
            cur_avg_rate = mean(cur_Robs(used_data));
            marg_gdist = sum(gsac_TB_dist,2);
            marg_sdist = sum(gsac_TB_dist);
            marg_gsacrate = sum(gsac_TB_dist.*gsac_TB_rate)./marg_sdist;
            marg_grate = sum(gsac_TB_dist.*gsac_TB_rate,2)./marg_gdist;
            gsacdep_info = nan(1,n_sac_bins);
            for tt = 1:n_sbins
                gsacdep_info(tt) = sum(gsac_TB_dist(:,tt).*gsac_TB_rate(:,tt).*log2(gsac_TB_rate(:,tt)/marg_gsacrate(tt)))/sum(gsac_TB_dist(:,tt));
            end
            gcumdist = cumsum(marg_gdist)/sum(marg_gdist);
            
            sacStimProc(cc).gsac_ov_TB_info = sum(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate))/cur_avg_rate;
            
            sacStimProc(cc).gsac_TB_avg_rate = marg_gsacrate;
            sacStimProc(cc).gsac_TB_info = gsacdep_info./marg_gsacrate;
            sacStimProc(cc).gsac_TB_gdist = marg_gdist;
            sacStimProc(cc).gsac_TB_grate = marg_grate;
            
            
            
            
            %% NOW DO INFO TIMING CALCS
            shuf_stim = all_shift_stimmat_up;
            shuf_stim(used_inds,:) = all_shift_stimmat_up(used_inds(randi(NT,NT,1)),:);
            shuf_X = create_time_embedding(shuf_stim,stim_params_us);
            shuf_X = shuf_X(used_inds(cc_uinds),use_kInds_up);
            
            sacGainMod = sacStimProc(cc).gsacGainMod;
            [gainLL,gain_pred_rate,G] = eval_sacgain_mod( sacGainMod, cur_Robs, all_Xmat_shift, cur_Xsac);
            gain_LLseq = cur_Robs.*log2(gain_pred_rate)-gain_pred_rate;
            
            [base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,cur_Robs,all_Xmat_shift);
            base_LLseq = cur_Robs.*log2(base_predrate)-base_predrate;
            
            [~,Tinds] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
            Tinds = Tinds(use_kInds_up);
            temp_sparams = NMMcreate_stim_params(1);
            tempmod = NMMinitialize_model(temp_sparams,1,{'lin'});
            
            info_backlag = round(0.3/dt);
            info_forlag = round(0.3/dt);
            info_slags = -info_backlag:info_forlag;
            
            cur_sac_start_inds = saccade_start_inds(big_sacs);
            cur_sac_stop_inds = saccade_stop_inds(big_sacs);
            
            rand_jit = randi(100,length(big_sacs),1);
            rcur_sac_start_inds = cur_sac_start_inds + rand_jit;
            rcur_sac_stop_inds = cur_sac_stop_inds + rand_jit;
            % rcur_sac_stop_inds = rcur_sac_start_inds;
            bad = find(rcur_sac_stop_inds > NT);
            rcur_sac_start_inds(bad) = [];
            rcur_sac_stop_inds(bad) = [];
            
            saccade_stop_trial_inds = all_trialvec(used_inds(cur_sac_stop_inds));
            rsaccade_stop_trial_inds = all_trialvec(used_inds(rcur_sac_stop_inds));
            
            Xsac_end = zeros(NT,length(info_slags));
            Xsac_rend = zeros(NT,length(info_slags));
            for ii = 1:length(info_slags)
                cur_sac_target = cur_sac_stop_inds + info_slags(ii);
                uu = find(cur_sac_target > 1 & cur_sac_target < NT);
                cur_sac_target = cur_sac_target(uu);
                cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_stop_trial_inds(uu)) = [];
                Xsac_end(cur_sac_target,ii) = 1;
                
                cur_sac_target = rcur_sac_stop_inds + info_slags(ii);
                uu = find(cur_sac_target > 1 & cur_sac_target < NT);
                cur_sac_target = cur_sac_target(uu);
                cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= rsaccade_stop_trial_inds(uu)) = [];
                Xsac_rend(cur_sac_target,ii) = 1;
            end
            
            Xsac_end = Xsac_end(cc_uinds,:);
            Xsac_rend = Xsac_rend(cc_uinds,:);
            
            %
            all_during = [];
            for ss = 1:length(cur_sac_start_inds)
                cur_during = (cur_sac_start_inds(ss)):(cur_sac_stop_inds(ss)-1);
                all_during = cat(2,all_during,cur_during);
            end
            is_during = false(length(all_t_axis),1);
            is_during(used_inds(all_during)) = true;
            is_during = repmat(is_during,[1 use_nPix_us]);
            is_during = create_time_embedding(is_during,NMMcreate_stim_params([flen use_nPix_us]));
            is_during = logical(is_during(used_inds(cc_uinds),:));
            
            all_during = [];
            for ss = 1:length(rcur_sac_start_inds)
                cur_during = (rcur_sac_start_inds(ss)):(rcur_sac_stop_inds(ss)-1);
                all_during = cat(2,all_during,cur_during);
            end
            ris_during = false(length(all_t_axis),1);
            ris_during(used_inds(all_during)) = true;
            ris_during = repmat(ris_during,[1 use_nPix_us]);
            ris_during = create_time_embedding(ris_during,NMMcreate_stim_params([flen use_nPix_us]));
            ris_during = logical(ris_during(used_inds(cc_uinds),:));
            
            %%
            for ii = 1:length(info_slags)
                cur_set = find(Xsac_end(:,ii) == 1);
                cur_X = all_Xmat_shift(cur_set,:);
                scramb_lags = find(Tinds-1 > info_slags(ii));
                cur_X(:,scramb_lags) = shuf_X(cur_set,scramb_lags);
                cur_X(is_during(cur_set,:)) = shuf_X(is_during(cur_set,:));
                
                [shuf_LL,shuf_pred_rate,Gshuff] = eval_sacgain_mod( sacGainMod, cur_Robs(cur_set), cur_X, cur_Xsac(cur_set,:));
                tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
                [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
                shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                LL_imp_before(ii) = nansum((gain_LLseq(cur_set) - shuf_LLseq))/sum(cur_Robs(cur_set));
                
                cur_set = find(Xsac_rend(:,ii) == 1);
                cur_X = all_Xmat_shift(cur_set,:);
                scramb_lags = find(Tinds-1 > info_slags(ii));
                cur_X(:,scramb_lags) = shuf_X(cur_set,scramb_lags);
                cur_X(is_during(cur_set,:)) = shuf_X(is_during(cur_set,:));
                
                [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),cur_X);
                Gshuff = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
                tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
                [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
                base_shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                base_LL_imp_before(ii) = nansum((base_LLseq(cur_set) - base_shuf_LLseq))/sum(cur_Robs(cur_set));
            end
            
            %%
            for ii = 1:length(info_slags)
                cur_set = find(Xsac_end(:,ii) == 1);
                cur_X = all_Xmat_shift(cur_set,:);
                scramb_lags = find((Tinds-1) <= info_slags(ii));
                cur_X(:,scramb_lags) = shuf_X(cur_set,scramb_lags);
                
                [shuf_LL,shuf_pred_rate,Gshuff] = eval_sacgain_mod( sacGainMod, cur_Robs(cur_set), cur_X, cur_Xsac(cur_set,:));
                tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
                [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
                shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                
                LL_imp_after(ii) = nansum((gain_LLseq(cur_set) - shuf_LLseq))/sum(cur_Robs(cur_set));
                
                sta_avg_rate(ii) = nanmean(cur_Robs(cur_set));
                
                cur_set = find(Xsac_rend(:,ii) == 1);
                cur_X = all_Xmat_shift(cur_set,:);
                scramb_lags = find((Tinds-1) <= info_slags(ii));
                cur_X(:,scramb_lags) = shuf_X(cur_set,scramb_lags);
                
                [base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),all_Xmat_shift(cur_set,:));
                [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),cur_X);
                Gshuff = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
                tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
                [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
                base_shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                base_LL_imp_after(ii) = nansum((base_LLseq(cur_set) - base_shuf_LLseq))/sum(cur_Robs(cur_set));
                
            end
            
            
            %%
            cur_X = all_Xmat_shift;
            cur_X(is_during) = shuf_X(is_during);
            cur_rX = all_Xmat_shift;
            cur_rX(ris_during) = shuf_X(ris_during);
            %
            for ii = 1:length(info_slags)
                cur_set = find(Xsac_end(:,ii) == 1);
                
                [shuf_LL,shuf_pred_rate,Gshuff] = eval_sacgain_mod( sacGainMod, cur_Robs(cur_set), cur_X(cur_set,:), cur_Xsac(cur_set,:));
                tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
                [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
                shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                
                LL_imp_during(ii) = nansum((gain_LLseq(cur_set) - shuf_LLseq))/sum(cur_Robs(cur_set));
                
                cur_set = find(Xsac_rend(:,ii) == 1);
                [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),cur_rX(cur_set,:));
                Gshuff = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
                tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
                [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
                base_shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                base_LL_imp_during(ii) = nansum((base_LLseq(cur_set) - base_shuf_LLseq))/sum(cur_Robs(cur_set));
                
            end
            
            %%
            sacStimProc(cc).info_slags = info_slags;
            sacStimProc(cc).gsac_Ltime_before = LL_imp_before;
            sacStimProc(cc).gsac_Ltime_before_control = base_LL_imp_before;
            sacStimProc(cc).gsac_Ltime_after = LL_imp_after;
            sacStimProc(cc).gsac_Ltime_after_control = base_LL_imp_after;
            sacStimProc(cc).gsac_Ltime_during = LL_imp_during;
            sacStimProc(cc).gsac_Ltime_during_control = base_LL_imp_during;
            

            %% FOR MSACS
        if fitMsacs
            any_sac_inds = find(any(Xmsac(cc_uinds,:) > 0,2));
            cur_Xsac = Xmsac(cc_uinds,:);
            
            %% FIT UPSTREAM STIM-MODULATION
            if fitUpstream
                fprintf('Estimating m-sac upstream filter\n');
                Xsac_mat = Xmsac(cc_uinds(any_sac_inds),:);
                sacGainMod = fit_sacgain_model(cur_GQM,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:),Xsac_mat,lambda_d2T,lambda_L2);
                sacStimProc(cc).msacGainMod = sacGainMod;
            end
            
            %% FIT POST-INTEGRATION GAIN
            fprintf('Fitting post-filter models\n');
            
            sac_reg_params = NMMcreate_reg_params('lambda_d2T',lambda_d2T,'lambda_L2',lambda_L2,'boundary_conds',[0 0 0]);
            
            if sum(stim_mod_signs == 1) > 0 && sum(stim_mod_signs == -1) > 0
                g_exc = sum(fgint(:,stim_mod_signs==1),2);
                g_inh = sum(fgint(:,stim_mod_signs==-1),2);
                Xsac_estim = bsxfun(@times,cur_Xsac,g_exc);
                Xsac_istim = bsxfun(@times,cur_Xsac,g_inh);
                clear tr_stim
                tr_stim{1} = [g_exc g_inh];
                tr_stim{2} = cur_Xsac;
                tr_stim{3} = Xsac_estim;
                tr_stim{4} = Xsac_istim;
                sac_stim_params(1) = NMMcreate_stim_params(2);
                sac_stim_params(2:4) = NMMcreate_stim_params([size(Xsac_estim,2)]);
                mod_signs = [1 1 1 1];
                Xtargets = [1 2 3 4];
                NL_types = {'lin','lin','lin','lin'};
                post_msac_EImod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
                post_msac_EImod.mods(1).reg_params = NMMcreate_reg_params();
                post_msac_EImod = NMMfit_filters(post_msac_EImod,cur_Robs,tr_stim,[],[],silent,optim_params);
                [post_msac_EImod_LL,~,post_EImod_predrate] = NMMmodel_eval(post_msac_EImod,cur_Robs(any_sac_inds),get_Xcell_tInds(tr_stim,any_sac_inds));
                sacStimProc(cc).msac_post_EImod = post_msac_EImod;
                sacStimProc(cc).msac_post_Egains = post_msac_EImod.mods(3).filtK;
                sacStimProc(cc).msac_post_Igains = post_msac_EImod.mods(4).filtK;
            end
            
%             g_tot = stimG - cur_GQM.spk_NL_params(1);
            g_tot = stimG;
            Xsac_tot = bsxfun(@times,cur_Xsac,g_tot);
            clear tr_stim
            tr_stim{1} = [g_tot];
            tr_stim{2} = cur_Xsac;
            tr_stim{3} = Xsac_tot;
            clear sac_stim_params
            sac_stim_params(1) = NMMcreate_stim_params(1);
            sac_stim_params(2:3) = NMMcreate_stim_params([size(Xsac_tot,2)]);
            mod_signs = [1 1 1];
            Xtargets = [1 2 3];
            NL_types = {'lin','lin','lin'};
            post_msac_Smod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
            post_msac_Smod = NMMfit_filters(post_msac_Smod,cur_Robs,tr_stim,[],[],silent);
            post_msac_Smod_LL = NMMmodel_eval(post_msac_Smod,cur_Robs(any_sac_inds),get_Xcell_tInds(tr_stim,any_sac_inds));
            [post_Smod_LL,~,post_Smod_predrate] = NMMmodel_eval( post_msac_Smod, cur_Robs, tr_stim);
            sacStimProc(cc).msac_post_singmod = post_msac_Smod;
            
            
            %% COMPUTE SUBSPAC SAC MODEL
            if fit_subMod
                fprintf('Estimating sac-dep subspace model\n');
                
                X{1} = cur_Xsac(any_sac_inds,:);
                X{2} = reshape(bsxfun(@times,cur_Xsac(any_sac_inds,:),reshape(ind_Gints(any_sac_inds,:),length(any_sac_inds),1,[])),length(any_sac_inds),[]);
                
                cur_stim_params(1) = NMMcreate_stim_params(length(slags));
                cur_stim_params(2) = NMMcreate_stim_params([length(slags) length(cur_GQM.mods)]);
                
                optim_params.optTol = 1e-4;
                optim_params.progTol = 1e-8;
                modSeq_d2T = 25;
                modSeq_L2 = 0;
                mod_signs = [1 1 1 1];
                Xtargs = [1 2 2 2];
                NL_types = {'lin','lin','quad','quad'};
                reg_params = NMMcreate_reg_params('lambda_d2T',[lambda_d2T repmat(modSeq_d2T,1,3)]','lambda_L2',[lambda_L2 repmat(modSeq_L2,1,3)]');
                init_mod = NMMinitialize_model(cur_stim_params,mod_signs,NL_types,reg_params,Xtargs);
                init_mod = NMMfit_filters(init_mod,cur_Robs(any_sac_inds),X,[],[],silent,optim_params);
                %         [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(init_mod,cur_Robs(any_sac_inds),X);
                %         vgint = var(gint);
                %         sgint = std(gint);
                %         init_mod = NMMadjust_regularization(init_mod,2:4,'lambda_d2T',modSeq_d2T./vgint(2:end));
                subspace_mod = NMMfit_filters(init_mod,cur_Robs(any_sac_inds),X,[],[],silent,optim_params);
                [subspace_LL,~,subspace_predrate] = NMMmodel_eval(subspace_mod,cur_Robs(any_sac_inds),X);
                
                cur_filts = reshape([subspace_mod.mods(2:end).filtK],[length(slags) length(cur_GQM.mods) 3]);
                stim_filts = reshape([cur_GQM.mods.filtK],[flen*use_nPix_us length(cur_GQM.mods)]);
                sacdep_filts = nan(3,length(slags),flen*use_nPix_us);
                for jj = 1:3
                    sacdep_filts(jj,:,:) = squeeze(cur_filts(:,:,jj))*stim_filts';
                end
                sacStimProc(cc).msac_phaseDep_subfilt = squeeze(sacdep_filts(1,:,:));
                sacStimProc(cc).msac_phaseInd_subfilt = squeeze(sqrt(sum(sacdep_filts(2:3,:,:).^2)));
            end
            
            %% COMPUTE MODEL-BASED INFORMATION
            fprintf('Computing saccade-triggered information\n');
            if fitUpstream
                [gainLL,gain_pred_rate] = eval_sacgain_mod( sacGainMod, cur_Robs, all_Xmat_shift, cur_Xsac);
            end
            [sac_fpost_info,sac_spost_info,sac_subpost_info,sac_info,sac_LL,sac_fpost_LL,sac_spost_LL,sac_subpost_LL,sac_nullLL,sac_Nspks,sac_avgrate] = deal(nan(length(slags),1));
            for ii = 1:length(slags)
                temp = find(cur_Xsac(:,ii) == 1);
                sac_avgrate(ii) = mean(cur_Robs(temp));
                cur_avg_rate = sac_avgrate(ii)*ones(size(temp));
                sac_nullLL(ii) = nansum(cur_Robs(temp).*log2(sac_avgrate(ii)) - sac_avgrate(ii));
                sac_Nspks(ii) = sum(cur_Robs(temp));
                
                temp2 = find(cur_Xsac(any_sac_inds,ii) == 1);
                if fitUpstream
                    sac_LL(ii) = nansum(cur_Robs(temp).*log2(gain_pred_rate(temp)) - gain_pred_rate(temp));
                    sac_info(ii) = nanmean(gain_pred_rate(temp).*log2(gain_pred_rate(temp)/mean(gain_pred_rate(temp))))/mean(gain_pred_rate(temp));
                end
                
%                 sac_fpost_LL(ii) = nansum(cur_Robs(temp).*log2(post_Fmod_predrate(temp)) - post_Fmod_predrate(temp));
%                 sac_fpost_info(ii) = nanmean(post_Fmod_predrate(temp).*log2(post_Fmod_predrate(temp)/mean(post_Fmod_predrate(temp))))/mean(post_Fmod_predrate(temp));
                sac_spost_LL(ii) = nansum(cur_Robs(temp).*log2(post_Smod_predrate(temp)) - post_Smod_predrate(temp));
                sac_spost_info(ii) = nanmean(post_Smod_predrate(temp).*log2(post_Smod_predrate(temp)/mean(post_Smod_predrate(temp))))/mean(post_Smod_predrate(temp));
                
                if fit_subMod
                    sac_subpost_LL(ii) = nansum(cur_Robs(temp).*log2(subspace_predrate(temp2)) - subspace_predrate(temp2));
                    sac_subpost_info(ii) = nanmean(subspace_predrate(temp2).*log2(subspace_predrate(temp2)/mean(subspace_predrate(temp2))))/mean(subspace_predrate(temp2));
                end
            end
            
            if fitUpstream
                sacStimProc(cc).msac_LLinfo = (sac_LL - sac_nullLL)./sac_Nspks;
                sacStimProc(cc).msac_ov_LLinfo = (gainLL-nullLL)/log(2);
                sacStimProc(cc).msac_ov_modinfo = mean(gain_pred_rate/mean(gain_pred_rate).*log2(gain_pred_rate/mean(gain_pred_rate)));
                sacStimProc(cc).msac_modinfo = sac_info;
            end
            
            sacStimProc(cc).msac_avg_rate = sac_avgrate;
            
%             sacStimProc(cc).msac_fpost_LLinfo = (sac_fpost_LL - sac_nullLL)./sac_Nspks;
%             sacStimProc(cc).msac_fpost_ov_LLinfo = (post_msac_Fmod_LL-nullLL)/log(2);
%             sacStimProc(cc).msac_fpost_modinfo = sac_fpost_info;
%             sacStimProc(cc).msac_fpost_ov_modinfo = mean(post_Fmod_predrate/mean(post_Fmod_predrate).*log2(post_Fmod_predrate/mean(post_Fmod_predrate)));
            sacStimProc(cc).msac_spost_LLinfo = (sac_spost_LL - sac_nullLL)./sac_Nspks;
            sacStimProc(cc).msac_spost_ov_LLinfo = (post_msac_Smod_LL-nullLL)/log(2);
            sacStimProc(cc).msac_spost_modinfo = sac_spost_info;
            sacStimProc(cc).msac_spost_ov_modinfo = mean(post_Smod_predrate/mean(post_Smod_predrate).*log2(post_Smod_predrate/mean(post_Smod_predrate)));
            
            if fit_subMod
                sacStimProc(cc).msac_sub_LLinfo = (sac_subpost_LL - sac_nullLL)./sac_Nspks;
                sacStimProc(cc).msac_sub_ov_LLinfo = (subspace_LL-nullLL)/log(2);
                sacStimProc(cc).msac_sub_modinfo = sac_subpost_info;
                sacStimProc(cc).msac_sub_ov_modinfo = mean(subspace_predrate/mean(subspace_predrate).*log2(subspace_predrate/mean(subspace_predrate)));
            end
            %% CREATE TENT_BASIS MODEL OF SACCADE-MODULATION
            fprintf('Estimating m-sac tent-basis model\n');
            
            Xtick = -(backlag-1/2):(1):(forlag+1/2);
            n_sbins = length(Xtick);
            [LL, penLL, pred_rate, G,gint,fgint] = NMMmodel_eval(cur_GQM, cur_Robs, all_Xmat_shift);
            g_tot = G - cur_GQM.spk_NL_params(1);
            
            cur_sac_starts = saccade_start_inds(big_sacs);
            cur_sac_stops = saccade_stop_inds(big_sacs);
            t_since_sac_start = nan(NT,1);
            for ii = 1:length(cur_sac_starts)
                prev_tstart = find(trial_start_inds <= cur_sac_starts(ii),1,'last');
                next_tstop = find(trial_end_inds >= cur_sac_starts(ii),1,'first');
                cur_inds = (cur_sac_starts(ii) - backlag):(cur_sac_starts(ii) + forlag);
                cur_uset = find(cur_inds > trial_start_inds(prev_tstart) & cur_inds < trial_end_inds(next_tstop));
                t_since_sac_start(cur_inds(cur_uset)) = slags(cur_uset);
            end
            
            TB_stim = [t_since_sac_start(cc_uinds) g_tot];
            Ytick = linspace(my_prctile(TB_stim(:,2),0.1),my_prctile(TB_stim(:,2),100-1),n_Gbins);
            TB = TentBasis2D(Xtick, Ytick);
            
            used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
                TB_stim(:,2) >= Ytick(1) & TB_stim(:,2) <= Ytick(end));
            [TB_Xmat,TB_counts] = TB.InputNL2D(TB_stim(used_data,:));
            L2_params = create_L2_params([],[1 n_sbins*n_Gbins],[n_sbins n_Gbins],2,3,[Inf Inf],[0.025 1]);
            TB_fitmod = regGLM_fit(TB_Xmat,cur_Robs(used_data),L2_params,TB_lambda,[],[],1);
            [LL, penLL, pred_rate, G] = regGLM_eval(TB_fitmod,cur_Robs(used_data),TB_Xmat);
            TB_K = reshape(TB_fitmod.K,n_sbins,n_Gbins)';
            bin_areas = TB.GetBinAreas();
            msac_TB_dist = TB_counts./bin_areas;
            msac_TB_dist = msac_TB_dist'/sum(msac_TB_dist(:));
            msac_TB_rate = log(1 + exp(TB_K + TB_fitmod.theta));
            sacStimProc(cc).msac_TB_rate = msac_TB_rate;
            
            %INFO CALS
            cur_avg_rate = mean(cur_Robs(used_data));
            marg_gdist = sum(msac_TB_dist,2);
            marg_sdist = sum(msac_TB_dist);
            marg_msacrate = sum(msac_TB_dist.*msac_TB_rate)./marg_sdist;
            marg_grate = sum(msac_TB_dist.*msac_TB_rate,2)./marg_gdist;
            msacdep_info = nan(1,n_sac_bins);
            for tt = 1:n_sbins
                msacdep_info(tt) = sum(msac_TB_dist(:,tt).*msac_TB_rate(:,tt).*log2(msac_TB_rate(:,tt)/marg_msacrate(tt)))/sum(msac_TB_dist(:,tt));
            end
            gcumdist = cumsum(marg_gdist)/sum(marg_gdist);
            
            sacStimProc(cc).msac_ov_TB_info = sum(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate))/cur_avg_rate;
            
            sacStimProc(cc).msac_TB_avg_rate = marg_msacrate;
            sacStimProc(cc).msac_TB_info = msacdep_info./marg_msacrate;
            sacStimProc(cc).msac_TB_gdist = marg_gdist;
            sacStimProc(cc).msac_TB_grate = marg_grate;
            
            
            
            
            
            
            
            
            
                        %% NOW DO INFO TIMING CALCS
            shuf_stim = all_shift_stimmat_up;
            shuf_stim(used_inds,:) = all_shift_stimmat_up(used_inds(randi(NT,NT,1)),:);
            shuf_X = create_time_embedding(shuf_stim,stim_params_us);
            shuf_X = shuf_X(used_inds(cc_uinds),use_kInds_up);
            
            sacGainMod = sacStimProc(cc).msacGainMod;
            [gainLL,gain_pred_rate,G] = eval_sacgain_mod( sacGainMod, cur_Robs, all_Xmat_shift, cur_Xsac);
            gain_LLseq = cur_Robs.*log2(gain_pred_rate)-gain_pred_rate;
            
            [base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,cur_Robs,all_Xmat_shift);
            base_LLseq = cur_Robs.*log2(base_predrate)-base_predrate;
            
            [~,Tinds] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
            Tinds = Tinds(use_kInds_up);
            temp_sparams = NMMcreate_stim_params(1);
            tempmod = NMMinitialize_model(temp_sparams,1,{'lin'});
            
            info_backlag = round(0.3/dt);
            info_forlag = round(0.3/dt);
            info_slags = -info_backlag:info_forlag;
            
            cur_sac_start_inds = saccade_start_inds(micro_sacs);
            cur_sac_stop_inds = saccade_stop_inds(micro_sacs);
            
            rand_jit = randi(100,length(micro_sacs),1);
            rcur_sac_start_inds = cur_sac_start_inds + rand_jit;
            rcur_sac_stop_inds = cur_sac_stop_inds + rand_jit;
            % rcur_sac_stop_inds = rcur_sac_start_inds;
            bad = find(rcur_sac_stop_inds > NT);
            rcur_sac_start_inds(bad) = [];
            rcur_sac_stop_inds(bad) = [];
            
            saccade_stop_trial_inds = all_trialvec(used_inds(cur_sac_stop_inds));
            rsaccade_stop_trial_inds = all_trialvec(used_inds(rcur_sac_stop_inds));
            
            Xsac_end = zeros(NT,length(info_slags));
            Xsac_rend = zeros(NT,length(info_slags));
            for ii = 1:length(info_slags)
                cur_sac_target = cur_sac_stop_inds + info_slags(ii);
                uu = find(cur_sac_target > 1 & cur_sac_target < NT);
                cur_sac_target = cur_sac_target(uu);
                cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_stop_trial_inds(uu)) = [];
                Xsac_end(cur_sac_target,ii) = 1;
                
                cur_sac_target = rcur_sac_stop_inds + info_slags(ii);
                uu = find(cur_sac_target > 1 & cur_sac_target < NT);
                cur_sac_target = cur_sac_target(uu);
                cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= rsaccade_stop_trial_inds(uu)) = [];
                Xsac_rend(cur_sac_target,ii) = 1;
            end
            
            Xsac_end = Xsac_end(cc_uinds,:);
            Xsac_rend = Xsac_rend(cc_uinds,:);
            
            %
            all_during = [];
            for ss = 1:length(cur_sac_start_inds)
                cur_during = (cur_sac_start_inds(ss)):(cur_sac_stop_inds(ss)-1);
                all_during = cat(2,all_during,cur_during);
            end
            is_during = false(length(all_t_axis),1);
            is_during(used_inds(all_during)) = true;
            is_during = repmat(is_during,[1 use_nPix_us]);
            is_during = create_time_embedding(is_during,NMMcreate_stim_params([flen use_nPix_us]));
            is_during = logical(is_during(used_inds(cc_uinds),:));
            
            all_during = [];
            for ss = 1:length(rcur_sac_start_inds)
                cur_during = (rcur_sac_start_inds(ss)):(rcur_sac_stop_inds(ss)-1);
                all_during = cat(2,all_during,cur_during);
            end
            ris_during = false(length(all_t_axis),1);
            ris_during(used_inds(all_during)) = true;
            ris_during = repmat(ris_during,[1 use_nPix_us]);
            ris_during = create_time_embedding(ris_during,NMMcreate_stim_params([flen use_nPix_us]));
            ris_during = logical(ris_during(used_inds(cc_uinds),:));
            
            %%
            for ii = 1:length(info_slags)
                cur_set = find(Xsac_end(:,ii) == 1);
                cur_X = all_Xmat_shift(cur_set,:);
                scramb_lags = find(Tinds-1 > info_slags(ii));
                cur_X(:,scramb_lags) = shuf_X(cur_set,scramb_lags);
                cur_X(is_during(cur_set,:)) = shuf_X(is_during(cur_set,:));
                
                [shuf_LL,shuf_pred_rate,Gshuff] = eval_sacgain_mod( sacGainMod, cur_Robs(cur_set), cur_X, cur_Xsac(cur_set,:));
                tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
                [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
                shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                LL_imp_before(ii) = nansum((gain_LLseq(cur_set) - shuf_LLseq))/sum(cur_Robs(cur_set));
                
                cur_set = find(Xsac_rend(:,ii) == 1);
                cur_X = all_Xmat_shift(cur_set,:);
                scramb_lags = find(Tinds-1 > info_slags(ii));
                cur_X(:,scramb_lags) = shuf_X(cur_set,scramb_lags);
                cur_X(is_during(cur_set,:)) = shuf_X(is_during(cur_set,:));
                
                [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),cur_X);
                Gshuff = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
                tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
                [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
                base_shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                base_LL_imp_before(ii) = nansum((base_LLseq(cur_set) - base_shuf_LLseq))/sum(cur_Robs(cur_set));
            end
            
            %%
            for ii = 1:length(info_slags)
                cur_set = find(Xsac_end(:,ii) == 1);
                cur_X = all_Xmat_shift(cur_set,:);
                scramb_lags = find((Tinds-1) <= info_slags(ii));
                cur_X(:,scramb_lags) = shuf_X(cur_set,scramb_lags);
                
                [shuf_LL,shuf_pred_rate,Gshuff] = eval_sacgain_mod( sacGainMod, cur_Robs(cur_set), cur_X, cur_Xsac(cur_set,:));
                tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
                [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
                shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                
                LL_imp_after(ii) = nansum((gain_LLseq(cur_set) - shuf_LLseq))/sum(cur_Robs(cur_set));
                
                sta_avg_rate(ii) = nanmean(cur_Robs(cur_set));
                
                cur_set = find(Xsac_rend(:,ii) == 1);
                cur_X = all_Xmat_shift(cur_set,:);
                scramb_lags = find((Tinds-1) <= info_slags(ii));
                cur_X(:,scramb_lags) = shuf_X(cur_set,scramb_lags);
                
                [base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),all_Xmat_shift(cur_set,:));
                [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),cur_X);
                Gshuff = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
                tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
                [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
                base_shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                base_LL_imp_after(ii) = nansum((base_LLseq(cur_set) - base_shuf_LLseq))/sum(cur_Robs(cur_set));
                
            end
            
            
            %%
            cur_X = all_Xmat_shift;
            cur_X(is_during) = shuf_X(is_during);
            cur_rX = all_Xmat_shift;
            cur_rX(ris_during) = shuf_X(ris_during);
            %
            for ii = 1:length(info_slags)
                cur_set = find(Xsac_end(:,ii) == 1);
                
                [shuf_LL,shuf_pred_rate,Gshuff] = eval_sacgain_mod( sacGainMod, cur_Robs(cur_set), cur_X(cur_set,:), cur_Xsac(cur_set,:));
                tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
                [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
                shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                
                LL_imp_during(ii) = nansum((gain_LLseq(cur_set) - shuf_LLseq))/sum(cur_Robs(cur_set));
                
                cur_set = find(Xsac_rend(:,ii) == 1);
                [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,cur_Robs(cur_set),cur_rX(cur_set,:));
                Gshuff = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
                tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
                [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
                base_shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                base_LL_imp_during(ii) = nansum((base_LLseq(cur_set) - base_shuf_LLseq))/sum(cur_Robs(cur_set));
                
            end
            
            %%
            sacStimProc(cc).info_slags = info_slags;
            sacStimProc(cc).msac_Ltime_before = LL_imp_before;
            sacStimProc(cc).msac_Ltime_before_control = base_LL_imp_before;
            sacStimProc(cc).msac_Ltime_after = LL_imp_after;
            sacStimProc(cc).msac_Ltime_after_control = base_LL_imp_after;
            sacStimProc(cc).msac_Ltime_during = LL_imp_during;
            sacStimProc(cc).msac_Ltime_during_control = base_LL_imp_during;
            

        end
    else
        sacStimProc(cc).used = false;
    end
end

%%
anal_dir = ['/home/james/Analysis/bruce/' Expt_name '/sac_mod/'];
fname = 'sacStimProc_v4';
% fname = 'sacStimProc_sta';
% fname = 'sacStimProc_mua';
if strcmp(rec_type,'UA') && bar_ori == 90
    fname = [fname '_vbars'];
end
if fit_unCor
    fname = [fname '_unCor'];
end
cd(anal_dir)
save(fname,'targs','slags','dt','sacStimProc');

