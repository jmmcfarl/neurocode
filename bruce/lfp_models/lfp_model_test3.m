%
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements//');
addpath('~/James_scripts/bruce/processing/');
addpath('~/James_scripts/bruce/saccade_modulation/');
addpath('~/James_scripts/TentBasis2D/');

global Expt_name bar_ori use_MUA


Expt_name = 'M277';
% Expt_name = 'G086';
use_MUA = false;
bar_ori = 0; %bar orientation to use (only for UA recs)


fit_unCor = false;
fit_subMod = false;
fitModSeq = false;
fitUpstream = false;
fitSTA = false;
fitMsacs = true;

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
        case 297
            if ~ismember(bar_ori,[0 90])
                error('M297 is either 0 or 90 deg bars');
            end
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
mod_data_dir = ['~/Analysis/bruce/' Expt_name '/models'];

% et_mod_data_name = 'full_eyetrack_initmods';
% et_anal_name = 'full_eyetrack';
et_mod_data_name = 'full_eyetrack_initmods_Rinit';
et_anal_name = 'full_eyetrack_Rinit';
mod_data_name = 'corrected_models2';

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

et_mod_data_name = [et_mod_data_name sprintf('_ori%d',bar_ori)];
et_anal_name = [et_anal_name sprintf('_ori%d',bar_ori)];
mod_data_name = [mod_data_name sprintf('_ori%d',bar_ori)];

%%

flen = 15;
spatial_usfac = 2;


%these recs have larger bar widths
if ismember(Expt_num,[287 289 294])
    use_nPix = 15;
    spatial_usfac = 4;
elseif ismember(Expt_num,[296 297])
    use_nPix = 22;
    spatial_usfac = 2;
end

min_trial_dur = 0.75;

stim_fs = 100; %in Hz
dt = 0.01;

Fr = 1;

backlag = round(0.3/dt);
forlag = round(0.5/dt);
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
    if all_mod_SU(ss) > 0
        su_probe_ind = find(Clust_data.SU_numbers == all_mod_SUnum(ss));
        Robs_mat(:,ss) = all_binned_sua(used_inds,su_probe_ind);
    elseif ~isnan(all_mod_SU(ss))
        Robs_mat(:,ss) = all_binned_mua(used_inds,ss);
    end
end

%% Create set of TR and XV trials
use_trials = unique(all_trialvec(used_inds));
if ~isnan(rpt_seed)
    rpt_trials = find(all_trial_Se == rpt_seed);
    use_trials(ismember(use_trials,rpt_trials)) = [];
end
nuse_trials = length(use_trials);
xv_frac = 0.2;
% if ~isnan(rpt_seed)
%     xv_trials = find(all_trial_Se==rpt_seed);
%     n_xv_trials = length(xv_trials);
% else
n_xv_trials = round(xv_frac*nuse_trials);
xv_trials = randperm(nuse_trials);
xv_trials(n_xv_trials+1:end) = [];
xv_trials = use_trials(xv_trials);
% end
tr_trials = setdiff(use_trials,xv_trials);
n_tr_trials = length(tr_trials);
fprintf('Initializing models with %d training trials and %d xval trials\n',n_tr_trials,n_xv_trials);

tr_inds = find(ismember(all_trialvec(used_inds),tr_trials));
xv_inds = find(ismember(all_trialvec(used_inds),xv_trials));
all_inds = union(tr_inds,xv_inds);

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
% if ~fit_unCor
%     for i=1:NT
%         best_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
%     end
% end
all_Xmat_shift = create_time_embedding(best_shift_stimmat_up,stim_params_us);
all_Xmat_shift = all_Xmat_shift(used_inds,use_kInds_up);

%% LOAD LFP DATA
cd(data_dir)
if strcmp(rec_type,'LP')
    Fs = 1000;
    dsf = 5;
    Fsd = Fs/dsf;
    niqf = Fs/2;
    new_niqf = Fsd/2;
    [bb,aa] = butter(2,0.8*new_niqf/niqf);
    use_lfps = 1:2:n_probes;
elseif strcmp(rec_type,'UA')
    Fs = 400.0032;
    dsf = 2;
    Fsd = Fs/dsf;
    niqf = Fs/2;
    new_niqf = Fsd/2;
    [bb,aa] = butter(2,0.8*new_niqf/niqf);
    use_lfps = 1:n_probes;
    %     use_lfps = [28 SU_probes];
end

%wavelet parameters
nwfreqs = 25;
min_freq = 1.5; max_freq = 50;
min_scale = 1/max_freq*Fsd;
max_scale = 1/min_freq*Fsd;
wavetype = 'cmor1-1';
scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
wfreqs = scal2frq(scales,wavetype,1/Fsd);


full_lfps = [];
full_cwts = [];
full_lfp_taxis = [];
cur_toffset = 0;
ublock_set = 1:length(cur_block_set);
if Expt_num == 275
    bad_blocks = 9;
    ublock_set(bad_blocks) = [];
end
for ee = ublock_set
    
    fprintf('Loading LFPs, Expt %d of %d\n',ee,length(cur_block_set));
    
    if strcmp(rec_type,'LP')
        fname = sprintf('lemM%dA.%d.lfp.mat',Expt_num,cur_block_set(ee));
        load(fname);
        
        tlens = arrayfun(@(X) length(X.ftime),LFP.Trials);
        bad_trials = find(tlens == 0);
        LFP.Trials(bad_trials) = [];
        n_trials(ee) = length(LFP.Trials);
        lfp_trial_starts = [LFP.Trials(:).ftime]/1e4;
        lfp_trial_ends = [LFP.Trials(:).End]/1e4;
        expt_lfp_t_axis = [];
        expt_lfps = [];
        for tt = 1:n_trials(ee)
            %         tt
            cur_npts = size(LFP.Trials(tt).LFP,1);
            cur_t_end(tt) = lfp_trial_starts(tt)+(cur_npts-1)/Fs;
            cur_t_axis = (lfp_trial_starts(tt):1/Fs:cur_t_end(tt)) + cur_toffset;
            
            if ~isempty(expt_lfp_t_axis)
                cur_sp = find(cur_t_axis > max(expt_lfp_t_axis),1,'first');
            else
                cur_sp = 1;
            end
            cur_t_axis = cur_t_axis(cur_sp:end);
            if length(cur_t_axis) > 50 & size(LFP.Trials(tt).LFP,2) == n_probes
                cur_LFP = double([LFP.Trials(tt).LFP]);
                cur_LFP = cur_LFP(cur_sp:end,use_lfps);
                if dsf > 1
                    cur_LFP = filtfilt(bb,aa,cur_LFP);
                    cur_LFP = downsample(cur_LFP,dsf);
                    cur_t_axis = downsample(cur_t_axis,dsf);
                end
                if size(cur_LFP,2) == length(use_lfps)
                    expt_lfp_t_axis = [expt_lfp_t_axis; cur_t_axis(:)];
                    expt_lfps = [expt_lfps; cur_LFP];
                end
            end
        end
        expt_cwt = nan(size(expt_lfps,1),length(wfreqs),length(use_lfps));
        for cc = 1:length(use_lfps)
            expt_cwt(:,:,cc) = cwt(expt_lfps(:,cc),scales,'cmor1-1')';
        end
    else
        lfp_fname = sprintf('Expt%d_LFP.mat',cur_block_set(ee));
        load(lfp_fname);
        
        cur_lfps = bsxfun(@times,double(lfp_mat(:,use_lfps)),lfp_int2V(use_lfps)');
        if dsf > 1
            cur_lfps = filtfilt(bb,aa,cur_lfps);
            expt_lfps = downsample(cur_lfps,dsf);
            expt_lfp_t_axis = downsample(lfp_t_ax',dsf);
        end
    end
    
    cur_uset = find(all_blockvec == ee);
    if ~isempty(cur_uset)
        uinds = find(expt_lfp_t_axis >= all_t_axis(cur_uset(1)) & expt_lfp_t_axis <= all_t_axis(cur_uset(end)));
        full_lfps = cat(1,full_lfps,expt_lfps(uinds,:));
        full_cwts = cat(1,full_cwts,expt_cwt(uinds,:,:));
        full_lfp_taxis = cat(1,full_lfp_taxis,expt_lfp_t_axis(uinds));
    end
    cur_toffset = trial_toffset(ee);
end

%%
interp_lfps = interp1(full_lfp_taxis,full_lfps,all_t_axis);
% interp_lfps = interp_lfps(used_inds,:);
interp_lfps = nanzscore(interp_lfps);

%%
interp_cwt_real = interp1(full_lfp_taxis,real(full_cwts),all_t_axis);
interp_cwt_imag = interp1(full_lfp_taxis,imag(full_cwts),all_t_axis);
interp_cwt_mag = sqrt(interp_cwt_real.^2 + interp_cwt_imag.^2);

interp_cwt_real = bsxfun(@rdivide,interp_cwt_real,nanstd(interp_cwt_mag));
interp_cwt_imag = bsxfun(@rdivide,interp_cwt_imag,nanstd(interp_cwt_mag));

%%
test_trials = find(all_trial_Se == rpt_seed);
Ntrials = length(test_trials);
beg_buff = 0; %number of bins from beginning of trial to exclude
end_buff = 0;

fullRobs_ms = bsxfun(@minus,Robs_mat,nanmean(Robs_mat(ismember(all_trialvec(used_inds),test_trials),:)));
Nframes = 375;
Nunits = size(Robs_mat,2);
tbt_Robs = nan(Ntrials,Nframes,Nunits);
tbt_LFP_real = nan(Ntrials,Nframes,length(wfreqs),length(use_lfps));
tbt_LFP_imag = nan(Ntrials,Nframes,length(wfreqs),length(use_lfps));
tbt_LFPs = nan(Ntrials,Nframes,length(use_lfps));
for ii = 1:Ntrials
   cur_inds = find(all_trialvec(used_inds) == test_trials(ii));
   tbt_Robs(ii,:,:) = fullRobs_ms(cur_inds,:);
   cur_inds = used_inds(all_trialvec(used_inds) == test_trials(ii));
   tbt_LFP_real(ii,:,:,:) = interp_cwt_real(cur_inds,:,:);   
   tbt_LFP_imag(ii,:,:,:) = interp_cwt_imag(cur_inds,:,:);   
   tbt_LFPs(ii,:,:) = interp_lfps(cur_inds,:);   
end

% tbt_Robs = bsxfun(@minus,tbt_Robs,nanmean(tbt_Robs));

tbt_LFP_amp = sqrt(tbt_LFP_real.^2 + tbt_LFP_imag.^2);
tbt_LFP_real = tbt_LFP_real./tbt_LFP_amp;
tbt_LFP_imag = tbt_LFP_imag./tbt_LFP_amp;

%%
cur_ch = 1;
xstart = 0.75;
xlen = 1;
trial_t_ax = (1:Nframes)*dt;
f1 = figure();
imagesc(trial_t_ax-xstart,1:Ntrials,squeeze(tbt_LFPs(:,:,cur_ch)));
caxis([-2.5 2.5]);
xlim([0 xlen]);
xlabel('Time (s)');
ylabel('Trial number');
colorbar


%%

[II,JJ] = meshgrid(1:length(test_trials));

clear all_LFP_metric all_rand_XC
for ww = 1:length(wfreqs)
    ww
    % ww =10;
    
    LFP_real = squeeze(tbt_LFP_real(:,:,ww,:));
    LFP_imag = squeeze(tbt_LFP_imag(:,:,ww,:));
    
    %compute bin edge locations for LFP metric
    all_LFP_dists = [];
    for tt = 1:Nframes
        cur_LFP_state = cat(2,squeeze(LFP_real(:,tt,:)),squeeze(LFP_imag(:,tt,:)));
        
        cur_Dmat = squareform(pdist(cur_LFP_state))/sqrt(length(use_lfps));
        cur_Dmat(JJ >= II) = nan;
        all_LFP_dists = cat(1,all_LFP_dists,cur_Dmat(~isnan(cur_Dmat)));
    end
    n_LFP_bins = 20;
    ED_bin_edges = prctile(all_LFP_dists,0:100/n_LFP_bins:100);
    ED_bin_centers = (ED_bin_edges(1:end-1)+ED_bin_edges(2:end))/2;
    
    
    LFP_metric_XC = zeros(length(ED_bin_centers),Nunits);
    LFP_metric_cnt = zeros(length(ED_bin_centers),Nunits);
    rand_XC = zeros(Nunits,1);
    rand_cnt = zeros(Nunits,1);
    
    LFP_metric_XC = nan(Nframes,n_LFP_bins,Nunits);
    rand_XC = nan(Nframes,Nunits);
    for tt = 1:Nframes
        cur_Robs = squeeze(tbt_Robs(:,tt,:));
        cur_LFP_state = cat(2,squeeze(LFP_real(:,tt,:)),squeeze(LFP_imag(:,tt,:)));
        
        cur_Dmat = squareform(pdist(cur_LFP_state))/sqrt(length(use_lfps));
        cur_Dmat(JJ == II) = nan;
        for jj = 1:length(ED_bin_centers)
            curset = find(cur_Dmat > ED_bin_edges(jj) & cur_Dmat <= ED_bin_edges(jj+1));
            %             temp = bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:));
            temp = cur_Robs(II(curset),:).*cur_Robs(JJ(curset),:);
%             LFP_metric_XC(jj,:) = LFP_metric_XC(jj,:,:) + nansum(temp,1);
            LFP_metric_XC(tt,jj,:) = nansum(temp,1);
%             LFP_metric_cnt(jj,:) = LFP_metric_cnt(jj,:,:) + sum(~isnan(temp),1);
            %         cur_XC(tt,jj,:,:) = (nanmean(bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:)),1));
        end
        curset = ~isnan(cur_Dmat);
        %         temp = bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:));
        temp = cur_Robs(II(curset),:).*cur_Robs(JJ(curset),:);
%         rand_XC = rand_XC + squeeze(nansum(temp,1))';
        rand_XC(tt,:) = squeeze(nansum(temp,1))';
%         rand_cnt = rand_cnt + squeeze(sum(~isnan(temp),1))';
        %     rand_XC(tt,:,:) = squeeze(nanmean(bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:)),1));
    end
    
%     LFP_metric_XC = LFP_metric_XC./LFP_metric_cnt;
%     rand_XC = rand_XC./rand_cnt;
LFP_metric_XC = squeeze(nanmean(LFP_metric_XC));
rand_XC = squeeze(nanmean(rand_XC))';

    tot_var = nanvar(reshape(tbt_Robs,[],Nunits));
    
    % normfac = sqrt(tot_var'*tot_var);
    normfac = tot_var;
    
    LFP_expl_cov = squeeze(LFP_metric_XC(1,:)) - rand_XC';
    LFP_expl_corr = LFP_expl_cov./normfac;
    
    all_LFP_metric(ww,:,:) = LFP_metric_XC;
    all_rand_XC(ww,:) = rand_XC;
end

%%
close all
for cc = 1:Nunits
    pcolor(1:n_LFP_bins,wfreqs,squeeze(all_LFP_metric(:,:,cc)));shading flat
    colorbar
    pause
    clf
end
%%

[II,JJ] = meshgrid(1:length(test_trials));

% clear all_LFP_metric all_rand_XC
% for ww = 1:length(wfreqs)
    ww
ww =  2;

LFP_real = squeeze(tbt_LFP_real(:,:,ww,:));
LFP_imag = squeeze(tbt_LFP_imag(:,:,ww,:));

%compute bin edge locations for LFP metric
all_LFP_dists = [];
for tt = 1:Nframes
    cur_LFP_state = cat(2,squeeze(LFP_real(:,tt,:)),squeeze(LFP_imag(:,tt,:)));
    
    cur_Dmat = squareform(pdist(cur_LFP_state))/sqrt(length(use_lfps));
    cur_Dmat(JJ >= II) = nan;
    all_LFP_dists = cat(1,all_LFP_dists,cur_Dmat(~isnan(cur_Dmat)));
end
n_LFP_bins = 20;
ED_bin_edges = prctile(all_LFP_dists,0:100/n_LFP_bins:100);
ED_bin_centers = (ED_bin_edges(1:end-1)+ED_bin_edges(2:end))/2;


LFP_metric_XC = zeros(length(ED_bin_centers),Nunits,Nunits);
LFP_metric_cnt = zeros(length(ED_bin_centers),Nunits,Nunits);
rand_XC = zeros(Nunits,Nunits);
rand_cnt = zeros(Nunits,Nunits);

for tt = 1:Nframes
    cur_Robs = squeeze(tbt_Robs(:,tt,:));
    cur_Robs2 = reshape(cur_Robs,[],1,Nunits);
    
    cur_LFP_state = cat(2,squeeze(LFP_real(:,tt,:)),squeeze(LFP_imag(:,tt,:)));
    
    cur_Dmat = squareform(pdist(cur_LFP_state))/sqrt(length(use_lfps));
    cur_Dmat(JJ == II) = nan;
    for jj = 1:length(ED_bin_centers)
        curset = find(cur_Dmat > ED_bin_edges(jj) & cur_Dmat <= ED_bin_edges(jj+1));
                    temp = bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:));
%         temp = cur_Robs(II(curset),:).*cur_Robs(JJ(curset),:);
        LFP_metric_XC(jj,:,:) = LFP_metric_XC(jj,:,:) + nansum(temp,1);
        LFP_metric_cnt(jj,:,:) = LFP_metric_cnt(jj,:,:) + sum(~isnan(temp),1);
    end
    curset = ~isnan(cur_Dmat);
            temp = bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:));
%     temp = cur_Robs(II(curset),:).*cur_Robs(JJ(curset),:);
%     rand_XC = rand_XC + squeeze(nansum(temp,1))';
%     rand_cnt = rand_cnt + squeeze(sum(~isnan(temp),1))';
    rand_XC = rand_XC + squeeze(nansum(temp,1));
    rand_cnt = rand_cnt + squeeze(sum(~isnan(temp),1));
    %     rand_XC(tt,:,:) = squeeze(nanmean(bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:)),1));
end

LFP_metric_XC = LFP_metric_XC./LFP_metric_cnt;
rand_XC = rand_XC./rand_cnt;

tot_var = nanvar(reshape(tbt_Robs,[],Nunits));

normfac = sqrt(tot_var'*tot_var);

% LFP_expl_cov = squeeze(LFP_metric_XC(1,:)) - rand_XC';
% LFP_expl_corr = LFP_expl_cov./normfac;
% LFP_expl_cov = squeeze(LFP_metric_XC(1,:,:)) - rand_XC;
LFP_expl_cov = squeeze(LFP_metric_XC(1,:,:));
LFP_expl_corr = LFP_expl_cov./normfac;

LFP_noise_cov = squeeze(LFP_metric_XC(1,:,:)) - rand_XC;
LFP_noise_corr = LFP_noise_cov./normfac;

stim_corr = rand_XC./normfac;

temp = reshape(tbt_Robs,[],Nunits);
noise_cov = squeeze(nanmean(bsxfun(@times,temp,reshape(temp,[],1,Nunits))));
corrected_noise_cov = noise_cov - LFP_expl_cov;
scorrected_noise_cov = noise_cov - rand_XC;
noise_cov = noise_cov./normfac;
corrected_noise_cov = corrected_noise_cov./normfac;
scorrected_noise_cov = scorrected_noise_cov./normfac;
% all_LFP_metric(ww,:,:) = LFP_metric_XC;
% all_rand_XC(ww,:) = rand_XC;

% end

%% compute xcorrs as a function of LFP state
[II,JJ] = meshgrid(1:length(test_trials));
cur_units = [13 30];
NcurUnits = length(cur_units);

max_tlag = 20; %max time lag for computing autocorrs
tlags = [-max_tlag:max_tlag];
Robs_shifted = nan(Ntrials,Nframes,NcurUnits,length(tlags));
for tt = 1:length(tlags)
    Robs_shifted(:,:,:,tt) = shift_matrix_Nd(tbt_Robs(:,:,cur_units),-tlags(tt),2);
end
% temp = reshape(Robs_shifted,[],NcurUnits,length(tlags));
% Robs_shifted = bsxfun(@minus,Robs_shifted,reshape(nanmean(temp),[1 1 NcurUnits length(tlags)]));

clear all_xc all_rxc
for ww = 1:length(wfreqs)
    ww
% ww =9;

LFP_real = squeeze(tbt_LFP_real(:,:,ww,:));
LFP_imag = squeeze(tbt_LFP_imag(:,:,ww,:));

%compute bin edge locations for LFP metric
all_LFP_dists = [];
for tt = 1:Nframes
    cur_LFP_state = cat(2,squeeze(LFP_real(:,tt,:)),squeeze(LFP_imag(:,tt,:)));
    
    cur_Dmat = squareform(pdist(cur_LFP_state))/sqrt(length(use_lfps));
    cur_Dmat(JJ >= II) = nan;
    all_LFP_dists = cat(1,all_LFP_dists,cur_Dmat(~isnan(cur_Dmat)));
end
n_LFP_bins = 15;
ED_bin_edges = prctile(all_LFP_dists,0:100/n_LFP_bins:100);
ED_bin_centers = (ED_bin_edges(1:end-1)+ED_bin_edges(2:end))/2;

LFP_metric_XC = zeros(length(ED_bin_centers),NcurUnits,NcurUnits,length(tlags));
LFP_metric_cnt = zeros(length(ED_bin_centers),NcurUnits,NcurUnits,length(tlags));
rand_XC = zeros(NcurUnits,NcurUnits,length(tlags));
rand_cnt = zeros(NcurUnits,NcurUnits,length(tlags));
same_XC = zeros(NcurUnits,NcurUnits,length(tlags));
same_cnt = zeros(NcurUnits,NcurUnits,length(tlags));

for tt = 1:Nframes
    cur_Robs = squeeze(tbt_Robs(:,tt,cur_units));
%     cur_Robs2 = reshape(cur_Robs,[],1,NcurUnits);
    cur_Robs2 = reshape(Robs_shifted(:,tt,:,:),[],1,NcurUnits,length(tlags));
    
    cur_LFP_state = cat(2,squeeze(LFP_real(:,tt,:)),squeeze(LFP_imag(:,tt,:)));
    
    cur_Dmat = squareform(pdist(cur_LFP_state))/sqrt(length(use_lfps));
    cur_Dmat(JJ == II) = nan;
    for jj = 1:length(ED_bin_centers)
        curset = find(cur_Dmat > ED_bin_edges(jj) & cur_Dmat <= ED_bin_edges(jj+1));
                    temp = bsxfun(@times,cur_Robs2(II(curset),:,:,:),cur_Robs(JJ(curset),:));
%         temp = cur_Robs(II(curset),:).*cur_Robs(JJ(curset),:);
        LFP_metric_XC(jj,:,:,:) = LFP_metric_XC(jj,:,:,:) + nansum(temp,1);
        LFP_metric_cnt(jj,:,:,:) = LFP_metric_cnt(jj,:,:,:) + sum(~isnan(temp),1);
    end
    curset = ~isnan(cur_Dmat);
            temp = bsxfun(@times,cur_Robs2(II(curset),:,:,:),cur_Robs(JJ(curset),:));
%     temp = cur_Robs(II(curset),:).*cur_Robs(JJ(curset),:);
%     rand_XC = rand_XC + squeeze(nansum(temp,1))';
%     rand_cnt = rand_cnt + squeeze(sum(~isnan(temp),1))';
    rand_XC = rand_XC + squeeze(nansum(temp,1));
    rand_cnt = rand_cnt + squeeze(sum(~isnan(temp),1));
    %     rand_XC(tt,:,:) = squeeze(nanmean(bsxfun(@times,cur_Robs2(II(curset),:,:),cur_Robs(JJ(curset),:)),1));

    curset = isnan(cur_Dmat);
    temp = bsxfun(@times,cur_Robs2(II(curset),:,:,:),cur_Robs(JJ(curset),:));
    same_XC = same_XC + squeeze(nansum(temp,1));
    same_cnt = same_cnt + squeeze(sum(~isnan(temp),1));
end

LFP_metric_XC = LFP_metric_XC./LFP_metric_cnt;
rand_XC = rand_XC./rand_cnt;
same_XC = same_XC./same_cnt;

% tot_var = nanvar(reshape(tbt_Robs(:,:,cur_units),[],NcurUnits));
% 
% normfac = sqrt(tot_var'*tot_var);
% 
% % LFP_expl_cov = squeeze(LFP_metric_XC(1,:)) - rand_XC';
% % LFP_expl_corr = LFP_expl_cov./normfac;
% LFP_expl_cov = squeeze(LFP_metric_XC(1,:,:)) - rand_XC;
% LFP_expl_corr = LFP_expl_cov./normfac;
% 
% stim_corr = rand_XC./normfac;

all_xc(ww,:,:) = squeeze(LFP_metric_XC(:,1,2,:));
all_rxc(ww,:) = squeeze(rand_XC(1,2,:));
all_sxc(ww,:) = squeeze(same_XC(1,2,:));
end

temp1 = reshape(Robs_shifted(:,:,2,:),[],1,length(tlags));
temp2 = reshape(tbt_Robs(:,:,cur_units(1)),[],1);
noise_xcor = squeeze(nanmean(bsxfun(@times,temp1,temp2)));
%%
%cc = 30 is interesting
% for cc = tr_set(all_mod_SU(tr_set) > 0)
    %%
    cc = 25
    
    cur_GQM = ModData(cc).rectGQM;
    cur_Robs = Robs_mat(:,cc);
    cur_Robs_us = Robs_mat_us(:,cc);
    
    cc_uinds = find(~isnan(cur_Robs));
    cc_tr_inds = cc_uinds(ismember(cc_uinds,tr_inds));
    cc_xv_inds = cc_uinds(ismember(cc_uinds,xv_inds));
    cc_all_inds = union(cc_tr_inds,cc_xv_inds);
    
    if cc > n_probes
        su_ind = find(SU_numbers == all_mod_SUnum(cc));
        [~,nearest_lfp] = min(abs(use_lfps-SU_probes(su_ind)));
    else
        [~,nearest_lfp] = min(abs(use_lfps - cc));
    end
    cc_uinds_us = find(~isnan(cur_Robs_us) & ~isnan(interp_lfps(:,nearest_lfp)));
    cc_tr_inds_us = cc_uinds_us(ismember(cc_uinds_us,tr_inds_us));
    cc_xv_inds_us = cc_uinds_us(ismember(cc_uinds_us,xv_inds_us));
    cc_all_inds_us = union(cc_tr_inds_us,cc_xv_inds_us);
    
    tr_cur_Robs_us = cur_Robs_us(cc_all_inds_us);
    tr_cur_Robs = cur_Robs(cc_all_inds);

    cur_uchs = 1:length(use_lfps);
    % cur_uchs(nearest_lfp) = [];
    
    %     cur_GQM = NMMfit_logexp_spkNL_withoffset(cur_GQM,cur_Robs,all_Xmat_shift(cc_uinds,:));
    cur_GQM = NMMfit_logexp_spkNL(cur_GQM,tr_cur_Robs,all_Xmat_shift(cc_all_inds,:));
    [~, ~, ~, tot_G,ind_Gints,fgint] = NMMmodel_eval(cur_GQM, [], all_Xmat_shift);
    modSigns = [cur_GQM.mods(:).sign];
    stimG = sum(bsxfun(@times,fgint,modSigns),2);
    excG = sum(fgint(:,modSigns==1),2);
    inhG = sum(fgint(:,modSigns==-1),2);
    
    g_tot = interp1(all_t_axis(used_inds),stimG,all_t_axis_us(used_inds_us(cc_all_inds_us)));
    g_tot(isnan(g_tot)) = 0;
    E_g = interp1(all_t_axis(used_inds),excG,all_t_axis_us(used_inds_us(cc_all_inds_us)));
    E_g(isnan(E_g)) = 0;
    I_g = interp1(all_t_axis(used_inds),inhG,all_t_axis_us(used_inds_us(cc_all_inds_us)));
    I_g(isnan(I_g)) = 0;
    Xcos_tot = bsxfun(@times,Xcos(cc_all_inds_us,:,:),g_tot);
    Xsin_tot = bsxfun(@times,Xsin(cc_all_inds_us,:,:),g_tot);
    Xcos_E = bsxfun(@times,Xcos(cc_all_inds_us,:,:),E_g);
    Xsin_E = bsxfun(@times,Xsin(cc_all_inds_us,:,:),E_g);
    Xcos_I = bsxfun(@times,Xcos(cc_all_inds_us,:,:),I_g);
    Xsin_I = bsxfun(@times,Xsin(cc_all_inds_us,:,:),I_g);
    
    % sac_stim_params(2:7) = NMMcreate_stim_params(length(cur_uchs));
    sac_stim_params(2:7) = NMMcreate_stim_params(1);
    sac_reg_params = NMMcreate_reg_params('lambda_d2T',1000);
    NL_types = repmat({'lin'},1,7);
    mod_signs = [1 1 1 1 1 1 1];
    Xtargets = [1 2 3 4 5 6 7];
    silent = 1;
    
%     for ww = 1:nwfreqs
%         ww
%         tr_stim{1} = [g_tot];
%         sac_stim_params(1) = NMMcreate_stim_params(size(tr_stim{1},2));
%         for hh = 1:length(cur_uchs)
%             tr_stim{2} = squeeze(Xcos(cc_uinds_us,cur_uchs(hh),ww));
%             tr_stim{3} = squeeze(Xsin(cc_uinds_us,cur_uchs(hh),ww));
%             %     tr_stim{4} = squeeze(Xcos_tot(:,cur_uchs,ww));
%             %     tr_stim{5} = squeeze(Xsin_tot(:,cur_uchs,ww));
%             tr_stim{4} = squeeze(Xcos_E(:,cur_uchs(hh),ww));
%             tr_stim{5} = squeeze(Xsin_E(:,cur_uchs(hh),ww));
%             tr_stim{6} = squeeze(Xcos_I(:,cur_uchs(hh),ww));
%             tr_stim{7} = squeeze(Xsin_I(:,cur_uchs(hh),ww));
%             gain_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
%             gain_mod.spk_NL_params = cur_GQM.spk_NL_params;
%             gain_mod.mods(1).reg_params = NMMcreate_reg_params();
%             gain_mod = NMMfit_filters(gain_mod,cur_Robs_us,tr_stim,[],[],silent);
%             [gain_LL,penLL, pred_rate, G, gint] = NMMmodel_eval(gain_mod,cur_Robs_us,tr_stim);
%             
%             all_filts = [gain_mod.mods(2:end).filtK];
%             
%             freq_ovgain(cc,ww,cur_uchs(hh)) = gain_mod.mods(1).filtK;
%             
%             freq_prefphase(cc,ww,cur_uchs(hh)) = atan2(all_filts(:,2),all_filts(:,1));
%             freq_prefamp(cc,ww,cur_uchs(hh)) = sqrt(sum(all_filts(:,1:2).^2,2));
%             
%             freq_gainphase(cc,ww,cur_uchs(hh)) = atan2(all_filts(:,4),all_filts(:,3));
%             freq_gainamp(cc,ww,cur_uchs(hh)) = sqrt(sum(all_filts(:,3:4).^2,2));
%             freq_Igainphase(cc,ww,cur_uchs(hh)) = atan2(all_filts(:,6),all_filts(:,5));
%             freq_Igainamp(cc,ww,cur_uchs(hh)) = sqrt(sum(all_filts(:,5:6).^2,2));
%             freq_LL(cc,ww) = gain_LL;
%             %     spkinds = convert_to_spikebins(cur_Robs_us);
%             %     freq_kappa(cc,ww) = circ_kappa(interp_phasegram(used_inds_us(cc_uinds_us(spkinds)),ww));
%             
%             ov_arate = mean(pred_rate);
%             freq_info(cc,ww) = mean(pred_rate.*log2(pred_rate/ov_arate))/ov_arate;
%         end
%     end
% end

%%
sac_stim_params(1:2) = NMMcreate_stim_params(1);
sac_stim_params(3:8) = NMMcreate_stim_params([length(cur_uchs) length(wfreqs)]);
sac_reg_params = NMMcreate_reg_params('lambda_d2T',1000,'lambda_L2',50);
% NL_types = repmat({'lin'},1,8);
% mod_signs = [1 1 1 1 1 1 1 1];
% Xtargets = [1 2 3 4 5 6 7 8];
NL_types = repmat({'lin'},1,4);
mod_signs = [1 1 1 1];
Xtargets = [1 2 3 4];
silent = 0;

tr_stim{1} = [E_g];
tr_stim{2} = [I_g];

tr_stim{3} = reshape(Xcos(cc_all_inds_us,:,:),[],length(cur_uchs)*length(wfreqs));
tr_stim{4} = reshape(Xsin(cc_all_inds_us,:,:),[],length(cur_uchs)*length(wfreqs));
tr_stim{5} = reshape(Xcos_E,[],length(cur_uchs)*length(wfreqs));
tr_stim{6} = reshape(Xsin_E,[],length(cur_uchs)*length(wfreqs));
tr_stim{7} = reshape(Xcos_I,[],length(cur_uchs)*length(wfreqs));
tr_stim{8} = reshape(Xsin_I,[],length(cur_uchs)*length(wfreqs));

gain_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
% gain_mod.spk_NL_params = cur_GQM.spk_NL_params;
gain_mod.mods(1).reg_params = NMMcreate_reg_params();
gain_mod.mods(2).reg_params = NMMcreate_reg_params();
% gain_mod.mods(1).filtK = 1;
% gain_mod.mods(2).filtK = 1;
gain_mod = NMMfit_filters(gain_mod,tr_cur_Robs_us,tr_stim,[],find(ismember(cc_all_inds_us,cc_tr_inds_us)),silent);
gain_mod = NMMfit_logexp_spkNL(gain_mod,Robs_mat_us(:,cc),tr_stim,[],find(ismember(cc_all_inds_us,cc_tr_inds_us)));
% [gain_LL,penLL, pred_rate, G, gint] = NMMeval_model(gain_mod,tr_cur_Robs_us,tr_stim);

% tempE = (gint(:,5)+gint(:,6))./E_g + gain_mod.mods(1).filtK;
% tempI = (gint(:,7)+gint(:,8))./I_g + gain_mod.mods(2).filtK;
% tempO = gint(:,3) + gint(:,4);

[gain_xvLL,null_xvLL, pred_xvrate, G, gint] = NMMeval_model(gain_mod,tr_cur_Robs_us,tr_stim,[],find(ismember(cc_all_inds_us,cc_xv_inds_us)));

[~,~, stim_pred_xvrate] = NMMeval_model(cur_GQM,tr_cur_Robs,all_Xmat_shift(cc_all_inds,:),[],find(ismember(cc_all_inds,cc_xv_inds)));
stim_pred_xvrate = interp1(all_t_axis(used_inds(cc_xv_inds)),stim_pred_xvrate,all_t_axis_us(used_inds_us(cc_xv_inds_us)))/(dt/new_dt);
stim_xvLL = nansum(tr_cur_Robs_us(ismember(cc_all_inds_us,cc_xv_inds_us)).*log(stim_pred_xvrate) - stim_pred_xvrate)/nansum(tr_cur_Robs_us(ismember(cc_all_inds_us,cc_xv_inds_us)));
%%
[~, ~, ~, tot_G,ind_Gints,fgint] = NMMmodel_eval(cur_GQM, [], all_Xmat_shift);
modSigns = [cur_GQM.mods(:).sign];
stimG = sum(bsxfun(@times,fgint,modSigns),2);
excG = sum(fgint(:,modSigns==1),2);
inhG = sum(fgint(:,modSigns==-1),2);

g_tot = interp1(all_t_axis(used_inds),stimG,all_t_axis_us(used_inds_us));
g_tot(isnan(g_tot)) = 0;
E_g = interp1(all_t_axis(used_inds),excG,all_t_axis_us(used_inds_us));
E_g(isnan(E_g)) = 0;
I_g = interp1(all_t_axis(used_inds),inhG,all_t_axis_us(used_inds_us));
I_g(isnan(I_g)) = 0;

Xcos_E = bsxfun(@times,Xcos,E_g);
Xsin_E = bsxfun(@times,Xsin,E_g);
Xcos_I = bsxfun(@times,Xcos,I_g);
Xsin_I = bsxfun(@times,Xsin,I_g);

tr_stim{1} = [E_g];
tr_stim{2} = [I_g];

tr_stim{3} = reshape(Xcos,[],length(cur_uchs)*length(wfreqs));
tr_stim{4} = reshape(Xsin,[],length(cur_uchs)*length(wfreqs));
tr_stim{5} = reshape(Xcos_E,[],length(cur_uchs)*length(wfreqs));
tr_stim{6} = reshape(Xsin_E,[],length(cur_uchs)*length(wfreqs));
tr_stim{7} = reshape(Xcos_I,[],length(cur_uchs)*length(wfreqs));
tr_stim{8} = reshape(Xsin_I,[],length(cur_uchs)*length(wfreqs));

rpt_trials = find(all_trial_Se == rpt_seed);
n_rpts = length(rpt_trials);
rpt_inds = find(ismember(all_trialvec(used_inds),rpt_trials));
rpt_inds_us = find(ismember(all_trialvec_us(used_inds_us),rpt_trials));

uu = rpt_inds_us(~isnan(Robs_mat_us(rpt_inds_us,cc)));
gain_mod_rpt = NMMfit_logexp_spkNL(gain_mod,Robs_mat_us(:,cc),tr_stim,[],uu);
[gain_LL,penLL,pred_rate,G,gint] = NMMmodel_eval(gain_mod_rpt,[],tr_stim);

uu = rpt_inds(~isnan(Robs_mat(rpt_inds,cc)));
cur_GQM_rpt = NMMfit_logexp_spkNL(cur_GQM,Robs_mat(:,cc),all_Xmat_shift,[],uu);
[gain_LL,penLL,stim_pred_rate] = NMMmodel_eval(cur_GQM_rpt,[],all_Xmat_shift);
stim_pred_rate = interp1(all_t_axis(used_inds),stim_pred_rate,all_t_axis_us(used_inds_us))/(dt/new_dt);

rpt_lfp_pred_rate = nan(n_rpts,750);
rpt_base_pred_rate = nan(n_rpts,750);
rpt_obs_rate = nan(n_rpts,750);
for ii = 1:n_rpts
    cur_tid = find(all_trialvec_us(used_inds_us) == rpt_trials(ii));
    rpt_lfp_pred_rate(ii,:) = pred_rate(cur_tid);
    rpt_base_pred_rate(ii,:) = stim_pred_rate(cur_tid);
    
    rpt_obs_rate(ii,:) = Robs_mat_us(cur_tid,cc);
end

lfp_resid = rpt_obs_rate - rpt_lfp_pred_rate;
tot_tbt_var = nanmean(nanvar(rpt_obs_rate));
lfp_tbt_var = nanmean(nanvar(lfp_resid));
stim_resid = rpt_obs_rate - rpt_base_pred_rate;
stim_tbt_var = nanmean(nanvar(stim_resid));
%%
%
%     Xcos_tot = bsxfun(@times,Xcos,g_tot);
%     Xsin_tot = bsxfun(@times,Xsin,g_tot);
%      Xcos_E = bsxfun(@times,Xcos,E_g);
%     Xsin_E = bsxfun(@times,Xsin,E_g);
%     Xcos_I = bsxfun(@times,Xcos,I_g);
%     Xsin_I = bsxfun(@times,Xsin,I_g);
%     clear tr_stim
%     tr_stim{2} = Xcos;
%     tr_stim{3} = Xsin;
%     tr_stim{4} = Xcos_E;
%     tr_stim{5} = Xsin_E;
%     tr_stim{6} = Xcos_I;
%     tr_stim{7} = Xsin_I;
%
%     sac_stim_params(1) = NMMcreate_stim_params(size(tr_stim{1},2));
%     sac_stim_params(2:7) = NMMcreate_stim_params(length(scales));
% %     sac_reg_params = NMMcreate_reg_params('lambda_d2T',10,'boundary_conds',[0 0 0]);
%     sac_reg_params = NMMcreate_reg_params('lambda_L2',5,'boundary_conds',[0 0 0]);
%     NL_types = {'lin','lin','lin','lin','lin','lin','lin'};
%
%     tr_stim{1} = [E_g I_g];
%     sac_stim_params(1) = NMMcreate_stim_params(size(tr_stim{1},2));
%     mod_signs = [1 1 1 1 1 1 1];
%     Xtargets = [1 2 3 4 5 6 7];
%     full_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
%     full_mod = NMMfit_filters(full_mod,cur_Robs_us(cc_tr_inds_us),get_Xcell_tInds(tr_stim,cc_tr_inds_us),[],[],silent);
%     full_xvLL = NMMmodel_eval(full_mod,cur_Robs_us(cc_xv_inds_us),get_Xcell_tInds(tr_stim,cc_xv_inds_us));
%
%
%      tr_stim{1} = g_tot;
%     tr_stim{4} = Xcos_tot;
%     tr_stim{5} = Xsin_tot;
%     sac_stim_params(1) = NMMcreate_stim_params(size(tr_stim{1},2));
%    mod_signs = [1 1 1 1 1];
%     Xtargets = [1 2 3 4 5];
%     gain_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
%     gain_mod = NMMfit_filters(gain_mod,cur_Robs_us(cc_tr_inds_us),get_Xcell_tInds(tr_stim,cc_tr_inds_us),[],[],silent);
%     gain_xvLL = NMMmodel_eval(gain_mod,cur_Robs_us(cc_xv_inds_us),get_Xcell_tInds(tr_stim,cc_xv_inds_us));
%      gain_LL = NMMmodel_eval(gain_mod,cur_Robs_us(cc_tr_inds_us),get_Xcell_tInds(tr_stim,cc_tr_inds_us));
%
%     mod_signs = [1 1 1];
%     Xtargets = [1 2 3];
%     off_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
%     off_mod = NMMfit_filters(off_mod,cur_Robs_us(cc_tr_inds_us),get_Xcell_tInds(tr_stim,cc_tr_inds_us),[],[],silent);
%     off_xvLL = NMMmodel_eval(off_mod,cur_Robs_us(cc_xv_inds_us),get_Xcell_tInds(tr_stim,cc_xv_inds_us));
%     off_LL = NMMmodel_eval(off_mod,cur_Robs_us(cc_tr_inds_us),get_Xcell_tInds(tr_stim,cc_tr_inds_us));
%     [~,~,off_predrate] = NMMmodel_eval(off_mod,cur_Robs_us,tr_stim);
%
% %     mod_signs = [1 1 1];
% %     Xtargets = [1 4 5];
% %     go_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
% %     go_mod = NMMfit_filters(go_mod,cur_Robs_us(cc_tr_inds_us),get_Xcell_tInds(tr_stim,cc_tr_inds_us),[],[],silent);
% %     go_xvLL = NMMmodel_eval(go_mod,cur_Robs_us(cc_xv_inds_us),get_Xcell_tInds(tr_stim,cc_xv_inds_us));
% %     go_LL = NMMmodel_eval(go_mod,cur_Robs_us(cc_tr_inds_us),get_Xcell_tInds(tr_stim,cc_tr_inds_us));
%
%     mod_signs = [1];
%     Xtargets = [1];
%     base_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
%     base_mod = NMMfit_filters(base_mod,cur_Robs_us(cc_tr_inds_us),get_Xcell_tInds(tr_stim,cc_tr_inds_us),[],[],silent);
%     base_xvLL = NMMmodel_eval(base_mod,cur_Robs_us(cc_xv_inds_us),get_Xcell_tInds(tr_stim,cc_xv_inds_us));
%     base_LL = NMMmodel_eval(base_mod,cur_Robs_us(cc_tr_inds_us),get_Xcell_tInds(tr_stim,cc_tr_inds_us));
%     [~,~,base_predrate] = NMMmodel_eval(base_mod,cur_Robs_us,tr_stim);
%
% %     mod_signs = [1 1];
% %     Xtargets = [2 3];
% %     lfp_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
% %     lfp_mod = NMMfit_filters(lfp_mod,cur_Robs(cc_tr_inds),get_Xcell_tInds(tr_stim,cc_tr_inds),[],[],silent);
% %     lfp_xvLL = NMMmodel_eval(lfp_mod,cur_Robs(cc_xv_inds),get_Xcell_tInds(tr_stim,cc_xv_inds));
% %     lfp_LL = NMMmodel_eval(lfp_mod,cur_Robs(cc_tr_inds),get_Xcell_tInds(tr_stim,cc_tr_inds));
% %     [~,~,lfp_predrate] = NMMmodel_eval(lfp_mod,cur_Robs,tr_stim);
%
%     null_pred = ones(size(cc_xv_inds_us))*mean(cur_Robs_us(cc_tr_inds_us));
%     null_xvLL = sum(cur_Robs_us(cc_xv_inds_us).*log(null_pred)-null_pred)/sum(cur_Robs_us(cc_xv_inds_us));
%     null_pred = ones(size(cc_tr_inds_us))*mean(cur_Robs_us(cc_tr_inds_us));
%     null_LL = sum(cur_Robs_us(cc_tr_inds_us).*log(null_pred)-null_pred)/sum(cur_Robs_us(cc_tr_inds_us));
%
%     all_filts = [off_mod.mods(2:end).filtK];
%     offset_amp = sqrt(all_filts(:,1).^2 + all_filts(:,2).^2);
%     all_filts = [gain_mod.mods(2:end).filtK];
%     gain_amp = sqrt(all_filts(:,3).^2 + all_filts(:,4).^2);
%     % Egain_amp = sqrt(all_filts(:,3).^2 + all_filts(:,4).^2);
%     % Igain_amp = sqrt(all_filts(:,5).^2 + all_filts(:,6).^2);
%
%     off_imp(cc) = off_xvLL - base_xvLL;
%     base_imp(cc) = base_xvLL - null_xvLL;
% %     lfp_imp(cc) = lfp_xvLL - null_xvLL;
%     gain_imp(cc) = gain_xvLL - base_xvLL;
%     full_imp(cc) = full_xvLL - base_xvLL;
%     off_rimp(cc) = off_LL - base_LL;
%     base_rimp(cc) = base_LL - null_LL;
%     gain_rimp(cc) = gain_LL - base_LL;
%
%     cur_off_xvLL(cc) = off_xvLL;
%     cur_base_xvLL(cc) = base_xvLL;
%     cur_null_xvLL(cc) = null_xvLL;
% %     cur_lfp_xvLL(cc) = lfp_xvLL;
%
%
%     X = [tr_stim{2} tr_stim{3}];
%     Xcov = cov(X);
%     full_filt = [off_mod.mods(2).filtK; off_mod.mods(3).filtK];
%     re_filt = Xcov*full_filt;
%     cond_dep_offset(cc,:) = sqrt(re_filt(1:nwfreqs).^2 + re_filt((nwfreqs+1):end).^2);
%     full_gfilt = [gain_mod.mods(4).filtK; gain_mod.mods(5).filtK];
%     re_gfilt = Xcov*full_gfilt;
%     cond_dep_gain(cc,:) = sqrt(re_gfilt(1:nwfreqs).^2 + re_gfilt((nwfreqs+1):end).^2);
%     full_filt = [gain_mod.mods(2).filtK; gain_mod.mods(3).filtK];
%     re_filt = Xcov*full_filt;
%     cond_dep_off_gmod(cc,:) = sqrt(re_filt(1:nwfreqs).^2 + re_filt((nwfreqs+1):end).^2);
%
%     full_filt = [full_mod.mods(2).filtK; full_mod.mods(3).filtK];
%     re_filt = Xcov*full_filt;
%     cond_fmod_off(cc,:) = sqrt(re_filt(1:nwfreqs).^2 + re_filt((nwfreqs+1):end).^2);
%     full_filt = [full_mod.mods(4).filtK; full_mod.mods(5).filtK];
%     re_filt = Xcov*full_filt;
%     cond_fmod_eg(cc,:) = sqrt(re_filt(1:nwfreqs).^2 + re_filt((nwfreqs+1):end).^2);
%     full_filt = [full_mod.mods(6).filtK; full_mod.mods(7).filtK];
%     re_filt = Xcov*full_filt;
%     cond_fmod_ig(cc,:) = sqrt(re_filt(1:nwfreqs).^2 + re_filt((nwfreqs+1):end).^2);

%%
%     spkinds = convert_to_spikebins(cur_Robs);
%     phase_lock(cc,:) = nan(1,length(wfreqs));
%     for ww = 1:length(wfreqs)
%         phase_lock(cc,ww) = circ_kappa(interp_phasegram(used_inds(spkinds),ww));
%     end


% end

%%
% phase_bin_edges = linspace(0,2*pi,30);
% phase_bin_cents = 0.5*phase_bin_edges(1:end-1) + 0.5*phase_bin_edges(2:end);
%
% phase_avg = nan(length(wfreqs),length(phase_bin_cents));
% phase_pred_avg = nan(length(wfreqs),length(phase_bin_cents));
% for ii = 1:length(phase_bin_cents)
%     for ww = 1:length(wfreqs)
%         curset = find(interp_phasegram(used_inds(cc_uinds),ww) >= phase_bin_edges(ii) & ...
%             interp_phasegram(used_inds(cc_uinds),ww) < phase_bin_edges(ii+1));
%         phase_avg(ww,ii) = mean(cur_Robs(curset));
%         phase_pred_avg(ww,ii) = mean(off_predrate(curset));
%     end
% end
