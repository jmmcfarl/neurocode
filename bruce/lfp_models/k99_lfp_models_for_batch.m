%
% clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements//');
addpath('~/James_scripts/bruce/processing/');
addpath('~/James_scripts/bruce/saccade_modulation/');
addpath('~/James_scripts/TentBasis2D/');

global Expt_name bar_ori use_MUA


% Expt_name = 'M296';
% Expt_name = 'G093';
use_MUA = false;
% bar_ori = 0; %bar orientation to use (only for UA recs)


fit_unCor = false;

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

if strcmp(Expt_name,'M297') && bar_ori == 0
   mod_data_name = strcat(mod_data_name,'_old2'); %changed the  
end
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
expt_sac_amp = nan(1,length(Expts));
expt_imback = nan(1,length(Expts));
included_type = false(1,length(Expts));
for ii = 1:length(Expts)
    if ~isempty(Expts{ii})
        expt_names{ii} = Expts{ii}.Header.expname;
        expt_dds(ii) = Expts{ii}.Stimvals.dd;
        expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
        expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
        expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
        expt_sac_amp(ii) = Expts{ii}.Stimvals.Fs;
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
if length(unique(expt_dds(cur_block_set))) > 1
    fprintf('Warning, multiple dds detected!\n');
    main_dds = mode(expt_dds(cur_block_set));
    cur_block_set(expt_dds(cur_block_set) ~= main_dds) = [];
end


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
%for G093 use only data where stripe width is AT LEAST 2 deg
if strcmp(Expt_name,'G093')
    un_wi_vals = unique(all_trial_wi);
    use_wi_trials = find(all_trial_wi >= un_wi_vals(2));
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
if isfield(Expts{cur_block_set(1)}.Header,'exptno')
    em_block_nums = cellfun(@(X) X.Header.exptno,Expts(cur_block_set),'uniformoutput',1); %block numbering for EM/LFP data sometimes isnt aligned with Expts struct
else
    em_block_nums = cur_block_set;
end

[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data_v2(all_t_axis,all_blockvec,em_block_nums,Expt_name,trial_toffset,good_coils);
% [all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);
% [all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data_v2(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset,good_coils);
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
if ~fit_unCor
    for i=1:NT
        best_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
    end
end
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
    use_lfps = 1:2:n_probes; %use every other probe for estimating full LFP filters
elseif strcmp(rec_type,'UA')
    Fs = 400.0032;
    dsf = 2;
    Fsd = Fs/dsf;
    niqf = Fs/2;
    new_niqf = Fsd/2;
    [bb,aa] = butter(4,0.8*new_niqf/niqf);
    use_lfps = SU_probes;
    %     use_lfps = [28 SU_probes];
end

%wavelet parameters
nwfreqs = 15;
min_freq = 1.5; max_freq = 40;
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
        fname = sprintf('lemM%dA.%d.lfp.mat',Expt_num,em_block_nums(ee));
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
        lfp_fname = sprintf('Expt%d_LFP.mat',em_block_nums(ee));
        load(lfp_fname);
        
        cur_lfps = bsxfun(@times,double(lfp_mat(:,use_lfps)),lfp_int2V(use_lfps)');
        if dsf > 1
            cur_lfps = filtfilt(bb,aa,cur_lfps);
            expt_lfps = downsample(cur_lfps,dsf);
            expt_lfp_t_axis = downsample(lfp_t_ax',dsf);
        end
           expt_cwt = nan(size(expt_lfps,1),length(wfreqs),length(use_lfps));
        for cc = 1:length(use_lfps)
            expt_cwt(:,:,cc) = cwt(expt_lfps(:,cc),scales,'cmor1-1')';
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

%% Interpolate raw LFPs onto t-axis
interp_lfps = interp1(full_lfp_taxis,full_lfps,all_t_axis);
interp_lfps = interp_lfps(used_inds,:);
interp_lfps = nanzscore(interp_lfps);

%% Interpolate CWT components
interp_cwt_real = interp1(full_lfp_taxis,real(full_cwts),all_t_axis);
interp_cwt_imag = interp1(full_lfp_taxis,imag(full_cwts),all_t_axis);
interp_cwt_mag = sqrt(interp_cwt_real.^2 + interp_cwt_imag.^2);

interp_cwt_real = bsxfun(@rdivide,interp_cwt_real,nanstd(interp_cwt_mag));
interp_cwt_imag = bsxfun(@rdivide,interp_cwt_imag,nanstd(interp_cwt_mag));

%% compute X-mat real and imaginary components for model-fitting
Xreal = interp_cwt_real(used_inds,:,:);
Ximag = interp_cwt_imag(used_inds,:,:);
clear interp_cwt*
Xmag = sqrt(Xreal.^2 + Ximag.^2);

%prevent some occasional nans from messing up models
Xreal(isnan(Xreal)) = 0;
Ximag(isnan(Ximag)) = 0;
Xmag(isnan(Xmag)) = 0;

%% FIT LFP PHASE-DEPENDENT MODELS
uwfreqs = find(wfreqs < 50); %use all freqs up to this values
addpath(genpath('~/James_scripts/chronux/spectral_analysis/'))
% for cc = 26
for cc = (n_probes+1):length(ModData)
    %%
    cc    
    cur_GQM = ModData(cc).rectGQM;
    cur_Robs = Robs_mat(:,cc);
    tot_nspks = nansum(cur_Robs);
    avg_rate = nanmean(cur_Robs)/dt;
    
    cc_uinds = find(~isnan(cur_Robs));
    cc_all_inds = cc_uinds(ismember(cc_uinds,all_inds));
    cur_tr_ind_set = find(ismember(cc_all_inds,tr_inds));
    cur_xv_ind_set = find(ismember(cc_all_inds,xv_inds));
    
    if ~isempty(cur_xv_ind_set)
    if cc > n_probes
        su_ind = find(SU_numbers == all_mod_SUnum(cc));
        [~,nearest_lfp] = min(abs(use_lfps-SU_probes(su_ind)));
    else
        [~,nearest_lfp] = min(abs(use_lfps - cc));
    end
        
    tr_cur_Robs = cur_Robs(cc_all_inds);
    
    %fit spk NL params of base stim-proc model
    cur_GQM = NMMfit_logexp_spkNL(cur_GQM,tr_cur_Robs,all_Xmat_shift(cc_all_inds,:));
    
    %compute outputs of Exc and Inh subunits
    [~, ~, ~, tot_G,ind_Gints,fgint] = NMMmodel_eval(cur_GQM, [], all_Xmat_shift);
    modSigns = [cur_GQM.mods(:).sign];
    g_tot = sum(bsxfun(@times,fgint,modSigns),2);
    E_g = sum(fgint(:,modSigns==1),2);
    I_g = sum(fgint(:,modSigns==-1),2);
    
    if strcmp(rec_type,'LP')
        uchs = 1:length(use_lfps);
    else
        uchs = nearest_lfp;
    end
    
    %compute product of LFP X_mat terms and Exc and Inh inputs.
    %     Xreal_tot = bsxfun(@times,Xreal(cc_all_inds,:,:),g_tot(cc_all_inds));
    %     Ximag_tot = bsxfun(@times,Ximag(cc_all_inds,:,:),g_tot(cc_all_inds));
    Xreal_E = bsxfun(@times,Xreal(cc_all_inds,uwfreqs,uchs),E_g(cc_all_inds));
    Ximag_E = bsxfun(@times,Ximag(cc_all_inds,uwfreqs,uchs),E_g(cc_all_inds));
    Xreal_I = bsxfun(@times,Xreal(cc_all_inds,uwfreqs,uchs),I_g(cc_all_inds));
    Ximag_I = bsxfun(@times,Ximag(cc_all_inds,uwfreqs,uchs),I_g(cc_all_inds));
    
    %just scalar E and I inputs
    tr_stim{1} = [E_g(cc_all_inds)];
    tr_stim{2} = [I_g(cc_all_inds)];
    sac_stim_params(1:2) = NMMcreate_stim_params(1);
    
    %X-mat 3 and 4 are the additive LFP terms
    tr_stim{3} = reshape(Xreal(cc_all_inds,uwfreqs,uchs),length(cc_all_inds),[]);
    tr_stim{4} = reshape(Ximag(cc_all_inds,uwfreqs,uchs),length(cc_all_inds),[]);
    
    %5-8 are the E and I gain terms
    tr_stim{5} = reshape(Xreal_E,length(cc_all_inds),[]);
    tr_stim{6} = reshape(Ximag_E,length(cc_all_inds),[]);
    tr_stim{7} = reshape(Xreal_I,length(cc_all_inds),[]);
    tr_stim{8} = reshape(Ximag_I,length(cc_all_inds),[]);
    
    poss_lambda = logspace(log10(10),log10(1e4),6);
%     poss_lambda = [10 100 1000 10000];
    lambda_xvLLs = nan(length(poss_lambda),1);
    lambda_models = cell(length(poss_lambda),1);
    for ll = 1:length(poss_lambda)
        fprintf('Poss lambda %d/%d\n',ll,length(poss_lambda));
        %these reg values have been ballparked but also roughly checked as
        %cross-validating well
        if strcmp(rec_type,'LP')
            sac_reg_params = NMMcreate_reg_params('lambda_d2XT',100,'lambda_L2',poss_lambda(ll)); %smoothing reg across freq and depth
        else
            sac_reg_params = NMMcreate_reg_params('lambda_d2T',100,'lambda_L2',poss_lambda(ll));
        end
        
        sac_stim_params(3:8) = NMMcreate_stim_params([length(uwfreqs),length(uchs)]);
        NL_types = repmat({'lin'},1,8);
        mod_signs = [1 1 1 1 1 1 1 1];
        Xtargets = [1 2 3 4 5 6 7 8];
        silent = 1;
        
        gain_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
        gain_mod.spk_NL_params = cur_GQM.spk_NL_params;
        gain_mod.mods(1).reg_params = NMMcreate_reg_params();
        gain_mod.mods(2).reg_params = NMMcreate_reg_params();
        gain_mod = NMMfit_filters(gain_mod,tr_cur_Robs,tr_stim,[],cur_tr_ind_set,silent);
        [gain_LL,nullLL, pred_rate, G, gint] = NMMeval_model(gain_mod,tr_cur_Robs,tr_stim,[],cur_tr_ind_set);
        
        %output of model on xv set
        [gain_xvLL,nullxvLL,~,~,gintxv] = NMMeval_model(gain_mod,tr_cur_Robs,tr_stim,[],cur_xv_ind_set);        
        lambda_LL(cc,ll) = gain_xvLL;
        lambda_xvLLs(ll) = gain_xvLL;
        lambda_models{ll} = gain_mod;
    end
  
    [gain_xvLL,best_reg_loc] = max(lambda_xvLLs);
    gain_mod = lambda_models{best_reg_loc};
    [gain_LL,nullLL, gain_rate, G, gint] = NMMeval_model(gain_mod,tr_cur_Robs,tr_stim,[],cur_tr_ind_set);
    
%     [gain_xvLL,nullxvLL,~,~,gintxv] = NMMeval_model(gain_mod,tr_cur_Robs,tr_stim,[],cur_xv_ind_set);
%     
%         %re-estimate the weights of each model component on the xval set
%     xvstim_params(1:8) = NMMcreate_stim_params(1);
%     clear cur_stim
%     for xx = 1:8
%         cur_stim{xx} = gintxv(:,xx);
%     end
%     xvgain_mod = NMMinitialize_model(xvstim_params,ones(1,8),NL_types,NMMcreate_reg_params(),1:8);
%     xvgain_mod = NMMfit_filters(xvgain_mod,tr_cur_Robs(cur_xv_ind_set),cur_stim,[],[],silent);
%     [~,~,~,~,gintxv2] = NMMeval_model(xvgain_mod,tr_cur_Robs(cur_xv_ind_set),cur_stim);

if strcmp(rec_type,'LP')
    sac_reg_params = NMMcreate_reg_params('lambda_d2XT',100,'lambda_L2',poss_lambda(best_reg_loc)); %smoothing reg across freq and depth
else
    sac_reg_params = NMMcreate_reg_params('lambda_d2T',100,'lambda_L2',poss_lambda(best_reg_loc));
end

    %for comparison fit a model that just has the additive LFP filter
    off_mod = NMMinitialize_model(sac_stim_params,[1 1 1 1],{'lin','lin','lin','lin'},sac_reg_params,[1 2 3 4]);
    off_mod.spk_NL_params = cur_GQM.spk_NL_params;
    off_mod.mods(1).reg_params = NMMcreate_reg_params();
    off_mod.mods(2).reg_params = NMMcreate_reg_params();
    off_mod = NMMfit_filters(off_mod,tr_cur_Robs,tr_stim,[],cur_tr_ind_set,silent);
    [off_xvLL] = NMMeval_model(off_mod,tr_cur_Robs,tr_stim,[],cur_xv_ind_set);

    %and a directly comparable 'base' stim-proc model
    base_mod = NMMinitialize_model(sac_stim_params(1:2),[1 1],{'lin','lin'},NMMcreate_reg_params(),[1 2]);
    base_mod = NMMfit_filters(base_mod,tr_cur_Robs,tr_stim,[],cur_tr_ind_set,silent);
    [base_xvLL] = NMMeval_model(base_mod,tr_cur_Robs,tr_stim,[],cur_xv_ind_set);
    [~,~,base_rate] = NMMeval_model(base_mod,tr_cur_Robs,tr_stim,[],cur_tr_ind_set);
    
    %these are all the LFP filters
    all_filts = [gain_mod.mods(3:end).filtK];
    
%     ov_arate = mean(pred_rate);
    
    pred_rate_diff = gain_rate - base_rate;
    
    %compute effective gains of the Exc and Inh inputs
    [gain_LL,nullLL, pred_rate, G, gint] = NMMeval_model(gain_mod,tr_cur_Robs,tr_stim);
    Eout = sum(gint(:,[1 5 6]),2);
    Iout = sum(gint(:,[2 7 8]),2);
    Oout = sum(gint(:,[3 4]),2);
    Egain = Eout./E_g(cc_all_inds);
    Igain = Iout./I_g(cc_all_inds);
    Egain = Egain - nanmean(Egain);
    Igain = Igain - nanmean(Igain);
%     Igain(isnan(Igain)) = 1;
%     Egain(isnan(Egain)) = 1;
    Offset = Oout - nanmean(Oout); 
    
%     %gains evaluated on xval set
%     Eoutxv = sum(gintxv(:,[1 5 6]),2);
%     Ioutxv = sum(gintxv(:,[2 7 8]),2);
%     Ooutxv = sum(gintxv(:,[3 4]),2);
%     Egainxv = Eoutxv./E_g(cc_all_inds(cur_xv_ind_set));
%     Igainxv = Ioutxv./I_g(cc_all_inds(cur_xv_ind_set));
% 
%     %gains with weights re-estimated on xval set
%     Eoutxv2 = sum(gintxv2(:,[1 5 6]),2);
%     Ioutxv2 = sum(gintxv2(:,[2 7 8]),2);
%     Ooutxv2 = sum(gintxv2(:,[3 4]),2);
%     Egainxv2 = Eoutxv2./E_g(cc_all_inds(cur_xv_ind_set));
%     Igainxv2 = Ioutxv2./I_g(cc_all_inds(cur_xv_ind_set));
    
    lfp_models(cc).Egain_SD = nanstd(Egain);
    lfp_models(cc).Igain_SD = nanstd(Igain);
    lfp_models(cc).gain_rate_diff_SD = nanstd(pred_rate_diff);
    lfp_models(cc).base_rate_SD = nanstd(base_rate);
    lfp_models(cc).base_avg_rate = nanmean(base_rate);
%     lfp_models(cc).Egain_SDxv = nanstd(Egainxv);
%     lfp_models(cc).Igain_SDxv = nanstd(Igainxv);
%     lfp_models(cc).Egain_SDxv2 = nanstd(Egainxv2);
%     lfp_models(cc).Igain_SDxv2 = nanstd(Igainxv2);
% 
    %% compute power spectra and coherence and xcorr of E and I gain signals
    
    Egain(isnan(Egain)) = 1;
    Igain(isnan(Igain)) = 1;
    Egain = Egain - nanmean(Egain);
    Igain = Igain - nanmean(Igain);
    Offset(isnan(Offset)) = nanmean(Offset);
    
    movingwin = [1.5 1.5];
    params.Fs = 1/dt;
    params.tapers = [3 5];
    cur_tvec = all_trialvec(used_inds(cc_all_inds));
    cur_tstart_inds = [1; find(cur_tvec(2:end) > cur_tvec(1:end-1))]; 
    cur_tend_inds = [cur_tstart_inds(2:end); length(cur_tvec)];
    sMarkers = [cur_tstart_inds cur_tend_inds];
    [Cmn,Phimn,Smn,Smm,f] = coherencyc_unequal_length_trials([Egain(:) Igain(:) Offset(:)], movingwin, params, sMarkers );
    
    [xc,lags] = xcov(-Igain(:),Egain(:),round(0.2/dt),'coeff');
    %%
    lfp_models(cc).pow_spectra = Smm;
    lfp_models(cc).spectra_f = f;
    lfp_models(cc).coherence = Cmn;
    lfp_models(cc).xcov = xc;
    lfp_models(cc).lags = lags*dt;
    lfp_models(cc).base_xvImp = (base_xvLL - nullxvLL)/log(2);
    lfp_models(cc).lfp_xvImp = (gain_xvLL - nullxvLL)/log(2);
    lfp_models(cc).off_xvImp = (off_xvLL - nullxvLL)/log(2);
    
    lfp_models(cc).ModData = ModData(cc);
    lfp_models(cc).best_reg = poss_lambda(best_reg_loc);
    lfp_models(cc).gain_mod = gain_mod;
    end
end

%% FOR FITTING MODELS TO LFP POWER
% clear tr_stim
% uwfreqs = find(wfreqs < 50);
% addpath(genpath('~/James_scripts/chronux/spectral_analysis/'))
% % for cc = 26
% for cc = (n_probes+1):length(ModData)
%     %%
%     cc
%     %     cc = 26
%     
%     cur_GQM = ModData(cc).rectGQM;
%     cur_Robs = Robs_mat(:,cc);
%     
%     cc_uinds = find(~isnan(cur_Robs));
% %     cc_tr_inds = cc_uinds(ismember(cc_uinds,tr_inds));
% %     cc_xv_inds = cc_uinds(ismember(cc_uinds,xv_inds));
% %     cc_all_inds = union(cc_tr_inds,cc_xv_inds);
%     cc_all_inds = cc_uinds(ismember(cc_uinds,all_inds));
%     cur_tr_ind_set = find(ismember(cc_all_inds,tr_inds));
%     cur_xv_ind_set = find(ismember(cc_all_inds,xv_inds));
%     
%     if ~isempty(cur_xv_ind_set)
%         
%     if cc > n_probes
%         su_ind = find(SU_numbers == all_mod_SUnum(cc));
%         [~,nearest_lfp] = min(abs(use_lfps-SU_probes(su_ind)));
%     else
%         [~,nearest_lfp] = min(abs(use_lfps - cc));
%     end
%     
%     if strcmp(rec_type,'LP')
%         uchs = 1:length(use_lfps);
%     sac_reg_params = NMMcreate_reg_params('lambda_d2XT',100,'lambda_L2',1000);
%     else
%         uchs = nearest_lfp;
%     sac_reg_params = NMMcreate_reg_params('lambda_d2T',100,'lambda_L2',1000);
%     end
%     
%     tr_cur_Robs = cur_Robs(cc_all_inds);
%     
%     cur_GQM = NMMfit_logexp_spkNL(cur_GQM,tr_cur_Robs,all_Xmat_shift(cc_all_inds,:));
%     [~, ~, ~, tot_G,ind_Gints,fgint] = NMMmodel_eval(cur_GQM, [], all_Xmat_shift);
%     modSigns = [cur_GQM.mods(:).sign];
%     g_tot = sum(bsxfun(@times,fgint,modSigns),2);
%     E_g = sum(fgint(:,modSigns==1),2);
%     I_g = sum(fgint(:,modSigns==-1),2);
%     
%     Xmag_E = bsxfun(@times,Xmag(cc_all_inds,uwfreqs,uchs),E_g(cc_all_inds));
%     Xmag_I = bsxfun(@times,Xmag(cc_all_inds,uwfreqs,uchs),I_g(cc_all_inds));
%     
%     sac_stim_params(3:5) = NMMcreate_stim_params([length(uwfreqs),length(uchs)]);
% %     sac_reg_params = NMMcreate_reg_params('lambda_d2XT',poss_lambda(ll),'lambda_L2',50);
% %     sac_reg_params = NMMcreate_reg_params('lambda_d2XT',100,'lambda_L2',poss_lambda(ll));
%     NL_types = repmat({'lin'},1,5);
%     mod_signs = [1 1 1 1 1];
%     Xtargets = [1 2 3 4 5];
%     silent = 1;
%     
%     tr_stim{1} = [E_g(cc_all_inds)];
%     tr_stim{2} = [I_g(cc_all_inds)];
%     sac_stim_params(1:2) = NMMcreate_stim_params(1);
%     
%     tr_stim{3} = reshape(Xmag(cc_all_inds,uwfreqs,uchs),length(cc_all_inds),[]);
%     tr_stim{4} = reshape(Xmag_E,length(cc_all_inds),[]);
%     tr_stim{5} = reshape(Xmag_I,length(cc_all_inds),[]);
%     
%     gain_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
%     gain_mod.spk_NL_params = cur_GQM.spk_NL_params;
%     gain_mod.mods(1).reg_params = NMMcreate_reg_params();
%     gain_mod.mods(2).reg_params = NMMcreate_reg_params();
%     gain_mod = NMMfit_filters(gain_mod,tr_cur_Robs,tr_stim,[],cur_tr_ind_set,silent);
%     [gain_LL,nullLL, pred_rate, G, gint] = NMMeval_model(gain_mod,tr_cur_Robs,tr_stim,[],cur_tr_ind_set);
%     [gain_xvLL,nullxvLL] = NMMeval_model(gain_mod,tr_cur_Robs,tr_stim,[],cur_xv_ind_set);
% %     lambda_LL(cc,ll) = gain_xvLL;
% %     end
%     
%     base_mod = NMMinitialize_model(sac_stim_params(1:2),[1 1],{'lin','lin'},NMMcreate_reg_params(),[1 2]);
%     base_mod = NMMfit_filters(base_mod,tr_cur_Robs,tr_stim,[],cur_tr_ind_set,silent);
%     [base_xvLL] = NMMeval_model(base_mod,tr_cur_Robs,tr_stim,[],cur_xv_ind_set);
%     
%     all_filts = [gain_mod.mods(3:end).filtK];
%     
%     ov_arate = mean(pred_rate);
%     
%     [gain_LL,nullLL, pred_rate, G, gint] = NMMeval_model(gain_mod,tr_cur_Robs,tr_stim);
%     Eout = sum(gint(:,[1 4]),2);
%     Iout = sum(gint(:,[2 5]),2);
%     Oout = sum(gint(:,[3]),2);
%     Egain = Eout./E_g(cc_all_inds);
%     Igain = Iout./I_g(cc_all_inds);
%     Egain = Egain - nanmean(Egain);
%     Igain = Igain - nanmean(Igain);
%     Igain(isnan(Igain)) = 1;
%     Egain(isnan(Egain)) = 1;
%     Offset = Oout - nanmean(Oout); 
%     %%
%     movingwin = [1.5 1.5];
%     params.Fs = 1/dt;
%     params.tapers = [3 5];
%     cur_tvec = all_trialvec(used_inds(cc_all_inds));
%     cur_tstart_inds = [1; find(cur_tvec(2:end) > cur_tvec(1:end-1))]; 
%     cur_tend_inds = [cur_tstart_inds(2:end); length(cur_tvec)];
%     sMarkers = [cur_tstart_inds cur_tend_inds];
%     [Cmn,Phimn,Smn,Smm,f] = coherencyc_unequal_length_trials([Egain(:) Igain(:) Offset(:)], movingwin, params, sMarkers );
%     
%     [xc,lags] = xcov(-Igain(:),Egain(:),round(0.2/dt),'coeff');
%     %%
%     lfp_ampmodels(cc).Egain_SD = nanstd(Egain);
%     lfp_ampmodels(cc).Igain_SD = nanstd(Igain);
%     lfp_ampmodels(cc).pow_spectra = Smm;
%     lfp_ampmodels(cc).spectra_f = f;
%     lfp_ampmodels(cc).coherence = Cmn;
%     lfp_ampmodels(cc).xcov = xc;
%     lfp_ampmodels(cc).lags = lags*dt;
%     lfp_ampmodels(cc).base_xvImp = (base_xvLL - nullxvLL)/log(2);
%     lfp_ampmodels(cc).lfp_xvImp = (gain_xvLL - nullxvLL)/log(2);
%     end
% end

%%
anal_dir = ['~/Analysis/bruce/' Expt_name '/lfp_models/'];
if ~exist(anal_dir)
    mkdir(anal_dir)
end
cd(anal_dir);

sname = 'lfp_models4';
sname = [sname sprintf('_ori%d',bar_ori)];

% save(sname,'targs','lfp_ampmodels','lfp_models');
save(sname,'targs','lfp_models');