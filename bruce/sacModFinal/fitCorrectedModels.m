clear all

addpath('~/James_scripts/bruce/processing/');
addpath('~/James_scripts/bruce/eye_tracking_improvements/');

global Expt_name bar_ori use_MUA fit_unCor

Expt_name = 'G085';
bar_ori = 90; %bar orientation to use (only for UA recs)
use_MUA = false;
fit_unCor = true; %fit models to uncorrected stim?

%%
%load in array RF position data
load ~/Data/bruce/general_array_data/array_pos_data.mat
interp_ecc = sqrt(interp_x.^2 + interp_y.^2);

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
save_dir = ['~/Analysis/bruce/' Expt_name '/models'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

% et_mod_data_name = 'full_eyetrack_initmods';
% et_anal_name = 'full_eyetrack';
et_mod_data_name = 'full_eyetrack_initmods_Rinit';
et_anal_name = 'full_eyetrack_Rinit';
save_name = 'corrected_models';

%if using coil info
if any(use_coils > 0)
    et_anal_name = [et_anal_name '_Cprior'];
end

% if strcmp(rec_type,'UA') && bar_ori == 90
%     et_mod_data_name = 'full_eyetrack_initmods_vbars';
%     et_anal_name = 'full_eyetrack_vbars';
%     save_name = 'corrected_models_vbars';
% end
et_mod_data_name = [et_mod_data_name sprintf('_ori%d',bar_ori)];
et_anal_name = [et_anal_name sprintf('_ori%d',bar_ori)];
save_name = [save_name sprintf('_ori%d',bar_ori)];

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

%if using trial-by-trial interleaving of expt conditions
if Expt_num >= 275
    is_TBT_expt = true;
else
    is_TBT_expt = false;
end

%%
xv_frac = 0.2; %fraction of trials to use for Xval
xv_type = 'uni';

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
%     case 296
%         full_nPix = 54;
end

%exclude data at beginning and end of each trial
beg_buffer = 0.2;
end_buffer = 0.05;
trial_dur = 4;

use_nPix_us = use_nPix*spatial_usfac;
klen_us = use_nPix_us*flen;

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

%% Spatial resolution
all_dws = cellfun(@(x) x.Stimvals.dw,Expts(cur_block_set));
base_sp_dx = mode(all_dws);
if length(unique(all_dws)) > 1
    fprintf('Warning, multiple different dws detected, using %.3f\n',base_sp_dx);
end
sp_dx = base_sp_dx/spatial_usfac/scale_fac; %model dx in deg
ov_RF_pos = Expts{cur_block_set(1)}.Stimvals.rf(1:2)/scale_fac;


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
all_trial_rptframes = [];

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
    cur_nrpt_frames = zeros(n_trials,1);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start'/1e4;
        n_frames = size(left_stim_mats{use_trials(tt)},1);
         if isfield(Expts{cur_block}.Trials(use_trials(tt)),'rptframes')
            cur_nrpt_frames(tt) = length(Expts{cur_block}.Trials(use_trials(tt)).rptframes);
        end
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
        all_trial_rptframes = [all_trial_rptframes; cur_nrpt_frames];

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

%select subset of pixels used for model fitting
buffer_pix = floor((full_nPix - use_nPix)/2);
[Xinds_up,~] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
cur_use_pix = (1/spatial_usfac:1/spatial_usfac:use_nPix) + buffer_pix;
use_kInds_up = find(ismember(Xinds_up(:),cur_use_pix));

stim_params_us = NMMcreate_stim_params([flen full_nPix_us],dt);
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

%% LOAD EYE-TRACKING DATA
cd(et_dir)
load(et_mod_data_name,'all_mod*');
load(et_anal_name,'drift*','it_*','et_tr_set','et_saccades');
tr_set = et_tr_set;

%% PROCESS EYE TRACKING DATA
cd(data_dir)

if isfield(Expts{cur_block_set(1)}.Header,'exptno')
    em_block_nums = cellfun(@(X) X.Header.exptno,Expts(cur_block_set),'uniformoutput',1); %block numbering for EM/LFP data sometimes isnt aligned with Expts struct
else
    em_block_nums = cur_block_set;
end

[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data_v2(all_t_axis,all_blockvec,em_block_nums,Expt_name,trial_toffset,good_coils);
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

if length(saccades) ~= length(et_saccades)
    fprintf('Using old blink detection algo!\n');
    [saccades,et_params] = detect_saccades_v2(corrected_eye_vals,all_eye_vals,all_eye_speed,all_eye_ts,et_params);
    is_blink = detect_blinks_old(all_eye_ts,all_eye_vals,saccades,et_params);
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
end

%% CREATE SACCADE PREDICTOR MATS
saccade_start_inds = find(ismember(used_inds,interp_sac_start_inds));
used_saccade_set = find(ismember(interp_sac_start_inds,used_inds));
%nearest index in the used data set of the saccade stop time
saccade_stop_inds = round(interp1(used_inds,1:length(used_inds),interp_sac_stop_inds(used_saccade_set)))';

saccade_stop_inds(isnan(saccade_stop_inds)) = length(used_inds);
used_is_blink = is_blink(used_saccade_set);

% sac_amps = [saccades(:).amplitude];
% is_micro = sac_amps(used_saccade_set) < 1;
% big_sacs = find(~is_micro & ~used_is_blink');
% micro_sacs = find(is_micro & ~used_is_blink');

saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));

%% DEFINE FIXATION PERIODS
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

%% Recon retinal stim for non LOO data
cur_fix_post_mean = squeeze(it_fix_post_mean(end,:));
cur_fix_post_std = squeeze(it_fix_post_std(end,:));
cur_drift_post_mean = squeeze(drift_post_mean(end,:));
cur_drift_post_std = squeeze(drift_post_std(end,:));
[fin_tot_corr,fin_tot_std] = construct_eye_position(cur_fix_post_mean,cur_fix_post_std,...
    cur_drift_post_mean,cur_drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);

fin_shift_cor = round(fin_tot_corr);

%RECOMPUTE XMAT
best_shift_stimmat_up = all_stimmat_up;
for i=1:NT
    best_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
end

if fit_unCor
    all_Xmat_unCor = create_time_embedding(all_stimmat_up,stim_params_us);
    all_Xmat_unCor = all_Xmat_unCor(used_inds,use_kInds_up);
end

%% COMBINE SUA AND MUA
n_chs = size(all_binned_mua,2) + size(all_binned_sua,2);

% make Robs_mat
poss_blocks = 1:(n_blocks-1);
Robs_mat = nan(length(used_inds),n_chs);
for ss = 1:n_chs
    if ss > n_probes
        su_probe_ind = ss-n_probes;
        Robs_mat(:,ss) = all_binned_sua(used_inds,su_probe_ind);
    else
        Robs_mat(:,ss) = all_binned_mua(used_inds,ss);
    end
end

%% IDENTIFY REPEAT TRIALS
%don't use repeat trials for training or cross-validation here.
if strcmp(xv_type,'rpt')
    rpt_trials = find(all_trial_Se==rpt_seed);
    rpt_inds = find(ismember(all_trialvec(used_inds),rpt_trials));
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

if strcmp(rec_type,'UA')
    probe_nums = [1:n_probes SU_probes];
    uset = find(~isnan(probe_nums));
    [all_ecc,all_Y,all_X] = deal(nan(1,length(probe_nums)));
    all_ecc(uset) = interp_ecc(probe_nums(uset));
    all_Y(uset) = interp_y(probe_nums(uset));
    all_X(uset) = interp_x(probe_nums(uset));
end

%%
% cd(save_dir)
% load(save_name);

%%
for cc = targs
    fprintf('Starting model fits for unit %d\n',cc);
    
    loo_cc = find(loo_set == cc);
    cc_uinds = find(~isnan(Robs_mat(:,cc)));
    
    if ~isempty(cc_uinds)
        use_trials = unique(all_trialvec(used_inds(cc_uinds)));
        use_trials(ismember(use_trials,rpt_trials)) = []; %DONT USE REPEAT TRIALS!
        
        nuse_trials = length(use_trials);
        n_xv_trials = round(xv_frac*nuse_trials);
        xv_trials = randperm(nuse_trials);
        xv_trials(n_xv_trials+1:end) = [];
        xv_trials = use_trials(xv_trials);
        tr_trials = setdiff(use_trials,xv_trials);
        n_tr_trials = length(tr_trials);
        fprintf('Initializing models with %d training trials and %d xval trials\n',n_tr_trials,n_xv_trials);
        
        cc_uinds(~ismember(all_trialvec(used_inds(cc_uinds)),use_trials)) = []; %eliminate repeat trials from the used indices
        cur_tr_inds = find(ismember(all_trialvec(used_inds(cc_uinds)),tr_trials));
        cur_xv_inds = find(ismember(all_trialvec(used_inds(cc_uinds)),xv_trials));
        cur_full_inds = find(ismember(all_trialvec(used_inds(cc_uinds)),use_trials));
        
        cur_Robs = Robs_mat(cc_uinds,cc);
        
        %% COMPUTE UNIT DATA
        unit_data.isLOO = ismember(cc,loo_set);
        unit_data.avg_rate = mean(cur_Robs)/dt;
        unit_data.tot_spikes = sum(cur_Robs);
        unit_data.N_used_samps = length(cur_Robs);
        used_blocks = unique(all_blockvec(used_inds(cc_uinds)));
        unit_data.n_used_blocks = length(used_blocks);
        block_rates = nan(unit_data.n_used_blocks,1);
        for ii = 1:unit_data.n_used_blocks
            block_rates(ii) = mean(cur_Robs(all_blockvec(used_inds(cc_uinds)) == used_blocks(ii)));
        end
        unit_data.rate_stability_cv = std(block_rates)/mean(block_rates);
        unit_data.block_rates = block_rates/dt;
        if cc > n_probes
            su_Ind = cc - n_probes;
            unit_data.SU_number = SU_numbers(su_Ind);
            unit_data.probe_number = SU_probes(su_Ind);
            unit_data.SU_Lratio = nanmean(Clust_data.SU_Lratios(su_Ind,used_blocks),2);
            unit_data.SU_isodist = nanmean(Clust_data.SU_isodists(su_Ind,used_blocks(Clust_data.SU_isoRel(su_Ind,used_blocks)==1)),2);
            unit_data.SU_refract = nanmean(Clust_data.SU_refract(su_Ind,used_blocks),2);
            unit_data.SU_dprime = nanmean(Clust_data.SU_dprimes(su_Ind,used_blocks),2);
            unit_data.dprime_stability_cv = nanstd(Clust_data.SU_dprimes(su_Ind,used_blocks),[],2)/nanmean(Clust_data.SU_dprimes(su_Ind,used_blocks),2);
        else
            unit_data.SU_number = nan;
            unit_data.probe_number = cc;
            unit_data.SU_Lratio = nan;
            unit_data.SU_isodist = nan;
            unit_data.SU_refract = nan;
            unit_data.SU_dprime = nan;
        end
        
        %% RECON RETINAL STIM
        if cc == targs(1) || cc > n_probes
            fprintf('Reconstructing retinal stim for unit %d\n',cc);
            if ismember(cc,loo_set)
                cur_fix_post_mean = squeeze(it_fix_post_mean_LOO(loo_cc,end,:));
                cur_fix_post_std = squeeze(it_fix_post_std_LOO(loo_cc,end,:));
                cur_drift_post_mean = squeeze(drift_post_mean_LOO(loo_cc,end,:));
                cur_drift_post_std = squeeze(drift_post_std_LOO(loo_cc,end,:));
                [fin_tot_corr,fin_tot_std] = construct_eye_position(cur_fix_post_mean,cur_fix_post_std,...
                    cur_drift_post_mean,cur_drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);
                
                fin_shift_cor = round(fin_tot_corr);
                
                %RECOMPUTE XMAT
                all_shift_stimmat_up = all_stimmat_up;
                for i=1:NT
                    all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
                end
                all_Xmat_shift = create_time_embedding(all_shift_stimmat_up,stim_params_us);
                all_Xmat_shift = all_Xmat_shift(used_inds(cc_uinds),use_kInds_up);
                
            else
                all_Xmat_shift = create_time_embedding(best_shift_stimmat_up,stim_params_us);
                all_Xmat_shift = all_Xmat_shift(used_inds(cc_uinds),use_kInds_up);
            end
        end
        
        %% FIT (eye-corrected) STIM-PROCESSING MODEL
        fprintf('Fitting stim models for unit %d\n',cc);
        
        max_Emods = 2;
        max_Imods = 2;
        
        silent = 1;
%         if mode(expt_dds(cur_block_set)) == 12
%             base_lambda_d2XT = 50;
%             base_lambda_L1 = 5;
%         else
%             base_lambda_d2XT = 200;
%             base_lambda_L1 = 10;
%         end
        if mode(expt_dds(cur_block_set)) == 12
            base_lambda_d2XT = 25;
            base_lambda_L1 = 5;
        else
            base_lambda_d2XT = 100;
            base_lambda_L1 = 10;
        end
        init_lambda_d2XT = 10;
        init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_lambda_d2XT,'boundary_conds',[0 0 0]);
        mod_stim_params(1) = NMMcreate_stim_params([flen use_nPix_us],dt);
        
        %start out with just a linear filter
        nEfilts = 0;
        nIfilts = 0;
        
        stim_mod_signs = [1 ones(1,nEfilts) -1*ones(1,nIfilts)];
        stim_Xtargs = [ones(1,length(stim_mod_signs))];
        stim_NL_types = [{'lin'} repmat({'quad'},1,nEfilts) repmat({'quad'},1,nIfilts)];
        
        init_mod = NMMinitialize_model( mod_stim_params, stim_mod_signs, stim_NL_types, init_reg_params,stim_Xtargs);
        init_mod = NMMfit_filters(init_mod,cur_Robs(cur_tr_inds),all_Xmat_shift(cur_tr_inds,:),[],[],silent);
        [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(init_mod,cur_Robs(cur_tr_inds),all_Xmat_shift(cur_tr_inds,:));
        vgint = var(gint);sgint = std(gint);
        init_mod = NMMadjust_regularization(init_mod,find(stim_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./vgint(stim_Xtargs==1));
        init_mod = NMMadjust_regularization(init_mod,find(stim_Xtargs==1),'lambda_L1',base_lambda_L1./sgint(stim_Xtargs==1));
        
        bestGQM = NMMfit_filters(init_mod,cur_Robs(cur_tr_inds),all_Xmat_shift(cur_tr_inds,:),[],[],silent);
        bestGQM_spkNL = NMMfit_logexp_spkNL(bestGQM,cur_Robs(cur_tr_inds),all_Xmat_shift(cur_tr_inds,:));
        if isempty(cur_xv_inds)
            error('Need xval trials');
        end
        [bestGQM_xvLL, ~, ~, ~, ~, ~, nullxvLL] = NMMmodel_eval(bestGQM_spkNL, cur_Robs(cur_xv_inds), all_Xmat_shift(cur_xv_inds,:));
        
        %now try adding up to max_Emods exc squared filters. Stop when xval
        %LL doesn't improve
        cur_imp = Inf;
        while nEfilts < max_Emods && cur_imp > 0
            cur_mod = bestGQM;
            cur_mod = NMMadd_NLinput(cur_mod,'quad',1,1,0.1*randn(flen*use_nPix_us,1),init_reg_params);
            nEfilts = nEfilts + 1;
            
            fprintf('Fitting model with LIN + %dE and %dI\n',nEfilts,nIfilts);
            cur_mod = NMMfit_filters(cur_mod,cur_Robs(cur_tr_inds),all_Xmat_shift(cur_tr_inds,:),[],[],silent);
            [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(cur_mod,cur_Robs(cur_tr_inds),all_Xmat_shift(cur_tr_inds,:));
            vgint = var(gint); sgint = std(gint);
            cur_mod = NMMadjust_regularization(cur_mod,length(cur_mod.mods),'lambda_d2XT',base_lambda_d2XT./vgint(end));
            cur_mod = NMMadjust_regularization(cur_mod,length(cur_mod.mods),'lambda_L1',base_lambda_L1./sgint(end));
            cur_mod = NMMfit_filters(cur_mod,cur_Robs(cur_tr_inds),all_Xmat_shift(cur_tr_inds,:),[],[],silent);
            cur_mod_spkNL = NMMfit_logexp_spkNL(cur_mod,cur_Robs(cur_tr_inds),all_Xmat_shift(cur_tr_inds,:));
            cur_mod_xvLL = NMMmodel_eval(cur_mod_spkNL, cur_Robs(cur_xv_inds), all_Xmat_shift(cur_xv_inds,:));
            
            cur_imp = cur_mod_xvLL - bestGQM_xvLL;
            if cur_imp > 0
                fprintf('Keeping new model\n');
                bestGQM = cur_mod;
                bestGQM_spkNL = cur_mod_spkNL;
                bestGQM_xvLL = cur_mod_xvLL;
            end
        end
        
        %now try adding up to max_Imods inh squared filters. Stop when xval
        %LL doesn't improve
        cur_imp = Inf;
        while nIfilts < max_Emods && cur_imp > 0
            cur_mod = bestGQM;
            cur_mod = NMMadd_NLinput(cur_mod,'quad',-1,1,0.1*randn(flen*use_nPix_us,1),init_reg_params);
            nIfilts = nIfilts + 1;
            
            fprintf('Fitting model with LIN + %dE and %dI\n',nEfilts,nIfilts);
            cur_mod = NMMfit_filters(cur_mod,cur_Robs(cur_tr_inds),all_Xmat_shift(cur_tr_inds,:),[],[],silent);
            [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(cur_mod,cur_Robs(cur_tr_inds),all_Xmat_shift(cur_tr_inds,:));
            vgint = var(gint); sgint = std(gint);
            cur_mod = NMMadjust_regularization(cur_mod,length(cur_mod.mods),'lambda_d2XT',base_lambda_d2XT./vgint(end));
            cur_mod = NMMadjust_regularization(cur_mod,length(cur_mod.mods),'lambda_L1',base_lambda_L1./sgint(end));
            cur_mod = NMMfit_filters(cur_mod,cur_Robs(cur_tr_inds),all_Xmat_shift(cur_tr_inds,:),[],[],silent);
            cur_mod_spkNL = NMMfit_logexp_spkNL(cur_mod,cur_Robs(cur_tr_inds),all_Xmat_shift(cur_tr_inds,:));
            cur_mod_xvLL = NMMmodel_eval(cur_mod_spkNL, cur_Robs(cur_xv_inds), all_Xmat_shift(cur_xv_inds,:));
            
            cur_imp = cur_mod_xvLL - bestGQM_xvLL;
            if cur_imp > 0
                fprintf('Keeping new model\n');
                bestGQM = cur_mod;
                bestGQM_spkNL = cur_mod_spkNL;
                bestGQM_xvLL = cur_mod_xvLL;
            end
        end
        
        bestGQM.xvLLimp = (bestGQM_xvLL-nullxvLL)/log(2);
        
        %% Fit a GQM with the linear filter split into two rectified subunits
        %convert linear filter into separate thresh-lin filters and refit
        rectGQM = bestGQM;
        rectGQM.mods = [rectGQM.mods; rectGQM.mods(1)];
        rectGQM.mods(end).sign = -1;
        rectGQM.mods(end).filtK = -rectGQM.mods(end).filtK; %initialize filter to be negative of the linear filter
        rectGQM.mods(1).NLtype = 'threshlin'; rectGQM.mods(end).NLtype = 'threshlin';
        rectGQM = NMMfit_filters(rectGQM,cur_Robs(cur_tr_inds),all_Xmat_shift(cur_tr_inds,:),[],[],silent);
        rectGQM_spkNL = NMMfit_logexp_spkNL(rectGQM,cur_Robs(cur_tr_inds),all_Xmat_shift(cur_tr_inds,:));
        [rectGQM_xvLL] = NMMmodel_eval(rectGQM_spkNL, cur_Robs(cur_xv_inds), all_Xmat_shift(cur_xv_inds,:));
        rectGQM.xvLLimp = (rectGQM_xvLL-nullxvLL)/log(2);
        
        %% Fit uncorrected model
        if fit_unCor
            bestGQM_unCor = bestGQM;
            bestGQM_unCor = NMMfit_filters(bestGQM_unCor,cur_Robs(cur_tr_inds),all_Xmat_unCor(cc_uinds(cur_tr_inds),:));
            cur_mod_spkNL = NMMfit_logexp_spkNL(bestGQM_unCor,cur_Robs(cur_tr_inds),all_Xmat_unCor(cc_uinds(cur_tr_inds),:));
            cur_xvLL = NMMmodel_eval(cur_mod_spkNL,cur_Robs(cur_xv_inds),all_Xmat_unCor(cc_uinds(cur_xv_inds),:));
            bestGQM_unCor.xvLLimp = (cur_xvLL - nullxvLL)/log(2);
            
            rectGQM_unCor = rectGQM;
            rectGQM_unCor = NMMfit_filters(rectGQM_unCor,cur_Robs(cur_tr_inds),all_Xmat_unCor(cc_uinds(cur_tr_inds),:));
            cur_mod_spkNL = NMMfit_logexp_spkNL(rectGQM_unCor,cur_Robs(cur_tr_inds),all_Xmat_unCor(cc_uinds(cur_tr_inds),:));
            cur_xvLL = NMMmodel_eval(cur_mod_spkNL,cur_Robs(cur_xv_inds),all_Xmat_unCor(cc_uinds(cur_xv_inds),:));
            rectGQM_unCor.xvLLimp = (cur_xvLL - nullxvLL)/log(2);
        end
        
        %% Refit model parameters on the full data set
        if ~isempty(cur_xv_inds)
            bestGQM = NMMfit_filters(bestGQM,cur_Robs(cur_full_inds),all_Xmat_shift(cur_full_inds,:),[],[],silent);
            bestGQM_spkNL = NMMfit_logexp_spkNL(bestGQM,cur_Robs(cur_full_inds),all_Xmat_shift(cur_full_inds,:));
            [bestGQM_LL, ~, ~, ~, ~, ~, nullLL] = NMMmodel_eval(bestGQM_spkNL, cur_Robs(cur_full_inds), all_Xmat_shift(cur_full_inds,:));
            bestGQM.LLimp = (bestGQM_LL - nullLL)/log(2);
            
            if fit_unCor
                bestGQM_unCor = NMMfit_filters(bestGQM_unCor,cur_Robs(cur_full_inds),all_Xmat_unCor(cc_uinds(cur_full_inds),:));
                bestGQM_unCor_spkNL = NMMfit_logexp_spkNL(bestGQM_unCor,cur_Robs(cur_full_inds),all_Xmat_unCor(cc_uinds(cur_full_inds),:));
                bestGQM_unCor.LLimp = (bestGQM_unCor_spkNL.LL_seq(end)-nullLL)/log(2);
            end
            rectGQM = NMMfit_filters(rectGQM,cur_Robs(cur_full_inds),all_Xmat_shift(cur_full_inds,:),[],[],silent);
            rectGQM_spkNL = NMMfit_logexp_spkNL(rectGQM,cur_Robs(cur_full_inds),all_Xmat_shift(cur_full_inds,:));
            rectGQM.LLimp = (rectGQM_spkNL.LL_seq(end)-nullLL)/log(2);
            
            if fit_unCor
                rectGQM_unCor = NMMfit_filters(rectGQM_unCor,cur_Robs(cur_full_inds),all_Xmat_unCor(cc_uinds(cur_full_inds),:));
                rectGQM_unCor_spkNL = NMMfit_logexp_spkNL(rectGQM_unCor,cur_Robs(cur_full_inds),all_Xmat_unCor(cc_uinds(cur_full_inds),:));
                rectGQM_unCor.LLimp = (rectGQM_unCor_spkNL.LL_seq(end)-nullLL)/log(2);
            end
        end
        
        %% Store data
        fullLL = nansum(cur_Robs(cur_full_inds).*log(cur_Robs(cur_full_inds)) - cur_Robs(cur_full_inds));
        fullxvLL = nansum(cur_Robs(cur_xv_inds).*log(cur_Robs(cur_xv_inds)) - cur_Robs(cur_xv_inds));
        unit_data.fullLL = fullLL;
        unit_data.nullLL = nullLL*nansum(cur_Robs(cur_full_inds));
        unit_data.fullxvLL = fullxvLL;
        unit_data.nullxvLL = nullxvLL*nansum(cur_Robs(cur_xv_inds));
        unit_data.tr_trials = tr_trials;
        unit_data.xv_trials = xv_trials;
        
        ModData(cc).unit_data = unit_data;
        ModData(cc).bestGQM = bestGQM;
        if fit_unCor
            ModData(cc).bestGQM_unCor = bestGQM_unCor;
        end
        ModData(cc).rectGQM = rectGQM;
        if fit_unCor
            ModData(cc).rectGQM_unCor = rectGQM_unCor;
        end
        
        %% CALCULATE TUNING PROPERTIES
        clear tune_props
        rectGQM = ModData(cc).rectGQM;
        
        stim_filters = [rectGQM.mods(:).filtK];
        stim_mod_signs = [rectGQM.mods(:).sign];
        mod_stim_params = rectGQM.stim_params;
        filt_data = get_filter_properties_v2(stim_filters,mod_stim_params,sp_dx);
        tune_props.filt_data = filt_data;
        
        [~,~,best_pred_rate,~,~,filt_outs] = NMMmodel_eval(rectGQM,cur_Robs,all_Xmat_shift);
        [~,~,rev_pred_rate] = NMMmodel_eval(rectGQM,cur_Robs,-all_Xmat_shift);
        tune_props.PRM = mean(abs(best_pred_rate - rev_pred_rate))/mean(best_pred_rate);
        
        %only use E/lin filters for these calculations
        non_supp_filts = find(stim_mod_signs ~= -1);
        filt_out_weights = std(filt_outs);
        filt_out_weights = filt_out_weights(non_supp_filts)/nansum(filt_out_weights(non_supp_filts));
        
        tune_props.RF_mean = nansum(filt_out_weights(non_supp_filts).*filt_data.gest(non_supp_filts,1)') - use_nPix_us*sp_dx/2 -sp_dx;
        tune_props.RF_sigma = nansum(filt_out_weights(non_supp_filts).*filt_data.gest(non_supp_filts,2)');
        tune_props.RF_gSF = nansum(filt_out_weights(non_supp_filts).*filt_data.gest(non_supp_filts,3)');
        tune_props.RF_dirsel = nansum(filt_out_weights(non_supp_filts).*filt_data.dir_selectivity(non_supp_filts)');
        tune_props.RF_FTF = nansum(filt_out_weights(non_supp_filts).*filt_data.FFt(non_supp_filts)');
        tune_props.RF_FSF = nansum(filt_out_weights(non_supp_filts).*filt_data.FFx(non_supp_filts)');
        
        rect_stim_filters = [rectGQM.mods([1 end]).filtK];
        avg_rect_filts = mean(rect_stim_filters);
        absavg_rect_filts = mean(abs(rect_stim_filters(:)));
        tune_props.net_phase_polarity = (avg_rect_filts(1) - avg_rect_filts(end))/absavg_rect_filts;
        
        screen_X =  -tune_props.RF_mean'.*sind(bar_ori) + ov_RF_pos(1);
        screen_Y =  tune_props.RF_mean'.*cosd(bar_ori) + ov_RF_pos(2);
        if strcmp(rec_type,'UA')
            if bar_ori == 0
                screen_X = interp_x(unit_data.probe_number);
            elseif bar_ori == 90
                screen_Y = interp_y(unit_data.probe_number);
            else
                error('Code doesnt work with this condition!');
            end
        end
        
        tune_props.screen_X = screen_X;
        tune_props.screen_Y = screen_Y;
        tune_props.RF_ecc = sqrt(screen_X.^2 + screen_Y.^2);
        
        ModData(cc).tune_props = tune_props;
 
        %% CALCULATE UNCORRECTED TUNING PROPERTIES
        clear tune_props_unCor
        if fit_unCor
            rectGQM_unCor = ModData(cc).rectGQM_unCor;
            
            stim_filters = [rectGQM_unCor.mods(:).filtK];
            stim_mod_signs = [rectGQM_unCor.mods(:).sign];
            mod_stim_params = rectGQM_unCor.stim_params;
            filt_data = get_filter_properties_v2(stim_filters,mod_stim_params,sp_dx);
            tune_props_unCor.filt_data = filt_data;
            
            [~,~,best_pred_rate,~,~,filt_outs] = NMMmodel_eval(rectGQM_unCor,cur_Robs,all_Xmat_shift);
            [~,~,rev_pred_rate] = NMMmodel_eval(rectGQM_unCor,cur_Robs,-all_Xmat_shift);
            tune_props_unCor.PRM = mean(abs(best_pred_rate - rev_pred_rate))/mean(best_pred_rate);
            
            %only use E/lin filters for these calculations
            non_supp_filts = find(stim_mod_signs ~= -1);
            filt_out_weights = std(filt_outs);
            filt_out_weights = filt_out_weights(non_supp_filts)/nansum(filt_out_weights(non_supp_filts));
            
            tune_props_unCor.RF_mean = filt_out_weights(non_supp_filts)*filt_data.gest(non_supp_filts,1) - use_nPix_us*sp_dx/2 -sp_dx;
            tune_props_unCor.RF_sigma = filt_out_weights(non_supp_filts)*filt_data.gest(non_supp_filts,2);
            tune_props_unCor.RF_gSF = filt_out_weights(non_supp_filts)*filt_data.gest(non_supp_filts,3);
            tune_props_unCor.RF_dirsel = filt_out_weights(non_supp_filts)*filt_data.dir_selectivity(non_supp_filts);
            tune_props_unCor.RF_FTF = filt_out_weights(non_supp_filts)*filt_data.FFt(non_supp_filts);
            tune_props_unCor.RF_FSF = filt_out_weights(non_supp_filts)*filt_data.FFx(non_supp_filts);
            
            rect_stim_filters = [rectGQM_unCor.mods([1 end]).filtK];
            avg_rect_filts = mean(rect_stim_filters);
            absavg_rect_filts = mean(abs(rect_stim_filters(:)));
            tune_props_unCor.net_phase_polarity = (avg_rect_filts(1) - avg_rect_filts(end))/absavg_rect_filts;
            
            screen_X =  -tune_props_unCor.RF_mean'.*sind(bar_ori) + ov_RF_pos(1);
            screen_Y =  tune_props_unCor.RF_mean'.*cosd(bar_ori) + ov_RF_pos(2);
            if strcmp(rec_type,'UA')
                if bar_ori == 0
                    screen_X = interp_x(unit_data.probe_number);
                elseif bar_ori == 90
                    screen_Y = interp_y(unit_data.probe_number);
                else
                    error('Code doesnt work with this condition!');
                end
            end
            
            tune_props_unCor.screen_X = screen_X;
            tune_props_unCor.screen_Y = screen_Y;
            tune_props_unCor.RF_ecc = sqrt(screen_X.^2 + screen_Y.^2);
            
            ModData(cc).tune_props_unCor = tune_props_unCor;
        end
    end
end

%%
modFitParams = struct('bar_ori',bar_ori,'use_MUA',use_MUA,'fit_uncor',fit_unCor,...
    'xv_frac',xv_frac,'xv_type',xv_type,'use_nPix_us',use_nPix_us,'flen',flen,'spatial_usfac',spatial_usfac,'sp_dx',sp_dx,'dt',dt,...
    'base_lambda_d2XT',base_lambda_d2XT,'base_lambda_L1',base_lambda_L1,...
    'init_lambda_d2XT',init_lambda_d2XT);

cd(save_dir)
save(save_name,'ModData','modFitParams');