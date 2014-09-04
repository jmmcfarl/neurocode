%
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements//');
addpath('~/James_scripts/bruce/processing/');
addpath('~/James_scripts/bruce/saccade_modulation/');
addpath('~/James_scripts/TentBasis2D/');

global Expt_name bar_ori use_MUA

% Expt_name = 'M297';
Expt_name = 'M296';
use_MUA = false;
bar_ori = 0; %bar orientation to use (only for UA recs)

mod_data_name = 'corrected_models2';

%%

micro_thresh = 1; %max amp of microsac (deg)
EP_bounds = 1;%eye position boundary (deg from central FP)
sac_burst_isi = 0.15;
max_gsac_dur = 0.1;

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

anal_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
if ~exist(anal_dir,'dir')
    mkdir(anal_dir)
end

et_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
mod_data_dir = ['~/Analysis/bruce/' Expt_name '/models'];

% et_mod_data_name = 'full_eyetrack_initmods';
% et_anal_name = 'full_eyetrack';
et_mod_data_name = 'full_eyetrack_initmods_Rinit';
et_anal_name = 'full_eyetrack_Rinit';

%if using coil info
if any(use_coils > 0)
    et_anal_name = [et_anal_name '_Cprior'];
end

et_mod_data_name = [et_mod_data_name sprintf('_ori%d',bar_ori)];
et_anal_name = [et_anal_name sprintf('_ori%d',bar_ori)];
mod_data_name = [mod_data_name sprintf('_ori%d',bar_ori)];

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
elseif ismember(Expt_num,[296 297])
    use_nPix = 22;
    spatial_usfac = 2;
end

min_trial_dur = 0.75;

stim_fs = 100; %in Hz
dt = 0.01;
new_dt = 0.0025;
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
        %     case 296
        %         full_nPix = 54;
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


if ~isempty(grayback_gs_expts) || ~isempty(imback_gs_expts)
    gsac_amp = unique(expt_sac_amp(cur_block_set([grayback_gs_expts; imback_gs_expts])));
else
    gsac_amp = unique(expt_sac_amp(cur_block_set));
end
if length(gsac_amp) > 1
    fprintf('Multiple guided sac amps detected!\n');
end
%minimum (parallel) amplitude for a guided saccade to be included in
%analysis
gsac_thresh = mean(gsac_amp)/2;

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
all_t_axis_new = [];
all_t_bin_edges_new = [];
all_bin_edge_pts_new = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_Se = [];
all_trial_wi = [];
all_trial_back = [];
all_trial_Ff = [];
all_trial_exvals = [];
all_trial_blk = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_trial_rptframes = [];
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
    
    if is_TBT_expt
        if isfield(Expts{cur_block}.Trials,'Bs')
            trial_back = strcmp('image',{Expts{cur_block}.Trials(:).Bs});
        else
            trial_back = nan(1,length(use_trials));
        end
        trial_back = trial_back(id_inds);
        
        if isfield(Expts{cur_block}.Trials,'exvals')
            exvals = reshape([Expts{cur_block}.Trials(:).exvals],length(Expts{cur_block}.Trials(1).exvals),[]);
            trial_exvals = exvals(:,id_inds)';
        else
            trial_exvals = nan(length(id_inds),3);
        end
        
        if isfield(Expts{cur_block}.Trials,'Ff')
            trial_Ff = [Expts{cur_block}.Trials(:).Ff];
            trial_Ff = trial_Ff(id_inds);
        else
            trial_Ff = nan(1,length(id_inds));
        end
        all_trial_back = cat(1,all_trial_back,trial_back(use_trials)');
        all_trial_Ff = cat(1,all_trial_Ff,trial_Ff(use_trials)');
        all_trial_exvals = cat(1,all_trial_exvals,trial_exvals(use_trials,:));
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
        
        cur_t_edges_new = cur_t_edges(1):new_dt:cur_t_edges(end);
        cur_t_axis_new = 0.5*cur_t_edges_new(1:end-1) + 0.5*cur_t_edges_new(2:end);
        
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
            
            all_t_axis_new = [all_t_axis_new; cur_t_axis_new' + cur_toffset];
            all_t_bin_edges_new = [all_t_bin_edges_new; cur_t_edges_new' + cur_toffset];
            all_bin_edge_pts_new = [all_bin_edge_pts_new; length(all_t_bin_edges_new)];
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

if is_TBT_expt
    if any(~isnan(all_trial_exvals(:,3)))
        fprintf('Using exvals to define trial-by-trial conditions\n');
        all_trial_Ff(all_trial_exvals(:,3) == 1) = 0; %these are sim sac trials
        all_trial_Ff(all_trial_exvals(:,3) > 1) = 70;
        all_trial_back(all_trial_exvals(:,3) == 2) = 0; %these are gray-back trials
        all_trial_back(ismember(all_trial_exvals(:,3),[1 3])) = 1;
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
[all_binned_mua,all_binned_sua,Clust_data,all_su_spk_times] = ...
    get_binned_spikes(cluster_dir,all_spk_times,all_clust_ids,all_spk_inds,...
    all_t_axis,all_t_bin_edges,all_bin_edge_pts,cur_block_set,all_blockvec,clust_params);
SU_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;

all_blockvec_new = round(interp1(all_t_axis,all_blockvec,all_t_axis_new));
[all_binned_mua_new,all_binned_sua_new] = ...
    get_binned_spikes(cluster_dir,all_spk_times,all_clust_ids,all_spk_inds,...
    all_t_axis_new,all_t_bin_edges_new,all_bin_edge_pts_new,cur_block_set,all_blockvec_new,clust_params);

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

%%
fprintf('Loading ET data\n');
cd(et_dir)
load(et_mod_data_name,'all_mod*');
load(et_anal_name,'drift*','it_*','et_tr_set','et_saccades');
tr_set = et_tr_set;

fprintf('Loading model fits\n');
cd(mod_data_dir)
load(mod_data_name);

%% PROCESS EYE TRACKING DATA
cd(data_dir)

if isfield(Expts{cur_block_set(1)}.Header,'exptno')
    em_block_nums = cellfun(@(X) X.Header.exptno,Expts(cur_block_set),'uniformoutput',1); %block numbering for EM/LFP data sometimes isnt aligned with Expts struct
else
    em_block_nums = cur_block_set;
end

% [all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);
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
saccades = saccades(used_saccade_set);

sac_amps = [saccades(:).amplitude];
sac_direction = [saccades(:).direction];
sac_durs = [saccades(:).duration];
sac_prepos = reshape([saccades(:).pre_pos],[],length(saccades));
sac_postpos = reshape([saccades(:).post_pos],[],length(saccades));
sac_deltaX = sac_postpos(1,:) - sac_prepos(1,:);

%find saccades that start or end out of window
out_bounds = abs(sac_prepos(2,:)) > EP_bounds | abs(sac_postpos(2,:)) > EP_bounds;

sacburst_set = find([saccades(:).isi] < sac_burst_isi | [saccades(:).next_isi] < sac_burst_isi);
micro_sacs = find([saccades(:).amplitude] < micro_thresh & ~used_is_blink' & ~out_bounds);

msac_bursts = micro_sacs(ismember(micro_sacs,sacburst_set));
micro_sacs(ismember(micro_sacs,sacburst_set)) = []; %eliminate microsacs that are part of a 'burst'

%guided saccades are those whose parallel component is large enough and
%that aren't blinks (and whose duration is not too long to be suspicious
big_sacs = find(abs(sac_deltaX) > gsac_thresh & ~used_is_blink' & ~out_bounds & sac_durs <= max_gsac_dur);

saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));

msac_thresh = prctile(sac_amps(micro_sacs),50);
big_msacs = micro_sacs(sac_amps(micro_sacs) > msac_thresh);
small_msacs = micro_sacs(sac_amps(micro_sacs) < msac_thresh);
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
%%
use_vec = zeros(size(all_t_axis));
use_vec(used_inds) = 1;
use_vec_new = round(interp1(all_t_axis,use_vec,all_t_axis_new));
use_vec_new(isnan(use_vec_new)) = 0;
use_vec_new = logical(use_vec_new);
%% Create set of TR trials
rpt_trials = find(all_trial_Se==rpt_seed);
n_rpt_trials = length(rpt_trials);

use_trials = unique(all_trialvec(used_inds));
use_trials(ismember(use_trials,rpt_trials)) = []; %DONT USE REPEAT TRIALS!
nuse_trials = length(use_trials);

xv_frac = 0.2;
n_xv_trials = round(xv_frac*nuse_trials);
xv_trials = randperm(nuse_trials);
xv_trials(n_xv_trials+1:end) = [];
xv_trials = use_trials(xv_trials);
tr_trials = setdiff(use_trials,xv_trials);
n_tr_trials = length(tr_trials);
fprintf('Initializing models with %d training trials and %d xval trials\n',n_tr_trials,n_xv_trials);

tr_inds = find(ismember(all_trialvec(used_inds),tr_trials));
xv_inds = find(ismember(all_trialvec(used_inds),xv_trials));
full_inds = find(ismember(all_trialvec(used_inds),use_trials));
%% CREATE SACCADE AND MICROSAC INDICATOR XMATS

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
%%
rpt_taxis = (1:round(trial_dur)/dt)*dt-dt/2;
rpt_taxis(rpt_taxis < beg_buffer) = [];
rpt_taxis(trial_dur - rpt_taxis < end_buffer) = [];
all_trial_dur = all_trial_end_times-all_trial_start_times;

rpt_trials = find(all_trial_Se == rpt_seed);
rpt_trials(all_trial_rptframes(rpt_trials) > 0) = []; %get rid of any repeat trials where there were repeat frames

n_rpts = length(rpt_trials);
all_rpt_inds = find(ismember(all_trialvec(used_inds),rpt_trials));
all_nonrpt_inds = find(~ismember(all_trialvec(used_inds),rpt_trials));

%%
full_psth = nan(n_rpts,length(rpt_taxis),length(targs));
cor_pred_rate = nan(length(all_rpt_inds),length(targs));
for cc = targs
    
    fprintf('Starting model fits for unit %d\n',cc);
    loo_cc = find(loo_set == cc); %index within the LOOXV set
    cc_uinds = full_inds(~isnan(Robs_mat(full_inds,cc))); %set of used indices where this unit was isolated
    cur_Robs = Robs_mat(cc_uinds,cc);
    
    if ~isempty(cc_uinds)
        fprintf('Unit %d of %d\n',cc,length(targs));
        
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
        
    else %otherwise use overall EP sequence
        all_Xmat_shift = create_time_embedding(best_shift_stimmat_up,stim_params_us);
        all_Xmat_shift = all_Xmat_shift(used_inds(cc_uinds),use_kInds_up);
    end
    
    all_used_rpt_inds = all_rpt_inds(~isnan(Robs_mat(all_rpt_inds,cc)));
    used_rpt_inds = find(ismember(all_rpt_inds,all_used_rpt_inds));
    
    %%
    su_num = find(SU_numbers == all_mod_SUnum(cc));
    all_spk_trials{cc} = [];
    all_rel_times{cc} = [];
    for ii = 1:length(rpt_trials)
        cur_inds = find(all_trialvec(used_inds) == rpt_trials(ii));
        cur_inds(size(full_psth,2)+1:end) = [];
        full_psth(ii,1:length(cur_inds),cc) = Robs_mat(cur_inds,cc);
        
        if ~isempty(su_num)
            cur_spk_set = find(all_su_spk_times{su_num} > all_trial_start_times(rpt_trials(ii)) & ...
                all_su_spk_times{su_num} < all_trial_end_times(rpt_trials(ii)));
            spk_times = all_su_spk_times{su_num}(cur_spk_set);
            all_spk_trials{cc} = [all_spk_trials{cc}; ones(size(cur_spk_set))*ii];
            all_rel_times{cc} = [all_rel_times{cc}; spk_times-all_trial_start_times(rpt_trials(ii))];
        end
    end
    
    %%    
    %CORRECTED AND UNCORRECTED MODEL FITS
    cur_cormod = ModData(cc).rectGQM;
    cur_Xtargs = [cur_cormod.mods(:).Xtarget];
    cur_cormod = NMMfit_logexp_spkNL(cur_cormod,Robs_mat(all_used_rpt_inds,cc),all_Xmat_shift(all_rpt_inds,:),[],[],[2; 3]);
    [corrLL, penLL, cor_pred_rate(:,cc),G,~,~,nullLL] = NMMmodel_eval(cur_cormod,Robs_mat(all_rpt_inds,cc),all_X);
    var_data.cor_mods(cc) = cur_cormod;
    
    neur_eye_pos(:,:,cc) = nan(length(rpt_taxis),length(rpt_trials));
    for ii = 1:length(rpt_trials)
        test_rpt_inds = find(all_trialvec(used_inds(all_rpt_inds)) == rpt_trials(ii));
        test_rpt_inds(length(rpt_taxis)+1:end) = [];
        neur_eye_pos(1:length(test_rpt_inds),ii,cc) = fin_tot_corr_neural(all_rpt_inds(test_rpt_inds));
    end
    
    %% calculate model-predicted EM variance during single-trial data
    %     %NOTE< SHOULD IMPLEMENT THIS IN A PARALLELIZED FASHION FOR SPEED>>>
    %     all_used_nonrpt_inds = all_nonrpt_inds(~isnan(Robs_mat(all_nonrpt_inds,cur_runit_ind)));
    %     used_nonrpt_inds = find(ismember(all_nonrpt_inds,all_used_nonrpt_inds));
    %
    %     all_X{1} = drift_Xmat(all_used_nonrpt_inds,:);
    %     cur_cormod = dit_mods_spkNL_LOO{cur_runit_ind,end}(cur_unit_ind);
    %     cur_cormod.mods(cur_Xtargs > 1) = [];
    %     cur_cormod = NMMfit_logexp_spkNL(cur_cormod,Robs_mat(used_nonrpt_inds,cur_runit_ind),all_X{1},[],[],[1 2]);
    % %     all_X{2} = Xblock(used_inds(all_nonrpt_inds),:);
    % %     all_X{3} = Xsac(all_nonrpt_inds,:);
    % %     all_X{4} = Xmsac(all_nonrpt_inds,:);
    %     rate_out = nan(n_shifts,length(all_nonrpt_inds));
    %     for ii = 1:n_shifts
    % %         fprintf('Shift %d of %d\n',ii,n_shifts);
    %         all_X{1} = all_Xmat_us(used_inds(all_nonrpt_inds),:)*shift_mat{ii};
    %         [~,~,rate_out(ii,:)] = NMMmodel_eval(cur_cormod,[],all_X{1});
    %     end
    %     disc_fin_tot_corr = round(fin_tot_corr/sp_dx);
    %     uu = find(disc_fin_tot_corr >= -max_shift & disc_fin_tot_corr <= max_shift);
    %     ep_dist = hist(disc_fin_tot_corr(uu),shifts);
    %     ep_dist = ep_dist/sum(ep_dist);
    %
    %     cond_mrate = sum(bsxfun(@times,rate_out(:,used_nonrpt_inds),ep_dist'));
    %     rate_out_ms = bsxfun(@minus,rate_out(:,used_nonrpt_inds),cond_mrate);
    %     cond_var = sum(bsxfun(@times,rate_out_ms.^2,ep_dist'));
    %     em_var = mean(cond_var);
    %
    %     ov_mrate = mean(Robs_mat(used_nonrpt_inds,cur_runit_ind));
    %     tot_var = mean(sum(bsxfun(@times,(rate_out(:,used_nonrpt_inds)-ov_mrate).^2,ep_dist')));
    %     var_data.st_mod_emvar(ss) = em_var;
    %     var_data.st_mod_psthvar(ss) = nanvar(cond_mrate);
    %     var_data.st_mod_totvar(ss) = tot_var;
    %     var_data.st_mod_emvfrac(ss) = em_var/tot_var;
    
end

%% Compute overall eye pos estimate
fin_fix_corr = nan(NT,1);
fin_fix_std = nan(NT,1);
fin_fix_corr(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
fin_fix_std(~isnan(fix_ids)) = it_fix_post_std(end,fix_ids(~isnan(fix_ids)));
fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);
fin_fix_corr = fin_fix_corr*sp_dx;
fin_fix_std = fin_fix_std*sp_dx;
fin_drift_corr = squeeze(drift_post_mean(end,:)*sp_dx);
fin_drift_std = squeeze(drift_post_std(end,:)*sp_dx);
fin_drift_corr_neural = fin_drift_corr;

for ii = 1:length(trial_start_inds)
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
    fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
end
fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT);
fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);

fin_tot_corr = fin_fix_corr + fin_drift_corr;
fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2);

for ii = 1:length(fix_start_inds)
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > sac_shift
        fin_drift_corr_neural(cur_inds(1:sac_shift)) = nan;
    else
        fin_drift_corr_neural(cur_inds) = nan;
    end
end

for ii = 1:length(fix_start_inds)-1
    cur_inds = fix_stop_inds(ii):fix_start_inds(ii+1);
    fin_drift_corr_neural(cur_inds) = nan;
end
fin_tot_corr_neural = fin_fix_corr + fin_drift_corr_neural;

tot_neur_eye_pos = nan(length(rpt_taxis),length(rpt_trials));
for ii = 1:length(rpt_trials)
    test_rpt_inds = find(all_trialvec(used_inds(all_rpt_inds)) == rpt_trials(ii));
    test_rpt_inds(length(rpt_taxis)+1:end) = [];
    tot_neur_eye_pos(1:length(test_rpt_inds),ii) = fin_tot_corr_neural(all_rpt_inds(test_rpt_inds));
end


%% COMPUTE TRIAL-BY-TRIAL PREDICTED RATES
tbt_cor_pred_rate = nan(length(rpt_taxis),length(rpt_trials),length(targs));
% tbt_corsac_pred_rate = nan(length(rpt_taxis),length(rpt_trials),length(tr_SUs));
% tbt_sac_pred_rate = nan(length(rpt_taxis),length(rpt_trials),length(tr_SUs));
% tbt_uncor_pred_rate = nan(length(rpt_taxis),length(rpt_trials),length(tr_SUs));
inf_eye_pos_t = nan(length(rpt_taxis),length(rpt_trials));
row_ind = nan(length(rpt_trials),1);
for ii = 1:length(rpt_trials)
    row_ind(ii) = ii;
    test_rpt_inds = find(all_trialvec(used_inds(all_rpt_inds)) == rpt_trials(ii));
    test_rpt_inds(length(rpt_taxis)+1:end) = [];
    tbt_cor_pred_rate(1:length(test_rpt_inds),ii,:) = cor_pred_rate(test_rpt_inds,:);
%     tbt_corsac_pred_rate(1:length(test_rpt_inds),ii,:) = corsac_pred_rate(test_rpt_inds,1:length(tr_SUs));
%     tbt_sac_pred_rate(1:length(test_rpt_inds),ii,:) = sac_pred_rate(test_rpt_inds,1:length(tr_SUs));
%     tbt_uncor_pred_rate(1:length(test_rpt_inds),ii,:) = uncor_pred_rate(test_rpt_inds,1:length(tr_SUs));    
    inf_eye_pos_t(1:length(test_rpt_inds),ii) = fin_tot_corr(all_rpt_inds(test_rpt_inds));
end

[sorted_inf_eyepos,inf_eyepos_ord] = sort(inf_eye_pos_t);

tr_sub_psth = bsxfun(@minus,full_psth,nanmean(full_psth,2));
avg_rates = nanmean(full_psth);
% avg_rates = nanmean(tr_sub_psth);
tbt_psth_mod_rate = permute(repmat(avg_rates,[length(rpt_trials),1,1]),[2 1 3]);

%% initialize parameters for spline models
use_eye_range = [-0.4 0.4];
n_pts = 200;
eye_xx = linspace(use_eye_range(1),use_eye_range(2),n_pts);

n_knots = 25;
knot_pts = linspace(-0.4,0.4,n_knots);

% all_trials = 1:n_rpts;
all_trials = 1:n_rpts;
xv_frac = 0.2;
xv_trials = randperm(length(all_trials));
xv_trials(round(xv_frac*n_rpts)+1:end) = [];
xv_trials = all_trials(xv_trials);
tr_trials = setdiff(all_trials,xv_trials);
[RR,TT] = meshgrid(1:length(rpt_taxis),1:n_rpts);

%% estimate spline models
cnt = 0;
eps = 0.1;

poss_lambdas = [1 5 10 50 100 500 1000 5e3];
% poss_lambdas = [1 2 5 10 20];
% for ll = 1:length(poss_lambdas)
ll = 3;
lambda = poss_lambdas(ll);
% lambda = 0.5;
% trial_avg_rates = squeeze(nanmean(full_psth(:,:,cc),2));

n_tr_chs = length(targs);
smoothed_eye_obs_fun = nan(length(eye_xx),length(rpt_taxis),n_tr_chs);
np_prate = nan(length(rpt_trials),length(rpt_taxis),n_tr_chs);
% for cc = 1:n_tr_chs
cc = 3
sp = NMMcreate_stim_params([1 n_knots]);
spr = NMMcreate_reg_params('lambda_d2X',lambda,'boundary_conds',[Inf Inf Inf]);

fprintf('Unit %d of %d\n',cc,length(n_tr_chs));
for ii = 1:length(rpt_taxis)
    %         fprintf('Time %d of %d\n',ii,length(rpt_taxis));
    xx = squeeze(neur_eye_pos(ii,:,cc))';
    xx(xx > use_eye_range(2) | xx < use_eye_range(1)) = nan;
    yy_obs = squeeze(full_psth(:,ii,cc));
    
    uu = all_trials(~isnan(xx(all_trials)) & ~isnan(yy_obs(all_trials)));
    %         knot_pts = prctile(xx(uu),linspace(0,100,n_knots));
    
    uu_tr_set = ismember(TT(uu,ii),tr_trials);
    uu_tr = uu(ismember(TT(uu,ii),tr_trials));
    uu_xv = uu(ismember(TT(uu,ii),xv_trials));
    if sum(yy_obs(uu_tr)>0) < 3
        np_prate(uu,ii,cc) = nan(size(uu));
        smoothed_eye_obs_fun(:,ii) = 0;
    else
        temp = fnval(csape(knot_pts,eye(length(knot_pts)),'var'),xx(uu));
        mod = NMMinitialize_model(sp,[1],{'lin'},spr,[],[],'logexp');
        mod = NMMfit_filters(mod,yy_obs(uu_tr),temp(:,uu_tr_set)');
        mod_spk = NMMfit_logexp_spkNL(mod,yy_obs(uu),temp',[],[],[3]);
        
        filtK = mod_spk.mods(1).filtK;
        %             filtK = mod.mods(1).filtK;
        pfun = csape(knot_pts,filtK,xx(uu),'var');
        fit_spk_nl = mod_spk.spk_NL_params;
        %             fit_spk_nl = mod.spk_NL_params;
        smoothed_eye_obs_fun(:,ii,cc) = fit_spk_nl(3)*log(1+exp(fit_spk_nl(2)*(fit_spk_nl(1) + fnval(pfun,eye_xx))));
        
        beta = fit_spk_nl(2)*(fit_spk_nl(1)+fnval(pfun,xx(uu)));
        too_big = beta > 50;
        cur_eyepos_p = fit_spk_nl(3)*log(1+exp(beta));
        cur_eyepos_p(too_big) = fit_spk_nl(3)*beta(too_big);
        %             cur_eyepos_p(cur_eyepos_p < eps) = eps;
        np_prate(uu,ii,cc) = cur_eyepos_p;
        
        out_bounds = find(eye_xx > max(xx(uu)) | eye_xx < min(xx(uu)));
        smoothed_eye_obs_fun(out_bounds,ii,cc) = nan;
        
        if cnt == 0
            close all
        end
        cnt = 1;
        
        %             [xs,ord] = sort(xx(uu));
        %             plot(eye_xx,squeeze(smoothed_eye_obs_fun(:,ii,cc)),'r'); hold on
        %             % plot(eye_xx,smooth_10000(:,ii),'m')
        %             plot(xs,smooth(yy_obs(uu(ord)),8,'loess'),'b-')
        %             % plot(xs,yy_obs(uu(ord)),'k.-')
        % %             plot(xx(uu),cur_eyepos_p,'g.')
        %             [xs,ord] = sort(xx(uu_xv));
        %             plot(xs,smooth(yy_obs(uu_xv(ord)),8,'loess'),'m-','linewidth',1)
        % %             plot(xx(uu_tr),yy_obs(uu_tr),'co')
        % plot(eye_xx,eyefun_pred_out_int(:,ii),'g','linewidth',2);
        %             %         [xs,ord] = sort(xx(uu));
        %             %         plot(xs,yy_rate(uu(ord)),'b.-')
        %             xlim(use_eye_range);
        %             pause
        %             clf
        
    end
end

% %     if xv_frac > 0
% %                 
%         cur_vec_obs_psth = reshape(squeeze(full_psth(:,:,cc)),[],1);
%         udata = find(ismember(reshape(TT,[],1),xv_trials));
%         vec_np_prate = reshape(squeeze(np_prate(:,:,cc)),[],1);
%         
% %         min_prate = n_knots/length(tr_trials);
% %         vec_np_prate(vec_np_prate < min_prate) = min_prate;
%         
% %         poss_eps = [0 0.01 0.05 0.1 0.15 0.2 0.3 0.4 0.5 1];
% %         for ee = 1:length(poss_eps)
% %             cur_nprate = vec_np_prate;
% %             cur_nprate(cur_nprate < poss_eps(ee)) = poss_eps(ee);
% %         espline_xvLL(ee) = nansum(cur_vec_obs_psth(udata).*log(cur_nprate(udata)) - cur_nprate(udata,:))./nansum(cur_vec_obs_psth(udata));
% %         end
%         spline_xvLL(ll,cc) = nansum(cur_vec_obs_psth(udata).*log(vec_np_prate(udata)) - vec_np_prate(udata,:))./nansum(cur_vec_obs_psth(udata));
% %     end
% %     
% %     opts = optimset('algorithm','active-set');
% %     [X,FVAL] = fmincon(@(x) offset_LL(x,vec_np_prate(udata),cur_vec_obs_psth(udata)),0.001,[],[],[],[],0,max(vec_np_prate(udata)),[],opts);
%      [X,FVAL] = fminbnd(@(x) offset_LL(x,vec_np_prate(udata),cur_vec_obs_psth(udata)),0,max(vec_np_prate));
%     spline_xvLL_min(ll,cc) = -FVAL/nansum(cur_vec_obs_psth(udata));
%      
%     np_psth_var(ll,cc) = nanvar(nanmean(np_prate(:,:,cc)));
%     np_em_var(ll,cc) = nanmean(nanvar(np_prate(:,:,cc)));
%      
% % opts.Display = 'off';
% % opts.GradObj = 'off';
% % opts.Algorithm = 'active-set';
% % initial_params = [100 100];
% % LB = [-Inf 0]; UB = [Inf Inf];
% % [fit_params,fval] = fmincon(@(K) spkNL_internal_LL(K,vec_np_prate(udata),cur_vec_obs_psth(udata)), initial_params,[],[],[],[],LB,UB,[],opts);
% 
% % kp = glmfit(vec_np_prate,cur_vec_obs_psth,'poisson');
% % glm_pvals = glmval(kp,vec_np_prate,'log');
% % glm_xvLL = nansum(cur_vec_obs_psth(udata).*log(glm_pvals(udata)) - glm_pvals(udata,:))./nansum(cur_vec_obs_psth(udata));
% 
% % end
% % end
% 
%%
% for cc = 1:length(tr_SUs);
% 
% n_bins = 20;
% all_xx = squeeze(neur_eye_pos(:,:,cc))';
% all_yy = squeeze(full_psth(:,:,cc));
% nnan = sum(~isnan(all_xx) & ~isnan(all_yy));
% nperbin = mean(nnan)/n_bins;
% ov_avg = nanmean(all_yy(:));
% prior_k = ov_avg*nperbin;
% prior_theta = 1/nperbin;
% post_theta = prior_theta./(nperbin*prior_theta + 1);
% 
% post_lambdas = nan(n_bins,length(rpt_taxis));
% for ii = 1:length(rpt_taxis)
% 
% xx = all_xx(:,ii);
% % xx = xx(randperm(length(xx)));
% xx(xx > use_eye_range(2) | xx < use_eye_range(1)) = nan;
% yy_obs = all_yy(:,ii);
% uu = all_trials(~isnan(xx(all_trials)) & ~isnan(yy_obs(all_trials)));
% 
% bin_edges = prctile(xx,linspace(0,100,n_bins+1));
% [n,ni] = histc(xx,bin_edges);
% n_spks = nan(n_bins,1);
% for bb = 1:n_bins
%     n_spks(bb) = sum(yy_obs(ni==bb));
% end
% post_k = prior_k + n_spks;
% 
% % post_lambdas(:,ii) = (post_k-1)*post_theta;
% post_lambdas(:,ii) = (post_k)*post_theta;
% end
% cur_EM_var(cc) = nanmean(nanvar(post_lambdas));
% cur_psth_var(cc) = nanvar(nanmean(post_lambdas));
% end
%% estimate model-rate as a function of time and eye pos.
% %generate shift matrices. Must be applied to the stimulus (not the filters)
% It = speye(flen);
% shift_mat = cell(n_shifts,1);
% for xx = 1:n_shifts
%     temp = spdiags( ones(full_nPix_us,1), -shifts(xx), full_nPix_us, full_nPix_us);
%     temp = kron(temp,It);
%     shift_mat{xx} = temp;
% end
% 
% cc = 21;
% srange = find(shifts*sp_dx >= use_eye_range(1) & shifts*sp_dx <= use_eye_range(2));
% upts = sum(~isnan(inf_eye_pos_t),2);
% [~,ut] = max(upts);
% test_rpt_inds = find(all_trialvec(used_inds(all_used_rpt_inds)) == rpt_trials(ut));
% eyefun_pred_out = nan(length(srange),length(rpt_taxis));
% uncor_X = all_Xmat_us(used_inds(all_used_rpt_inds),:);
% for ss = 1:length(srange)
%     fprintf('%d of %d\n',ss,length(srange));
%     shift_data = uncor_X*shift_mat{srange(ss)}';
%     cur_cormod = var_data.cor_mods(cc);
%     [~, ~, temp_prate] = NMMmodel_eval(cur_cormod,[],shift_data(test_rpt_inds,use_kInds_up));
%     eyefun_pred_out(ss,:) = temp_prate;
% end
% 
% eyefun_pred_out_int = interp1(-shifts(srange)*sp_dx,eyefun_pred_out,eye_xx');
% 
%%
for cc = 1:length(tr_SUs)
n_nearest = 40;
all_xx = squeeze(neur_eye_pos(:,:,cc))';
all_yy = squeeze(full_psth(:,:,cc));
% all_yy = squeeze(tr_sub_psth(:,:,cc));

clear cur_var cur_mean cur_rvar cur_rmean
[all_xx_sort,all_xx_sort_inds] = sort(abs(all_xx));
for ii = 1:length(rpt_taxis)
    uset = all_xx_sort_inds(1:n_nearest,ii);
    uset(isnan(all_xx(uset,ii))) = [];
    if length(uset) < n_nearest
        cur_var(ii) = nan;
        cur_mean(ii) = nan;
    else
        cur_var(ii) = var(all_yy(uset,ii));
        cur_mean(ii) = mean(all_yy(uset,ii));
    end
    
    poss = find(~isnan(all_xx(:,ii)));
    if length(poss) > n_nearest
        rperm = randperm(length(poss));
        uset = poss(rperm(1:n_nearest));
        cur_rvar(ii) = var(all_yy(uset,ii));
        cur_rmean(ii) = mean(all_yy(uset,ii));
    else
        cur_rvar(ii) = nan;
        cur_rmean(ii) = nan;
    end
end
EP_weight_FF(cc) = nansum(cur_var)./nansum(cur_mean);
PSTH_weight_FF(cc) = nansum(cur_rvar)./nansum(cur_rmean);
EP_weight_FF2(cc) = nanmean(cur_var./cur_mean);
PSTH_weight_FF2(cc) = nanmean(cur_rvar./cur_rmean);
tot_FF(cc) = nanvar(all_yy(:))/nanmean(all_yy(:));
end
%% %vectorize trial-to-trial firing rates

vec_obs_psth = reshape(permute(tr_sub_psth,[2 1 3]),[],length(tr_SUs));
% vec_spline_pred_rate = reshape(permute(np_prate,[2 1 3]),[],length(tr_SUs));
vec_psth_mod_rate = reshape(tbt_psth_mod_rate,[],length(tr_SUs));
vec_cor_pred_rate = reshape(tbt_cor_pred_rate,[],length(tr_SUs));
vec_uncor_pred_rate = reshape(tbt_uncor_pred_rate,[],length(tr_SUs));
vec_corsac_pred_rate = reshape(tbt_corsac_pred_rate,[],length(tr_SUs));
vec_sac_pred_rate = reshape(tbt_sac_pred_rate,[],length(tr_SUs));
% vec_sc_mod_pred = reshape(tbt_sc_mod_pred,[],length(tr_SUs));
%% Compute model predictions of trial-to-trial variability
tbt_rate_ms = bsxfun(@minus,tbt_cor_pred_rate,reshape(nanmean(vec_cor_pred_rate),[1 1 length(tr_SUs)]));
tbt_pred_diff = bsxfun(@minus,tbt_rate_ms,nanmean(tbt_rate_ms,2));

mod_stimlock_var = squeeze(nanvar(nanmean(tbt_rate_ms,2)));
mod_tot_var = nanvar(reshape(tbt_rate_ms,[],length(tr_SUs)));
mod_EM_var = nanvar(reshape(tbt_pred_diff,[],length(tr_SUs)));
mod_EM_var_frac = mod_EM_var./mod_tot_var;

var_data.mod_EM_var_frac = mod_EM_var_frac;
var_data.mod_tot_var = mod_tot_var;
var_data.mod_stimlock_var = mod_stimlock_var;

% %for spline model
% var_data.spline_em_var = squeeze(nanmean(nanvar(np_prate)));
% var_data.spline_psth_var = squeeze(nanvar(nanmean(np_prate)));
% var_data.spline_tot_var = nanvar(vec_spline_pred_rate);
% var_data.spline_EM_var_frac = var_data.spline_em_var./var_data.spline_tot_var';

%% create eye-sorted shift-predictor matrices
[II,JJ] = meshgrid(1:n_rpts);
Tinds = [II(:) JJ(:)];

use_neur_eye_pos = tot_neur_eye_pos;
% use_neur_eye_pos = squeeze(neur_eye_pos(:,:,26));

[a,b] = sort(use_neur_eye_pos,2);
% b(isnan(a)) = nan;

rand_eyepos = randn(size(use_neur_eye_pos));
rand_eyepos(isnan(use_neur_eye_pos)) = nan;
[arand,brand] = sort(rand_eyepos,2);

not_usable = isnan(use_neur_eye_pos');

b = b';a = a';
arand = arand'; brand = brand';
btshuff = b;

tr_sub_psth_ec = tr_sub_psth;
for ii = 1:length(tr_SUs)
    temp = squeeze(tr_sub_psth_ec(:,:,ii));
    temp(not_usable) = nan;
%     temp = temp - nanmean(temp(:)); %subtract off mean spike counts
    tr_sub_psth_ec(:,:,ii) = temp;
end

max_adiff = sp_dx*2;
adiff = [zeros(1,length(rpt_taxis)); diff(a)];

shift_size = 1;

full_psth_eyeperm1 = tr_sub_psth_ec;
full_psth_eyeperm2 = tr_sub_psth_ec;

% full_psth_rperm1 = tr_sub_psth_ec;
% full_psth_rperm2 = tr_sub_psth_ec;
for ii = 1:length(rpt_taxis)
    temp = squeeze(tr_sub_psth_ec(:,ii,:));
    tempr = temp(b(:,ii),:);
    full_psth_eyeperm1(:,ii,:) = tempr;

    temp = squeeze(tr_sub_psth_ec(:,ii,:));
    tempr(1:end-shift_size,:) = temp(b((1+shift_size):end,ii),:);
    tempr((end-shift_size+1):end,:) = nan;
    tempr(adiff(2:end,ii) > max_adiff,:) = nan;
    full_psth_eyeperm2(:,ii,:) = tempr;

%     temp = squeeze(tr_sub_psth_ec(:,ii,:));
%     tempr = temp(brand(:,ii),:);
%     full_psth_rperm1(:,ii,:) = tempr;
% 
%     temp = squeeze(tr_sub_psth_ec(:,ii,:));
%     tempr(1:end-1,:) = temp(brand(2:end,ii),:);
%     tempr(adiff(2:end,ii) > max_adiff,:) = nan;
%     full_psth_rperm2(:,ii,:) = tempr;
end

vec_obs_eyeperm1 = reshape(permute(full_psth_eyeperm1,[2 1 3]),[],length(tr_SUs));
vec_obs_eyeperm2 = reshape(permute(full_psth_eyeperm2,[2 1 3]),[],length(tr_SUs));

%%
% use_avgs = squeeze(nanmean(full_psth_rperm1));
% use_psth_var = nanvar(use_avgs);
% use_ntrials = squeeze(nanmean(sum(~isnan(full_psth_rperm1))));
% use_noise_var = squeeze(nanmean(nanvar(full_psth_rperm1)))./use_ntrials;
% var_data.np_cor_psth_var = use_psth_var - use_noise_var';

% avg_eyeperm = nanmean(full_psth_eyeperm1);
% full_psth_eyeperm1 = bsxfun(@minus,full_psth_eyeperm1,avg_eyeperm);
% avg_eyeperm = nanmean(full_psth_eyeperm2);
% full_psth_eyeperm2 = bsxfun(@minus,full_psth_eyeperm2,avg_eyeperm);


% avg_rperm = nanmean(full_psth_rperm1);
% full_psth_rperm1 = bsxfun(@minus,full_psth_rperm1,avg_rperm);
% avg_rperm = nanmean(full_psth_rperm2);
% full_psth_rperm2 = bsxfun(@minus,full_psth_rperm2,avg_rperm);

% vec_obs_rperm1 = reshape(permute(full_psth_rperm1,[2 1 3]),[],length(tr_SUs));
% vec_obs_rperm2 = reshape(permute(full_psth_rperm2,[2 1 3]),[],length(tr_SUs));

%%
var_data.shift_em_var = nanmean(vec_obs_eyeperm1.*vec_obs_eyeperm2) - nanmean(vec_obs_eyeperm1).*nanmean(vec_obs_eyeperm2);
var_data.psth_var = nanvar(vec_psth_mod_rate);
n_utrials = squeeze(mean(sum(~isnan(tr_sub_psth_ec))));
var_data.psth_noise_var = squeeze(nanmean(nanvar(tr_sub_psth_ec,[],2)))./n_utrials;
var_data.psth_var_cor = var_data.psth_var.*(n_utrials'./(n_utrials'-1)) - var_data.psth_noise_var';

var_data.shift_em_varfrac = 1-var_data.psth_var_cor./var_data.shift_em_var;

var_data.tot_var = nanvar(vec_obs_psth);
var_data.tot_mean = nanmean(reshape(full_psth,[],length(tr_SUs)));
var_data.tot_var2 = nanvar(reshape(full_psth,[],length(tr_SUs)));
corr_sc = 1./(1-mod_EM_var_frac);

var_data.probe_nums = tr_set;

%%

maxlag = round(0.1/dt);
lags = (-maxlag:maxlag);
Cmat_obs = nan(length(tr_SUs),length(tr_SUs),2*maxlag+1);
Cmat_mod = Cmat_obs;
% Cmat_spline = Cmat_obs;
Cmat_sac = Cmat_obs;
Cmat_eye = Cmat_obs;
Cmat_psth = Cmat_obs;
Cmat_psthr = Cmat_obs;
for ii = 1:length(tr_SUs)
%     ii
    for jj = 1:length(tr_SUs)
        if jj > ii
            uset = ~isnan(vec_obs_psth(:,ii)) & ~isnan(vec_obs_psth(:,jj));
%             usets = ~isnan(vec_spline_pred_rate(:,ii)) & ~isnan(vec_spline_pred_rate(:,jj));
            uset2 = ~isnan(vec_obs_eyeperm1(:,ii)) & ~isnan(vec_obs_eyeperm2(:,jj));
%             uset2r = ~isnan(vec_obs_rperm1(:,ii)) & ~isnan(vec_obs_rperm2(:,jj));
            if sum(uset) >0
                Cmat_obs(ii,jj,:) = xcov(vec_obs_psth(uset,ii),vec_obs_psth(uset,jj),maxlag,'unbiased');
                Cmat_mod(ii,jj,:) = xcov(vec_cor_pred_rate(uset,ii),vec_cor_pred_rate(uset,jj),maxlag,'unbiased');
                Cmat_sac(ii,jj,:) = xcov(vec_sac_pred_rate(uset,ii),vec_sac_pred_rate(uset,jj),maxlag,'unbiased');
                Cmat_psth(ii,jj,:) = xcov(vec_psth_mod_rate(uset,ii),vec_psth_mod_rate(uset,jj),maxlag,'unbiased');
            end
%             if sum(usets) > 0
%                 Cmat_spline(ii,jj,:) = xcov(vec_spline_pred_rate(usets,ii),vec_spline_pred_rate(usets,jj),maxlag,'unbiased');
%             end
            if sum(uset2) > 0
                Cmat_eye(ii,jj,:) = xcov(vec_obs_eyeperm1(uset2,ii),vec_obs_eyeperm2(uset2,jj),maxlag,'unbiased');
            else
                Cmat_eye(ii,jj,:) = nan(2*maxlag+1,1);
            end            
%             if sum(uset2r) > 0
%                 Cmat_psthr(ii,jj,:) = xcov(vec_obs_rperm1(uset2r,ii),vec_obs_rperm2(uset2r,jj),maxlag,'unbiased');
%             else
%                 Cmat_psthr(ii,jj,:) = nan(2*maxlag+1,1);
%             end
        elseif jj < ii
            Cmat_obs(ii,jj,:) = Cmat_obs(jj,ii,:);
            Cmat_mod(ii,jj,:) = Cmat_mod(jj,ii,:);
%             Cmat_spline(ii,jj,:) = Cmat_spline(jj,ii,:);
            Cmat_sac(ii,jj,:) = Cmat_sac(jj,ii,:);
            Cmat_eye(ii,jj,:) = Cmat_eye(jj,ii,:);
            Cmat_psth(ii,jj,:) = Cmat_psth(jj,ii,:);
            Cmat_psthr(ii,jj,:) = Cmat_psthr(jj,ii,:);
        end
    end
end

Cmat_eye_res = Cmat_obs - Cmat_eye;
% Cmat_spline_res = Cmat_obs - Cmat_spline;
Cmat_psth_res = Cmat_obs - Cmat_psth;
Cmat_mod_res = Cmat_obs - Cmat_mod;

% em_var_mat = bsxfun(@plus,cor_EM_var,cor_EM_var');
% tot_var_mat = bsxfun(@plus,cor_tot_var,cor_tot_var');
psth_var_frac = 1-mod_EM_var_frac;
% em_var_mat = bsxfun(@plus,cor_EM_var,cor_EM_var');
% tot_var_mat = bsxfun(@plus,cor_tot_var,cor_tot_var');
corr_sc_mat = bsxfun(@times,sqrt(psth_var_frac)',sqrt(psth_var_frac));
Cmat_psth_res2 = Cmat_obs - bsxfun(@rdivide,Cmat_psth,corr_sc_mat);

%%

% tot_var_stim = em_var + psth_var;
% tot_var_stim = em_var;
% tot_sc_mat_stim = sqrt(tot_var_stim'*tot_var_stim);
% var_data.EM_covar_frac = bsxfun(@rdivide,Cmat_eye,tot_sc_mat_stim);
% var_data.psth_covar_frac = bsxfun(@rdivide,Cmat_psth,tot_sc_mat_stim);


for ii = 1:length(tr_SUs)
    neighbs = [ii-1 ii+1];
    neighbs(neighbs > length(tr_SUs) | neighbs < 1) = [];
    Cmat_obs(ii,neighbs,maxlag+1) = nan;
%     Cmat_sac(ii,neighbs,maxlag+1) = nan;
    Cmat_eye_res(ii,neighbs,maxlag+1) = nan;
    Cmat_psth_res(ii,neighbs,maxlag+1) = nan;
    Cmat_mod_res(ii,neighbs,maxlag+1) = nan;
%     Cmat_spline_res(ii,neighbs,maxlag+1) = nan;
end

for ii = 1:length(tr_SUs)
    temp = squeeze(Cmat_obs(ii,:,:));
    temp2 = squeeze(Cmat_psth_res(ii,:,:));
    temp3 = squeeze(Cmat_eye_res(ii,:,:));
    temp4 = squeeze(Cmat_psth_res2(ii,:,:));
    if sum(~isnan(temp(:,maxlag+1))) > 3
        beta_psth(ii) = regress(temp2(:,maxlag+1),temp(:,maxlag+1));
        beta_psth2(ii) = regress(temp4(:,maxlag+1),temp(:,maxlag+1));
        beta_eye(ii) = regress(temp3(:,maxlag+1),temp(:,maxlag+1));
        beta_psth_sur(ii) = regress(reshape(temp2(:,[maxlag maxlag+2]),2*length(tr_SUs),[]),reshape(temp(:,[maxlag maxlag+2]),2*length(tr_SUs),[]));
        beta_psth_sur2(ii) = regress(reshape(temp4(:,[maxlag maxlag+2]),2*length(tr_SUs),[]),reshape(temp(:,[maxlag maxlag+2]),2*length(tr_SUs),[]));
        beta_eye_sur(ii) = regress(reshape(temp3(:,[maxlag maxlag+2]),2*length(tr_SUs),[]),reshape(temp(:,[maxlag maxlag+2]),2*length(tr_SUs),[]));
    else
        beta_psth(ii) = nan;
        beta_psth2(ii) = nan;
        beta_eye(ii) = nan;
        beta_psth_sur(ii) = nan;
        beta_psth_sur2(ii) = nan;
        beta_eye_sur(ii) = nan;
    end
    avg_cor_psth(ii) = nanmean(temp2(:,maxlag+1));
    avg_cor_psth2(ii) = nanmean(temp4(:,maxlag+1));
    avg_cor_eye(ii) = nanmean(temp3(:,maxlag+1));
    avg_cor_obs(ii) = nanmean(temp(:,maxlag+1));
    
    avg_cor_psth_sur(ii) = nanmean(reshape(temp2(:,[maxlag maxlag+2]),2*length(tr_SUs),[]));
    avg_cor_psth_sur2(ii) = nanmean(reshape(temp4(:,[maxlag maxlag+2]),2*length(tr_SUs),[]));
    avg_cor_eye_sur(ii) = nanmean(reshape(temp3(:,[maxlag maxlag+2]),2*length(tr_SUs),[]));
    avg_cor_obs_sur(ii) = nanmean(reshape(temp(:,[maxlag maxlag+2]),2*length(tr_SUs),[]));
end

tot_var_rates = nanvar(vec_obs_psth);
tot_sc_mat = sqrt(tot_var_rates'*tot_var_rates);

Cmat_obs = bsxfun(@rdivide,Cmat_obs,tot_sc_mat);
Cmat_mod = bsxfun(@rdivide,Cmat_mod,tot_sc_mat);
% Cmat_spline = bsxfun(@rdivide,Cmat_spline,tot_sc_mat);
Cmat_sac = bsxfun(@rdivide,Cmat_sac,tot_sc_mat);
Cmat_mod_res = bsxfun(@rdivide,Cmat_mod_res,tot_sc_mat);
% Cmat_spline_res = bsxfun(@rdivide,Cmat_spline_res,tot_sc_mat);
Cmat_eye = bsxfun(@rdivide,Cmat_eye,tot_sc_mat);
Cmat_eye_res = bsxfun(@rdivide,Cmat_eye_res,tot_sc_mat);
Cmat_psth_res = bsxfun(@rdivide,Cmat_psth_res,tot_sc_mat);
Cmat_psth_res2 = bsxfun(@rdivide,Cmat_psth_res2,tot_sc_mat);
Cmat_psth = bsxfun(@rdivide,Cmat_psth,tot_sc_mat);

var_data.cc_lags = lags*dt;
var_data.cc_maxlag = maxlag;
var_data.avg_cc_obs = avg_cor_obs;
var_data.avg_cc_eye = avg_cor_eye;
var_data.avg_cc_psth = avg_cor_psth;
var_data.avg_cc_psth2 = avg_cor_psth2;
var_data.avg_cc_obs_sur = avg_cor_obs_sur;
var_data.avg_cc_eye_sur = avg_cor_eye_sur;
var_data.avg_cc_psth_sur = avg_cor_psth_sur;
var_data.avg_cc_psth_sur2 = avg_cor_psth_sur2;
var_data.beta_cc_psth = beta_psth;
var_data.beta_cc_psth2 = beta_psth2;
var_data.beta_cc_eye = beta_eye;
var_data.beta_cc_psth_sur = beta_psth_sur;
var_data.beta_cc_psth_sur2 = beta_psth_sur2;
var_data.beta_cc_eye_sur = beta_eye_sur;

var_data.Cmat_obs = Cmat_obs;
var_data.Cmat_mod = Cmat_mod;
var_data.Cmat_sac = Cmat_sac;
var_data.Cmat_eye = Cmat_eye;
var_data.Cmat_psth = Cmat_psth;
var_data.Cmat_mod_res = Cmat_mod_res;
var_data.Cmat_psth_res = Cmat_psth_res;
var_data.Cmat_psth_res2 = Cmat_psth_res2;
var_data.Cmat_eye_res = Cmat_eye_res;

%%
[RR,TT] = meshgrid(1:length(rpt_taxis),1:n_rpts);
utrials = 1:n_rpts;
% utrials = xv_trials;
% utrials = tr_trials;

udata = find(ismember(reshape(TT',[],1),utrials));

neur_avg_rates = nanmean(reshape(full_psth,[],length(tr_SUs))); %overall avg rate for each neuron
neur_tot_spikes = nansum(reshape(full_psth,[],length(tr_SUs)));

vec_obs_psth = reshape(permute(full_psth,[2 1 3]),[],length(tr_SUs));
% tr_sub_psth = bsxfun(@rdivide,full_psth,nanmean(full_psth,2)); %divide out trial-by-trial diff in avg rate
% vec_obs_psth = reshape(permute(tr_sub_psth,[2 1 3]),[],length(tr_SUs));
% vec_obs_psth = bsxfun(@times,vec_obs_psth,neur_avg_rates); %multiply by overall avg rate

%output of psth model
ov_avg_rates = nanmean(full_psth(utrials,:,:));
% ov_avg_rates = nanmean(full_psth(utrials,:,:));
tbt_psth_mod_rate = permute(repmat(ov_avg_rates,[length(rpt_trials),1,1]),[2 1 3]);

% %vectorized model predictions
% vec_cor_pred_rate = reshape(tbt_cor_pred_rate,[],length(tr_SUs));
% vec_uncor_pred_rate = reshape(tbt_uncor_pred_rate,[],length(tr_SUs));
vec_psth_mod_rate = reshape(tbt_psth_mod_rate,[],length(tr_SUs));
vec_null_prate = repmat(neur_avg_rates,size(vec_obs_psth,1),1);

% vec_cor_pred_rate = bsxfun(@rdivide,vec_cor_pred_rate,nanmean(vec_cor_pred_rate));
% vec_uncor_pred_rate = bsxfun(@rdivide,vec_uncor_pred_rate,nanmean(vec_uncor_pred_rate));
% vec_psth_mod_rate = bsxfun(@rdivide,vec_psth_mod_rate,nanmean(vec_psth_mod_rate));

full_LL = nansum(vec_obs_psth(udata,:).*log(vec_obs_psth(udata,:)) - vec_obs_psth(udata,:))./nansum(vec_obs_psth(udata,:));
cormod_LL = nansum(vec_obs_psth(udata,:).*log(vec_cor_pred_rate(udata,:)) - vec_cor_pred_rate(udata,:))./nansum(vec_obs_psth(udata,:));
uncormod_LL = nansum(vec_obs_psth(udata,:).*log(vec_uncor_pred_rate(udata,:)) - vec_uncor_pred_rate(udata,:))./nansum(vec_obs_psth(udata,:));
psthmod_LL = nansum(vec_obs_psth(udata,:).*log(vec_psth_mod_rate(udata,:)) - vec_psth_mod_rate(udata,:))./nansum(vec_obs_psth(udata,:));
null_LL = nansum(vec_obs_psth(udata,:).*log(vec_null_prate(udata,:)) - vec_null_prate(udata,:))./nansum(vec_obs_psth(udata,:));

% %compute all LLs
% for ii = 1:length(tr_SUs)
%     u2 = udata(~isnan(vec_spline_pred_rate(udata,ii)));
%  
%     full_mLL(ii) = nansum(vec_obs_psth(u2,ii).*log(vec_obs_psth(u2,ii)) - vec_obs_psth(u2,ii))./nansum(vec_obs_psth(u2,ii));
%     cormod_mLL(ii) = nansum(vec_obs_psth(u2,ii).*log(vec_cor_pred_rate(u2,ii)) - vec_cor_pred_rate(u2,ii))./nansum(vec_obs_psth(u2,ii));
%     uncormod_mLL(ii) = nansum(vec_obs_psth(u2,ii).*log(vec_uncor_pred_rate(u2,ii)) - vec_uncor_pred_rate(u2,ii))./nansum(vec_obs_psth(u2,ii));
%     psthmod_mLL(ii) = nansum(vec_obs_psth(u2,ii).*log(vec_psth_mod_rate(u2,ii)) - vec_psth_mod_rate(u2,ii))./nansum(vec_obs_psth(u2,ii));
% %     spline_mLL(ii) = nansum(vec_obs_psth(u2,ii).*log(vec_spline_pred_rate(u2,ii)) - vec_spline_pred_rate(u2,ii))./nansum(vec_obs_psth(u2,ii));
%     null_mLL(ii) = nansum(vec_obs_psth(u2,ii).*log(vec_null_prate(u2,ii)) - vec_null_prate(u2,ii))./nansum(vec_obs_psth(u2,ii));
% end

%comptue all R1
uncor_r2 = 1-(full_LL-uncormod_LL)./(full_LL-null_LL);
cor_r2 = 1-(full_LL-cormod_LL)./(full_LL-null_LL);
psth_r2 = 1-(full_LL-psthmod_LL)./(full_LL-null_LL);
% spline_r2 = 1-(full_LL - spline_LL)./(full_LL-null_LL);
% uncor_mr2 = 1-(full_mLL-uncormod_mLL)./(full_mLL-null_mLL);
% cor_mr2 = 1-(full_mLL-cormod_mLL)./(full_mLL-null_mLL);
% psth_mr2 = 1-(full_mLL-psthmod_mLL)./(full_mLL-null_mLL);
% % spline_mr2 = 1-(full_mLL - spline_mLL)./(full_mLL-null_mLL);

var_data.avg_rates = neur_avg_rates;
var_data.tot_spikes = neur_tot_spikes;
var_data.uncor_r2 = uncor_r2;
var_data.cor_r2 = cor_r2;
var_data.psth_r2 = psth_r2;
% var_data.spline_r2 = spline_r2;
var_data.SU_nums = tr_set(tr_SUs);
var_data.n_rpts = n_rpts;
%%
var_anal_name = 'variability_anal5';
var_anal_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
cd(var_anal_dir);
save(var_anal_name,'var_data');

%%
% su_pt = find(tr_set(tr_SUs) > n_probes,1);
% % close all
% h = figure();
% for ii = 1:length(tr_SUs)
% 
% % ii = 3;
% 
%    subplot(2,4,1)
%    imagescnan(lags*dt,1:length(tr_SUs),squeeze(Cmat_obs(ii,:,:)));
%    ca = caxis(); cam = 0.75*max(abs(ca)); caxis([-cam cam]);
%    yl = ylim();
%    line(lags([maxlag+1 maxlag+1])*dt,yl,'color','w');
%    line(lags([1 end])*dt,[su_pt su_pt]-0.5)
%    title('Measured');
%      xlabel('Time lag (s)');
%    ylabel('Unit');
%   
%    subplot(2,4,2)
%    imagescnan(lags*dt,1:length(tr_SUs),squeeze(Cmat_psth_res(ii,:,:)));
%    caxis([-cam cam]);
%    line(lags([maxlag+1 maxlag+1])*dt,yl,'color','w');
%    line(lags([1 end])*dt,[su_pt su_pt]-0.5)
%    title('PSTH-corrected');
%     xlabel('Time lag (s)');
%    ylabel('Unit');
%     
% 
%    subplot(2,4,3)
%    imagescnan(lags*dt,1:length(tr_SUs),squeeze(Cmat_eye_res(ii,:,:)));
%    caxis([-cam cam]);
%    line(lags([maxlag+1 maxlag+1])*dt,yl,'color','w');
%    line(lags([1 end])*dt,[su_pt su_pt]-0.5)
%    title('Eyepos-corrected');
%     xlabel('Time lag (s)');
%    ylabel('Unit');
% 
%       subplot(2,4,4)
% %    imagescnan(lags*dt,1:length(tr_SUs),squeeze(Cmat_mod_res(ii,:,:)));
%    imagescnan(lags*dt,1:length(tr_SUs),squeeze(Cmat_psth_res2(ii,:,:)));
%    line(lags([maxlag+1 maxlag+1])*dt,yl,'color','w');
%    line(lags([1 end])*dt,[su_pt su_pt]-0.5)
%    caxis([-cam cam]);
%    title('EM-PSTH-corrected');
%     xlabel('Time lag (s)');
%    ylabel('Unit');
% 
%    subplot(2,4,5)
%    imagescnan(lags*dt,1:length(tr_SUs),squeeze(Cmat_spline(ii,:,:)));
% %    ca = caxis(); cam = 0.75*max(abs(ca)); caxis([-cam cam]);
%    caxis([-cam cam]);
%    line(lags([maxlag+1 maxlag+1])*dt,yl,'color','w');
%    line(lags([1 end])*dt,[su_pt su_pt]-0.5)
%    title('Spline-pred');
%     xlabel('Time lag (s)');
%    ylabel('Unit');
% 
%    subplot(2,4,6)
%    imagescnan(lags*dt,1:length(tr_SUs),squeeze(Cmat_psth(ii,:,:)));
% %    ca = caxis(); cam = 0.75*max(abs(ca)); caxis([-cam cam]);
%    caxis([-cam cam]);
%    line(lags([maxlag+1 maxlag+1])*dt,yl,'color','w');
%    line(lags([1 end])*dt,[su_pt su_pt]-0.5)
%    title('PSTH-pred');
%     xlabel('Time lag (s)');
%    ylabel('Unit');
% 
%    subplot(2,4,7)
%    imagescnan(lags*dt,1:length(tr_SUs),squeeze(Cmat_mod(ii,:,:)));
% %    ca = caxis(); cam = 0.75*max(abs(ca)); caxis([-cam cam]);
%    caxis([-cam cam]);
%    line(lags([maxlag+1 maxlag+1])*dt,yl,'color','w');
%    line(lags([1 end])*dt,[su_pt su_pt]-0.5)
%    title('Model-pred');
%     xlabel('Time lag (s)');
%    ylabel('Unit');
% 
%    
%    subplot(2,4,8)
%    imagescnan(lags*dt,1:length(tr_SUs),squeeze(Cmat_spline_res(ii,:,:)));
% %    ca = caxis(); cam = 0.75*max(abs(ca)); caxis([-cam cam]);
%    caxis([-cam cam]);
%    line(lags([maxlag+1 maxlag+1])*dt,yl,'color','w');
%    line(lags([1 end])*dt,[su_pt su_pt]-0.5)
%    title('Spline-res');
%     xlabel('Time lag (s)');
%    ylabel('Unit');
% 
%    pause
%    clf
% 
% % out_dir = '/home/james/Analysis/bruce/variability/';
% % fname = sprintf('examp_xcorr_C%d.pdf',ii);
% % fig_width = 10;
% % rel_height = 0.5;
% % figufy(h);
% % exportfig(h,[out_dir fname],'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% 
% 
% end
% 
%% PRINT RATE VS EYE POS PLOT
% % usm = 1;
% 
% out_dir = '/home/james/Analysis/bruce/variability/';
% fname = sprintf('examp_eye_rate_map_C%d.pdf',cc);
% xl = [1 2];
% ca = [0.2 0.9];
% 
% h = figure();
% subplot(3,1,1)
% imagescnan(rpt_taxis,eye_xx,squeeze(smoothed_eye_obs_fun(:,:,cc)));
% xlim(xl);
% xlabel('Time (s)');
% ylabel('Eye position (deg)');
% caxis(ca);
% 
% subplot(3,1,2)
% imagescnan(rpt_taxis,eye_xx,eyefun_pred_out_int);
% xlim(xl);
% caxis(ca);
% xlabel('Time (s)');
% ylabel('Eye position (deg)');
% 
% subplot(3,1,3)
% plot(rpt_taxis,nanmean(full_psth(:,:,cc)));
% hold on
% plot(rpt_taxis,nanmean(tbt_cor_pred_rate(:,:,cc),2),'r')
% xlabel('Time (s)');
% ylabel('Rate (Hz)');
% xlim(xl);
% ylim(ca)
% 
% fig_width = 3;
% rel_height = 2.1;
% figufy(h);
% exportfig(h,[out_dir fname],'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%% SORTED VS UNSORTED PSTHS
% close all
% cc = 31;
% out_dir = '/home/james/Analysis/bruce/variability/';
% fname = sprintf('sorted_PSTH_C%d.pdf',cc);
% xl = [1. 3];
% ca = [0 2.5];
% 
% h = figure();
% subplot(2,1,1)
% imagescnan(rpt_taxis,1:n_rpts,squeeze(full_psth(:,:,cc)));
% xlim(xl);
% xlabel('Time (s)');
% ylabel('Trial');
% caxis(ca);
% 
% subplot(2,1,2)
% imagescnan(rpt_taxis,1:n_rpts,squeeze(full_psth_eyeperm1(:,:,cc)));
% xlim(xl);
% caxis(ca);
% xlabel('Time (s)');
% ylabel('Trial');
% caxis(ca);
% 
% fig_width = 3;
% rel_height = 1.5;
% 
% figufy(h);
% exportfig(h,[out_dir fname],'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

%%
% xl = [-0.025 0.2];
% out_dir = '/home/james/Analysis/bruce/variability/';
% fname = 'Correlation_scatter.pdf';
% 
% h = figure();
% subplot(2,1,1)
% plot(Cmat_obs(:,:,maxlag+1),Cmat_psth(:,:,maxlag+1),'r.')
% hold on
% plot(Cmat_obs(:,:,maxlag+1),Cmat_psth(:,:,maxlag+1)./corr_sc_mat,'k.')
% line(xl,xl);
% xlim(xl);ylim(xl);
% xlabel('Observed correlation');
% ylabel('PSTH correlation');
% 
% xl = [-0.015 0.1];
% subplot(2,1,2)
% plot(mean(Cmat_obs(:,:,[maxlag maxlag+2]),3),mean(Cmat_psth(:,:,[maxlag maxlag+2]),3),'r.')
% hold on
% plot(mean(Cmat_obs(:,:,[maxlag maxlag+2]),3),mean(Cmat_psth(:,:,[maxlag maxlag+2]),3)./corr_sc_mat,'k.')
% xlabel('Observed correlation');
% ylabel('PSTH correlation');
% line(xl,xl);
% xlim(xl);ylim(xl);
% 
% fig_width = 4;
% rel_height = 1.5;
% 
% figufy(h);
% exportfig(h,[out_dir fname],'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% 
%%
% cc = 26;
% 
% out_dir = '/home/james/Analysis/bruce/variability/';
% fname = sprintf('examp_raster_E%d_C%d.pdf',Expt_num,cc);
% 
% close all
% xr = [0.5 1.5];
% % xr = [0.8 1.8];
% max_t = 45;
% 
% ep = find(rpt_taxis > xr(2),1);
% utrials = find(~any(isnan(full_psth(:,1:ep,cc)),2)); %only use trials that are complete up to xr(2)
% n_utrials = length(utrials);
% 
% lw = 0.5;
% theight = 0.75;
% 
% f1 = figure(); 
% subplot(3,1,[1 2]); hold on
% for nn = 1:length(utrials)
%     cur_spikes = find(all_spk_trials{cc} == utrials(nn));
%     for ii = 1:length(cur_spikes)
%         line(all_rel_times{cc}(cur_spikes([ii ii])),nn + [0 theight],'color','k','linewidth',lw);
%     end    
% end
% xlim(xr);
% ut = min(n_utrials,max_t);
% ylim([1 ut+1]);
% xlabel('Time (s)');
% ylabel('Trial');
% 
% rpt_ut = linspace(rpt_taxis(1),rpt_taxis(end),1e3);
% cur_psth = squeeze(nanmean(full_psth(utrials(1:ut),:,cc)))/dt;
% int_psth = interp1(rpt_taxis,cur_psth,rpt_ut,'pchip');
% tuse = find(rpt_taxis >= xr(1) & rpt_taxis <= xr(2));
% subplot(3,1,3);hold on
% % plot(rpt_taxis,squeeze(nanmean(full_psth(utrials(1:ut),:,cc)))/dt,'k')
% plot(rpt_ut,int_psth,'k')
% xlim(xr);
% xlabel('Time (s)');
% ylabel('Firing rate (Hz)');
% mm = max(cur_psth(tuse));
% ylim([0 mm+10]);
% 
% fig_width = 4;
% rel_height = 0.8;
% 
% % figufy(f1);
% % exportfig(f1,[out_dir fname],'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);

