%
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements//');
addpath('~/James_scripts/bruce/processing/');
addpath('~/James_scripts/bruce/saccade_modulation/');
addpath('~/James_scripts/TentBasis2D/');

global Expt_name bar_ori use_MUA

% Expt_name = 'M294';
Expt_name = 'M275';
use_MUA = false;
bar_ori = 90; %bar orientation to use (only for UA recs)

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
sac_buff = round(0.05/dt);
in_sac_inds = zeros(NT,1);
for ii = 1:sac_buff+1
    cur_inds = saccade_stop_inds + (ii-1);
    cur_inds(cur_inds > NT) = [];
    in_sac_inds(cur_inds) = 1;
end
for ii = 1:length(saccade_start_inds)
    cur_inds = saccade_start_inds(ii):saccade_stop_inds(ii);
    in_sac_inds(cur_inds) = 1;
end

in_sac_inds = logical(in_sac_inds);

blink_inds = find(used_is_blink);
in_blink_inds = zeros(NT,1);
for ii = 1:length(blink_inds)
    cur_inds = saccade_start_inds(blink_inds(ii)):saccade_stop_inds(blink_inds(ii));
    in_blink_inds(cur_inds) = 1;
end
in_blink_inds = logical(in_blink_inds);

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
all_Xmat = create_time_embedding(all_stimmat_up,stim_params_us);
all_Xmat = all_Xmat(used_inds(all_rpt_inds),use_kInds_up);

%% Recon retinal stim for non LOO data
cur_fix_post_mean = squeeze(it_fix_post_mean(end,:));
cur_fix_post_std = squeeze(it_fix_post_std(end,:));
cur_drift_post_mean = squeeze(drift_post_mean(end,:));
cur_drift_post_std = squeeze(drift_post_std(end,:));
[fin_tot_corr,fin_tot_std] = construct_eye_position(cur_fix_post_mean,cur_fix_post_std,...
    cur_drift_post_mean,cur_drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);

if length(used_inds) ~= length(cur_drift_post_mean)
    error('ET data mismatch!');
end

fin_shift_cor = round(fin_tot_corr);
fin_shift_cor(isnan(fin_shift_cor)) = 0;

%RECOMPUTE XMAT
best_shift_stimmat_up = all_stimmat_up;
for i=1:NT
    best_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
end

all_Xmat_shift = create_time_embedding(best_shift_stimmat_up,stim_params_us);
all_Xmat_shift = all_Xmat_shift(used_inds,use_kInds_up);

%%
all_EP = nan(n_rpts,length(rpt_taxis),length(targs));
full_psth = nan(n_rpts,length(rpt_taxis),length(targs));
mod_prates = nan(n_rpts,length(rpt_taxis),length(targs));
% uncor_prates = nan(n_rpts,length(rpt_taxis),length(targs));
all_fin_tot_corr = nan(length(all_rpt_inds),length(targs));
for ss = 1:length(targs)
    cc = targs(ss)
    cur_Robs = Robs_mat(:,cc);
    cc_uinds = find(~isnan(cur_Robs(all_rpt_inds)));
    if ~isempty(cc_uinds)
        %% RECONSTRUCT LOO STIM
        loo_cc = find(loo_set == cc); %index within the LOOXV set
        
%         cur_fix_post_mean = squeeze(it_fix_post_mean_LOO(loo_cc,end,:));
%         cur_fix_post_std = squeeze(it_fix_post_std_LOO(loo_cc,end,:));
%         cur_drift_post_mean = squeeze(drift_post_mean_LOO(loo_cc,end,:));
%         cur_drift_post_std = squeeze(drift_post_std_LOO(loo_cc,end,:));
%         [fin_tot_corr,fin_tot_std] = construct_eye_position(cur_fix_post_mean,cur_fix_post_std,...
%             cur_drift_post_mean,cur_drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);
%         
%         all_fin_tot_corr(:,ss) = fin_tot_corr(all_rpt_inds);
%         fin_shift_cor = round(fin_tot_corr);
%         
%         %RECOMPUTE XMAT
%         all_shift_stimmat_up = all_stimmat_up;
%         for i=1:NT
%             all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
%         end
%         all_Xmat_shift = create_time_embedding(all_shift_stimmat_up,stim_params_us);
%         all_Xmat_shift = all_Xmat_shift(used_inds(all_rpt_inds),use_kInds_up);
        
        %%
        cor_rGQM = ModData(cc).rectGQM;
        cor_rGQM = NMMfit_logexp_spkNL(cor_rGQM,cur_Robs(all_rpt_inds(cc_uinds)),all_Xmat_shift(all_rpt_inds(cc_uinds),:));
        [~,~,cor_prate] = NMMmodel_eval(cor_rGQM,cur_Robs(all_rpt_inds),all_Xmat_shift(all_rpt_inds,:));
        
        %     uncor_rGQM = ModData(cc).rectGQM_unCor;
        %     uncor_rGQM = NMMfit_logexp_spkNL(uncor_rGQM,cur_Robs(all_rpt_inds(cc_uinds)),all_Xmat(cc_uinds,:));
        %     [~,~,uncor_prate] = NMMmodel_eval(uncor_rGQM,cur_Robs(all_rpt_inds),all_Xmat);
        
        %%
        for ii = 1:length(rpt_trials)
            cur_inds = find(all_trialvec(used_inds) == rpt_trials(ii));
            cur_inds(size(full_psth,2)+1:end) = [];
            full_psth(ii,1:length(cur_inds),ss) = Robs_mat(cur_inds,cc);
            
            %         all_EP(ii,1:length(cur_inds),cc) = fin_tot_corr(cur_inds);
            
            cur_rinds = find(ismember(all_rpt_inds,cur_inds));
            mod_prates(ii,1:length(cur_rinds),ss) = cor_prate(cur_rinds);
            %         uncor_prates(ii,1:length(cur_rinds),ss) = uncor_prate(cur_rinds);
        end
    end
end
%%

back_look = 7;
back_kern = zeros(back_look*2+1,1);
back_kern(1:back_look+1) = 1;
back_kern = flipud(back_kern/sum(back_kern));

rpt_EP = fin_tot_corr(all_rpt_inds);
rpt_EPnan = rpt_EP;
rpt_EPnan(in_sac_inds(all_rpt_inds)) = nan;
rpt_EPnan(in_blink_inds(all_rpt_inds)) = nan;

sp = NMMcreate_stim_params(back_look);
rpt_EP_emb = create_time_embedding(rpt_EPnan(:),sp);


in_sac = false(n_rpts,length(rpt_taxis));
full_EP = nan(n_rpts,length(rpt_taxis));
full_EP_emb = nan(n_rpts,length(rpt_taxis),back_look);
full_EP_filt = nan(n_rpts,length(rpt_taxis));
for ii = 1:length(rpt_trials)
    cur_inds = find(all_trialvec(used_inds(all_rpt_inds)) == rpt_trials(ii));
    in_sac(ii,1:length(cur_inds)) = in_sac_inds(all_rpt_inds(cur_inds)) | in_blink_inds(all_rpt_inds(cur_inds));
    full_EP(ii,1:length(cur_inds)) = rpt_EP(cur_inds);
    full_EP_emb(ii,1:length(cur_inds),:) = rpt_EP_emb(cur_inds,:);
    
    full_EP_filt(ii,1:length(cur_inds)) = conv(rpt_EP(cur_inds),back_kern,'same');
end
full_EP(in_sac) = nan;
full_EP_filt(in_sac) = nan;


%% LOAD LFP DATA
cd(data_dir)
if strcmp(rec_type,'LP')
    Fs = 1000;
    dsf = 6;
    Fsd = Fs/dsf;
    niqf = Fs/2;
    new_niqf = Fsd/2;
    [bb,aa] = butter(2,0.8*new_niqf/niqf);
    use_lfps = 1:2:n_probes;
elseif strcmp(rec_type,'UA')
    Fs = 400.0032;
    dsf = 4;
    Fsd = Fs/dsf;
    niqf = Fs/2;
    new_niqf = Fsd/2;
    [bb,aa] = butter(2,0.8*new_niqf/niqf);
    use_lfps = 1:n_probes;
    %     use_lfps = [28 SU_probes];
end

full_lfps = [];
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
        full_lfp_taxis = cat(1,full_lfp_taxis,expt_lfp_t_axis(uinds));
    end
    cur_toffset = trial_toffset(ee);
end

interp_lfps = interp1(full_lfp_taxis,full_lfps,all_t_axis);
% interp_lfps(ismember(all_blockvec,bad_blocks),:) = nan;
interp_lfps = interp_lfps(used_inds,:);
interp_lfps = nanzscore(interp_lfps);

%%
tbt_LFPs = nan(length(rpt_trials),length(rpt_taxis),length(use_lfps));
for ii = 1:length(rpt_trials)
    cur_inds = find(all_trialvec(used_inds(all_rpt_inds)) == rpt_trials(ii));
    tbt_LFPs(ii,1:length(cur_inds),:) = interp_lfps(cur_inds,:);
end
[a,b] = sort(full_EP_filt);
sorted_LFPs = nan(size(tbt_LFPs));
for tt = 1:length(rpt_taxis)
    sorted_LFPs(:,tt,:) = tbt_LFPs(b(:,tt),tt,:);
    sorted_LFPs(isnan(a(:,tt)),tt,:) = nan;
end

%%
silent = 1;
for cc = 1:length(use_lfps);
    cur_LFP = interp_lfps(:,cc);
    Xmat = abs(all_Xmat_shift);
    
    fit_inds = find(ismember(all_trialvec(used_inds),use_trials));
    mod_signs = [1];
    Xtargets = [1];
    lfp_stim_params = NMMcreate_stim_params([flen use_nPix_us]);
    
    d2XT = 20;
    reg_params = NMMcreate_reg_params('lambda_d2XT',d2XT);
    lfp_mod = NMMinitialize_model(lfp_stim_params,mod_signs,{'lin'},reg_params,Xtargets,[],'linear');
    lfp_mod = NMMfit_filters(lfp_mod,cur_LFP(fit_inds),Xmat(fit_inds,:),[],[],silent);
    
    [~,~,predlfp] = NMMmodel_eval(lfp_mod,cur_LFP,Xmat);
    lfp_expvar(cc) = var(predlfp);
    
    stimpred_LFP(cc,:) = predlfp;
end
residual_LFPs = interp_lfps - stimpred_LFP';

interp_stimpred = interp1(all_t_axis(used_inds),stimpred_LFP',full_lfp_taxis);
interp_stimpred(isnan(interp_stimpred)) = 0;
full_residual_LFPs = full_lfps - interp_stimpred;
%%
maxl = round(1/dt);
all_LFP_TAs = nan(length(targs),maxl*2+1,length(use_lfps));
for ss = 1:length(targs)
    cc = targs(ss)
  cur_Robs = Robs_mat(:,cc);
  cc_uinds = find(~isnan(cur_Robs));
  
  spk_inds = convert_to_spikebins(cur_Robs(cc_uinds));
  [cur_TA,cur_lags] = get_event_trig_avg_v3(interp_lfps(cc_uinds,:),spk_inds,maxl,maxl);
  
  all_LFP_TAs(ss,:,:) = cur_TA;
    
end

%%
nwfreqs = 15;
min_freq = 1; max_freq = 40;
min_scale = 1/max_freq*Fsd;
max_scale = 1/min_freq*Fsd;
scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);

all_ampgram = nan(length(all_t_axis),length(use_lfps),length(wfreqs));
all_phasegram = nan(length(all_t_axis),length(use_lfps),length(wfreqs));
for cc = 1:length(use_lfps)
    cc
%     cur_LFP = full_lfps(:,cc);
    cur_LFP = full_residual_LFPs(:,cc);
    cur_cwt = cwt(cur_LFP,scales,'cmor1-1');
    cur_phasegram = angle(cur_cwt);
    cur_ampgram = abs(cur_cwt);
    interp_phasegram = mod(interp1(full_lfp_taxis,unwrap(cur_phasegram'+pi),all_t_axis),2*pi)-pi;
    interp_ampgram = interp1(full_lfp_taxis,cur_ampgram',all_t_axis);
    
    all_ampgram(:,cc,:) = interp_ampgram;
    all_phasegram(:,cc,:) = interp_phasegram;
end

all_ampgram = all_ampgram(used_inds,:,:);
all_phasegram = all_phasegram(used_inds,:,:);

X{1} = reshape(cos(all_phasegram).*all_ampgram,[],length(wfreqs)*length(use_lfps));
X{2} = reshape(sin(all_phasegram).*all_ampgram,[],length(wfreqs)*length(use_lfps));

%%
for ss = 1:length(targs)
    cc = targs(ss)
    
    if cc > n_probes
        su_ind = find(SU_numbers == all_mod_SUnum(cc));
        [~,nearest_lfp] = min(abs(use_lfps-SU_probes(su_ind)));
    else
        [~,nearest_lfp] = min(abs(use_lfps - cc));
    end
    cur_LFP = full_lfps(:,nearest_lfp);
    
    cur_Robs = Robs_mat(:,cc);
    cc_uinds = find(~isnan(cur_Robs) & ~isnan(interp_lfps(:,nearest_lfp)));
    cc_uinds(~ismember(all_trialvec(used_inds(cc_uinds)),use_trials)) = [];
    if ~isempty(cc_uinds)
    %%
    cor_rGQM = ModData(cc).rectGQM;
    cor_rGQM = NMMfit_logexp_spkNL(cor_rGQM,cur_Robs(cc_uinds),all_Xmat_shift(cc_uinds,:));
    [~,~,cor_prate,G] = NMMmodel_eval(cor_rGQM,cur_Robs,all_Xmat_shift);
    stim_G = G - cor_rGQM.spk_NL_params(1);
    X{3} = stim_G;
    %%
    
    %     nwfreqs = 15;
    %     min_freq = 1; max_freq = 40;
    %     min_scale = 1/max_freq*Fsd;
    %     max_scale = 1/min_freq*Fsd;
    %     scales = logspace(log10(min_scale),log10(max_scale),nwfreqs);
    %     wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
    %
    %     cur_cwt = cwt(cur_LFP,scales,'cmor1-1');
    %     cur_phasegram = angle(cur_cwt);
    %     cur_ampgram = abs(cur_cwt);
    %     interp_phasegram = mod(interp1(full_lfp_taxis,unwrap(cur_phasegram'),all_t_axis),2*pi);
    %     interp_ampgram = interp1(full_lfp_taxis,cur_ampgram',all_t_axis);
    %     interp_ampgram = bsxfun(@rdivide,interp_ampgram,nanstd(interp_ampgram));
    %
    %     AXcos = cos(interp_phasegram(used_inds(cc_uinds),:)).*interp_ampgram(used_inds(cc_uinds),:);
    %     AXsin = sin(interp_phasegram(used_inds(cc_uinds),:)).*interp_ampgram(used_inds(cc_uinds),:);
    %     Xcos = cos(interp_phasegram(used_inds(cc_uinds),:));
    %     Xsin = sin(interp_phasegram(used_inds(cc_uinds),:));
    
    %%
    
    silent = 1;
    NL_types = repmat({'lin'},1,3);
%     mod_signs = [1 1 1];
%     Xtargets = [1 2 3];
    mod_signs = [1 1];
    Xtargets = [1 2];
    sac_stim_params(1:2) = NMMcreate_stim_params([length(wfreqs) length(use_lfps)]);
    sac_stim_params(3) = NMMcreate_stim_params(1);
    %     tr_stim{1} = AXcos;
    %     tr_stim{2} = AXsin;
    %     tr_stim{1} = Xcos;
    %     tr_stim{2} = Xsin;
    
    d2T = 10;
    sac_reg_params = NMMcreate_reg_params('lambda_d2T',d2T);
    lfp_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
    if length(mod_signs) == 3
    lfp_mod.mods(3).reg_params = NMMcreate_reg_params();
    end
    lfp_mod = NMMfit_filters(lfp_mod,cur_Robs(cc_uinds),get_Xcell_tInds(X,cc_uinds),[],[],silent);
    
    %%
    %     AXcos = cos(interp_phasegram(used_inds,:)).*interp_ampgram(used_inds,:);
    %     AXsin = sin(interp_phasegram(used_inds,:)).*interp_ampgram(used_inds,:);
    %     tr_stim{1} = AXcos;
    %     tr_stim{2} = AXsin;
    [~,~,all_lfp_prate,lfp_G,~,lfp_fgint] = NMMmodel_eval(lfp_mod,[],X);
    lfp_G = sum(lfp_fgint(:,1:2),2);
    
    cc_rpt_inds = all_rpt_inds(~isnan(Robs_mat(all_rpt_inds,cc)));
    [rpt_LL,~,~,~,~,~,rpt_nullLL] = NMMmodel_eval(lfp_mod,Robs_mat(cc_rpt_inds,cc),get_Xcell_tInds(X,cc_rpt_inds));
    
    %%
    %     max_lag = 40;
    %     lag_vec = -max_lag:max_lag;
    %     tot_lags = length(lag_vec);
    %     cur_LFP = interp_lfps(cc_uinds,:);
    %     %otherwise loop over lags and manually shift the stim matrix
    %     Npix = 24;
    %     Xmat = zeros(length(cc_uinds), Npix*tot_lags);
    %     for n = 1:tot_lags
    %         Xmat(:,n-1+(1:tot_lags:(Npix*tot_lags))) = shift_mat_zpad( cur_LFP, lag_vec(n), 1);
    %     end
    %
    %     stim_params = NMMcreate_stim_params([tot_lags Npix]);
    %     reg_params = NMMcreate_reg_params('lambda_L2',5);
    %
    %       lfp_mod = NMMinitialize_model(stim_params,1,{'lin'},reg_params,1);
    %     lfp_mod = NMMfit_filters(lfp_mod,cur_Robs(cc_uinds),Xmat,[],[],silent);
    
    %%
    tbt_LFP_ratepred = nan(length(rpt_trials),length(rpt_taxis));
    tbt_stim_ratepred = nan(length(rpt_trials),length(rpt_taxis));
    tbt_LFP_G = nan(length(rpt_trials),length(rpt_taxis));
    for ii = 1:length(rpt_trials)
        cur_inds = find(all_trialvec(used_inds(all_rpt_inds)) == rpt_trials(ii));
        tbt_LFP_ratepred(ii,1:length(cur_inds)) = all_lfp_prate(all_rpt_inds(cur_inds));
        tbt_stim_ratepred(ii,1:length(cur_inds)) = cor_prate(all_rpt_inds(cur_inds));
        tbt_LFP_G(ii,1:length(cur_inds)) = lfp_G(all_rpt_inds(cur_inds));
    end
    % [a,b] = sort(full_EP_filt);
    [a,b] = sort(tbt_LFP_G);
    sorted_LFP_pred = nan(size(tbt_LFP_ratepred));
    sorted_LFP_psth = nan(size(tbt_LFP_ratepred));
    for tt = 1:length(rpt_taxis)
        sorted_LFP_pred(:,tt) = tbt_LFP_ratepred(b(:,tt),tt,:);
        sorted_LFP_pred(isnan(a(:,tt)),tt) = nan;
        sorted_LFP_psth(:,tt) = full_psth(b(:,tt),tt,find(targs==cc));
        sorted_LFP_psth(isnan(a(:,tt)),tt) = nan;
    end
    
    lfp_mod_var(ss) = nanvar(tbt_LFP_ratepred(:));
    lfp_mod_psth(ss,:) = nanmean(tbt_LFP_ratepred);
    %%
    end
end







%%
rpt_avg_rates = nanmean(reshape(full_psth,[],length(targs)));
rpt_avg_prates = nanmean(reshape(mod_prates,[],length(targs)));
full_psth_ms = bsxfun(@minus,full_psth,reshape(rpt_avg_rates,1,1,[]));
trial_avg_rates = nanmean(full_psth_ms,2);
full_psth_ms = bsxfun(@minus,full_psth_ms,trial_avg_rates);
mod_prates_ms = bsxfun(@minus,mod_prates,reshape(rpt_avg_prates,1,1,[]));
all_psths = squeeze(nanmean(full_psth_ms));

%%
ep_naninds = isnan(full_EP(:));
full_psth_ms = reshape(full_psth_ms,[],length(targs));
full_psth_ms(ep_naninds,:) = nan;
full_psth_ms = reshape(full_psth_ms,n_rpts,[],length(targs));

mod_prates_ms = reshape(mod_prates_ms,[],length(targs));
mod_prates_ms(ep_naninds,:) = nan;
mod_prates_ms = reshape(mod_prates_ms,n_rpts,[],length(targs));

%%
[a,b] = sort(full_EP_filt);
sorted_mod_prates = nan(size(mod_prates_ms));
sorted_full_psth = nan(size(full_psth_ms));
for tt = 1:length(rpt_taxis)
    sorted_mod_prates(:,tt,:) = mod_prates_ms(b(:,tt),tt,:);
    sorted_mod_prates(isnan(a(:,tt)),tt,:) = nan;
    
    sorted_full_psth(:,tt,:) = full_psth_ms(b(:,tt),tt,:);
    sorted_full_psth(isnan(a(:,tt)),tt,:) = nan;
end
