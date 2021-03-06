%
% clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements//');
addpath('~/James_scripts/bruce/processing/');
addpath('~/James_scripts/bruce/saccade_modulation/');
addpath('~/James_scripts/TentBasis2D/');

global Expt_name bar_ori use_MUA


% Expt_name = 'G093';
% use_MUA = false;
% bar_ori = 0; %bar orientation to use (only for UA recs)

mod_data_name = 'corrected_models2';
base_fname = 'sac_info_timing_noXV3';

include_bursts = 0;

if include_bursts
    base_fname = [base_fname '_withbursts'];
end

%%
lambda_off = 4;
lambda_gain = 3;

Nrpts = 5; %number of random jitterings of saccade times to perform
sac_jit_amp = 200; %in units of dt

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
Fr = 1;

backlag = round(0.1/dt);
forlag = round(0.3/dt);
slags = -backlag:forlag;
n_sac_bins = length(slags);


info_backlag = round(0.3/dt);
info_forlag = round(0.3/dt);
info_slags = -info_backlag:info_forlag;


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
%that aren't blinks
sac_durs = [saccades(:).duration];
big_sacs = find(abs(sac_deltaX) > gsac_thresh & ~used_is_blink' & ~out_bounds & sac_durs <= max_gsac_dur);

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
fname = 'sacStimProcFin_noXV';
if include_bursts
    fname = [fname '_withbursts'];
end
fname = [fname sprintf('_ori%d',bar_ori)];

anal_dir = ['~/Analysis/bruce/' Expt_name '/FINsac_mod/'];
cd(anal_dir)
load(fname);
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
for i=1:NT
    best_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
end

%%
un_trials = unique(all_trialvec);
n_trials = length(un_trials);
all_ttill_end = nan(size(all_tsince_start));
for ii = 1:n_trials
    cur_set = find(all_trialvec == un_trials(ii));
    all_ttill_end(cur_set) = (length(cur_set):-1:1)*dt;
end
%% MAIN ANALYSIS LOOP
cd(anal_dir)
silent = 1;

for cc = targs
    
    fprintf('Starting model fits for unit %d\n',cc);
    loo_cc = find(loo_set == cc); %index within the LOOXV set
    cc_uinds = full_inds(~isnan(Robs_mat(full_inds,cc))); %set of used indices where this unit was isolated
    cur_tr_inds = find(ismember(cc_uinds,tr_inds)); %subset for training
    cur_xv_inds = find(ismember(cc_uinds,xv_inds)); %subset for XV
    
    cur_Robs = Robs_mat(cc_uinds,cc);
    
    %%
    if ~isempty(cc_uinds)
        %%
        cur_GQM = ModData(cc).rectGQM;
        sacStimProc(cc).ModData = ModData(cc);
        sacStimProc(cc).used = true;
        
        fprintf('Reconstructing retinal stim for unit %d\n',cc);
        if ismember(cc,loo_set) %if unit is member of LOOXV set, use its unique EP sequence
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
        
        %%
        %lag index values in Xmat
        [~,Tinds] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
        Tinds = Tinds(use_kInds_up);
        
        cur_sac_start_inds = saccade_start_inds(big_sacs);
        cur_sac_stop_inds = saccade_stop_inds(big_sacs);
        is_insac = false(NT,1);
        for ii = 1:length(cur_sac_start_inds)
            is_insac(cur_sac_start_inds(ii):cur_sac_stop_inds(ii)) = true;
        end
        
        cur_msac_start_inds = saccade_start_inds(micro_sacs);
        cur_msac_stop_inds = saccade_stop_inds(micro_sacs);
        is_inmsac = false(NT,1);
        for ii = 1:length(cur_msac_start_inds)
            is_inmsac(cur_msac_start_inds(ii):cur_msac_start_inds(ii)) = true;
        end
        is_in_anysac = is_insac | is_inmsac;
        
        %         %realign the detected end points of saccades to better match
        %         %downslope of eye speed
        %         min_eye_speed = 10;
        %         for ii = 1:length(cur_sac_start_inds)
        %             if interp_eye_speed(used_inds(cur_sac_stop_inds(ii))) < min_eye_speed
        %                 cur_sac_stop_inds(ii) = 1+find(interp_eye_speed(used_inds(1:cur_sac_stop_inds(ii))) >= min_eye_speed,1,'last');
        %             end
        %         end
        
        %%
        saccade_stop_trial_inds = all_trialvec(used_inds(cur_sac_stop_inds));
        saccade_start_trial_inds = all_trialvec(used_inds(cur_sac_start_inds));
        cur_sac_durs = cur_sac_stop_inds - cur_sac_start_inds;
        sacStimProc(cc).avg_gsac_dur = mean(cur_sac_durs);
        sacStimProc(cc).prc_gsac_durs = prctile(cur_sac_durs,[25 50 75]);
        
        Xsac_end = zeros(NT,length(info_slags));
        for ii = 1:length(info_slags)
            cur_sac_target = cur_sac_stop_inds + info_slags(ii);
            uu = find(cur_sac_target > 1 & cur_sac_target < NT);
            cur_sac_target = cur_sac_target(uu);
            cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_stop_trial_inds(uu)) = [];
            Xsac_end(cur_sac_target,ii) = 1;
        end
        Xsac_end = Xsac_end(cc_uinds,:);
        
        %store duration of previous saccade
        prev_sac_dur = ones(NT,1)*cur_sac_durs(end);
        for ii = 1:length(cur_sac_start_inds)-1
            cur_inds = cur_sac_start_inds(ii):cur_sac_start_inds(ii+1);
            next_tstop = trial_end_inds(find(trial_end_inds >= cur_sac_start_inds(ii),1));
            cur_inds(cur_inds > next_tstop) = [];
            prev_sac_dur(cur_inds) = cur_sac_durs(ii);
        end
        
        Xsac_start = zeros(NT,length(info_slags));
        for ii = 1:length(info_slags)
            cur_sac_target = cur_sac_start_inds + info_slags(ii);
            uu = find(cur_sac_target > 1 & cur_sac_target < NT);
            cur_sac_target = cur_sac_target(uu);
            cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_start_trial_inds(uu)) = [];
            Xsac_start(cur_sac_target,ii) = 1;
        end
        Xsac_start = Xsac_start(cc_uinds,:);
        
        %%
        cur_Xsac = Xsac(cc_uinds,:);
        %         any_sac_inds = find(any(cur_Xsac > 0,2));
        any_sac_inds = find(any(Xsac_end > 0,2));
        if ~isempty(any_sac_inds) && ~isempty(sacStimProc(cc).gsac_avg_rate)
            
            cur_GQM = NMMfit_logexp_spkNL(cur_GQM,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:));
            cur_GQM = NMMfit_scale(cur_GQM,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:));
            
            stim_mod_signs = [cur_GQM.mods(:).sign];
            [LL,~,~,~,filt_outs,fgint,nullLL] = NMMmodel_eval(cur_GQM,cur_Robs,all_Xmat_shift);
            fgint = bsxfun(@times,fgint,stim_mod_signs);
            stimG = sum(fgint,2);
            
            tr_stim{1} = stimG;
            tr_stim{2} = cur_Xsac; %saccade timing indicator matrix
            tr_stim{3} = bsxfun(@times,cur_Xsac,stimG);
            
            %% NOW DO INFO TIMING CALCS
            %create a shuffled version of the stimulus, where different frames
            %(within the original used_inds set) are randomly shuffled with
            %replacement
            shuf_stim = all_shift_stimmat_up;
            shuf_stim(used_inds,:) = all_shift_stimmat_up(used_inds(randi(NT,NT,1)),:);
            shuf_X = create_time_embedding(shuf_stim,stim_params_us);
            shuf_X = shuf_X(used_inds(cc_uinds),use_kInds_up);
            
            %%
            
            n_pregain_regs = length(sacStimProc(cc).gsacPreGainMod);
            %             for ll = 1:n_pregain_regs
            ll = 3;
            fprintf('Computing info timing for preGain Reg %d/%d\n',ll,n_pregain_regs);
            sacGainMod = sacStimProc(cc).gsacPreGainMod{ll};
            
            %initialize simple 1-parameter GLM for gain/offset fitting
            temp_sparams = NMMcreate_stim_params(1);
            tempmod = NMMinitialize_model(temp_sparams,1,{'lin'});
            
            stim_mod_signs = [sacGainMod.stim_mod.mods(:).sign];
            
            [sac_avg_rate,sac_nspks,LL_imp_before,LL_imp_after,LL_imp_during] = deal(nan(length(info_slags),1));
            for ii = 1:length(info_slags)
                cur_set = find(Xsac_end(:,ii) == 1);
                cur_X = all_Xmat_shift(cur_set,:);
                
                sac_avg_rate(ii) = mean(cur_Robs(cur_set));
                sac_nspks(ii) = sum(cur_Robs(cur_set));
                if sac_avg_rate(ii) > 0 %so long as there is at least one spike at this latency
                    
                    %get the G-sig from the sacMod, fit a GLM to it,
                    %and compute the unshuffled LL
                    [unshuf_LL,unshuf_pred_rate,G,fG] = eval_pre_gainmodel( sacGainMod, cur_Robs(cur_set), cur_X, cur_Xsac(cur_set,:));
                    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),G);
                    tempmod = NMMfit_logexp_spkNL(tempmod,cur_Robs(cur_set),G);
                    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),G);
                    unshuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                    null_pred = mean(cur_Robs(cur_set));
                    null_LLseq = cur_Robs(cur_set).*log2(null_pred)-null_pred;
                    unshuf_LLimp(ii) = nansum(unshuf_LLseq-null_LLseq)/sum(cur_Robs(cur_set));
                    
                    %find Xmat indices before the current sacs and
                    %shuffle these
                    cur_Tinds = bsxfun(@plus,-prev_sac_dur(cur_set),Tinds'-1);
                    scramb_lags = find(cur_Tinds > info_slags(ii));
                    cur_shufX = shuf_X(cur_set,:);
                    scr_X = cur_X;
                    scr_X(scramb_lags) = cur_shufX(scramb_lags);
                    %compute scrambled G-sig and compute the shuffled
                    %LL of a GLM-fit
                    [shuf_LL,shuf_pred_rate,Gshuff,fGshuff] = eval_pre_gainmodel( sacGainMod, cur_Robs(cur_set), scr_X, cur_Xsac(cur_set,:));
                    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
                     tempmod = NMMfit_logexp_spkNL(tempmod,cur_Robs(cur_set),Gshuff);
                   [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
                    shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                    %compute LL improvement (per spike) from unshuffling
                    LL_imp_before(ii) = nansum(unshuf_LLseq - shuf_LLseq)/sum(cur_Robs(cur_set));
                    
                    %find set of Xmat indices after fix onset and
                    %scramble these
                    scramb_lags = find((Tinds-1) < info_slags(ii));
                    scr_X = cur_X;
                    scr_X(:,scramb_lags) = shuf_X(cur_set,scramb_lags);
                    [shuf_LL,shuf_pred_rate,Gshuff,fGshuff] = eval_pre_gainmodel( sacGainMod, cur_Robs(cur_set), scr_X, cur_Xsac(cur_set,:));
                    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
                     tempmod = NMMfit_logexp_spkNL(tempmod,cur_Robs(cur_set),Gshuff);
                   [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
                    shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                    LL_imp_after(ii) = nansum(unshuf_LLseq - shuf_LLseq)/sum(cur_Robs(cur_set));
                    
                    poss_scramb_lags = find((Tinds-1) >= info_slags(ii)); %set of lags before fixation onset
                    
                    %set of lags after saccade onset
                    cur_Tinds = bsxfun(@plus,-prev_sac_dur(cur_set),Tinds'-1);
                    scramb_lags = find(cur_Tinds <= info_slags(ii));
                    
                    %these are the set of lags during the saccade
                    [scramb_i,scramb_j] = ind2sub([length(cur_set) length(Tinds)],scramb_lags);
                    scramb_lags = scramb_lags(ismember(scramb_j,poss_scramb_lags));
                    
                    scr_X = cur_X;
                    cur_shufX = shuf_X(cur_set,:);
                    scr_X(scramb_lags) = cur_shufX(scramb_lags);
                    [shuf_LL,shuf_pred_rate,Gshuff,fGshuff] = eval_pre_gainmodel( sacGainMod, cur_Robs(cur_set), scr_X, cur_Xsac(cur_set,:));
                    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
                     tempmod = NMMfit_logexp_spkNL(tempmod,cur_Robs(cur_set),Gshuff);
                   [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
                    shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                    LL_imp_during(ii) = nansum(unshuf_LLseq - shuf_LLseq)/sum(cur_Robs(cur_set));
                end
            end
            
            sacInfoTiming(cc).gsac_info_before(ll,:) = LL_imp_before;
            sacInfoTiming(cc).gsac_info_after(ll,:) = LL_imp_after;
            sacInfoTiming(cc).gsac_info_during(ll,:) = LL_imp_during;
            sacInfoTiming(cc).gsac_avg_rate = sac_avg_rate;
            %             end
            %%
            
            [base_LL_imp_before,base_LL_imp_after,base_LL_imp_during,base_unshuf_LLimp] = deal(nan(Nrpts,length(info_slags)));
            for rr = 1:Nrpts
                fprintf('Scrambling with rand saclocs, iter %d of %d\n',rr,Nrpts);
                
                rand_jit = randi(sac_jit_amp,length(cur_sac_start_inds),1);
                % rand_jit = -ones(length(cur_sac_start_inds),1)*30;
                rcur_sac_start_inds = cur_sac_start_inds + rand_jit;
                rcur_sac_stop_inds = cur_sac_stop_inds + rand_jit;
                bad = find(rcur_sac_stop_inds > NT);
                rcur_sac_start_inds(bad) = [];
                rcur_sac_stop_inds(bad) = [];
               
                %if the sim fixation starts too close to the end of a
                %trial, eliminate it
                bad = find(all_ttill_end(used_inds(rcur_sac_stop_inds)) <= 0.3);
                rcur_sac_start_inds(bad) = [];
                rcur_sac_stop_inds(bad) = [];

%                 %eliminate any sim sacs that include any real sac indices
%                 to_elim = [];
%                 for ii = 1:length(rcur_sac_start_inds)
%                     cur_set = rcur_sac_start_inds(ii):rcur_sac_stop_inds(ii);
%                     %                     if any(is_insac(cur_set))
%                     if any(is_in_anysac(cur_set))
%                         to_elim = cat(1,to_elim,ii);
%                     end
%                 end
%                 rcur_sac_start_inds(to_elim) = [];
%                 rcur_sac_stop_inds(to_elim) = [];
                
                rsaccade_stop_trial_inds = all_trialvec(used_inds(rcur_sac_stop_inds));
                
                Xsac_rend = zeros(NT,length(info_slags));
                for ii = 1:length(info_slags)
                    cur_sac_target = rcur_sac_stop_inds + info_slags(ii);
                    uu = find(cur_sac_target > 1 & cur_sac_target < NT);
                    cur_sac_target = cur_sac_target(uu);
                    cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= rsaccade_stop_trial_inds(uu)) = [];
                    Xsac_rend(cur_sac_target,ii) = 1;
                end
                Xsac_rend = Xsac_rend(cc_uinds,:);
                
                cur_sac_durs = rcur_sac_stop_inds - rcur_sac_start_inds;
                %store duration of previous saccade
                rprev_sac_dur = ones(NT,1)*cur_sac_durs(end);
                for ii = 1:length(rcur_sac_start_inds)-1
                    cur_inds = rcur_sac_start_inds(ii):rcur_sac_start_inds(ii+1);
                    next_tstop = trial_end_inds(find(trial_end_inds >= rcur_sac_start_inds(ii),1));
                    cur_inds(cur_inds > next_tstop) = [];
                    rprev_sac_dur(cur_inds) = cur_sac_durs(ii);
                end
                
                Xsac_end = Xsac_rend;
                prev_sac_dur = rprev_sac_dur;
                %                 rand_shift = randi(sac_jit_amp,1,1);
                rand_shift = 0;
                rand_Robs = circshift(cur_Robs,rand_shift);
                rand_X = circshift(all_Xmat_shift,rand_shift);
                
                for ii = 1:length(info_slags)
                    
                    cur_set = find(Xsac_end(:,ii) == 1);
                    
                    if sum(rand_Robs(cur_set)) > 0
                        cur_X = rand_X(cur_set,:);
                        
                        [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,rand_Robs(cur_set),cur_X);
                        G = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
                        tempmod = NMMfit_filters(tempmod,rand_Robs(cur_set),G);
                        tempmod = NMMfit_logexp_spkNL(tempmod,cur_Robs(cur_set),G);
                        [~,~,tempprate] = NMMmodel_eval(tempmod,rand_Robs(cur_set),G);
                        %                         tempprate = base_shuff_prate;
                        base_unshuf_LLseq = rand_Robs(cur_set).*log2(tempprate)-tempprate;
                        null_pred = mean(rand_Robs(cur_set));
                        null_LLseq = rand_Robs(cur_set).*log2(null_pred)-null_pred;
                        base_unshuf_LLimp(rr,ii) = nansum(base_unshuf_LLseq-null_LLseq)/sum(rand_Robs(cur_set));
                        
                        cur_Tinds = bsxfun(@plus,-prev_sac_dur(cur_set),Tinds'-1);
                        scramb_lags = find(cur_Tinds > info_slags(ii));
                        scr_X = cur_X;
                        cur_shufX = shuf_X(cur_set,:);
                        scr_X(scramb_lags) = cur_shufX(scramb_lags);
                        
                        [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,rand_Robs(cur_set),scr_X);
                        Gshuff = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
                        tempmod = NMMfit_filters(tempmod,rand_Robs(cur_set),Gshuff);
                        tempmod = NMMfit_logexp_spkNL(tempmod,cur_Robs(cur_set),Gshuff);
                        [~,~,tempprate] = NMMmodel_eval(tempmod,rand_Robs(cur_set),Gshuff);
                        base_shuf_LLseq = rand_Robs(cur_set).*log2(tempprate)-tempprate;
                        base_LL_imp_before(rr,ii) = nansum(base_unshuf_LLseq - base_shuf_LLseq)/sum(rand_Robs(cur_set));
                        
                        
                        scramb_lags = find((Tinds-1) < info_slags(ii));
                        scr_X = cur_X;
                        scr_X(:,scramb_lags) = shuf_X(cur_set,scramb_lags);
                        
                        [base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,rand_Robs(cur_set),all_Xmat_shift(cur_set,:));
                        [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,rand_Robs(cur_set),scr_X);
                        Gshuff = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
                        tempmod = NMMfit_filters(tempmod,rand_Robs(cur_set),Gshuff);
                        tempmod = NMMfit_logexp_spkNL(tempmod,cur_Robs(cur_set),Gshuff);
                        [~,~,tempprate] = NMMmodel_eval(tempmod,rand_Robs(cur_set),Gshuff);
                        base_shuf_LLseq = rand_Robs(cur_set).*log2(tempprate)-tempprate;
                        base_LL_imp_after(rr,ii) = nansum(base_unshuf_LLseq - base_shuf_LLseq)/sum(rand_Robs(cur_set));
                        
                        
                        poss_scramb_lags = find((Tinds-1) >= info_slags(ii));
                        cur_Tinds = bsxfun(@plus,-prev_sac_dur(cur_set),Tinds'-1);
                        scramb_lags = find(cur_Tinds <= info_slags(ii));
                        [scramb_i,scramb_j] = ind2sub([length(cur_set) length(Tinds)],scramb_lags);
                        scramb_lags = scramb_lags(ismember(scramb_j,poss_scramb_lags));
                        scr_X = cur_X;
                        cur_shufX = shuf_X(cur_set,:);
                        scr_X(scramb_lags) = cur_shufX(scramb_lags);
                        
                        [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,rand_Robs(cur_set),scr_X);
                        Gshuff = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
                        tempmod = NMMfit_filters(tempmod,rand_Robs(cur_set),Gshuff);
                        tempmod = NMMfit_logexp_spkNL(tempmod,cur_Robs(cur_set),Gshuff);
                        [~,~,tempprate] = NMMmodel_eval(tempmod,rand_Robs(cur_set),Gshuff);
                        base_shuf_LLseq = rand_Robs(cur_set).*log2(tempprate)-tempprate;
                        base_LL_imp_during(rr,ii) = nansum(base_unshuf_LLseq - base_shuf_LLseq)/sum(rand_Robs(cur_set));
                    end
                end
            end
            sacInfoTiming(cc).gsac_base_info_before = nanmean(base_LL_imp_before,1);
            sacInfoTiming(cc).gsac_base_info_after = nanmean(base_LL_imp_after,1);
            sacInfoTiming(cc).gsac_base_info_during = nanmean(base_LL_imp_during,1);
            sacInfoTiming(cc).gsac_base_info_unshuff = nanmean(base_unshuf_LLimp,1);
        end
        
        
        %% MICROSACS
        %lag index values in Xmat
        cur_sac_start_inds = saccade_start_inds(micro_sacs);
        cur_sac_stop_inds = saccade_stop_inds(micro_sacs);
        
        %         %realign the detected end points of saccades to better match
        %         %downslope of eye speed
        %         min_eye_speed = 3;
        %         for ii = 1:length(cur_sac_start_inds)
        %             if interp_eye_speed(used_inds(cur_sac_stop_inds(ii))) < min_eye_speed
        %                 cur_sac_stop_inds(ii) = 1+find(interp_eye_speed(used_inds(1:cur_sac_stop_inds(ii))) >= min_eye_speed,1,'last');
        %             end
        %         end
        
        %%
        saccade_stop_trial_inds = all_trialvec(used_inds(cur_sac_stop_inds));
        saccade_start_trial_inds = all_trialvec(used_inds(cur_sac_start_inds));
        cur_sac_durs = cur_sac_stop_inds - cur_sac_start_inds;
        sacStimProc(cc).avg_msac_dur = mean(cur_sac_durs);
        sacStimProc(cc).prc_msac_durs = prctile(cur_sac_durs,[25 50 75]);
        
        Xsac_end = zeros(NT,length(info_slags));
        for ii = 1:length(info_slags)
            cur_sac_target = cur_sac_stop_inds + info_slags(ii);
            uu = find(cur_sac_target > 1 & cur_sac_target < NT);
            cur_sac_target = cur_sac_target(uu);
            cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_stop_trial_inds(uu)) = [];
            Xsac_end(cur_sac_target,ii) = 1;
        end
        Xsac_end = Xsac_end(cc_uinds,:);
        
        %store duration of previous saccade
        prev_sac_dur = ones(NT,1)*cur_sac_durs(end);
        for ii = 1:length(cur_sac_start_inds)-1
            cur_inds = cur_sac_start_inds(ii):cur_sac_start_inds(ii+1);
            next_tstop = trial_end_inds(find(trial_end_inds >= cur_sac_start_inds(ii),1));
            cur_inds(cur_inds > next_tstop) = [];
            prev_sac_dur(cur_inds) = cur_sac_durs(ii);
        end
        
        Xsac_start = zeros(NT,length(info_slags));
        for ii = 1:length(info_slags)
            cur_sac_target = cur_sac_start_inds + info_slags(ii);
            uu = find(cur_sac_target > 1 & cur_sac_target < NT);
            cur_sac_target = cur_sac_target(uu);
            cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_start_trial_inds(uu)) = [];
            Xsac_start(cur_sac_target,ii) = 1;
        end
        Xsac_start = Xsac_start(cc_uinds,:);
        
        
        cur_Xsac = Xmsac(cc_uinds,:);
        any_sac_inds = find(any(Xsac_end > 0,2));
        
        %%
        if ~isempty(any_sac_inds) && ~isempty(sacStimProc(cc).msac_avg_rate)
            
            cur_GQM = NMMfit_logexp_spkNL(cur_GQM,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:));
            cur_GQM = NMMfit_scale(cur_GQM,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:));
            
            stim_mod_signs = [cur_GQM.mods(:).sign];
            [~,~,~,~,filt_outs,fgint] = NMMmodel_eval(cur_GQM,cur_Robs,all_Xmat_shift);
            fgint = bsxfun(@times,fgint,stim_mod_signs);
            stimG = sum(fgint,2);
            
            tr_stim{1} = stimG;
            tr_stim{2} = cur_Xsac; %saccade timing indicator matrix
            tr_stim{3} = bsxfun(@times,cur_Xsac,stimG);
            
            
            %% NOW DO INFO TIMING CALCS
            %create a shuffled version of the stimulus, where different frames
            %(within the original used_inds set) are randomly shuffled with
            %replacement
            shuf_stim = all_shift_stimmat_up;
            shuf_stim(used_inds,:) = all_shift_stimmat_up(used_inds(randi(NT,NT,1)),:);
            shuf_X = create_time_embedding(shuf_stim,stim_params_us);
            shuf_X = shuf_X(used_inds(cc_uinds),use_kInds_up);
            
            %%
            n_pregain_regs = length(sacStimProc(cc).msacPreGainMod);
            %             for ll = 1:n_pregain_regs
            ll = 3;
            fprintf('Computing info timing for preGain Reg %d/%d\n',ll,n_pregain_regs);
            sacGainMod = sacStimProc(cc).msacPreGainMod{ll};
            
            %initialize simple 1-parameter GLM for gain/offset fitting
            temp_sparams = NMMcreate_stim_params(1);
            tempmod = NMMinitialize_model(temp_sparams,1,{'lin'});
            
            stim_mod_signs = [sacGainMod.stim_mod.mods(:).sign];
            
            [sac_avg_rate,LL_imp_before,LL_imp_after,LL_imp_during] = deal(nan(length(info_slags),1));
            for ii = 1:length(info_slags)
                cur_set = find(Xsac_end(:,ii) == 1);
                cur_X = all_Xmat_shift(cur_set,:);
                
                sac_avg_rate(ii) = mean(cur_Robs(cur_set));
                if sac_avg_rate(ii) > 0
                    [unshuf_LL,unshuf_pred_rate,G,fG] = eval_pre_gainmodel( sacGainMod, cur_Robs(cur_set), cur_X, cur_Xsac(cur_set,:));
                    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),G);
                    tempmod = NMMfit_logexp_spkNL(tempmod,cur_Robs(cur_set),G);
                    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),G);
                    unshuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                    
                    cur_Tinds = bsxfun(@plus,-prev_sac_dur(cur_set),Tinds'-1);
                    scramb_lags = find(cur_Tinds > info_slags(ii));
                    cur_shufX = shuf_X(cur_set,:);
                    scr_X = cur_X;
                    scr_X(scramb_lags) = cur_shufX(scramb_lags);
                    [shuf_LL,shuf_pred_rate,Gshuff,fGshuff] = eval_pre_gainmodel( sacGainMod, cur_Robs(cur_set), scr_X, cur_Xsac(cur_set,:));
                    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
                    tempmod = NMMfit_logexp_spkNL(tempmod,cur_Robs(cur_set),Gshuff);
                    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
                    shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                    LL_imp_before(ii) = nansum(unshuf_LLseq - shuf_LLseq)/sum(cur_Robs(cur_set));
                    
                    scramb_lags = find((Tinds-1) < info_slags(ii));
                    scr_X = cur_X;
                    scr_X(:,scramb_lags) = shuf_X(cur_set,scramb_lags);
                    [shuf_LL,shuf_pred_rate,Gshuff,fGshuff] = eval_pre_gainmodel( sacGainMod, cur_Robs(cur_set), scr_X, cur_Xsac(cur_set,:));
                    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
                     tempmod = NMMfit_logexp_spkNL(tempmod,cur_Robs(cur_set),Gshuff);
                   [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
                    shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                    LL_imp_after(ii) = nansum(unshuf_LLseq - shuf_LLseq)/sum(cur_Robs(cur_set));
                    
                    poss_scramb_lags = find((Tinds-1) >= info_slags(ii)); %set of lags before fixation onset
                    
                    %set of lags after saccade onset
                    cur_Tinds = bsxfun(@plus,-prev_sac_dur(cur_set),Tinds'-1);
                    scramb_lags = find(cur_Tinds <= info_slags(ii));
                    
                    [scramb_i,scramb_j] = ind2sub([length(cur_set) length(Tinds)],scramb_lags);
                    scramb_lags = scramb_lags(ismember(scramb_j,poss_scramb_lags));
                    
                    scr_X = cur_X;
                    cur_shufX = shuf_X(cur_set,:);
                    scr_X(scramb_lags) = cur_shufX(scramb_lags);
                    [shuf_LL,shuf_pred_rate,Gshuff,fGshuff] = eval_pre_gainmodel( sacGainMod, cur_Robs(cur_set), scr_X, cur_Xsac(cur_set,:));
                    tempmod = NMMfit_filters(tempmod,cur_Robs(cur_set),Gshuff);
                    tempmod = NMMfit_logexp_spkNL(tempmod,cur_Robs(cur_set),Gshuff);
                    [~,~,tempprate] = NMMmodel_eval(tempmod,cur_Robs(cur_set),Gshuff);
                    shuf_LLseq = cur_Robs(cur_set).*log2(tempprate)-tempprate;
                    LL_imp_during(ii) = nansum(unshuf_LLseq - shuf_LLseq)/sum(cur_Robs(cur_set));
                end
            end
            
            sacInfoTiming(cc).msac_info_before(ll,:) = LL_imp_before;
            sacInfoTiming(cc).msac_info_after(ll,:) = LL_imp_after;
            sacInfoTiming(cc).msac_info_during(ll,:) = LL_imp_during;
            sacInfoTiming(cc).msac_avg_rate = sac_avg_rate;
            %             end
            %%
            [base_LL_imp_before,base_LL_imp_after,base_LL_imp_during,base_unshuf_LLimp] = deal(nan(Nrpts,length(info_slags)));
            for rr = 1:Nrpts
                fprintf('Scrambling with rand saclocs, iter %d of %d\n',rr,Nrpts);
                
                rand_jit = randi(sac_jit_amp,length(cur_sac_start_inds),1);
                rcur_sac_start_inds = cur_sac_start_inds + rand_jit;
                rcur_sac_stop_inds = cur_sac_stop_inds + rand_jit;
                bad = find(rcur_sac_stop_inds > NT);
                rcur_sac_start_inds(bad) = [];
                rcur_sac_stop_inds(bad) = [];
          
                %if the sim fixation starts too close to the end of a
                %trial, eliminate it
                bad = find(all_ttill_end(used_inds(rcur_sac_stop_inds)) <= 0.3);
                rcur_sac_start_inds(bad) = [];
                rcur_sac_stop_inds(bad) = [];

                %eliminate any sim sacs that include any real sac indices
                to_elim = [];
                for ii = 1:length(rcur_sac_start_inds)
                    cur_set = rcur_sac_start_inds(ii):rcur_sac_stop_inds(ii);
                    %                     if any(is_insac(cur_set))
                    if any(is_in_anysac(cur_set))
                        to_elim = cat(1,to_elim,ii);
                    end
                end
                rcur_sac_start_inds(to_elim) = [];
                rcur_sac_stop_inds(to_elim) = [];

                
                rsaccade_stop_trial_inds = all_trialvec(used_inds(rcur_sac_stop_inds));
                
                Xsac_rend = zeros(NT,length(info_slags));
                for ii = 1:length(info_slags)
                    cur_sac_target = rcur_sac_stop_inds + info_slags(ii);
                    uu = find(cur_sac_target > 1 & cur_sac_target < NT);
                    cur_sac_target = cur_sac_target(uu);
                    cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= rsaccade_stop_trial_inds(uu)) = [];
                    Xsac_rend(cur_sac_target,ii) = 1;
                end
                Xsac_rend = Xsac_rend(cc_uinds,:);
                
                cur_sac_durs = rcur_sac_stop_inds - rcur_sac_start_inds;
                %store duration of previous saccade
                rprev_sac_dur = ones(NT,1)*cur_sac_durs(end);
                for ii = 1:length(rcur_sac_start_inds)-1
                    cur_inds = rcur_sac_start_inds(ii):rcur_sac_start_inds(ii+1);
                    next_tstop = trial_end_inds(find(trial_end_inds >= rcur_sac_start_inds(ii),1));
                    cur_inds(cur_inds > next_tstop) = [];
                    rprev_sac_dur(cur_inds) = cur_sac_durs(ii);
                end
                
                Xsac_end = Xsac_rend;
                prev_sac_dur = rprev_sac_dur;
                %                 rand_shift = randi(sac_jit_amp,1,1);
                rand_shift = 0;
                rand_Robs = circshift(cur_Robs,rand_shift);
                rand_X = circshift(all_Xmat_shift,rand_shift);
                
                for ii = 1:length(info_slags)
                    
                    cur_set = find(Xsac_end(:,ii) == 1);
                    
                    if sum(rand_Robs(cur_set)) > 0
                        %                         cur_X = all_Xmat_shift(cur_set,:);
                        cur_X = rand_X(cur_set,:);
                        
                        [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,rand_Robs(cur_set),cur_X);
                        G = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
                        tempmod = NMMfit_filters(tempmod,rand_Robs(cur_set),G);
                      tempmod = NMMfit_logexp_spkNL(tempmod,cur_Robs(cur_set),G);
                      [~,~,tempprate] = NMMmodel_eval(tempmod,rand_Robs(cur_set),G);
                        base_unshuf_LLseq = rand_Robs(cur_set).*log2(tempprate)-tempprate;
                        null_pred = mean(rand_Robs(cur_set));
                        null_LLseq = rand_Robs(cur_set).*log2(null_pred)-null_pred;
                        base_unshuf_LLimp(rr,ii) = nansum(base_unshuf_LLseq-null_LLseq)/sum(rand_Robs(cur_set));
                        
                        cur_Tinds = bsxfun(@plus,-prev_sac_dur(cur_set),Tinds'-1);
                        scramb_lags = find(cur_Tinds > info_slags(ii));
                        scr_X = cur_X;
                        cur_shufX = shuf_X(cur_set,:);
                        scr_X(scramb_lags) = cur_shufX(scramb_lags);
                        
                        [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,rand_Robs(cur_set),scr_X);
                        Gshuff = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
                        tempmod = NMMfit_filters(tempmod,rand_Robs(cur_set),Gshuff);
                     tempmod = NMMfit_logexp_spkNL(tempmod,cur_Robs(cur_set),Gshuff);
                       [~,~,tempprate] = NMMmodel_eval(tempmod,rand_Robs(cur_set),Gshuff);
                        base_shuf_LLseq = rand_Robs(cur_set).*log2(tempprate)-tempprate;
                        base_LL_imp_before(rr,ii) = nansum(base_unshuf_LLseq - base_shuf_LLseq)/sum(rand_Robs(cur_set));
                        
                        
                        scramb_lags = find((Tinds-1) < info_slags(ii));
                        scr_X = cur_X;
                        scr_X(:,scramb_lags) = shuf_X(cur_set,scramb_lags);
                        
                        [base_LL,~,base_predrate] = NMMmodel_eval(cur_GQM,rand_Robs(cur_set),all_Xmat_shift(cur_set,:));
                        [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,rand_Robs(cur_set),scr_X);
                        Gshuff = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
                        tempmod = NMMfit_filters(tempmod,rand_Robs(cur_set),Gshuff);
                     tempmod = NMMfit_logexp_spkNL(tempmod,cur_Robs(cur_set),Gshuff);
                       [~,~,tempprate] = NMMmodel_eval(tempmod,rand_Robs(cur_set),Gshuff);
                        base_shuf_LLseq = rand_Robs(cur_set).*log2(tempprate)-tempprate;
                        base_LL_imp_after(rr,ii) = nansum(base_unshuf_LLseq - base_shuf_LLseq)/sum(rand_Robs(cur_set));
                        
                        
                        poss_scramb_lags = find((Tinds-1) >= info_slags(ii));
                        cur_Tinds = bsxfun(@plus,-prev_sac_dur(cur_set),Tinds'-1);
                        scramb_lags = find(cur_Tinds <= info_slags(ii));
                        [scramb_i,scramb_j] = ind2sub([length(cur_set) length(Tinds)],scramb_lags);
                        scramb_lags = scramb_lags(ismember(scramb_j,poss_scramb_lags));
                        scr_X = cur_X;
                        cur_shufX = shuf_X(cur_set,:);
                        scr_X(scramb_lags) = cur_shufX(scramb_lags);
                        
                        [base_shufLL,~,base_shuff_prate,~,~,temp_fgint] = NMMmodel_eval(cur_GQM,rand_Robs(cur_set),scr_X);
                        Gshuff = sum(bsxfun(@times,temp_fgint,[cur_GQM.mods(:).sign]),2);
                        tempmod = NMMfit_filters(tempmod,rand_Robs(cur_set),Gshuff);
                      tempmod = NMMfit_logexp_spkNL(tempmod,cur_Robs(cur_set),Gshuff);
                      [~,~,tempprate] = NMMmodel_eval(tempmod,rand_Robs(cur_set),Gshuff);
                        base_shuf_LLseq = rand_Robs(cur_set).*log2(tempprate)-tempprate;
                        base_LL_imp_during(rr,ii) = nansum(base_unshuf_LLseq - base_shuf_LLseq)/sum(rand_Robs(cur_set));
                    end
                end
            end
            sacInfoTiming(cc).msac_base_info_before = nanmean(base_LL_imp_before,1);
            sacInfoTiming(cc).msac_base_info_after = nanmean(base_LL_imp_after,1);
            sacInfoTiming(cc).msac_base_info_during = nanmean(base_LL_imp_during,1);
            sacInfoTiming(cc).msac_base_info_unshuff = nanmean(base_unshuf_LLimp,1);
        end
    end
    sacInfoTiming(cc).lag_axis = info_slags*dt;
    
end
%%
cd(anal_dir)
fname = [base_fname sprintf('_ori%d',bar_ori)];
save(fname,'sacInfoTiming','info_slags','dt');