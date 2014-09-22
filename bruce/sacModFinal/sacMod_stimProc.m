
%
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements//');
addpath('~/James_scripts/bruce/processing/');
addpath('~/James_scripts/bruce/saccade_modulation/');
addpath('~/James_scripts/bruce/sacModFinal/');
addpath('~/James_scripts/TentBasis2D/');

global Expt_name bar_ori use_MUA

% Expt_name = 'M297';
Expt_name = 'G093';
use_MUA = false;
bar_ori = 0; %bar orientation to use (only for UA recs)


fit_unCor = false;
fit_subMod = true;
fitUpstream = true;
fitSTA = true;
fitMsacs = true;
fit_msacUpstream = true;
fitFullPostMod = true;

include_bursts = 0;

sname = 'sacStimProcFin';
if include_bursts
    sname = [sname '_withbursts'];
end

mod_data_name = 'corrected_models2';

%%
poss_gain_d2T = logspace(log10(1),log10(1e3),8); %range of d2T reg values for post-gain models
poss_gain_L2 = [0 1 5 10]; %range of L2 reg values 
poss_pre_d2T = logspace(log10(1),log10(200),6); %range of d2T reg values for pre-gain models
poss_sub_d2T = logspace(log10(10),log10(1e4),6); %range of d2T reg values for subspace models
poss_TB_lambdas = logspace(log10(0.1),log10(50),6); %range of d2T reg values for TB models
% poss_gain_d2T = 50;
% poss_pre_d2T = logspace(log10(1),log10(500),5);
% poss_gain_L2 = 0.5;
% poss_sub_d2T = 100;
% poss_TB_lambdas = 5;

%reg parameters for full model (separate gain kernels for each subunit)
fullMod_d2T = 5;
fullMod_L2 = 5;

n_Gbins = 35; %number of bins for TB model
G_lambdas = 100; %d2G reg parameter for TB model (smoothness in g-dimension)

micro_thresh = 1; %max amp of microsac (deg)
EP_bounds = 1;%eye position boundary (deg from central FP)
sac_burst_isi = 0.15; %minimum inter-saccade interval for micros (to eliminate 'bursts')
max_gsac_dur = 0.1; %maximum saccade duration before we call it a likely blink

xv_frac = 0.2; %fraction of trials to use for cross-validation

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

%% create up-sampled stimulus
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
%dont use data at beginning and end of trial
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer);

%for G093 use only data where stripe width is AT LEAST 2 deg
if strcmp(Expt_name,'G093')
    un_wi_vals = unique(all_trial_wi);
    use_wi_trials = find(all_trial_wi >= un_wi_vals(2));
    used_inds = used_inds(ismember(all_trialvec(used_inds),use_wi_trials));
end

NT = length(used_inds);

%% load in ET data and model fits
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

%% PROCESS DETECTED SACCADES
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
if ~include_bursts
    micro_sacs(ismember(micro_sacs,sacburst_set)) = []; %eliminate microsacs that are part of a 'burst'
end

%guided saccades are those whose parallel component is large enough and
%that aren't blinks (and whose duration is not too long to be suspicious
big_sacs = find(abs(sac_deltaX) > gsac_thresh & ~used_is_blink' & ~out_bounds & sac_durs <= max_gsac_dur);

saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));

msac_thresh = prctile(sac_amps(micro_sacs),50);
big_msacs = micro_sacs(sac_amps(micro_sacs) > msac_thresh);
small_msacs = micro_sacs(sac_amps(micro_sacs) < msac_thresh);

%% DEFINE FIXATION DATA (FOR RECONSTRUCTING EYE POSITION SEQUENCE)
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

%% Compile Binned Spike Vectors
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

%% Identify indices from non-repeat trials
rpt_trials = find(all_trial_Se==rpt_seed);
n_rpt_trials = length(rpt_trials);

use_trials = unique(all_trialvec(used_inds));
use_trials(ismember(use_trials,rpt_trials)) = []; %DONT USE REPEAT TRIALS!
nuse_trials = length(use_trials);
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
    if ~fit_unCor
        best_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
    end
end

%% MAIN ANALYSIS LOOP
cd(anal_dir)
silent = 1;
for cc = targs
    
    fprintf('Starting model fits for unit %d\n',cc);
    loo_cc = find(loo_set == cc); %index within the LOOXV set
    cc_uinds = full_inds(~isnan(Robs_mat(full_inds,cc))); %set of used indices where this unit was isolated
    
    use_trials = unique(all_trialvec(used_inds(cc_uinds))); %unique trials for this unit
    use_trials(ismember(use_trials,rpt_trials)) = []; %DONT USE REPEAT TRIALS!
    
    nuse_trials = length(use_trials);
    n_xv_trials = round(xv_frac*nuse_trials);
    xv_trials = randperm(nuse_trials);
    xv_trials(n_xv_trials+1:end) = [];
    xv_trials = use_trials(xv_trials);
    tr_trials = setdiff(use_trials,xv_trials);
    n_tr_trials = length(tr_trials);
    
    cur_tr_inds = find(ismember(all_trialvec(used_inds(cc_uinds)),tr_trials));
    cur_xv_inds = find(ismember(all_trialvec(used_inds(cc_uinds)),xv_trials));
    
    cur_Robs = Robs_mat(cc_uinds,cc);
    
    %%
    if ~isempty(cc_uinds)
        
        if fit_unCor %if not using ET corrections
            cur_rGQM = ModData(cc).rectGQM_unCor;
        else
            cur_rGQM = ModData(cc).rectGQM;
        end
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
            if ~fit_unCor
                for i=1:NT
                    all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-fin_shift_cor(i),2);
                end
            end
            all_Xmat_shift = create_time_embedding(all_shift_stimmat_up,stim_params_us);
            all_Xmat_shift = all_Xmat_shift(used_inds(cc_uinds),use_kInds_up);
            
        else %otherwise use overall EP sequence
            all_Xmat_shift = create_time_embedding(best_shift_stimmat_up,stim_params_us);
            all_Xmat_shift = all_Xmat_shift(used_inds(cc_uinds),use_kInds_up);
        end
                        
        %% FOR GSACS
        cur_Xsac = Xsac(cc_uinds,:); %saccade indicator Xmat
                        
        %only use indices within lagrange of a saccade
        any_sac_inds = find(any(cur_Xsac > 0,2));
        tr_sac_inds = cur_tr_inds(ismember(cur_tr_inds,any_sac_inds));
        xv_sac_inds = cur_xv_inds(ismember(cur_xv_inds,any_sac_inds));        
               
        %%
        if length(any_sac_inds) > 1e4
            %% Fit spk NL params and refit scale of each filter using target data (within trange of sacs)
            cur_rGQM = NMMfit_logexp_spkNL(cur_rGQM,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:));
            cur_rGQM = NMMfit_scale(cur_rGQM,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:));
            
            stim_mod_signs = [cur_rGQM.mods(:).sign];
            [~,~,basemod_pred_rate,~,filt_outs,fgint] = NMMmodel_eval(cur_rGQM,cur_Robs,all_Xmat_shift);
            fgint = bsxfun(@times,fgint,stim_mod_signs);
            stimG = sum(fgint,2);
            norm_stimG = zscore(stimG);
            
            sacStimProc(cc).gsac_ovavg_rate = mean(cur_Robs(any_sac_inds));
                        
            %% FOR SIMPLE POST_GAIN MODEL, SCAN RANGE OF L2s AND SELECT BEST USING XVAL LL
            
            clear tr_stim
            tr_stim{1} = [stimG]; %scalar-valued generating signal
            tr_stim{2} = cur_Xsac; %saccade timing indicator matrix
            tr_stim{3} = bsxfun(@times,cur_Xsac,stimG); %product of generating signal and saccade timing matrix
            clear sac_stim_params
            sac_stim_params(1) = NMMcreate_stim_params(1);
            sac_stim_params(2:3) = NMMcreate_stim_params(length(slags));
            mod_signs = [1 1 1];
            Xtargets = [1 2 3];
            NL_types = {'lin','lin','lin'};
            
            L2_gain_xvLL = nan(length(poss_gain_d2T),length(poss_gain_L2));
            for jj = 1:length(poss_gain_d2T)
                for ii = 1:length(poss_gain_L2)
                    sac_reg_params = NMMcreate_reg_params('lambda_d2T',poss_gain_d2T(jj),'lambda_L2',poss_gain_L2(ii),'boundary_conds',[0 0 0]);
                    cur_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
                    cur_mod.mods(1).reg_params = NMMcreate_reg_params();
                    cur_mod.spk_NL_params = cur_rGQM.spk_NL_params;
                    cur_mod = NMMfit_filters(cur_mod,cur_Robs(tr_sac_inds),get_Xcell_tInds(tr_stim,tr_sac_inds),[],[],silent);
                    L2_gain_xvLL(jj,ii) = NMMmodel_eval(cur_mod,cur_Robs(xv_sac_inds),get_Xcell_tInds(tr_stim,xv_sac_inds));
                end
            end
            
            %store values for optimal regularization
            sacStimProc(cc).gsac_spost_xvLL = L2_gain_xvLL;
            [~,optloc] = max(L2_gain_xvLL(:));
            [optloc_x,optloc_y] = ind2sub([length(poss_gain_d2T) length(poss_gain_L2)],optloc);
            opt_d2T = poss_gain_d2T(optloc_x);
            opt_L2 = poss_gain_L2(optloc_y);
            sacStimProc(cc).gsac_optL2 = opt_L2;
            sacStimProc(cc).gsac_optd2T = opt_d2T;
            
            
            sacMod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,[],Xtargets);
            sacMod.spk_NL_params = cur_rGQM.spk_NL_params;
            [sacMod] = sacMod_scan_regularization(sacMod,cur_Robs,tr_stim,tr_sac_inds,xv_sac_inds,poss_gain_d2T,poss_gain_L2);
            
            %% FIT POST-INTEGRATION GAIN USING OPTIMAL REGULARIZATION
            fprintf('Fitting post-filter models\n');
            sac_reg_params = NMMcreate_reg_params('lambda_d2T',opt_d2T,'lambda_L2',opt_L2,'boundary_conds',[0 0 0]);
            
            %if there are both significant E and I filters, fit model with sep
            %E and I gain kernels
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
                sac_stim_params(2:4) = NMMcreate_stim_params(length(slags));
                mod_signs = [1 1 1 1];
                Xtargets = [1 2 3 4];
                NL_types = {'lin','lin','lin','lin'};
                post_gsac_EImod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
                post_gsac_EImod.mods(1).reg_params = NMMcreate_reg_params();
                post_gsac_EImod.spk_NL_params = cur_rGQM.spk_NL_params;
                post_gsac_EImod = NMMfit_filters(post_gsac_EImod,cur_Robs(any_sac_inds),get_Xcell_tInds(tr_stim,any_sac_inds),[],[],silent);
                
                [post_gsac_EImod_LL,~,post_EImod_predrate] = NMMmodel_eval(post_gsac_EImod,cur_Robs,tr_stim);
                sacStimProc(cc).gsac_post_EImod = post_gsac_EImod;
                sacStimProc(cc).gsac_post_Egains = post_gsac_EImod.mods(3).filtK;
                sacStimProc(cc).gsac_post_Igains = post_gsac_EImod.mods(4).filtK;
            else
                sacStimProc(cc).gsac_post_EImod = nan;
                sacStimProc(cc).gsac_post_Egains = nan;
                sacStimProc(cc).gsac_post_Igains = nan;
            end
            
            [EI_xc,xc_lags] = xcov(g_exc,g_inh,flen,'coeff');
            sacStimProc(cc).EI_xc = EI_xc;
            sacStimProc(cc).EI_xc_lags = xc_lags;
            
            %now fit single post-gain model
            Xsac_tot = bsxfun(@times,cur_Xsac,stimG);
            clear tr_stim
            tr_stim{1} = [stimG];
            tr_stim{2} = cur_Xsac;
            tr_stim{3} = Xsac_tot;
            clear sac_stim_params
            sac_stim_params(1) = NMMcreate_stim_params(1);
            sac_stim_params(2:3) = NMMcreate_stim_params(size(Xsac_tot,2));
            mod_signs = [1 1 1];
            Xtargets = [1 2 3];
            NL_types = {'lin','lin','lin'};
            post_gsac_Smod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
            post_gsac_Smod.mods(1).reg_params = NMMcreate_reg_params();
            post_gsac_Smod.spk_NL_params = cur_rGQM.spk_NL_params;
            post_gsac_Smod = NMMfit_filters(post_gsac_Smod,cur_Robs(any_sac_inds),get_Xcell_tInds(tr_stim,any_sac_inds),[],[],silent);
            [post_gsac_Smod_LL,~,post_Smod_predrate,~,~,~,nullLL] = NMMmodel_eval(post_gsac_Smod,cur_Robs(any_sac_inds),get_Xcell_tInds(tr_stim,any_sac_inds));
            
            %overall info (evaluated on anysac inds)
            sacStimProc(cc).gsac_spost_ov_modinfo = mean(post_Smod_predrate/mean(post_Smod_predrate).*log2(post_Smod_predrate/mean(post_Smod_predrate)));
            
            [~,~,post_Smod_predrate] = NMMmodel_eval( post_gsac_Smod, cur_Robs, tr_stim);
            sacStimProc(cc).gsac_post_singmod = post_gsac_Smod;
            
            %now compute model with unique post-sac gain for each subunit
            %use pre-specified reg params for the full models
            sac_reg_params = NMMcreate_reg_params('lambda_d2T',fullMod_d2T,'lambda_L2',fullMod_L2,'boundary_conds',[0 0 0]);
            if fitFullPostMod
                Xsac_tot = bsxfun(@times,cur_Xsac,reshape(fgint,[],1,size(fgint,2)));
                clear Ftr_stim
                Ftr_stim{1} = [fgint];
                Ftr_stim{2} = cur_Xsac;
                Ftr_stim{3} = reshape(Xsac_tot,size(fgint,1),[]);
                clear sac_stim_params
                sac_stim_params(1) = NMMcreate_stim_params(size(fgint,2));
                sac_stim_params(2) = NMMcreate_stim_params(length(slags));
                sac_stim_params(3) = NMMcreate_stim_params([length(slags) size(fgint,2)]);
                mod_signs = [1 1 1];
                Xtargets = [1 2 3];
                NL_types = {'lin','lin','lin'};
                post_gsac_Fmod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
                post_gsac_Fmod.mods(1).reg_params = NMMcreate_reg_params();
                post_gsac_Fmod.spk_NL_params = cur_rGQM.spk_NL_params;
                post_gsac_Fmod = NMMfit_filters(post_gsac_Fmod,cur_Robs(any_sac_inds),get_Xcell_tInds(Ftr_stim,any_sac_inds),[],[],silent);
                [post_gsac_Fmod_LL,~,post_Fmod_predrate] = NMMmodel_eval(post_gsac_Fmod,cur_Robs(any_sac_inds),get_Xcell_tInds(Ftr_stim,any_sac_inds));
                
                %overall info (evaluated on anysac inds)
                sacStimProc(cc).gsac_Fpost_ov_modinfo = mean(post_Fmod_predrate/mean(post_Fmod_predrate).*log2(post_Fmod_predrate/mean(post_Fmod_predrate)));
                
                [~,~,post_Fmod_predrate] = NMMmodel_eval( post_gsac_Fmod, cur_Robs, Ftr_stim);
                sacStimProc(cc).gsac_post_Fmod = post_gsac_Fmod;
            end
 
            fullSacMod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,[],Xtargets);
            fullSacMod.spk_NL_params = cur_rGQM.spk_NL_params;
            [fullSacMod] = sacMod_scan_regularization(fullSacMod,cur_Robs,Ftr_stim,tr_sac_inds,xv_sac_inds,poss_gain_d2T,poss_gain_L2);
           

            %% FIT UPSTREAM STIM-MODULATION
            if fitUpstream
                fprintf('Fitting upstream saccade kernel\n');
                Xsac_mat = cur_Xsac(any_sac_inds,:);
                cur_tr_inds = find(ismember(any_sac_inds,tr_sac_inds));
                cur_xv_inds = find(ismember(any_sac_inds,xv_sac_inds));
                [preGainMod] = fit_pre_gainmodel(cur_rGQM,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:),Xsac_mat,poss_pre_d2T,opt_L2,cur_tr_inds,cur_xv_inds);
                sacStimProc(cc).gsacPreGainMod = preGainMod;
                [preLL,pre_pred_rate] = eval_pre_gainmodel( preGainMod, cur_Robs(any_sac_inds), all_Xmat_shift(any_sac_inds,:), cur_Xsac(any_sac_inds,:));
                sacStimProc(cc).gsac_pre_ov_modinfo = mean(pre_pred_rate/mean(pre_pred_rate).*log2(pre_pred_rate/mean(pre_pred_rate)));
                [~,pre_pred_rate] = eval_pre_gainmodel( preGainMod, cur_Robs, all_Xmat_shift, cur_Xsac);
             
                %                 [sacGainMod,sacGainOnlyMod] = fit_prepost_gainmodel(cur_rGQM,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:),Xsac_mat,opt_d2T,opt_L2,[],maxIter);
%                 sacStimProc(cc).gsacGainMod = sacGainMod;
%                 sacStimProc(cc).gsacGainOnlyMod = sacGainOnlyMod;
%                 [gainLL,gain_pred_rate] = eval_prepost_gainmodel( sacGainOnlyMod, cur_Robs(any_sac_inds), all_Xmat_shift(any_sac_inds,:), cur_Xsac(any_sac_inds,:));
%                 sacStimProc(cc).gsac_ov_modinfo = mean(gain_pred_rate/mean(gain_pred_rate).*log2(gain_pred_rate/mean(gain_pred_rate)));
%                 [~,gain_pred_rate] = eval_prepost_gainmodel( sacGainOnlyMod, cur_Robs, all_Xmat_shift, cur_Xsac);
            end
            %         sacGainMod = sacStimProc(cc).gsacGainMod;
            
            %% COMPUTE SUBSPAC SAC MODEL
            if fit_subMod
                fprintf('Estimating sac-dep subspace model\n');
                
                sub_Xsac = cur_Xsac; %create new Xsac mat that will not have any 'double sacs' present
                net_Xsac = sum(sub_Xsac,2); 
                multi_inds = find(net_Xsac > 1); %times when there were more than one sac present in the lag window
                %for these times drop the sac with the longer absolute lag
                for ii = 1:length(multi_inds)
                    cur_set = find(sub_Xsac(multi_inds(ii),:) == 1);
                    [~,larger] = min(abs(slags(cur_set)));
                    sub_Xsac(multi_inds(ii),cur_set(larger)) = 0;
                end
                
                X{1} = sub_Xsac;
                X{2} = reshape(bsxfun(@times,sub_Xsac,reshape(filt_outs,length(cc_uinds),1,[])),length(cc_uinds),[]);
                
                cur_stim_params(1) = NMMcreate_stim_params(length(slags));
                cur_stim_params(2) = NMMcreate_stim_params([length(slags) length(cur_rGQM.mods)]);
                
                %initialize model with same number (and type) of subunits, but
                %where the filter coefs can mix within the subspace
                mod_signs = [1 cur_rGQM.mods(:).sign];
                Xtargs = [1 2*ones(1,length(cur_rGQM.mods))];
                NL_types = cat(2,{'lin'},{cur_rGQM.mods(:).NLtype});
                modSeq_L2 = 0;
                
                %cycle over a range of possible d2T reg lambdas
                subspace_xvLL = nan(length(poss_sub_d2T),1);
                clear all_sub_mods
                for jj = 1:length(poss_sub_d2T)
                    modSeq_d2T = poss_sub_d2T(jj);
                    fprintf('Xval %d of %d\n',jj,length(poss_sub_d2T));
                    %use opt d2T lambda from gain/offset model for the offset
                    %kernel
                    reg_params = NMMcreate_reg_params('lambda_d2T',[opt_d2T repmat(modSeq_d2T,1,length(mod_signs)-1)]',...
                        'lambda_L2',[0 repmat(modSeq_L2,1,length(mod_signs)-1)]','boundary_conds',repmat([Inf Inf Inf],length(mod_signs),1));
                    init_mod = NMMinitialize_model(cur_stim_params,mod_signs,NL_types,reg_params,Xtargs);
                    init_mod.mods(1).reg_params.boundary_conds = [0 Inf Inf]; %0 boundary on offset term
                    
                    %set initial filters to have no mixing
                    for ii = 1:length(cur_rGQM.mods)
                        init_filt = zeros(length(slags),length(cur_rGQM.mods));
                        init_filt(:,ii) = 1;
                        init_mod.mods(ii+1).filtK = init_filt(:);
                    end
                    init_mod.spk_NL_params = cur_rGQM.spk_NL_params;
                    init_mod.mods(1).filtK = post_gsac_Smod.mods(2).filtK; %set initial offset filter
                    subspace_mod = NMMfit_filters(init_mod,cur_Robs(tr_sac_inds),get_Xcell_tInds(X,tr_sac_inds),[],[],silent);
                    subspace_xvLL(jj) = NMMmodel_eval(subspace_mod,cur_Robs(xv_sac_inds),get_Xcell_tInds(X,xv_sac_inds));
                    all_sub_mods(jj) = subspace_mod;
                end
                %select best lambda using xvLL
                [~,optloc] = max(subspace_xvLL);
                subspace_optL2 = poss_sub_d2T(optloc);
                sacStimProc(cc).gsac_subspace_optL2 = subspace_optL2;
                sacStimProc(cc).gsac_subspace_xvLL = subspace_xvLL;
                
                %now initialize model with optimal d2T lambda and fit to all
                %used data
                reg_params = NMMcreate_reg_params('lambda_d2T',[opt_d2T repmat(subspace_optL2,1,length(mod_signs)-1)]',...
                    'boundary_conds',repmat([Inf Inf Inf],length(mod_signs),1));
                init_mod = NMMinitialize_model(cur_stim_params,mod_signs,NL_types,reg_params,Xtargs);
                init_mod.spk_NL_params = cur_rGQM.spk_NL_params;
                init_mod.mods(1).reg_params.boundary_conds = [0 Inf Inf];
                for ii = 1:length(init_mod.mods)-1
                    init_filt = zeros(length(slags),length(cur_rGQM.mods));
                    init_filt(:,ii) = 1;
                    init_mod.mods(ii+1).filtK = init_filt(:);
                end
                init_mod.mods(1).filtK = post_gsac_Smod.mods(2).filtK;
                init_mod.spk_NL_params(1) = post_gsac_Smod.spk_NL_params(1);
                subspace_mod = NMMfit_filters(init_mod,cur_Robs(any_sac_inds),get_Xcell_tInds(X,any_sac_inds),[],[],silent);
                
                %evaluate overall info on anysac inds
                [subspace_LL,~,subspace_predrate,subG,subintG,sub_fgint] = NMMmodel_eval(subspace_mod,cur_Robs(any_sac_inds),get_Xcell_tInds(X,any_sac_inds));
                sacStimProc(cc).gsac_sub_ov_modinfo = mean(subspace_predrate/mean(subspace_predrate).*log2(subspace_predrate/mean(subspace_predrate)));
                
                sacStimProc(cc).gsac_submod = subspace_mod;
                
                [~,~,subspace_predrate] = NMMmodel_eval(subspace_mod,cur_Robs,X);
                
                %                 %extract the subspace filters
                %                 sub_efilts = find(Xtargs == 2 & mod_signs == 1);
                %                 cur_filts = reshape([subspace_mod.mods(sub_efilts).filtK],[length(slags) length(cur_rGQM.mods) length(sub_efilts)]);
                %                 stim_filts = reshape([cur_rGQM.mods.filtK],[flen*use_nPix_us length(cur_rGQM.mods)]);
                %                 sacdep_filts = nan(length(sub_efilts),length(slags),flen*use_nPix_us);
                %                 for jj = 1:length(sub_efilts)
                %                     sacdep_filts(jj,:,:) = squeeze(cur_filts(:,:,jj))*stim_filts';
                %                 end
                %                 sacStimProc(cc).gsac_phaseDep_subfilt = squeeze(sacdep_filts(1,:,:));
                %                 sacStimProc(cc).gsac_phaseInd_subfilt = squeeze(sqrt(sum(sacdep_filts(2:end,:,:).^2)));
            end
            
            %% CREATE TENT_BASIS MODEL OF SACCADE-MODULATION
            fprintf('Estimating tent-basis model\n');
            xbuff = 3; %add a lag-axis buffer to the model to minimize the impact of any boundary effects (with smoothing)
            
            Xtick = -(backlag+xbuff+1/2):(1):(forlag+xbuff+1/2);
            n_sbins = length(Xtick);
            
            %compute a single-valued 'time-since-saccade' parameter
            cur_sac_starts = saccade_start_inds(big_sacs);
            cur_sac_stops = saccade_stop_inds(big_sacs);
            t_since_sac_start = nan(NT,1);
            temp_cnt = zeros(NT,1);
            for ii = 1:length(cur_sac_starts)
                prev_tstart = find(trial_start_inds <= cur_sac_starts(ii),1,'last');
                next_tstop = find(trial_end_inds >= cur_sac_starts(ii),1,'first');
                cur_inds = (cur_sac_starts(ii) - backlag - xbuff):(cur_sac_starts(ii) + forlag + xbuff);
                cur_uset = find(cur_inds > trial_start_inds(prev_tstart) & cur_inds < trial_end_inds(next_tstop));
                %             t_since_sac_start(cur_inds(cur_uset)) = slags(cur_uset);
                t_since_sac_start(cur_inds(cur_uset)) = Xtick(cur_uset)+0.5;
                temp_cnt(cur_inds(cur_uset)) = temp_cnt(cur_inds(cur_uset)) + 1;
            end
            
            %initialize 2D TB data using t-since-sac and normalized G
            TB_stim = [t_since_sac_start(cc_uinds) norm_stimG];
            
            %set G-bins based on prctiles
            %         Ytick = linspace(my_prctile(TB_stim(any_sac_inds,2),1),my_prctile(TB_stim(any_sac_inds,2),99),n_Gbins);
            Ytick = my_prctile(TB_stim(any_sac_inds,2),linspace(0.5,99.5,n_Gbins));
            
            %in some cases there are an excess of G==0 values, which makes
            %some adjacent bins identical. Need to cut these bin edges.
            nd_bins = find(diff(Ytick) <= 0);
            Ytick(nd_bins) = [];
            cur_nGbins = length(Ytick);
            
            %initialize TBs
            TB = TentBasis2D(Xtick, Ytick);
            
            %this is data within range of the TB centers
            TB_used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
                TB_stim(:,2) >= Ytick(1) & TB_stim(:,2) <= Ytick(end));
            
            double_sac_inds = find(temp_cnt(cc_uinds) > 1);
            TB_used_data(ismember(TB_used_data,double_sac_inds)) = [];
            
            udata_tr = find(ismember(TB_used_data,cur_tr_inds));
            udata_xv = find(ismember(TB_used_data,cur_xv_inds));
            
            %process data with TBs
            [TB_Xmat,TB_counts] = TB.InputNL2D(TB_stim(TB_used_data,:));
            
            TB_xvLL = nan(1,length(poss_TB_lambdas));
            for ll = 1:length(poss_TB_lambdas)
                fprintf('Fitting TB gsac model lambda %d/%d\n',ll,length(poss_TB_lambdas));
                cur_sac_lambda = poss_TB_lambdas(ll);
                
                %fit a penalized GLM on the TB outputs
                %            L2_params = create_L2_params([],[1 n_sbins*n_Gbins],[n_sbins n_Gbins],2,3,[Inf Inf],[TB_relsac_pen 1]);
                %this gives fixed reg strength in G-direction, but
                %cur_sac_lambda in T-direction
                L2_params = create_L2_params([],[1 n_sbins*cur_nGbins],[n_sbins cur_nGbins],2,3,[Inf Inf],[1 1/cur_sac_lambda*G_lambdas]);
                TB_fitmod = regGLM_fit(TB_Xmat(udata_tr,:),cur_Robs(TB_used_data(udata_tr)),L2_params,cur_sac_lambda,[],[],silent);
                [TB_xvLL(ll)] = regGLM_eval(TB_fitmod,cur_Robs(TB_used_data(udata_xv)),TB_Xmat(udata_xv,:));
            end
            
            %select optimal regularization lambda
            [~,optloc] = max(TB_xvLL);
            TB_optL2 = poss_TB_lambdas(optloc);
            sacStimProc(cc).gsac_TB_optL2 = TB_optL2;
            sacStimProc(cc).gsac_TB_xvLL = TB_xvLL;
            
            %fit a penalized GLM on the TB outputs using optimal lambda
            %         L2_params = create_L2_params([],[1 n_sbins*n_Gbins],[n_sbins n_Gbins],2,3,[Inf Inf],[TB_relsac_pen 1]);
            L2_params = create_L2_params([],[1 n_sbins*cur_nGbins],[n_sbins cur_nGbins],2,3,[Inf Inf],[1 1/TB_optL2*G_lambdas]);
            TB_fitmod = regGLM_fit(TB_Xmat,cur_Robs(TB_used_data),L2_params,TB_optL2,[],[],silent);
            [LL, penLL, TB_pred_rate, G] = regGLM_eval(TB_fitmod,cur_Robs(TB_used_data),TB_Xmat);
            cur = mean(cur_Robs(TB_used_data));
            TB_nullLL = nansum(cur_Robs(TB_used_data).*log2(cur) - cur)/sum(cur_Robs(TB_used_data));
            TB_LL = nansum(cur_Robs(TB_used_data).*log2(TB_pred_rate) - TB_pred_rate)/sum(cur_Robs(TB_used_data));
            
            %compute output of TB model
            TB_K = reshape(TB_fitmod.K,n_sbins,cur_nGbins)';
            bin_areas = TB.GetBinAreas();
            gsac_TB_dist = TB_counts./bin_areas;
            gsac_TB_dist = gsac_TB_dist'/sum(gsac_TB_dist(:));
            gsac_TB_rate = log(1 + exp(TB_K + TB_fitmod.theta));
            sacStimProc(cc).gsac_TB_rate = gsac_TB_rate(:,xbuff+1:end-xbuff);
            
            %INFO CALS
            cur_avg_rate = mean(cur_Robs(TB_used_data));
            marg_gdist = trapz(gsac_TB_dist,2);
            marg_sdist = trapz(gsac_TB_dist);
            marg_gsacrate = trapz(gsac_TB_dist.*gsac_TB_rate)./marg_sdist;
            marg_grate = trapz(gsac_TB_dist.*gsac_TB_rate,2)./marg_gdist;
            gsacdep_info = nan(1,n_sac_bins);
            for tt = 1:n_sbins
                gsacdep_info(tt) = trapz(gsac_TB_dist(:,tt).*gsac_TB_rate(:,tt)/marg_gsacrate(tt).*log2(gsac_TB_rate(:,tt)/marg_gsacrate(tt)))/trapz(gsac_TB_dist(:,tt));
            end
            
%             sacStimProc(cc).gsac_ov_TB_info = sum(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate))/cur_avg_rate;
            sacStimProc(cc).gsac_ov_TB_info = trapz(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate))/cur_avg_rate;
            
            sacStimProc(cc).gsac_TB_avg_rate = marg_gsacrate(xbuff+1:end-xbuff);
            sacStimProc(cc).gsac_TB_info = gsacdep_info(xbuff+1:end-xbuff);
            sacStimProc(cc).gsac_TB_gdist = marg_gdist;
            sacStimProc(cc).gsac_TB_grate = marg_grate;
            sacStimProc(cc).gsac_TB_lagX = Xtick(xbuff+1:end-xbuff);
            sacStimProc(cc).gsac_TB_gX = Ytick;
            
            dYtick = median(diff(Ytick));
            equi_space_gX = linspace(Ytick(1)-dYtick,Ytick(end)+dYtick,51);
            equi_space_gdist = histc(TB_stim(TB_used_data,2),equi_space_gX);
            sacStimProc(cc).gsac_equi_space_gdist = equi_space_gdist(1:end-1)/sum(equi_space_gdist);
            sacStimProc(cc).gsac_equi_space_gX = 0.5*equi_space_gX(1:end-1) + 0.5*equi_space_gX(2:end);
            
            [sac_offset,sac_gain] = deal(nan(length(sacStimProc(cc).gsac_TB_lagX),1));
            for ii = 1:length(sacStimProc(cc).gsac_TB_lagX)
               temp = lscov([ones(cur_nGbins,1) marg_grate],sacStimProc(cc).gsac_TB_rate(:,ii),marg_gdist);
               sac_offset(ii) = temp(1); sac_gain(ii) = temp(2);
            end
            sacStimProc(cc).gsac_TB_gain = sac_gain;
            sacStimProc(cc).gsac_TB_offset = sac_offset;
            
            
            %% COMPUTE MODEL-BASED INFORMATION
            [sac_spost_info,sac_fpost_info,sac_subpost_info,sac_pre_info,sac_spost_offset,sac_spost_gain,sac_TB_modinfo] = deal(nan(length(slags),1));
            
            for ii = 1:length(slags)
                temp = find(cur_Xsac(:,ii) == 1);
                
                rr = regress(post_Smod_predrate(temp),[ones(length(temp),1) basemod_pred_rate(temp)]);
                sac_spost_offset(ii) = rr(1);
                sac_spost_gain(ii) = rr(2);
                
                %compute LL and info for upstream model
                if fitUpstream
                    sac_pre_info(ii) = nanmean(pre_pred_rate(temp).*log2(pre_pred_rate(temp)/mean(pre_pred_rate(temp))))/mean(pre_pred_rate(temp));
                end
                
                %compute LL and info for gain/offset model
                sac_spost_info(ii) = nanmean(post_Smod_predrate(temp).*log2(post_Smod_predrate(temp)/mean(post_Smod_predrate(temp))))/mean(post_Smod_predrate(temp));
                
                if fitFullPostMod
                    sac_fpost_info(ii) = nanmean(post_Fmod_predrate(temp).*log2(post_Fmod_predrate(temp)/mean(post_Fmod_predrate(temp))))/mean(post_Fmod_predrate(temp));
                end
                
                temp2 = find(ismember(TB_used_data,temp));
                sac_TB_modinfo(ii) = nanmean(TB_pred_rate(temp2).*log2(TB_pred_rate(temp2)/mean(TB_pred_rate(temp2))))/mean(TB_pred_rate(temp2));

                %compute LL and info for subpsace model
                if fit_subMod
                    temp = find(sub_Xsac(:,ii) == 1);
                    sac_subpost_info(ii) = nanmean(subspace_predrate(temp).*log2(subspace_predrate(temp)/mean(subspace_predrate(temp))))/mean(subspace_predrate(temp));
                end
            end
            sacStimProc(cc).gsac_spost_modinfo = sac_spost_info;
            sacStimProc(cc).gsac_spost_offset = sac_spost_offset;
            sacStimProc(cc).gsac_spost_gain = sac_spost_gain;
            
            sacStimProc(cc).gsac_TB_modinfo = sac_TB_modinfo;
            if fitFullPostMod
                sacStimProc(cc).gsac_fpost_modinfo = sac_fpost_info;
            end
            
            if fitUpstream
                sacStimProc(cc).gsac_pre_modinfo = sac_pre_info;
            end
            if fit_subMod
                sacStimProc(cc).gsac_sub_modinfo = sac_subpost_info;
            end
            
            %% Get LL-based infos
            [sac_pre_LL,sac_fpost_LL,sac_spost_LL,sac_subpost_LL,sac_TB_LL,sac_nullLL,sac_Nspks,sac_avgrate] = deal(nan(length(slags),1));
            
            norm_stimE = nanzscore(g_exc);
            norm_stimI = nanzscore(g_inh);
            
            spk_cond_G = nan(length(slags),1);
            spk_cond_E = nan(length(slags),1);
            spk_cond_I = nan(length(slags),1);
            for ii = 1:length(slags)
                temp = find(cur_Xsac(:,ii) == 1);
                sac_avgrate(ii) = mean(cur_Robs(temp));
                
                %compute nullLL for data at this latency
                sac_nullLL(ii) = nansum(cur_Robs(temp).*log2(sac_avgrate(ii)) - sac_avgrate(ii));
                sac_Nspks(ii) = sum(cur_Robs(temp));
                
                %compute LL and info for upstream model
                if fitUpstream
                    sac_pre_LL(ii) = nansum(cur_Robs(temp).*log2(pre_pred_rate(temp)) - pre_pred_rate(temp));
                end
                
                %compute LL and info for gain/offset model
                sac_spost_LL(ii) = nansum(cur_Robs(temp).*log2(post_Smod_predrate(temp)) - post_Smod_predrate(temp));

                if fitFullPostMod
                sac_fpost_LL(ii) = nansum(cur_Robs(temp).*log2(post_Fmod_predrate(temp)) - post_Fmod_predrate(temp));                
                end
                %spike-weighted average of G (normalized)
                spk_cond_G(ii) = sum(cur_Robs(temp).*norm_stimG(temp))/sum(cur_Robs(temp));
                spk_cond_E(ii) = sum(cur_Robs(temp).*norm_stimE(temp))/sum(cur_Robs(temp));
                spk_cond_I(ii) = sum(cur_Robs(temp).*norm_stimI(temp))/sum(cur_Robs(temp));
                                
                %compute LL and info for subpsace model
                if fit_subMod
                    subtemp = temp(sub_Xsac(temp,ii) == 1);
                    sac_subpost_LL(ii) = nansum(cur_Robs(subtemp).*log2(subspace_predrate(subtemp)) - subspace_predrate(subtemp));
                    cur = mean(cur_Robs(subtemp));
                    cur_nullLL = nansum(cur_Robs(subtemp).*log2(cur) - cur);
                    cur_nspk = nansum(cur_Robs(subtemp));
                    sac_subpost_LL(ii) = (sac_subpost_LL(ii) - cur_nullLL)/cur_nspk;
                end
                
                [lia,lib] = ismember(temp,TB_used_data);
                sac_TB_LL(ii) = nansum(cur_Robs(temp(lia)).*log2(TB_pred_rate(lib(lia))) - TB_pred_rate(lib(lia)));
                cur = mean(cur_Robs(temp(lia)));
                cur_nullLL = nansum(cur_Robs(temp(lia)).*log2(cur) - cur);
                cur_nspks = nansum(cur_Robs(temp(lia)));
                sac_TB_LL(ii) = (sac_TB_LL(ii) - cur_nullLL)/cur_nspks;
            end
            
            %store spike-weighted G data
            sacStimProc(cc).gsac_spkCondG = spk_cond_G;
            sacStimProc(cc).gsac_spkCondE = spk_cond_E;
            sacStimProc(cc).gsac_spkCondI = spk_cond_I;
            ov_spkCondG = sum(cur_Robs(any_sac_inds).*norm_stimG(any_sac_inds))/sum(cur_Robs(any_sac_inds));
            ov_spkCondE = sum(cur_Robs(any_sac_inds).*norm_stimE(any_sac_inds))/sum(cur_Robs(any_sac_inds));
            ov_spkCondI = sum(cur_Robs(any_sac_inds).*norm_stimI(any_sac_inds))/sum(cur_Robs(any_sac_inds));
            sacStimProc(cc).gsac_ovspkCondG = ov_spkCondG;
            sacStimProc(cc).gsac_ovspkCondE = ov_spkCondE;
            sacStimProc(cc).gsac_ovspkCondI = ov_spkCondI;
            
            sacStimProc(cc).gsac_avg_rate = sac_avgrate;
            
            %store gain/offset model info calcs
            sacStimProc(cc).gsac_spost_LLinfo = (sac_spost_LL - sac_nullLL)./sac_Nspks;
            sacStimProc(cc).gsac_spost_ov_LLinfo = (post_gsac_Smod_LL-nullLL)/log(2);
            if fitFullPostMod
            sacStimProc(cc).gsac_fpost_LLinfo = (sac_fpost_LL - sac_nullLL)./sac_Nspks;
            sacStimProc(cc).gsac_fpost_ov_LLinfo = (post_gsac_Fmod_LL-nullLL)/log(2);
            end
            %store upstream model info calcs
            if fitUpstream
                sacStimProc(cc).gsac_pre_LLinfo = (sac_pre_LL - sac_nullLL)./sac_Nspks;
                sacStimProc(cc).gsac_pre_ov_LLinfo = (preGainMod.LL-nullLL)/log(2);
            end
            
            sacStimProc(cc).gsac_TB_LLinfo = sac_TB_LL;
            sacStimProc(cc).gsac_TB_ov_LLinfo = TB_LL-TB_nullLL;
            
            %store subspace model info calcs
            if fit_subMod
                sacStimProc(cc).gsac_sub_LLinfo = sac_subpost_LL;
                sacStimProc(cc).gsac_sub_ov_LLinfo = (subspace_LL-nullLL)/log(2);
            end
            
            %% COMPUTE SAC-CONDITIONAL STAS
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
                
                sacStimProc(cc).ov_phaseDep_sta = reshape(sum(bsxfun(@times,all_Xmat_shift,cur_Robs))/sum(cur_Robs)-mean(all_Xmat_shift),flen,[]);
                sacStimProc(cc).ov_phaseInd_sta = reshape(sum(bsxfun(@times,abs(all_Xmat_shift),cur_Robs))/sum(cur_Robs)-mean(abs(all_Xmat_shift)),flen,[]);
            end
            
            
        end
        %% FOR MSACS
        if fitMsacs
            cur_Xsac = Xmsac(cc_uinds,:); %saccade indicator Xmat
            
            %only use indices during guided saccade expts here
            any_sac_inds = find(any(cur_Xsac > 0,2));
            tr_sac_inds = cur_tr_inds(ismember(cur_tr_inds,any_sac_inds));
            xv_sac_inds = cur_xv_inds(ismember(cur_xv_inds,any_sac_inds));
            
            %%
            if length(any_sac_inds) > 1e4
                %% Fit spk NL params and refit scale of each filter using target data (within trange of sacs)
                if fit_unCor
                cur_rGQM = ModData(cc).rectGQM_unCor;    
                else
                cur_rGQM = ModData(cc).rectGQM;
                end
                
                cur_rGQM = NMMfit_logexp_spkNL(cur_rGQM,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:));
                cur_rGQM = NMMfit_scale(cur_rGQM,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:));
                
                stim_mod_signs = [cur_rGQM.mods(:).sign];
                [~,~,basemod_pred_rate,~,filt_outs,fgint] = NMMmodel_eval(cur_rGQM,cur_Robs,all_Xmat_shift);
                fgint = bsxfun(@times,fgint,stim_mod_signs);
                stimG = sum(fgint,2);
                norm_stimG = zscore(stimG);
                
                sacStimProc(cc).msac_ovavg_rate = mean(cur_Robs(any_sac_inds));
                %% FOR SIMPLE POST_GAIN MODEL, SCAN RANGE OF L2s AND SELECT BEST USING XVAL LL
                
%                 Xsac_tot = bsxfun(@times,cur_Xsac,stimG);
%                 clear tr_stim
%                 tr_stim{1} = [stimG];
%                 tr_stim{2} = cur_Xsac;
%                 tr_stim{3} = Xsac_tot;
%                 clear sac_stim_params
%                 sac_stim_params(1) = NMMcreate_stim_params(1);
%                 sac_stim_params(2:3) = NMMcreate_stim_params([size(Xsac_tot,2)]);
%                 mod_signs = [1 1 1];
%                 Xtargets = [1 2 3];
%                 NL_types = {'lin','lin','lin'};
%                 
%                 L2_gain_xvLL = nan(length(poss_gain_d2T),length(poss_gain_L2));
%                 L2_gain_LL = nan(length(poss_gain_d2T),length(poss_gain_L2));
%                 for jj = 1:length(poss_gain_d2T)
%                     for ii = 1:length(poss_gain_L2)
%                         sac_reg_params = NMMcreate_reg_params('lambda_d2T',poss_gain_d2T(jj),'lambda_L2',poss_gain_L2(ii),'boundary_conds',[0 0 0]);
%                         cur_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
%                         cur_mod.mods(1).reg_params = NMMcreate_reg_params();
%                         cur_mod.spk_NL_params = cur_rGQM.spk_NL_params;
%                         cur_mod = NMMfit_filters(cur_mod,cur_Robs(tr_sac_inds),get_Xcell_tInds(tr_stim,tr_sac_inds),[],[],silent);
%                         L2_gain_xvLL(jj,ii) = NMMmodel_eval(cur_mod,cur_Robs(xv_sac_inds),get_Xcell_tInds(tr_stim,xv_sac_inds));
%                         L2_gain_LL(jj,ii) = NMMmodel_eval(cur_mod,cur_Robs(tr_sac_inds),get_Xcell_tInds(tr_stim,tr_sac_inds));
%                     end
%                 end
%                 
%                 sacStimProc(cc).msac_spost_xvLL = L2_gain_xvLL;
%                 [~,optloc] = max(L2_gain_xvLL(:));
%                 [optloc_x,optloc_y] = ind2sub([length(poss_gain_d2T) length(poss_gain_L2)],optloc);
%                 opt_d2T = poss_gain_d2T(optloc_x);
%                 opt_L2 = poss_gain_L2(optloc_y);
%                 sacStimProc(cc).msac_optd2T = opt_d2T;
%                 sacStimProc(cc).msac_optL2 = opt_L2;
                
                %% FIT POST-INTEGRATION GAIN
                fprintf('Fitting post-filter models\n');
                sac_reg_params = NMMcreate_reg_params('lambda_d2T',opt_d2T,'lambda_L2',opt_L2,'boundary_conds',[0 0 0]);
                
%                 %if there are both significant E and I filters, fit model with sep
%                 %E and I gain kernels
%                 if sum(stim_mod_signs == 1) > 0 && sum(stim_mod_signs == -1) > 0
%                     g_exc = sum(fgint(:,stim_mod_signs==1),2);
%                     g_inh = sum(fgint(:,stim_mod_signs==-1),2);
%                     Xsac_estim = bsxfun(@times,cur_Xsac,g_exc);
%                     Xsac_istim = bsxfun(@times,cur_Xsac,g_inh);
%                     clear tr_stim
%                     tr_stim{1} = [g_exc g_inh];
%                     tr_stim{2} = cur_Xsac;
%                     tr_stim{3} = Xsac_estim;
%                     tr_stim{4} = Xsac_istim;
%                     sac_stim_params(1) = NMMcreate_stim_params(2);
%                     sac_stim_params(2:4) = NMMcreate_stim_params(size(cur_Xsac,2));
%                     mod_signs = [1 1 1 1];
%                     Xtargets = [1 2 3 4];
%                     NL_types = {'lin','lin','lin','lin'};
%                     post_msac_EImod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
%                     post_msac_EImod.mods(1).reg_params = NMMcreate_reg_params();
%                     post_msac_EImod.spk_NL_params = cur_rGQM.spk_NL_params;
%                     post_msac_EImod = NMMfit_filters(post_msac_EImod,cur_Robs(any_sac_inds),get_Xcell_tInds(tr_stim,any_sac_inds),[],[],silent);
%                     
%                     [post_msac_EImod_LL,~,post_EImod_predrate] = NMMmodel_eval(post_msac_EImod,cur_Robs,tr_stim);
%                     sacStimProc(cc).msac_post_EImod = post_msac_EImod;
%                     sacStimProc(cc).msac_post_Egains = post_msac_EImod.mods(3).filtK;
%                     sacStimProc(cc).msac_post_Igains = post_msac_EImod.mods(4).filtK;
%                 else
%                     sacStimProc(cc).msac_post_EImod = nan;
%                     sacStimProc(cc).msac_post_Egains = nan;
%                     sacStimProc(cc).msac_post_Igains = nan;
%                 end
                
                Xsac_tot = bsxfun(@times,cur_Xsac,stimG);
                clear tr_stim
                tr_stim{1} = [stimG];
                tr_stim{2} = cur_Xsac;
                tr_stim{3} = Xsac_tot;
                clear sac_stim_params
                sac_stim_params(1) = NMMcreate_stim_params(1);
                sac_stim_params(2:3) = NMMcreate_stim_params(size(Xsac_tot,2));
                mod_signs = [1 1 1];
                Xtargets = [1 2 3];
                NL_types = {'lin','lin','lin'};
                post_msac_Smod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
                post_msac_Smod.mods(1).reg_params = NMMcreate_reg_params();
                post_msac_Smod.spk_NL_params = cur_rGQM.spk_NL_params;
                post_msac_Smod = NMMfit_filters(post_msac_Smod,cur_Robs(any_sac_inds),get_Xcell_tInds(tr_stim,any_sac_inds),[],[],silent);
                [post_msac_Smod_LL,~,post_Smod_predrate,~,~,~,nullLL] = NMMmodel_eval(post_msac_Smod,cur_Robs(any_sac_inds),get_Xcell_tInds(tr_stim,any_sac_inds));
                
                %overall info (evaluated on anysac inds)
                sacStimProc(cc).msac_spost_ov_modinfo = mean(post_Smod_predrate/mean(post_Smod_predrate).*log2(post_Smod_predrate/mean(post_Smod_predrate)));
                
                [~,~,post_Smod_predrate] = NMMmodel_eval( post_msac_Smod, cur_Robs, tr_stim);
                sacStimProc(cc).msac_post_singmod = post_msac_Smod;
                
                %% FIT UPSTREAM STIM-MODULATION
                if fit_msacUpstream
                    fprintf('Fitting upstream saccade kernel\n');
                    Xsac_mat = cur_Xsac(any_sac_inds,:);
%                     cur_poss_d2T = sacStimProc(cc).gsacPreGainMod.opt_d2T;
%                     cur_poss_L2 = sacStimProc(cc).gsacPreGainMod.opt_L2;
                    cur_tr_inds = find(ismember(any_sac_inds,tr_sac_inds));
                    cur_xv_inds = find(ismember(any_sac_inds,xv_sac_inds));
                    [preGainMod] = fit_pre_gainmodel(cur_rGQM,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:),Xsac_mat,poss_pre_d2T,opt_L2,cur_tr_inds,cur_xv_inds);
                    sacStimProc(cc).msacPreGainMod = preGainMod;
                    [preLL,pre_pred_rate] = eval_pre_gainmodel( preGainMod, cur_Robs(any_sac_inds), all_Xmat_shift(any_sac_inds,:), cur_Xsac(any_sac_inds,:));
                    sacStimProc(cc).msac_pre_ov_modinfo = mean(pre_pred_rate/mean(pre_pred_rate).*log2(pre_pred_rate/mean(pre_pred_rate)));
                    [~,pre_pred_rate] = eval_pre_gainmodel( preGainMod, cur_Robs, all_Xmat_shift, cur_Xsac);
                    
                    %                     Xsac_mat = cur_Xsac(any_sac_inds,:);
%                     [sacGainMod,sacGainOnlyMod] = fit_prepost_gainmodel(cur_rGQM,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:),Xsac_mat,opt_d2T,opt_L2,[],maxIter);
%                     sacStimProc(cc).msacGainMod = sacGainMod;
%                     sacStimProc(cc).msacGainOnlyMod = sacGainOnlyMod;
%                     [gainLL,gain_pred_rate] = eval_prepost_gainmodel( sacGainOnlyMod, cur_Robs(any_sac_inds), all_Xmat_shift(any_sac_inds,:), cur_Xsac(any_sac_inds,:));
%                     sacStimProc(cc).msac_ov_modinfo = mean(gain_pred_rate/mean(gain_pred_rate).*log2(gain_pred_rate/mean(gain_pred_rate)));
%                     [~,gain_pred_rate] = eval_prepost_gainmodel( sacGainOnlyMod, cur_Robs, all_Xmat_shift, cur_Xsac);
                end
                
                %% COMPUTE MODEL-BASED INFORMATION
                [sac_spost_info,sac_pre_info,sac_spost_offset,sac_spost_gain] = deal(nan(length(slags),1));
                for ii = 1:length(slags)
                    temp = find(cur_Xsac(:,ii) == 1);
                    
                    rr = regress(post_Smod_predrate(temp),[ones(length(temp),1) basemod_pred_rate(temp)]);
                    sac_spost_offset(ii) = rr(1);
                    sac_spost_gain(ii) = rr(2);
                    
                    %compute LL and info for upstream model
                    if fit_msacUpstream
                        sac_pre_info(ii) = nanmean(pre_pred_rate(temp).*log2(pre_pred_rate(temp)/mean(pre_pred_rate(temp))))/mean(pre_pred_rate(temp));
                    end
                    
                    %compute LL and info for gain/offset model
                    sac_spost_info(ii) = nanmean(post_Smod_predrate(temp).*log2(post_Smod_predrate(temp)/mean(post_Smod_predrate(temp))))/mean(post_Smod_predrate(temp));
                end
                
                sacStimProc(cc).msac_spost_modinfo = sac_spost_info;
                sacStimProc(cc).msac_spost_offset = sac_spost_offset;
                sacStimProc(cc).msac_spost_gain = sac_spost_gain;
                if fit_msacUpstream
                    sacStimProc(cc).msac_pre_modinfo = sac_pre_info;
                end
                
                %% Get LL-based infos
                [sac_pre_LL,sac_fpost_LL,sac_spost_LL,sac_nullLL,sac_Nspks,sac_avgrate] = deal(nan(length(slags),1));
                
                spk_cond_G = nan(length(slags),1);
                spk_cond_lagG = nan(length(slags),flen);
                for ii = 1:length(slags)
                    temp = find(cur_Xsac(:,ii) == 1);
                    sac_avgrate(ii) = mean(cur_Robs(temp));
                    
                    %compute nullLL for data at this latency
                    sac_nullLL(ii) = nansum(cur_Robs(temp).*log2(sac_avgrate(ii)) - sac_avgrate(ii));
                    sac_Nspks(ii) = sum(cur_Robs(temp));
                    
                    %compute LL and info for upstream model
                    if fit_msacUpstream
                        sac_pre_LL(ii) = nansum(cur_Robs(temp).*log2(pre_pred_rate(temp)) - pre_pred_rate(temp));
                    end
                    
                    %compute LL and info for gain/offset model
                    sac_spost_LL(ii) = nansum(cur_Robs(temp).*log2(post_Smod_predrate(temp)) - post_Smod_predrate(temp));
                    
                    %spike-weighted average of G (normalized)
                    spk_cond_G(ii) = sum(cur_Robs(temp).*norm_stimG(temp))/sum(cur_Robs(temp));
                    
                end
                
                sacStimProc(cc).msac_spkCondG = spk_cond_G;
                ov_spkCondG = sum(cur_Robs(any_sac_inds).*norm_stimG(any_sac_inds))/sum(cur_Robs(any_sac_inds));
                sacStimProc(cc).msac_ovspkCondG = ov_spkCondG;
                                
                sacStimProc(cc).msac_avg_rate = sac_avgrate;
                
                %store gain/offset model info calcs
                sacStimProc(cc).msac_spost_LLinfo = (sac_spost_LL - sac_nullLL)./sac_Nspks;
                sacStimProc(cc).msac_spost_ov_LLinfo = (post_msac_Smod_LL-nullLL)/log(2);
                
                %store upstream model info calcs
                if fit_msacUpstream
                    sacStimProc(cc).msac_pre_LLinfo = (sac_pre_LL - sac_nullLL)./sac_Nspks;
                    sacStimProc(cc).msac_pre_ov_LLinfo = (preGainMod.LL-nullLL)/log(2);
                end
                
                %% COMPUTE SAC-CONDITIONAL STAS
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
                    sacStimProc(cc).msac_phaseDep_sta = reshape(lag_dep_sta,length(slags),flen, use_nPix_us);
                    sacStimProc(cc).msac_phaseInd_sta = reshape(lag_dep_asta,length(slags),flen,use_nPix_us);
                    
                    sacStimProc(cc).ov_phaseDep_sta = reshape(sum(bsxfun(@times,all_Xmat_shift,cur_Robs))/sum(cur_Robs)-mean(all_Xmat_shift),flen,[]);
                    sacStimProc(cc).ov_phaseInd_sta = reshape(sum(bsxfun(@times,abs(all_Xmat_shift),cur_Robs))/sum(cur_Robs)-mean(abs(all_Xmat_shift)),flen,[]);
                end
                
            end
        end
    else
        sacStimProc(cc).used = false;
    end
end

%%
anal_dir = ['/home/james/Analysis/bruce/' Expt_name '/FINsac_mod/'];

sname = [sname sprintf('_ori%d',bar_ori)];
if fit_unCor
    sname = [sname '_unCor'];
end
cd(anal_dir)
save(sname,'targs','slags','dt','sacStimProc');

