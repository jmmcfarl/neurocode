
%
% clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements//');
addpath('~/James_scripts/bruce/processing/');
addpath('~/James_scripts/bruce/saccade_modulation/');
addpath('~/James_scripts/bruce/sacModFinal/');
addpath('~/James_scripts/TentBasis2D/');

global Expt_name bar_ori use_MUA

% Expt_name = 'G093';
% use_MUA = false;
% bar_ori = 0; %bar orientation to use (only for UA recs)

fit_unCor = false;
fit_subMod = true;
fitUpstream = true;
fitSTA = true;
fitMsacs = true;
fitFullPostMod = true;

include_bursts = 0;

% sname = 'sacStimDelay_noXV';
sname = 'sacStimDelay_noXV2';
if include_bursts
    sname = [sname '_withbursts'];
end

mod_data_name = 'corrected_models2';
base_sname = 'sacStimProcFin_noXV'; %base datat file name for full sacmod analysis

%%
% poss_gain_d2T = logspace(log10(1),log10(1e3),8); %range of d2T reg values for post-gain models
% poss_full_gain_d2T = 10;
% poss_gain_L2 = [0 logspace(log10(1),log10(50),4)]; %range of L2 reg values
% poss_pre_d2T = logspace(log10(1),log10(1e3),8); %range of d2T reg values for pre-gain models
% poss_sub_d2T = logspace(log10(10),log10(1e3),4); %range of d2T reg values for subspace models
% poss_TB_lambdas = logspace(log10(0.1),log10(500),8); %range of d2T reg values for TB models

poss_off_d2T = [1 10 50 100 1000];
cent_off_d2T = 100;
poss_gain_d2T = [1 10 50 100 1000]; %range of d2T reg values for post-gain models
poss_full_gain_d2T = [1 10 50 100 1000];
poss_gain_L2 = [0 1 5 10]; %range of L2 reg values
poss_pre_d2T = [1 10 50 100 1000]; %range of d2T reg values for pre-gain models
poss_sub_d2T = [1 10 50 100 1000]; %range of d2T reg values for subspace models
poss_TB_lambdas = [0.1 1 5 10 100]; %range of d2T reg values for TB models

dt = 0.01;
backlag = round(0.1/dt);
forlag = round(0.3/dt);
slags = -backlag:forlag;
n_sac_bins = length(slags);

TB_params.xbuff = 3;
TB_params.n_Gbins = 35; %number of bins for TB model
TB_params.G_lambdas = 100; %d2G reg parameter for TB model (smoothness in g-dimension)
TB_params.backlag = backlag;
TB_params.forlag = forlag;

micro_thresh = 1; %max amp of microsac (deg)
EP_bounds = 1;%eye position boundary (deg from central FP)
sac_burst_isi = 0.15; %minimum inter-saccade interval for micros (to eliminate 'bursts')
max_gsac_dur = 0.1; %maximum saccade duration before we call it a likely blink

xv_frac = 0.2; %fraction of trials to use for cross-validation

post_lambda_off = 4; %offset d2T lambda for post-model
pre_lambda = 3; %this is just the gain kernel lambda, the offset kernel lambda is fixed equal to post-model (4)
post_lambda_gain = 3; %gain d2T for post model

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
% et_anal_name = 'full_eyetrack_Rinit_Cprior_ori0_old2.mat';
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

%% LOAD IN FULL SACMOD DATA STRUCTURE
cd(anal_dir)
data_name = [base_sname sprintf('_ori%d',bar_ori)];
load(data_name,'sacStimProc')

%% MAIN ANALYSIS LOOP
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
    
    sacDelay(cc).xv_trials = xv_trials;
    sacDelay(cc).tr_trials = tr_trials;
    
    base_tr_inds = find(ismember(all_trialvec(used_inds(cc_uinds)),tr_trials));
    base_xv_inds = find(ismember(all_trialvec(used_inds(cc_uinds)),xv_trials));
    
    cur_Robs = Robs_mat(cc_uinds,cc);
    
    %%
    if ~isempty(cc_uinds)
        
        if fit_unCor %if not using ET corrections
            cur_rGQM = ModData(cc).rectGQM_unCor;
        else
            cur_rGQM = ModData(cc).rectGQM;
        end
        sacDelay(cc).ModData = ModData(cc);
        sacDelay(cc).used = true;
        
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
        tr_sac_inds = base_tr_inds(ismember(base_tr_inds,any_sac_inds));
        xv_sac_inds = base_xv_inds(ismember(base_xv_inds,any_sac_inds));
        
        %%
        if length(any_sac_inds) > 1e4
            
            %% Fit spk NL params and refit scale of each filter using target data (within trange of sacs)
            cur_rGQM = NMMfit_logexp_spkNL(cur_rGQM,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:));
            cur_rGQM = NMMfit_scale(cur_rGQM,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:));
            
            stim_mod_signs = [cur_rGQM.mods(:).sign];
            [~,~,basemod_pred_rate,~,filt_outs,fgint] = NMMmodel_eval(cur_rGQM,cur_Robs,all_Xmat_shift);
            fgint = bsxfun(@times,fgint,stim_mod_signs);
            
            sacDelay(cc).gsac_ovavg_rate = mean(cur_Robs(any_sac_inds));
            [baseLL,~,~,~,~,~,base_nullLL] = NMMmodel_eval(cur_rGQM,cur_Robs(any_sac_inds),all_Xmat_shift(any_sac_inds,:));
            sacDelay(cc).gsac_base_LLimp = (baseLL-base_nullLL)/log(2);
            
            %% CREATE SIMULATED SPIKE TRAINS USING THE PRE AND POST-FILTERING MODELS
            post_mod = sacStimProc(cc).gsac_post_mod{post_lambda_off,post_lambda_gain};
            pre_gain_mod = sacStimProc(cc).gsacPreGainMod{pre_lambda};
            
            %get predicted rate from the post-filtering model
            stimG = sum(fgint,2);
            trX{1} = stimG;
            trX{2} = cur_Xsac;
            trX{3} = bsxfun(@times,cur_Xsac,stimG);
            [~,~,sim_prate] = NMMmodel_eval(post_mod,cur_Robs,trX);
            
            %get predicted rate from the pre-filtering model
            [~,pre_pred_rate] = eval_pre_gainmodel( pre_gain_mod, cur_Robs, all_Xmat_shift, cur_Xsac);
            
            %get simulated binned spike responses
            postsim_spikes = poissrnd(sim_prate);
            presim_spikes = poissrnd(pre_pred_rate);
            
            %%
            [~,Tinds] = meshgrid(1:use_nPix_us,1:flen);
            cur_mod_signs = [cur_rGQM.mods(:).sign];
            cur_NL_types = {cur_rGQM.mods(:).NLtype};
            cur_stim_params = NMMcreate_stim_params([1 use_nPix_us]);
            all_filts = [cur_rGQM.mods(:).filtK];
            
            %loop over all stimulus latencies and fit gain kernel models
            %using only that latency
            [gain_filts,offset_filts,out_gain,out_offset] = deal(nan(flen,length(slags)));
            [presim_gain_filts,presim_offset_filts,postsim_gain_filts,postsim_offset_filts] = deal(nan(flen,length(slags)));
            single_mod_filts = nan(flen,use_nPix_us,length(cur_mod_signs));
            [base_gweights, single_gSD,single_rCV] = deal(nan(flen,1));
            [all_flen_stas,all_high_avgs,all_low_avgs,all_hq_avgs,all_lq_avgs] = deal(nan(flen,length(slags)));
            [sac_gains,sac_offsets] = deal(nan(flen,length(slags)));
            for ff = 1:flen
                fprintf('fitting single-latency model %d/%d\n',ff,flen);
                cur_kInds = find(Tinds(:) == ff); %stimulus indices with this time lag
                
                %                 cur_filts = all_filts(cur_kInds,:);
                %                 single_mod = NMMinitialize_model(cur_stim_params,cur_mod_signs,cur_NL_types,[],ones(length(cur_mod_signs),1));
                %                 for ii = 1:length(cur_mod_signs)
                %                     single_mod.mods(ii).filtK = cur_filts(:,ii); %initialize filters to be slices from the original ST filters
                %                     single_mod.mods(ii).reg_params.lambda_d2X = cur_rGQM.mods(ii).reg_params.lambda_d2XT;
                %                     single_mod.mods(ii).reg_params.lambda_L1 = cur_rGQM.mods(ii).reg_params.lambda_L1/10; %decrease sparsity penalty so we dont get so many 0-filters
                %                 end
                %                 single_mod = NMMfit_filters(single_mod,cur_Robs,all_Xmat_shift(:,cur_kInds),[],[],silent);
                %                 [~,~,basemod_pred_rate,cur_gint,~,fgint] = NMMeval_model(single_mod,cur_Robs,all_Xmat_shift(:,cur_kInds));
                %                 single_mod_filts(ff,:,:) = [single_mod.mods(:).filtK];
                %
                
                %initialize filters of a model using only stimulus latency
                %ff
                cur_filts = all_filts(cur_kInds,:); %these are indices of stim elements at latency ff
                init_filts = cell(length(cur_mod_signs),1);
                for ii = 1:length(cur_mod_signs)
                    init_filts{ii} = cur_filts(:,ii);
                end
                single_mod = NMMinitialize_model(cur_stim_params,cur_mod_signs,cur_NL_types,[],ones(length(cur_mod_signs),1),init_filts);
                single_mod_filts(ff,:,:) = [single_mod.mods(:).filtK]; %store stim filters
                [~, ~, ~, ~, ~,fgint] = NMMmodel_eval(single_mod,cur_Robs,all_Xmat_shift(:,cur_kInds)); %get output of each filter
                temp_sp(1:length(cur_mod_signs)) = NMMcreate_stim_params([1 1]); %initialize stim param struct (treating the scalar valued filter outputs as the stim)
                
                %make output of each filter into a cell array
                fout = cell(length(cur_mod_signs),1);
                for ii = 1:length(cur_mod_signs)
                    fout{ii} = fgint(:,ii);
                end
                %initialize a model carrying weights on each stim filter output
                single_mod = NMMinitialize_model(temp_sp,cur_mod_signs,repmat({'lin'},length(cur_mod_signs),1),[],1:length(cur_mod_signs));
                single_mod.spk_NL_params = cur_rGQM.spk_NL_params; %get estimated spk NL function
                single_mod = NMMfit_filters(single_mod,cur_Robs,fout,[],any_sac_inds,silent); %fit scalar weights on the output of each stim filter
                [~,~,basemod_pred_rate,cur_gint,~,fgint] = NMMeval_model(single_mod,cur_Robs,fout); %get outputs of this reweighted model
                
                fgint = bsxfun(@times,fgint,cur_mod_signs); %sign-corrected output of each filter slice
                
                %output of stim model
                stimG = sum(fgint,2); %total output
                single_gSD(ff) = std(stimG); %SD of generating signal as measure of modulation by stimuli at this latency
                single_rCV(ff) = std(basemod_pred_rate)/mean(basemod_pred_rate); %CV of firing rate as a measure of modulation by stimuli at this latency

                %% model-free analysis of saccade modulation.
                normStimG = nanzscore(stimG); %normalize generating signal
                mdpt = nanmedian(normStimG); %median of gen signal
                quarts = prctile(normStimG,[25 75]); %quartiles
                [slag_stas,high_avg,low_avg,hq_avg,lq_avg] = deal(nan(length(slags),1));
                ov_sac_rate = nan(length(slags),1);
                for ss = 1:length(slags)
                    cur_tp = find(cur_Xsac(:,ss) == 1); %all times when there was a saccade at latency tau
                    slag_stas(ss) = mean(cur_Robs(cur_tp).*normStimG(cur_tp)); %saccade trig avg gen signal
                    ov_sac_rate(ss) = mean(cur_Robs(cur_tp)); %avg rate
                    
                    high_pts = cur_tp(normStimG(cur_tp) > mdpt); %good stims at this latency
                    low_pts = cur_tp(normStimG(cur_tp) < mdpt); %bad stims at this latency
                    high_avg(ss) = mean(cur_Robs(high_pts)); %avg rate with good stims
                    low_avg(ss) = mean(cur_Robs(low_pts)); %avg rate with bad stims
                    uq_pts = cur_tp(normStimG(cur_tp) > quarts(2)); %good stims at this latency
                    lq_pts = cur_tp(normStimG(cur_tp) < quarts(1)); %bad stims at this latency
                    hq_avg(ss) = mean(cur_Robs(uq_pts)); %avg rate with good stims
                    lq_avg(ss) = mean(cur_Robs(lq_pts)); %avg rate with bad stims
                end
                all_flen_stas(ff,:) = slag_stas;
                all_high_avgs(ff,:) = high_avg;
                all_low_avgs(ff,:) = low_avg;
                all_hq_avgs(ff,:) = hq_avg;
                all_lq_avgs(ff,:) = lq_avg;
                
                
                %% GET PIXEL TRIG AVGS
                sac_trg_avg = nan(length(slags),1);
                [white_avg,black_avg] = deal(nan(length(cur_kInds),1));
                [white_savg,black_savg] = deal(nan(length(cur_kInds),length(slags)));
                for ii = 1:length(cur_kInds)
                    curwhite = all_Xmat_shift(:,cur_kInds(ii)) == 1;
                    curblack = all_Xmat_shift(:,cur_kInds(ii)) == -1;
                    white_avg(ii) = mean(cur_Robs(curwhite));
                    black_avg(ii) = mean(cur_Robs(curblack));
                    
                    for ss = 1:length(slags)
                        cur_tp = find(cur_Xsac(:,ss) == 1); %all times when there was a saccade at latency tau
                        ccwhite = cur_tp(curwhite(cur_tp));
                        ccblack = cur_tp(curblack(cur_tp));
                        white_savg(ii,ss) = mean(cur_Robs(ccwhite));
                        black_savg(ii,ss) = mean(cur_Robs(ccblack));
                        
                    end
                end
                sacDelay(cc).white_savg(ff,:,:) = white_savg;
                sacDelay(cc).black_savg(ff,:,:) = black_savg;
                sacDelay(cc).white_avg(ff,:) = white_avg;
                sacDelay(cc).black_avg(ff,:) = black_avg;
                
                %% FIT GAIN/OFFSET MODEL
                n_slags = size(cur_Xsac,2);
                n_gains = size(stimG,2);
                
                sac_stim_params(1) = NMMcreate_stim_params(n_gains);
                sac_stim_params(2) = NMMcreate_stim_params(n_slags);
                
                tr_stim{1} = stimG; %total output of stim model
                tr_stim{2} = cur_Xsac; %saccade timing indicator matrix
                tr_stim{3} = reshape(bsxfun(@times,cur_Xsac,reshape(stimG,[],1,n_gains)), size(stimG,1),[]); %product of saccade indicator matrix and gen signal
                sac_stim_params(3) = NMMcreate_stim_params([n_slags n_gains]);
                
                mod_signs = [1 1 1];
                Xtargets = [1 2 3];
                NL_types = {'lin','lin','lin'};
                silent = 1;
                %set lower tolerances to ensure convergence, even for
                %latencies where likelihood is very flat
                optim_params.optTol = 1e-8;
                optim_params.progTol = 1e-11;
                sac_reg_params = NMMcreate_reg_params('boundary_conds',repmat([Inf 0 0],length(mod_signs),1));
                
                %initialize sacmod filter to zeros. This ensures that we
                %get the exact same estimates every time, even for the
                %estimates for lags with minimal stim response (where the
                %likelihood surface is extremely flat).
                clear init_filts
                init_filts{1} = 1;
                init_filts{2} = zeros(length(slags),1);
                init_filts{3} = zeros(length(slags),1);
                
                %estimates for the real data
                cur_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets,init_filts);
                cur_mod.mods(1).filtK(:) = 1; %initialize base gains to 1
                cur_mod.spk_NL_params = cur_rGQM.spk_NL_params;
                cur_mod = NMMadjust_regularization(cur_mod,[2],'lambda_d2T',cent_off_d2T); %temporal smoothness reg for offset filter
                cur_mod = NMMadjust_regularization(cur_mod,[3],'lambda_d2T',10,'lambda_L2',1); %temporal smoothness reg for gain filter
                cur_mod = NMMfit_filters(cur_mod,cur_Robs,tr_stim,[],any_sac_inds,silent,optim_params,[],[2 3]); %estimate saccade filters
                gain_filts(ff,:) = cur_mod.mods(3).filtK;
                offset_filts(ff,:) = cur_mod.mods(2).filtK;
                base_gweights(ff) = cur_mod.mods(1).filtK; 
                %                 modmods{ff} = cur_mod;
                
                %fit linear regression based gain/offset parameters
                [~,~,pred_rate] = NMMeval_model(cur_mod,cur_Robs,tr_stim,[],any_sac_inds);
                for ss = 1:length(slags)
                    cur_tp = find(cur_Xsac(any_sac_inds,ss) == 1); %all times when there was a saccade at latency tau
                    
                    rr = regress(pred_rate(cur_tp),[ones(length(cur_tp),1) basemod_pred_rate(any_sac_inds(cur_tp))]);
                    sac_gains(ff,ss) = rr(2);
                    sac_offsets(ff,ss) = rr(1);
                end
                
                %                 %now for post-filtering sim spikes
                %                 cur_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets,init_filts);
                %                 cur_mod.mods(1).filtK(:) = 1; %initialize base gains to 1
                %                 cur_mod.spk_NL_params = cur_rGQM.spk_NL_params;
                %                 cur_mod = NMMadjust_regularization(cur_mod,[2],'lambda_d2T',cent_off_d2T); %temporal smoothness reg for offset filter
                %                 cur_mod = NMMadjust_regularization(cur_mod,[3],'lambda_d2T',10,'lambda_L2',1); %temporal smoothness reg for gain filter
                %                 cur_mod = NMMfit_filters(cur_mod,postsim_spikes,tr_stim,[],any_sac_inds,silent,optim_params,[],[2 3]); %estimate saccade filters
                %                 postsim_gain_filts(ff,:) = cur_mod.mods(3).filtK;
                %                 postsim_offset_filts(ff,:) = cur_mod.mods(2).filtK;
                %
                %                 %now for pre-filtering sim spikes
                %                 cur_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets,init_filts);
                %                 cur_mod.mods(1).filtK(:) = 1; %initialize base gains to 1
                %                 cur_mod.spk_NL_params = cur_rGQM.spk_NL_params;
                %                 cur_mod = NMMadjust_regularization(cur_mod,[2],'lambda_d2T',cent_off_d2T); %temporal smoothness reg for offset filter
                %                 cur_mod = NMMadjust_regularization(cur_mod,[3],'lambda_d2T',10,'lambda_L2',1); %temporal smoothness reg for gain filter
                %                 cur_mod = NMMfit_filters(cur_mod,presim_spikes,tr_stim,[],any_sac_inds,silent,optim_params,[],[2 3]); %estimate saccade filters
                %                 presim_gain_filts(ff,:) = cur_mod.mods(3).filtK;
                %                 presim_offset_filts(ff,:) = cur_mod.mods(2).filtK;
                
            end
            %% STORE RESULTS
            sacDelay(cc).single_mod_filts = single_mod_filts;
            sacDelay(cc).gain_filts = gain_filts;
            sacDelay(cc).offset_filts = offset_filts;
            sacDelay(cc).single_gSD = single_gSD; %SD of gen signal
            sacDelay(cc).single_rCV = single_rCV; %CV of predicted firing rate
            sacDelay(cc).base_gweights = base_gweights;
            
            sacDelay(cc).sac_gains = sac_gains;
            sacDelay(cc).sac_offsets = sac_offsets;
            
            sacDelay(cc).flen_stas = all_flen_stas;
            sacDelay(cc).high_avgs = all_high_avgs;
            sacDelay(cc).low_avgs = all_low_avgs;
            sacDelay(cc).hq_avgs = all_hq_avgs;
            sacDelay(cc).lq_avgs = all_lq_avgs;
            sacDelay(cc).raw_sac_rate = ov_sac_rate;
            
            %             sacDelay(cc).presim_gain_filts = presim_gain_filts;
            %             sacDelay(cc).presim_offset_filts = presim_offset_filts;
            %             sacDelay(cc).postsim_gain_filts = postsim_gain_filts;
            %             sacDelay(cc).postsim_offset_filts = postsim_offset_filts;
        end
    else
        sacDelay(cc).used = false;
    end
end

%%
anal_dir = ['/home/james/Analysis/bruce/' Expt_name '/FINsac_mod/'];

sname = [sname sprintf('_ori%d',bar_ori)];
if fit_unCor
    sname = [sname '_unCor'];
end
cd(anal_dir)
save(sname,'targs','slags','dt','sacDelay');

