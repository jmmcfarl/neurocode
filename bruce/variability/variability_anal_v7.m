
clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements//');
addpath('~/James_scripts/bruce/processing/');
addpath('~/James_scripts/bruce/saccade_modulation/');
addpath('~/James_scripts/TentBasis2D/');

global Expt_name bar_ori use_MUA

Expt_name = 'M296';
use_MUA = false;
bar_ori = 45; %bar orientation to use (only for UA recs)

mod_data_name = 'corrected_models2';
ep_dist_bin_edges = linspace(-1,1,100);

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

% if Expt_num >= 280
%     data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
% else
    data_dir = ['~/Data/bruce/' Expt_name];
% end

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
et_hres_anal_name = [et_anal_name '_hres'];

et_mod_data_name = [et_mod_data_name sprintf('_ori%d',bar_ori)];
et_anal_name = [et_anal_name sprintf('_ori%d',bar_ori)];
mod_data_name = [mod_data_name sprintf('_ori%d',bar_ori)];
et_hres_anal_name = [et_hres_anal_name sprintf('_ori%d',bar_ori)];

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
        
%         cur_t_edges_up = cur_t_edges(1):up_dt:cur_t_edges(end);
%         cur_t_axis_up = 0.5*cur_t_edges_up(1:end-1) + 0.5*cur_t_edges_up(2:end);
        
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
[all_binned_mua,all_binned_sua,Clust_data,all_su_spk_times,~,all_mu_spk_times] = ...
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

%% Create set of TR trials
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

%% CREATE TIME AXIS FOR REPEAT TRIALS
rpt_taxis = (1:round(trial_dur)/dt)*dt-dt/2;
rpt_taxis(rpt_taxis < beg_buffer) = [];
rpt_taxis(trial_dur - rpt_taxis < end_buffer) = [];

all_trial_dur = all_trial_end_times-all_trial_start_times;

rpt_trials = find(all_trial_Se == rpt_seed);
rpt_trials(all_trial_rptframes(rpt_trials) > 0) = []; %get rid of any repeat trials where there were repeat frames

n_rpts = length(rpt_trials);
all_rpt_inds = find(ismember(all_trialvec(used_inds),rpt_trials));
all_nonrpt_inds = find(~ismember(all_trialvec(used_inds),rpt_trials));

%% Recon retinal stim with high res ET data
cd(et_dir)
load(et_hres_anal_name)
[hres_fin_tot_corr,hres_tot_std] = construct_eye_position(best_fix_cor,best_fix_std,...
    drift_post_mean,drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);
hres_sp_dx = et_params.sp_dx;
hres_ep = hres_fin_tot_corr*hres_sp_dx;

hres_fin_tot_corr_LOO = nan(length(loo_set),NT);
for ss = 1:length(loo_set)
[hres_fin_tot_corr_LOO(ss,:)] = construct_eye_position(best_fix_cor,best_fix_std,...
    drift_post_mean_LOO(ss,:),drift_post_std_LOO(ss,:),fix_ids,trial_start_inds,trial_end_inds,sac_shift);    
end
hres_fin_tot_corr_LOO = hres_fin_tot_corr_LOO*hres_sp_dx;


%% DEFINE IN-SAC AND IN-BUFF INDS
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

%% PROCESS REPEAT EYE POSITION DATA

back_look = 7; %look back this many time steps to parse EP trajectories

%make an anticausal filter for processing EP history
back_kern = zeros(back_look*2+1,1);
back_kern(1:back_look+1) = 1;
back_kern = flipud(back_kern/sum(back_kern));

%eye position during repeats
rpt_EP = hres_ep(all_rpt_inds);

rpt_EPnan = rpt_EP;
rpt_EPnan(in_sac_inds(all_rpt_inds)) = nan;
rpt_EPnan(in_blink_inds(all_rpt_inds)) = nan;

%initialize a time-embedded version of the EPs
sp = NMMcreate_stim_params(back_look);
rpt_EP_emb = create_time_embedding(rpt_EPnan(:),sp);

rpt_EP_emb_LOO = nan(length(loo_set),size(rpt_EP_emb,1),back_look);
for ll = 1:length(loo_set)
    temp = hres_fin_tot_corr_LOO(ll,all_rpt_inds)';
    temp(isnan(rpt_EPnan)) = nan;
   rpt_EP_emb_LOO(ll,:,:,:) = create_time_embedding(temp,sp);
end

%make TBT matrices of EP data
in_sac = false(n_rpts,length(rpt_taxis));
full_EP = nan(n_rpts,length(rpt_taxis));
full_EP_emb = nan(n_rpts,length(rpt_taxis),back_look);
full_EP_emb_LOO = nan(length(loo_set),n_rpts,length(rpt_taxis),back_look);
for ii = 1:length(rpt_trials)
    cur_inds = find(all_trialvec(used_inds(all_rpt_inds)) == rpt_trials(ii));
    in_sac(ii,1:length(cur_inds)) = in_sac_inds(all_rpt_inds(cur_inds)) | in_blink_inds(all_rpt_inds(cur_inds));
    full_EP(ii,1:length(cur_inds)) = rpt_EP(cur_inds);
    full_EP_emb(ii,1:length(cur_inds),:) = rpt_EP_emb(cur_inds,:);    

    for ll = 1:length(loo_set)
        full_EP_emb_LOO(ll,ii,1:length(cur_inds),:) = (rpt_EP_emb_LOO(ll,cur_inds,:));    
    end
end

%nan out EP data either in blinks or sacs
full_EP(in_sac) = nan;
ep_naninds = isnan(full_EP(:)); %find within-sac inds

%% bin spike data at desired dt res and interpolate ET data onto this up-sampled t-axis
up_dt = 0.005;
new_tlen = ceil(length(rpt_taxis)*dt/up_dt);
rpt_bs_up = nan(n_rpts,new_tlen,n_chs);
rpt_EP = nan(n_rpts,new_tlen);
rpt_EP_emb = nan(n_rpts,new_tlen,back_look);
rpt_EP_emb_LOO = nan(length(loo_set),n_rpts,new_tlen,back_look);
for ii = 1:length(rpt_trials)
    
    cur_st = all_trial_start_times(rpt_trials(ii)) + beg_buffer;
    cur_et = all_trial_end_times(rpt_trials(ii)) - end_buffer;
    cur_up_tedges = cur_st:up_dt:cur_et;
    
    if length(cur_up_tedges) > (new_tlen + 1)
        cur_up_tedges = cur_up_tedges(1:(new_tlen+1));
    elseif length(cur_up_tedges) < (new_tlen + 1)
       cur_up_tedges = [cur_up_tedges cur_up_tedges(end) + up_dt];
    end
    cur_up_tcents = 0.5*cur_up_tedges(1:end-1) + 0.5*cur_up_tedges(2:end);

    cur_inds = find(all_trialvec(used_inds) == rpt_trials(ii));
    cur_tax = all_t_axis(used_inds(cur_inds));

    cur_block_id = unique(all_blockvec(used_inds(cur_inds)));
    if length(cur_block_id) ~= 1
        error('!')
    end

    cur_up_EP = interp1(cur_tax,full_EP(ii,:),cur_up_tcents);
    cur_up_EP_emb = interp1(cur_tax,squeeze(full_EP_emb(ii,:,:)),cur_up_tcents);
    cur_up_EP(1) = cur_up_EP(2); cur_up_EP(end) = cur_up_EP(end-1);
    cur_up_EP_emb(1,:) = cur_up_EP_emb(2,:); cur_up_EP_emb(end,:) = cur_up_EP_emb(end-1,:);
    rpt_EP(ii,:) = cur_up_EP;
    rpt_EP_emb(ii,:,:) = cur_up_EP_emb;
    
    for ll = 1:length(loo_set)
        cur_up_EP_emb = interp1(cur_tax,squeeze(full_EP_emb_LOO(ll,ii,:,:)),cur_up_tcents);
        cur_up_EP_emb(1,:) = cur_up_EP_emb(2,:); cur_up_EP_emb(end,:) = cur_up_EP_emb(end-1,:);
        rpt_EP_emb_LOO(ll,ii,:,:) = cur_up_EP_emb;
    end
    
    for cc = 1:n_probes
        cur_bs = histc(all_mu_spk_times{cc},cur_up_tedges);
        rpt_bs_up(ii,:,cc) = cur_bs(1:end-1);
    end
    for cc = 1:length(all_su_spk_times)
        if ~isnan(Clust_data.SU_block_probes(cc,cur_block_id))
        cur_bs = histc(all_su_spk_times{cc},cur_up_tedges);
        rpt_bs_up(ii,:,cc+n_probes) = cur_bs(1:end-1);
        end 
    end
end

ep_naninds_up = find(isnan(rpt_EP(:)));
%% SUBTRACT OUT MEAN RATES
%nan out within-sac (or within blink) times in psth and mod prates
rpt_bs_up = reshape(rpt_bs_up,[],n_chs);
rpt_bs_up(ep_naninds_up,:) = nan;
rpt_bs_up = reshape(rpt_bs_up,n_rpts,[],n_chs);

psths_raw = squeeze(nanmean(rpt_bs_up)); %raw psths

%compute a version of bs matrix with trial-avgs subtracted out
rpt_bs_up_ms = rpt_bs_up;
trial_avg_rates = nanmean(rpt_bs_up_ms,2);%compute trial avg rates
rpt_bs_up_ms = bsxfun(@minus,rpt_bs_up_ms,trial_avg_rates); %subtract off trial-mean rates

trial_avg_rates = squeeze(trial_avg_rates);

%compute psth of tavg-subtracted
psths_ms = squeeze(nanmean(rpt_bs_up_ms));

%subtract off overall mean rate from rpt_bs
rpt_avg_rates = nanmean(reshape(rpt_bs_up,[],n_chs)); %mean rates during rpts (after nanning out)
rpt_bs_up = bsxfun(@minus,rpt_bs_up,reshape(rpt_avg_rates,1,1,[])); %subtract off overall mean rate
rpt_avg_rates_ms = nanmean(reshape(rpt_bs_up_ms,[],n_chs)); %mean rates during rpts (after nanning out)
rpt_bs_up_ms = bsxfun(@minus,rpt_bs_up_ms,reshape(rpt_avg_rates_ms,1,1,[])); %subtract off overall mean rate

%% COMPUTE VARIANCES
raw_psth_var = nanvar(psths_raw);
ms_psth_var = nanvar(psths_ms);

unfolded_psth_raw = reshape(permute(psths_raw,[2 1 3]),[],n_chs);
raw_resp_var = nanvar(unfolded_psth_raw); %total variance of binned spike data

unfolded_psth_ms = reshape(permute(psths_ms,[2 1 3]),[],n_chs);
ms_resp_var = nanvar(unfolded_psth_ms); %total variance of binned spike data

raw_direct_noise_vars = squeeze(nanmean(nanvar(rpt_bs_up)));
ms_direct_noise_vars = squeeze(nanmean(nanvar(rpt_bs_up_ms)));

n_utrials = squeeze(mean(sum(~isnan(rpt_bs_up))));
raw_psth_var_cor = raw_psth_var.*(n_utrials'./(n_utrials'-1)) - raw_direct_noise_vars'./n_utrials';
ms_psth_var_cor = ms_psth_var.*(n_utrials'./(n_utrials'-1)) - ms_direct_noise_vars'./n_utrials';

trial_avg_var = squeeze(nanvar(trial_avg_rates)); %variance of trial-avg rates

%%
% for ss = 1:n_chs
%     Rpt_Data(ss).mod_resp_var = mod_tot_var(ss);
%     Rpt_Data(ss).tot_resp_var = tot_resp_var(ss);
%     Rpt_Data(ss).dir_noise_var = direct_noise_vars(ss);
%     Rpt_Data(ss).ms_resp_var = resp_var(ss);
%     Rpt_Data(ss).emp_psth = all_psths(:,ss);
%     Rpt_Data(ss).psth_var = psth_var(ss);
%     Rpt_Data(ss).psth_var_cor = psth_var_cor(ss);
%     Rpt_Data(ss).rpt_avg_rate = rpt_avg_rates(ss);
%     Rpt_Data(ss).n_utrials = n_utrials(ss);
%     uset = ~isnan(unfolded_psth_ms(:,ss));
%     Rpt_Data(ss).rpt_ep_sd = robust_std_dev(hres_ep(uset));
% end


%% SORT TBT EY EPOSITION FOR VISUALIZATION PURPOSES
% [a,b] = sort(rpt_EP);
% sorted_full_psth = nan(size(rpt_bs_up_ms));
% for tt = 1:new_tlen
%     sorted_full_psth(:,tt,:) = rpt_bs_up_ms(b(:,tt),tt,:);
%     sorted_full_psth(isnan(a(:,tt)),tt,:) = nan;
% end
%% COMPUTE VARIANCES FOR BINNED DELTA ES FOR ALL CELLS (using full recon et)
maxlag_ED = 0.15;
ED_space = 0.0025;
ED_bin_edges = 0:ED_space:maxlag_ED;
ED_bin_centers = (ED_bin_edges(1:end-1)+ED_bin_edges(2:end))/2;
[II,JJ] = meshgrid(1:n_rpts);

cur_XC = nan(new_tlen,length(ED_bin_centers),n_chs);
cur_cnt = zeros(length(ED_bin_centers),n_chs);
rand_XC = nan(new_tlen,n_chs);
% cur_XC = nan(new_tlen,length(ED_bin_centers),length(targs));
% cur_cnt = zeros(length(ED_bin_centers),length(targs));
% rand_XC = nan(new_tlen,length(targs));
for tt = 1:new_tlen
%     Y1 = squeeze(rpt_bs_up(:,tt,:));
    Y1 = squeeze(rpt_bs_up_ms(:,tt,:));
%     Y1 = squeeze(rpt_bs_up_ms(:,tt,targs));

    cur_Dmat = abs(squareform(pdist(squeeze(rpt_EP_emb(:,tt,:)))))/sqrt(back_look);
    cur_Dmat(logical(eye(n_rpts))) = nan;
    for jj = 1:length(ED_bin_centers)
        curset = find(cur_Dmat > ED_bin_edges(jj) & cur_Dmat <= ED_bin_edges(jj+1));
        cur_XC(tt,jj,:) = squeeze(nanmean(bsxfun(@times,Y1(II(curset),:),Y1(JJ(curset),:)),1));
        cur_cnt(jj,:) = cur_cnt(jj,:) + sum(~isnan(Y1(II(curset),:)));
    end
    curset = ~isnan(cur_Dmat);
    rand_XC(tt,:) = squeeze(nanmean(bsxfun(@times,Y1(II(curset),:),Y1(JJ(curset),:)),1));
end
new_psth_var = nanmean(rand_XC);

var_ep_binned = squeeze(nanmean(cur_XC));
all_relprobs = bsxfun(@rdivide,cur_cnt,sum(cur_cnt));

%% FIT SPLINES TO VAR VS ET FUNCTIONS (again using full et recon)
spline_spacing = 0.015;
spline_knots = spline_spacing:spline_spacing:maxlag_ED;
spline_eval = 0:0.005:maxlag_ED;

nboots = 0;

var_spline_ZPT = nan(n_chs,1);
var_spline_funs = nan(n_chs,length(spline_eval));
var_spline_ZPT_boot = nan(n_chs,nboots);
for ii = 1:n_chs
    x = ED_bin_centers;
    y = squeeze(var_ep_binned(:,ii));
    bad = find(isnan(y));  x(bad) = []; y(bad) = [];
    ss = fnxtr(csape(spline_knots,y(:).'/fnval(fnxtr(csape(spline_knots,eye(length(spline_knots)),'var')),x(:).'),'var'));
    var_spline_ZPT(ii) = fnval(ss,0);
    var_spline_funs(ii,:) = fnval(ss,spline_eval);
    
    npts = length(x);
    for nn = 1:nboots
        rsmp = randi(npts,npts,1);
        ss = fnxtr(csape(spline_knots,y(rsmp).'/fnval(fnxtr(csape(spline_knots,eye(length(spline_knots)),'var')),x(rsmp).'),'var'));
        var_spline_ZPT_boot(ii,nn) = fnval(ss,0);
    end
end

% psth_var_frac = ms_psth_var'./var_spline_ZPT;
new_psth_var_frac = new_psth_var'./var_spline_ZPT;
psth_var_frac_cor = ms_psth_var_cor'./var_spline_ZPT;

new_psth_var_frac_boot = bsxfun(@times,(var_spline_ZPT_boot.^(-1)),new_psth_var');
new_psth_var_frac_CIs = prctile(new_psth_var_frac_boot',[5 95]);

%% COMPUTE AUTOCOVARIANCES FOR BINNED DELTA ES FOR ALL CELLS (USING TARGS)

%create a set of temporally shifted versions of the spiking data
max_tlag = round(0.1/up_dt); %max time lag for computing autocorrs
tlags = [0:max_tlag];

zpad_len = length(tlags);
rpt_bs_up_zpad = cat(2,zeros(n_rpts,zpad_len,n_chs),rpt_bs_up_ms,zeros(n_rpts,zpad_len,n_chs));
rpt_EP_emb_zpad = cat(2,zeros(n_rpts,zpad_len,back_look),rpt_EP_emb,zeros(n_rpts,zpad_len,back_look));

% covar_ep_binned = zeros(length(ED_bin_centers),length(targs),length(tlags));
% covar_rand = zeros(length(targs),length(tlags));
% for cc = 1:length(targs)
%     fprintf('Computing autocovariances with unit %d\n',targs(cc));
    cur_XC = nan(new_tlen,length(ED_bin_centers),length(tlags),length(targs));
    rand_XC = nan(new_tlen,length(tlags),length(targs));
    obs_XC = nan(new_tlen,length(tlags),length(targs));
    for tt = 1:new_tlen
        fprintf('Time %d of %d\n',tt,new_tlen);
        for kk = 1:length(tlags)
            Y1 = squeeze(rpt_bs_up_ms(:,tt,targs));
            Y2 = squeeze(rpt_bs_up_zpad(:,tt+zpad_len+tlags(kk),targs));
            
%             cur_Dmat = abs(squareform(pdist(squeeze(rpt_EP_emb(:,tt,:)))))/sqrt(back_look);
            cur_Dmat = abs(pdist2(squeeze(rpt_EP_emb(:,tt,:)),squeeze(rpt_EP_emb_zpad(:,tt+zpad_len+tlags(kk),:))))/sqrt(back_look);
            cur_Dmat(logical(eye(n_rpts))) = nan;
            for jj = 1:length(ED_bin_centers)
                curset = find(cur_Dmat > ED_bin_edges(jj) & cur_Dmat <= ED_bin_edges(jj+1));
%                 cur_XC(tt,jj,:) = squeeze(nanmean(bsxfun(@times,Y1(II(curset),:),Y2(JJ(curset))),1));
                cur_XC(tt,jj,kk,:) = squeeze(nanmean(Y1(II(curset),:).*Y2(JJ(curset),:)));
            end
            curset = find(~isnan(cur_Dmat));
%             rand_XC(tt,:,:) = squeeze(nanmean(bsxfun(@times,Y1(II(curset),:),Y2(JJ(curset))),1));
            rand_XC(tt,kk,:) = squeeze(nanmean(Y1(II(curset),:).*Y2(JJ(curset),:)));
            obs_XC(tt,kk,:) = squeeze(nanmean(Y1(II(curset),:).*Y2(II(curset),:)));
        end
    end
    covar_ep_binned = squeeze(nanmean(cur_XC));
    covar_rand = squeeze(nanmean(rand_XC));
    covar_obs = squeeze(nanmean(obs_XC));
% end

acvar_spline_ZPT = nan(length(tlags),length(targs));
for cc = 1:length(targs)
    for ii = 1:length(tlags)
        x = ED_bin_centers;
        y = squeeze(covar_ep_binned(:,ii,cc));
        bad = find(isnan(y));  x(bad) = []; y(bad) = [];
        ss = fnxtr(csape(spline_knots,y(:).'/fnval(fnxtr(csape(spline_knots,eye(length(spline_knots)),'var')),x(:).'),'var'));
        acvar_spline_ZPT(ii,cc) = fnval(ss,0);
    end
end

tlag_corr_fac = (new_tlen - tlags)/new_tlen;
covar_rand = bsxfun(@rdivide,covar_rand,tlag_corr_fac');
acvar_spline_ZPT = bsxfun(@rdivide,acvar_spline_ZPT,tlag_corr_fac');

int_acvar_cor = cumtrapz(tlags,acvar_spline_ZPT);
int_acvar_raw = cumtrapz(tlags,covar_rand);
int_acvar_obs = cumtrapz(tlags,covar_obs);

ialpha_acvar = int_acvar_cor./int_acvar_raw;
% new_psth_var_frac = new_psth_var'./var_spline_ZPT;

% tbt_acovs = nan(n_rpts,length(tlags),length(targs));
% for ii = 1:n_rpts
%     for cc = 1:length(targs)
%         Y1 = squeeze(rpt_bs_up_ms(ii,:,targs(cc)));
%     end
% end
%% VIEW SPLINE FITS
close all
for ii = 1:n_chs
    fprintf('Unit %d\n',ii);
   plot(ED_bin_centers,var_ep_binned(:,ii),'.');
   hold on
   plot(spline_eval,var_spline_funs(ii,:),'r')
%    line(spline_eval([1 end]),psth_var([ii ii]),'color','m');
   line(spline_eval([1 end]),new_psth_var([ii ii]),'color','m');
%    line(spline_eval([1 end])*sp_dx,new_psth_var(ii)/new_mod_var_frac2(ii) + [0 0],'color','g');
   line(spline_eval([1 end]),[0 0],'color','k')
   pause
   clf
end

%% ESTIMATE FANO FACTOR CORRECTIONS USING NEAREST ET (FULL CORRECTION)
rpt_bs_up_ps = bsxfun(@minus,rpt_bs_up,nanmean(rpt_bs_up)); %subtract off raw psths
actual_tvars = squeeze(nanvar(rpt_bs_up,[],1)); %direct est of across-trial variance at each time point

%for doing adjacent-trial correction
adj1 = 1:(n_rpts-1);
adj2 = 2:n_rpts;

timewise_vars = nan(new_tlen,n_chs);
adj_vars = nan(new_tlen,n_chs);
for tt = 1:new_tlen
     Y1 = squeeze(rpt_bs_up_ps(:,tt,:));

     cur_Dmat = abs(squareform(pdist(squeeze(rpt_EP_emb(:,tt,:)))))/sqrt(back_look);
     cur_Dmat(logical(eye(n_rpts))) = nan;
          
     [closest_sep,closest_pair] = nanmin(cur_Dmat);
     uset = ~isnan(closest_sep);
    
     timewise_vars(tt,:) = squeeze(nanmean(bsxfun(@times,Y1(uset,:),Y1(closest_pair(uset),:)),1));
     adj_vars(tt,:) = squeeze(nanmean(bsxfun(@times,Y1(adj1,:),Y1(adj2,:)),1));   
end
corrected_noise_vars = actual_tvars - timewise_vars;
adj_corrected_noise_vars = actual_tvars - adj_vars;

[base_FF,cor_FF,adj_FF] = deal(nan(n_chs,1));
for ii = 1:n_chs
   base_FF(ii) = regress(actual_tvars(:,ii),psths_raw(:,ii));
   cor_FF(ii) = regress(corrected_noise_vars(:,ii),psths_raw(:,ii));
   adj_FF(ii) = regress(adj_corrected_noise_vars(:,ii),psths_raw(:,ii));
end

%% RE-ESTIMATE SINGLE-NEURON ET CORRECTIONS USING LOO ET SIGNALS

loo_spline_ZPT = nan(n_chs,1);
loo_cor_FF = nan(n_chs,1);
for ll = 1:length(loo_set)
    
    fprintf('LOO %d of %d\n',ll,length(loo_set));    
    cur_XC = nan(new_tlen,length(ED_bin_centers));
    cur_time_vars = nan(new_tlen,1);
    for tt = 1:new_tlen
        Y1 = squeeze(rpt_bs_up_ms(:,tt,loo_set(ll)));
        cur_Dmat = abs(squareform(pdist(squeeze(rpt_EP_emb_LOO(ll,:,tt,:)))))/sqrt(back_look);
        cur_Dmat(logical(eye(n_rpts))) = nan;
        cur_Dmat(exclude_trials,:) = nan;
        for jj = 1:length(ED_bin_centers)
            curset = find(cur_Dmat > ED_bin_edges(jj) & cur_Dmat <= ED_bin_edges(jj+1));
            cur_XC(tt,jj) = nanmean(Y1(II(curset)).*Y1(JJ(curset)));
        end
              
        Y1 = squeeze(rpt_bs_up_ps(:,tt,loo_set(ll)));
        [closest_sep,closest_pair] = nanmin(cur_Dmat);
        uset = ~isnan(closest_sep);
        cur_time_vars(tt) = nanmean(Y1(uset).*Y1(closest_pair(uset)));
        
    end
    loo_var_ep_binned = squeeze(nanmean(cur_XC));
    
    x = ED_bin_centers;
    y = loo_var_ep_binned;
    bad = find(isnan(y));  x(bad) = []; y(bad) = [];
    ss = fnxtr(csape(spline_knots,y(:).'/fnval(fnxtr(csape(spline_knots,eye(length(spline_knots)),'var')),x(:).'),'var'));
    loo_spline_ZPT(loo_set(ll)) = fnval(ss,0);

    cur_corrected_noise_vars = actual_tvars(:,loo_set(ll)) - cur_time_vars;
%     loo_cor_FF(loo_set(ll)) = all_psths_raw(:,loo_set(ll))\cur_corrected_noise_vars;
    loo_cor_FF(loo_set(ll)) = regress(cur_corrected_noise_vars,psths_raw(:,loo_set(ll)));
end

loo_psth_var_frac = new_psth_var'./loo_spline_ZPT;

%%
for ss = 1:n_chs
    Rpt_Data(ss).ep_binned_respvar = var_ep_binned(:,ss);
    Rpt_Data(ss).ep_binned_modvar = var_ep_mod(:,ss);
    Rpt_Data(ss).rand_psth_var = new_psth_var(ss);
    Rpt_Data(ss).rand_mod_var = new_mod_var(ss);
    Rpt_Data(ss).ep_spline_respvar = var_spline_funs(ss,:);
    Rpt_Data(ss).spline_resp_ZPT = var_spline_ZPT(ss);
    Rpt_Data(ss).spline_mod_ZPT = var_spline_mZPT(ss);
    
    Rpt_Data(ss).base_FF = base_FF(ss);
    Rpt_Data(ss).cor_FF = cor_FF(ss);
    Rpt_Data(ss).loo_spline_ZPT = loo_spline_ZPT(ss);
    Rpt_Data(ss).loo_cor_FF = loo_cor_FF(ss);
end


%% ESTIMATE FULL EP-BINNED XCOVs
%create a set of temporally shifted versions of the spiking data
max_tlag = 10; %max time lag for computing autocorrs
tlags = [-max_tlag:max_tlag];
full_psth_shifted = nan(n_rpts,length(rpt_taxis),n_chs,length(tlags));
for tt = 1:length(tlags)
    full_psth_shifted(:,:,:,tt) = shift_matrix_Nd(full_psth_ms,-tlags(tt),2);
end

covar_ep_binned = zeros(length(ED_bin_centers),n_chs,n_chs,length(tlags));
covar_rand = zeros(n_chs,n_chs,length(tlags));
for cc = 1:length(targs)
    fprintf('Computing covariances with unit %d\n',targs(cc));
    cur_XC = nan(length(rpt_taxis),length(ED_bin_centers),n_chs,length(tlags));
    rand_XC = nan(length(rpt_taxis),n_chs,length(tlags));
    for tt = 1:length(rpt_taxis)
        Y1 = squeeze(full_psth_shifted(:,tt,:,:));
        Y2 = squeeze(full_psth_ms(:,tt,targs(cc)));
        
        %         cur_Dmat = abs(squareform(pdist(full_EP(:,tt))));
        %         cur_Dmat(logical(eye(n_rpts))) = nan;
        cur_Dmat = abs(squareform(pdist(squeeze(full_EP_emb(:,tt,:)))))/sqrt(back_look);
        cur_Dmat(logical(eye(n_rpts))) = nan;
        for jj = 1:length(ED_bin_centers)
            curset = find(cur_Dmat > ED_bin_edges(jj) & cur_Dmat <= ED_bin_edges(jj+1));
            cur_XC(tt,jj,:,:) = squeeze(nanmean(bsxfun(@times,Y1(II(curset),:,:,:),Y2(JJ(curset),:,:)),1));
        end
        curset = find(~isnan(cur_Dmat));
        rand_XC(tt,:,:) = squeeze(nanmean(bsxfun(@times,Y1(II(curset),:,:,:),Y2(JJ(curset),:,:)),1));
    end
    covar_ep_binned(:,:,targs(cc),:) = squeeze(nanmean(cur_XC));
    covar_rand(:,targs(cc),:) = squeeze(nanmean(rand_XC));
end

%% SPLINE FIT FOR XCOVS
covar_spline_ZPT = nan(n_chs,n_chs,length(tlags));
for ii = 1:n_chs
    for jj = 1:length(targs)
        for kk = 1:length(tlags)
            x = ED_bin_centers;
            y = squeeze(covar_ep_binned(:,ii,targs(jj),kk));
            bad = find(isnan(y));  x(bad) = []; y(bad) = [];
            ss = fnxtr(csape(spline_knots,y(:).'/fnval(fnxtr(csape(spline_knots,eye(length(spline_knots)),'var')),x(:).'),'var'));
            covar_spline_ZPT(ii,targs(jj),kk) = fnval(ss,0);
        end
    end
end

%% ESTIMATE EMPIRICAL XCOV MATS, AS WELL AS PSTH XCOVS

psth_xcov = nan(n_chs,n_chs,2*max_tlag+1);
obs_xcov = nan(n_chs,n_chs,2*max_tlag+1);
trial_xcov = nan(n_chs,n_chs);
for ii = 1:n_chs
    for jj = 1:n_chs
        trial_xcov(ii,jj) = nanmean(trial_avg_rates(:,ii).*trial_avg_rates(:,jj));
        psth_xcov(ii,jj,:) = xcov(all_psths(:,ii),all_psths(:,jj),max_tlag,'biased');
        uset = find(~isnan(unfolded_psth_ms(:,ii)) & ~isnan(unfolded_psth_ms(:,jj)));
        if ~isempty(uset)
            obs_xcov(ii,jj,:) = xcov(unfolded_psth_ms(uset,ii),unfolded_psth_ms(uset,jj),max_tlag,'biased');
        end
    end
end
trial_vars = diag(trial_xcov);
%%
%nan-out elements of xcov mats corresponding to nearby MU channels for the
%SUs
zlag = find(tlags == 0);
nan_locs = false(n_chs,n_chs,length(tlags));
for ii = 1:n_chs
    %for MUs, nan-out the zero-point and surrounding channels
    if ii <= n_probes
        surr_units = [ii-1 ii+1];
        surr_units(surr_units < 1 | surr_units > n_probes) = [];
        nan_locs(ii,surr_units,zlag) = true;
        nan_locs(ii,ii,:) = true;
    else
        temp = find(targs == ii);
        su_pr = SU_probes(temp);
        if ~isnan(su_pr)
            surr_units = su_pr + [-1 0 1];
            surr_units(surr_units < 1 | surr_units > n_probes) = [];
            nan_locs(ii,surr_units,zlag) = true;
            nan_locs(surr_units,ii,zlag) = true;
            nan_locs(ii,ii,:) = true;
        end
    end
end

%%
est_sig_var = var_spline_ZPT;
est_sig_var(est_sig_var < psth_var') = psth_var(est_sig_var < psth_var');
est_sig_var = est_sig_var';

est_noise_var = resp_var' - var_spline_ZPT;
est_noise_var(est_noise_var < direct_noise_vars) = direct_noise_vars(est_noise_var < direct_noise_vars);
est_noise_var = est_noise_var';

poss_cnt_wins = [0:max_tlag];
psth_cnt_cov = nan(n_chs,n_chs,length(poss_cnt_wins));
obs_cnt_cov = nan(n_chs,n_chs,length(poss_cnt_wins));
spline_cnt_cov = nan(n_chs,n_chs,length(poss_cnt_wins));
for cw = 1:length(poss_cnt_wins)
    ulags = find(abs(tlags) <= poss_cnt_wins(cw));
    psth_cnt_cov(:,:,cw) = squeeze(sum(psth_xcov(:,:,ulags),3));
    obs_cnt_cov(:,:,cw) = squeeze(sum(obs_xcov(:,:,ulags),3));
    spline_cnt_cov(:,:,cw) = squeeze(sum(covar_spline_ZPT(:,:,ulags),3));
end

psth_cnt_var = nan(n_chs,length(poss_cnt_wins));
obs_cnt_var = nan(n_chs,length(poss_cnt_wins));
psth_noise_var = nan(n_chs,length(poss_cnt_wins));
for ss = 1:n_chs
    psth_cnt_var(ss,:) = squeeze(psth_cnt_cov(ss,ss,:));
    obs_cnt_var(ss,:) = squeeze(obs_cnt_cov(ss,ss,:));
end
psth_noise_cnt_var = obs_cnt_var - psth_cnt_var;
%%
for ss = 1:length(targs)
    Rpt_Data(targs(ss)).rand_xcov = squeeze(covar_rand(:,targs(ss),:));
    Rpt_Data(targs(ss)).psth_xcov = squeeze(psth_xcov(:,targs(ss),:));
    Rpt_Data(targs(ss)).emp_xcov = squeeze(obs_xcov(:,targs(ss),:));
    Rpt_Data(targs(ss)).spline_xcov = squeeze(covar_spline_ZPT(:,targs(ss),:));
    Rpt_Data(targs(ss)).tot_varnorm = sqrt(resp_var(targs(ss))*resp_var');
    Rpt_Data(targs(ss)).spline_varnorm = sqrt(est_sig_var(targs(ss))*est_sig_var');
    Rpt_Data(targs(ss)).spline_noise_varnorm = sqrt(est_noise_var(targs(ss))*est_noise_var');
    Rpt_Data(targs(ss)).psth_varnorm = sqrt(psth_var(targs(ss))*psth_var');
    Rpt_Data(targs(ss)).noise_varnorm = sqrt(direct_noise_vars(targs(ss))*direct_noise_vars);
    Rpt_Data(targs(ss)).geom_mean = sqrt(rpt_avg_rates(targs(ss))*rpt_avg_rates');
    Rpt_Data(targs(ss)).nan_locs = squeeze(nan_locs(targs(ss),:,:));
    
    cur_trial_corrs = trial_xcov(targs(ss),:)./sqrt(trial_vars(targs(ss)).*trial_vars');
    Rpt_Data(targs(ss)).trial_corrs = cur_trial_corrs;
    
    Rpt_Data(targs(ss)).otherunit_nutrials = n_utrials(:);
    Rpt_Data(targs(ss)).otherunit_avgrates = rpt_avg_rates(:);
    
    Rpt_Data(targs(ss)).psth_cnt_cov = squeeze(psth_cnt_cov(:,targs(ss),:));
    Rpt_Data(targs(ss)).obs_cnt_cov = squeeze(obs_cnt_cov(:,targs(ss),:));
    Rpt_Data(targs(ss)).spline_cnt_cov = squeeze(spline_cnt_cov(:,targs(ss),:));
    Rpt_Data(targs(ss)).psth_cnt_varnorm = sqrt(bsxfun(@times,psth_cnt_var,psth_cnt_var(targs(ss),:)));
    Rpt_Data(targs(ss)).obs_cnt_varnorm = sqrt(bsxfun(@times,obs_cnt_var,obs_cnt_var(targs(ss),:)));
    Rpt_Data(targs(ss)).psthnoise_cnt_varnorm = sqrt(bsxfun(@times,psth_noise_cnt_var,psth_noise_cnt_var(targs(ss),:)));
end

all_data.rand_xcov = covar_rand;
all_data.psth_xcov = psth_xcov;
all_data.obs_xcov = obs_xcov;
all_data.spline_xcov = covar_spline_ZPT;
all_data.nan_locs = nan_locs;
all_data.tot_var = resp_var;
all_data.dir_noisevar = direct_noise_vars;
all_data.avg_rates = rpt_avg_rates;

%%
anal_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
if ~exist(anal_dir)
    mkdir(anal_dir)
end
cd(anal_dir);

sname = 'rpt_variability_analysis4';
sname = [sname sprintf('_ori%d',bar_ori)];

save(sname,'targs','Rpt_Data','ED_*','spline_*','tlags','rpt_ep*','poss_cnt_wins','all_data');