clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

Expt_num = 91;
Expt_name = sprintf('G%.3d',Expt_num);
data_dir = ['~/Data/bruce/' Expt_name];
% data_dir = ['/Volumes/james/Data/bruce/' Expt_name];
cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final/'];
% anal_dir = ['/Volumes/james/Analysis/bruce/' Expt_name '/ET_final/'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
% cluster_dir = ['/Volumes/james/Analysis/bruce/' Expt_name '/clustering'];

mod_data_name = 'dualori_eyecorr_mods2';
anal_name = 'dualori_eyecorr5';

recompute_init_mods = 0;
use_measured_pos = 0;
use_sac_kerns = 1;
use_LOOXV = 0;
use_coils = [0 0]; %[L R]

if any(use_coils > 0)
    anal_name = [anal_name '_Cprior'];
end
if use_measured_pos == 1
    mod_data_name = [mod_data_name '_Cinit'];
    anal_name = [anal_name '_Cinit'];
end

%dont fit stim models using these blocks
if Expt_num == 86
    ignore_blocks = [16 17 28 30]; %G086
elseif Expt_num == 87
    ignore_blocks = [15];
elseif Expt_num == 93
    ignore_blocks = [28];
else
    ignore_blocks = [];
end

%%
xv_frac = 0.2;

flen = 12;
use_nPix = 16;

n_fix_inf_it = 3; %3
n_drift_inf_it = 1; %3

fix_prior_sigma = sqrt(2*0.15^2);
fix_noise_sigma = sqrt(2*0.1^2);
drift_noise_sigma = sqrt(2*0.004^2);
drift_prior_sigma = sqrt(2*0.004^2); %.004 may be best here
drift_jump_sigma = sqrt(2*0.075^2); %0.05 start
drift_dsf = 3;

min_trial_dur = 0.75;

spatial_usfac = 2;

%%
eps = -1e3;

stim_fs = 100; %in Hz
dt = 0.01;
full_nPix = 36;
Fr = 1;

beg_buffer = 0.2;
end_buffer = 0.05;
trial_dur = 4;

n_probes = 96;

use_right_eye = false;

n_use_blocks = Inf;

use_nPix_us = use_nPix*spatial_usfac;
klen_us = use_nPix_us*flen;

sac_backlag = round(0.2/dt);
sac_forlag = round(0.3/dt);
sac_bincents = -sac_backlag:sac_forlag;
n_sac_bins = length(sac_bincents);

sp_dx = 0.0565/spatial_usfac;
max_shift = round(15*spatial_usfac);
% max_shift = round(6*spatial_usfac);
dshift = 1;

max_Dshift = round(8*spatial_usfac);
% max_Dshift = round(4*spatial_usfac);

%%
include_expts = {'rls.orRC'};
for ii = 1:length(Expts)
    expt_names{ii} = Expts{ii}.Header.expname;
    expt_dds(ii) = Expts{ii}.Stimvals.dd;
    expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
    expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
    expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
    expt_imback(ii) = isfield(Expts{ii}.Trials,'imi');
    included_type(ii) = any(strcmp(expt_names{ii},include_expts));
end
expt_has_ds(isnan(expt_has_ds)) = 0;
expt_has_ds(expt_has_ds == -1) = 0;
expt_binoc(isnan(expt_binoc)) = 0;
% cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori & expt_sac_dir == bar_ori);
cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1);


cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];

n_blocks = length(cur_block_set);

if length(cur_block_set) > n_use_blocks
    cur_block_set = cur_block_set(1:n_use_blocks);
end

n_blocks = length(cur_block_set);
%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_stim_times = [];
all_stim_mat = [];
all_stim_ori = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);

trial_toffset = zeros(length(cur_block_set),1);
cur_spkind_offset = 0;
cur_toffset = 0;
for ee = 1:n_blocks;
    fprintf('Expt %d Block %d of %d;  UNMATCHED EXPT TYPE\n',Expt_num,ee,n_blocks);
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
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)');
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    if strcmp(Expt_name,'G093')
        trial_wi = [Expts{cur_block}.Trials(:).wi];
        trial_wi = trial_wi(id_inds);
        all_trial_wi = cat(1,all_trial_wi,trial_wi(use_trials)');
    end
    
    fname = sprintf('%s/stims/Expt%d_stim',data_dir,cur_block);
    load(fname);
    buffer_pix = floor((expt_npix(cur_block) - full_nPix)/2);
    cur_use_pix = (1:full_nPix) + buffer_pix;
    
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start/1e4;
        cur_stim_ori = Expts{cur_block}.Trials(use_trials(tt)).or;
        n_frames = size(left_stim_mats{use_trials(tt)},1);
        if n_frames > 0
            if length(cur_stim_times) == 1
                cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
                cur_stim_times(cur_stim_times > trial_end_times(use_trials(tt))) = [];
            end
        end
        cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > min_trial_dur/dt
            use_frames = min(length(cur_stim_times),n_frames);
            cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            cur_stim_ori = cur_stim_ori(1:use_frames);
            %             bar_Xmat = create_time_embedding(cur_stim_mat,stim_params);
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times];
            all_t_axis = [all_t_axis; cur_t_axis];
            all_stim_mat = [all_stim_mat; cur_stim_mat];
            all_stim_ori = [all_stim_ori; cur_stim_ori];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
    end
    trial_cnt = trial_cnt + n_trials;
end

%%
full_nPix_us = spatial_usfac*full_nPix;
if spatial_usfac == 2
    all_stimmat_up = zeros(size(all_stim_mat,1),full_nPix_us);
    for ii = 1:size(all_stim_mat,2)
        all_stimmat_up(:,2*(ii-1)+1) = all_stim_mat(:,ii);
        all_stimmat_up(:,2*(ii-1)+2) = all_stim_mat(:,ii);
    end
elseif spatial_usfac == 1
    all_stimmat_up = all_stim_mat;
else
    error('Unsupported spatial Up-sample factor!');
end

hor_stim_inds = find(all_stim_ori == 0);
ver_stim_inds = find(all_stim_ori == 90);

full_stim_mat_ver = zeros(size(all_stim_mat));
full_stim_mat_hor = zeros(size(all_stim_mat));
full_stim_mat_hor(hor_stim_inds,:) = all_stim_mat(hor_stim_inds,:);
full_stim_mat_ver(ver_stim_inds,:) = all_stim_mat(ver_stim_inds,:);

full_stim_mat_ver_up = zeros(size(all_stimmat_up));
full_stim_mat_hor_up = zeros(size(all_stimmat_up));
full_stim_mat_hor_up(hor_stim_inds,:) = all_stimmat_up(hor_stim_inds,:);
full_stim_mat_ver_up(ver_stim_inds,:) = all_stimmat_up(ver_stim_inds,:);

stim_params = NIMcreate_stim_params([flen full_nPix],dt);
stim_params_us = NMMcreate_stim_params([flen full_nPix_us],dt);

all_hor_Xmat = create_time_embedding(full_stim_mat_hor,stim_params);
all_ver_Xmat = create_time_embedding(full_stim_mat_ver,stim_params);
all_Xmat = [all_ver_Xmat all_hor_Xmat];
clear all_ver_Xmat all_hor_Xmat

all_hor_Xmat_us = create_time_embedding(full_stim_mat_hor_up,stim_params_us);
all_ver_Xmat_us = create_time_embedding(full_stim_mat_ver_up,stim_params_us);
all_Xmat_us = [all_ver_Xmat_us all_hor_Xmat_us];

%% select submatrix with central pixels
[Xinds,Tinds] = meshgrid(1:full_nPix,1:flen);
Xinds = cat(2,Xinds,Xinds);
buffer_pix = floor((full_nPix - use_nPix)/2);
cur_use_pix = (1:use_nPix) + buffer_pix;
use_kInds = find(ismember(Xinds(:),cur_use_pix));

[Xinds_up,Tinds] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
Xinds_up = cat(2,Xinds_up,Xinds_up);
if spatial_usfac > 1
    use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1)-1/spatial_usfac & Xinds_up(:) <= cur_use_pix(end));
else
    use_kInds_up = use_kInds;
end
use_kInds_back = find(ismember(Xinds_up(use_kInds_up),cur_use_pix));

%% BIN SPIKES FOR MU AND SU
% LOAD REFCLUSTERS
fname = [cluster_dir '/final_cluster.mat'];
if exist(fname,'file')
    load(fname);
    SU_numbers = unique(SU_ID_mat(~isnan(SU_ID_mat)));
    for ii = 1:length(SU_numbers)
        SU_tot_nblocks = sum(SU_ID_mat(:) == SU_numbers(ii));
    end
    fprintf('%d SUs Clustered\n',length(SU_numbers));
    
else
    disp('No Cluster data found.');
end

%for SU probes
fprintf('Using %d SUs\n',length(SU_numbers));
all_binned_sua = nan(length(all_t_axis),length(SU_numbers));
cur_used_blocks = find(ismember(SU_target_blocks,cur_block_set));
SU_ID_mat = SU_ID_mat(cur_used_blocks,:);
[CC,BB] = meshgrid(1:length(SU_clust_data),1:length(cur_used_blocks));

su_probes = nan(1,length(SU_numbers));
all_su_spk_times = cell(length(SU_numbers),1);
all_su_spk_inds = cell(length(SU_numbers),1);
for ss = 1:length(SU_numbers)
    used_clust_set = unique(CC(SU_ID_mat==SU_numbers(ss))); %set of clusters used to capture this SU
    SU_block_probes(ss,:) = nan(1,length(cur_block_set));
    cur_su_spk_times = [];
    cur_su_spk_inds = [];
    cur_blocks = [];
    if ~isempty(used_clust_set)
        for cc = 1:length(used_clust_set)
            cur_clust = used_clust_set(cc);
            cur_probe = SU_clust_data(cur_clust).probe_num;
            cur_clust_label = SU_clust_data(cur_clust).cluster_label;
            cur_blocks = [cur_blocks find(SU_ID_mat(:,cur_clust) == SU_numbers(ss))];
            SU_block_probes(ss,cur_blocks) = cur_probe;
            
            all_su_inds = all_clust_ids{cur_probe} == cur_clust_label;
            cur_su_spk_times = all_spk_times{cur_probe}(all_su_inds);
            cur_su_spk_inds = all_spk_inds{cur_probe}(all_su_inds);
            spk_block_inds = round(interp1(all_t_axis,all_blockvec,cur_su_spk_times));
            cur_su_spk_times = cur_su_spk_times(ismember(spk_block_inds,cur_blocks));
            cur_su_spk_inds = cur_su_spk_inds(ismember(spk_block_inds,cur_blocks));
            
            all_su_spk_times{ss} = cat(1,all_su_spk_times{ss},cur_su_spk_times(:));
            all_su_spk_inds{ss} = cat(1,all_su_spk_inds{ss},cur_su_spk_inds(:));
        end
        
        cur_suahist = histc(all_su_spk_times{ss},all_t_bin_edges);
        cur_suahist(all_bin_edge_pts) = [];
        cur_id_set = ismember(all_blockvec,cur_blocks);
        all_binned_sua(cur_id_set,ss) = cur_suahist(cur_id_set);
        su_probes(ss) = mode(SU_block_probes(ss,~isnan(SU_block_probes(ss,:))));
    end
end

%for only-MU probes
all_binned_mua = nan(length(all_t_axis),n_probes);
for cc = 1:n_probes
    
    %this is the set of blocks where this probe had an SU, and the
    %correspodning SU numbers
    cur_set = find(SU_block_probes == cc);
    if ~isempty(cur_set)
        [cur_SS,cur_BB] = ind2sub([length(SU_numbers) length(cur_used_blocks)],cur_set);
    else
        cur_SS = [];
    end
    unique_su_nums = unique(cur_SS);
    cur_mua_inds = find(all_clust_ids{cc} >= 1);
    
    %remove spikes from isolated SUs from the MU
    for ss = 1:length(unique_su_nums)
        cur_mua_inds(ismember(all_spk_inds{cc}(cur_mua_inds),all_su_spk_inds{unique_su_nums(ss)})) = [];
    end
    
    [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc}(cur_mua_inds),all_t_bin_edges);
    cur_spkhist(all_bin_edge_pts) = [];
    all_binned_mua(:,cc) = cur_spkhist;
end


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
trial_toffset = zeros(length(cur_block_set),1);
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

%compute corrected eye data in bar-oriented frame
lin_correction = false;
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,zeros(length(cur_block_set),1),used_inds,lin_correction);

[saccades,et_params] = detect_saccades(corrected_eye_vals,all_eye_speed,all_eye_ts,et_params);

par_thresh = 1.25;
orth_thresh = 1.25;
[out_of_range] = detect_bad_fixation(corrected_eye_vals_interp,all_trialvec,used_inds,par_thresh,orth_thresh);
fract_out = length(out_of_range)/length(used_inds);
fprintf('Eliminating %.4f of data out of window\n',fract_out);
used_inds(ismember(used_inds,out_of_range)) = [];
NT = length(used_inds);

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
interp_sac_start_inds(bad_sacs) = [];
interp_sac_stop_inds(bad_sacs) = [];

%% CREATE SACCADE PREDICTOR MATS
saccade_start_inds = find(ismember(used_inds,interp_sac_start_inds));
used_saccade_set = find(ismember(interp_sac_start_inds,used_inds));
%nearest index in the used data set of the saccade stop time
saccade_stop_inds = round(interp1(used_inds,1:length(used_inds),interp_sac_stop_inds(used_saccade_set)))';

saccade_stop_inds(isnan(saccade_stop_inds)) = length(used_inds);

sac_amps = [saccades(:).amplitude];
is_micro = sac_amps(used_saccade_set) < 1;
big_sacs = find(~is_micro);
micro_sacs = find(is_micro);

saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));

if use_sac_kerns
    Xsac = zeros(NT,n_sac_bins);
    Xmsac = zeros(NT,n_sac_bins);
    for ii = 1:n_sac_bins
        cur_sac_target = saccade_start_inds(big_sacs) + sac_bincents(ii);
        uu = find(cur_sac_target > 1 & cur_sac_target < NT);
        cur_sac_target = cur_sac_target(uu);
        cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_trial_inds(big_sacs(uu))) = [];
        Xsac(cur_sac_target,ii) = 1;
        
        cur_sac_target = saccade_start_inds(micro_sacs) + sac_bincents(ii);
        uu = find(cur_sac_target > 1 & cur_sac_target < NT);
        cur_sac_target = cur_sac_target(uu);
        cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_trial_inds(micro_sacs(uu))) = [];
        Xmsac(cur_sac_target,ii) = 1;
    end
end

%% DEFINE FIXATION POINTS
trial_start_inds = [1; find(diff(all_trialvec(used_inds)) ~= 0) + 1];
trial_end_inds = [find(diff(all_trialvec(used_inds)) ~= 0); NT];

fix_start_inds = sort([trial_start_inds; saccade_stop_inds]);
fix_stop_inds = sort([trial_end_inds; saccade_start_inds]);
fix_durs = fix_stop_inds-fix_start_inds;
fix_start_inds(fix_durs<=0) = [];
fix_stop_inds(fix_durs <= 0) = [];
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

n_trials = length(trial_start_inds);
trial_ids = nan(NT,1);
for ii = 1:n_trials
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    trial_ids(cur_inds) = ii;
end


%% Create set of TR and XV trials
use_trials = unique(all_trialvec(used_inds));
nuse_trials = length(use_trials);

n_xv_trials = round(xv_frac*nuse_trials);
xv_trials = randperm(nuse_trials);
xv_trials(n_xv_trials+1:end) = [];
xv_trials = use_trials(xv_trials);
tr_trials = setdiff(use_trials,xv_trials);

tr_inds = find(ismember(all_trialvec(used_inds),tr_trials));
xv_inds = find(ismember(all_trialvec(used_inds),xv_trials));

full_inds = sort([tr_inds; xv_inds]);

all_Xmat = all_Xmat(used_inds,:);
all_Xmat_us = all_Xmat_us(used_inds,:);
all_ver_Xmat_us = all_ver_Xmat_us(used_inds,:);
all_hor_Xmat_us = all_hor_Xmat_us(used_inds,:);

%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITS
init_stim_params = NMMcreate_stim_params([flen 2*use_nPix],dt);
init_stim_params(2) = NMMcreate_stim_params([n_blocks 1],dt);

fin_stim_params(1) = NMMcreate_stim_params([flen 2*use_nPix_us],dt);
fin_stim_params(2) = NMMcreate_stim_params([n_blocks 1],dt);
if use_sac_kerns
    fin_stim_params(3) = NMMcreate_stim_params(n_sac_bins,dt);
    fin_stim_params(4) = NMMcreate_stim_params(n_sac_bins,dt);
end

null_stim_params = fin_stim_params(2:end);

block_L2 = 1;
silent = 1;
sac_d2t = 100;

base_lambda_custom = 50;
base_lambda_L1 = 10;

%create L2 mat for left eye
L2_params = create_L2_params([],[1 flen*use_nPix],[flen use_nPix],2,3,[Inf 0]);
L2_mat = generate_L2_mat(L2_params,2*flen*use_nPix);
%add L2 mat for right eye
L2_params = create_L2_params([],[flen*use_nPix+1 2*flen*use_nPix],[flen use_nPix],2,3,[Inf 0]);
L2_mat = L2_mat + generate_L2_mat(L2_params,2*flen*use_nPix);

%create L2 mat for left eye
L2_params = create_L2_params([],[1 flen*use_nPix_us],[flen use_nPix_us],2,3,[Inf 0]);
L2_mat_us = generate_L2_mat(L2_params,2*flen*use_nPix_us);
%add L2 mat for right eye
L2_params = create_L2_params([],[flen*use_nPix_us+1 2*flen*use_nPix_us],[flen use_nPix_us],2,3,[Inf 0]);
L2_mat_us = L2_mat_us + generate_L2_mat(L2_params,2*flen*use_nPix_us);

init_optim_p.optTol = 1e-5; init_optim_p.progTol = 1e-9;

sac_reg_params = NMMcreate_reg_params('lambda_d2T',sac_d2t);
if use_sac_kerns
    null_reg_params = NMMcreate_reg_params('lambda_d2T',[0; sac_d2t; sac_d2t],'lambda_L2',[block_L2; 0; 0]);
else
    null_reg_params = NMMcreate_reg_params('lambda_L2',block_L2);
end

n_squared_filts = 2;
mod_signs = ones(1,n_squared_filts+2);
%if using more than 2 quadratic filters, take the additional ones as
%suppressive
if n_squared_filts > 2
    mod_signs(n_squared_filts+1) = -1;
end
NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts) {'lin'}];
init_custom = [ones(n_squared_filts+1,1); 0;];
init_L2 = [zeros(n_squared_filts+1,1); block_L2];
init_reg_params = NMMcreate_reg_params('lambda_custom',init_custom,'lambda_L2',init_L2);
init_Xtargs = [ones(n_squared_filts+1,1); 2];

init_filts = cell(length(mod_signs),1);
cd(anal_dir);

if ~exist(['./' mod_data_name '.mat'],'file') || recompute_init_mods == 1
    tot_nUnits = length(su_probes) + n_probes;
    all_mod_SU = zeros(tot_nUnits,1);
    all_mos_SUnum = zeros(tot_nUnits,1);
    for ss = 1:n_probes;
        fprintf('Computing base LLs for MU %d of %d\n',ss,n_probes);
        cur_tr_inds = tr_inds(~isnan(all_binned_mua(used_inds(tr_inds),ss)));
        cur_xv_inds = xv_inds(~isnan(all_binned_mua(used_inds(xv_inds),ss)));
        tr_NT = length(cur_tr_inds);
        xv_NT = length(cur_xv_inds);
        if ~isempty(cur_tr_inds)
            Robs = all_binned_mua(used_inds(cur_tr_inds),ss);
            Robsxv = all_binned_mua(used_inds(cur_xv_inds),ss);
            
            tr_X{1} = all_Xmat((cur_tr_inds),use_kInds);
            tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
            xv_X{1} = all_Xmat((cur_xv_inds),use_kInds);
            xv_X{2} = Xblock(used_inds(cur_xv_inds),:);
            
            if use_sac_kerns
                tr_X{3} = Xsac(cur_tr_inds,:);
                tr_X{4} = Xmsac(cur_tr_inds,:);
                xv_X{3} = Xsac(cur_xv_inds,:);
                xv_X{4} = Xmsac(cur_xv_inds,:);
                null_mod = NMMinitialize_model(null_stim_params,[1 1 1],repmat({'lin'},1,3),null_reg_params,[1 2 3]);
            else
                null_mod = NMMinitialize_model(null_stim_params,[1],{'lin'},null_reg_params,[1]);
            end
            
            null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));
            
            gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
            gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent,init_optim_p,L2_mat);
            
            tr_X{1} = all_Xmat_us((cur_tr_inds),use_kInds_up);
            xv_X{1} = all_Xmat_us((cur_xv_inds),use_kInds_up);
            %spatial up-sampling of filter estimates
            base_filts = reshape([gqm1.mods(find(init_Xtargs == 1)).filtK],[flen use_nPix 2 n_squared_filts+1]);
            if spatial_usfac == 2
                base_filts_up = zeros(flen,use_nPix_us,2,n_squared_filts+1);
                for ii = 1:use_nPix
                    base_filts_up(:,2*(ii-1)+1,:,:) = 0.5*base_filts(:,ii,:,:);
                    base_filts_up(:,2*(ii-1)+2,:,:) = 0.5*base_filts(:,ii,:,:);
                end
            elseif spatial_usfac == 1
                base_filts_up = base_filts;
            else
                error('unsupported')
            end
            base_filts_up = reshape(base_filts_up,2*use_nPix_us*flen,n_squared_filts+1);
            
            init_filts{end} = gqm1.mods(find(init_Xtargs==2)).filtK;
            for ii = 1:n_squared_filts+1
                init_filts{ii} = base_filts_up(:,ii);
            end
            gqm2 = NMMinitialize_model(fin_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs,init_filts);
            gqm2.spk_NL_params(1) = gqm1.spk_NL_params(1);
            [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm2,Robs,tr_X);
            gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_custom',base_lambda_custom./var(gint)');
            gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
            if use_sac_kerns
                gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,3,zeros(n_sac_bins,1),sac_reg_params); %sac term
                gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,4,zeros(n_sac_bins,1),sac_reg_params); %nsac term
            end
            gqm2 = NMMfit_filters(gqm2,Robs,tr_X,[],[],silent,[],L2_mat_us);
            
            xvLL = NMMmodel_eval(gqm2,Robsxv,xv_X);
            null_xvLL(ss) = NMMmodel_eval(null_mod,Robsxv,xv_X(2:end));
            all_mod_xvLLimp(ss) = (xvLL-null_xvLL(ss))/log(2);
            
            %now refit model using all (usable) data
            if xv_frac > 0
                cur_tr_inds = sort([cur_tr_inds; cur_xv_inds]);
                Robs = all_binned_mua(used_inds(cur_tr_inds),ss);
                tr_X{1} = all_Xmat_us((cur_tr_inds),use_kInds_up);
                tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
                if use_sac_kerns
                    tr_X{3} = Xsac(cur_tr_inds,:);
                    tr_X{4} = Xmsac(cur_tr_inds,:);
                end
            end
            null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));
            all_nullmod(ss) = null_mod;
            
            gqm2 = NMMfit_filters(gqm2,Robs,tr_X,[],[],silent,[],L2_mat_us);
            all_mod_fits(ss) = gqm2;
            all_mod_fits_withspkNL(ss) = NMMfit_logexp_spkNL(gqm2,Robs,tr_X);
            
            [LL,~,pred_rate] = NMMmodel_eval(all_mod_fits_withspkNL(ss),Robs,tr_X);
            [null_LL(ss),~,null_prate] = NMMmodel_eval(null_mod,Robs,tr_X(2:end));
            [all_mod_R2(ss),all_mod_dev(ss),all_null_dev(ss)] = pseudo_r2(Robs,pred_rate,null_prate);
            all_mod_LLimp(ss) = (LL-null_LL(ss))/log(2);
            
            if xv_frac > 0
                fprintf('xvLL imp %.4f\n',all_mod_xvLLimp(ss));
            else
                fprintf('LL imp %.4f\n',all_mod_LLimp(ss));
            end
        end
    end
    
    for ss = 1:length(su_probes);
        fprintf('Computing base LLs for SU %d of %d\n',ss,length(su_probes));
        all_mod_SU(ss+n_probes) = su_probes(ss);
        all_mod_SUnum(ss+n_probes) = SU_numbers(ss);
        cur_tr_inds = tr_inds(~isnan(all_binned_sua(used_inds(tr_inds),ss)));
        cur_xv_inds = xv_inds(~isnan(all_binned_sua(used_inds(xv_inds),ss)));
        tr_NT = length(cur_tr_inds);
        xv_NT = length(cur_xv_inds);
        if ~isempty(cur_tr_inds)
            Robs = all_binned_sua(used_inds(cur_tr_inds),ss);
            Robsxv = all_binned_sua(used_inds(cur_xv_inds),ss);
            
            tr_X{1} = all_Xmat((cur_tr_inds),use_kInds);
            tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
            xv_X{1} = all_Xmat((cur_xv_inds),use_kInds);
            xv_X{2} = Xblock(used_inds(cur_xv_inds),:);
            if use_sac_kerns
                tr_X{3} = Xsac(cur_tr_inds,:);
                tr_X{4} = Xmsac(cur_tr_inds,:);
                xv_X{3} = Xsac(cur_xv_inds,:);
                xv_X{4} = Xmsac(cur_xv_inds,:);
                null_mod = NMMinitialize_model(null_stim_params,[1 1 1],repmat({'lin'},1,3),null_reg_params,[1 2 3]);
            else
                null_mod = NMMinitialize_model(null_stim_params,[1],{'lin'},null_reg_params,[1]);
            end
            null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));
            
            gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
            gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent,init_optim_p,L2_mat);
            
            tr_X{1} = all_Xmat_us((cur_tr_inds),use_kInds_up);
            xv_X{1} = all_Xmat_us((cur_xv_inds),use_kInds_up);
            %spatial up-sampling of filter estimates
            base_filts = reshape([gqm1.mods(find(init_Xtargs == 1)).filtK],[flen use_nPix 2 n_squared_filts+1]);
            if spatial_usfac == 2
                base_filts_up = zeros(flen,use_nPix_us,2,n_squared_filts+1);
                for ii = 1:use_nPix
                    base_filts_up(:,2*(ii-1)+1,:,:) = 0.5*base_filts(:,ii,:,:);
                    base_filts_up(:,2*(ii-1)+2,:,:) = 0.5*base_filts(:,ii,:,:);
                end
            elseif spatial_usfac == 1
                base_filts_up = base_filts;
            else
                error('unsupported')
            end
            base_filts_up = reshape(base_filts_up,2*use_nPix_us*flen,n_squared_filts+1);
            
            init_filts{end} = gqm1.mods(find(init_Xtargs==2)).filtK;
            for ii = 1:n_squared_filts+1
                init_filts{ii} = base_filts_up(:,ii);
            end
            gqm2 = NMMinitialize_model(fin_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs,init_filts);
            gqm2.spk_NL_params(1) = gqm1.spk_NL_params(1);
            [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm2,Robs,tr_X);
            gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_custom',base_lambda_custom./var(gint)');
            gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
            if use_sac_kerns
                gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,3,zeros(n_sac_bins,1),sac_reg_params); %sac term
                gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,4,zeros(n_sac_bins,1),sac_reg_params); %nsac term
            end
            gqm2 = NMMfit_filters(gqm2,Robs,tr_X,[],[],silent,[],L2_mat_us);
            
            xvLL = NMMmodel_eval(gqm2,Robsxv,xv_X);
            null_xvLL(ss+n_probes) = NMMmodel_eval(null_mod,Robsxv,xv_X(2:end));
            all_mod_xvLLimp(ss+n_probes) = (xvLL-null_xvLL(ss+n_probes))/log(2);
            
            %now refit model using all (usable) data
            if xv_frac > 0
                cur_tr_inds = sort([cur_tr_inds; cur_xv_inds]);
                Robs = all_binned_sua(used_inds(cur_tr_inds),ss);
                tr_X{1} = all_Xmat_us((cur_tr_inds),use_kInds_up);
                tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
                if use_sac_kerns
                    tr_X{3} = Xsac(cur_tr_inds,:);
                    tr_X{4} = Xmsac(cur_tr_inds,:);
                end
            end
            null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));
            all_nullmod(ss+n_probes) = null_mod;
            
            gqm2 = NMMfit_filters(gqm2,Robs,tr_X,[],[],silent,[],L2_mat_us);
            all_mod_fits(ss+n_probes) = gqm2;
            all_mod_fits_withspkNL(ss+n_probes) = NMMfit_logexp_spkNL(gqm2,Robs,tr_X);
            
            [LL,~,pred_rate] = NMMmodel_eval(all_mod_fits_withspkNL(ss+n_probes),Robs,tr_X);
            [null_LL(ss+n_probes),~,null_prate] = NMMmodel_eval(null_mod,Robs,tr_X(2:end));
            [all_mod_R2(ss+n_probes),all_mod_dev(ss+n_probes),all_null_dev(ss+n_probes)] = pseudo_r2(Robs,pred_rate,null_prate);
            all_mod_LLimp(ss+n_probes) = (LL-null_LL(ss+n_probes))/log(2);
            if xv_frac > 0
                fprintf('xvLL imp %.4f\n',all_mod_xvLLimp(ss+n_probes));
            else
                fprintf('LL imp %.4f\n',all_mod_LLimp(ss+n_probes));
            end
        end
    end
    save(mod_data_name,'all_mod*','all_nullmod','su_probes','null_xvLL','null_LL','*_trials');
    clear tr_X xv_X
else
    fprintf('Loading pre-computed initial models\n');
    load(mod_data_name);
end

%% SELECT USABLE UNITS AND make Robs_mat
if xv_frac == 0
    LL_imp_thresh = 5e-3;
    usable_units = find(all_mod_LLimp >= LL_imp_thresh);
else
    LL_imp_thresh = 0;
    usable_units = find(all_mod_xvLLimp > LL_imp_thresh);
end

n_used_sus = sum(all_mod_SU(usable_units) ~= 0);
n_used_mus = sum(all_mod_SU(usable_units) == 0);
fprintf('Using %d SUs and %d MUs for analysis\n',n_used_sus,n_used_mus);
tr_set = usable_units;
n_tr_chs = length(tr_set);

% make Robs_mat
poss_blocks = 1:(n_blocks-1);
Robs_mat = nan(length(used_inds),n_tr_chs);
for ss = 1:n_tr_chs
    if all_mod_SU(tr_set(ss)) > 0
        su_probe_ind = find(SU_numbers == all_mod_SUnum(tr_set(ss)));
        Robs_mat(:,ss) = all_binned_sua(used_inds,su_probe_ind);
    else
        Robs_mat(:,ss) = all_binned_mua(used_inds,tr_set(ss));
    end
end

%don't use separate xv set for eye-tracking
tr_inds = full_inds;

%%
[Xinds_up,Tinds] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
use_kInds_up_single = find(Xinds_up(:) >= cur_use_pix(1)-1/spatial_usfac & Xinds_up(:) <= cur_use_pix(end));

x_shifts = -max_shift:dshift:max_shift;
y_shifts = -max_shift:dshift:max_shift;
n_xshifts = length(x_shifts);

[Xsh,Ysh] = meshgrid(x_shifts,y_shifts);
[Xsh_inds,Ysh_inds] = meshgrid(1:length(x_shifts),1:length(y_shifts));
SH = [Xsh(:) Ysh(:)];
SH_inds = [Xsh_inds(:) Ysh_inds(:)];
n_shifts = size(SH,1);
zero_frame = find(SH(:,1) == 0 & SH(:,2) == 0);

%generate shift matrices. Must be applied to the stimulus (not the filters)
It = speye(flen);
shift_mat = cell(n_xshifts,1);
for xx = 1:n_xshifts
    temp = spdiags( ones(full_nPix_us,1), -x_shifts(xx), full_nPix_us, full_nPix_us);
    temp = kron(temp,It);
    shift_mat{xx} = temp(:,use_kInds_up_single);
end

%shift matrices for images (non-time-embedded)
max_Tshift = max_shift;
Tshifts = -max_Tshift:dshift:max_Tshift;
n_Tshifts = length(Tshifts);
Tshift_mat = cell(n_Tshifts,1);
for xx = 1:n_Tshifts
    Tshift_mat{xx} = spdiags( ones(full_nPix_us,1), -Tshifts(xx), full_nPix_us, full_nPix_us);
end

D_shifts = -max_Dshift:dshift:max_Dshift;

[DXsh,DYsh] = meshgrid(D_shifts,D_shifts);
[DXsh_inds,DYsh_inds] = meshgrid(1:length(D_shifts),1:length(D_shifts));
DSH = [DXsh(:) DYsh(:)];
DSH_inds = [DXsh_inds(:) DYsh_inds(:)];
n_Dshifts = size(DSH,1);
%generate shift matrices. Must be applied to the stimulus (not the filters)
It = speye(flen);
Dshift_mat = cell(length(D_shifts),1);
for xx = 1:length(D_shifts)
    temp = spdiags( ones(full_nPix_us,1), -D_shifts(xx), full_nPix_us, full_nPix_us);
    temp = kron(temp,It);
    Dshift_mat{xx} = temp(:,use_kInds_up_single);
end


%%
fix_prior = -sum((SH*sp_dx).^2,2)/(2*fix_prior_sigma^2);
fix_prior = bsxfun(@minus,fix_prior,logsumexp(fix_prior)); %normalize

drift_jump_prior = -sum((DSH*sp_dx).^2,2)/(2*(drift_jump_sigma)^2);
drift_jump_prior = bsxfun(@minus,drift_jump_prior,logsumexp(drift_jump_prior)); %normalize

base_lA = pdist2(DSH*sp_dx,DSH*sp_dx);
base_lA = -base_lA.^2/(2*(drift_prior_sigma*drift_dsf)^2);
base_lA = bsxfun(@minus,base_lA,logsumexp(base_lA)); %normalize

%%
warning('OFF','stats:statrobustfit:IterationLimit');

measured_eyepos = corrected_eye_vals_interp;
max_sim_pos = max_shift*sp_dx;
% measured_eyepos = tanh(measured_eyepos/max_sim_pos)*max_sim_pos;

measured_eyepos(isnan(measured_eyepos)) = 0;

%smooth out fast transients in eye signal
eye_smooth_sig = round(0.025/dt);
interp_inds = [];
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > eye_smooth_sig*5;
        measured_eyepos(used_inds(cur_inds),1) = jmm_smooth_1d_cor(measured_eyepos(used_inds(cur_inds),1),eye_smooth_sig,2);
        measured_eyepos(used_inds(cur_inds),2) = jmm_smooth_1d_cor(measured_eyepos(used_inds(cur_inds),2),eye_smooth_sig,2);
        measured_eyepos(used_inds(cur_inds),3) = jmm_smooth_1d_cor(measured_eyepos(used_inds(cur_inds),3),eye_smooth_sig,2);
        measured_eyepos(used_inds(cur_inds),4) = jmm_smooth_1d_cor(measured_eyepos(used_inds(cur_inds),4),eye_smooth_sig,2);
    end
    interp_inds = [interp_inds; cur_inds'];
end
interp_inds = unique(interp_inds);
measured_eyepos(used_inds,:) = interp1(used_inds(interp_inds),measured_eyepos(used_inds(interp_inds),:),used_inds);
measured_eyepos(isnan(measured_eyepos)) = 0;

measured_eyepos = measured_eyepos(used_inds,:);

measured_drift = nan(length(used_inds),4);
measured_fix_avg = nan(n_fixs,4);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    measured_fix_avg(ii,:) = median(measured_eyepos(cur_inds,:));
    measured_drift(cur_inds,:) = bsxfun(@minus,measured_eyepos(cur_inds,:),measured_fix_avg(ii,:));
end

back_measured_drift = nan(size(measured_drift));
%back-project measured drift
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    cur_pinds = cur_inds + sac_shift;
    if length(cur_inds) > sac_shift
        back_measured_drift(cur_inds(sac_shift+1:end),:) = measured_drift(cur_inds(1:(end-sac_shift)),:);
    end
end
measured_drift = back_measured_drift;
clear back_measured_drift

measured_drift = [zeros(1,4); diff(measured_drift)];

%if using both coils
if use_coils(1) == 1 && use_coils(2) == 1
    
    measured_fix_deltas = nan(n_fixs,2);
    measured_fix_deltas(2:end,1) = mean(diff(measured_fix_avg(:,[1 3])),2);
    measured_fix_deltas(2:end,2) = mean(diff(measured_fix_avg(:,[2 4])),2);
    %     fix_noise_sigma = fix_noise_sigma/sqrt(2); %noise variance decreases by factor of sqrt(2) with 2 independent measures
    %     post_drift_var = 1/(2/drift_noise_sigma^2 + 1/drift_prior_sigma^2);
    %     post_mean_drift = post_drift_var*2*mean(measured_drift,2)/drift_noise_sigma^2;
    
    %if using left coil only
elseif use_coils(1) == 1
    measured_fix_deltas = nan(n_fixs,2);
    measured_fix_deltas(2:end,1) = diff(measured_fix_avg(:,1));
    measured_fix_deltas(2:end,2) = diff(measured_fix_avg(:,2));
    
    post_drift_var = 1/(1/drift_noise_sigma^2 + 1/drift_prior_sigma^2);
    post_mean_drift(:,1) = post_drift_var*measured_drift(:,1)/drift_noise_sigma^2;
    post_mean_drift(:,2) = post_drift_var*measured_drift(:,2)/drift_noise_sigma^2;
    
    %if using right coil only
elseif use_coils(2) == 1
    %     measured_fix_deltas = nan(n_fixs,1);
    %     measured_fix_deltas(2:end) = diff(measured_fix_avg(:,2));
    %     post_drift_var = 1/(1/drift_noise_sigma^2 + 1/drift_prior_sigma^2);
    %     post_mean_drift = post_drift_var*measured_drift(:,2)/drift_noise_sigma^2;
    
else
    post_drift_var = drift_prior_sigma^2;
    post_mean_drift = zeros(NT,2);
    measured_fix_deltas = nan(n_fixs,2);
end

measured_fix_deltas(isnan(measured_fix_deltas)) = 0;
post_drift_sigma = sqrt(post_drift_var);
% post_mean_drift(isnan(post_mean_drift)) = 0;


%%
su_inds = find(all_mod_SU(tr_set) > 0);
cd(anal_dir);
load(anal_name)
%% ITERATE FIXATION-BASED CORRECTIONS


%back-project saccade-times
all_fix_post_mean_cor = nan(NT,2);
all_fix_post_mean_cor(~isnan(fix_ids),:) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)),:);
all_fix_post_mean_cor = interp1(find(~isnan(fix_ids)),all_fix_post_mean_cor(~isnan(fix_ids),:),1:NT);

%%

all_fix_post_mean_cor = round(all_fix_post_mean_cor);

all_shift_stimmat_ver = full_stim_mat_ver_up;
all_shift_stimmat_hor = full_stim_mat_hor_up;
for i=1:NT
    if all_stim_ori(used_inds(i)) == 0
        %             all_shift_stimmat_hor(used_inds(i),:) = all_shift_stimmat_hor(used_inds(i),:)*Tshift_mat{all_fix_post_mean_inds(i,2)};
        all_shift_stimmat_hor(used_inds(i),:) = shift_matrix_Nd(all_shift_stimmat_hor(used_inds(i),:),-all_fix_post_mean_cor(i,2),2);
    else
        %             all_shift_stimmat_ver(used_inds(i),:) = all_shift_stimmat_ver(used_inds(i),:)*Tshift_mat{all_fix_post_mean_inds(i,1)};
        all_shift_stimmat_ver(used_inds(i),:) = shift_matrix_Nd(all_shift_stimmat_ver(used_inds(i),:),-all_fix_post_mean_cor(i,1),2);
    end
end

shift_hor_Xmat_us = create_time_embedding(all_shift_stimmat_hor,stim_params_us);
shift_ver_Xmat_us = create_time_embedding(all_shift_stimmat_ver,stim_params_us);
cur_X{1} = [shift_ver_Xmat_us shift_hor_Xmat_us];
cur_X{1} = cur_X{1}(used_inds,use_kInds_up);
clear shift_hor_Xmat_us shift_ver_Xmat_us

cur_X{2} = Xblock(used_inds,:);
if use_sac_kerns
    cur_X{3} = Xsac;
    cur_X{4} = Xmsac;
end

%%
all_ori_X = zeros(length(all_t_axis),2);
all_ori_X(all_stim_ori == 0,1) = 1;
all_ori_X(all_stim_ori == 90,2) = 1;
ori_params = NMMcreate_stim_params([flen 2]);
all_ori_Xmat = create_time_embedding(all_ori_X,ori_params);
all_ori_Xmat = all_ori_Xmat(used_inds,:);

%%


for cc = 1:length(tr_set)
    cur_Robs = Robs_mat(:,cc);
    init_mod = NMMinitialize_model(ori_params,[1],{'lin'});
    init_mod = NMMfit_filters(init_mod,cur_Robs,all_ori_Xmat);
    ori_mod(cc) = init_mod;
    [ori_LL,~,~,~,~,~,nullLL] = NMMmodel_eval(init_mod,cur_Robs,all_ori_Xmat);
    ori_LLimp(cc) = (ori_LL-nullLL)/log(2);
end

%%
for cc = 1:length(tr_set);
    cc
    cur_Robs = Robs_mat(:,cc);
    [ori_LL,~,~,~,~,fgint,nullLL] = NMMmodel_eval(ori_mod(cc),cur_Robs,all_ori_Xmat);
    g_tot = fgint;
    any_sac_inds = find(any(Xmsac > 0,2));
    lambda_d2T = 10;
    lambda_L2 = 5;
    sac_reg_params = NMMcreate_reg_params('lambda_d2T',lambda_d2T,'lambda_L2',lambda_L2,'boundary_conds',[0 0 0]);
    
    Xsac_tot = bsxfun(@times,Xmsac,g_tot);
    clear tr_stim
    tr_stim{1} = [g_tot];
    tr_stim{2} = Xmsac;
    tr_stim{3} = Xsac_tot;
    clear sac_stim_params
    sac_stim_params(1) = NMMcreate_stim_params(1);
    sac_stim_params(2:3) = NMMcreate_stim_params([size(Xsac_tot,2)]);
    mod_signs = [1 1 1];
    Xtargets = [1 2 3];
    NL_types = {'lin','lin','lin'};
    post_gsac_Smod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
    post_gsac_Smod.mods(1).reg_params = NMMcreate_reg_params();
    post_gsac_Smod = NMMfit_filters(post_gsac_Smod,cur_Robs(any_sac_inds),get_Xcell_tInds(tr_stim,any_sac_inds),[],[],silent);
    post_gsac_Smod_LL = NMMmodel_eval(post_gsac_Smod,cur_Robs(any_sac_inds),get_Xcell_tInds(tr_stim,any_sac_inds));
    [post_Smod_LL,~,post_Smod_predrate,~,~,~,nullLL] = NMMmodel_eval( post_gsac_Smod, cur_Robs, tr_stim);
    sacStimProc(cc).gsac_post_singmod = post_gsac_Smod;
    ori_gainkern(cc,:) = post_gsac_Smod.mods(3).filtK;
    ori_offkern(cc,:) = post_gsac_Smod.mods(2).filtK;
    
%     post_gsac_Smod.mods(3).filtK(:) = 0;
%     rrange = linspace(0.5,1.2,50);
    for ii = 1:n_sac_bins
       test_Xmsac = zeros(size(Xmsac));
       test_Xmsac(:,ii) = 1;
        test_Xsac_tot = bsxfun(@times,test_Xmsac,g_tot);
       tr_stim{2} = test_Xmsac;
       tr_stim{3} = test_Xsac_tot;
       [~,~,temp_predrate] = NMMmodel_eval(post_gsac_Smod,cur_Robs,tr_stim);
       aprate = mean(temp_predrate);
       temp_info(ii) = mean(temp_predrate.*log2(temp_predrate/aprate))/aprate;
%        temp_rhist(ii,:) = hist(temp_predrate,rrange);
    end
    sac_ori_newinfo(cc,:) = temp_info;
    
    %%
                fprintf('Estimating tent-basis model\n');
            n_Gbins = 15;
            Xtick = -(sac_backlag-1/2):(1):(sac_forlag+1/2);
            n_sbins = length(Xtick);
            
            cur_sac_starts = saccade_start_inds(micro_sacs);
            cur_sac_stops = saccade_stop_inds(micro_sacs);
            t_since_sac_start = nan(NT,1);
            for ii = 1:length(cur_sac_starts)
                prev_tstart = find(trial_start_inds <= cur_sac_starts(ii),1,'last');
                next_tstop = find(trial_end_inds >= cur_sac_starts(ii),1,'first');
                cur_inds = (cur_sac_starts(ii) - sac_backlag):(cur_sac_starts(ii) + sac_forlag);
                cur_uset = find(cur_inds > trial_start_inds(prev_tstart) & cur_inds < trial_end_inds(next_tstop));
                t_since_sac_start(cur_inds(cur_uset)) = sac_bincents(cur_uset);
            end
            
            TB_lambda = 10;
            TB_stim = [t_since_sac_start g_tot];
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
            
            gsac_ov_TB_info(cc) = sum(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate))/cur_avg_rate;
            
            gsac_TB_avg_rate(cc,:) = marg_gsacrate;
            gsac_TB_info(cc,:) = gsacdep_info./marg_gsacrate;
            sacStimProc(cc).gsac_TB_gdist = marg_gdist;
            sacStimProc(cc).gsac_TB_grate = marg_grate;
            
            
    %%
%     filts = reshape([ori_mod(cc).mods(1).filtK],[flen 2]);
%     filt_tkern = squeeze(std(filts,[],2));
%     [~,tloc] = max(filt_tkern);
%     [~,Tinds] = meshgrid(1:2,1:flen);
%     cur_ori_stim = all_ori_Xmat(:,Tinds(:) == tloc);
    for ii = 1:n_sac_bins
        temp = find(Xmsac(:,ii) == 1);
        sac_avgrate(cc,ii) = mean(cur_Robs(temp));
        cur_avg_rate = sac_avgrate(ii)*ones(size(temp));
        sac_nullLL(ii) = nansum(cur_Robs(temp).*log2(sac_avgrate(ii)) - sac_avgrate(ii));
        sac_Nspks(ii) = sum(cur_Robs(temp));
        
        sac_spost_LL(cc,ii) = nansum(cur_Robs(temp).*log2(post_Smod_predrate(temp)) - post_Smod_predrate(temp));
        sac_spost_info(cc,ii) = nanmean(post_Smod_predrate(temp).*log2(post_Smod_predrate(temp)/mean(post_Smod_predrate(temp))))/mean(post_Smod_predrate(temp));
        sac_spost_trigavg(cc,ii) = mean(post_Smod_predrate(temp));
        
%         set1 = temp(cur_ori_stim(temp,1) == 1);
%         set2 = temp(cur_ori_stim(temp,2) == 1);
%         r1 = mean(cur_Robs(set1));
%         r2 = mean(cur_Robs(set2));
%         ovr = mean(cur_Robs(temp));
%         a_r1(ii) = r1; a_r2(ii) = r2; a_ovr(ii) = ovr;
%         p1 = length(set1)/length(temp);
%         p2 = 1-p1;
%         brute_info(cc,ii) = p1*r1/ovr*log2(r1/ovr) + p2*r2/ovr*log2(r2/ovr);
    end
    gsac_spost_ov_modinfo(cc) = mean(post_Smod_predrate(any_sac_inds)/mean(post_Smod_predrate(any_sac_inds)).*log2(post_Smod_predrate(any_sac_inds)/mean(post_Smod_predrate(any_sac_inds))));
    
end

gsac_ori_info = bsxfun(@rdivide,sac_spost_info,gsac_spost_ov_modinfo');
new_avg_ori_info = mean(sac_ori_newinfo,2);
gsac_newori_info = bsxfun(@rdivide,sac_ori_newinfo,new_avg_ori_info);
gsac_TB_rinfo = bsxfun(@rdivide,gsac_TB_info,gsac_ov_TB_info');
%%
for cc = 1:length(tr_set);
    cc
    cur_Robs = Robs_mat(:,cc);
    cur_mod = it_mods{end}(tr_set(cc));
    cur_mod.mods(4:end) = [];
    %     [ori_LL,~,~,~,~,fgint,nullLL] = NMMmodel_eval(ori_mod(cc),cur_Robs,all_ori_Xmat);
    %     g_tot = fgint;
    [~,~,~,~,~,fgint,nullLL] = NMMmodel_eval(cur_mod,cur_Robs,cur_X);
    g_tot = sum(fgint,2);
    
    any_sac_inds = find(any(Xmsac > 0,2));
    lambda_d2T = 10;
    lambda_L2 = 5;
    sac_reg_params = NMMcreate_reg_params('lambda_d2T',lambda_d2T,'lambda_L2',lambda_L2,'boundary_conds',[0 0 0]);
    
    Xsac_tot = bsxfun(@times,Xmsac,g_tot);
    clear tr_stim
    tr_stim{1} = [g_tot];
    tr_stim{2} = Xmsac;
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
    post_gsac_Smod_LL = NMMmodel_eval(post_gsac_Smod,cur_Robs(any_sac_inds),get_Xcell_tInds(tr_stim,any_sac_inds));
    [post_Smod_LL,~,post_Smod_predrate,~,~,~,nullLL] = NMMmodel_eval( post_gsac_Smod, cur_Robs, tr_stim);
    sacStimProc(cc).gsac_post_Xmod = post_gsac_Smod;
    
    %%
    for ii = 1:n_sac_bins
        temp = find(Xmsac(:,ii) == 1);
        sac_avgrate(cc,ii) = mean(cur_Robs(temp));
        cur_avg_rate = sac_avgrate(ii)*ones(size(temp));
        sac_nullLL(ii) = nansum(cur_Robs(temp).*log2(sac_avgrate(ii)) - sac_avgrate(ii));
        sac_Nspks(ii) = sum(cur_Robs(temp));
        
        temp2 = find(Xmsac(any_sac_inds,ii) == 1);
        
        sac_spost_XLL(cc,ii) = nansum(cur_Robs(temp).*log2(post_Smod_predrate(temp)) - post_Smod_predrate(temp));
        sac_spost_Xinfo(cc,ii) = nanmean(post_Smod_predrate(temp).*log2(post_Smod_predrate(temp)/mean(post_Smod_predrate(temp))))/mean(post_Smod_predrate(temp));
        
        
    end
    gsac_spost_ov_Xmodinfo(cc) = mean(post_Smod_predrate(any_sac_inds)/mean(post_Smod_predrate(any_sac_inds)).*log2(post_Smod_predrate(any_sac_inds)/mean(post_Smod_predrate(any_sac_inds))));
    
    %%
    
end

gsac_stim_info = bsxfun(@rdivide,sac_spost_Xinfo,gsac_spost_ov_Xmodinfo');

%%
ov_avg_rates = mean(Robs_mat);
rel_sactrg_rate = bsxfun(@rdivide,sac_avgrate,ov_avg_rates');

%%
figure
for cc = 1:94
    subplot(3,1,1)
    plot(sac_bincents*dt,ori_offkern(cc,:))
    hold on
    plot(sac_bincents*dt,ori_gainkern(cc,:),'r')
    subplot(3,1,2)
    plot(sac_bincents*dt,sac_spost_info(cc,:))
    subplot(3,1,3)
    plot(sac_bincents*dt,sac_avgrate(cc,:))
    pause
    clf
end
    