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
old_anal_name = 'dualori_eyecorr5';
anal_name = 'dualori_eyecorr_highres';
hrmod_data_name = 'dualori_eyecorr_mods_hres';

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

n_drift_inf_it = 1; %3

drift_noise_sigma = sqrt(2*0.004^2);
drift_prior_sigma = sqrt(2*0.004^2); 
drift_jump_sigma = sqrt(2*0.075^2); 
% drift_noise_sigma = 0.003;
% drift_prior_sigma = 0.004; %.004 may be best here
% drift_jump_sigma = 0.075; %0.05 start
drift_dsf = 2;

min_trial_dur = 0.75;

spatial_usfac = 4;
old_usfac = 2;

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

sac_backlag = round(0.05/dt);
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
if spatial_usfac > 1
    all_stimmat_up = zeros(size(all_stim_mat,1),full_nPix_us);
    for ii = 1:size(all_stim_mat,2)
        for jj = 1:spatial_usfac
        all_stimmat_up(:,spatial_usfac*(ii-1)+jj) = all_stim_mat(:,ii);
        end
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


%% select submatrix with central pixels
[Xinds,Tinds] = meshgrid(1:full_nPix,1:flen);
Xinds = cat(2,Xinds,Xinds);
buffer_pix = floor((full_nPix - use_nPix)/2);
cur_use_pix = (1:use_nPix) + buffer_pix;
use_kInds = find(ismember(Xinds(:),cur_use_pix));

[Xinds_up,Tinds] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
Xinds_up = cat(2,Xinds_up,Xinds_up);
if spatial_usfac > 1
    cnt = 0;
    use_kInds_up = [];
    while length(use_kInds_up) < use_nPix_us*2*flen
        use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1)-cnt/spatial_usfac & Xinds_up(:) <= cur_use_pix(end)+cnt/spatial_usfac);
        if length(use_kInds_up) < use_nPix_us*2*flen
            use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1)-(cnt+1)/spatial_usfac & Xinds_up(:) <= cur_use_pix(end) + cnt/spatial_usfac);
        end
            cnt = cnt + 1;
    end
else
    use_kInds_up = use_kInds;
end
use_kInds_back = find(ismember(Xinds_up(use_kInds_up),cur_use_pix));

[Xinds_ups,Tinds] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
if spatial_usfac > 1
    cnt = 0;
    use_kInds_up_single = [];
    while length(use_kInds_up_single) < use_nPix_us*flen
        use_kInds_up_single = find(Xinds_ups(:) >= cur_use_pix(1)-cnt/spatial_usfac & Xinds_ups(:) <= cur_use_pix(end)+cnt/spatial_usfac);
        if length(use_kInds_up_single) < use_nPix_us*flen
        use_kInds_up_single = find(Xinds_ups(:) >= cur_use_pix(1)-(cnt+1)/spatial_usfac & Xinds_ups(:) <= cur_use_pix(end) + cnt/spatial_usfac);
        end
    cnt = cnt + 1;    
    end
end

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
%%
cd(anal_dir);

fprintf('Loading pre-computed initial models\n');
load(mod_data_name);
load(old_anal_name,'dit_mods','tr_set','it_fix_post_mean','drift_post_mean','it_fix_post_std','drift_post_std','it_LL*','dit_LL*');
old_best_mods = dit_mods{end};

best_fix_cor = squeeze(it_fix_post_mean(end,:,:))*spatial_usfac/old_usfac;
best_fix_std = squeeze(it_fix_post_std(end,:,:))*spatial_usfac/old_usfac;
best_drift_cor = squeeze(drift_post_mean(end,:,:))*spatial_usfac/old_usfac;
best_drift_std = squeeze(drift_post_std(end,:,:))*spatial_usfac/old_usfac;
best_LLimp = it_LLimp(end,:);
best_dLLimp = dit_LLimp(end,:);

clear it_fix_post_mean drift_post_mean
clear all_mod_fits

%% SELECT USABLE UNITS AND make Robs_mat
n_tr_chs = length(tr_set);

% make Robs_mat
poss_blocks = 1:(n_blocks-1);
Robs_mat = nan(length(used_inds),n_tr_chs);
for ss = 1:n_tr_chs
    if all_mod_SU(tr_set(ss)) > 0
%         su_probe_ind = find(su_probes == all_mod_SU(tr_set(ss)));
        su_probe_ind = find(SU_numbers == all_mod_SUnum(tr_set(ss)));
        Robs_mat(:,ss) = all_binned_sua(used_inds,su_probe_ind);
    else
        Robs_mat(:,ss) = all_binned_mua(used_inds,tr_set(ss));
    end
end

%don't use separate xv set for eye-tracking
tr_inds = full_inds;

%%
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

%% construct drift-corrected X-mat
fix_post_cor = nan(NT,2);
fix_post_cor(~isnan(fix_ids),:) = best_fix_cor(fix_ids(~isnan(fix_ids)),:);
fix_post_cor = interp1(find(~isnan(fix_ids)),fix_post_cor(~isnan(fix_ids),:),1:NT);
fix_post_cor(isnan(fix_post_cor)) = 0;
%back project drift (within-fixation) by sac_shift
drift_post_cor = best_drift_cor;

for ii = 1:length(trial_start_inds)
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    drift_post_cor(cur_inds(1:end-sac_shift),:) = drift_post_cor(cur_inds(sac_shift+1:end),:);
end
drift_post_cor = interp1(find(~isnan(fix_ids)),drift_post_cor(~isnan(fix_ids),:),1:NT);

drift_post_cor(isnan(drift_post_cor)) = 0;

all_post_cor = round((fix_post_cor+drift_post_cor));

all_shift_stimmat_ver = full_stim_mat_ver_up;
all_shift_stimmat_hor = full_stim_mat_hor_up;
for i=1:NT
    if all_stim_ori(used_inds(i)) == 0
        all_shift_stimmat_hor(used_inds(i),:) = shift_matrix_Nd(all_shift_stimmat_hor(used_inds(i),:),-all_post_cor(i,2),2);
    else
        all_shift_stimmat_ver(used_inds(i),:) = shift_matrix_Nd(all_shift_stimmat_ver(used_inds(i),:),-all_post_cor(i,1),2);
    end
end

shift_hor_Xmat_us = create_time_embedding(all_shift_stimmat_hor,stim_params_us);
shift_ver_Xmat_us = create_time_embedding(all_shift_stimmat_ver,stim_params_us);
all_Xmat_us = [shift_ver_Xmat_us shift_hor_Xmat_us];
clear shift_*Xmat_us

%%
cd(anal_dir);

if ~exist(['./' hrmod_data_name '.mat'],'file') || recompute_init_mods == 1
    for ss = 1:n_tr_chs
        fprintf('Computing base LLs for MU %d of %d\n',ss,n_probes);
        cur_tr_inds = tr_inds(~isnan(Robs_mat(tr_inds,ss)));
        tr_NT = length(cur_tr_inds);
        if ~isempty(cur_tr_inds)
            Robs = Robs_mat(cur_tr_inds,ss);
            null_mod = all_nullmod(tr_set(ss));
            
            tr_X{1} = all_Xmat_us(used_inds(cur_tr_inds),use_kInds_up);
            tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
            
            if use_sac_kerns
                tr_X{3} = Xsac(cur_tr_inds,:);
                tr_X{4} = Xmsac(cur_tr_inds,:);
            end
            
            %spatial up-sampling of filter estimates
            base_filts = reshape([old_best_mods(tr_set(ss)).mods(find(init_Xtargs == 1)).filtK],[flen use_nPix*4 n_squared_filts+1]);
            base_filts_up = zeros(flen,use_nPix_us*2,n_squared_filts+1);
            for ii = 1:use_nPix*4
                base_filts_up(:,2*(ii-1)+1,:) = 0.5*base_filts(:,ii,:);
                base_filts_up(:,2*(ii-1)+2,:) = 0.5*base_filts(:,ii,:);
            end
            base_filts_up = reshape(base_filts_up,use_nPix_us*flen*2,n_squared_filts+1);
            
            init_filts{end} = old_best_mods(tr_set(ss)).mods(find(init_Xtargs==2)).filtK;
            for ii = 1:n_squared_filts+1
                init_filts{ii} = base_filts_up(:,ii);
            end
            gqm2 = NMMinitialize_model(fin_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs,init_filts);
            [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm2,Robs,tr_X);
            gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_custom',base_lambda_custom./var(gint)');
            gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
            if use_sac_kerns
                gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,3,zeros(n_sac_bins,1),sac_reg_params); %sac term
                gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,4,zeros(n_sac_bins,1),sac_reg_params); %nsac term
            end
            all_mod_fits(tr_set(ss)) = gqm2;
            all_mod_fits(tr_set(ss)).spk_NL_params(1) = old_best_mods(tr_set(ss)).spk_NL_params(1);
            all_mod_fits(tr_set(ss)) = NMMfit_filters(all_mod_fits(tr_set(ss)),Robs,tr_X,[],[],silent,[],L2_mat_us);
            
            all_mod_fits_withspkNL(tr_set(ss)) = NMMfit_logexp_spkNL(all_mod_fits(tr_set(ss)),Robs,tr_X);
            
            [LL,~,pred_rate] = NMMmodel_eval(all_mod_fits_withspkNL(tr_set(ss)),Robs,tr_X);
            [null_LL(tr_set(ss)),~,null_prate] = NMMmodel_eval(null_mod,Robs,tr_X(2:end));
            all_mod_LLimp(tr_set(ss)) = (LL-null_LL(tr_set(ss)))/log(2);
        else
            all_mod_LLimp(tr_set(ss)) = nan;
        end
    end
    save(hrmod_data_name,'all_mod*','null_LL');
    
else
    fprintf('Loading pre-computed initial models\n');
    load(hrmod_data_name);
end

%%
clear all_Xmat_us
su_inds = find(all_mod_SU(tr_set) > 0);

%%
%back-project saccade-times
all_fix_post_mean_cor = nan(NT,2);
all_fix_post_mean_cor(~isnan(fix_ids),:) = best_fix_cor(fix_ids(~isnan(fix_ids)),:);
all_fix_post_mean_cor = interp1(find(~isnan(fix_ids)),all_fix_post_mean_cor(~isnan(fix_ids),:),1:NT);
all_fix_post_mean_cor(isnan(all_fix_post_mean_cor)) = 0;

all_fix_post_mean_cor = round(all_fix_post_mean_cor);

all_shift_stimmat_ver = full_stim_mat_ver_up;
all_shift_stimmat_hor = full_stim_mat_hor_up;
for i=1:NT
    if all_stim_ori(used_inds(i)) == 0
        all_shift_stimmat_hor(used_inds(i),:) = shift_matrix_Nd(all_shift_stimmat_hor(used_inds(i),:),-all_fix_post_mean_cor(i,2),2);
    else
        all_shift_stimmat_ver(used_inds(i),:) = shift_matrix_Nd(all_shift_stimmat_ver(used_inds(i),:),-all_fix_post_mean_cor(i,1),2);
    end
end

shift_hor_Xmat_us = create_time_embedding(all_shift_stimmat_hor,stim_params_us);
shift_ver_Xmat_us = create_time_embedding(all_shift_stimmat_ver,stim_params_us);
shift_hor_Xmat_us = shift_hor_Xmat_us(used_inds,:);
shift_ver_Xmat_us = shift_ver_Xmat_us(used_inds,:);
% clear shift_*Xmat_us

%% ITERATE FIXATION-BASED CORRECTIONS

dit_mods{1} = all_mod_fits;
dit_mods_spkNL{1} = all_mod_fits_withspkNL;
dit_LLimp(1,:) = all_mod_LLimp;
dit_R2(1,:) = all_mod_R2;
for nn = 1:n_drift_inf_it
    
    %% PREPROCESS MODEL COMPONENTS
    filt_bank = zeros(n_tr_chs,2*klen_us,n_squared_filts+1);
    lin_kerns = nan(n_tr_chs,n_blocks);
    if use_sac_kerns
        sac_kerns = nan(n_tr_chs,n_sac_bins);
        msac_kerns = nan(n_tr_chs,n_sac_bins);
    end
    mod_spkNL_params = nan(n_tr_chs,3);
    for ss = 1:n_tr_chs
        cur_Xtargs = [dit_mods{nn}(tr_set(ss)).mods(:).Xtarget];
        cur_k = [dit_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 1).filtK];
        n_used_filts = size(cur_k,2);
        filt_bank(ss,:,1:n_used_filts) = cur_k;
        mod_spkNL_params(ss,:) = dit_mods_spkNL{nn}(tr_set(ss)).spk_NL_params;
        lin_kerns(ss,:) = dit_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 2).filtK;
        if use_sac_kerns
            sac_kerns(ss,:) = dit_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 3).filtK;
            msac_kerns(ss,:) = dit_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 4).filtK;
        end
    end
    filt_bank = permute(filt_bank,[2 1 3]);
    
    %indicator predictions
    block_out = Xblock(used_inds,:)*lin_kerns';
    if use_sac_kerns
        sac_out = Xsac*sac_kerns';
        msac_out = Xmsac*msac_kerns';
    end
    
    %% ESTIMATE LL for each shift in each stimulus frame
    
    frame_LLs = nan(NT,n_Dshifts);
    for xx = 1:n_Dshifts
        fprintf('Dshift %d of %d\n',xx,n_Dshifts);
        cur_stim_shift = [shift_ver_Xmat_us*Dshift_mat{DSH_inds(xx,1)} shift_hor_Xmat_us*Dshift_mat{DSH_inds(xx,2)}];
        
        %outputs of stimulus models at current X-matrix shift
        gfuns = ones(NT,n_tr_chs);
        gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,1)');
        gfuns = gfuns + cur_stim_shift*squeeze(filt_bank(:,:,1));
        for ff = 2:(n_squared_filts+1)
            gfuns = gfuns + mod_signs(ff)*(cur_stim_shift*squeeze(filt_bank(:,:,ff))).^2;
        end
        
        %add contributions from extra lin kernels
        gfuns = gfuns + block_out;
        if use_sac_kerns
            gfuns = gfuns + sac_out + msac_out;
        end
        
        %incorporate beta
        gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,2)');
        
        %handle numerical overflow with log(1+exp)
        too_large = gfuns > 50;
        pred_rate = log(1+exp(gfuns));
        pred_rate(too_large) = gfuns(too_large);
        
        %incorporate alpha
        pred_rate = bsxfun(@times,pred_rate,mod_spkNL_params(:,3)');
        
        %enforce min predicted rate
        pred_rate(pred_rate < 1e-50) = 1e-50;
        
        frame_LLs(:,xx) = squeeze(nansum(Robs_mat.*log(pred_rate) - pred_rate,2));
    end
    
    %% INFER DRIFT CORRECTIONS
    edrift_jump_prior = exp(drift_jump_prior)';
    ebase_lA = exp(base_lA);
    
    for ff = 1:n_fixs
        if mod(ff,10)==0
            fprintf('Fixation %d of %d\n',ff,n_fixs);
        end
        
        tset = find(pfix_ids==ff)';
        ntset = length(tset);
        if ntset > drift_dsf
            nt_pts = ceil(ntset/drift_dsf);
            tset_inds = 1+floor((0:(ntset-1))/drift_dsf);
            alpha=zeros(nt_pts,n_Dshifts);
            beta = zeros(nt_pts,n_Dshifts);
            scale = zeros(nt_pts,1);
            
            tpt_loc = ceil(1:drift_dsf:nt_pts*drift_dsf);
            tpt_loc(end) = ntset;
            
%             cur_drift_mean = post_mean_drift(tset,:);
            cur_LL_set = frame_LLs(tset,:);
            if mod(ntset,drift_dsf) ~= 0
                dangling_pts = nt_pts*drift_dsf-ntset;
                cur_LL_set = cat(1,cur_LL_set,zeros(dangling_pts,n_Dshifts));
%                 cur_drift_mean = cat(1,cur_drift_mean,zeros(dangling_pts,2));
            end
            cur_LL_set = reshape(cur_LL_set,[drift_dsf nt_pts n_Dshifts]);
            cur_LL_set = exp(squeeze(sum(cur_LL_set,1)));
            
            %compute rescaled forward messages
            alpha(1,:)= edrift_jump_prior.*cur_LL_set(1,:);
            scale(1)=sum(alpha(1,:));
            alpha(1,:)=alpha(1,:)/scale(1);
            for t=2:nt_pts
                alpha(t,:)=(alpha(t-1,:)*ebase_lA).*cur_LL_set(t,:);
                scale(t)=sum(alpha(t,:));
                alpha(t,:)=alpha(t,:)/scale(t);
            end
            
            %compute rescaled backward messages
            beta(nt_pts,:)=ones(1,n_Dshifts)/scale(nt_pts);
            for t=nt_pts-1:-1:1
                beta(t,:)=(beta(t+1,:).*cur_LL_set(t+1,:))*(ebase_lA')/scale(t);
            end
            
            temp_gamma = alpha.*beta;
%             temp_gamma = bsxfun(@rdivide,temp_gamma,sum(temp_gamma,2));
            int_gamma = interp1(tpt_loc,log(temp_gamma),1:ntset);
            int_gamma = bsxfun(@minus,int_gamma,logsumexp(int_gamma,2));    
            int_gamma = exp(int_gamma);
            
            cur_dpost_mean = nan(ntset,2);
            cur_dpost_std = nan(ntset,2);
            cur_dpost_mean(:,1) = sum(bsxfun(@times,int_gamma,DSH(:,1)'),2);
            cur_dpost_mean(:,2) = sum(bsxfun(@times,int_gamma,DSH(:,2)'),2);
            
            cur_xdiff = bsxfun(@minus,squeeze(cur_dpost_mean(:,1))',DSH(:,1)).^2;
            cur_dpost_std(:,1) = sqrt(sum(cur_xdiff'.*int_gamma,2));
            cur_ydiff = bsxfun(@minus,squeeze(cur_dpost_mean(:,2))',DSH(:,2)).^2;
            cur_dpost_std(:,2) = sqrt(sum(cur_ydiff'.*int_gamma,2));
           
            drift_post_mean(nn,tset,:) = cur_dpost_mean;
            drift_post_std(nn,tset,:) = cur_dpost_std;
            
        end
    end
    uset = find(~isnan(drift_post_mean(nn,:,1)));
    drift_post_mean(nn,:,1) = interp1(uset,drift_post_mean(nn,uset,1),1:NT);
    drift_post_mean(nn,:,2) = interp1(uset,drift_post_mean(nn,uset,2),1:NT);
    drift_post_std(nn,:,1) = interp1(uset,drift_post_std(nn,uset,1),1:NT);
    drift_post_std(nn,:,2) = interp1(uset,drift_post_std(nn,uset,2),1:NT);
   
    
%     lgamma = nan(NT,n_Dshifts);
%     for ff = 1:n_fixs
%         if mod(ff,100)==0
%             fprintf('Fixation %d of %d\n',ff,n_fixs);
%         end
%         
%         tset = find(fix_ids==ff)';
%         tset = find(pfix_ids==ff)';
%         ntset = length(tset);
%         if ntset > drift_dsf
%         nt_pts = ceil(ntset/drift_dsf);
%         tset_inds = 1+floor((0:(ntset-1))/drift_dsf);
%         talpha=zeros(nt_pts,n_Dshifts);
%         tbeta = zeros(nt_pts,n_Dshifts);
%         
%         tpt_loc = ceil(1:drift_dsf:nt_pts*drift_dsf);
%         tpt_loc(end) = ntset;
%         
%         cur_drift_mean = post_mean_drift(tset,:);
%         cur_LL_set = frame_LLs(tset,:);
%         if mod(ntset,drift_dsf) ~= 0
%             dangling_pts = nt_pts*drift_dsf-ntset;
%             cur_LL_set = cat(1,cur_LL_set,zeros(dangling_pts,n_Dshifts));
%             cur_drift_mean = cat(1,cur_drift_mean,zeros(dangling_pts,2));
%         end
%         cur_LL_set = reshape(cur_LL_set,[drift_dsf nt_pts n_Dshifts]);
%         cur_LL_set = squeeze(sum(cur_LL_set,1));
%         
%         cur_drift_mean = drift_dsf*reshape(mean(reshape(cur_drift_mean,[drift_dsf nt_pts 2])),nt_pts,2);
%         
%         if all(use_coils==0)
%             cur_lA = repmat(base_lA,[1 1 nt_pts]);
%         else
%             cur_lA = nan(n_Dshifts,n_Dshifts,nt_pts);
%             for iii = 1:nt_pts
%                 %                 gprobs = diff(erf((Dshift_dx_edges+cur_drift_mean(iii))/(post_drift_sigma*sqrt(2)*drift_dsf)));
%                 %                 cur_lA(:,:,iii) = toeplitz(gprobs(n_Dshifts:end),gprobs(n_Dshifts:-1:1));
%                 %
%                 
%                 cdist = pdist2(DSH*sp_dx,[DSH(:,1)*sp_dx - cur_drift_mean(iii,1) DSH(:,2)*sp_dx - cur_drift_mean(iii,2)]);
%                 cur_lA(:,:,iii) = -cdist.^2/(2*(post_drift_sigma*drift_dsf)^2);
%             end
%             %             cur_lA = log(bsxfun(@rdivide,cur_lA,sum(cur_lA,2)));
%             %             cur_lA(cur_lA < eps) = eps;
%             cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2));
%         end
%         
%         talpha(1,:) = drift_jump_prior' + cur_LL_set(1,:);
%         for t = 2:nt_pts
%             talpha(t,:) = logmulexp(talpha(t-1,:),cur_lA(:,:,t)) + cur_LL_set(t,:);
%         end
%         
%         tbeta(end,:)=log(ones(1,n_Dshifts));
%         for t = (nt_pts-1):-1:1
%             lf1 = tbeta(t+1,:) + cur_LL_set(t+1,:);
%             tbeta(t,:) = logmulexp(lf1,cur_lA(:,:,t)');
%         end
%         temp_gamma = talpha + tbeta;
%         int_gamma = interp1(tpt_loc,temp_gamma,1:ntset);
%         int_gamma = bsxfun(@minus,int_gamma,logsumexp(int_gamma,2));
%         int_gamma = exp(int_gamma);
%         
%         cur_dpost_mean2 = nan(ntset,2);
%         cur_dpost_std2 = nan(ntset,2);
%         cur_dpost_mean2(:,1) = sum(bsxfun(@times,int_gamma,DSH(:,1)'),2);
%         cur_dpost_mean2(:,2) = sum(bsxfun(@times,int_gamma,DSH(:,2)'),2);
%         
%         cur_xdiff = bsxfun(@minus,squeeze(cur_dpost_mean2(:,1))',DSH(:,1)).^2;
%         cur_dpost_std2(:,1) = sqrt(sum(cur_xdiff'.*int_gamma,2));
%         cur_ydiff = bsxfun(@minus,squeeze(cur_dpost_mean2(:,2))',DSH(:,2)).^2;
%         cur_dpost_std2(:,2) = sqrt(sum(cur_ydiff'.*int_gamma,2));
        
        %
        %         if drift_dsf > 1
%             if nt_pts > 1
%                 int_gamma = interp1(tpt_loc,temp_gamma,1:ntset);
%                 lgamma(tset,:) = int_gamma;
%             else
%                 lgamma(tset,:) = repmat(temp_gamma,ntset,1);
%             end
%         else
%             lgamma(tset,:) = temp_gamma;
%         end
%         end
%     end
%     lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
%     gamma = exp(lgamma);
%     
%     drift_post_mean(nn,:,1) = sum(bsxfun(@times,gamma,DSH(:,1)'),2);
%     drift_post_mean(nn,:,2) = sum(bsxfun(@times,gamma,DSH(:,2)'),2);
% %     drift_post_mean_inds(:,1) = sum(bsxfun(@times,gamma,DSH_inds(:,1)'),2);
% %     drift_post_mean_inds(:,2) = sum(bsxfun(@times,gamma,DSH_inds(:,2)'),2);
%     
%     cur_xdiff = bsxfun(@minus,squeeze(drift_post_mean(nn,:,1)),DSH(:,1)).^2';
%     drift_post_std(nn,:,1) = sqrt(sum(cur_xdiff.*gamma,2));
%     cur_ydiff = bsxfun(@minus,squeeze(drift_post_mean(nn,:,2)),DSH(:,2)).^2';
%     drift_post_std(nn,:,2) = sqrt(sum(cur_ydiff.*gamma,2));
%     
%     [~,drift_post_max] = max(gamma,[],2);
%     
%     
% %     drift_post_mean(:,1) = interp1(find(~isnan(fix_ids)),drift_post_mean(~isnan(fix_ids),1),1:NT);
% %     drift_post_mean(:,2) = interp1(find(~isnan(fix_ids)),drift_post_mean(~isnan(fix_ids),2),1:NT);
% %     drift_post_std(:,1) = interp1(find(~isnan(fix_ids)),drift_post_std(~isnan(fix_ids),1),1:NT);
% %     drift_post_std(:,2) = interp1(find(~isnan(fix_ids)),drift_post_std(~isnan(fix_ids),2),1:NT);
%     drift_post_mean(nn,:,1) = interp1(find(~isnan(pfix_ids)),drift_post_mean(nn,~isnan(pfix_ids),1),1:NT);
%     drift_post_mean(nn,:,2) = interp1(find(~isnan(pfix_ids)),drift_post_mean(nn,~isnan(pfix_ids),2),1:NT);
%     drift_post_std(nn,:,1) = interp1(find(~isnan(pfix_ids)),drift_post_std(nn,~isnan(pfix_ids),1),1:NT);
%     drift_post_std(nn,:,2) = interp1(find(~isnan(pfix_ids)),drift_post_std(nn,~isnan(pfix_ids),2),1:NT);
%    
%     %     drift_post_mean_inds(:,1) = round(interp1(find(~isnan(fix_ids)),drift_post_mean_inds(~isnan(fix_ids),1),1:NT));
%     %     drift_post_mean_inds(:,2) = round(interp1(find(~isnan(fix_ids)),drift_post_mean_inds(~isnan(fix_ids),2),1:NT));

    %% construct drift-corrected X-mat
    fix_post_cor = nan(NT,2);
    fix_post_cor(~isnan(fix_ids),:) = best_fix_cor(fix_ids(~isnan(fix_ids)),:);
    fix_post_cor = interp1(find(~isnan(fix_ids)),fix_post_cor(~isnan(fix_ids),:),1:NT);
    fix_post_cor(isnan(fix_post_cor)) = 0;
%     drift_post_cor = squeeze(drift_post_mean(nn,:,:));
    drift_post_cor = squeeze(best_drift_cor);

    for ii = 1:length(trial_start_inds)
        cur_inds = trial_start_inds(ii):trial_end_inds(ii);
        drift_post_cor(cur_inds(1:end-sac_shift),:) = drift_post_cor(cur_inds(sac_shift+1:end),:);
    end
    drift_post_cor = interp1(find(~isnan(fix_ids)),drift_post_cor(~isnan(fix_ids),:),1:NT);
    drift_post_cor(isnan(drift_post_cor)) = 0;
    
    all_dshift_stimmat_hor = all_shift_stimmat_hor;
    all_dshift_stimmat_ver = all_shift_stimmat_ver;
    for i = 1:NT
        if all_stim_ori(used_inds(i)) == 0
            all_dshift_stimmat_hor(used_inds(i),:) = shift_matrix_Nd(all_dshift_stimmat_hor(used_inds(i),:),-round(drift_post_cor(i,2)),2);
        elseif all_stim_ori(used_inds(i)) == 90
            all_dshift_stimmat_ver(used_inds(i),:) = shift_matrix_Nd(all_dshift_stimmat_ver(used_inds(i),:),-round(drift_post_cor(i,1)),2);
        end
    end
    
    dshift_hor_Xmat_us = create_time_embedding(all_dshift_stimmat_hor,stim_params_us);
    dshift_ver_Xmat_us = create_time_embedding(all_dshift_stimmat_ver,stim_params_us);
    dshift_hor_Xmat_us = dshift_hor_Xmat_us(used_inds,:);
    dshift_ver_Xmat_us = dshift_ver_Xmat_us(used_inds,:);
    
%     cur_X{1} = [dshift_ver_Xmat_us(:,use_kInds_up_single) dshift_hor_Xmat_us(:,use_kInds_up_single)];
    cur_X{1} = [dshift_ver_Xmat_us dshift_hor_Xmat_us];
    cur_X{1} = cur_X{1}(:,use_kInds_up);
    clear dshift_*_Xmat_us
    
        cur_X{2} = Xblock(used_inds,:);
    if use_sac_kerns
        cur_X{3} = Xsac;
        cur_X{4} = Xmsac;
    end

    
    %% REFIT ALL CELLS
    silent = 1;
    for ss = 1:length(tr_set)
        cur_cell = tr_set(ss);
        fprintf('Refitting model for tr cell %d of %d\n',ss,length(tr_set));
        cur_unit_ind = find(tr_set == cur_cell);
        cur_tr_uset = tr_inds(~isnan(Robs_mat(tr_inds,cur_unit_ind)));
        
        tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
        
        dit_mods{nn+1}(cur_cell) = dit_mods{nn}(cur_cell);
        dit_mods{nn+1}(cur_cell) = NMMfit_filters(dit_mods{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),...
            tr_X,[],[],silent,[],L2_mat_us); %fit stimulus filters
        
        dit_mods_spkNL{nn+1}(cur_cell) = NMMfit_logexp_spkNL(dit_mods{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
        
        [newLL,~,new_prate] = NMMmodel_eval(dit_mods_spkNL{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
        [~,~,null_prate] = NMMmodel_eval(all_nullmod(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X(2:end));
        [dit_R2(nn+1,cur_cell),dit_dev(nn+1,cur_cell)] = pseudo_r2(Robs_mat(cur_tr_uset,cur_unit_ind),new_prate,null_prate);
        
        dit_LLimp(nn+1,cur_cell) = (newLL - null_LL(cur_cell))/log(2);
        
        fprintf('Original: %.5f  Prev: %.5f  New: %.5f\n',all_mod_LLimp(cur_cell),dit_LLimp(nn,cur_cell),dit_LLimp(nn+1,cur_cell));
    end
    
end

%% SAVE EYE-TRACKING RESULTS
cd(anal_dir);
save(anal_name,'it_*','fix_*','tr_set','dit_*','drift_*')

%%
fin_fix_corr = nan(NT,2);
fin_fix_std = nan(NT,2);
fin_fix_corr(~isnan(fix_ids),:) = best_fix_cor(fix_ids(~isnan(fix_ids)),:);
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids),:),1:NT);
fin_fix_std(~isnan(fix_ids),:) = best_fix_std(fix_ids(~isnan(fix_ids)),:);
fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids),:),1:NT);

fin_fix_corr = fin_fix_corr*sp_dx;
fin_fix_std = fin_fix_std*sp_dx;

% fin_drift_corr = squeeze(drift_post_mean)*sp_dx;
% fin_drift_std = squeeze(drift_post_std)*sp_dx;

% fin_drift_corr = squeeze(best_drift_cor)*sp_dx;
% fin_drift_std = squeeze(best_drift_std)*sp_dx;

fin_drift_corr = squeeze(drift_post_mean(end,:,:))*sp_dx;
fin_drift_std = squeeze(drift_post_std(end,:,:))*sp_dx;

min_fix_dur = 0.15;
fix_inds = [];
long_fix_inds = [];
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
%     if length(cur_inds) > sac_shift
%         fin_drift_corr(cur_inds(1:end-sac_shift+1),:) = fin_drift_corr(cur_inds(sac_shift:end),:);
%         fin_drift_std(cur_inds(1:end-sac_shift+1),:) = fin_drift_std(cur_inds(sac_shift:end),:);
%     end
    fix_inds = [fix_inds cur_inds];
    if length(cur_inds)*dt >= min_fix_dur
        long_fix_inds = [long_fix_inds cur_inds];
    end
end

for ii = 1:length(trial_start_inds)
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    fin_drift_corr(cur_inds(1:end-sac_shift),:) = fin_drift_corr(cur_inds(sac_shift+1:end),:);
    fin_drift_std(cur_inds(1:end-sac_shift),:) = fin_drift_std(cur_inds(sac_shift+1:end),:);
end
fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids),:),1:NT);
fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids),:),1:NT);

fin_tot_corr = fin_fix_corr + fin_drift_corr;
fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2);
% fin_tot_corr = interp1(find(~isnan(fix_ids)),fin_tot_corr(~isnan(fix_ids)),1:NT);
% fin_tot_std = interp1(find(~isnan(fix_ids)),fin_tot_std(~isnan(fix_ids)),1:NT);
    
%%
measured_seq = corrected_eye_vals_interp(used_inds,:);

sac_buff_inds = 2;
min_fix_dur = 0.1;

saccade_prefix = fix_ids(saccade_start_inds);
saccade_postfix = fix_ids(saccade_stop_inds);
saccade_prefix_dur = nan(size(saccade_start_inds));
saccade_postfix_dur = nan(size(saccade_start_inds));
saccade_prefix_dur(~isnan(saccade_prefix)) = fix_durs(saccade_prefix(~isnan(saccade_prefix)));
saccade_postfix_dur(~isnan(saccade_postfix)) = fix_durs(saccade_postfix(~isnan(saccade_postfix)));

too_short = find(saccade_prefix_dur  < min_fix_dur | saccade_postfix_dur < min_fix_dur);
long_enough = setdiff(1:length(saccade_start_inds),too_short);

start_pts = saccade_start_inds(long_enough) - sac_buff_inds;
end_pts = saccade_stop_inds(long_enough) + sac_buff_inds;
start_pts(start_pts < 1) = 1; end_pts(end_pts > length(used_inds)) = length(used_inds);
m_pre_pos = measured_seq(start_pts,:);
m_post_pos = measured_seq(end_pts,:);
m_delta_pos = (m_post_pos - m_pre_pos);


inferred_pre_pos = fin_tot_corr(start_pts,:);
inferred_post_pos = fin_tot_corr(end_pts,:);
inferred_delta_pos = (inferred_post_pos - inferred_pre_pos);

saccade_blocks = all_blockvec(used_inds(saccade_start_inds));

use_micros = ismember(long_enough',micro_sacs) & ~isnan(m_delta_pos(:,1)) & ~isnan(inferred_delta_pos(:,1));
use_nonmicros = ~ismember(long_enough',micro_sacs) & ~isnan(m_delta_pos(:,1)) & ~isnan(inferred_delta_pos(:,1));


[msac_corrs,msac_pvals] = corr(m_delta_pos(use_micros,:),inferred_delta_pos(use_micros,:),'type','spearman');

%%
measured_seq = corrected_eye_vals_interp(used_inds,:);

min_fix_dur = 0.15;
inferred_drift = nan(size(fin_tot_corr));
measured_drift = nan(length(fin_tot_corr),4);
inferred_fix_avg = nan(n_fixs,2);
measured_fix_avg = nan(n_fixs,4);
fix_inds = [];
long_fix_inds = [];
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    fix_inds = [fix_inds cur_inds];
    if length(cur_inds)*dt >= min_fix_dur
        long_fix_inds = [long_fix_inds cur_inds];
    end
    cur_inf = fin_tot_corr(cur_inds,:);
    inferred_fix_avg(ii,:) = median(fin_tot_corr(cur_inds,:));
    inferred_drift(cur_inds,:) = bsxfun(@minus,cur_inf,inferred_fix_avg(ii,:));
    
    measured_fix_avg(ii,:) = median(measured_seq(cur_inds,:));
    measured_drift(cur_inds,:) = bsxfun(@minus,measured_seq(cur_inds,:),measured_fix_avg(ii,:));
end

ufix_inds = fix_inds(~isnan(fin_tot_corr(fix_inds,1)));
[tot_corrs,tot_pvals] = corr(measured_seq(ufix_inds,:),fin_tot_corr(ufix_inds,:),'type','spearman');

ufix_inds = long_fix_inds(~isnan(inferred_drift(long_fix_inds,1)));
[drift_corrs,drift_pvals] = corr(measured_drift(ufix_inds,:),inferred_drift(ufix_inds,:),'type','spearman');

%%
fig_dir = ['/home/james/Analysis/bruce/ET_final/'];

yr = [-0.4 0.4];
print_on = false;
fig_width = 4.86; %3.27 4.86 6.83
rel_height = 1.2;

close all
H = figure();
n_trials = length(unique(all_trialvec));
for tt = 1:n_trials
%     for tt = [96 137 154 179 376 409]
%     for tt = [96 217 249 320 366 400 410 453 485 507 510 610 671]
% for tt = [217 671]
    uu = find(all_trialvec(used_inds) == tt);
    if ~isempty(uu)
        bt = all_t_axis(used_inds(uu(1)));
        et = all_t_axis(used_inds(uu(end)));
        dur = et-bt;
        if dur > 3
            cur_sac_inds = find(ismember(saccade_start_inds,uu));
            rel_sac_start_times = all_t_axis(used_inds(saccade_start_inds(cur_sac_inds))) - bt;
            rel_sac_end_times = all_t_axis(used_inds(saccade_stop_inds(cur_sac_inds))) - bt;
            
            plot_dim = 2;
            subplot(2,1,1)
            hold on
            h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_tot_corr(uu,plot_dim),fin_tot_std(uu,plot_dim),{'color','k'});
%             h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,old_corr(uu,plot_dim),old_std(uu,plot_dim),{'color','m'});
            
            l1=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),plot_dim),'r','linewidth',2);
%             h4=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),plot_dim+2)-median(corrected_eye_vals_interp(used_inds(uu),plot_dim+2)),'color',[0.2 0.8 0.2],'linewidth',2);
%             plot(all_t_axis(used_inds(uu))-bt,smooth_eyepos(used_inds(uu),plot_dim),'r--','linewidth',2);
%             plot(all_t_axis(used_inds(uu))-bt,smooth_eyepos(used_inds(uu),plot_dim+2)-median(smooth_eyepos(used_inds(uu),plot_dim+2)),'--','color',[0.2 0.8 0.2],'linewidth',2);
            
            xlim([0 dur]);
            ylim(yr);
            
            yl = ylim();
            for ii = 1:length(rel_sac_start_times)
                line(rel_sac_start_times([ii ii]),yl,'color','g','linestyle','--');
%                 line(rel_sac_end_times([ii ii]),yl,'color','r','linestyle','--');
            end
            
            xlabel('Time (s)','fontsize',10);
            ylabel('Orthoganol position (deg)','fontsize',10);
            title(sprintf('Trial %d',tt));
            set(gca,'fontsize',8,'fontname','arial');
%             fillPage(gcf,'papersize',[8 5]);
            

            plot_dim = 1;
            subplot(2,1,2)
            hold on
%             h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_fix_corr(uu,plot_dim),fin_fix_std(uu,plot_dim),{'color','m'});
            h2=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_tot_corr(uu,plot_dim),fin_tot_std(uu,plot_dim),{'color','k'});
            
            l2=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),plot_dim),'r','linewidth',2);
%             h4=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),plot_dim+2)-median(corrected_eye_vals_interp(used_inds(uu),plot_dim+2)),'color',[0.2 0.8 0.2],'linewidth',2);
%             plot(all_t_axis(used_inds(uu))-bt,smooth_eyepos(used_inds(uu),plot_dim),'r--','linewidth',2);
%             plot(all_t_axis(used_inds(uu))-bt,smooth_eyepos(used_inds(uu),plot_dim+2)-median(smooth_eyepos(used_inds(uu),plot_dim+2)),'--','color',[0.2 0.8 0.2],'linewidth',2);
            
            xlim([0 dur]);
            ylim(yr);
        
            yl = ylim();
            for ii = 1:length(rel_sac_start_times)
                line(rel_sac_start_times([ii ii]),yl,'color','g','linestyle','--');
%                 line(rel_sac_end_times([ii ii]),yl,'color','r','linestyle','--');
            end
            
            xlabel('Time (s)','fontsize',10);
            ylabel('Orthoganol position (deg)','fontsize',10);
            title(sprintf('Trial %d',tt));
            
            if print_on
                delete(h1.edge);
                delete(h2.edge);
                figufy(H);
                fname = [fig_dir sprintf('Dualori_example_trace_%s_T%d',Expt_name,tt)];
                exportfig(H,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
                close(H);
            else
                pause
                clf(H);
            end
        end
    end
end

%%
close all
f1 = figure();
f2 = figure();
for ss = 1:length(tr_set)
    % sbeg = find(all_mod_SU(tr_set) > 0,1);
    % for ss = sbeg:length(tr_set)
    ss
    init_mod = all_mod_fits(tr_set(ss));
    xtargs = [init_mod.mods(:).Xtarget];
    kmat = [init_mod.mods(xtargs == 1).filtK];
    figure(f1); clf
    subplot(2,2,1)
    imagesc(reshape(kmat(:,1),flen,use_nPix_us));
    ca = max(abs(kmat(:,1))); caxis([-ca ca]);
    for ii = 1:(size(kmat,2)-1)
        subplot(2,2,2+ii)
        imagesc(reshape(kmat(:,ii+1),flen,use_nPix_us));
        ca = max(abs(kmat(:,ii+1))); caxis([-ca ca]);
    end
    colormap(gray)
    
    fin_mod = dit_mods{end}(tr_set(ss));
    xtargs = [fin_mod.mods(:).Xtarget];
    kmat = [fin_mod.mods(xtargs == 1).filtK];
    figure(f2); clf
    subplot(2,2,1)
    imagesc(reshape(kmat(:,1),flen,use_nPix_us));
    ca = max(abs(kmat(:,1))); caxis([-ca ca]);
    for ii = 1:(size(kmat,2)-1)
        subplot(2,2,2+ii)
        imagesc(reshape(kmat(:,ii+1),flen,use_nPix_us));
        ca = max(abs(kmat(:,ii+1))); caxis([-ca ca]);
    end
    colormap(gray)
    
    fprintf('Cell %d of %d\n',ss,length(tr_set));
    fprintf('Original: %.4f  Fin: %.4f\n',all_mod_LLimp(tr_set(ss)),dit_LLimp(end,tr_set(ss)));
    pause
end

