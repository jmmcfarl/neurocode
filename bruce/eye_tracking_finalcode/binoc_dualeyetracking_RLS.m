clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

Expt_num = 86;
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

bar_ori = 90;
use_eye = 'joint';

if bar_ori == 0
    mod_data_name = 'eyecorr_hbar_mods';
    joint_anal_name = 'eyecorr_hbar';
% base_eyecorr_name = 'monoc_eyecorr_hbar';
base_eyecorr_name = 'joint_eyecorr_hbar';
else
    mod_data_name = 'eyecorr_vbar_mods';
    joint_anal_name = 'eyecorr_vbar';
%     base_eyecorr_name = 'monoc_eyecorr_vbar';
    base_eyecorr_name = 'joint_eyecorr_vbar';
end
mod_data_name = strcat(use_eye,'_',mod_data_name);
joint_anal_name = strcat(use_eye,'_',joint_anal_name);
joint_anal_name = [joint_anal_name '_v2'];

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

fix_prior_sigma = 0.15;
fix_noise_sigma = 0.1;
drift_noise_sigma = 0.003;
drift_prior_sigma = 0.005;
drift_jump_sigma = 0.075; %0.05 start
drift_dsf = 3;

min_trial_dur = 0.75;

spatial_usfac = 2;

%%
eps = -1e3;

stim_fs = 100; %in Hz
dt = 0.01;
full_nPix = 36;
stim_params = NIMcreate_stim_params([flen full_nPix*2],dt);
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
dshift = 1;
shifts = -max_shift:dshift:max_shift;
n_shifts = length(shifts);
zero_frame = find(shifts==0);
temp_dx = [-2*max_shift:dshift:2*max_shift];
shift_dx_edges = [temp_dx*sp_dx-dshift*sp_dx/2 temp_dx(end)*sp_dx+dshift*sp_dx/2];
ovshift_dx_edges = [shifts*sp_dx-dshift*sp_dx/2 shifts(end)*sp_dx+dshift*sp_dx/2];

max_Dshift = round(8*spatial_usfac);
Dshifts = -max_Dshift:dshift:max_Dshift;
n_Dshifts = length(Dshifts);
zero_Dframe = find(Dshifts==0);
temp_dx = [-2*max_Dshift:dshift:2*max_Dshift];
Dshift_dx_edges = [temp_dx*sp_dx-dshift*sp_dx/2 temp_dx(end)*sp_dx+dshift*sp_dx/2];
ovDshift_dx_edges = [Dshifts*sp_dx-dshift*sp_dx/2 Dshifts(end)*sp_dx+dshift*sp_dx/2];

max_Tshift = max_shift + max_Dshift;
Tshifts = -max_Tshift:dshift:max_Tshift;

%% load overall su data
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

%%
if strcmp(Expt_name,'G093')
    include_expts = {'rls.FaXwi','rls.FaXwiXimi'};
else
    include_expts = {'rls.Fa', 'rls.FaXimi'};
end
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
% cur_block_set = find(included_type & expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori);
cur_block_set = find(included_type & expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori & expt_dds == 12);

cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];

n_blocks = length(cur_block_set);

sim_sac_expts = find(~expt_has_ds(cur_block_set));
imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');

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
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_wi = [];
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
            cur_stim_mat = cat(2,double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix)),...
                double(right_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix)));
            if ~isempty(all_stim_times)
                if any(cur_stim_times < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times];
            all_t_axis = [all_t_axis; cur_t_axis];
            all_stim_mat = [all_stim_mat; cur_stim_mat];
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
    all_stimmat_up = zeros(size(all_stim_mat,1),full_nPix_us*2);
    for ii = 1:size(all_stim_mat,2)
        all_stimmat_up(:,2*(ii-1)+1) = all_stim_mat(:,ii);
        all_stimmat_up(:,2*(ii-1)+2) = all_stim_mat(:,ii);
    end
elseif spatial_usfac == 1
    all_stimmat_up = all_stim_mat;
else
    error('Unsupported spatial Up-sample factor!');
end
stim_params_us = NMMcreate_stim_params([flen full_nPix_us*2],dt);
all_Xmat_us = create_time_embedding(all_stimmat_up,stim_params_us);

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
if strcmp(Expt_name,'G093')
    un_wi_vals = unique(all_trial_wi);
    use_wi_trials = find(all_trial_wi == un_wi_vals(2));
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
trial_toffset = zeros(length(cur_block_set),1);
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

%compute corrected eye data in bar-oriented frame
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,expt_bar_ori(cur_block_set),used_inds);

[saccades,et_params] = detect_saccades(corrected_eye_vals,all_eye_speed,all_eye_ts,et_params);

par_thresh = 4;
orth_thresh = 1;
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
for i = 1:length(saccade_start_inds)
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
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
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

%%
cd(anal_dir)
load(mod_data_name,'all_mod_SU*');
load(joint_anal_name,'dit_mods','et_tr_set');
tr_set = et_tr_set;
joint_mods = dit_mods{end};
load(base_eyecorr_name,'it_fix*','drift_post*');

fin_fix_corr = nan(NT,1);
fin_fix_std = nan(NT,1);
fin_fix_corr(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
fin_fix_std(~isnan(fix_ids)) = it_fix_post_std(end,fix_ids(~isnan(fix_ids)));
fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);
fin_drift_corr = drift_post_mean(end,:);
fin_drift_std = drift_post_std(end,:);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > sac_shift
        fin_drift_corr(cur_inds(1:end-sac_shift+1)) = fin_drift_corr(cur_inds(sac_shift:end));
        fin_drift_std(cur_inds(1:end-sac_shift+1)) = fin_drift_std(cur_inds(sac_shift:end));
    end
end

fin_tot_corr = round(fin_fix_corr + fin_drift_corr);

%%
left_eye_Sinds = 1:(full_nPix_us);
right_eye_Sinds = (full_nPix_us+1):(2*full_nPix_us);

all_shift_stimmat_left = all_stimmat_up(:,left_eye_Sinds);
all_shift_stimmat_right = all_stimmat_up(:,right_eye_Sinds);
for i=1:NT
    all_shift_stimmat_left(used_inds(i),:) = shift_matrix_Nd(all_shift_stimmat_left(used_inds(i),:),-fin_tot_corr(i),2);
    all_shift_stimmat_right(used_inds(i),:) = shift_matrix_Nd(all_shift_stimmat_right(used_inds(i),:),-fin_tot_corr(i),2);
end
shift_Xmat = create_time_embedding([all_shift_stimmat_left all_shift_stimmat_right],stim_params_us);

%%
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

%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITS

fin_stim_params(1) = NMMcreate_stim_params([flen use_nPix_us*2],dt);
fin_stim_params(2) = NMMcreate_stim_params([n_blocks 1],dt);
if use_sac_kerns
    fin_stim_params(3) = NMMcreate_stim_params(n_sac_bins,dt);
    fin_stim_params(4) = NMMcreate_stim_params(n_sac_bins,dt);
end

null_stim_params = fin_stim_params(2:end);

block_L2 = 1;
silent = 1;
sac_d2t = 100;

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


base_lambda_custom = 20;
base_lambda_L1 = 4;

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
init_d2XT = [ones(n_squared_filts+1,1); 0;];
init_L2 = [zeros(n_squared_filts+1,1); block_L2];
init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_d2XT,'lambda_L2',init_L2);
init_Xtargs = [ones(n_squared_filts+1,1); 2];

init_filts = cell(length(mod_signs),1);
cd(anal_dir);

for ss = 1:n_tr_chs;
    fprintf('Computing base LLs for MU %d of %d\n',ss,n_tr_chs);
    cur_tr_inds = tr_inds(~isnan(Robs_mat(tr_inds,ss)));
    tr_NT = length(cur_tr_inds);
    Robs = Robs_mat(cur_tr_inds,ss);
    
    tr_X{1} = shift_Xmat(used_inds(cur_tr_inds),use_kInds_up);
    tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
    if use_sac_kerns
        tr_X{3} = Xsac(cur_tr_inds,:);
        tr_X{4} = Xmsac(cur_tr_inds,:);
        null_mod = NMMinitialize_model(null_stim_params,[1 1 1],repmat({'lin'},1,3),null_reg_params,[1 2 3]);
    else
        null_mod = NMMinitialize_model(null_stim_params,[1],{'lin'},null_reg_params,[1]);
    end
    
    null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));
    
    all_mod_fits(ss) = joint_mods(tr_set(ss));
    all_mod_fits(ss).stim_params = fin_stim_params;
    all_mod_fits(ss).mods(n_squared_filts+2).filtK = zeros(n_blocks,1);
    all_mod_fits(ss) = NMMfit_filters(all_mod_fits(ss),Robs,tr_X,[],(n_squared_filts+2):length(all_mod_fits(ss).mods),silent,[],L2_mat);
    all_mod_fits_withspkNL(ss) = NMMfit_logexp_spkNL(all_mod_fits(ss),Robs,tr_X);
    
    null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));
    all_nullmod(ss) = null_mod;
       
    [LL,~,pred_rate] = NMMmodel_eval(all_mod_fits_withspkNL(ss),Robs,tr_X);
    [null_LL(ss),~,null_prate] = NMMmodel_eval(null_mod,Robs,tr_X(2:end));
    [all_mod_R2(ss),all_mod_dev(ss),all_null_dev(ss)] = pseudo_r2(Robs,pred_rate,null_prate);
    all_mod_LLimp(ss) = (LL-null_LL(ss))/log(2);
end
clear shift_Xmat

%% DEFINE POINTS TO ALLOW RAPID CHANGES (TRIAL STARTS AFTER SACCADES)
%push the effects of saccades forward in time
sac_shift = round(0.05/dt);
pfix_start_inds = fix_start_inds;
pfix_stop_inds = fix_stop_inds;
for i = 1:length(saccade_start_inds)
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
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
end

trial_ids = nan(NT,1);
for ii = 1:n_trials
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    trial_ids(cur_inds) = ii;
end

%%
left_eye_inds = 1:(flen*full_nPix_us);
right_eye_inds = (flen*full_nPix_us+1):(2*flen*full_nPix_us);

[Xinds_up,Tinds] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
use_kInds_up_single = find(Xinds_up(:) >= cur_use_pix(1)-1/spatial_usfac & Xinds_up(:) <= cur_use_pix(end));

x_shifts = -max_Dshift:dshift:max_Dshift;
y_shifts = -max_Dshift:dshift:max_Dshift;
% y_shifts = 0;
n_xshifts = length(x_shifts);

[Xsh,Ysh] = meshgrid(x_shifts,y_shifts);
[Xsh_inds,Ysh_inds] = meshgrid(1:length(x_shifts),1:length(y_shifts));
SH = [Xsh(:) Ysh(:)];
SH_inds = [Xsh_inds(:) Ysh_inds(:)];
n_shifts = size(SH,1);
zero_frame = find(SH(:,1) == 0 & SH(:,2) == 0);

%generate shift matrices. Must be applied to the stimulus (not the filters)
It = speye(flen);
xshift_mat = cell(n_xshifts,1);
for xx = 1:n_xshifts
    temp = spdiags( ones(full_nPix_us,1), -x_shifts(xx), full_nPix_us, full_nPix_us);
    temp = kron(temp,It);
%     xshift_mat{xx} = temp(:,use_kInds_up(1:flen*use_nPix_us));
    xshift_mat{xx} = temp;
end

shift_mat = cell(n_shifts,1);
for xx = 1:n_shifts
       tempF = spalloc(2*full_nPix_us*flen,2*full_nPix_us*flen,2*full_nPix_us*flen); 
       tempF(left_eye_inds,left_eye_inds) = xshift_mat{SH_inds(xx,1)};
       tempF(right_eye_inds,right_eye_inds) = xshift_mat{SH_inds(xx,2)};
%        tempF(right_eye_inds,right_eye_inds) = xshift_mat{SH_inds(xx,1)};
       shift_mat{xx} = tempF(:,use_kInds_up);
end


% x_Dshifts = -max_Dshift:dshift:max_Dshift;
% y_Dshifts = -max_Dshift:dshift:max_Dshift;
% n_xDshifts = length(x_Dshifts);
% [XDsh,YDsh] = meshgrid(x_Dshifts,y_Dshifts);
% [XDsh_inds,YDsh_inds] = meshgrid(1:length(x_Dshifts),1:length(y_Dshifts));
% DSH = [XDsh(:) YDsh(:)];
% DSH_inds = [XDsh_inds(:) YDsh_inds(:)];
% Dshift_mat = cell(n_xDshifts,1);
% for xx = 1:n_Dshifts
%     temp = spdiags( ones(full_nPix_us,1), -x_Dshifts(xx), full_nPix_us, full_nPix_us);
%     temp = kron(temp,It);
%     Dshift_mat{xx} = temp(:,use_kInds_up(1:flen*use_nPix_us));
% end

%%
warning('OFF','stats:statrobustfit:IterationLimit');

measured_eyepos = [corrected_eye_vals_interp(:,2) corrected_eye_vals_interp(:,4)];
max_sim_pos = max_shift*sp_dx;
measured_eyepos = tanh(measured_eyepos/max_sim_pos)*max_sim_pos;

measured_eyepos(isnan(measured_eyepos)) = 0;

%smooth out fast transients in eye signal
eye_smooth_sig = round(0.025/dt);
interp_inds = [];
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > eye_smooth_sig*5;
        measured_eyepos(used_inds(cur_inds),1) = jmm_smooth_1d_cor(measured_eyepos(used_inds(cur_inds),1),eye_smooth_sig,2);
        measured_eyepos(used_inds(cur_inds),2) = jmm_smooth_1d_cor(measured_eyepos(used_inds(cur_inds),2),eye_smooth_sig,2);
    end
    interp_inds = [interp_inds; cur_inds'];
end
interp_inds = unique(interp_inds);
measured_eyepos(used_inds,:) = interp1(used_inds(interp_inds),measured_eyepos(used_inds(interp_inds),:),used_inds);
measured_eyepos(isnan(measured_eyepos)) = 0;

measured_eyepos = measured_eyepos(used_inds,:);

measured_drift = nan(length(used_inds),2);
measured_fix_avg = nan(n_fixs,2);
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
        back_measured_drift(cur_inds(sac_shift+1:end),:) = measured_drift(cur_inds(1:(end-sac_shift:end)),:);
    end
end
measured_drift = back_measured_drift;
clear back_measured_drift

measured_drift = [zeros(1,2); diff(measured_drift)];

%if using both coils
if use_coils(1) == 1 && use_coils(2) == 1
    
    measured_fix_deltas = nan(n_fixs,1);
    measured_fix_deltas(2:end) = mean(diff(measured_fix_avg),2);
    fix_noise_sigma = fix_noise_sigma/sqrt(2); %noise variance decreases by factor of sqrt(2) with 2 independent measures
    post_drift_var = 1/(2/drift_noise_sigma^2 + 1/drift_prior_sigma^2);
    post_mean_drift = post_drift_var*2*mean(measured_drift,2)/drift_noise_sigma^2;
    
    %if using left coil only
elseif use_coils(1) == 1
    measured_fix_deltas = nan(n_fixs,1);
    measured_fix_deltas(2:end) = diff(measured_fix_avg(:,1));
    post_drift_var = 1/(1/drift_noise_sigma^2 + 1/drift_prior_sigma^2);
    post_mean_drift = post_drift_var*measured_drift(:,1)/drift_noise_sigma^2;
    
    %if using right coil only
elseif use_coils(2) == 1
    measured_fix_deltas = nan(n_fixs,1);
    measured_fix_deltas(2:end) = diff(measured_fix_avg(:,2));
    post_drift_var = 1/(1/drift_noise_sigma^2 + 1/drift_prior_sigma^2);
    post_mean_drift = post_drift_var*measured_drift(:,2)/drift_noise_sigma^2;
    
else
    post_drift_var = drift_prior_sigma^2;
    post_mean_drift = zeros(NT,1);
end

post_drift_sigma = sqrt(post_drift_var);
post_mean_drift(isnan(post_mean_drift)) = 0;

clear drift_post_*

%%
fix_prior_sigma = 0.2;
drift_jump_sigma = 0.2;
drift_prior_sigma = 0.0045;

fix_prior = -sum((SH*sp_dx).^2,2)/(2*fix_prior_sigma^2);
fix_prior = bsxfun(@minus,fix_prior,logsumexp(fix_prior)); %normalize

drift_jump_prior = -sum((SH*sp_dx).^2,2)/(2*(drift_jump_sigma)^2);
drift_jump_prior = bsxfun(@minus,drift_jump_prior,logsumexp(drift_jump_prior)); %normalize

base_lA = pdist2(SH*sp_dx,SH*sp_dx);
base_lA = -base_lA.^2/(2*(drift_prior_sigma*drift_dsf)^2);
base_lA = bsxfun(@minus,base_lA,logsumexp(base_lA)); %normalize


%%fix_prior = -sum((SH*sp_dx).^2,2)/(2*fix_prior_sigma^2);
fix_prior = bsxfun(@minus,fix_prior,logsumexp(fix_prior)); %normalize

su_inds = find(all_mod_SU(tr_set) > 0);

Sused_blocks = 1:length(cur_block_set);
subind_set = find(ismember(all_blockvec(used_inds),Sused_blocks));
used_subinds = used_inds(subind_set);
sNT = length(used_subinds);
used_fix_set = find(ismember(fix_start_inds,subind_set));
%% ITERATE FIXATION-BASED CORRECTIONS
it_mods{1} = all_mod_fits;
it_mods_spkNL{1} = all_mod_fits_withspkNL;
it_LLimp(1,:) = all_mod_LLimp;
it_R2(1,:) = all_mod_R2;
it_dev(1,:) = all_mod_dev;
for nn = 1:n_fix_inf_it
    fprintf('Inferring fixation corrections, iter %d of %d\n',nn,n_fix_inf_it);
    
    %% PREPROCESS MODEL COMPONENTS
    filt_bank = zeros(n_tr_chs,2*klen_us,n_squared_filts+1);
    lin_kerns = nan(n_tr_chs,n_blocks);
    if use_sac_kerns
        sac_kerns = nan(n_tr_chs,n_sac_bins);
        msac_kerns = nan(n_tr_chs,n_sac_bins);
    end
    mod_spkNL_params = nan(n_tr_chs,3);
    for ss = 1:n_tr_chs
        cur_Xtargs = [it_mods{nn}(ss).mods(:).Xtarget];
        cur_k = [it_mods{nn}(ss).mods(cur_Xtargs == 1).filtK];
        n_used_filts = size(cur_k,2);
        filt_bank(ss,:,1:n_used_filts) = cur_k;
        mod_spkNL_params(ss,:) = it_mods_spkNL{nn}(ss).spk_NL_params;
        lin_kerns(ss,:) = it_mods{nn}(ss).mods(cur_Xtargs == 2).filtK;
        if use_sac_kerns
            sac_kerns(ss,:) = it_mods{nn}(ss).mods(cur_Xtargs == 3).filtK;
            msac_kerns(ss,:) = it_mods{nn}(ss).mods(cur_Xtargs == 4).filtK;
        end
    end
    filt_bank = permute(filt_bank,[2 1 3]);
    
    %indicator predictions
    block_out = Xblock(used_inds(subind_set),:)*lin_kerns';
    if use_sac_kerns
        sac_out = Xsac(subind_set,:)*sac_kerns';
        msac_out = Xmsac(subind_set,:)*msac_kerns';
    end
    
    %% ESTIMATE LL for each shift in each stimulus frame
%     cur_Xmat = shift_Xmat(used_subinds,:);
    cur_Xmat = all_Xmat_us(used_subinds,:);
    %precompute LL at all shifts for all units
    frame_LLs = nan(sNT,n_shifts);
    for xx = 1:n_shifts
        fprintf('Shift %d of %d\n',xx,n_shifts);
        cur_stim_shift = cur_Xmat*shift_mat{xx};
        
        %outputs of stimulus models at current X-matrix shift
        gfuns = ones(sNT,n_tr_chs);
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
        
        frame_LLs(:,xx) = squeeze(nansum(Robs_mat(subind_set,:).*log(pred_rate) - pred_rate,2));
    end
    
        %% INFER DRIFT CORRECTIONS
        
        n_Sfixs = length(used_fix_set);
        lgamma = nan(sNT,n_shifts);
        for ff = 1:n_Sfixs
%         if mod(ff,100)==0
            fprintf('Fixation %d of %d\n',ff,n_Sfixs);
%         end
        
        tset = find(fix_ids(subind_set)==used_fix_set(ff))';
        ntset = length(tset);
        nt_pts = ceil(ntset/drift_dsf);
        tset_inds = 1+floor((0:(ntset-1))/drift_dsf);
        talpha=zeros(nt_pts,n_shifts);
        tbeta = zeros(nt_pts,n_shifts);
        
        tpt_loc = ceil(1:drift_dsf:nt_pts*drift_dsf);
        tpt_loc(end) = ntset;
        
        cur_drift_mean = post_mean_drift(tset);
        cur_LL_set = frame_LLs(tset,:);
        if mod(ntset,drift_dsf) ~= 0
            dangling_pts = nt_pts*drift_dsf-ntset;
            cur_LL_set = cat(1,cur_LL_set,zeros(dangling_pts,n_shifts));
            cur_drift_mean = cat(1,cur_drift_mean,zeros(dangling_pts,1));
        end
        cur_LL_set = reshape(cur_LL_set,[drift_dsf nt_pts n_shifts]);
        cur_LL_set = squeeze(sum(cur_LL_set,1));
        
        cur_drift_mean = drift_dsf*mean(reshape(cur_drift_mean,[drift_dsf nt_pts]));
        
        if all(use_coils==0)
            cur_lA = repmat(base_lA,[1 1 nt_pts]);
        else
            cur_lA = nan(n_shifts,n_shifts,nt_pts);
            for iii = 1:nt_pts
                cdist = pdist2(shifts'*sp_dx + cur_drift_mean(iii),shifts'*sp_dx);
                cur_lA(:,:,iii) = -cdist.^2/(2*(post_drift_sigma*drift_dsf)^2);
            end
            cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2));
        end
        
        talpha(1,:) = drift_jump_prior' + cur_LL_set(1,:);
        for t = 2:nt_pts
            talpha(t,:) = logmulexp(talpha(t-1,:),cur_lA(:,:,t)) + cur_LL_set(t,:);
        end
        
        tbeta(end,:)=log(ones(1,n_shifts));
        for t = (nt_pts-1):-1:1
            lf1 = tbeta(t+1,:) + cur_LL_set(t+1,:);
            tbeta(t,:) = logmulexp(lf1,cur_lA(:,:,t)');
        end
        temp_gamma = talpha + tbeta;
        
        if drift_dsf > 1
            if nt_pts > 1
                int_gamma = interp1(tpt_loc,temp_gamma,1:ntset);
                lgamma(tset,:) = int_gamma;
            else
                lgamma(tset,:) = repmat(temp_gamma,ntset,1);
            end
        else
            lgamma(tset,:) = temp_gamma;
        end
    end
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    gamma = exp(lgamma);
    
        drift_post_mean(:,1) = sum(bsxfun(@times,gamma,SH(:,1)'),2);
    drift_post_mean(:,2) = sum(bsxfun(@times,gamma,SH(:,2)'),2);

        cur_xdiff = bsxfun(@minus,drift_post_mean(:,1)',SH(:,1)).^2';
    drift_post_std(:,1) = sqrt(sum(cur_xdiff.*gamma,2));
    cur_ydiff = bsxfun(@minus,drift_post_mean(:,2)',SH(:,2)).^2';
    drift_post_std(:,2) = sqrt(sum(cur_ydiff.*gamma,2));

    drift_post_mean(:,1) = interp1(find(~isnan(fix_ids(subind_set))),drift_post_mean(~isnan(fix_ids(subind_set)),1),1:sNT);
    drift_post_mean(:,2) = interp1(find(~isnan(fix_ids(subind_set))),drift_post_mean(~isnan(fix_ids(subind_set)),2),1:sNT);
    drift_post_std(:,1) = interp1(find(~isnan(fix_ids(subind_set))),drift_post_std(~isnan(fix_ids(subind_set)),1),1:sNT);
    drift_post_std(:,2) = interp1(find(~isnan(fix_ids(subind_set))),drift_post_std(~isnan(fix_ids(subind_set)),2),1:sNT);
    
    %%
    
    all_shift_stimmat_left = all_stimmat_up(:,left_eye_Sinds);
    all_shift_stimmat_right = all_stimmat_up(:,right_eye_Sinds);
    for i = 1:NT
        all_shift_stimmat_left(used_inds(i),:) = shift_matrix_Nd(all_shift_stimmat_left(used_inds(i),:),-round(drift_post_mean(i,1)),2);
        all_shift_stimmat_right(used_inds(i),:) = shift_matrix_Nd(all_shift_stimmat_right(used_inds(i),:),-round(drift_post_mean(i,2)),2);
    end
 shift_Xmat = create_time_embedding([all_shift_stimmat_left all_shift_stimmat_right],stim_params_us);
   shift_Xmat = shift_Xmat(used_inds,:);
%     all_fix_post_mean_cor = round(all_fix_post_mean_cor) + max_Tshift + 1;
%     all_shift_stimmat_up = all_stimmat_up;
%     for i=1:NT
%         all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_fix_post_mean_cor(i)};
%     end
%     all_Xmat_up_fixcor = create_time_embedding(all_shift_stimmat_up,stim_params_us);
%     all_Xmat_up_fixcor = all_Xmat_up_fixcor(used_inds,:);
    
    %% REFIT ALL CELLS
    cur_X{1} = shift_Xmat(:,use_kInds_up);
    cur_X{2} = Xblock(used_inds,:);
    if use_sac_kerns
        cur_X{3} = Xsac;
        cur_X{4} = Xmsac;
    end
    
    silent = 1;
    for ss = 1:length(tr_set)
        cur_cell = tr_set(ss);
        fprintf('Refitting model for tr cell %d of %d\n',ss,length(tr_set));
        cur_unit_ind = find(tr_set == cur_cell);
        cur_tr_uset = tr_inds(~isnan(Robs_mat(tr_inds,cur_unit_ind)));
        
        tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
        
        it_mods{nn+1}(cur_unit_ind) = it_mods{nn}(cur_unit_ind);
        it_mods{nn+1}(cur_unit_ind) = NMMfit_filters(it_mods{nn+1}(cur_unit_ind),Robs_mat(cur_tr_uset,cur_unit_ind),...
            tr_X,[],[],silent); %fit stimulus filters
        
        %refit spk NL
        it_mods_spkNL{nn+1}(cur_unit_ind) = NMMfit_logexp_spkNL(it_mods{nn+1}(cur_unit_ind),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
        
        [newLL,~,new_prate] = NMMmodel_eval(it_mods_spkNL{nn+1}(cur_unit_ind),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
        [~,~,null_prate] = NMMmodel_eval(all_nullmod(cur_unit_ind),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X(2:end));
        [it_R2(nn+1,cur_unit_ind),it_dev(nn+1,cur_unit_ind)] = pseudo_r2(Robs_mat(cur_tr_uset,cur_unit_ind),new_prate,null_prate);
        
        it_LLimp(nn+1,cur_unit_ind) = (newLL - null_LL(cur_unit_ind))/log(2);
        if nn > 1
            fprintf('Original: %.4f  Prev: %.4f  New: %.4f\n',all_mod_LLimp(cur_unit_ind),it_LLimp(nn,cur_unit_ind),it_LLimp(nn+1,cur_unit_ind));
        else
            fprintf('Original: %.4f  New: %.4f\n',all_mod_LLimp(cur_unit_ind),it_LLimp(nn+1,cur_unit_ind));
        end
    end
    
end

%% NOW INFER DRIFT CORRECTIONS

dit_mods{1} = it_mods{n_fix_inf_it+1};
dit_mods_spkNL{1} = it_mods_spkNL{n_fix_inf_it+1};
dit_LLimp(1,:) = it_LLimp(n_fix_inf_it+1,:);
dit_R2(1,:) = it_R2(n_fix_inf_it+1,:);
dit_dev(1,:) = it_dev(n_fix_inf_it+1,:);
if use_LOOXV
    for xv = 1:length(su_inds)
        dit_mods_LOO{xv,1} = it_mods_LOO{xv,n_fix_inf_it+1};
        dit_mods_spkNL_LOO{xv,1} = it_mods_spkNL_LOO{xv,n_fix_inf_it+1};
    end
    dit_LLimp_LOO(:,1,:) = it_LLimp_LOO(:,n_fix_inf_it+1,:);
    dit_R2_LOO(:,1,:) = it_R2_LOO(:,n_fix_inf_it+1,:);
    dit_dev_LOO(:,1,:) = it_dev_LOO(:,n_fix_inf_it+1,:);
end
for nn = 1:n_drift_inf_it
    fprintf('Inferring drift corrections, iter %d of %d\n',nn,n_drift_inf_it);
    if use_LOOXV
        %re-create all-cell fixation-corrected Xmat
        all_fix_post_mean_cor(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
        all_fix_post_mean_cor = interp1(find(~isnan(fix_ids)),all_fix_post_mean_cor(~isnan(fix_ids)),1:NT);
        all_fix_post_mean_cor(isnan(all_fix_post_mean_cor)) = 0;
        all_fix_post_mean_cor = round(all_fix_post_mean_cor) + max_Tshift + 1;
        for i=1:NT
            all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_fix_post_mean_cor(i)};
        end
        all_Xmat_up_fixcor = create_time_embedding(all_shift_stimmat_up,stim_params_us);
        all_Xmat_up_fixcor = all_Xmat_up_fixcor(used_inds,:);
    end
    
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
    for xx = 1:length(Dshifts)
        fprintf('Dshift %d of %d\n',xx,n_Dshifts);
        cur_stim_shift = all_Xmat_up_fixcor*Dshift_mat{xx};
        
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
    
    lgamma = nan(NT,n_Dshifts);
    for ff = 1:n_fixs
        if mod(ff,100)==0
            fprintf('Fixation %d of %d\n',ff,n_fixs);
        end
        
        tset = find(fix_ids==ff)';
        ntset = length(tset);
        nt_pts = ceil(ntset/drift_dsf);
        tset_inds = 1+floor((0:(ntset-1))/drift_dsf);
        talpha=zeros(nt_pts,n_Dshifts);
        tbeta = zeros(nt_pts,n_Dshifts);
        
        tpt_loc = ceil(1:drift_dsf:nt_pts*drift_dsf);
        tpt_loc(end) = ntset;
        
        cur_drift_mean = post_mean_drift(tset);
        cur_LL_set = frame_LLs(tset,:);
        if mod(ntset,drift_dsf) ~= 0
            dangling_pts = nt_pts*drift_dsf-ntset;
            cur_LL_set = cat(1,cur_LL_set,zeros(dangling_pts,n_Dshifts));
            cur_drift_mean = cat(1,cur_drift_mean,zeros(dangling_pts,1));
        end
        cur_LL_set = reshape(cur_LL_set,[drift_dsf nt_pts n_Dshifts]);
        cur_LL_set = squeeze(sum(cur_LL_set,1));
        
        cur_drift_mean = drift_dsf*mean(reshape(cur_drift_mean,[drift_dsf nt_pts]));
        
        if all(use_coils==0)
            cur_lA = repmat(base_lA,[1 1 nt_pts]);
        else
            cur_lA = nan(n_Dshifts,n_Dshifts,nt_pts);
            for iii = 1:nt_pts
                %                                 gprobs = diff(erf((Dshift_dx_edges+cur_drift_mean(iii))/(post_drift_sigma*sqrt(2)*drift_dsf)));
                %                                 cur_lA(:,:,iii) = toeplitz(gprobs(n_Dshifts:end),gprobs(n_Dshifts:-1:1));
                cdist = pdist2(Dshifts'*sp_dx + cur_drift_mean(iii),Dshifts'*sp_dx);
                cur_lA(:,:,iii) = -cdist.^2/(2*(post_drift_sigma*drift_dsf)^2);
            end
            %             cur_lA = log(bsxfun(@rdivide,cur_lA,sum(cur_lA,2)));
            %             cur_lA(cur_lA < eps) = eps;
            cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2));
        end
        
        talpha(1,:) = drift_jump_prior + cur_LL_set(1,:);
        for t = 2:nt_pts
            talpha(t,:) = logmulexp(talpha(t-1,:),cur_lA(:,:,t)) + cur_LL_set(t,:);
        end
        
        tbeta(end,:)=log(ones(1,n_Dshifts));
        for t = (nt_pts-1):-1:1
            lf1 = tbeta(t+1,:) + cur_LL_set(t+1,:);
            tbeta(t,:) = logmulexp(lf1,cur_lA(:,:,t)');
        end
        temp_gamma = talpha + tbeta;
        
        if drift_dsf > 1
            if nt_pts > 1
                int_gamma = interp1(tpt_loc,temp_gamma,1:ntset);
                lgamma(tset,:) = int_gamma;
            else
                lgamma(tset,:) = repmat(temp_gamma,ntset,1);
            end
        else
            lgamma(tset,:) = temp_gamma;
        end
    end
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    gamma = exp(lgamma);
    drift_post_mean(nn,:) = sum(bsxfun(@times,gamma,Dshifts),2);
    drift_post_std(nn,:) = sqrt(sum(bsxfun(@times,gamma,Dshifts.^2),2) - drift_post_mean(nn,:)'.^2);
    
    drift_post_mean(nn,:) = interp1(find(~isnan(fix_ids)),drift_post_mean(nn,~isnan(fix_ids)),1:NT);
    drift_post_std(nn,:) = interp1(find(~isnan(fix_ids)),drift_post_std(nn,~isnan(fix_ids)),1:NT);
    drift_post_mean(nn,isnan(drift_post_mean(nn,:))) = 0;
    
    %% construct drift-corrected X-mat
    fix_post_cor = nan(NT,1);
    fix_post_cor(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
    fix_post_cor = interp1(find(~isnan(fix_ids)),fix_post_cor(~isnan(fix_ids)),1:NT);
    fix_post_cor(isnan(fix_post_cor)) = 0;
    %back project drift (within-fixation) by sac_shift
    drift_post_cor = squeeze(drift_post_mean(nn,:));
    for ii = 1:n_fixs
        cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
        if length(cur_inds) > sac_shift
            drift_post_cor(cur_inds(1:end-sac_shift+1)) = drift_post_cor(cur_inds(sac_shift:end));
        end
    end
    drift_post_cor(isnan(drift_post_cor)) = 0;
    
    all_post_cor = round((fix_post_cor+drift_post_cor)) + max_Tshift + 1;
    
    %RECOMPUTE XMAT
    for i=1:NT
        all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_post_cor(i)};
    end
    cur_X{1} = create_time_embedding(all_shift_stimmat_up,stim_params_us);
    cur_X{1} = cur_X{1}(used_inds,use_kInds_up);
    
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
            tr_X,[],[],silent); %fit stimulus filters
        
        dit_mods_spkNL{nn+1}(cur_cell) = NMMfit_logexp_spkNL(dit_mods{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
        
        [newLL,~,new_prate] = NMMmodel_eval(dit_mods_spkNL{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
        [~,~,null_prate] = NMMmodel_eval(all_nullmod(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X(2:end));
        [dit_R2(nn+1,cur_cell),dit_dev(nn+1,cur_cell)] = pseudo_r2(Robs_mat(cur_tr_uset,cur_unit_ind),new_prate,null_prate);
        
        dit_LLimp(nn+1,cur_cell) = (newLL - null_LL(cur_cell))/log(2);
        
        fprintf('Original: %.5f  Prev: %.5f  New: %.5f\n',all_mod_LLimp(cur_cell),dit_LLimp(nn,cur_cell),dit_LLimp(nn+1,cur_cell));
    end
    
    %%
    if use_LOOXV == 1
        for xv = 1:length(su_inds)
            fprintf('Inferring drift corrections, XV %d of %d\n',xv,length(su_inds));
            
            %re-create LOOXV fixation-corrected Xmat
            all_fix_post_mean_cor(~isnan(fix_ids)) = squeeze(it_fix_post_mean_LOO(xv,end,fix_ids(~isnan(fix_ids))));
            all_fix_post_mean_cor = interp1(find(~isnan(fix_ids)),all_fix_post_mean_cor(~isnan(fix_ids)),1:NT);
            all_fix_post_mean_cor(isnan(all_fix_post_mean_cor)) = 0;
            all_fix_post_mean_cor = round(all_fix_post_mean_cor) + max_Tshift + 1;
            for i=1:NT
                all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_fix_post_mean_cor(i)};
            end
            all_Xmat_up_fixcor = create_time_embedding(all_shift_stimmat_up,stim_params_us);
            all_Xmat_up_fixcor = all_Xmat_up_fixcor(used_inds,:);
            
            %% PREPROCESS MODEL COMPONENTS
            cur_uset = setdiff(1:n_tr_chs,su_inds(xv));
            n_uset = length(cur_uset);
            filt_bank = zeros(n_uset,klen_us,n_squared_filts+1);
            lin_kerns = nan(n_uset,n_blocks);
            if use_sac_kerns
                sac_kerns = nan(n_uset,n_sac_bins);
                msac_kerns = nan(n_uset,n_sac_bins);
            end
            mod_spkNL_params = nan(n_uset,3);
            for ss = 1:n_uset
                cur_Xtargs = [dit_mods_LOO{xv,nn}(tr_set(cur_uset(ss))).mods(:).Xtarget];
                cur_k = [dit_mods_LOO{xv,nn}(tr_set(cur_uset(ss))).mods(cur_Xtargs == 1).filtK];
                n_used_filts = size(cur_k,2);
                filt_bank(ss,:,1:n_used_filts) = cur_k;
                mod_spkNL_params(ss,:) = dit_mods_spkNL_LOO{xv,nn}(tr_set(cur_uset(ss))).spk_NL_params;
                lin_kerns(ss,:) = dit_mods_LOO{xv,nn}(tr_set(cur_uset(ss))).mods(cur_Xtargs == 2).filtK;
                if use_sac_kerns
                    sac_kerns(ss,:) = dit_mods_LOO{xv,nn}(tr_set(cur_uset(ss))).mods(cur_Xtargs == 3).filtK;
                    msac_kerns(ss,:) = dit_mods_LOO{xv,nn}(tr_set(cur_uset(ss))).mods(cur_Xtargs == 4).filtK;
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
            for xx = 1:length(Dshifts)
                fprintf('Dshift %d of %d\n',xx,n_Dshifts);
                cur_stim_shift = all_Xmat_up_fixcor*Dshift_mat{xx};
                
                %outputs of stimulus models at current X-matrix shift
                gfuns = ones(NT,n_uset);
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
                
                frame_LLs(:,xx) = squeeze(nansum(Robs_mat(:,cur_uset).*log(pred_rate) - pred_rate,2));
            end
            
            %% INFER MICRO-SAC SEQUENCE
            lgamma = nan(NT,n_Dshifts);
            for ff = 1:n_fixs
                if mod(ff,100)==0
                    fprintf('Fixation %d of %d\n',ff,n_fixs);
                end
                
                tset = find(fix_ids==ff)';
                ntset = length(tset);
                nt_pts = ceil(ntset/drift_dsf);
                tset_inds = 1+floor((0:(ntset-1))/drift_dsf);
                talpha=zeros(nt_pts,n_Dshifts);
                tbeta = zeros(nt_pts,n_Dshifts);
                
                tpt_loc = ceil(1:drift_dsf:nt_pts*drift_dsf);
                tpt_loc(end) = ntset;
                
                cur_drift_post = post_mean_drift(tset);
                cur_LL_set = frame_LLs(tset,:);
                if mod(ntset,drift_dsf) ~= 0
                    dangling_pts = nt_pts*drift_dsf-ntset;
                    cur_LL_set = cat(1,cur_LL_set,zeros(dangling_pts,n_Dshifts));
                    cur_drift_post = cat(1,cur_drift_post,zeros(dangling_pts,1));
                end
                cur_LL_set = reshape(cur_LL_set,[drift_dsf nt_pts n_Dshifts]);
                cur_LL_set = squeeze(sum(cur_LL_set,1));
                
                cur_drift_post = drift_dsf*mean(reshape(cur_drift_post,[drift_dsf nt_pts]));
                
                if all(use_coils==0)
                    cur_lA = repmat(base_lA,[1 1 nt_pts]);
                else
                    cur_lA = nan(n_Dshifts,n_Dshifts,nt_pts);
                    for iii = 1:nt_pts
                        %                                 gprobs = diff(erf((Dshift_dx_edges+cur_drift_mean(iii))/(post_drift_sigma*sqrt(2)*drift_dsf)));
                        %                                 cur_lA(:,:,iii) = toeplitz(gprobs(n_Dshifts:end),gprobs(n_Dshifts:-1:1));
                        cdist = pdist2(Dshifts'*sp_dx + cur_drift_mean(iii),Dshifts'*sp_dx);
                        cur_lA(:,:,iii) = -cdist.^2/(2*(post_drift_sigma*drift_dsf)^2);
                    end
                    %             cur_lA = log(bsxfun(@rdivide,cur_lA,sum(cur_lA,2)));
                    %             cur_lA(cur_lA < eps) = eps;
                    cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2));
                end
                
                talpha(1,:) = drift_jump_prior + cur_LL_set(1,:);
                for t = 2:nt_pts
                    talpha(t,:) = logmulexp(talpha(t-1,:),cur_lA(:,:,t)) + cur_LL_set(t,:);
                end
                
                tbeta(end,:)=log(ones(1,n_Dshifts));
                for t = (nt_pts-1):-1:1
                    lf1 = tbeta(t+1,:) + cur_LL_set(t+1,:);
                    tbeta(t,:) = logmulexp(lf1,cur_lA(:,:,t)');
                end
                temp_gamma = talpha + tbeta;
                
                if drift_dsf > 1
                    if nt_pts > 1
                        int_gamma = interp1(tpt_loc,temp_gamma,1:ntset);
                        lgamma(tset,:) = int_gamma;
                    else
                        lgamma(tset,:) = repmat(temp_gamma,ntset,1);
                    end
                else
                    lgamma(tset,:) = temp_gamma;
                end
            end
            lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
            
            gamma = exp(lgamma);
            drift_post_mean_LOO(xv,nn,:) = sum(bsxfun(@times,gamma,Dshifts),2);
            drift_post_std_LOO(xv,nn,:) = sqrt(sum(bsxfun(@times,gamma,Dshifts.^2),2) - squeeze(drift_post_mean_LOO(xv,nn,:)).^2);
            drift_post_mean_LOO(xv,nn,:) = interp1(find(~isnan(fix_ids)),squeeze(drift_post_mean_LOO(xv,nn,~isnan(fix_ids))),1:NT);
            drift_post_std_LOO(xv,nn,:) = interp1(find(~isnan(fix_ids)),squeeze(drift_post_std_LOO(xv,nn,~isnan(fix_ids))),1:NT);
            drift_post_mean_cor = squeeze(drift_post_mean_LOO(xv,nn,:));
            drift_post_mean_cor(isnan(drift_post_mean_cor)) = 0;
            
            %% construct drift-corrected X-mat
            fix_post_cor = nan(NT,1);
            fix_post_cor(~isnan(fix_ids)) = squeeze(it_fix_post_mean_LOO(xv,end,fix_ids(~isnan(fix_ids))));
            fix_post_cor = interp1(find(~isnan(fix_ids)),fix_post_cor(~isnan(fix_ids)),1:NT);
            fix_post_cor(isnan(fix_post_cor)) = 0;
            
            %back project drift (within-fixation) by sac_shift
            drift_post_cor = squeeze(drift_post_mean_LOO(xv,nn,:));
            for ii = 1:n_fixs
                cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
                if length(cur_inds) > sac_shift
                    drift_post_cor(cur_inds(1:end-sac_shift+1)) = drift_post_cor(cur_inds(sac_shift:end));
                end
            end
            drift_post_cor(isnan(drift_post_cor)) = 0;
            
            all_post_cor = round(fix_post_cor+drift_post_cor') + max_Tshift + 1;
            
            %RECOMPUTE XMAT
            for i=1:NT
                all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_post_cor(i)};
            end
            cur_X{1} = create_time_embedding(all_shift_stimmat_up,stim_params_us);
            cur_X{1} = cur_X{1}(used_inds,use_kInds_up);
            
            %% REFIT XV CELLS
            silent = 1;
            for ss = 1:length(tr_set)
                cur_cell = tr_set(ss);
                fprintf('Drift LOO %d/%d, Refitting model %d of %d\n',xv,length(su_inds),ss,length(tr_set));
                cur_unit_ind = find(tr_set == cur_cell);
                cur_tr_uset = tr_inds(~isnan(Robs_mat(tr_inds,cur_unit_ind)));
                
                tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
                
                dit_mods_LOO{xv,nn+1}(cur_cell) = dit_mods{nn+1}(cur_cell);
                dit_mods_LOO{xv,nn+1}(cur_cell) = NMMfit_filters(dit_mods_LOO{xv,nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),...
                    tr_X,[],[],silent); %fit stimulus filters
                
                dit_mods_spkNL_LOO{xv,nn+1}(cur_cell) = NMMfit_logexp_spkNL(dit_mods_LOO{xv,nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
                
                [newLL,~,new_prate] = NMMmodel_eval(dit_mods_spkNL_LOO{xv,nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
                [~,~,null_prate] = NMMmodel_eval(all_nullmod(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X(2:end));
                [dit_R2_LOO(xv,nn+1,cur_cell),dit_dev_LOO(xv,nn+1,cur_cell)] = pseudo_r2(Robs_mat(cur_tr_uset,cur_unit_ind),new_prate,null_prate);
                
                dit_LLimp_LOO(xv,nn+1,cur_cell) = (newLL - null_LL(cur_cell))/log(2);
                
                %                 fprintf('Original: %.5f  Prev: %.5f  New: %.5f\n',all_mod_LLimp(cur_cell),dit_LLimp(nn,cur_cell),dit_LLimp(nn+1,cur_cell));
            end
        end
    end
end

%% SAVE EYE-TRACKING RESULTS
et_params = struct('beg_buffer',beg_buffer,'end_buffer',end_buffer,'min_trial_dur',min_trial_dur,'bar_ori',bar_ori,...
    'use_nPix',use_nPix,'flen',flen,'dt',dt,'drift_jump_sigma',drift_jump_sigma,'drift_prior_sigma',drift_prior_sigma,...
    'fix_prior_sigma',fix_prior_sigma,'fix_noise_sigma',fix_noise_sigma,'drift_noise_sigma',drift_noise_sigma,...
    'drift_dsf',drift_dsf,'n_fix_inf_it',n_fix_inf_it,'n_drift_inf_it',n_drift_inf_it,'use_sac_kerns',use_sac_kerns,'shifts',shifts,...
    'use_measured_pos',use_measured_pos,'sac_bincents',sac_bincents,'spatial_usfac',spatial_usfac,'sac_shift',sac_shift,'use_coils',use_coils);

et_used_inds = used_inds;
et_tr_set = tr_set;
et_tr_trials = tr_trials;
et_xv_trials = xv_trials;
et_saccade_inds = saccade_start_inds;
cd(anal_dir);
save(anal_name,'it_*','drift_post_*','fix_ids','dit_*','et_used_inds','et_tr_set','et_saccade_inds','et_params','et_tr_trials','et_xv_trials');

%%
fin_fix_corr = nan(NT,1);
fin_fix_std = nan(NT,1);
fin_fix_corr(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
fin_fix_std(~isnan(fix_ids)) = it_fix_post_std(end,fix_ids(~isnan(fix_ids)));
fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);

fin_fix_corr = fin_fix_corr*sp_dx;
fin_fix_std = fin_fix_std*sp_dx;

fin_drift_corr = drift_post_mean(end,:)*sp_dx;
fin_drift_std = drift_post_std(end,:)*sp_dx;
min_fix_dur = 0.15;
fix_inds = [];
long_fix_inds = [];
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > sac_shift
        fin_drift_corr(cur_inds(1:end-sac_shift+1)) = fin_drift_corr(cur_inds(sac_shift:end));
        fin_drift_std(cur_inds(1:end-sac_shift+1)) = fin_drift_std(cur_inds(sac_shift:end));
    end
    fix_inds = [fix_inds cur_inds];
    if length(cur_inds)*dt >= min_fix_dur
        long_fix_inds = [long_fix_inds cur_inds];
    end
end

fin_tot_corr = fin_fix_corr + fin_drift_corr;
fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2);
% fin_tot_corr = interp1(find(~isnan(fix_ids)),fin_tot_corr(~isnan(fix_ids)),1:NT);
% fin_tot_std = interp1(find(~isnan(fix_ids)),fin_tot_std(~isnan(fix_ids)),1:NT);

if use_measured_pos==1
    fin_tot_corr = fin_tot_corr + init_eyepos(used_inds)';
end
%%
