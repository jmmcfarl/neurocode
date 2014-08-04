clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

% dir_prefix = '/home/james';
dir_prefix = '/Volumes/james';

%run on 86 [0, 90]; 
Expt_num = 88;    
Expt_name = sprintf('G%.3d',Expt_num);
data_dir = [dir_prefix '/Data/bruce/' Expt_name];
cd(data_dir);
fig_dir = [dir_prefix '/Analysis/bruce/ET_final/'];

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = [dir_prefix '/Analysis/bruce/' Expt_name '/ET_final/'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

bar_ori = 0;
use_eye = 'left';

if bar_ori == 0
    mod_data_name = 'eyecorr_hbar_mods3';
    old_anal_name = 'eyecorr_hbar3';
    anal_name = 'eyecorr_hbar_highres';
    hrmod_data_name = 'eyecorr_hbar_mods_hres';
else
    mod_data_name = 'eyecorr_vbar_mods3';
    old_anal_name = 'eyecorr_vbar3';
    anal_name = 'eyecorr_vbar_highres';
    hrmod_data_name = 'eyecorr_vbar_mods_hres';
end
mod_data_name = strcat(use_eye,'_',mod_data_name);
 old_anal_name = strcat(use_eye,'_', old_anal_name);
anal_name = strcat(use_eye,'_',anal_name);
 hrmod_data_name = strcat(use_eye,'_', hrmod_data_name);


recompute_init_mods = 0;
use_sac_kerns = 1;
use_LOOXV = 0;
use_coils = [0 0]; %[L R]


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

% use training set defined from monocular expt
if bar_ori == 0
base_name = 'monoc_eyecorr_hbar';
else
    base_name = 'monoc_eyecorr_vbar';
end
cd(anal_dir);
temp = load(base_name,'et_tr_set');
tr_set = temp.et_tr_set;
%%
xv_frac = 0;

flen = 12;
use_nPix = 16;
n_drift_inf_it = 1; %2

drift_noise_sigma = 0.003;
drift_prior_sigma = 0.004; %.004 may be best here
drift_jump_sigma = 0.1; %0.05 start
drift_dsf = 2;

min_trial_dur = 0.75;

spatial_usfac = 4;
old_usfac = 2;

%%
eps = -1e3;

stim_fs = 100; %in Hz
dt = 0.01;
full_nPix = 36;
stim_params = NIMcreate_stim_params([flen full_nPix],dt);
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
% temp_dx = [-2*max_shift:dshift:2*max_shift];
% shift_dx_edges = [temp_dx*sp_dx-dshift*sp_dx/2 temp_dx(end)*sp_dx+dshift*sp_dx/2];
% ovshift_dx_edges = [shifts*sp_dx-dshift*sp_dx/2 shifts(end)*sp_dx+dshift*sp_dx/2];

max_Dshift = round(8*spatial_usfac);
Dshifts = -max_Dshift:dshift:max_Dshift;
n_Dshifts = length(Dshifts);
zero_Dframe = find(Dshifts==0);
% temp_dx = [-2*max_Dshift:dshift:2*max_Dshift];
% Dshift_dx_edges = [temp_dx*sp_dx-dshift*sp_dx/2 temp_dx(end)*sp_dx+dshift*sp_dx/2];
% ovDshift_dx_edges = [Dshifts*sp_dx-dshift*sp_dx/2 Dshifts(end)*sp_dx+dshift*sp_dx/2];

max_Tshift = max_shift + max_Dshift;
Tshifts = -max_Tshift:dshift:max_Tshift;

%% load overall su data
% LOAD REFCLUSTERS
cluster_dir = [dir_prefix '/Analysis/bruce/' Expt_name '/clustering'];
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
% all_Xmat = [];
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
            if strcmp(use_eye,'left')
            cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            elseif strcmp(use_eye,'right')
            cur_stim_mat = double(right_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            else
                error('Unsupported eye specified');
            end
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
    all_stimmat_up = zeros(size(all_stim_mat,1),full_nPix_us);
    for ii = 1:size(all_stim_mat,2)
        all_stimmat_up(:,2*(ii-1)+1) = all_stim_mat(:,ii);
        all_stimmat_up(:,2*(ii-1)+2) = all_stim_mat(:,ii);
    end
elseif spatial_usfac == 4
     all_stimmat_up = zeros(size(all_stim_mat,1),full_nPix_us);
    for ii = 1:size(all_stim_mat,2)
        all_stimmat_up(:,4*(ii-1)+1) = all_stim_mat(:,ii);
        all_stimmat_up(:,4*(ii-1)+2) = all_stim_mat(:,ii);
        all_stimmat_up(:,4*(ii-1)+3) = all_stim_mat(:,ii);
        all_stimmat_up(:,4*(ii-1)+4) = all_stim_mat(:,ii);
    end   
elseif spatial_usfac == 1
    all_stimmat_up = all_stim_mat;
else
    error('Unsupported spatial Up-sample factor!');
end
stim_params_us = NMMcreate_stim_params([flen full_nPix_us],dt);
all_Xmat = create_time_embedding(all_stim_mat,stim_params);
all_Xmat_us = create_time_embedding(all_stimmat_up,stim_params_us);

%% select submatrix with central pixels
[Xinds,Tinds] = meshgrid(1:full_nPix,1:flen);
buffer_pix = floor((full_nPix - use_nPix)/2);
cur_use_pix = (1:use_nPix) + buffer_pix;
use_kInds = find(ismember(Xinds(:),cur_use_pix));

[Xinds_up,Tinds] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
if spatial_usfac > 1
%     use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1)-1/spatial_usfac & Xinds_up(:) <= cur_use_pix(end));
    use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1)-2/spatial_usfac & Xinds_up(:) <= cur_use_pix(end) + 1/spatial_usfac);
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
fix_ids = nan(NT,1);
pfix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
    cur_inds = pfix_start_inds(ii):(pfix_stop_inds(ii));
    pfix_ids(cur_inds) = ii;
end

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
init_stim_params = NMMcreate_stim_params([flen use_nPix],dt);
init_stim_params(2) = NMMcreate_stim_params([n_blocks 1],dt);

fin_stim_params(1) = NMMcreate_stim_params([flen use_nPix_us],dt);
fin_stim_params(2) = NMMcreate_stim_params([n_blocks 1],dt);
if use_sac_kerns
    fin_stim_params(3) = NMMcreate_stim_params(n_sac_bins,dt);
    fin_stim_params(4) = NMMcreate_stim_params(n_sac_bins,dt);
end

null_stim_params = fin_stim_params(2:end);

block_L2 = 1;
silent = 1;
sac_d2t = 100;

if expt_dds(cur_block_set(1)) == 67
    base_lambda_d2XT = 60;
    base_lambda_L1 = 10;
else
    base_lambda_d2XT = 10;
    base_lambda_L1 = 4;
end
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

%%
cd(anal_dir);

fprintf('Loading pre-computed initial models\n');
load(mod_data_name);
load(old_anal_name,'dit_mods','et_tr_set','it_fix_post_mean','drift_post_mean','it_fix_post_std');
old_best_mods = dit_mods{end};
tr_set = et_tr_set;

best_fix_cor = it_fix_post_mean(end,:)*spatial_usfac/old_usfac;
best_fix_std = it_fix_post_std(end,:)*spatial_usfac/old_usfac;
best_drift_cor = drift_post_mean(end,:)*spatial_usfac/old_usfac;

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

%% DEFINE POINTS TO ALLOW RAPID CHANGES (TRIAL STARTS AFTER SACCADES)

%generate shift matrices. Must be applied to the stimulus (not the filters)
It = speye(flen);
shift_mat = cell(n_shifts,1);
for xx = 1:n_shifts
    temp = spdiags( ones(full_nPix_us,1), -shifts(xx), full_nPix_us, full_nPix_us);
    temp = kron(temp,It);
    shift_mat{xx} = temp(:,use_kInds_up);
end
Dshift_mat = cell(n_Dshifts,1);
for xx = 1:n_Dshifts
    temp = spdiags( ones(full_nPix_us,1), -Dshifts(xx), full_nPix_us, full_nPix_us);
    temp = kron(temp,It);
    Dshift_mat{xx} = temp(:,use_kInds_up);
end

%shift matrices for images (non-time-embedded)    
n_Tshifts = length(Tshifts);
Tshift_mat = cell(n_Tshifts,1);
for xx = 1:n_Tshifts
    Tshift_mat{xx} = spdiags( ones(full_nPix_us,1), -Tshifts(xx), full_nPix_us, full_nPix_us);
end

%%
warning('OFF','stats:statrobustfit:IterationLimit');

measured_eyepos = [corrected_eye_vals_interp(:,2) corrected_eye_vals_interp(:,4)];
max_sim_pos = max_shift*sp_dx;

% measured_eyepos = tanh(measured_eyepos/max_sim_pos)*max_sim_pos;
measured_eyepos(measured_eyepos > max_sim_pos) = max_sim_pos;
measured_eyepos(measured_eyepos < -max_sim_pos) = -max_sim_pos;

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
        back_measured_drift(cur_inds(sac_shift+1:end),:) = measured_drift(cur_inds(1:(end-sac_shift)),:);
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
% post_mean_drift(isnan(post_mean_drift)) = 0;

% % fix_prior = diff(erf(ovshift_dx_edges/(fix_prior_sigma*sqrt(2))));
% % fix_prior = log(fix_prior/sum(fix_prior));
% % fix_prior(fix_prior < eps) = eps;
% fix_prior = -(shifts*sp_dx).^2./(2*fix_prior_sigma^2);
% fix_prior = fix_prior - logsumexp(fix_prior);

% %overall prior on shifts
% temp_edges = [Dshifts*sp_dx-dshift*sp_dx/2 Dshifts(end)*sp_dx + dshift*sp_dx/2];
% drift_jump_prior = diff(erf(temp_edges/(sqrt(2)*drift_jump_sigma)));
% drift_jump_prior = log(drift_jump_prior/sum(drift_jump_prior));
% drift_jump_prior(drift_jump_prior < eps) = eps;
drift_jump_prior = -(Dshifts*sp_dx).^2./(2*drift_jump_sigma^2);
drift_jump_prior = drift_jump_prior - logsumexp(drift_jump_prior);

% gprobs = diff(erf((Dshift_dx_edges)/(post_drift_sigma*sqrt(2)*drift_dsf)));
% base_lA = toeplitz(gprobs(n_Dshifts:end),gprobs(n_Dshifts:-1:1));
% base_lA = log(bsxfun(@rdivide,base_lA,sum(base_lA,2)));
% base_lA(base_lA < eps) = eps;
cdist = squareform(pdist(Dshifts'*sp_dx));
base_lA = -cdist.^2/(2*(post_drift_sigma*drift_dsf)^2);
base_lA = bsxfun(@minus,base_lA,logsumexp(base_lA,2)); %normalize


%% construct drift-corrected X-mat
fix_post_cor = nan(NT,1);
fix_post_cor(~isnan(fix_ids)) = best_fix_cor(fix_ids(~isnan(fix_ids)));
fix_post_cor = interp1(find(~isnan(fix_ids)),fix_post_cor(~isnan(fix_ids)),1:NT);
fix_post_cor(isnan(fix_post_cor)) = 0;
%back project drift (within-fixation) by sac_shift
drift_post_cor = best_drift_cor;

    for ii = 1:length(trial_start_inds)
        cur_inds = trial_start_inds(ii):trial_end_inds(ii);
        drift_post_cor(cur_inds(1:end-sac_shift)) = drift_post_cor(cur_inds(sac_shift+1:end));
    end
    drift_post_cor = interp1(find(~isnan(fix_ids)),drift_post_cor(~isnan(fix_ids)),1:NT);

% for ii = 1:n_fixs
%     cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
%     if length(cur_inds) > sac_shift
%         drift_post_cor(cur_inds(1:end-sac_shift+1)) = drift_post_cor(cur_inds(sac_shift:end));
%     end
% end
drift_post_cor(isnan(drift_post_cor)) = 0;

all_post_cor = round((fix_post_cor+drift_post_cor)) + max_Tshift + 1;

%RECOMPUTE XMAT
all_shift_stimmat_up = all_stimmat_up;
for i=1:NT
    all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_post_cor(i)};
end
all_Xmat_cor = create_time_embedding(all_shift_stimmat_up,stim_params_us);

%%
if ~exist(['./' hrmod_data_name '.mat'],'file') || recompute_init_mods == 1

for ss = 1:n_tr_chs
    fprintf('Computing base LLs for Unit %d of %d\n',ss,n_tr_chs);
    cur_tr_inds = tr_inds(~isnan(Robs_mat(tr_inds,ss)));
    tr_NT = length(cur_tr_inds);
    Robs = Robs_mat(cur_tr_inds,ss);
    null_mod = all_nullmod(tr_set(ss));
    if ~isempty(cur_tr_inds) && nansum(Robs) > 0
        
        tr_X{1} = all_Xmat_cor(used_inds(cur_tr_inds),use_kInds_up);
        tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
        
        if use_sac_kerns
            tr_X{3} = Xsac(cur_tr_inds,:);
            tr_X{4} = Xmsac(cur_tr_inds,:);
        end
        
        %spatial up-sampling of filter estimates
        base_filts = reshape([old_best_mods(tr_set(ss)).mods(find(init_Xtargs == 1)).filtK],[flen use_nPix*2 n_squared_filts+1]);
        base_filts_up = zeros(flen,use_nPix_us,n_squared_filts+1);
        for ii = 1:use_nPix*2
            base_filts_up(:,2*(ii-1)+1,:) = 0.5*base_filts(:,ii,:);
            base_filts_up(:,2*(ii-1)+2,:) = 0.5*base_filts(:,ii,:);
        end
        base_filts_up = reshape(base_filts_up,use_nPix_us*flen,n_squared_filts+1);
        
        init_filts{end} = old_best_mods(tr_set(ss)).mods(find(init_Xtargs==2)).filtK;
        for ii = 1:n_squared_filts+1
            init_filts{ii} = base_filts_up(:,ii);
        end
        gqm2 = NMMinitialize_model(fin_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs,init_filts);
        [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm2,Robs,tr_X);
        gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
        gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
        if use_sac_kerns
            gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,3,zeros(n_sac_bins,1),sac_reg_params); %sac term
            gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,4,zeros(n_sac_bins,1),sac_reg_params); %nsac term
        end
        all_mod_fits(tr_set(ss)) = gqm2;
        all_mod_fits(tr_set(ss)).spk_NL_params(1) = old_best_mods(tr_set(ss)).spk_NL_params(1);
        all_mod_fits(tr_set(ss)) = NMMfit_filters(all_mod_fits(tr_set(ss)),Robs,tr_X);
        
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
clear all_Xmat_cor
su_inds = find(all_mod_SU(tr_set) > 0);

%%
%back-project saccade-times
all_fix_post_mean_cor = nan(NT,1);
all_fix_post_mean_cor(~isnan(fix_ids)) = best_fix_cor(fix_ids(~isnan(fix_ids)));
all_fix_post_mean_cor = interp1(find(~isnan(fix_ids)),all_fix_post_mean_cor(~isnan(fix_ids)),1:NT);
all_fix_post_mean_cor(isnan(all_fix_post_mean_cor)) = 0;

all_fix_post_mean_cor = round(all_fix_post_mean_cor) + max_Tshift + 1;
all_shift_stimmat_up = all_stimmat_up;
for i=1:NT
    all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_fix_post_mean_cor(i)};
end
all_Xmat_up_fixcor = create_time_embedding(all_shift_stimmat_up,stim_params_us);
all_Xmat_up_fixcor = all_Xmat_up_fixcor(used_inds,:);

%% NOW INFER DRIFT CORRECTIONS

dit_mods{1} = all_mod_fits;
dit_mods_spkNL{1} = all_mod_fits_withspkNL;
dit_LLimp(1,:) = all_mod_LLimp;
dit_R2(1,:) = all_mod_R2;
dit_dev(1,:) = all_mod_dev;
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
        all_fix_post_mean_cor(~isnan(fix_ids)) = best_fix_cor(fix_ids(~isnan(fix_ids)));
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
    filt_bank = zeros(n_tr_chs,klen_us,n_squared_filts+1);
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
        
%         tset = find(fix_ids==ff)';
        tset = find(pfix_ids==ff)';
        ntset = length(tset);
        if ntset > drift_dsf
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
            cur_drift_mean = cat(1,cur_drift_mean,nan(dangling_pts,1));
        end
        cur_LL_set = reshape(cur_LL_set,[drift_dsf nt_pts n_Dshifts]);
        cur_LL_set = squeeze(sum(cur_LL_set,1));
        
        cur_drift_mean = drift_dsf*nanmean(reshape(cur_drift_mean,[drift_dsf nt_pts]));
        
        if all(use_coils==0)
            cur_lA = repmat(base_lA,[1 1 nt_pts]);
        else
            cur_lA = nan(n_Dshifts,n_Dshifts,nt_pts);
            for iii = 1:nt_pts
                if ~isnan(cur_drift_mean(iii))
%                 gprobs = diff(erf((Dshift_dx_edges+cur_drift_mean(iii))/(post_drift_sigma*sqrt(2)*drift_dsf)));
%                 cur_lA(:,:,iii) = toeplitz(gprobs(n_Dshifts:end),gprobs(n_Dshifts:-1:1));
                    cdist = pdist2(Dshifts'*sp_dx + cur_drift_mean(iii),Dshifts'*sp_dx);
                    cur_lA(:,:,iii) = -cdist.^2/(2*(post_drift_sigma*drift_dsf)^2);
                else
                   cur_lA(:,:,iii) = base_lA; 
                end
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
            tbeta(t,:) = logmulexp(lf1,cur_lA(:,:,t+1)');
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
    end
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    gamma = exp(lgamma);
    drift_post_mean(nn,:) = sum(bsxfun(@times,gamma,Dshifts),2);
    drift_post_std(nn,:) = sqrt(sum(bsxfun(@times,gamma,Dshifts.^2),2) - drift_post_mean(nn,:)'.^2);
    
    drift_post_mean(nn,:) = interp1(find(~isnan(pfix_ids)),drift_post_mean(nn,~isnan(pfix_ids)),1:NT);
    drift_post_std(nn,:) = interp1(find(~isnan(pfix_ids)),drift_post_std(nn,~isnan(pfix_ids)),1:NT);
%     drift_post_mean(nn,:) = interp1(find(~isnan(fix_ids)),drift_post_mean(nn,~isnan(fix_ids)),1:NT);
%     drift_post_std(nn,:) = interp1(find(~isnan(fix_ids)),drift_post_std(nn,~isnan(fix_ids)),1:NT);
    drift_post_mean(nn,isnan(drift_post_mean(nn,:))) = 0;
     
    %% construct drift-corrected X-mat
%     fix_post_cor = nan(NT,1);
%     fix_post_cor(~isnan(fix_ids)) = best_fix_cor(fix_ids(~isnan(fix_ids)));
%     fix_post_cor = interp1(find(~isnan(fix_ids)),fix_post_cor(~isnan(fix_ids)),1:NT);
%     fix_post_cor(isnan(fix_post_cor)) = 0;
%     drift_post_cor = squeeze(drift_post_mean(nn,:));
% 
%     for ii = 1:length(trial_start_inds)
%         cur_inds = trial_start_inds(ii):trial_end_inds(ii);
%         drift_post_cor(cur_inds(1:end-sac_shift)) = drift_post_cor(cur_inds(sac_shift+1:end));
%     end
%     drift_post_cor = interp1(find(~isnan(fix_ids)),drift_post_cor(~isnan(fix_ids)),1:NT);
% 
% %     back project drift (within-fixation) by sac_shift
% %     for ii = 1:n_fixs
% %         cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
% %         if length(cur_inds) > sac_shift
% %             drift_post_cor(cur_inds(1:end-sac_shift+1)) = drift_post_cor(cur_inds(sac_shift:end));
% %         end
% %     end
%     drift_post_cor(isnan(drift_post_cor)) = 0;
%     
%     all_post_cor = round((fix_post_cor+drift_post_cor)) + max_Tshift + 1;
%     
%     %RECOMPUTE XMAT
%     for i=1:NT
%         all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_post_cor(i)};
%     end
%     cur_X{1} = create_time_embedding(all_shift_stimmat_up,stim_params_us);
%     cur_X{1} = cur_X{1}(used_inds,use_kInds_up);
%     
%     cur_X{2} = Xblock(used_inds,:);
%     if use_sac_kerns
%         cur_X{3} = Xsac;
%         cur_X{4} = Xmsac;
%     end

    %% REFIT ALL CELLS    
%     silent = 1;
%     for ss = 1:length(tr_set)
%         cur_cell = tr_set(ss);
%         fprintf('Refitting model for tr cell %d of %d\n',ss,length(tr_set));
%         cur_unit_ind = find(tr_set == cur_cell);
%         cur_tr_uset = tr_inds(~isnan(Robs_mat(tr_inds,cur_unit_ind)));
%         
%         tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
%         
%         dit_mods{nn+1}(cur_cell) = dit_mods{nn}(cur_cell);
%         dit_mods{nn+1}(cur_cell) = NMMfit_filters(dit_mods{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),...
%             tr_X,[],[],silent); %fit stimulus filters
%         
%         dit_mods_spkNL{nn+1}(cur_cell) = NMMfit_logexp_spkNL(dit_mods{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
%         
%         [newLL,~,new_prate] = NMMmodel_eval(dit_mods_spkNL{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
%         [~,~,null_prate] = NMMmodel_eval(all_nullmod(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X(2:end));
%         [dit_R2(nn+1,cur_cell),dit_dev(nn+1,cur_cell)] = pseudo_r2(Robs_mat(cur_tr_uset,cur_unit_ind),new_prate,null_prate);
%          
%         dit_LLimp(nn+1,cur_cell) = (newLL - null_LL(cur_cell))/log(2);
%         
%         fprintf('Original: %.5f  Prev: %.5f  New: %.5f\n',all_mod_LLimp(cur_cell),dit_LLimp(nn,cur_cell),dit_LLimp(nn+1,cur_cell));
%     end
end

%% SAVE EYE-TRACKING RESULTS
et_params = struct('beg_buffer',beg_buffer,'end_buffer',end_buffer,'min_trial_dur',min_trial_dur,'bar_ori',bar_ori,...
    'use_nPix',use_nPix,'flen',flen,'dt',dt,'drift_jump_sigma',drift_jump_sigma,'drift_prior_sigma',drift_prior_sigma,'drift_noise_sigma',drift_noise_sigma,...
    'drift_dsf',drift_dsf,'n_drift_inf_it',n_drift_inf_it,'use_sac_kerns',use_sac_kerns,'shifts',shifts,...
    'sac_bincents',sac_bincents,'spatial_usfac',spatial_usfac,'sac_shift',sac_shift,'use_coils',use_coils);

et_used_inds = used_inds;
et_tr_set = tr_set;
et_tr_trials = tr_trials;
et_xv_trials = xv_trials;
et_saccade_inds = saccade_start_inds;
cd(anal_dir);
save(anal_name,'best_fix_*','drift_post_*','fix_ids','dit_*','et_used_inds','et_tr_set','et_saccade_inds','et_params','et_tr_trials','et_xv_trials');

%%
fin_fix_corr = nan(NT,1);
fin_fix_std = nan(NT,1);
fin_fix_corr(~isnan(fix_ids)) = best_fix_cor(fix_ids(~isnan(fix_ids)));
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
fin_fix_std(~isnan(fix_ids)) = best_fix_std(fix_ids(~isnan(fix_ids)));
fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);
% fin_fix_corr(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)))*2;
% fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
% fin_fix_std(~isnan(fix_ids)) = it_fix_post_std(end,fix_ids(~isnan(fix_ids)))*2;
% fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);

fin_fix_corr = fin_fix_corr*sp_dx;
fin_fix_std = fin_fix_std*sp_dx;

fin_drift_corr = drift_post_mean(end,:)*sp_dx;
fin_drift_std = drift_post_std(end,:)*sp_dx;
min_fix_dur = 0.15;
fix_inds = [];
long_fix_inds = [];
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
%     if length(cur_inds) > sac_shift
% %         fin_drift_corr(cur_inds(1:end-sac_shift+1)) = fin_drift_corr(cur_inds(sac_shift:end));
% %         fin_drift_std(cur_inds(1:end-sac_shift+1)) = fin_drift_std(cur_inds(sac_shift:end));
%         fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
%         fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
%     end
    fix_inds = [fix_inds cur_inds];
    if length(cur_inds)*dt >= min_fix_dur
        long_fix_inds = [long_fix_inds cur_inds];
    end
end

n_trials = length(trial_start_inds);
for ii = 1:n_trials
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
    fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
end

%within-saccade estimates use linear interpolation
fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT);
fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);

fin_tot_corr = fin_fix_corr + fin_drift_corr;
fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2);
% fin_tot_corr = interp1(find(~isnan(fix_ids)),fin_tot_corr(~isnan(fix_ids)),1:NT);
% fin_tot_std = interp1(find(~isnan(fix_ids)),fin_tot_std(~isnan(fix_ids)),1:NT);


%%
smooth_eyepos = corrected_eye_vals_interp;

%smooth out fast transients in eye signal
eye_smooth_sig = round(0.025/dt);
interp_inds = [];
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > eye_smooth_sig*5;
        smooth_eyepos(used_inds(cur_inds),1) = jmm_smooth_1d_cor(smooth_eyepos(used_inds(cur_inds),1),eye_smooth_sig,2);
        smooth_eyepos(used_inds(cur_inds),2) = jmm_smooth_1d_cor(smooth_eyepos(used_inds(cur_inds),2),eye_smooth_sig,2);
    end
    interp_inds = [interp_inds; cur_inds'];
end
interp_inds = unique(interp_inds);
smooth_eyepos(used_inds,:) = interp1(used_inds(interp_inds),smooth_eyepos(used_inds(interp_inds),:),used_inds);

%%

measured_seqL = corrected_eye_vals_interp(used_inds,2);
measured_seqR = corrected_eye_vals_interp(used_inds,4);
% measured_seqL = smooth_eyepos(used_inds,2);
% measured_seqR = smooth_eyepos(used_inds,4);

min_fix_dur = 0.15;
inferred_drift = nan(size(fin_tot_corr));
measured_driftL = nan(size(fin_tot_corr));
measured_driftR = nan(size(fin_tot_corr));
inferred_fix_avg = nan(n_fixs,1);
measured_fix_avgL = nan(n_fixs,1);
measured_fix_avgR = nan(n_fixs,1);
fix_inds = [];
long_fix_inds = [];
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    cur_inf = fin_tot_corr(cur_inds);
    inferred_fix_avg(ii) = median(fin_tot_corr(cur_inds));
    inferred_drift(cur_inds) = cur_inf - inferred_fix_avg(ii);
        measured_fix_avgL(ii) = median(measured_seqL(cur_inds));
        measured_fix_avgR(ii) = median(measured_seqR(cur_inds));
    measured_driftL(cur_inds) = measured_seqL(cur_inds) - measured_fix_avgL(ii);
    measured_driftR(cur_inds) = measured_seqR(cur_inds) - measured_fix_avgR(ii);
    fix_inds = [fix_inds; cur_inds(:)];
        %         inferred_drift_avg(cur_inds) = mean(drift_post_mean(end,cur_inds))*sp_dx;
    if length(cur_inds)*dt >= min_fix_dur
        long_fix_inds = [long_fix_inds; cur_inds(:)];
    end
end

u = find(~isnan(measured_driftL) & ~isnan(inferred_drift));
[drift_corrs,drif_pvals] = corr([measured_driftL(u)' measured_driftR(u)'],inferred_drift(u)','type','spearman');
u = find(~isnan(measured_fix_avgL) & ~isnan(inferred_fix_avg));
[fix_corrs,fix_pvals] = corr([measured_fix_avgL(u) measured_fix_avgR(u)],inferred_fix_avg(u),'type','spearman');
[tot_corrs,tot_pvals] = corr([measured_seqL measured_seqR],fin_tot_corr','type','spearman');
[tot_corrs_pear,tot_pvals_pear] = corr([measured_seqL measured_seqR],fin_tot_corr','type','pearson');

%%
close all

print_on = false;
fig_width = 4.86; %3.27 4.86 6.83
rel_height = 0.6;

yl = [-0.4 0.4];

n_trials = length(unique(all_trialvec));
H = figure();
for tt = 1:n_trials
% for tt = [30 31 74 113 121 179 223 239 436 522 556 563 630 658 688 736 746 757 759 913 930 949] %G086
% for tt = [179] %G086
    uu = find(all_trialvec(used_inds) == tt);
    if ~isempty(uu)
        bt = all_t_axis(used_inds(uu(1)));
        et = all_t_axis(used_inds(uu(end)));
        dur = et-bt;
        if dur > 3.5
           cur_sac_inds = find(ismember(saccade_start_inds,uu));
            rel_sac_start_times = all_t_axis(used_inds(saccade_start_inds(cur_sac_inds))) - bt;
            rel_sac_end_times = all_t_axis(used_inds(saccade_stop_inds(cur_sac_inds))) - bt;            
 
            hold on
            h2=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_fix_corr(uu),fin_fix_std(uu),{'color','m'});
            h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_tot_corr(uu),fin_tot_std(uu),{'color','k'});
            h3=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),2),'r','linewidth',2);
            line([0 dur],[0 0],'color','k','linestyle','--');
            
            
                        cur_fix_starts = find(ismember(fix_start_inds,uu));
            plot(all_t_axis(used_inds(fix_start_inds(cur_fix_starts)))-bt,zeros(size(cur_fix_starts)),'ko');
                        cur_fix_stops = find(ismember(fix_stop_inds,uu));
            plot(all_t_axis(used_inds(fix_stop_inds(cur_fix_stops)))-bt,zeros(size(cur_fix_stops)),'mo');
            plot(all_t_axis(used_inds(uu))-bt,fin_drift_corr(1,uu),'g');
            plot(all_t_axis(used_inds(uu))-bt,fin_drift_std(1,uu),'c');

            xlim([0 dur]);
            ylim(yl);
            xlabel('Time (s)','fontsize',10);
            ylabel('Orthoganol position (deg)','fontsize',10);
            title(sprintf('Trial %d',tt));
            xlabel('Time (s)');
            ylabel('Orthoganol position (deg)');
            title(sprintf('Trial %d',tt));
            for ii = 1:length(rel_sac_start_times)
                line(rel_sac_start_times([ii ii]),yl,'color',[0.2 0.8 0.2],'linestyle','--');
                %                 line(rel_sac_end_times([ii ii]),yl,'color',[0.2 0.8 0.2],'linestyle','--');
            end
                delete(h1.edge);
            
            if print_on
                delete(h1.edge);
                figufy(H);
                fname = [fig_dir sprintf('Example_loosetrace_%s_T%d',Expt_name,tt)];
                exportfig(H,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
                close(H);
            else
                pause
                clf(H);
            end
        end
    end
end

