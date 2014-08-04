clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

%run on 86 [0, 90];
Expt_num = 86;
Expt_name = sprintf('G%.3d',Expt_num);
data_dir = ['~/Data/bruce/' Expt_name];
cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final/'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

bar_ori = 0;
sim_gain = 1;

if sim_gain >= 1
    gname = num2str(sim_gain);
elseif sim_gain == 0.5
    gname = 'p5';
elseif sim_gain == 0.25
    gname = 'p25';
else
    error('unsupported');
end
sim_data_name = ['monoc_eyecorr_hbar_simdata_g' gname]';
% mod_data_name = ['monoc_eyecorr_hbar_sim4_mods_g' gname];
% data_name = ['monoc_eyecorr_hbar_sim4_g' gname];
% data_name_hres = ['monoc_eyecorr_hbar_hres_sim4_g' gname];
mod_data_name = ['monoc_eyecorr_hbar_sim4_mods_g' gname '_nosac'];
data_name = ['monoc_eyecorr_hbar_sim4_g' gname '_nosac'];
data_name_hres = ['monoc_eyecorr_hbar_hres_sim4_g' gname '_nosac'];

do_highres_est = true;
use_sac_kerns = 0;
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

recompute_init_mods = 1;

switch sim_gain
    case 4
        n_fix_inf_it = 5; 
        n_drift_inf_it = 3;
        n_HR_inf_it = 1;
        hr_spatial_usfac = 4;
spatial_usfac = 2;
    case 2
        n_fix_inf_it = 4; 
        n_drift_inf_it = 2;
        n_HR_inf_it = 1;
        hr_spatial_usfac = 4;
spatial_usfac = 2;
    case 1
        n_fix_inf_it = 3; 
        n_drift_inf_it = 1; 
        n_HR_inf_it = 1;
        hr_spatial_usfac = 4;
spatial_usfac = 2;
   case 0.5
        n_fix_inf_it = 2; 
        n_drift_inf_it = 1;
        n_HR_inf_it = 1;
        hr_spatial_usfac = 8;
spatial_usfac = 4;
    case 0.25
        n_fix_inf_it = 2; 
        n_drift_inf_it = 1; 
        n_HR_inf_it = 1;
        hr_spatial_usfac = 8;
spatial_usfac = 4;
end
%%
xv_frac = 0;

flen = 12;
use_nPix = 16;

fix_prior_sigma = 0.15*sim_gain;
fix_noise_sigma = 0.1*sim_gain;
drift_noise_sigma = 0.003*sim_gain;
drift_prior_sigma = 0.004*sim_gain; %.004 may be best here
drift_jump_sigma = 0.075*sim_gain; %0.05 start
% drift_jump_sigma = 0.15*sim_gain; %0.05 start
drift_dsf = 3;

% %apply ao minimum threshold to the drift prior sigma so lattice jumps are
% %possible.
% drift_prior_sigma = max(drift_prior_sigma,0.003);

min_trial_dur = 0.75;


%%
eps = -1e3;

stim_fs = 100; %in Hz
dt = 0.01;
full_nPix = 36;
if sim_gain == 8
    full_nPix = full_nPix + 18;
end

% stim_params = NIMcreate_stim_params([flen full_nPix],dt);
Fr = 1;

beg_buffer = 0.2;
beg_buffer_new = 0.26; %simulation based on inferred eye positions, so there is no inferred eye position for the first 50ms of each trial
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
max_shift = round(16*spatial_usfac*sim_gain);
dshift = 1;
shifts = -max_shift:dshift:max_shift;
n_shifts = length(shifts);
zero_frame = find(shifts==0);
temp_dx = [-2*max_shift:dshift:2*max_shift];
shift_dx_edges = [temp_dx*sp_dx-dshift*sp_dx/2 temp_dx(end)*sp_dx+dshift*sp_dx/2];
ovshift_dx_edges = [shifts*sp_dx-dshift*sp_dx/2 shifts(end)*sp_dx+dshift*sp_dx/2];

max_Dshift = round(10*spatial_usfac*sim_gain);
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
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
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
% cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori & expt_sac_dir == bar_ori);
cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori);


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
            
                cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            if ~isempty(all_stim_times)
                if any(cur_stim_times < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times];
            all_t_axis = [all_t_axis; cur_t_axis];
            %             all_Xmat = [all_Xmat; bar_Xmat];
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
cd(anal_dir)
load(sim_data_name,'extra_stim');
cd(data_dir);

if sim_gain > 1
    full_nPix = full_nPix + sim_gain*8;
end
stim_params = NIMcreate_stim_params([flen full_nPix],dt);
% if sim_gain > 1
%     extra_pix = full_nPix - 36;
% 
%     extra_stim = zeros(length(all_t_axis),extra_pix);
%     rand_mat = rand(length(all_t_axis),extra_pix);
%     extra_stim(rand_mat <= 0.06) = -1;
%     extra_stim(rand_mat >= 0.94) = 1;
% else
%     extra_stim = [];
% end
all_stim_mat = cat(2,all_stim_mat,extra_stim);

%%
full_nPix_us = spatial_usfac*full_nPix;
if spatial_usfac > 1
    all_stimmat_up = zeros(size(all_stim_mat,1),full_nPix_us);
    for ii = 1:size(all_stim_mat,2)
        for jj = 1:spatial_usfac
        all_stimmat_up(:,spatial_usfac*(ii-1)+jj) = all_stim_mat(:,ii);
        end
    end
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
use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1) & Xinds_up(:) <= cur_use_pix(end));
cnt = 0;
while length(use_kInds_up) < use_nPix_us*flen
    use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1)-cnt/spatial_usfac & Xinds_up(:) <= cur_use_pix(end) + cnt/spatial_usfac);
    if length(use_kInds_up) < use_nPix_us*flen
        use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1)-(cnt+1)/spatial_usfac & Xinds_up(:) <= cur_use_pix(end) + cnt/spatial_usfac);
    end
    cnt = cnt + 1;
end

%% DEFINE DATA USED FOR ANALYSIS
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer);
used_inds_new = find(all_tsince_start >= beg_buffer_new & (trial_dur-all_tsince_start) >= end_buffer);

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

sim_uinds = find(ismember(used_inds,used_inds_new));
used_inds = used_inds(sim_uinds);

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

%% DEFINE POINTS TO ALLOW RAPID CHANGES (TRIAL STARTS AFTER SACCADES)
n_trials = length(trial_start_inds);
trial_ids = nan(NT,1);
for ii = 1:n_trials
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    trial_ids(cur_inds) = ii;
end

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

% fix_prior = diff(erf(ovshift_dx_edges/(fix_prior_sigma*sqrt(2))));
% fix_prior = log(fix_prior/sum(fix_prior));
% fix_prior(fix_prior < eps) = eps;

fix_prior = -(shifts*sp_dx).^2./(2*fix_prior_sigma^2);
fix_prior = fix_prior - logsumexp(fix_prior);

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


%% INCORPORATE INFERRED EYE-POSITIONS AND MAKE SIMULATED ERRORS
fprintf('Loading simulated data\n');
cd(anal_dir)

load(sim_data_name);
tr_set = sim_tr_set;
n_tr_chs = length(tr_set);

Robs_mat = Robs_mat(sim_uinds,:);
Robs_xv = Robs_xv(sim_uinds,:);
sim_prate_mat = sim_prate_mat(sim_uinds,:);

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
cd(anal_dir);
if ~exist(['./' mod_data_name '.mat'],'file') || recompute_init_mods == 1
    
    for ss = 1:n_tr_chs;
        fprintf('Computing base LLs for U %d of %d\n',ss,n_tr_chs);
        cur_tr_inds = tr_inds(~isnan(Robs_mat(tr_inds,ss)));
        tr_NT = length(cur_tr_inds);
        Robs = Robs_mat(cur_tr_inds,ss);
        
        tr_X{1} = all_Xmat(used_inds(cur_tr_inds),use_kInds);
        tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
        
        if use_sac_kerns
            tr_X{3} = Xsac(cur_tr_inds,:);
            tr_X{4} = Xmsac(cur_tr_inds,:);
            null_mod = NMMinitialize_model(null_stim_params,[1 1 1],repmat({'lin'},1,3),null_reg_params,[1 2 3]);
        else
            null_mod = NMMinitialize_model(null_stim_params,[1],{'lin'},null_reg_params,[1]);
        end
        
        null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));
        
        gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
        gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent,init_optim_p);
        
        tr_X{1} = all_Xmat_us(used_inds(cur_tr_inds),use_kInds_up);
        %spatial up-sampling of filter estimates
        base_filts = reshape([gqm1.mods(find(init_Xtargs == 1)).filtK],[flen use_nPix n_squared_filts+1]);
        if spatial_usfac > 1
            base_filts_up = zeros(flen,use_nPix_us,n_squared_filts+1);
            for ii = 1:use_nPix
                for jj = 1:spatial_usfac
                    base_filts_up(:,spatial_usfac*(ii-1)+jj,:) = 0.5*base_filts(:,ii,:);
                end
            end
        elseif spatial_usfac == 1
            base_filts_up = base_filts;
        else
            error('unsupported')
        end
        base_filts_up = reshape(base_filts_up,use_nPix_us*flen,n_squared_filts+1);
        
        init_filts{end} = gqm1.mods(find(init_Xtargs==2)).filtK;
        for ii = 1:n_squared_filts+1
            init_filts{ii} = base_filts_up(:,ii);
        end
        gqm2 = NMMinitialize_model(fin_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs,init_filts);
        gqm2.spk_NL_params(1) = gqm1.spk_NL_params(1);
        [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm2,Robs,tr_X);
        gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
        gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
        if use_sac_kerns
            gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,3,zeros(n_sac_bins,1),sac_reg_params); %sac term
            gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,4,zeros(n_sac_bins,1),sac_reg_params); %nsac term
        end
        gqm2 = NMMfit_filters(gqm2,Robs,tr_X,[],[],silent);
        
        null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));
        all_nullmod(ss) = null_mod;
        
        all_mod_fits(ss) = gqm2;
        all_mod_fits_withspkNL(ss) = NMMfit_logexp_spkNL(gqm2,Robs,tr_X);
        
        [LL,~,pred_rate] = NMMmodel_eval(all_mod_fits_withspkNL(ss),Robs,tr_X);
        [null_LL(ss),~,null_prate] = NMMmodel_eval(null_mod,Robs,tr_X(2:end));
        [all_mod_R2(ss),all_mod_dev(ss),all_null_dev(ss)] = pseudo_r2(Robs,pred_rate,null_prate);
        all_mod_LLimp(ss) = (LL-null_LL(ss))/log(2);
        
        xv_nullmod(ss) = NMMfit_filters(all_nullmod(ss),Robs_xv(cur_tr_inds,ss),tr_X(2:end));
        [null_xvLL(ss),~,null_prate] = NMMmodel_eval(xv_nullmod(ss),Robs_mat(cur_tr_inds,ss),tr_X(2:end));
        
        [newxvLL,~,new_prate] = NMMmodel_eval(all_mod_fits_withspkNL(ss),Robs_xv(cur_tr_inds,ss),tr_X);
        all_mod_xvLLimp(ss) = (newxvLL - null_xvLL(ss))/log(2);
        all_mod_xvR2(ss) = pseudo_r2(Robs_xv(cur_tr_inds,ss),new_prate,null_prate);
        
    end
    save(mod_data_name,'all_mod*','all_nullmod','null_LL','*_trials','xv_nullmod','null_xvLL');
    
else
    fprintf('Loading pre-computed initial models\n');
    load(mod_data_name);
end

%% ITERATE FIXATION-BASED CORRECTIONS

it_mods{1} = all_mod_fits;
it_mods_spkNL{1} = all_mod_fits_withspkNL;
it_LLimp(1,:) = all_mod_LLimp;
it_R2(1,:) = all_mod_R2;
it_dev(1,:) = all_mod_dev;
it_xvR2(1,:) = all_mod_xvR2;
it_xvLLimp(1,:) = all_mod_xvLLimp;
it_fix_sigma(1) = fix_prior_sigma;
for nn = 1:n_fix_inf_it
    fprintf('Inferring fixation corrections, iter %d of %d\n',nn,n_fix_inf_it);
    
    %% PREPROCESS MODEL COMPONENTS
    filt_bank = zeros(n_tr_chs,klen_us,n_squared_filts+1);
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
    block_out = Xblock(used_inds,:)*lin_kerns';
    if use_sac_kerns
        sac_out = Xsac*sac_kerns';
        msac_out = Xmsac*msac_kerns';
    end
    
    %% ESTIMATE LL for each shift in each stimulus frame
    cur_Xmat = all_Xmat_us(used_inds,:);
    %precompute LL at all shifts for all units
    frame_LLs = nan(NT,n_shifts);
    for xx = 1:length(shifts)
        fprintf('Shift %d of %d\n',xx,n_shifts);
        cur_stim_shift = cur_Xmat*shift_mat{xx};
        
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
    
    %% INFER MICRO-SAC SEQUENCE
    fix_LLs = nan(n_fixs,n_shifts);
    for ii = 1:n_fixs
        cur_inds = pfix_start_inds(ii):pfix_stop_inds(ii);
        fix_LLs(ii,:) = sum(frame_LLs(cur_inds,:));
    end
    
    if all(use_coils==0)
        lgamma = bsxfun(@plus,fix_LLs,fix_prior);
        lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    else
        lalpha=zeros(n_fixs,n_shifts);
        lbeta = zeros(n_fixs,n_shifts);
        %compute forward messages
        lalpha(1,:) = fix_prior + fix_LLs(1,:);
        for t=2:n_fixs
            gprobs = diff(erf((shift_dx_edges + measured_fix_deltas(t))/(fix_noise_sigma*sqrt(2))));
            cur_lA = log(toeplitz(gprobs(n_shifts:end)',gprobs(n_shifts:-1:1)'));
            cur_lA(isinf(cur_lA)) = eps;
            cur_lA = bsxfun(@plus,cur_lA,fix_prior);
            cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2));
            lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + fix_LLs(t,:);
        end
        
        %compute backward messages
        lbeta(n_fixs,:)=log(ones(1,n_shifts));
        for t=n_fixs-1:-1:1
            gprobs = diff(erf((shift_dx_edges + measured_fix_deltas(t+1))/(fix_noise_sigma*sqrt(2))));
            cur_lA = log(toeplitz(gprobs(n_shifts:end)',gprobs(n_shifts:-1:1)'));
            cur_lA(isinf(cur_lA)) = eps;
            cur_lA = bsxfun(@plus,cur_lA,fix_prior);
            cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2));
            
            lf1 = lbeta(t+1,:) + fix_LLs(t+1,:);
            lbeta(t,:) = logmulexp(lf1,cur_lA');
        end
        lgamma= lalpha + lbeta;
        lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    end
    
    gamma = exp(lgamma);
    
    it_fix_sigma(nn) = sqrt(mean(sum(bsxfun(@times,gamma,shifts*sp_dx).^2,2)));
    
    it_fix_post_mean(nn,:) = sum(bsxfun(@times,gamma,shifts),2);
    cur_diff = bsxfun(@minus,it_fix_post_mean(nn,:)',shifts).^2;
    it_fix_post_std(nn,:) = sqrt(sum(cur_diff.*gamma,2));
    
    %back-project saccade-times
    all_fix_post_mean_cor = nan(NT,1);
    all_fix_post_mean_cor(~isnan(fix_ids)) = it_fix_post_mean(nn,fix_ids(~isnan(fix_ids)));
    all_fix_post_mean_cor = interp1(find(~isnan(fix_ids)),all_fix_post_mean_cor(~isnan(fix_ids)),1:NT);
    all_fix_post_mean_cor(isnan(all_fix_post_mean_cor)) = 0;
    %%
    all_fix_post_mean_cor = round(all_fix_post_mean_cor) + max_Tshift + 1;
    all_shift_stimmat_up = all_stimmat_up;
    for i=1:NT
        all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_fix_post_mean_cor(i)};
    end
    all_Xmat_up_fixcor = create_time_embedding(all_shift_stimmat_up,stim_params_us);
    all_Xmat_up_fixcor = all_Xmat_up_fixcor(used_inds,:);
    
    %% REFIT ALL CELLS
    cur_X{1} = all_Xmat_up_fixcor(:,use_kInds_up);
%     cur_X{1} = all_Xmat_us(used_inds,use_kInds_up);
    cur_X{2} = Xblock(used_inds,:);
    if use_sac_kerns
        cur_X{3} = Xsac;
        cur_X{4} = Xmsac;
    end
    
    silent = 1;
    for ss = 1:length(tr_set)
        fprintf('Refitting model for tr cell %d of %d\n',ss,length(tr_set));
        cur_tr_uset = tr_inds(~isnan(Robs_mat(tr_inds,ss)));
        
        tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
        
        it_mods{nn+1}(ss) = it_mods{nn}(ss);
        it_mods{nn+1}(ss) = NMMfit_filters(it_mods{nn+1}(ss),Robs_mat(cur_tr_uset,ss),...
            tr_X,[],[],silent); %fit stimulus filters
        
        %refit spk NL
        it_mods_spkNL{nn+1}(ss) = NMMfit_logexp_spkNL(it_mods{nn+1}(ss),Robs_mat(cur_tr_uset,ss),tr_X);
        
        [newLL,~,new_prate] = NMMmodel_eval(it_mods_spkNL{nn+1}(ss),Robs_mat(cur_tr_uset,ss),tr_X);
        [~,~,null_prate] = NMMmodel_eval(all_nullmod(ss),Robs_mat(cur_tr_uset,ss),tr_X(2:end));
        [it_R2(nn+1,ss),it_dev(nn+1,ss)] = pseudo_r2(Robs_mat(cur_tr_uset,ss),new_prate,null_prate);
        
        it_LLimp(nn+1,ss) = (newLL - null_LL(ss))/log(2);
        if nn > 1
            fprintf('Original: %.4f  Prev: %.4f  New: %.4f\n',all_mod_LLimp(ss),it_LLimp(nn,ss),it_LLimp(nn+1,ss));
        else
            fprintf('Original: %.4f  New: %.4f\n',all_mod_LLimp(ss),it_LLimp(nn+1,ss));
        end
                
        [newxvLL,~,new_prate] = NMMmodel_eval(it_mods_spkNL{nn+1}(ss),Robs_xv(cur_tr_uset,ss),tr_X);
        it_xvLLimp(nn+1,ss) = (newxvLL - null_xvLL(ss))/log(2);
        [it_xvR2(nn+1,ss)] = pseudo_r2(Robs_xv(cur_tr_uset,ss),new_prate,null_prate);
        
    end
end

%% NOW INFER DRIFT CORRECTIONS
if n_fix_inf_it == 0
    all_Xmat_up_fixcor = all_Xmat_us(used_inds,:);
    it_fix_post_mean = zeros(1,n_fixs);
end
dit_mods{1} = it_mods{n_fix_inf_it+1};
dit_mods_spkNL{1} = it_mods_spkNL{n_fix_inf_it+1};
dit_LLimp(1,:) = it_LLimp(n_fix_inf_it+1,:);
dit_R2(1,:) = it_R2(n_fix_inf_it+1,:);
dit_dev(1,:) = it_dev(n_fix_inf_it+1,:);
dit_xvR2(1,:) = it_xvR2(n_fix_inf_it+1,:);
dit_xvLLimp(1,:) = it_xvLLimp(n_fix_inf_it+1,:);

% load('monoc_eyecorr_hbar2.mat','dit_mods','dit_mods_spkNL');
% true_mods = dit_mods{end}(tr_set);
% true_mods_spkNL = dit_mods_spkNL{end}(tr_set);
% dit_mods{1} = true_mods;
% dit_mods_spkNL{1} = true_mods_spkNL;
for nn = 1:n_drift_inf_it
    fprintf('Inferring drift corrections, iter %d of %d\n',nn,n_drift_inf_it);
    %% PREPROCESS MODEL COMPONENTS
    filt_bank = zeros(n_tr_chs,klen_us,n_squared_filts+1);
    lin_kerns = nan(n_tr_chs,n_blocks);
    if use_sac_kerns
        sac_kerns = nan(n_tr_chs,n_sac_bins);
        msac_kerns = nan(n_tr_chs,n_sac_bins);
    end
    mod_spkNL_params = nan(n_tr_chs,3);
    for ss = 1:n_tr_chs
        cur_Xtargs = [dit_mods{nn}(ss).mods(:).Xtarget];
        cur_k = [dit_mods{nn}(ss).mods(cur_Xtargs == 1).filtK];
        n_used_filts = size(cur_k,2);
        filt_bank(ss,:,1:n_used_filts) = cur_k;
        mod_spkNL_params(ss,:) = dit_mods_spkNL{nn}(ss).spk_NL_params;
        lin_kerns(ss,:) = dit_mods{nn}(ss).mods(cur_Xtargs == 2).filtK;
        if use_sac_kerns
            sac_kerns(ss,:) = dit_mods{nn}(ss).mods(cur_Xtargs == 3).filtK;
            msac_kerns(ss,:) = dit_mods{nn}(ss).mods(cur_Xtargs == 4).filtK;
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
                    gprobs = diff(erf((Dshift_dx_edges+cur_drift_mean(iii))/(post_drift_sigma*sqrt(2)*drift_dsf)));
                    cur_lA(:,:,iii) = toeplitz(gprobs(n_Dshifts:end),gprobs(n_Dshifts:-1:1));
                end
                cur_lA = log(bsxfun(@rdivide,cur_lA,sum(cur_lA,2)));
                cur_lA(cur_lA < eps) = eps;
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
    end
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    gamma = exp(lgamma);
    drift_post_mean(nn,:) = sum(bsxfun(@times,gamma,Dshifts),2);
    drift_post_std(nn,:) = sqrt(sum(bsxfun(@times,gamma,Dshifts.^2),2) - drift_post_mean(nn,:)'.^2);
    
    %     drift_post_mean(nn,:) = interp1(find(~isnan(fix_ids)),drift_post_mean(nn,~isnan(fix_ids)),1:NT);
    %     drift_post_std(nn,:) = interp1(find(~isnan(fix_ids)),drift_post_std(nn,~isnan(fix_ids)),1:NT);
    drift_post_mean(nn,:) = interp1(find(~isnan(pfix_ids)),drift_post_mean(nn,~isnan(pfix_ids)),1:NT);
    drift_post_std(nn,:) = interp1(find(~isnan(pfix_ids)),drift_post_std(nn,~isnan(pfix_ids)),1:NT);
    drift_post_mean(nn,isnan(drift_post_mean(nn,:))) = 0;
    
    %% construct drift-corrected X-mat
    fix_post_cor = nan(NT,1);
    fix_post_cor(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
    fix_post_cor = interp1(find(~isnan(fix_ids)),fix_post_cor(~isnan(fix_ids)),1:NT);
    fix_post_cor(isnan(fix_post_cor)) = 0;
    drift_post_cor = squeeze(drift_post_mean(nn,:));
    
    for ii = 1:length(trial_start_inds)
        cur_inds = trial_start_inds(ii):trial_end_inds(ii);
        drift_post_cor(cur_inds(1:end-sac_shift)) = drift_post_cor(cur_inds(sac_shift+1:end));
    end
    drift_post_cor = interp1(find(~isnan(fix_ids)),drift_post_cor(~isnan(fix_ids)),1:NT);
    
    %     %back project drift (within-fixation) by sac_shift
    %     for ii = 1:n_fixs
    %         cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    %         if length(cur_inds) > sac_shift
    %             drift_post_cor(cur_inds(1:end-sac_shift+1)) = drift_post_cor(cur_inds(sac_shift:end));
    %         end
    %     end
    
    drift_post_cor(isnan(drift_post_cor)) = 0;
    
    all_post_cor = round((fix_post_cor+drift_post_cor)) + max_Tshift + 1;
    
    %RECOMPUTE XMAT
    all_shift_stimmat_up = all_stimmat_up;
    for i=1:NT
        all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_post_cor(i)};
    end
    cur_X{1} = create_time_embedding(all_shift_stimmat_up,stim_params_us);
    cur_X{1} = cur_X{1}(used_inds,use_kInds_up);
    cur_X{2} = Xblock(used_inds,:);
    if use_sac_kerns
        cur_X{3} = Xsac;
        cur_X{4} = Xmsac;
    end
    
    %% REFIT ALL CELLS
    silent = 1;
    for ss = 1:length(tr_set)
        fprintf('Refitting model for tr cell %d of %d\n',ss,length(tr_set));
        cur_tr_uset = tr_inds(~isnan(Robs_mat(tr_inds,ss)));
        
        tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
        
        dit_mods{nn+1}(ss) = dit_mods{nn}(ss);
        dit_mods{nn+1}(ss) = NMMfit_filters(dit_mods{nn+1}(ss),Robs_mat(cur_tr_uset,ss),...
            tr_X,[],[],silent); %fit stimulus filters
        
        dit_mods_spkNL{nn+1}(ss) = NMMfit_logexp_spkNL(dit_mods{nn+1}(ss),Robs_mat(cur_tr_uset,ss),tr_X);
        
        [newLL,~,new_prate] = NMMmodel_eval(dit_mods_spkNL{nn+1}(ss),Robs_mat(cur_tr_uset,ss),tr_X);
        [~,~,null_prate] = NMMmodel_eval(all_nullmod(ss),Robs_mat(cur_tr_uset,ss),tr_X(2:end));
        [dit_R2(nn+1,ss),dit_dev(nn+1,ss)] = pseudo_r2(Robs_mat(cur_tr_uset,ss),new_prate,null_prate);
        
        dit_LLimp(nn+1,ss) = (newLL - null_LL(ss))/log(2);
        
        fprintf('Original: %.5f  Prev: %.5f  New: %.5f\n',all_mod_LLimp(ss),dit_LLimp(nn,ss),dit_LLimp(nn+1,ss));
        
        [newxvLL,~,new_prate] = NMMmodel_eval(dit_mods_spkNL{nn+1}(ss),Robs_xv(cur_tr_uset,ss),tr_X);
        dit_xvLLimp(nn+1,ss) = (newxvLL - null_xvLL(ss))/log(2);
        [dit_xvR2(nn+1,ss)] = pseudo_r2(Robs_xv(cur_tr_uset,ss),new_prate,null_prate);

    end
end

%%
cd(anal_dir);
save(data_name,'sim_gain','dit_*','it_*','all_mod_*','tr_set','drift_*','it_fix_*');

%%
if do_highres_est
     old_usfac = spatial_usfac;
   spatial_usfac = hr_spatial_usfac;
    new_usratio = spatial_usfac/old_usfac;
    
        full_nPix_us = spatial_usfac*full_nPix;

        drift_dsf = 2;

    use_nPix_us = use_nPix*spatial_usfac;
    klen_us = use_nPix_us*flen;
    sp_dx = 0.0565/spatial_usfac;
    max_shift = round(15*spatial_usfac*sim_gain);
    dshift = 1;
    shifts = -max_shift:dshift:max_shift;
    n_shifts = length(shifts);
    zero_frame = find(shifts==0);
    temp_dx = [-2*max_shift:dshift:2*max_shift];
    shift_dx_edges = [temp_dx*sp_dx-dshift*sp_dx/2 temp_dx(end)*sp_dx+dshift*sp_dx/2];
    ovshift_dx_edges = [shifts*sp_dx-dshift*sp_dx/2 shifts(end)*sp_dx+dshift*sp_dx/2];
    
    max_Dshift = round(8*spatial_usfac*sim_gain);
    Dshifts = -max_Dshift:dshift:max_Dshift;
    n_Dshifts = length(Dshifts);
    zero_Dframe = find(Dshifts==0);
    temp_dx = [-2*max_Dshift:dshift:2*max_Dshift];
    Dshift_dx_edges = [temp_dx*sp_dx-dshift*sp_dx/2 temp_dx(end)*sp_dx+dshift*sp_dx/2];
    ovDshift_dx_edges = [Dshifts*sp_dx-dshift*sp_dx/2 Dshifts(end)*sp_dx+dshift*sp_dx/2];
    
    max_Tshift = max_shift + max_Dshift;
    Tshifts = -max_Tshift:dshift:max_Tshift;
    
    
    %%
    all_stimmat_up = zeros(size(all_stim_mat,1),full_nPix_us);
    for ii = 1:size(all_stim_mat,2)
        for jj = 1:spatial_usfac
            all_stimmat_up(:,spatial_usfac*(ii-1) + jj) = all_stim_mat(:,ii);
        end
    end
    stim_params_us = NMMcreate_stim_params([flen full_nPix_us],dt);
    all_Xmat_us = create_time_embedding(all_stimmat_up,stim_params_us);
    
    [Xinds_up,Tinds] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
    use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1)-1/spatial_usfac & Xinds_up(:) <= cur_use_pix(end) + 1/spatial_usfac);
    if length(use_kInds_up) < use_nPix_us*flen
        use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1)-2/spatial_usfac & Xinds_up(:) <= cur_use_pix(end) + 1/spatial_usfac);
    end
    cnt = 2;
    while length(use_kInds_up) < use_nPix_us*flen
        use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1)-cnt/spatial_usfac & Xinds_up(:) <= cur_use_pix(end) + cnt/spatial_usfac);
        if length(use_kInds_up) < use_nPix_us*flen
            use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1)-(cnt+1)/spatial_usfac & Xinds_up(:) <= cur_use_pix(end) + cnt/spatial_usfac);
        end
        cnt = cnt + 1;
    end
    %%
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
    
    drift_jump_prior = -(Dshifts*sp_dx).^2./(2*drift_jump_sigma^2);
    drift_jump_prior = drift_jump_prior - logsumexp(drift_jump_prior);
    
    cdist = squareform(pdist(Dshifts'*sp_dx));
    base_lA = -cdist.^2/(2*(post_drift_sigma*drift_dsf)^2);
    base_lA = bsxfun(@minus,base_lA,logsumexp(base_lA,2)); %normalize
    
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

    %% construct drift-corrected X-mat
    old_best_mods = dit_mods{end};
    best_fix_cor = it_fix_post_mean(end,:)*spatial_usfac/old_usfac;
    best_fix_std = it_fix_post_std(end,:)*spatial_usfac/old_usfac;
    best_drift_cor = drift_post_mean(end,:)*spatial_usfac/old_usfac;
    
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
    
    drift_post_cor(isnan(drift_post_cor)) = 0;
    
    all_post_cor = round((fix_post_cor+drift_post_cor)) + max_Tshift + 1;
    
    %RECOMPUTE XMAT
    all_shift_stimmat_up = all_stimmat_up;
    for i=1:NT
        all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_post_cor(i)};
    end
    all_Xmat_cor = create_time_embedding(all_shift_stimmat_up,stim_params_us);
    
    %%
    fin_stim_params(1) = NMMcreate_stim_params([flen use_nPix_us],dt);
    
    for ss = 1:n_tr_chs
        fprintf('Computing base LLs for Unit %d of %d\n',ss,n_tr_chs);
        cur_tr_inds = tr_inds(~isnan(Robs_mat(tr_inds,ss)));
        tr_NT = length(cur_tr_inds);
        Robs = Robs_mat(cur_tr_inds,ss);
        null_mod = all_nullmod(ss);
        if ~isempty(cur_tr_inds) && nansum(Robs) > 0
            
            tr_X{1} = all_Xmat_cor(used_inds(cur_tr_inds),use_kInds_up);
            tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
            
            if use_sac_kerns
                tr_X{3} = Xsac(cur_tr_inds,:);
                tr_X{4} = Xmsac(cur_tr_inds,:);
            end
            
            %spatial up-sampling of filter estimates
            base_filts = reshape([old_best_mods(ss).mods(find(init_Xtargs == 1)).filtK],[flen use_nPix*old_usfac n_squared_filts+1]);
            base_filts_up = zeros(flen,use_nPix_us,n_squared_filts+1);
            for ii = 1:use_nPix*2
                for jj = 1:new_usratio
                    base_filts_up(:,new_usratio*(ii-1)+jj,:) = 1/new_usratio*base_filts(:,ii,:);
                end
            end
            base_filts_up = reshape(base_filts_up,use_nPix_us*flen,n_squared_filts+1);
            
            init_filts{end} = old_best_mods(ss).mods(find(init_Xtargs==2)).filtK;
            for ii = 1:n_squared_filts+1
                init_filts{ii} = base_filts_up(:,ii);
            end
            gqm2 = NMMinitialize_model(fin_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs,init_filts);
            [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm2,Robs,tr_X);
            vgint = var(gint)'; sgint = std(gint)';
            vgint(vgint == 0) = 1; sgint(sgint == 0) = 1;
            gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./vgint);
            gqm2 = NMMadjust_regularization(gqm2,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./sgint);
            if use_sac_kerns
                gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,3,zeros(n_sac_bins,1),sac_reg_params); %sac term
                gqm2 = NMMadd_NLinput(gqm2,{'lin'},1,4,zeros(n_sac_bins,1),sac_reg_params); %nsac term
            end
            HR_mod_fits(ss) = gqm2;
            HR_mod_fits(ss).spk_NL_params(1) = old_best_mods(ss).spk_NL_params(1);
            HR_mod_fits(ss) = NMMfit_filters(HR_mod_fits(ss),Robs,tr_X,[],[],1);
            
            HR_mod_fits_withspkNL(ss) = NMMfit_logexp_spkNL(HR_mod_fits(ss),Robs,tr_X);
            
            [LL,~,pred_rate] = NMMmodel_eval(HR_mod_fits_withspkNL(ss),Robs,tr_X);
            [null_LL(ss),~,null_prate] = NMMmodel_eval(null_mod,Robs,tr_X(2:end));
            HR_mod_LLimp(ss) = (LL-null_LL(ss))/log(2);
            fprintf('Prev %.4f  New: %.4f\n',dit_LLimp(end,ss),HR_mod_LLimp(ss));
        else
            HR_mod_LLimp(ss) = nan;
        end
        
    end
    
    clear all_Xmat_cor
    
    
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
    n_HR_inf_it = 1;
    HRit_mods{1} = HR_mod_fits;
    HRit_mods_spkNL{1} = HR_mod_fits_withspkNL;
    HRit_LLimp(1,:) = HR_mod_LLimp;
    for nn = 1:n_HR_inf_it
        fprintf('Inferring drift corrections, iter %d of %d\n',nn,n_HR_inf_it);
        
        %% PREPROCESS MODEL COMPONENTS
        filt_bank = zeros(n_tr_chs,klen_us,n_squared_filts+1);
        lin_kerns = nan(n_tr_chs,n_blocks);
        if use_sac_kerns
            sac_kerns = nan(n_tr_chs,n_sac_bins);
            msac_kerns = nan(n_tr_chs,n_sac_bins);
        end
        mod_spkNL_params = nan(n_tr_chs,3);
        for ss = 1:n_tr_chs
            cur_Xtargs = [HRit_mods{nn}(ss).mods(:).Xtarget];
            cur_k = [HRit_mods{nn}(ss).mods(cur_Xtargs == 1).filtK];
            n_used_filts = size(cur_k,2);
            filt_bank(ss,:,1:n_used_filts) = cur_k;
            mod_spkNL_params(ss,:) = HRit_mods_spkNL{nn}(ss).spk_NL_params;
            lin_kerns(ss,:) = HRit_mods{nn}(ss).mods(cur_Xtargs == 2).filtK;
            if use_sac_kerns
                sac_kerns(ss,:) = HRit_mods{nn}(ss).mods(cur_Xtargs == 3).filtK;
                msac_kerns(ss,:) = HRit_mods{nn}(ss).mods(cur_Xtargs == 4).filtK;
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
                        gprobs = diff(erf((Dshift_dx_edges+cur_drift_mean(iii))/(post_drift_sigma*sqrt(2)*drift_dsf)));
                        cur_lA(:,:,iii) = toeplitz(gprobs(n_Dshifts:end),gprobs(n_Dshifts:-1:1));
                    end
                    cur_lA = log(bsxfun(@rdivide,cur_lA,sum(cur_lA,2)));
                    cur_lA(cur_lA < eps) = eps;
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
        end
        lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
        gamma = exp(lgamma);
        HRdrift_post_mean(nn,:) = sum(bsxfun(@times,gamma,Dshifts),2);
        HRdrift_post_std(nn,:) = sqrt(sum(bsxfun(@times,gamma,Dshifts.^2),2) - HRdrift_post_mean(nn,:)'.^2);
        
        HRdrift_post_mean(nn,:) = interp1(find(~isnan(pfix_ids)),HRdrift_post_mean(nn,~isnan(pfix_ids)),1:NT);
        HRdrift_post_std(nn,:) = interp1(find(~isnan(pfix_ids)),HRdrift_post_std(nn,~isnan(pfix_ids)),1:NT);
        HRdrift_post_mean(nn,isnan(HRdrift_post_mean(nn,:))) = 0;
        
        %% construct drift-corrected X-mat
%         fix_post_cor = nan(NT,1);
%         fix_post_cor(~isnan(fix_ids)) = best_fix_cor(fix_ids(~isnan(fix_ids)));
%         fix_post_cor = interp1(find(~isnan(fix_ids)),fix_post_cor(~isnan(fix_ids)),1:NT);
%         fix_post_cor(isnan(fix_post_cor)) = 0;
%         drift_post_cor = squeeze(HRdrift_post_mean(nn,:));
%         
%         for ii = 1:length(trial_start_inds)
%             cur_inds = trial_start_inds(ii):trial_end_inds(ii);
%             drift_post_cor(cur_inds(1:end-sac_shift)) = drift_post_cor(cur_inds(sac_shift+1:end));
%         end
%         drift_post_cor = interp1(find(~isnan(fix_ids)),drift_post_cor(~isnan(fix_ids)),1:NT);
%         
%         drift_post_cor(isnan(drift_post_cor)) = 0;
%         
%         all_post_cor = round((fix_post_cor+drift_post_cor)) + max_Tshift + 1;
%         
%         %RECOMPUTE XMAT
%         for i=1:NT
%             all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_post_cor(i)};
%         end
%         cur_X{1} = create_time_embedding(all_shift_stimmat_up,stim_params_us);
%         cur_X{1} = cur_X{1}(used_inds,use_kInds_up);
%         
%         cur_X{2} = Xblock(used_inds,:);
%         if use_sac_kerns
%             cur_X{3} = Xsac;
%             cur_X{4} = Xmsac;
%         end
        
        %% REFIT ALL CELLS
%         silent = 1;
%         for ss = 1:length(tr_set)
%             cur_cell = ss;
%             fprintf('Refitting model for tr cell %d of %d\n',ss,length(tr_set));
%             cur_tr_uset = tr_inds(~isnan(Robs_mat(tr_inds,ss)));
%             
%             tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
%             
%             HRit_mods{nn+1}(cur_cell) = HRit_mods{nn}(cur_cell);
%             HRit_mods{nn+1}(cur_cell) = NMMfit_filters(HRit_mods{nn+1}(cur_cell),Robs_mat(cur_tr_uset,ss),...
%                 tr_X,[],[],silent); %fit stimulus filters
%             
%             HRit_mods_spkNL{nn+1}(cur_cell) = NMMfit_logexp_spkNL(HRit_mods{nn+1}(cur_cell),Robs_mat(cur_tr_uset,ss),tr_X);
%             
%             [newLL,~,new_prate] = NMMmodel_eval(HRit_mods_spkNL{nn+1}(cur_cell),Robs_mat(cur_tr_uset,ss),tr_X);
%             [~,~,null_prate] = NMMmodel_eval(all_nullmod(cur_cell),Robs_mat(cur_tr_uset,ss),tr_X(2:end));
%             
%             HRit_LLimp(nn+1,cur_cell) = (newLL - null_LL(cur_cell))/log(2);
%             
%             fprintf('Original: %.5f  Prev: %.5f  New: %.5f\n',all_mod_LLimp(cur_cell),HRit_LLimp(nn,cur_cell),HRit_LLimp(nn+1,cur_cell));
%             
%             [newxvLL,~,new_prate] = NMMmodel_eval(HRit_mods_spkNL{nn+1}(cur_cell),Robs_xv(cur_tr_uset,ss),tr_X);
%             HRit_xvLLimp(nn+1,cur_cell) = (newxvLL - null_xvLL(cur_cell))/log(2);
%             [HRit_xvR2(nn+1,cur_cell)] = pseudo_r2(Robs_xv(cur_tr_uset,ss),new_prate,null_prate);
% 
%         end
    end
    
    %%
    cd(anal_dir)
    save(data_name_hres,'HRit*','HRdrift*','HR_mod_*','sp_dx','n_*inf*');
end

%%
    %%
    fin_fix_corr = nan(NT,1);
    fin_fix_std = nan(NT,1);
    fin_fix_corr(~isnan(fix_ids)) = best_fix_cor(fix_ids(~isnan(fix_ids)));
    fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
    fin_fix_std(~isnan(fix_ids)) = best_fix_std(fix_ids(~isnan(fix_ids)));
    fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);
    
    fin_fix_corr = fin_fix_corr*sp_dx;
    fin_fix_std = fin_fix_std*sp_dx;
    
    fin_drift_corr = HRdrift_post_mean(end,:)*sp_dx;
    fin_drift_std = HRdrift_post_std(end,:)*sp_dx;
    
    fix_inds = [];
    for ii = 1:length(trial_start_inds)
        cur_inds = trial_start_inds(ii):trial_end_inds(ii);
        fix_inds = [fix_inds cur_inds];
        fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
        fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
    end
    fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT);
    fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);
    
    fin_tot_corr = fin_fix_corr + fin_drift_corr;
    fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2);

    fin_tot_corr = fin_tot_corr - nanmedian(fin_tot_corr - sim_eyepos(used_inds)');

