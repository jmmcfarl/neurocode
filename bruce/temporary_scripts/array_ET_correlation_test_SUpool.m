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

if bar_ori == 0
    true_moddata_name = 'monoc_eyecorr_hbar_mods';
    true_eyedata_name = 'monoc_eyecorr_hbar';
else
    true_moddata_name = 'monoc_eyecorr_vbar_mods';
    true_eyedata_name = 'monoc_eyecorr_vbar';
end

use_sac_kerns = 1;
use_coils = [0 0]; %[L R]

if any(use_coils > 0)
    anal_name = [anal_name '_Cprior'];
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

use_smooth_eyepos = true;
use_med_sub = true;
%%
xv_frac = 0;

flen = 12;
use_nPix = 16;

n_fix_inf_it = 4; %4
n_drift_inf_it = 4; %2

fix_prior_sigma = 0.15;
fix_noise_sigma = 0.1;
drift_noise_sigma = 0.003;
drift_prior_sigma = 0.003; %.004 may be best here
drift_jump_sigma = 0.1; %0.05 start
drift_dsf = 3;

min_trial_dur = 0.75;

spatial_usfac = 2;

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
max_shift = round(16*spatial_usfac);
dshift = 1;
shifts = -max_shift:dshift:max_shift;
n_shifts = length(shifts);
zero_frame = find(shifts==0);
temp_dx = [-2*max_shift:dshift:2*max_shift];
shift_dx_edges = [temp_dx*sp_dx-dshift*sp_dx/2 temp_dx(end)*sp_dx+dshift*sp_dx/2];
ovshift_dx_edges = [shifts*sp_dx-dshift*sp_dx/2 shifts(end)*sp_dx+dshift*sp_dx/2];

max_Dshift = round(10*spatial_usfac);
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
    extra_pix = full_nPix - expt_npix(cur_block);
    
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
stim_params_us = NMMcreate_stim_params([flen full_nPix_us],dt);

%% select submatrix with central pixels
[Xinds,Tinds] = meshgrid(1:full_nPix,1:flen);
buffer_pix = floor((full_nPix - use_nPix)/2);
cur_use_pix = (1:use_nPix) + buffer_pix;
use_kInds = find(ismember(Xinds(:),cur_use_pix));

[Xinds_up,Tinds] = meshgrid(1/spatial_usfac:1/spatial_usfac:full_nPix,1:flen);
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

fix_prior = diff(erf(ovshift_dx_edges/(sqrt(2)*fix_prior_sigma)));
fix_prior = log(fix_prior/sum(fix_prior));
eps = -1e3;
fix_prior(fix_prior < eps) = eps;

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


%% INCORPORATE INFERRED EYE-POSITIONS AND MAKE SIMULATED ERRORS
fprintf('Incorporating inferred eye-positions\n');
cd(anal_dir)

load(true_eyedata_name,'dit_mods*','et_tr_set','drift_*','it_fix_*');
true_mods = dit_mods_spkNL{end};
tr_set = et_tr_set;
n_tr_chs = length(tr_set);
load(true_moddata_name,'all_mod_SU*');
true_mod_SU = all_mod_SU;
true_mod_SUnum = all_mod_SUnum;

clear dit_mods

NT = length(used_inds);
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
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > sac_shift
        fin_drift_corr(cur_inds(1:end-sac_shift+1)) = fin_drift_corr(cur_inds(sac_shift:end));
        fin_drift_std(cur_inds(1:end-sac_shift+1)) = fin_drift_std(cur_inds(sac_shift:end));
    end
end

sim_eyepos = zeros(length(all_t_axis),1);
sim_eyepos(used_inds) = fin_fix_corr + fin_drift_corr;

%maximum initial corrections
max_sim_pos = max_shift*sp_dx;
sim_eyepos(sim_eyepos > max_sim_pos) = max_sim_pos;
sim_eyepos(sim_eyepos < - max_sim_pos) = -max_sim_pos;
sim_eyepos_us_rnd = round(sim_eyepos/sp_dx);


%%
cd ~/Analysis/bruce/ET_final/
load jbe_all_vert_units
comb_mods = [all_su_data(used_sus).after_mod];
comb_weights = [weight_mean(used_sus)];
comb_expts = [all_su_data(used_sus).expt_num];

load jbe_all_hori_units
comb_mods = [comb_mods all_su_data(used_sus).after_mod];
comb_weights = [comb_weights; weight_mean(used_sus)];
comb_expts = [comb_expts all_su_data(used_sus).expt_num];

n_tr_chs = length(comb_mods);

nPix = use_nPix_us;
stim_params = NMMcreate_stim_params([flen nPix],0.01);
n_frames = 1e5;

rand_mat = rand(n_frames,nPix);

dds = 0.67;
rand_stim = zeros(n_frames,nPix);
rand_stim(rand_mat <= dds/2) = -1;
rand_stim(rand_mat >= 1-(dds/2)) = 1;
test_stim_dense = create_time_embedding(rand_stim,stim_params);

dds = 0.12;
rand_stim = zeros(n_frames,nPix);
rand_stim(rand_mat <= dds/2) = -1;
rand_stim(rand_mat >= 1-(dds/2)) = 1;
test_stim_sparse = create_time_embedding(rand_stim,stim_params);

test_prate_out = nan(n_frames,n_tr_chs);
for ss = 1:n_tr_chs
    avg_block = mean(comb_mods(ss).mods(4).filtK);
    comb_mods(ss).spk_NL_params(1) = comb_mods(ss).spk_NL_params(1) + avg_block;
    comb_mods(ss).mods(4:end) = [];
    
    if ismember(comb_expts(ss),[91 92 93])
        [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(comb_mods(ss),[],test_stim_dense);
        dense_out = std(gint);
        
        [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(comb_mods(ss),[],test_stim_sparse);
        sparse_out = std(gint);
        for jj = 1:3
            comb_mods(ss).mods(jj).filtK = comb_mods(ss).mods(jj).filtK*dense_out(jj)/sparse_out(jj);
        end
    end
    [LL, penLL, pred_rate] = NMMmodel_eval(comb_mods(ss),[],test_stim_sparse);
    test_prate_out(:,ss) = pred_rate;
end

stim_corrmat = nancorr(test_prate_out);
%%
all_trial_durs = (trial_end_inds - trial_start_inds)*dt;
poss_trials = find(all_trial_durs == 3.74);
utrial = poss_trials(2);

cur_inds = find(all_trialvec == utrial);

n_rpt_trials = 200;
used_rpt_set = poss_trials(randperm(length(poss_trials)));
used_rpt_set(n_rpt_trials+1:end) = [];
used_rpt_set = sort(used_rpt_set);

for ii = 1:n_rpt_trials
    cur_dur(ii) = sum(all_trialvec(used_inds) == used_rpt_set(ii));
end
bad = find(cur_dur < 375);
used_rpt_set(bad) = [];
n_rpt_trials = length(used_rpt_set);

rpt_inds = find(ismember(all_trialvec,used_rpt_set));

rpt_inds_isused = ismember(rpt_inds,used_inds);
used_rpt_indset = find(ismember(rpt_inds,used_inds));
used_ind_rptset = find(ismember(used_inds,rpt_inds));

used_rpt_inds = used_inds(ismember(used_inds,rpt_inds));
NT = length(used_rpt_inds);

%%
rpt_stim = all_stimmat_up(cur_inds,:);
all_rpt_stimmat = nan(401,n_rpt_trials,size(all_stimmat_up,2));
for ii = 1:n_rpt_trials
   cur_set = find(all_trialvec == used_rpt_set(ii));
   all_rpt_stimmat(1:length(cur_set),ii,:) = rpt_stim(1:length(cur_set),:);
end
all_rpt_stimmat = reshape(all_rpt_stimmat,401*n_rpt_trials,size(all_stimmat_up,2));
bad_t = find(any(isnan(all_rpt_stimmat),2));
all_rpt_stimmat(bad_t,:) = [];

sim_gain = 1;
sim_eyepos = zeros(length(all_t_axis),1);
sim_eyepos(used_inds) = fin_fix_corr + fin_drift_corr;
sim_eyepos = round(sim_eyepos*sim_gain/sp_dx);

all_rpt_eyepos = sim_eyepos(used_ind_rptset);
all_rpt_eyepos(all_rpt_eyepos > max_shift) = max_shift;
all_rpt_eyepos(all_rpt_eyepos < -max_shift) = -max_shift;

%%
all_stimmat_true = all_rpt_stimmat;
for ii = 1:NT
    all_stimmat_true(used_rpt_indset(ii),:) = shift_matrix_Nd(all_stimmat_true(used_rpt_indset(ii),:),-all_rpt_eyepos(ii),2);
end

all_Xmat_true = create_time_embedding(all_stimmat_true,stim_params_us);


%% SIMULATE SPIKING RESPONSES
fprintf('Creating simulated spike trains\n');
% sim_Xblock = Xblock;
% sim_Xblock(:,1) = 1; sim_Xblock(:,2:end) = 0;

tr_X{1} = all_Xmat_true(used_rpt_indset,use_kInds_up);
% tr_X{2} = sim_Xblock(used_rpt_inds,:);
% tr_X{3} = Xsac(used_ind_rptset,:);
% tr_X{4} = Xmsac(used_ind_rptset,:);

% make Robs_mat
prate_mat = nan(NT,n_tr_chs);
for ss = 1:n_tr_chs
%     cur_mod = true_mods(tr_set(ss));
%     cur_mod.mods([cur_mod.mods.Xtarget] > 2) = [];
%     cur_mod = full_mod_fitspkNL(ss);
%     cur_mod = true_mods2(ss);
    cur_mod = comb_mods(ss);
    [~, ~, pred_rate] = NMMmodel_eval(cur_mod,[], tr_X);
    prate_mat(:,ss) = pred_rate;
end
Robs_mat = poissrnd(prate_mat);


%%
prate_tmat = nan(375,n_rpt_trials,n_tr_chs);
stim_tmat = nan(375,n_rpt_trials,72);
Robs_tmat = nan(375,n_rpt_trials,n_tr_chs);
for ii = 1:n_rpt_trials
   cur_set = find(all_trialvec(used_rpt_inds) == used_rpt_set(ii));
   prate_tmat(1:length(cur_set),ii,:) = prate_mat(cur_set,:);
   Robs_tmat(1:length(cur_set),ii,:) = Robs_mat(cur_set,:);
   stim_tmat(1:length(cur_set),ii,:) = all_stimmat_true(used_rpt_indset(cur_set),:);
end

em_var = squeeze(mean(var(prate_tmat,[],2)));
tot_var = var(prate_mat);
em_var_frac = em_var./tot_var';
avg_var = var(squeeze(mean(prate_tmat,2)));



%%
dsfac = 100;
temp = Robs_tmat;
temp = reshape(temp(1:floor(375/dsfac)*dsfac,:,:) ,[dsfac floor(375/dsfac) n_rpt_trials n_tr_chs]);
Robs_tmat_DS = squeeze(nansum(temp));
Robs_mat_DS = reshape(Robs_tmat_DS,floor(375/dsfac)*n_rpt_trials,n_tr_chs);

%%
n_shuffs = 25;
shuff_corrs = zeros(n_tr_chs);
shuff_pcorrs = zeros(n_tr_chs);
shuff_pcorrs_DS = zeros(n_tr_chs);
for ii = 1:n_shuffs
    ii
    cur_cpred_rate = prate_tmat;
    cur_obs_rate = Robs_tmat;
    cur_obs_rate_DS = Robs_tmat_DS;
    for jj = 1:n_tr_chs
        cur_tord = randperm(n_rpt_trials);
        cur_cpred_rate(:,:,jj) = prate_tmat(:,cur_tord,jj);
        cur_obs_rate(:,:,jj) = cur_obs_rate(:,cur_tord,jj);
        cur_obs_rate_DS(:,:,jj) = cur_obs_rate_DS(:,cur_tord,jj);
    end
   shuff_corr_pred_rate_reshaped = reshape(cur_cpred_rate,n_rpt_trials*375,n_tr_chs);
   cur_corr = corr(shuff_corr_pred_rate_reshaped);
   shuff_corrs = shuff_corrs + cur_corr;
  
   shuff_psth_reshaped = reshape(cur_obs_rate,n_rpt_trials*375,n_tr_chs);
   cur_corr = corr(shuff_psth_reshaped);  shuff_pcorrs = shuff_pcorrs + cur_corr;
   shuff_pcorrs = shuff_pcorrs + cur_corr;

   shuff_psth_reshaped_DS = reshape(cur_obs_rate_DS,n_rpt_trials*floor(375/dsfac),n_tr_chs);
   cur_corr = corr(shuff_psth_reshaped_DS);
   shuff_pcorrs_DS = shuff_pcorrs_DS + cur_corr;
end
shuff_corrs = shuff_corrs/n_shuffs;
shuff_pcorrs = shuff_pcorrs/n_shuffs;
shuff_pcorrs_DS = shuff_pcorrs_DS/n_shuffs;

corrected_rate_corrs = corr(prate_mat) - shuff_corrs;
corrected_cnt_corrs = corr(Robs_mat)-shuff_pcorrs;
corrected_cnt_corrs_DS = corr(Robs_mat_DS)-shuff_pcorrs_DS;

corrected_rate_corrs(logical(eye(n_tr_chs))) = nan;
corrected_cnt_corrs(logical(eye(n_tr_chs))) = nan;
corrected_cnt_corrs_DS(logical(eye(n_tr_chs))) = nan;

[~,int_yord] = sort(comb_weights);

% stim_corrs = 

%%
fig_width = 3.27; %3.27 4.86 6.83
rel_height = 0.7;
fig_dir = '/home/james/Analysis/bruce/ET_final/';

figure
h = imagescnan(corrected_rate_corrs(int_yord,int_yord));
colorbar;
figufy(h);
box on
% caxis([-0.75 0.75])
% fname = [fig_dir 'sim_noisecorr_simp.pdf'];
caxis([-0.45 0.45])
fname = [fig_dir 'sim_noisecorr_SUs.pdf'];
exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close;

%%
fig_width = 3.27; %3.27 4.86 6.83o
rel_height = 0.8;

[I,J] = meshgrid(1:n_tr_chs);
p1 = [42 68];
p2 = [42 40];
b = robustfit(stim_corrmat(:),corrected_rate_corrs(:));
h = figure; hold on
plot(stim_corrmat(:),corrected_rate_corrs(:),'b.');
hold on
plot(stim_corrmat(p1(1),p1(2)),corrected_rate_corrs(p1(1),p1(2)),'ro','markersize',10,'linewidth',2)
plot(stim_corrmat(p1(1),p2(2)),corrected_rate_corrs(p2(1),p2(2)),'ro','markersize',10,'linewidth',2)

xlabel('Stimulus correlation');
ylabel('EM-induced noise correlation');
xlim([-0.6 1]); ylim([-0.6 1]);
xl = xlim();yl = ylim();
line([0 0],yl,'color','b','linestyle','--');
line(xl,[0 0],'color','b','linestyle','--');
line([-0.5 1],[-0.5 1],'color','k','linestyle','--');
xx = linspace(xl(1),xl(end),100);
% plot(xx,b(1) + b(2)*xx,'r--');

figufy(h);
fname = [fig_dir 'sim_noisecorr_SUs_scatter2.pdf'];
exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close;


%%

flen = 12;
dt = dt;

lag_axis = (0:flen-1)*dt*1e3;
pix_axis = (1:use_nPix_us)*sp_dx - use_nPix_us/2*sp_dx;
n_squared_filts = 2;
disp_pix = find(pix_axis >= -0.4 & pix_axis <= 0.4);

m1_filts = [comb_mods(p1(1)).mods(:).filtK];
max_vals1 = 0.9*max(abs(m1_filts));
m1_filts = reshape(m1_filts,[flen use_nPix_us n_squared_filts+1]);
m2_filts = [comb_mods(p1(2)).mods(:).filtK];
max_vals2 = 0.9*max(abs(m2_filts));
m2_filts = reshape(m2_filts,[flen use_nPix_us n_squared_filts+1]);
m3_filts = [comb_mods(p2(2)).mods(:).filtK];
max_vals3 = 0.9*max(abs(m3_filts));
m3_filts = reshape(m3_filts,[flen use_nPix_us n_squared_filts+1]);

h1 = figure();
for ii = 1:3;
subplot(3,3,3*(ii-1)+1)
imagesc(pix_axis(disp_pix),lag_axis,squeeze(m1_filts(:,disp_pix,ii))); caxis([-max_vals1(ii) max_vals1(ii)]);
colormap(gray)
xlabel('Position (deg)');
ylabel('Lag (ms)');

subplot(3,3,3*(ii-1)+2)
imagesc(pix_axis(disp_pix),lag_axis,squeeze(m2_filts(:,disp_pix,ii))); caxis([-max_vals2(ii) max_vals2(ii)]);
colormap(gray)
xlabel('Position (deg)');
ylabel('Lag (ms)');

subplot(3,3,3*(ii-1)+3)
imagesc(pix_axis(disp_pix),lag_axis,squeeze(m3_filts(:,disp_pix,ii))); caxis([-max_vals3(ii) max_vals3(ii)]);
colormap(gray)
xlabel('Position (deg)');
ylabel('Lag (ms)');

end

fig_width = 5; %3.27 4.86 6.83o
rel_height = 0.8;


figufy(h1);
fname = [fig_dir 'all_mod_compares.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h1);

