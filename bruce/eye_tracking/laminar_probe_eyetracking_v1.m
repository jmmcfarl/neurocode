clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

Expt_num = 266;
bar_ori = 80; %266
n_probes = 24;

Expt_name = sprintf('M%d',Expt_num);
data_dir = ['~/Data/bruce/' Expt_name];
cd(data_dir);

load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

mod_data_name = 'monoc_eyecorr_LP_mods';
anal_name = 'monoc_eyecorr_LP_ET';

recompute_init_mods = 0;
use_measured_pos = 0;
use_sac_kerns = 1;
use_LOOXV = 1;

ignore_blocks = [];
ignore_chs = 24;

%%
xv_frac = 0;

n_trial_inf_it = 0;
n_fix_inf_it = 5;

fix_prior_sigma = 0.3;
fix_noise_sigma = 100;

stim_fs = 100; %in Hz
dt = 0.01;
flen = 12;
use_nPix = 24;
full_nPix = 36;
stim_params = NIMcreate_stim_params([flen full_nPix],dt);
Fr = 1;
min_trial_dur = 0.75;

beg_buffer = 0.2;
end_buffer = 0.05;
trial_dur = 4;

use_right_eye = false;

n_use_blocks = Inf;

spatial_usfac = 2;
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

%%
include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaXFaXFs'};
for ii = 1:length(Expts)
    if ~isempty(Expts{ii})
        expt_names{ii} = Expts{ii}.Header.expname;
        expt_dds(ii) = Expts{ii}.Stimvals.dd;
        expt_bar_ori(ii) = Expts{ii}.Stimvals.or;
        expt_sac_dir(ii) = mod(Expts{ii}.Stimvals.Fa,180);
        expt_Fr(ii) = Expts{ii}.Stimvals.Fr;
        expt_imback(ii) = isfield(Expts{ii}.Trials,'imi');
        included_type(ii) = any(strcmp(expt_names{ii},include_expts));
    end
end
expt_has_ds(isnan(expt_has_ds)) = 0;
expt_has_ds(expt_has_ds == -1) = 0;
expt_binoc(isnan(expt_binoc)) = 0;
cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori);

if exist('./fin_aclust_data.mat','file')
    load('./fin_aclust_data.mat');
    [n_aclust_expts,n_aclust_probes] = size(autoclust);
else
    disp('No fin_aclust_data found.');
    autoclust = [];
    n_aclust_expts = 0; n_aclust_probes = 0;
end

n_blocks = length(cur_block_set);

sim_sac_expts = find(~expt_has_ds(cur_block_set));
imback_gs_expts = find(expt_has_ds(cur_block_set) & expt_imback(cur_block_set)');
grayback_gs_expts = find(expt_has_ds(cur_block_set) & ~expt_imback(cur_block_set)');

if length(cur_block_set) > n_use_blocks
    cur_block_set = cur_block_set(1:n_use_blocks);
end
n_blocks = length(cur_block_set);


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
all_clust_ids = cell(n_probes,1);
trial_toffset = zeros(length(cur_block_set),1);
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
            %             bar_Xmat = create_time_embedding(cur_stim_mat,stim_params);
            
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
    trial_toffset(ee) = all_t_bin_edges(end);
    cur_toffset = trial_toffset(ee);
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
all_Xmat = create_time_embedding(all_stim_mat,stim_params);
all_Xmat_us = create_time_embedding(all_stimmat_up,stim_params_us);

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
all_binned_mua = nan(length(all_t_axis),n_probes);
%for only-MU probes
for cc = 1:n_probes
    cur_mua_inds = all_clust_ids{cc} > 0;
    [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc}(cur_mua_inds),all_t_bin_edges);
    cur_spkhist(all_bin_edge_pts) = [];
    all_binned_mua(:,cc) = cur_spkhist;
end

%for SU probes
fprintf('Using %d SUs\n',length(SU_numbers));
all_binned_sua = nan(length(all_t_axis),length(SU_numbers));
cur_used_blocks = find(ismember(SU_target_blocks,cur_block_set));
SU_ID_mat = SU_ID_mat(cur_used_blocks,:);
[CC,BB] = meshgrid(1:length(SU_clust_data),1:length(cur_used_blocks));

su_probes = nan(1,length(SU_numbers));
for ss = 1:length(SU_numbers)
    used_clust_set = unique(CC(SU_ID_mat==ss)); %set of clusters used to capture this SU
    SU_block_probes(ss,:) = nan(1,length(cur_block_set));
    for cc = 1:length(used_clust_set)
        cur_clust = used_clust_set(cc);
        cur_probe = SU_clust_data(cur_clust).probe_num;
        cur_clust_label = SU_clust_data(cur_clust).cluster_label;
        cur_blocks = find(SU_ID_mat(:,cur_clust) == ss);
        SU_block_probes(ss,cur_blocks) = cur_probe;
        
        all_su_inds = all_clust_ids{cur_probe} == cur_clust_label;
        all_su_spk_times = all_spk_times{cur_probe}(all_su_inds);
        spk_block_inds = round(interp1(all_t_axis,all_blockvec,all_su_spk_times));
        all_su_spk_times = all_su_spk_times(ismember(spk_block_inds,cur_blocks));
        
        cur_suahist = histc(all_su_spk_times,all_t_bin_edges);
        cur_suahist(all_bin_edge_pts) = [];
        cur_id_set = ismember(all_blockvec,cur_blocks);
        all_binned_sua(cur_id_set,ss) = cur_suahist(cur_id_set);
    end
    su_probes(ss) = mode(SU_block_probes(ss,~isnan(SU_block_probes(ss,:))));
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
fix_start_inds(fix_durs==0) = [];
fix_stop_inds(fix_durs == 0) = [];
n_fixs = length(fix_start_inds);

%% INCORPORATE MEASURED EYE-POSITIONS
if use_measured_pos
    fprintf('Incorporating measured eye-corrections\n');
    
%     init_eyepos = corrected_eye_vals_interp(:,2);
    init_eyepos = 0.5*corrected_eye_vals_interp(:,2) + 0.5*corrected_eye_vals_interp(:,4);
    max_sim_pos = max_shift*sp_dx;
    init_eyepos = tanh(init_eyepos/max_sim_pos)*max_sim_pos;
    
    init_eyepos(isnan(init_eyepos)) = 0;
    
    %smooth out fast transients in eye signal
    eye_smooth_sig = round(0.025/dt);
    interp_inds = [];
    for ii = 1:n_fixs
        cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
        if length(cur_inds) > eye_smooth_sig*5;
            init_eyepos(used_inds(cur_inds)) = jmm_smooth_1d_cor(init_eyepos(used_inds(cur_inds)),eye_smooth_sig,2);
        end
        interp_inds = [interp_inds; cur_inds'];
    end
    interp_inds = unique(interp_inds);
    init_eyepos(used_inds) = interp1(used_inds(interp_inds),init_eyepos(used_inds(interp_inds)),used_inds);
    
    init_eyepos_rnd = round(init_eyepos/sp_dx);
    init_eyepos_ds = round(init_eyepos_rnd/spatial_usfac);
    for ii = 1:NT
        all_stimmat_up(used_inds(ii),:) = shift_matrix_Nd(all_stimmat_up(used_inds(ii),:),-init_eyepos_rnd(used_inds(ii)),2);
        all_stim_mat(used_inds(ii),:) = shift_matrix_Nd(all_stim_mat(used_inds(ii),:),-init_eyepos_ds(used_inds(ii)),2);
    end
    all_Xmat_us = create_time_embedding(all_stimmat_up,stim_params_us);
    all_Xmat = create_time_embedding(all_stim_mat,stim_params);
end

%% Create set of TR and XV trials
use_trials = unique(all_trialvec(used_inds));
nuse_trials = length(use_trials);

poss_tr_trials = find(~ismember(cur_block_set(all_trial_blocknums),ignore_blocks));
nposs_tr_trials = length(poss_tr_trials);
n_xv_trials = round(xv_frac*nposs_tr_trials);
xv_trials = randperm(nposs_tr_trials);
xv_trials(n_xv_trials+1:end) = [];
xv_trials = poss_tr_trials(xv_trials);
tr_trials = setdiff(poss_tr_trials,xv_trials);

tr_inds = find(ismember(all_trialvec(used_inds),tr_trials));
xv_inds = find(ismember(all_trialvec(used_inds),xv_trials));

tr_blocks = find(~ismember(cur_block_set,ignore_blocks));

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
    tot_nUnits = length(su_probes) + n_probes;
    all_mod_SU = zeros(tot_nUnits,1);
    for ss = 1:n_probes;
        fprintf('Computing base LLs for MU %d of %d\n',ss,n_probes);
        Robs = all_binned_mua(used_inds(tr_inds),ss);
        Robsxv = all_binned_mua(used_inds(xv_inds),ss);
        
        tr_X{1} = all_Xmat(used_inds(tr_inds),use_kInds);
        tr_X{2} = Xblock(used_inds(tr_inds),:);
        xv_X{1} = all_Xmat(used_inds(xv_inds),use_kInds);
        xv_X{2} = Xblock(used_inds(xv_inds),:);
        
        if use_sac_kerns
            tr_X{3} = Xsac(tr_inds,:);
            tr_X{4} = Xmsac(tr_inds,:);
            xv_X{3} = Xsac(xv_inds,:);
            xv_X{4} = Xmsac(xv_inds,:);
            null_mod = NMMinitialize_model(null_stim_params,[1 1 1],repmat({'lin'},1,3),null_reg_params,[1 2 3]);
        else
            null_mod = NMMinitialize_model(null_stim_params,[1],{'lin'},null_reg_params,[1]);
        end
        
        null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));
        all_nullmod(ss) = null_mod;
        
        gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
        gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent,init_optim_p);
        
        tr_X{1} = all_Xmat_us(used_inds(tr_inds),use_kInds_up);
        xv_X{1} = all_Xmat_us(used_inds(xv_inds),use_kInds_up);
        %spatial up-sampling of filter estimates
        base_filts = reshape([gqm1.mods(find(init_Xtargs == 1)).filtK],[flen use_nPix n_squared_filts+1]);
        if spatial_usfac == 2
            base_filts_up = zeros(flen,use_nPix_us,n_squared_filts+1);
            for ii = 1:use_nPix
                base_filts_up(:,2*(ii-1)+1,:) = 0.5*base_filts(:,ii,:);
                base_filts_up(:,2*(ii-1)+2,:) = 0.5*base_filts(:,ii,:);
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
        
        all_mod_fits(ss) = gqm2;
        all_mod_fits_withspkNL(ss) = NMMfit_logexp_spkNL(gqm2,Robs,tr_X);
        
        LL = NMMmodel_eval(all_mod_fits_withspkNL(ss),Robs,tr_X);
        null_LL(ss) = NMMmodel_eval(null_mod,Robs,tr_X(2:end));
        all_mod_LLimp(ss) = (LL-null_LL(ss))/log(2);
        
        truexvLL = NMMmodel_eval(all_mod_fits_withspkNL(ss),Robsxv,xv_X);
        null_xvLL(ss) = NMMmodel_eval(null_mod,Robsxv,xv_X(2:end));
        all_mod_truexvLLimp(ss) = (truexvLL-null_xvLL(ss))/log(2);
    end
    
    for ss = 1:length(su_probes);
        fprintf('Computing base LLs for SU %d of %d\n',ss,length(su_probes));
        all_mod_SU(ss+n_probes) = su_probes(ss);
        cur_used_blocks = tr_blocks(~isnan(SU_block_probes(ss,tr_blocks)));
        cur_poss_inds = find(ismember(all_blockvec,cur_used_blocks));
        cur_tr_inds = tr_inds(ismember(used_inds(tr_inds),cur_poss_inds));
        cur_xv_inds = xv_inds(ismember(used_inds(xv_inds),cur_poss_inds));
        tr_NT = length(cur_tr_inds);
        xv_NT = length(cur_xv_inds);
        if ~isempty(cur_tr_inds)
            Robs = all_binned_sua(used_inds(cur_tr_inds),ss);
            Robsxv = all_binned_sua(used_inds(cur_xv_inds),ss);
            
            tr_X{1} = all_Xmat(used_inds(cur_tr_inds),use_kInds);
            tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
            xv_X{1} = all_Xmat(used_inds(cur_xv_inds),use_kInds);
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
            all_nullmod(ss+n_probes) = null_mod;
            
            gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
            gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent,init_optim_p);
            
            tr_X{1} = all_Xmat_us(used_inds(cur_tr_inds),use_kInds_up);
            xv_X{1} = all_Xmat_us(used_inds(cur_xv_inds),use_kInds_up);
            %spatial up-sampling of filter estimates
            base_filts = reshape([gqm1.mods(find(init_Xtargs == 1)).filtK],[flen use_nPix n_squared_filts+1]);
            if spatial_usfac == 2
                base_filts_up = zeros(flen,use_nPix_us,n_squared_filts+1);
                for ii = 1:use_nPix
                    base_filts_up(:,2*(ii-1)+1,:) = 0.5*base_filts(:,ii,:);
                    base_filts_up(:,2*(ii-1)+2,:) = 0.5*base_filts(:,ii,:);
                end
            elseif spatial_usfac == 1
                base_filts_up = base_filts;
            else
                error('unsupported')
            end
            base_filts_up = reshape(base_filts_up,use_nPix_us*flen,	n_squared_filts+1);
            
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
            
            all_mod_fits(ss+n_probes) = gqm2;
            
            all_mod_fits_withspkNL(ss+n_probes) = NMMfit_logexp_spkNL(gqm2,Robs,tr_X);
            
            LL = NMMmodel_eval(all_mod_fits_withspkNL(ss+n_probes),Robs,tr_X);
            null_LL(ss+n_probes) = NMMmodel_eval(null_mod,Robs,tr_X(2:end));
            all_mod_LLimp(ss+n_probes) = (LL-null_LL(ss+n_probes))/log(2);
            
            truexvLL = NMMmodel_eval(all_mod_fits_withspkNL(ss+n_probes),Robsxv,xv_X);
            null_xvLL(ss+n_probes) = NMMmodel_eval(null_mod,Robsxv,xv_X(2:end));
            all_mod_truexvLLimp(ss+n_probes) = (truexvLL-null_xvLL(ss+n_probes))/log(2);
        end
    end
    if xv_frac > 0
        save(mod_data_name,'all_mod*','all_nullmod','SU_*','null_xvLL','null_LL','*_trials');
    else
        save(mod_data_name,'all_mod*','all_nullmod','SU_*','null_xvLL','null_LL','*_trials');
    end
    
else
    fprintf('Loading pre-computed initial models\n');
    load(mod_data_name);
    
    tr_inds = find(ismember(all_trialvec(used_inds),tr_trials));
    xv_inds = find(ismember(all_trialvec(used_inds),xv_trials));
end

%% SELECT USABLE UNITS AND make Robs_mat
if xv_frac == 0
    LL_imp_thresh = 2.5e-3;
    usable_units = find(all_mod_LLimp >= LL_imp_thresh);
else
    LL_imp_thresh = 0;
    usable_units = find(all_mod_truexvLLimp > LL_imp_thresh);
end

%exclude MU probes which also have a used SU (MAY WANT TO DO THIS ON A
%BLOCK-BY-BLOCK BASIS
su_inds = find(all_mod_SU > 0);
used_sus = find(ismember(su_inds,usable_units));
mu_copies = su_probes(used_sus);
usable_units(ismember(usable_units,mu_copies)) = [];

usable_units(ismember(usable_units,ignore_chs)) = [];
usable_units(ismember(all_mod_SU(usable_units),ignore_chs)) = [];

n_used_sus = sum(all_mod_SU(usable_units) ~= 0);
n_used_mus = sum(all_mod_SU(usable_units) == 0);
fprintf('Using %d SUs and %d MUs for analysis\n',n_used_sus,n_used_mus);
tr_set = usable_units;
n_tr_chs = length(tr_set);

% make Robs_mat
poss_blocks = 1:(n_blocks-1);
Robs_mat = nan(length(used_inds),n_tr_chs);
for ss = 1:n_tr_chs
    if all_mod_SU(tr_set(ss)) == 0
        Robs_mat(:,ss) = all_binned_mua(used_inds,tr_set(ss));
    else
        su_ind = find(su_probes == all_mod_SU(tr_set(ss)));
        cur_used_blocks = find(~isnan(SU_block_probes(su_ind,:)));
        cur_use = find(ismember(all_blockvec(used_inds),cur_used_blocks));
        Robs_mat(cur_use,ss) = all_binned_sua(used_inds(cur_use),su_ind);
    end
end

%% DIAGNOSTICS
n_trials = length(trial_start_inds);
% cur_X{1} = all_Xmat_us(used_inds,use_kInds_up);
% cur_X{2} = Xblock(used_inds,:);
% if use_sac_kerns
%     cur_X{3} = Xsac;
%     cur_X{4} = Xmsac;
% end
% full_prates = nan(NT,n_tr_chs);
% % full_xvprates = nan(NT,n_tr_chs);
% full_nullrates = nan(NT,n_tr_chs);
% for cc = 1:n_tr_chs
%     [~,~,full_prates(:,cc)] = NMMmodel_eval(all_mod_fits_withspkNL(tr_set(cc)),Robs_mat(:,cc),cur_X);
%     %     [~,~,full_xvprates(:,cc)] = NMMmodel_eval(all_xvmod_fits(tr_set(cc)),Robs_mat(:,cc),cur_X);
%     [~,~,full_nullrates(:,cc)] = NMMmodel_eval(all_nullmod(tr_set(cc)),Robs_mat(:,cc),cur_X(2:end));
% end
% full_modLL = Robs_mat.*log(full_prates) - full_prates;
% full_nullLL = Robs_mat.*log(full_nullrates) - full_nullrates;
% full_LLimp = full_modLL-full_nullLL;
% 
% trial_blocknums = nan(nuse_trials,1);
% trial_LLimp = nan(nuse_trials,n_tr_chs);
% trial_meanrate = nan(nuse_trials,n_tr_chs);
% trial_nspks = nan(nuse_trials,n_tr_chs);
% trial_durs = nan(nuse_trials,1);
% trial_eyeerr = nan(nuse_trials,1);
% for tt = 1:nuse_trials
%     cur_used_inds = find(all_trialvec(used_inds) == use_trials(tt));
%     trial_LLimp(tt,:) = nanmean(full_LLimp(cur_used_inds,:));
%     trial_meanrate(tt,:) = nanmean(Robs_mat(cur_used_inds,:));
%     trial_nspks(tt,:) = nansum(Robs_mat(cur_used_inds,:));
%     trial_blocknums(tt) = unique(all_blockvec(used_inds(cur_used_inds)));
%     trial_durs(tt) = length(cur_used_inds);
%     trial_eyeerr(tt) = mean(corrected_eye_vals_interp(used_inds(cur_used_inds),2));
% end
% 
% block_LLimp = nan(n_blocks,n_tr_chs);
% block_meanrate = nan(n_blocks,n_tr_chs);
% block_nspks = nan(n_blocks,n_tr_chs);
% for tt = 1:n_blocks
%     cur_used_inds = find(all_blockvec(used_inds) == tt);
%     block_LLimp(tt,:) = nanmean(full_LLimp(cur_used_inds,:));
%     block_meanrate(tt,:) = nanmean(Robs_mat(cur_used_inds,:));
%     block_nspks(tt,:) = nansum(Robs_mat(cur_used_inds,:));
% end
% 
% tr_trialset = find(ismember(use_trials,tr_trials));
% xv_trialset = find(ismember(use_trials,xv_trials));

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

% leps_prior = -(shifts*sp_dx).^2/(2*fix_prior_sigma^2);
% leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
%overall prior on shifts
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

%%
warning('OFF','stats:statrobustfit:IterationLimit');
% measured_seq = corrected_eye_vals_interp(used_inds,2);
measured_seq = 0.5*corrected_eye_vals_interp(used_inds,2) + 0.5*corrected_eye_vals_interp(used_inds,4);
measured_drift = nan(size(used_inds));
measured_fix_avg = nan(n_fixs,1);
measured_drift_avg = zeros(n_fixs,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    measured_fix_avg(ii) = median(measured_seq(cur_inds));
    measured_drift(cur_inds) = measured_seq(cur_inds) - measured_fix_avg(ii);
    if length(cur_inds) > 10
        temp = robustfit(1:length(cur_inds),measured_seq(cur_inds));
        measured_drift_avg(ii) = temp(2);
    end
end

post_fix_var = 1/(1/fix_noise_sigma^2 + 1/fix_prior_sigma^2);
post_fix_sigma = sqrt(post_fix_var);

fix_prior = diff(erf(ovshift_dx_edges/(fix_prior_sigma*sqrt(2))));
fix_prior = log(fix_prior/sum(fix_prior));
fix_prior(fix_prior < eps) = eps;

measured_fix_deltas = nan(n_fixs,1);
measured_fix_deltas(2:end) = diff(measured_fix_avg);
post_mean_delta = post_fix_var*measured_fix_deltas/fix_noise_sigma^2;


%% ITERATE FIXATION-BASED CORRECTIONS
su_inds = find(all_mod_SU(tr_set) > 0);

Tit_mods{1} = all_mod_fits;
Tit_mods_spkNL{1} = all_mod_fits_withspkNL;
Tit_LLimp(1,:) = all_mod_LLimp;
Tit_xvLLimp(1,:) = all_mod_truexvLLimp;
it_trial_sigma(1) = fix_prior_sigma;
if use_LOOXV
    Tit_mods_LOO{1} = all_mod_fits;
    Tit_mods_spkNL_LOO{1} = all_mod_fits_withspkNL;
    Tit_LLimp_LOO(1,:) = all_mod_LLimp;
    Tit_xvLLimp_LOO(1,:) = all_mod_truexvLLimp;
end
for nn = 1:n_trial_inf_it
    fprintf('Inferring trial corrections, iter %d of %d\n',nn,n_trial_inf_it);
    
    %% PREPROCESS MODEL COMPONENTS
    filt_bank = zeros(n_tr_chs,klen_us,n_squared_filts+1);
    lin_kerns = nan(n_tr_chs,n_blocks);
    if use_sac_kerns
        sac_kerns = nan(n_tr_chs,n_sac_bins);
        msac_kerns = nan(n_tr_chs,n_sac_bins);
    end
    mod_spkNL_params = nan(n_tr_chs,3);
    for ss = 1:n_tr_chs
        cur_Xtargs = [Tit_mods{nn}(tr_set(ss)).mods(:).Xtarget];
        cur_k = [Tit_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 1).filtK];
        n_used_filts = size(cur_k,2);
        filt_bank(ss,:,1:n_used_filts) = cur_k;
        mod_spkNL_params(ss,:) = Tit_mods_spkNL{nn}(tr_set(ss)).spk_NL_params;
        lin_kerns(ss,:) = Tit_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 2).filtK;
        if use_sac_kerns
            sac_kerns(ss,:) = Tit_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 3).filtK;
            msac_kerns(ss,:) = Tit_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 4).filtK;
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
    all_frame_LLs = nan(NT,n_tr_chs,n_shifts);
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
        
        all_frame_LLs(:,:,xx) = Robs_mat.*log(pred_rate) - pred_rate;
    end
    frame_LLs = squeeze(nansum(all_frame_LLs,2));
    
    %% INFER MICRO-SAC SEQUENCE
    trial_LLs = nan(n_trials,n_shifts);
    for ii = 1:n_trials
        cur_inds = trial_start_inds(ii):trial_end_inds(ii);
        trial_LLs(ii,:) = sum(frame_LLs(cur_inds,:));
    end
    
    lPost = bsxfun(@plus,trial_LLs,fix_prior);
    lPost = bsxfun(@minus,lPost,logsumexp(lPost,2));
    Post = exp(lPost);
    it_trial_sigma(nn) = sqrt(mean(sum(bsxfun(@times,Post,shifts*sp_dx).^2,2)));
    
    it_trial_post_mean(nn,:) = sum(bsxfun(@times,Post,shifts),2);
    cur_diff = bsxfun(@minus,it_trial_post_mean(nn,:)',shifts).^2;
    it_trial_post_std(nn,:) = sqrt(sum(cur_diff.*Post,2));
    
    %back-project saccade-times
    all_trial_post_mean_cor = it_trial_post_mean(nn,trial_ids);
    
    %% RECOMPUTE XMAT
    all_shift_stimmat_up = all_stimmat_up;
    for i=1:NT
        all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-round(all_trial_post_mean_cor(i)),2);
    end
    all_Xmat_up_trialcor = create_time_embedding(all_shift_stimmat_up,stim_params_us);
    all_Xmat_up_trialcor = all_Xmat_up_trialcor(used_inds,:);
    
    %% REFIT ALL CELLS
    cur_X{1} = all_Xmat_up_trialcor(:,use_kInds_up);
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
        cur_xv_uset = xv_inds(~isnan(Robs_mat(xv_inds,cur_unit_ind)));
        
        tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
        xv_X = get_Xcell_tInds(cur_X,cur_xv_uset);
        
        Tit_mods{nn+1}(cur_cell) = Tit_mods{nn}(cur_cell);
        Tit_mods{nn+1}(cur_cell) = NMMfit_filters(Tit_mods{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),...
            tr_X,[],[],silent); %fit stimulus filters
        
        %refit spk NL
        Tit_mods_spkNL{nn+1}(cur_cell) = NMMfit_logexp_spkNL(Tit_mods{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
        
        newLL = NMMmodel_eval(Tit_mods_spkNL{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
        truenewxvLL = NMMmodel_eval(Tit_mods_spkNL{nn+1}(cur_cell),Robs_mat(cur_xv_uset,cur_unit_ind),xv_X);
        
        Tit_LLimp(nn+1,cur_cell) = (newLL - null_LL(cur_cell))/log(2);
        Tit_xvLLimp(nn+1,cur_cell) = (truenewxvLL - null_xvLL(cur_cell))/log(2);
        if nn > 1
            fprintf('Original: %.4f  Prev: %.4f  New: %.4f\n',all_mod_LLimp(cur_cell),Tit_LLimp(nn,cur_cell),Tit_LLimp(nn+1,cur_cell));
        else
            fprintf('Original: %.4f  New: %.4f\n',all_mod_LLimp(cur_cell),Tit_LLimp(nn+1,cur_cell));
        end
    end
    
    if use_LOOXV == 1
        for xv = 1:length(su_inds)
            fprintf('Inferring drift corrections, XV %d of %d\n',xv,length(su_inds));
            
            %% PREPROCESS MODEL COMPONENTS
            cur_xv_cell = tr_set(su_inds(xv));
            cur_tr_cells = setdiff(tr_set,cur_xv_cell);
            cur_tr_cinds = find(ismember(tr_set,cur_tr_cells));

            frame_LLs = squeeze(nansum(all_frame_LLs(:,cur_tr_cinds,:),2));
            
            %% INFER MICRO-SAC SEQUENCE            
            trial_LLs = nan(n_trials,n_shifts);
            for ii = 1:n_trials
                cur_inds = trial_start_inds(ii):trial_end_inds(ii);
                trial_LLs(ii,:) = sum(frame_LLs(cur_inds,:));
            end
            
            lPost = bsxfun(@plus,trial_LLs,fix_prior);
            lPost = bsxfun(@minus,lPost,logsumexp(lPost,2));
            Post = exp(lPost);
            it_trial_sigma_LOO(nn,xv) = sqrt(mean(sum(bsxfun(@times,Post,shifts*sp_dx).^2,2)));
            
            it_trial_post_mean_LOO(xv,nn,:) = sum(bsxfun(@times,Post,shifts),2);
            cur_diff = bsxfun(@minus,squeeze(it_trial_post_mean_LOO(xv,nn,:)),shifts).^2;
            it_trial_post_std_LOO(xv,nn,:) = sqrt(sum(cur_diff.*Post,2));
            
            %back-project saccade-times
            all_trial_post_mean_cor = squeeze(it_trial_post_mean_LOO(xv,nn,trial_ids));
            
            %% RECOMPUTE XMAT
            all_shift_stimmat_up = all_stimmat_up;
            for i=1:NT
                all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-round(all_trial_post_mean_cor(i)),2);
            end
            cur_Xmat = create_time_embedding(all_shift_stimmat_up,stim_params_us);
            cur_Xmat = cur_Xmat(used_inds,:);
            
            %% REFIT XV CELLS
            cur_X{1} = cur_Xmat(:,use_kInds_up);
            
            silent = 1;
            cur_cell = cur_xv_cell;
            fprintf('Refitting model \n');
            cur_unit_ind = find(tr_set == cur_cell);
            cur_tr_uset = tr_inds(~isnan(Robs_mat(tr_inds,cur_unit_ind)));
            cur_xv_uset = xv_inds(~isnan(Robs_mat(xv_inds,cur_unit_ind)));
            
            tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
            xv_X = get_Xcell_tInds(cur_X,cur_xv_uset);
            
            Tit_mods_LOO{nn+1}(cur_cell) = Tit_mods_LOO{nn}(cur_cell);
            Tit_mods_LOO{nn+1}(cur_cell) = NMMfit_filters(Tit_mods_LOO{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),...
                tr_X,[],[],silent); %fit stimulus filters
            
            %refit spk NL
            Tit_mods_spkNL_LOO{nn+1}(cur_cell) = NMMfit_logexp_spkNL(Tit_mods_LOO{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
            
            newLL = NMMmodel_eval(Tit_mods_spkNL_LOO{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
            newxvLL = NMMmodel_eval(Tit_mods_spkNL_LOO{nn+1}(cur_cell),Robs_mat(cur_xv_uset,cur_unit_ind),xv_X);
            
            Tit_LLimp_LOO(nn+1,cur_cell) = (newLL - null_LL(cur_cell))/log(2);
            Tit_xvLLimp_LOO(nn+1,cur_cell) = (newxvLL - null_xvLL(cur_cell))/log(2);
        end
    end
end
%%
if n_trial_inf_it == 0
    all_Xmat_up_trialcor = all_Xmat_us(used_inds,:);
    it_trial_post_mean = zeros(1,n_trials);
    it_trial_post_std = zeros(1,n_trials);
    all_trial_post_mean_cor = zeros(NT,1);
    if use_LOOXV
        it_trial_post_mean_LOO = zeros(length(su_probes),1,n_trials);
        it_trial_post_std_LOO = zeros(length(su_probes),1,n_trials);
    end
end
%% ITERATE FIXATION-BASED CORRECTIONS
su_inds = find(all_mod_SU(tr_set) > 0);

it_mods{1} = Tit_mods{n_trial_inf_it+1};
it_mods_spkNL{1} = Tit_mods_spkNL{n_trial_inf_it+1};
it_LLimp(1,:) = Tit_LLimp(n_trial_inf_it+1,:);
it_xvLLimp(1,:) = Tit_xvLLimp(n_trial_inf_it+1,:);
it_fix_sigma(1) = fix_prior_sigma;
if use_LOOXV
    it_mods_LOO{1} = Tit_mods_LOO{n_trial_inf_it+1};
    it_mods_spkNL_LOO{1} = Tit_mods_spkNL_LOO{n_trial_inf_it+1};
    it_LLimp_LOO(1,:) = Tit_LLimp_LOO(n_trial_inf_it+1,:);
    it_xvLLimp_LOO(1,:) = Tit_xvLLimp_LOO(n_trial_inf_it+1,:);
end
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
        cur_Xtargs = [it_mods{nn}(tr_set(ss)).mods(:).Xtarget];
        cur_k = [it_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 1).filtK];
        n_used_filts = size(cur_k,2);
        filt_bank(ss,:,1:n_used_filts) = cur_k;
        mod_spkNL_params(ss,:) = it_mods_spkNL{nn}(tr_set(ss)).spk_NL_params;
        lin_kerns(ss,:) = it_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 2).filtK;
        if use_sac_kerns
            sac_kerns(ss,:) = it_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 3).filtK;
            msac_kerns(ss,:) = it_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 4).filtK;
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
    %precompute LL at all shifts for all units
    all_frame_LLs = nan(NT,n_tr_chs,n_shifts);
    for xx = 1:length(shifts)
        fprintf('Shift %d of %d\n',xx,n_shifts);
        cur_stim_shift = all_Xmat_up_trialcor*shift_mat{xx};
        
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
        
        all_frame_LLs(:,:,xx) = Robs_mat.*log(pred_rate) - pred_rate;
    end
    frame_LLs = squeeze(nansum(all_frame_LLs,2));
    
    %% INFER MICRO-SAC SEQUENCE
    fix_LLs = nan(n_fixs,n_shifts);
    for ii = 1:n_fixs
        cur_inds = pfix_start_inds(ii):pfix_stop_inds(ii);
        fix_LLs(ii,:) = sum(frame_LLs(cur_inds,:));
    end
    
    lalpha=zeros(n_fixs,n_shifts);
    lbeta = zeros(n_fixs,n_shifts);
    %compute rescaled forward messages
    lalpha(1,:) = fix_prior + fix_LLs(1,:);
    for t=2:n_fixs
        gprobs = diff(erf((shift_dx_edges + post_mean_delta(t))/(post_fix_sigma*sqrt(2))));
        cur_lA = toeplitz(gprobs(n_shifts:end)',gprobs(n_shifts:-1:1)');
        cur_lA = log(bsxfun(@rdivide,cur_lA,sum(cur_lA,2)));
        cur_lA(isnan(cur_lA)) = eps;
        cur_lA(cur_lA < eps) = eps;
        lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + fix_LLs(t,:);
    end
    
    %compute rescaled backward messages
    lbeta(n_fixs,:)=log(ones(1,n_shifts));
    for t=n_fixs-1:-1:1
         gprobs = diff(erf((shift_dx_edges + post_mean_delta(t+1))/(post_fix_sigma*sqrt(2))));
        cur_lA = toeplitz(gprobs(n_shifts:end),gprobs(n_shifts:-1:1));
        cur_lA = log(bsxfun(@rdivide,cur_lA,sum(cur_lA,2)));
        cur_lA(isnan(cur_lA)) = eps;
        cur_lA(cur_lA < eps) = eps;
        lf1 = lbeta(t+1,:) + fix_LLs(t+1,:);
        lbeta(t,:) = logmulexp(lf1,cur_lA');
    end
    lgamma= lalpha + lbeta;
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    gamma = exp(lgamma); 
    
    it_fix_post_mean(nn,:) = sum(bsxfun(@times,gamma,shifts),2);
    cur_diff = bsxfun(@minus,it_fix_post_mean(nn,:)',shifts).^2;
    it_fix_post_std(nn,:) = sqrt(sum(cur_diff.*gamma,2));
    
    %back-project saccade-times
    all_fix_post_mean_cor = nan(NT,1);
    all_fix_post_mean_cor(~isnan(fix_ids)) = it_fix_post_mean(nn,fix_ids(~isnan(fix_ids)));
    all_fix_post_mean_cor = interp1(find(~isnan(fix_ids)),all_fix_post_mean_cor(~isnan(fix_ids)),1:NT);
    
    %back-project saccade-times
    all_trial_post_mean_cor = squeeze(it_trial_post_mean(end,trial_ids));
    all_fix_post_mean_cor = all_fix_post_mean_cor + all_trial_post_mean_cor;
    
    %% RECOMPUTE XMAT
    all_shift_stimmat_up = all_stimmat_up;
    for i=1:NT
        all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-round(all_fix_post_mean_cor(i)),2);
    end
    all_Xmat_up_fixcor = create_time_embedding(all_shift_stimmat_up,stim_params_us);
    all_Xmat_up_fixcor = all_Xmat_up_fixcor(used_inds,:);
    
    %% REFIT ALL CELLS
    cur_X{1} = all_Xmat_up_fixcor(:,use_kInds_up);
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
        cur_xv_uset = xv_inds(~isnan(Robs_mat(xv_inds,cur_unit_ind)));
        
        tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
        xv_X = get_Xcell_tInds(cur_X,cur_xv_uset);
        
        it_mods{nn+1}(cur_cell) = it_mods{nn}(cur_cell);
        it_mods{nn+1}(cur_cell) = NMMfit_filters(it_mods{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),...
            tr_X,[],[],silent); %fit stimulus filters
        
        %refit spk NL
        it_mods_spkNL{nn+1}(cur_cell) = NMMfit_logexp_spkNL(it_mods{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
        
        newLL = NMMmodel_eval(it_mods_spkNL{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
        truenewxvLL = NMMmodel_eval(it_mods_spkNL{nn+1}(cur_cell),Robs_mat(cur_xv_uset,cur_unit_ind),xv_X);
        
        it_LLimp(nn+1,cur_cell) = (newLL - null_LL(cur_cell))/log(2);
        it_xvLLimp(nn+1,cur_cell) = (truenewxvLL - null_xvLL(cur_cell))/log(2);
        if nn > 1
            fprintf('Original: %.4f  Prev: %.4f  New: %.4f\n',all_mod_LLimp(cur_cell),it_LLimp(nn,cur_cell),it_LLimp(nn+1,cur_cell));
        else
            fprintf('Original: %.4f  New: %.4f\n',all_mod_LLimp(cur_cell),it_LLimp(nn+1,cur_cell));
        end
    end
    
    %%
    if use_LOOXV == 1
        for xv = 1:length(su_inds)
            fprintf('Inferring drift corrections, XV %d of %d\n',xv,length(su_inds));
            
            %% PREPROCESS MODEL COMPONENTS
            cur_xv_cell = tr_set(su_inds(xv));
            cur_tr_cells = setdiff(tr_set,cur_xv_cell);
            cur_tr_cinds = find(ismember(tr_set,cur_tr_cells));
            frame_LLs = squeeze(nansum(all_frame_LLs(:,cur_tr_cinds,:),2));
            
            %% INFER MICRO-SAC SEQUENCE            
            fix_LLs = nan(n_fixs,n_shifts);
            for ii = 1:n_fixs
                cur_inds = pfix_start_inds(ii):pfix_stop_inds(ii);
                fix_LLs(ii,:) = sum(frame_LLs(cur_inds,:));
            end
            
            lalpha=zeros(n_fixs,n_shifts);
            lbeta = zeros(n_fixs,n_shifts);
            %compute rescaled forward messages
            lalpha(1,:) = fix_prior + fix_LLs(1,:);
            for t=2:n_fixs
                gprobs = diff(erf((shift_dx_edges + post_mean_delta(t))/(post_fix_sigma*sqrt(2))));
                cur_lA = toeplitz(gprobs(n_shifts:end)',gprobs(n_shifts:-1:1)');
                cur_lA = log(bsxfun(@rdivide,cur_lA,sum(cur_lA,2)));
                cur_lA(isnan(cur_lA)) = eps;
                cur_lA(cur_lA < eps) = eps;
                lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + fix_LLs(t,:);
            end
            
            %compute rescaled backward messages
            lbeta(n_fixs,:)=log(ones(1,n_shifts));
            for t=n_fixs-1:-1:1
                gprobs = diff(erf((shift_dx_edges + post_mean_delta(t+1))/(post_fix_sigma*sqrt(2))));
                cur_lA = toeplitz(gprobs(n_shifts:end),gprobs(n_shifts:-1:1));
                cur_lA = log(bsxfun(@rdivide,cur_lA,sum(cur_lA,2)));
                cur_lA(isnan(cur_lA)) = eps;
                cur_lA(cur_lA < eps) = eps;
                lf1 = lbeta(t+1,:) + fix_LLs(t+1,:);
                lbeta(t,:) = logmulexp(lf1,cur_lA');
            end
            lgamma= lalpha + lbeta;
            lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
            gamma = exp(lgamma);
            
            it_fix_post_mean_LOO(xv,nn,:) = sum(bsxfun(@times,gamma,shifts),2);
            cur_diff = bsxfun(@minus,squeeze(it_fix_post_mean_LOO(xv,nn,:)),shifts).^2;
            it_fix_post_std_LOO(xv,nn,:) = sqrt(sum(cur_diff.*gamma,2));
            
            %back-project saccade-times
            all_fix_post_mean_cor = nan(NT,1);
            all_fix_post_mean_cor(~isnan(fix_ids)) = squeeze(it_fix_post_mean_LOO(xv,nn,fix_ids(~isnan(fix_ids))));
            all_fix_post_mean_cor = interp1(find(~isnan(fix_ids)),all_fix_post_mean_cor(~isnan(fix_ids)),1:NT);
            
            %back-project saccade-times
            all_trial_post_mean_cor = squeeze(it_trial_post_mean(end,trial_ids));
            all_fix_post_mean_cor = all_fix_post_mean_cor + all_trial_post_mean_cor;

            %% RECOMPUTE XMAT
            all_shift_stimmat_up = all_stimmat_up;
            for i=1:NT
                all_shift_stimmat_up(used_inds(i),:) = shift_matrix_Nd(all_stimmat_up(used_inds(i),:),-round(all_fix_post_mean_cor(i)),2);
            end
            all_Xmat_up_fixcor_LOO = create_time_embedding(all_shift_stimmat_up,stim_params_us);
            all_Xmat_up_fixcor_LOO = all_Xmat_up_fixcor_LOO(used_inds,:);
            
            %% REFIT XV CELLS
            cur_X{1} = all_Xmat_up_fixcor_LOO(:,use_kInds_up);
            
            silent = 1;
            cur_cell = cur_xv_cell;
            fprintf('Refitting model \n');
            cur_unit_ind = find(tr_set == cur_cell);
            cur_tr_uset = tr_inds(~isnan(Robs_mat(tr_inds,cur_unit_ind)));
            cur_xv_uset = xv_inds(~isnan(Robs_mat(xv_inds,cur_unit_ind)));
            
            tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
            xv_X = get_Xcell_tInds(cur_X,cur_xv_uset);
            
            it_mods_LOO{nn+1}(cur_cell) = it_mods_LOO{nn}(cur_cell);
            it_mods_LOO{nn+1}(cur_cell) = NMMfit_filters(it_mods_LOO{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),...
                tr_X,[],[],silent); %fit stimulus filters
            
            %refit spk NL
            it_mods_spkNL_LOO{nn+1}(cur_cell) = NMMfit_logexp_spkNL(it_mods_LOO{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
            
            newLL = NMMmodel_eval(it_mods_spkNL_LOO{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
            newxvLL = NMMmodel_eval(it_mods_spkNL_LOO{nn+1}(cur_cell),Robs_mat(cur_xv_uset,cur_unit_ind),xv_X);
            
            it_LLimp_LOO(nn+1,cur_cell) = (newLL - null_LL(cur_cell))/log(2);
            it_xvLLimp_LOO(nn+1,cur_cell) = (newxvLL - null_xvLL(cur_cell))/log(2);
        end
    end
end

%% SAVE EYE-TRACKING RESULTS
et_params = struct('beg_buffer',beg_buffer,'end_buffer',end_buffer,'min_trial_dur',min_trial_dur,'bar_ori',bar_ori,...
    'use_nPix',use_nPix,'flen',flen,'dt',dt,...
    'fix_prior_sigma',fix_prior_sigma,'n_fix_inf_it',n_fix_inf_it,'n_trial_inf_it',n_trial_inf_it,'use_sac_kerns',use_sac_kerns,'shifts',shifts,...
    'use_measured_pos',use_measured_pos,'sac_bincents',sac_bincents,'spatial_usfac',spatial_usfac,'sac_shift',sac_shift);

et_used_inds = used_inds;
et_tr_set = tr_set;
et_tr_trials = tr_trials;
et_xv_trials = xv_trials;
et_saccade_inds = saccade_start_inds;
cd(anal_dir);
save(anal_name,'it_*','Tit_*','fix_ids','trial_ids','et_used_inds','et_tr_set','et_saccade_inds','et_params','et_tr_trials','et_xv_trials');

%%
usable_inds = find(~ismember(all_blockvec(used_inds),ignore_blocks));

fin_fix_corr = nan(NT,1);
fin_fix_std = nan(NT,1);
fin_fix_corr(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
fin_fix_std(~isnan(fix_ids)) = it_fix_post_std(end,fix_ids(~isnan(fix_ids)));
fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);

fin_fix_corr = fin_fix_corr*sp_dx;
fin_fix_std = fin_fix_std*sp_dx;

fin_trial_corr = it_trial_post_mean(end,trial_ids);
fin_trial_std = it_trial_post_std(end,trial_ids);
fin_trial_corr = fin_trial_corr*sp_dx;
fin_trial_std = fin_trial_std*sp_dx;

fin_tot_corr = fin_fix_corr + fin_trial_corr;
fin_tot_std = sqrt(fin_fix_std.^2 + fin_trial_std.^2);

if use_measured_pos==1
    fin_tot_corr = fin_tot_corr + init_eyepos(used_inds)';
end

%%

measured_seqL = corrected_eye_vals_interp(used_inds,2);
measured_seqR = corrected_eye_vals_interp(used_inds,4);

min_fix_dur = 0.15;
inferred_drift = nan(size(fin_tot_corr));
measured_driftL = nan(size(fin_tot_corr));
measured_driftR = nan(size(fin_tot_corr));
inferred_fix_avg = nan(n_fixs,1);
measured_fix_avgL = nan(n_fixs,1);
measured_fix_avgR = nan(n_fixs,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds)*dt >= min_fix_dur & ismember(used_inds(cur_inds(1)),usable_inds)
        cur_inf = fin_tot_corr(cur_inds);
        inferred_fix_avg(ii) = median(fin_tot_corr(cur_inds));
        inferred_drift(cur_inds) = cur_inf - inferred_fix_avg(ii);
        
        measured_fix_avgL(ii) = median(measured_seqL(cur_inds));
        measured_fix_avgR(ii) = median(measured_seqR(cur_inds));
        measured_driftL(cur_inds) = measured_seqL(cur_inds) - measured_fix_avgL(ii);
        measured_driftR(cur_inds) = measured_seqR(cur_inds) - measured_fix_avgR(ii);
    end
end

u = find(~isnan(measured_driftL) & ~isnan(inferred_drift));
[drift_corrs,drif_pvals] = corr([measured_driftL(u)' measured_driftR(u)'],inferred_drift(u)','type','spearman');
u = find(~isnan(measured_fix_avgL) & ~isnan(inferred_fix_avg));
[fix_corrs,fix_pvals] = corr([measured_fix_avgL(u) measured_fix_avgR(u)],inferred_fix_avg(u),'type','spearman');
[tot_corrs,tot_pvals] = corr([measured_seqL(usable_inds) measured_seqR(usable_inds)],fin_tot_corr(usable_inds)','type','spearman');
[tot_corrs_pear,tot_pvals_pear] = corr([measured_seqL(usable_inds) measured_seqR(usable_inds)],fin_tot_corr(usable_inds)','type','pearson');

%% FOR SACCADES
use_trials = unique(all_trialvec(used_inds));
nuse_trials = length(use_trials);
trial_blocknums = nan(nuse_trials,1);
for tt = 1:nuse_trials
    cur_used_inds = find(all_trialvec(used_inds) == use_trials(tt));
    trial_blocknums(tt) = unique(all_blockvec(used_inds(cur_used_inds)));
end

saccade_blocks = cur_block_set(all_blockvec(used_inds(saccade_start_inds)));
trial_blocks = cur_block_set(trial_blocknums);

saccade_blocks = cur_block_set(all_blockvec(used_inds(saccade_start_inds)));
use_ss = find(saccade_start_inds > 5 & saccade_start_inds < (NT-10) & ~ismember(saccade_blocks',ignore_blocks));

inferred_pre_pos = fin_tot_corr(saccade_start_inds(use_ss)-2);
inferred_post_pos = fin_tot_corr(saccade_stop_inds(use_ss) + 2);
inferred_delta_pos = (inferred_post_pos - inferred_pre_pos);
inferred_pre_pos_fo = fin_fix_corr(saccade_start_inds(use_ss)-2);
inferred_post_pos_fo = fin_fix_corr(saccade_stop_inds(use_ss) + 2);
inferred_delta_pos_fo = (inferred_post_pos_fo - inferred_pre_pos_fo);

measured_pre_pos = corrected_eye_vals_interp(used_inds(saccade_start_inds(use_ss)-2),[2 4]);
measured_post_pos = corrected_eye_vals_interp(used_inds(saccade_stop_inds(use_ss)+2),[2 4]);
measured_delta_pos = (measured_post_pos - measured_pre_pos);
measured_delta_posA = mean(measured_delta_pos,2);

u = find(~isnan(measured_delta_posA) & ~isnan(inferred_delta_pos'));
[sac_corrs,sac_pvals] = corr(measured_delta_pos(u,:),inferred_delta_pos(u)','type','spearman');

u_micro = find(~isnan(measured_delta_posA) & ~isnan(inferred_delta_pos') & use_micros);
[msac_corrs,msac_pvals] = corr(measured_delta_pos(u_micro),inferred_delta_pos(u_micro)','type','spearman');

u = find(~isnan(measured_delta_posA) & ~isnan(inferred_delta_pos') & use_nonmicros);
[nmsac_corrs,nmsac_pvals] = corr(measured_delta_pos(u),inferred_delta_pos(u)','type','spearman');


%% QUANTIFY LL IMPROVEMENTS
sus = tr_set(all_mod_SU(tr_set) > 0);
mus = tr_set(all_mod_SU(tr_set) == 0);

it_LLimp_LOO(1,:) = it_LLimp(1,:);
it_xvLLimp_LOO(1,:) = it_xvLLimp(1,:);

fix_LL_imp = bsxfun(@rdivide,it_LLimp,it_LLimp(1,:));
fix_xvLL_imp = bsxfun(@rdivide,it_xvLLimp,it_xvLLimp(1,:));
full_LL_imp = bsxfun(@rdivide,dit_LLimp,it_LLimp(1,:));
full_xvLL_imp = bsxfun(@rdivide,dit_xvLLimp,it_xvLLimp(1,:));

fix_LL_imp_LOO = bsxfun(@rdivide,it_LLimp_LOO(:,sus),it_LLimp_LOO(1,sus));
full_LL_imp_LOO = bsxfun(@rdivide,dit_LLimp_LOO(:,sus),it_LLimp_LOO(1,sus));

% figure
% plot(it_xvLLimp(1,mus),it_xvLLimp(end,mus),'k.','markersize',8)
% hold on
% plot(it_xvLLimp(1,sus),it_xvLLimp(end,sus),'r.','markersize',8)
% xlim([0 0.6]); ylim([0 0.6])
% line([0 0.6],[0 0.6],'color','k')
% xlabel('Initial xvLL (bits/spk)','fontsize',12);
% ylabel('Final xvLL (bits/spk)','fontsize',12);
% legend('MU','SU','Location','Southeast');
% box off
% set(gca,'fontname','arial','fontsize',10);
% fillPage(gcf,'papersize',[4 4]);
%
% figure
% plot(dit_xvLLimp(1,mus),dit_xvLLimp(end,mus),'k.','markersize',8)
% hold on
% plot(dit_xvLLimp(1,sus),dit_xvLLimp(end,sus),'r.','markersize',8)
% xlim([0 0.75]); ylim([0 0.75])
% line([0 0.75],[0 0.75],'color','k')
% xlabel('Initial xvLL (bits/spk)','fontsize',12);
% ylabel('Final xvLL (bits/spk)','fontsize',12);
% legend('MU','SU','Location','Southeast');
% box off
% set(gca,'fontname','arial','fontsize',10);
% fillPage(gcf,'papersize',[4 4]);

figure;hold on
errorbar(1:(n_fix_inf_it+1),nanmean(fix_LL_imp(:,mus),2),nanstd(fix_LL_imp(:,mus),[],2)/sqrt(length(mus)));
errorbar(1:(n_fix_inf_it+1),nanmean(fix_LL_imp(:,sus),2),nanstd(fix_LL_imp(:,sus),[],2)/sqrt(length(sus)),'r');

errorbar((n_fix_inf_it+2):(n_fix_inf_it+n_drift_inf_it+2),nanmean(full_LL_imp(:,mus),2),nanstd(full_LL_imp(:,mus),[],2)/sqrt(length(mus)));
errorbar((n_fix_inf_it+2):(n_fix_inf_it+n_drift_inf_it+2),nanmean(full_LL_imp(:,sus),2),nanstd(full_LL_imp(:,sus),[],2)/sqrt(length(sus)),'r');
xlabel('Iterations','fontsize',12);
ylabel('LL ratio','fontsize',12);
box off
set(gca,'fontname','arial','fontsize',10);
fillPage(gcf,'papersize',[4 4]);

figure;hold on
errorbar(1:(n_fix_inf_it+1),nanmean(fix_LL_imp(:,sus),2),nanstd(fix_LL_imp(:,sus),[],2)/sqrt(length(sus)),'k');
errorbar((n_fix_inf_it+2):(n_fix_inf_it+n_drift_inf_it+2),nanmean(full_LL_imp(:,sus),2),nanstd(full_LL_imp(:,sus),[],2)/sqrt(length(sus)),'k');
errorbar(1:(n_fix_inf_it+1),nanmean(fix_LL_imp_LOO,2),nanstd(fix_LL_imp_LOO,[],2)/sqrt(length(sus)),'r');
errorbar((n_fix_inf_it+2):(n_fix_inf_it+n_drift_inf_it+2),nanmean(full_LL_imp_LOO,2),nanstd(full_LL_imp_LOO,[],2)/sqrt(length(sus)),'r');
xlabel('Iterations','fontsize',12);
ylabel('LL ratio','fontsize',12);
box off
set(gca,'fontname','arial','fontsize',10);
fillPage(gcf,'papersize',[4 4]);

%%
close all
n_trials = length(unique(all_trialvec));
for tt = 1:n_trials
    % for tt = [96 137 154 179 376 409]
    uu = find(all_trialvec(used_inds) == tt);
    if ~isempty(uu)
        bt = all_t_axis(used_inds(uu(1)));
        et = all_t_axis(used_inds(uu(end)));
        dur = et-bt;
        if dur > 3.5
            hold on
%             h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_fix_corr(uu),fin_fix_std(uu),{'color','m'});
            h2=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_tot_corr(uu),fin_tot_std(uu),{'color','k'});
            h3=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),2),'r','linewidth',2);
            h4=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),4),'g','linewidth',2);
            h5=plot(all_t_axis(used_inds(uu))-bt,0.5*corrected_eye_vals_interp(used_inds(uu),2) +...
                0.5*corrected_eye_vals_interp(used_inds(uu),4),'c','linewidth',2);
            if use_measured_pos==1
                plot(all_t_axis(used_inds(uu))-bt,init_eyepos(used_inds(uu)),'c','linewidth',2)
            end
            %             plot(all_t_axis(used_inds(uu))-bt,nanmean(Robs_mat(uu,:),2)/5,'k');
            
            %             legend([h1.mainLine h2.mainLine h3 h4],{'Fixation corrections','Drift corrections','Left eye','Right eye'})
            xlim([0 dur]);
            ylim([-0.5 0.5]);
            xlabel('Time (s)','fontsize',10);
            ylabel('Orthoganol position (deg)','fontsize',10);
            title(sprintf('Trial %d',tt));
            set(gca,'fontsize',8,'fontname','arial');
            fillPage(gcf,'papersize',[8 5]);
            pause
            clf
        end
    end
end

%%
close all
f1 = figure();
f2 = figure();
% for ss = 1:length(tr_set)
sbeg = find(all_mod_SU(tr_set) > 0,1);
for ss = sbeg:length(tr_set)
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

