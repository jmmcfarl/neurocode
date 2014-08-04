clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

Expt_num = 86;
Expt_name = sprintf('G0%d',Expt_num);
data_dir = ['~/Data/bruce/' Expt_name];
cd(data_dir);

bar_ori = 0;

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

anal_name_PF = 'monoc_eyecorr_hbar_Fin2_loose';
anal_name_ET = 'monoc_eyecorr_hbar_Fin2ET_loose';

recompute_init_mods = 0;
use_measured_pos = 0;
use_sac_kerns = 1;

%dont fit stim models using these blocks
if Expt_num == 86
ignore_blocks = [16 17 28 30]; %G086
else
    ignore_blocks = [];
end
%%
xv_frac = 0.2;

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

%% LOAD OVERALL SU DATA
load ~/Analysis/bruce/summary_analysis/su_data.mat

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

if exist('./fin_aclust_data.mat','file')
    load('./fin_aclust_data.mat');
    [n_aclust_expts,n_aclust_probes] = size(autoclust);
else
    disp('No fin_aclust_data found.');
    autoclust = [];
    n_aclust_expts = 0; n_aclust_probes = 0;
end

%LIST OF BLOCKS TO EXCLUDE
if strcmp(Expt_name,'G086')
    %     cur_block_set(cur_block_set == 17) = []; %fails to align stim rc files.
    %     cur_block_set(cur_block_set == 28) = []; %clearly problem with the stimulus but don't know what it is!
end
if strcmp(Expt_name,'G087')
    cur_block_set(cur_block_set == 15) = []; %only 6 trials and causes problems
end
if strcmp(Expt_name,'G087')
    cur_block_set(cur_block_set == 15) = []; %only 6 trials and causes problems
end
if strcmp(Expt_name,'G093')
    cur_block_set(cur_block_set ==  28) = []; %only 6 trials and causes problems
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
cur_expt_id = find(su_data.expt_nums == Expt_num);
su_probes = find(su_data.is_su(cur_expt_id,:));
mua_probes = setdiff(1:96,su_probes); %probes with ONLY MU
aclust_probenums = [autoclust(cur_block_set(1),:).probenum];
autoclust = autoclust(:,ismember(aclust_probenums,su_probes));

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
all_spk_times = cell(96,1);
all_clust_ids = cell(length(su_probes),1);
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
    fname = sprintf('Expt%dClusterTimesDetails.mat',cur_block);
    load(fname);
    for cc = 1:96
        all_spk_times{cc} = cat(1,all_spk_times{cc},ClusterDetails{cc}.t');
        cur_ind = find(su_probes == cc);
        if ~isempty(cur_ind)
            all_clust_ids{cur_ind} = cat(1,all_clust_ids{cur_ind},autoclust(cur_block_set(ee),cur_ind).idx(:));
        end
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
%             cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            %             bar_Xmat = create_time_embedding(cur_stim_mat,stim_params);
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times];
            all_t_axis = [all_t_axis; cur_t_axis];
            %             all_Xmat = [all_Xmat; bar_Xmat];
%             all_stim_mat = [all_stim_mat; cur_stim_mat];
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
fix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
end

%% COMPUTE INITIAL EYE SEQ

init_eyepos = corrected_eye_vals_interp(:,2);
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

usable_inds = find(~ismember(all_blockvec(used_inds),ignore_blocks));

%%
cd(anal_dir);
load(anal_name_PF,'drift*','it_*','et_params','dit_*');

PF_fin_fix_corr = nan(NT,1);
PF_fin_fix_std = nan(NT,1);
PF_fin_fix_corr(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
PF_fin_fix_corr = interp1(find(~isnan(fix_ids)),PF_fin_fix_corr(~isnan(fix_ids)),1:NT);
PF_fin_fix_std(~isnan(fix_ids)) = it_fix_post_std(end,fix_ids(~isnan(fix_ids)));
PF_fin_fix_std = interp1(find(~isnan(fix_ids)),PF_fin_fix_std(~isnan(fix_ids)),1:NT);

PF_fin_fix_corr = PF_fin_fix_corr*sp_dx;
PF_fin_fix_std = PF_fin_fix_std*sp_dx;

PF_fin_drift_corr = drift_post_mean(end,:)*sp_dx;
PF_fin_drift_std = drift_post_std(end,:)*sp_dx;
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > sac_shift
        PF_fin_drift_corr(cur_inds(1:end-sac_shift+1)) = PF_fin_drift_corr(cur_inds(sac_shift:end));
        PF_fin_drift_std(cur_inds(1:end-sac_shift+1)) = PF_fin_drift_std(cur_inds(sac_shift:end));
    end
end

PF_fin_tot_corr = PF_fin_fix_corr + PF_fin_drift_corr;
PF_fin_tot_std = sqrt(PF_fin_fix_std.^2 + PF_fin_drift_std.^2);

PF_init_LLimp = it_LLimp(1,:);
PF_fin_LLimp = dit_LLimp(end,:);
%%
load(anal_name_ET,'drift*','it_*','et_params','dit_*');

ET_fin_fix_corr = nan(NT,1);
ET_fin_fix_std = nan(NT,1);
ET_fin_fix_corr(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
ET_fin_fix_corr = interp1(find(~isnan(fix_ids)),ET_fin_fix_corr(~isnan(fix_ids)),1:NT);
ET_fin_fix_std(~isnan(fix_ids)) = it_fix_post_std(end,fix_ids(~isnan(fix_ids)));
ET_fin_fix_std = interp1(find(~isnan(fix_ids)),ET_fin_fix_std(~isnan(fix_ids)),1:NT);

ET_fin_fix_corr = ET_fin_fix_corr*sp_dx;
ET_fin_fix_std = ET_fin_fix_std*sp_dx;

ET_fin_drift_corr = drift_post_mean(end,:)*sp_dx;
ET_fin_drift_std = drift_post_std(end,:)*sp_dx;
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > sac_shift
        ET_fin_drift_corr(cur_inds(1:end-sac_shift+1)) = ET_fin_drift_corr(cur_inds(sac_shift:end));
        ET_fin_drift_std(cur_inds(1:end-sac_shift+1)) = ET_fin_drift_std(cur_inds(sac_shift:end));
    end
end

ET_fin_tot_corr = ET_fin_fix_corr + ET_fin_drift_corr;
ET_fin_tot_std = sqrt(ET_fin_fix_std.^2 + ET_fin_drift_std.^2);

ET_fin_tot_corr_raw = ET_fin_tot_corr;
ET_fin_tot_corr = ET_fin_tot_corr + init_eyepos(used_inds)';

ET_init_LLimp = it_LLimp(1,:);
ET_fin_LLimp = dit_LLimp(end,:);

%%
measured_seqL = corrected_eye_vals_interp(used_inds,2);
measured_seqR = corrected_eye_vals_interp(used_inds,4);

min_fix_dur = 0.1;
PF_inferred_drift = nan(size(ET_fin_tot_corr));
ET_inferred_drift = nan(size(ET_fin_tot_corr));
measured_driftL = nan(size(ET_fin_tot_corr));
measured_driftR = nan(size(ET_fin_tot_corr));
PF_inferred_fix_avg = nan(n_fixs,1);
ET_inferred_fix_avg = nan(n_fixs,1);
measured_fix_avgL = nan(n_fixs,1);
measured_fix_avgR = nan(n_fixs,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds)*dt >= min_fix_dur & ismember(used_inds(cur_inds(1)),usable_inds)
        PF_cur_inf = PF_fin_tot_corr(cur_inds);
        ET_cur_inf = ET_fin_tot_corr(cur_inds);
        PF_inferred_fix_avg(ii) = median(PF_fin_tot_corr(cur_inds));
        ET_inferred_fix_avg(ii) = median(ET_fin_tot_corr(cur_inds));
        PF_inferred_drift(cur_inds) = PF_cur_inf - PF_inferred_fix_avg(ii);
        ET_inferred_drift(cur_inds) = ET_cur_inf - ET_inferred_fix_avg(ii);
        
        measured_fix_avgL(ii) = median(measured_seqL(cur_inds));
        measured_fix_avgR(ii) = median(measured_seqR(cur_inds));
        measured_driftL(cur_inds) = measured_seqL(cur_inds) - measured_fix_avgL(ii);
        measured_driftR(cur_inds) = measured_seqR(cur_inds) - measured_fix_avgR(ii);
     end
end

u = find(~isnan(measured_driftL) & ~isnan(ET_inferred_drift));
[ETdrift_corrs,drif_pvals] = corr([measured_driftL(u)' measured_driftR(u)'],ET_inferred_drift(u)','type','spearman');
[PFdrift_corrs,drif_pvals] = corr([measured_driftL(u)' measured_driftR(u)'],PF_inferred_drift(u)','type','spearman');
[ET_PFdrift_corrs,drif_pvals] = corr(ET_inferred_drift(u)',PF_inferred_drift(u)','type','spearman');
u = find(~isnan(measured_fix_avgL) & ~isnan(ET_inferred_fix_avg));
[ET_fix_corrs,fix_pvals] = corr([measured_fix_avgL(u) measured_fix_avgR(u)],ET_inferred_fix_avg(u),'type','spearman');
[PF_fix_corrs,fix_pvals] = corr([measured_fix_avgL(u) measured_fix_avgR(u)],PF_inferred_fix_avg(u),'type','spearman');
[ET_PF_fix_corrs,fix_pvals] = corr(ET_inferred_fix_avg(u),PF_inferred_fix_avg(u),'type','spearman');
[ET_tot_corrs,tot_pvals] = corr([measured_seqL(usable_inds) measured_seqR(usable_inds)],ET_fin_tot_corr(usable_inds)','type','spearman');
[PF_tot_corrs,tot_pvals] = corr([measured_seqL(usable_inds) measured_seqR(usable_inds)],PF_fin_tot_corr(usable_inds)','type','spearman');
[ET_PF_tot_corrs,tot_pvals] = corr(ET_fin_tot_corr(usable_inds)',PF_fin_tot_corr(usable_inds)','type','spearman');

%% FOR SACCADES
use_ss = find(saccade_start_inds > 5 & saccade_start_inds < (NT-10) & ~ismember(saccade_blocks',ignore_blocks));
use_micros = (ismember(use_ss,micro_sacs));
use_nonmicros = (~ismember(use_ss,micro_sacs));

saccade_blocks = cur_block_set(all_blockvec(used_inds(saccade_start_inds)));

inferred_pre_pos = ET_fin_tot_corr(saccade_start_inds(use_ss)-2);
inferred_post_pos = ET_fin_tot_corr(saccade_start_inds(use_ss) + 5);
ET_inferred_delta_pos = (inferred_post_pos - inferred_pre_pos);

inferred_pre_pos = PF_fin_tot_corr(saccade_start_inds(use_ss)-2);
inferred_post_pos = PF_fin_tot_corr(saccade_start_inds(use_ss) + 5);
PF_inferred_delta_pos = (inferred_post_pos - inferred_pre_pos)';

u = find(~isnan(ET_inferred_delta_pos) & ~isnan(PF_inferred_delta_pos) & use_nonmicros);
[nmsac_corrs,nmsac_pvals] = corr(ET_inferred_delta_pos(u),PF_inferred_delta_pos(u),'type','spearman');

u = find(~isnan(ET_inferred_delta_pos) & ~isnan(PF_inferred_delta_pos) & use_micros);
[msac_corrs,msac_pvals] = corr(ET_inferred_delta_pos(u),PF_inferred_delta_pos(u),'type','spearman');

%% MICROSACCADE AMP COMPARISON
cd ~/Analysis/bruce/summary_analysis/eyetrack_figs/

sqrt_scale = true;
n_bins = 40;

close all

figure
if sqrt_scale
    DensityPlot_jmm(ET_inferred_delta_pos,PF_inferred_delta_pos,'ynormal','xrange',[-0.5 0.5],'yrange',[-0.5 0.5],'sd',[3 3],'sqrtsc');
else
    DensityPlot_jmm(ET_inferred_delta_pos,PF_inferred_delta_pos,'ynormal','xrange',[-0.5 0.5],'yrange',[-0.5 0.5],'sd',[3 3]);
end
line([-0.5 0.5],[-0.5 0.5],'color','w');
xlabel('Inferred amplitude (deg)','fontsize',12);
ylabel('Measured amplitude (deg)','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
set(gca,'xtick',-0.5:0.1:0.5,'ytick',-0.5:0.1:0.5);
fillPage(gcf,'papersize',[5 5]);
fname = sprintf('Msac_dens_ori%d',bar_ori);
axis square

b_ax = linspace(-0.5,0.5,n_bins);
nx_m = histc(ET_inferred_delta_pos(use_micros),b_ax);
nx_m = nx_m/sum(nx_m);
nx_nm = histc(ET_inferred_delta_pos(use_nonmicros),b_ax);
nx_nm = nx_nm/sum(nx_nm);
figure; hold on
stairs(b_ax,nx_m,'b');
% stairs(b_ax,nx_nm,'r');
yl = ylim();
line([0 0],yl,'color','k')
xlim(b_ax([1 end])); set(gca,'xtick',[],'ytick',[]);
box off;
title('Inferred position','fontsize',10);
fillPage(gcf,'papersize',[5 5]);
fname = sprintf('Msac_infmarg_ori%d',bar_ori);


nx_m = histc(PF_inferred_delta_pos(use_micros),b_ax);
nx_m = nx_m/sum(nx_m);
nx_nm = histc(PF_inferred_delta_pos(use_nonmicros),b_ax);
nx_nm = nx_nm/sum(nx_nm);
figure; hold on
stairs(b_ax,nx_m,'b');
% stairs(b_ax,nx_nm,'r');
yl = ylim();
line([0 0],yl,'color','k')
xlim(b_ax([1 end])); set(gca,'xtick',[],'ytick',[]);
box off;
title('Measured position','fontsize',10);
fname = sprintf('Msac_measmarg_ori%d',bar_ori);


nx_m = histc(ET_inferred_delta_pos(use_micros) - PF_inferred_delta_pos(use_micros),b_ax);
nx_m = nx_m/sum(nx_m);
nx_nm = histc(ET_inferred_delta_pos(use_nonmicros) - PF_inferred_delta_pos(use_nonmicros),b_ax);
nx_nm = nx_nm/sum(nx_nm);
figure; hold on
stairs(b_ax,nx_m,'b');
% stairs(b_ax,nx_nm,'r');
yl = ylim();
line([0 0],yl,'color','k')
xlim(b_ax([1 end])); set(gca,'xtick',[],'ytick',[]);
box off;
title('Difference','fontsize',10);
fname = sprintf('Msac_diff_ori%d',bar_ori);

%% FIXATION BASED COMPARISON
close all

figure
if sqrt_scale
    DensityPlot_jmm(inferred_fix_avg_ET,inferred_fix_avg_PF,'ynormal','xrange',[-0.5 0.5],'yrange',[-0.5 0.5],'sd',[3 3],'sqrtsc');
else
    DensityPlot_jmm(inferred_fix_avg_ET,inferred_fix_avg_PF,'ynormal','xrange',[-0.5 0.5],'yrange',[-0.5 0.5],'sd',[3 3]);    
end
line([-0.5 0.5],[-0.5 0.5],'color','w');
xlabel('Inferred position (deg)','fontsize',12);
ylabel('Measured position (deg)','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
set(gca,'xtick',-0.5:0.1:0.5,'ytick',-0.5:0.1:0.5);
axis square
fillPage(gcf,'papersize',[5 5]);
fname = sprintf('Fix_dens_ori%d',bar_ori);

b_ax = linspace(-0.5,0.5,n_bins);
nx = histc(inferred_fix_avg_ET,b_ax);
nx = nx/sum(nx);
figure; hold on
stairs(b_ax,nx,'b');
yl = ylim();
line([0 0],yl,'color','k')
xlim(b_ax([1 end])); set(gca,'xtick',[],'ytick',[]);
box off;
title('ET Init','fontsize',10);
fname = sprintf('Fix_infmarg_ori%d',bar_ori);

nx = histc(inferred_fix_avg_PF,b_ax);
nx = nx/sum(nx);
figure; hold on
stairs(b_ax,nx,'b');
yl = ylim();
line([0 0],yl,'color','k')
xlim(b_ax([1 end])); set(gca,'xtick',[],'ytick',[]);
box off;
title('PF Init','fontsize',10);
fname = sprintf('Fix_measmarg_ori%d',bar_ori);

nx = histc(inferred_fix_avg_ET - inferred_fix_avg_PF,b_ax);
nx = nx/sum(nx);
figure; hold on
stairs(b_ax,nx,'b');
yl = ylim();
line([0 0],yl,'color','k')
xlim(b_ax([1 end])); set(gca,'xtick',[],'ytick',[]);
box off;
title('Difference','fontsize',10);
fname = sprintf('Fix_diff_ori%d',bar_ori);

%% DRIFT COMPARISON

close all
min_damp = 0.005;

figure
uset = find(abs(inferred_drift_ET) >= min_damp & abs(inferred_drift_PF) >= min_damp);
[h,D] = DensityPlot_jmm(inferred_drift_ET(uset),inferred_drift_PF(uset),'ynormal','xrange',[-0.1 0.1],'yrange',[-0.1 0.1],'sd',[3 3]);
line([-0.1 0.1],[-0.1 0.1],'color','w');
to_nan = abs(D.x) < min_damp | abs(D.y) < min_damp;
D.z(to_nan) = nan;
if sqrt_scale
    imagescnan(D.x,D.y,sqrt(D.z),'NanColor','w'); set(gca,'ydir','normal');
else
    imagescnan(D.x,D.y,D.z,'NanColor','w'); set(gca,'ydir','normal');
end
line([-0.1 0.1],[-0.1 0.1],'color','w');
xlabel('Inferred position (deg)','fontsize',10);
ylabel('Measured position (deg)','fontsize',10);
set(gca,'fontsize',10,'fontname','arial');
set(gca,'xtick',-0.1:0.02:0.1,'ytick',-0.1:0.02:0.1);
axis square
fillPage(gcf,'papersize',[5 5]);
fname = sprintf('Drift_dens_ori%d',bar_ori);

b_ax = linspace(-0.1,0.1,n_bins);
nx = histc(inferred_drift_ET,b_ax);
nx = nx/sum(nx);
figure; hold on
stairs(b_ax,nx,'b');
yl = ylim();
line([0 0],yl,'color','k')
xlim(b_ax([1 end])); set(gca,'xtick',[],'ytick',[]);
box off;
title('ET Init','fontsize',10);
fname = sprintf('Drift_infmarg_ori%d',bar_ori);

nx = histc(inferred_drift_PF,b_ax);
nx = nx/sum(nx);
figure; hold on
stairs(b_ax,nx,'b');
yl = ylim();
line([0 0],yl,'color','k')
xlim(b_ax([1 end])); set(gca,'xtick',[],'ytick',[]);
box off;
title('PF Init','fontsize',10);
fname = sprintf('Drift_measmarg_ori%d',bar_ori);

nx = histc(inferred_drift_ET - inferred_drift_PF,b_ax);
nx = nx/sum(nx);
figure; hold on
stairs(b_ax,nx,'b');
yl = ylim();
line([0 0],yl,'color','k')
xlim(b_ax([1 end])); set(gca,'xtick',[],'ytick',[]);
box off;
title('Difference','fontsize',10);
fname = sprintf('Drift_diff_ori%d',bar_ori);

%% TOTAL COMPARISON
close all

figure
if sqrt_scale
    DensityPlot_jmm(ET_fin_tot_corr,PF_fin_tot_corr','ynormal','xrange',[-0.5 0.5],'yrange',[-0.5 0.5],'sd',[3 3],'sqrtsc');
else
    DensityPlot_jmm(ET_fin_tot_corr,PF_fin_tot_corr','ynormal','xrange',[-0.5 0.5],'yrange',[-0.5 0.5],'sd',[3 3]);    
end
line([-0.5 0.5],[-0.5 0.5],'color','w');
xlabel('Inferred position (deg)','fontsize',12);
ylabel('Measured position (deg)','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
set(gca,'xtick',-0.5:0.1:0.5,'ytick',-0.5:0.1:0.5);
axis square
fillPage(gcf,'papersize',[5 5]);
fname = sprintf('Fix_dens_ori%d',bar_ori);

b_ax = linspace(-0.5,0.5,n_bins);
nx = histc(ET_fin_tot_corr,b_ax);
nx = nx/sum(nx);
figure; hold on
stairs(b_ax,nx,'b');
yl = ylim();
line([0 0],yl,'color','k')
xlim(b_ax([1 end])); set(gca,'xtick',[],'ytick',[]);
box off;
title('ET Init','fontsize',10);
fname = sprintf('Fix_infmarg_ori%d',bar_ori);

nx = histc(PF_fin_tot_corr',b_ax);
nx = nx/sum(nx);
figure; hold on
stairs(b_ax,nx,'b');
yl = ylim();
line([0 0],yl,'color','k')
xlim(b_ax([1 end])); set(gca,'xtick',[],'ytick',[]);
box off;
title('PF Init','fontsize',10);
fname = sprintf('Fix_measmarg_ori%d',bar_ori);

nx = histc(ET_fin_tot_corr - PF_fin_tot_corr',b_ax);
nx = nx/sum(nx);
figure; hold on
stairs(b_ax,nx,'b');
yl = ylim();
line([0 0],yl,'color','k')
xlim(b_ax([1 end])); set(gca,'xtick',[],'ytick',[]);
box off;
title('Difference','fontsize',10);
fname = sprintf('Fix_diff_ori%d',bar_ori);

%%
n_bins = 60;
b_ax = linspace(-0.5,0.5,n_bins);
nx = histc(measured_seqL,b_ax);
nx = nx/max(nx);
ny = histc(PF_fin_tot_corr-ET_fin_tot_corr',b_ax);
ny = ny/max(ny);
figure; hold on
stairs(b_ax,nx,'b','linewidth',2);
stairs(b_ax,ny,'r','linewidth',2);
legend('Before','After');
xlabel('Position difference (deg)','fontsize',10);
ylabel('Relative frequency (AU)','fontsize',10);
set(gca,'fontname','arial','fontsize',8);
set(gca,'xtick',[-0.5:0.1:0.5]);
fillPage(gcf,'papersize',[4 4]);
fname = 'PF_ET_init_diff';

%%

close all
n_trials = length(unique(all_trialvec));
%328
for tt = 250:n_trials
    uu = find(all_trialvec(used_inds) == tt);
    if ~isempty(uu)
        bt = all_t_axis(used_inds(uu(1)));
        et = all_t_axis(used_inds(uu(end)));
        dur = et-bt;
        if dur > 3.5
            hold on
            h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,ET_fin_tot_corr(uu),ET_fin_tot_std(uu),{'color','b'});
            h2=shadedErrorBar(all_t_axis(used_inds(uu))-bt,PF_fin_tot_corr(uu),PF_fin_tot_std(uu),{'color','k'});
%             h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,ET_fin_fix_corr(uu),ET_fin_fix_std(uu),{'color','b'});
%             h2=shadedErrorBar(all_t_axis(used_inds(uu))-bt,PF_fin_fix_corr(uu),PF_fin_fix_std(uu),{'color','k'});
            if bar_ori == 0
                h3=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),2),'r','linewidth',2);
%                 h4=plot(all_t_axis(used_inds(uu))-bt,corrected_interp_eyevals(used_inds(uu),4)-median(corrected_interp_eyevals(used_inds(uu),4)),'color',[0.2 0.8 0.2],'linewidth',2);
            else
                h3=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),1),'r','linewidth',2);
%                 h4=plot(all_t_axis(used_inds(uu))-bt,corrected_interp_eyevals(used_inds(uu),3)-median(corrected_interp_eyevals(used_inds(uu),3)),'color',[0.2 0.8 0.2],'linewidth',2);
            end
            
                        plot(all_t_axis(used_inds(uu))-bt,init_eyepos(used_inds(uu)),'c','linewidth',2);
            
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