clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

%run on 86 [0, 90];
Expt_num = 86;
Expt_name = sprintf('G%.3d',Expt_num);
data_dir = ['~/Data/bruce/' Expt_name];
fig_dir = '/home/james/Analysis/bruce/ET_final/';
cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final/'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

bar_ori = 0;
use_inf_sim = false;

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

poss_sim_gains = [0.25 0.5 1 2 4 8];
% poss_sim_gains = [3 1 0.25 0.5 3];
n_reps_per_gain = 1;


flen = 12;
use_nPix = 16;

n_fix_inf_it = 5; %4
n_drift_inf_it = 2; %2

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
            all_stim_times = [all_stim_times; cur_stim_times];
            all_t_axis = [all_t_axis; cur_t_axis];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
    end
    trial_cnt = trial_cnt + n_trials;
end

%% DEFINE DATA USED FOR ANALYSIS
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer);
NT = length(used_inds);

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

%%
cd(anal_dir)
load(true_eyedata_name,'drift_*','it_fix_*');

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
best_est_corr = fin_fix_corr + fin_drift_corr;

%%
for sg = 1:length(poss_sim_gains)
    
    sim_gain = poss_sim_gains(sg);
    
    %%
    if ~use_inf_sim
        sim_eyepos = corrected_eye_vals_interp(:,2);
        sim_eyepos = sim_eyepos*sim_gain;
        sim_eyepos(isnan(sim_eyepos)) = 0;
        
        %smooth out fast transients in eye signal
        eye_smooth_sig = round(0.025/dt);
        interp_inds = [];
        for ii = 1:n_fixs
            cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
            if length(cur_inds) > eye_smooth_sig*5;
                sim_eyepos(used_inds(cur_inds)) = jmm_smooth_1d_cor(sim_eyepos(used_inds(cur_inds)),eye_smooth_sig,2);
            end
            interp_inds = [interp_inds; cur_inds'];
        end
        interp_inds = unique(interp_inds);
        sim_eyepos(used_inds) = interp1(used_inds(interp_inds),sim_eyepos(used_inds(interp_inds)),used_inds);
    
%         for ii = 1:length(cur_block_set)
%             cur_inds = find(all_blockvec == ii);
%             sim_eyepos(cur_inds,:) = bsxfun(@minus,sim_eyepos(cur_inds,:),median(sim_eyepos(cur_inds,:)));
%         end

        %subtract off within-trial median
        for ii = 1:length(trial_start_inds)
            cur_inds = trial_start_inds(ii):trial_end_inds(ii);
            sim_eyepos(used_inds(cur_inds),:) = bsxfun(@minus,sim_eyepos(used_inds(cur_inds),:),median(sim_eyepos(used_inds(cur_inds),:)));
        end

        %maximum initial corrections
        if sim_gain >= 4
            max_sim_pos = max_shift*2*sp_dx;
        else
            max_sim_pos = max_shift*sp_dx;
        end
        sim_eyepos(sim_eyepos > max_sim_pos) = max_sim_pos;
        sim_eyepos(sim_eyepos < - max_sim_pos) = -max_sim_pos;
        
    else
        sim_eyepos = zeros(length(all_t_axis),1);
        sim_eyepos(used_inds) = best_est_corr;
        sim_eyepos = sim_eyepos*sim_gain;
        
        %maximum initial corrections
        if sim_gain >= 4
            max_sim_pos = max_shift*2*sp_dx;
        else
            max_sim_pos = max_shift*sp_dx;
        end
        sim_eyepos(sim_eyepos > max_sim_pos) = max_sim_pos;
        sim_eyepos(sim_eyepos < - max_sim_pos) = -max_sim_pos;
        sim_eyepos_us_rnd = round(sim_eyepos/sp_dx);
    end
    
    all_sim_eyepos(sg,:) = sim_eyepos(used_inds,:);
    
    %%
    for rpt = 1:n_reps_per_gain
        if use_inf_sim
            data_name = [anal_dir sprintf('ETsim_hinfEP_it%d_gain%d',rpt,sg)];
        else
            data_name = [anal_dir sprintf('ETsim_it%d_gain%d',rpt,sg)];
%             data_name = [anal_dir sprintf('ETsim2_it%d_gain%d',rpt,sg)];
%             data_name = [anal_dir sprintf('ETsim3_it%d_gain%d',rpt,sg)];
        end
        fprintf('Loading %s\n',data_name);
        load(data_name,'drift_*','it_fix*','it_R2','dit_R2');
        
        all_init_R2(sg,:,rpt) = it_R2(1,:);
        all_fin_R2(sg,:,rpt) = dit_R2(end,:);
        
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
        
        fin_tot_corr = fin_fix_corr + fin_drift_corr;
        fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2);
        
        all_inf_eyepos(sg,:,rpt) = fin_tot_corr;
        all_std_eyepos(sg,:,rpt) = fin_tot_std;
        cur_eyepos_err = squeeze(all_inf_eyepos(sg,:,rpt) - all_sim_eyepos(sg,:));
        all_eyepos_err(sg,:,rpt) = all_inf_eyepos(sg,:,rpt) - all_sim_eyepos(sg,:);
        all_inf_rms(sg,rpt) = robust_std_dev(cur_eyepos_err);
        all_err_std(sg,rpt) = robust_std_dev(abs(cur_eyepos_err));
        all_inf_SD(sg,rpt) = std(cur_eyepos_err);
        all_inf_prec(sg,rpt) = nanmedian(fin_tot_std);
        all_inf_prec_prc(sg,rpt,:) = prctile(fin_tot_std,[25 50 75]);
        all_inf_err_prc(sg,rpt,:) = prctile(abs(cur_eyepos_err),[25 50 75]);
%         cur_unc = 1./fin_tot_std;
%         cur_unc = cur_unc/sum(cur_unc);
%         wSD = sqrt(sum(cur_unc.*cur_eyepos_err.^2));
%         all_inf_wSD(sg,rpt) = wSD;
    end
    clear it_fix_* drift_post*
end

%% USE ROBUST SD ERR
% load('~/Analysis/bruce/ET_final/G086_hbar_mua_RFwidth.mat');
load('~/Analysis/bruce/ET_final/G086_hbar_mua_RFwidth2.mat');

uset = find(poss_sim_gains <= 4);

avg_rf_width = mean(all_RF_widths); %2sigma, avg across all units from G086, horizontal
std_rf_width = std(all_RF_widths);
prc_rf_widths = prctile(all_RF_widths,[25 50 75]);

fig_width = 3.27;
rel_height = 0.8;

all_inf_errSD_deg = bsxfun(@rdivide,all_err_std,poss_sim_gains');

all_inf_prc_deg = bsxfun(@rdivide,all_inf_err_prc,poss_sim_gains');
all_inf_prc_deg = squeeze(mean(all_inf_prc_deg,2));
avg_rf_widths_deg = avg_rf_width./poss_sim_gains(uset);
std_rf_widths_deg = std_rf_width./poss_sim_gains(uset);
prc_rf_widths_deg = bsxfun(@rdivide,prc_rf_widths,poss_sim_gains(uset)');

h=figure;
eh = errorbarxy(prc_rf_widths_deg(:,2),all_inf_prc_deg(uset,2),prc_rf_widths_deg(:,2)-prc_rf_widths_deg(:,1),prc_rf_widths_deg(:,3)-prc_rf_widths_deg(:,2),...
    all_inf_prc_deg(uset,2)-all_inf_prc_deg(uset,1),all_inf_prc_deg(uset,3)-all_inf_prc_deg(uset,2),{'k','k','k'});
% eh = errorbarxy(avg_rf_widths_deg,all_inf_prc_deg(uset,2),std_rf_widths_deg,all_inf_errSD_deg(uset),{'k','k','k'});
set(eh(1),'linewidth',1.5);
set(eh(2:end),'linewidth',0.75);
plot(prc_rf_widths_deg(:,2),all_inf_prc_deg(uset,2),'ko','linewidth',1,'markersize',4)
xlabel('Average RF width');
ylabel('Error (deg)');
xlim([0 0.6])
set(gca,'xscale','log');
xlim([.025 0.6]);

% h=figure;
% plot(avg_rf_width./poss_sim_gains(uset),all_inf_rms_deg(uset),'ko-','linewidth',1.5,'markersize',4);
% xlabel('Average RF width');
% ylabel('Error (deg)');
% yl = ylim();
% ylim([0 yl(2)]);
figufy(h);
fname = [fig_dir 'eyeErr_vs_RFwidth2.pdf'];
% exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close;

%% PRECISION VS RF WIDTH
load('~/Analysis/bruce/ET_final/G086_hbar_mua_RFwidth2.mat');

uset = find(poss_sim_gains <= 4);

avg_rf_width = mean(all_RF_widths); %2sigma, avg across all units from G086, horizontal
std_rf_width = std(all_RF_widths);
prc_rf_widths = prctile(all_RF_widths,[25 50 75]);

fig_width = 3.27;
rel_height = 0.8;

all_inf_prec_prc_deg = bsxfun(@rdivide,all_inf_prec_prc,poss_sim_gains');
all_inf_prec_prc_deg = squeeze(mean(all_inf_prec_prc_deg,2));
avg_rf_widths_deg = avg_rf_width./poss_sim_gains(uset);
std_rf_widths_deg = std_rf_width./poss_sim_gains(uset);
prc_rf_widths_deg = bsxfun(@rdivide,prc_rf_widths,poss_sim_gains(uset)');

h=figure;
eh = errorbarxy(prc_rf_widths_deg(:,2),all_inf_prec_prc_deg(uset,2),prc_rf_widths_deg(:,2)-prc_rf_widths_deg(:,1),prc_rf_widths_deg(:,3)-prc_rf_widths_deg(:,2),...
    all_inf_prec_prc_deg(uset,2)-all_inf_prec_prc_deg(uset,1),all_inf_prec_prc_deg(uset,3)-all_inf_prec_prc_deg(uset,2),{'k','k','k'});
set(eh(1),'linewidth',1.5);
set(eh(2:end),'linewidth',0.75);
plot(prc_rf_widths_deg(:,2),all_inf_prc_deg(uset,2),'ko','linewidth',1,'markersize',4)
xlabel('Average RF width');
ylabel('Error (deg)');
xlim([0 0.6])
set(gca,'xscale','log');
xlim([.025 0.6]);

% h=figure;
% plot(avg_rf_width./poss_sim_gains(uset),all_inf_rms_deg(uset),'ko-','linewidth',1.5,'markersize',4);
% xlabel('Average RF width');
% ylabel('Error (deg)');
% yl = ylim();
% ylim([0 yl(2)]);
figufy(h);
fname = [fig_dir 'eyeErr_vs_RFwidth2.pdf'];

%% CUM ERROR DISTS VS RF WIDTH

uset = find(poss_sim_gains <= 4);

n_pts = 200;
xx = linspace(0,0.1,n_pts);
% xx = linspace(0,0.5,n_pts);
ydist = nan(length(uset),size(all_eyepos_err,3),n_pts);
for cc = 1:length(uset)
    for rr = 1:size(all_eyepos_err,3)
        yobs = abs(squeeze(all_eyepos_err(cc,:,rr)));
        yobs = yobs/poss_sim_gains(cc);
        temp = histc(yobs,xx);
%         temp = cumsum(temp)/sum(temp);
        temp = cumsum(temp)/length(yobs);
        ydist(cc,rr,:) = temp;
    end
end
yobs = abs(squeeze(all_sim_eyepos(poss_sim_gains==1,:)));
yobs = yobs/poss_sim_gains(poss_sim_gains==1);
temp = histc(yobs,xx);
% temp = cumsum(temp)/sum(temp);
temp = cumsum(temp)/length(yobs);
base_err = temp;

cmap = jet(length(uset));
h1=figure; hold on
for ii = 1:length(uset)
    if poss_sim_gains(ii) == 1
    plot(xx,squeeze(nanmean(ydist(ii,:,:),2)),'color',cmap(ii,:),'linewidth',2);
    else
    plot(xx,squeeze(nanmean(ydist(ii,:,:),2)),'color',cmap(ii,:),'linewidth',1.5);
    end
end
plot(xx,base_err,'k','linewidth',2);
xlabel('Error magnitude (deg)');
ylabel('Cumulative probability');

fig_width = 3.27;
rel_height = 0.8;

% figufy(h1);
% fname = [fig_dir 'errDist_vs_RFwidth.pdf'];
% exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close;



%% 
all_R2_imp = all_fin_R2./all_init_R2;
ravg_R2_imp = squeeze(mean(all_R2_imp,3));
mean_R2_imp = mean(ravg_R2_imp,2);
std_R2_imp = std(ravg_R2_imp,[],2);
prc_R2_imp = prctile(all_R2_imp',[25 50 75])';

fig_width = 3.27;
rel_height = 0.8;

uset = find(poss_sim_gains <= 4);

h2 = figure(); hold on
% eh = errorbarxy(avg_rf_width./poss_sim_gains(uset),mean_R2_imp(uset),std_rf_width./poss_sim_gains(uset),std_R2_imp(uset),{'k','k','k'});
eh = errorbarxy(prc_rf_widths_deg(:,2),prc_R2_imp(uset,2),prc_rf_widths_deg(:,2)-prc_rf_widths_deg(:,1),prc_rf_widths_deg(:,3)-prc_rf_widths_deg(:,2),...
    prc_R2_imp(uset,2)-prc_R2_imp(uset,1),prc_R2_imp(uset,3)-prc_R2_imp(uset,2),{'k','k','k'});
set(eh(1),'linewidth',1.5);
set(eh(2:end),'linewidth',0.75);
% plot(avg_rf_width./poss_sim_gains(uset),mean_R2_imp(uset),'ko','linewidth',1,'markersize',4)
plot(prc_rf_widths_deg(:,2),prc_R2_imp(uset,2),'ko','linewidth',1,'markersize',4)
xlabel('RF width');
ylabel('Model improvement');
xlim([0 0.6]);
set(gca,'xscale','log');
xlim([.025 0.6]);

figufy(h2);
fname = [fig_dir 'modimp_vs_RFwidth.pdf'];
exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close;

%%
cd(anal_dir);
load(true_moddata_name,'all_mod*');
load(true_eyedata_name,'et_tr_set','dit_mods');
su_inds = find(all_mod_SU(et_tr_set) > 0);
cur_su = 7; 

true_mod = dit_mods{end}(et_tr_set(su_inds(cur_su)));

sg = 5;
data_name = [anal_dir sprintf('ETsim_it%d_gain%d',1,sg)];
fprintf('Loading %s\n',data_name);
load(data_name,'drift_*','it_fix*','it_R2','dit_R2','it_mods','dit_mods');

pre_mod = it_mods{1}(su_inds(cur_su));
post_mod = dit_mods{end}(su_inds(cur_su));

flen = 12;
dt = 0.01;
use_nPix_us = use_nPix*spatial_usfac;
sp_dx = 0.0565/spatial_usfac/poss_sim_gains(sg);

lag_axis = (0:flen-1)*dt*1e3;
pix_axis = (1:use_nPix_us)*sp_dx - use_nPix_us/2*sp_dx;

Xtargs = [post_mod.mods(:).Xtarget];
n_stimfilts = sum(Xtargs == 1);
n_squared_filts = sum(Xtargs == 1) - 1;
true_filts = [true_mod.mods((Xtargs == 1)).filtK];
cor_filts = [post_mod.mods((Xtargs == 1)).filtK];
uncor_filts = [pre_mod.mods((Xtargs == 1)).filtK];
max_vals = max(abs(cor_filts));
uc_max_vals = max(abs(uncor_filts));
true_max_vals = max(abs(true_filts));

true_filts = reshape(true_filts,[flen use_nPix_us n_squared_filts+1]);
cor_filts = reshape(cor_filts,[flen use_nPix_us n_squared_filts+1]);
uncor_filts = reshape(uncor_filts,[flen use_nPix_us n_squared_filts+1]);

true_temp_profiles = squeeze(var(true_filts,[],2));
cor_temp_profiles = squeeze(var(cor_filts,[],2));
uncor_temp_profiles = squeeze(var(uncor_filts,[],2));
[~,best_lags_true] = max(true_temp_profiles);
[~,best_lags_cor] = max(cor_temp_profiles);
[~,best_lags_uncor] = max(uncor_temp_profiles);

% best_lags_cor = [5 5 5];

for i = 1:n_stimfilts
    spatial_profiles_true(:,i) = squeeze(true_filts(best_lags_true(i),:,i));
    spatial_profiles_cor(:,i) = squeeze(cor_filts(best_lags_true(i),:,i));
    spatial_profiles_uncor(:,i) = squeeze(uncor_filts(best_lags_true(i),:,i));
end

disp_pix = find(pix_axis >= -0.5 & pix_axis <= 0.5);

h=figure;
for ii = 1:n_stimfilts
    subplot(n_stimfilts,3,(ii-1)*3+1);
    imagesc(pix_axis(disp_pix),lag_axis,squeeze(true_filts(:,disp_pix,ii))); 
    caxis([-true_max_vals(ii) true_max_vals(ii)]);
%     line(pix_axis([1 end]),lag_axis(best_lags_true([ii ii])),'color','r');
    xlabel('Position (deg)');
    ylabel('Lag (ms)');
    set(gca,'xtick',[-0.1:0.1:0.1]);
    set(gca,'ydir','normal');
    %     xlim([-0.5 0.5]);

    subplot(n_stimfilts,3,(ii-1)*3+2);
    imagesc(pix_axis(disp_pix),lag_axis,squeeze(uncor_filts(:,disp_pix,ii))); 
%     caxis([-true_max_vals(ii) true_max_vals(ii)]);
    caxis([-uc_max_vals(ii) uc_max_vals(ii)]);
%     line(pix_axis([1 end]),lag_axis(best_lags_true([ii ii])),'color','r');
    xlabel('Position (deg)');
    ylabel('Lag (ms)');
    set(gca,'xtick',[-0.1:0.1:0.1]);
    set(gca,'ydir','normal');
%     xlim([-0.5 0.5]);

subplot(n_stimfilts,3,(ii-1)*3+3);
    imagesc(pix_axis(disp_pix),lag_axis,squeeze(cor_filts(:,disp_pix,ii))); 
%     caxis([-true_max_vals(ii) true_max_vals(ii)]);
    caxis([-max_vals(ii) max_vals(ii)]);
%     line(pix_axis([1 end]),lag_axis(best_lags_true([ii ii])),'color','r');
    xlabel('Position (deg)');
    ylabel('Lag (ms)');
    set(gca,'xtick',[-0.1:0.1:0.1]);
    set(gca,'ydir','normal');
    %     xlim([-0.5 0.5]);
    
end
colormap(gray);
subplot(n_stimfilts,3,1);
title('True');
subplot(n_stimfilts,3,2);
title('Uncorrected');
subplot(n_stimfilts,3,3);
title('Corrected');

%%
fig_width = 4.86; %3.27 4.86 6.83
rel_height = 0.8;

figufy(h);
fname = [fig_dir 'RFsim_exampmod.pdf'];
exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close;
