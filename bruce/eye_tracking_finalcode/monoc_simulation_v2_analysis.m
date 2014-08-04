clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');
fig_dir = ['/home/james' '/Analysis/bruce/ET_final/'];

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
% sim_data_name = ['monoc_eyecorr_hbar_simdata_g' gname]';

do_highres_est = true;
use_sac_kerns = 1;
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

%%
xv_frac = 0;

flen = 12;
use_nPix = 16;

drift_dsf = 3;

% %apply ao minimum threshold to the drift prior sigma so lattice jumps are
% %possible.
% drift_prior_sigma = max(drift_prior_sigma,0.003);

min_trial_dur = 0.75;

spatial_usfac = 2;

%%
eps = -1e3;

stim_fs = 100; %in Hz
dt = 0.01;
full_nPix = 36;
% if sim_gain == 8
%     full_nPix = full_nPix + 18;
% end

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
            
            %                 cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
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

%% DEFINE DATA USED FOR ANALYSIS
used_inds = find(all_tsince_start >= beg_buffer & (trial_dur-all_tsince_start) >= end_buffer);
used_inds_new = find(all_tsince_start >= beg_buffer_new & (trial_dur-all_tsince_start) >= end_buffer);
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
n_trials = length(trial_start_inds);
trial_ids = nan(NT,1);
for ii = 1:n_trials
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    trial_ids(cur_inds) = ii;
end

%%
cd(anal_dir);
base_simdata = 'monoc_eyecorr_hbar_simdata_g1.mat';
load(base_simdata,'sim_eyepos');
%%
poss_gains = [0.25 0.5 1 2 4];

for sg = 1:length(poss_gains)
    fprintf('Sim gain %d/%d\n',sg,length(poss_gains));
    sim_gain = poss_gains(sg);
    if sim_gain >= 1
        gname = num2str(sim_gain);
    elseif sim_gain == 0.5
        gname = 'p5';
    elseif sim_gain == 0.25
        gname = 'p25';
    else
        error('unsupported');
    end
    
    switch sim_gain
        case 4
            hr_spatial_usfac = 4;
            spatial_usfac = 2;
        case 2
            hr_spatial_usfac = 4;
            spatial_usfac = 2;
        case 1
            hr_spatial_usfac = 4;
            spatial_usfac = 2;
        case 0.5
            hr_spatial_usfac = 8;
            spatial_usfac = 4;
        case 0.25
            hr_spatial_usfac = 8;
            spatial_usfac = 4;
    end
    
    cur_simdata = ['monoc_eyecorr_hbar_simdata_g' gname];
    temp = load(cur_simdata,'sim_eyepos');
    cur_sim_eyepos(sg,:) = temp.sim_eyepos/sim_gain;
%     if sim_gain > 1
%         data_name_hres = ['monoc_eyecorr_hbar_hres_sim3_g' gname];
%         data_name = ['monoc_eyecorr_hbar_sim3_g' gname];
%     else
%         data_name_hres = ['monoc_eyecorr_hbar_hres_sim2_g' gname];
%         data_name = ['monoc_eyecorr_hbar_sim2_g' gname];
%     end
        data_name_hres = ['monoc_eyecorr_hbar_hres_sim4_g' gname];
        data_name = ['monoc_eyecorr_hbar_sim4_g' gname];
mod_data_name = ['monoc_eyecorr_hbar_sim4_mods_g' gname];
    load(data_name_hres,'HRdrift*')
    load(data_name,'it_fix*','drift_*','dit_xvR2','tr_set','dit_R2','it_*R2','dit_*LLimp','it_*LLimp')
    load(mod_data_name,'all_mod_*R2','all_mod_*LLimp');
    new_usratio = hr_spatial_usfac/spatial_usfac;
    sp_dx = 0.0565/hr_spatial_usfac;
    best_fix_cor = it_fix_post_mean(end,:)*hr_spatial_usfac/spatial_usfac;
    best_fix_std = it_fix_post_std(end,:)*hr_spatial_usfac/spatial_usfac;
    
    cur_newR2dat_name = ['simdata_newR2_g' gname];
    load(cur_newR2dat_name);
    sim_init_cdR2(sg,:) = init_cdR2;
    sim_fin_cdR2(sg,:) = fin_cdR2;
%     sim_LLimp_fac = HRit_xvLLimp
    
    sim_xv_R2_imp(sg,:) = dit_xvR2(end,:)./all_mod_xvR2;
    sim_xv_R2_init(sg,:) = all_mod_xvR2;
    sim_xv_R2_fin(sg,:) = dit_xvR2(end,:);
    
    sim_R2_imp(sg,:) = dit_R2(end,:)./all_mod_R2;
    sim_R2_init(sg,:) = all_mod_R2;
    sim_R2_fin(sg,:) = dit_R2(end,:);

    sim_xv_LL_init(sg,:) = all_mod_xvLLimp;
    sim_xv_LL_fin(sg,:) = dit_xvLLimp(end,:);
    sim_LL_init(sg,:) = all_mod_LLimp;
    sim_LL_fin(sg,:) = dit_LLimp(end,:);

    sim_dit_R2{sg} = [median(it_R2,2); median(dit_R2,2)];
    sim_dit_xR2{sg} = [median(it_xvR2,2); median(dit_xvR2,2)];
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
    %%
    sim_tot_inf(sg,:) = fin_tot_corr/sim_gain;
    sim_tot_std(sg,:) = fin_tot_std/sim_gain;
    
    sim_errs(sg,:) = sim_tot_inf(sg,:) - cur_sim_eyepos(sg,used_inds);
end

%%
sim_unc = nanmedian(sim_tot_std,2);
sim_errSD = robust_std_dev(sim_errs');

%% LOAD THE SET OF TEMPLATE MODELS AND COMPUTE THEIR FILTER PROPERTIES
% load('monoc_eyecorr_hbar_fullXV','dit_mods_LOO');
load('monoc_eyecorr_hbar2','dit_mods');
true_mods = dit_mods{end}(tr_set);
 sp_dx = 0.0565/2;
for ii = 1:length(tr_set)
    ii
%     cur_mod = dit_mods_LOO{ii,end}(tr_set(ii));
    cur_mod = dit_mods{end}(tr_set(ii));
    Xtargs = [cur_mod.mods(:).Xtarget];
    cor_filters = [cur_mod.mods(Xtargs==1).filtK];
    cor_stim_params = cur_mod.stim_params(1);
    after_filt_data(ii) = get_filter_properties(cor_filters,cor_stim_params,sp_dx);
end
%% COMPUTE SIMPLICITY
n_frames = 1e5;
flen = 12;
weight_std = nan(length(tr_set),1);
best_std = nan(length(tr_set),1);

after_bsq_gests = reshape([after_filt_data(:).bsq_gest],[6 length(tr_set)]);
after_sq_gests = reshape([after_filt_data(:).sq_gest],[6 length(tr_set)]);
after_lin_gests = reshape([after_filt_data(:).lin_gest],[6 length(tr_set)]);

sp = true_mods(1).stim_params(1);
dx = sp_dx;
dds = 0.12;
nPix = sp.stim_dims(2);
rand_stim = zeros(n_frames,nPix);
rand_mat = rand(n_frames,nPix);
rand_stim(rand_mat <= dds/2) = -1;
rand_stim(rand_mat >= 1-(dds/2)) = 1;
stim_params = NMMcreate_stim_params([flen nPix],0.01);
test_stim = create_time_embedding(rand_stim,stim_params);

for ii = 1:length(tr_set)
    ii
%     cur_mod = dit_mods_LOO{ii,end}(tr_set(ii));
    cur_mod = dit_mods{end}(tr_set(ii));
    Xtargs = [cur_mod.mods(:).Xtarget];
    cor_filters = [cur_mod.mods(Xtargs==1).filtK];
    filt_out = test_stim*cor_filters;
    filt_out(:,2:end) = filt_out(:,2:end).^2;
    filt_out_var = std(filt_out);
    filt_out_var = filt_out_var/nansum(filt_out_var);
    simp(ii) = filt_out_var(1);
    
    all_stds = [after_lin_gests(2,ii) after_sq_gests(2,ii) ...
        after_bsq_gests(2,ii)];
    cur_ford = after_filt_data(ii).ford;
    weight_std(ii) = nansum(all_stds.*filt_out_var(cur_ford));
    [~,best_filt] = max(filt_out_var);
    best_std(ii) = all_stds(best_filt);
end


%% PLOT ET ERROR AS A FUNCTION OF AVG RF WIDTH
% load('~/Analysis/bruce/ET_final/G086_hbar_mua_RFwidth2.mat');
all_RF_widths = weight_std*2; %use RF-width estimate from the template models

avg_rf_width = mean(all_RF_widths); %2sigma, avg across all units from G086, horizontal
std_rf_width = std(all_RF_widths);
prc_rf_widths = prctile(all_RF_widths,[25 50 75]);
prc_rf_widths_deg = bsxfun(@rdivide,prc_rf_widths,poss_gains');

all_err_prc = prctile(abs(sim_errs)',[25 50 75]);

h=figure;
eh = errorbarxy(prc_rf_widths_deg(:,2),all_err_prc(2,:),prc_rf_widths_deg(:,2)-prc_rf_widths_deg(:,1),prc_rf_widths_deg(:,3)-prc_rf_widths_deg(:,2),...
    all_err_prc(2,:)-all_err_prc(1,:),all_err_prc(3,:)-all_err_prc(2,:),{'k','k','k'});
set(eh(1),'linewidth',1.5);
set(eh(2:end),'linewidth',0.75);
plot(prc_rf_widths_deg(:,2),all_err_prc(2,:),'ko','linewidth',1,'markersize',4)
% plot(prc_rf_widths_deg(:,2),sim_unc,'ro-');
xlabel('Average RF width');
ylabel('Error (deg)');
% xlim([0 0.6])
set(gca,'xscale','log');
xlim([.025 0.6]);
line([0.025 0.6],[1/60 1/60],'color','g','linestyle','--');

%%
fig_width = 3.27;
rel_height = 0.8;
figufy(h);
fname = [fig_dir 'eyeErr_vs_RFwidth4.pdf'];
exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close;

%% PLOT MODEL IMPROVEMENT VS RF WIDTH
all_R2_imp = prctile(sim_R2_imp',[25 50 75]);
all_R2_init = prctile(sim_R2_init',[25 50 75]);
all_R2_fin = prctile(sim_R2_fin',[25 50 75]);

h2 = figure;
eh = errorbarxy(prc_rf_widths_deg(:,2),all_R2_imp(2,:),prc_rf_widths_deg(:,2)-prc_rf_widths_deg(:,1),prc_rf_widths_deg(:,3)-prc_rf_widths_deg(:,2),...
    all_R2_imp(2,:)-all_R2_imp(1,:),all_R2_imp(3,:)-all_R2_imp(2,:),{'k','k','k'});
set(eh(1),'linewidth',1.5);
set(eh(2:end),'linewidth',0.75);
plot(prc_rf_widths_deg(:,2),all_R2_imp(2,:),'ko','linewidth',1,'markersize',4)
xlabel('Average RF width');
ylabel('R2 fold-improvement');
set(gca,'xscale','log');
set(gca,'yscale','log')
xlim([.025 0.6]);
yl = ylim();
ylim([1 11])
plot(all_RF_widths(101)/4,sim_R2_imp(5,101),'ro')

%%
fig_width = 3.27;
rel_height = 0.8;
figufy(h2);
fname = [fig_dir 'R2imp_vs_RFwidth4.pdf'];
exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h2);

%%
use_sg = 3;
sim_gain = poss_gains(use_sg);
switch sim_gain
    case 4
        hr_spatial_usfac = 4;
        spatial_usfac = 2;
    case 2
        hr_spatial_usfac = 4;
        spatial_usfac = 2;
    case 1
        hr_spatial_usfac = 4;
        spatial_usfac = 2;
    case 0.5
        hr_spatial_usfac = 8;
        spatial_usfac = 4;
    case 0.25
        hr_spatial_usfac = 8;
        spatial_usfac = 4;
end

max_shift = round(16*hr_spatial_usfac*sim_gain);
sp_dx = 0.0565/hr_spatial_usfac;
min_unc = 0;
est_unc = sim_tot_std(use_sg,:);
est_unc(est_unc < min_unc) = min_unc;
not_fix_inds = setdiff(1:length(used_inds),fix_inds);
err_m = sim_tot_inf(use_sg,:) - sim_eyepos(used_inds)';
err_z = err_m./est_unc;
err_z(not_fix_inds) = nan;

%points when true eye position was outside decodable range (very rare),
%exclude these
rid = find(abs(sim_eyepos(used_inds)) > max_shift*sp_dx);
err_z(rid) = nan;
err_m(rid) = nan;

%% CREATE PLOT OF ERROR DIST (NORMALIZED AND ABSOLUTE)
n_bins = 100;
% err_mag_b = linspace(0,0.04,n_bins);
% err_magz_b = linspace(0,4,n_bins);
err_mag_b = linspace(-0.04,0.04,n_bins);
err_magz_b = linspace(-4,4,n_bins);

% err_std_est = robust_std_dev(err_m(fix_inds));
err_std_est = nanstd(err_m(fix_inds));

% err_cdist = histc(abs(err_m(fix_inds)),err_mag_b);
% err_cdist = cumsum(err_cdist)/sum(~isnan(err_m(fix_inds)));
err_cdist = histc(err_m(fix_inds),err_mag_b);
err_cdist = err_cdist/sum(err_cdist);

unc_cdist = histc(est_unc(fix_inds),err_mag_b);
unc_cdist = unc_cdist/max(unc_cdist);
% unc_cdist = cumsum(unc_cdist)/sum(~isnan(est_unc(fix_inds)));

% temp = randn(length(est_unc),1).*est_unc';
% norm_cdf = histc(abs(temp),err_mag_b);
% norm_cdf = cumsum(norm_cdf)/sum(~isnan(temp));

% errz_cdist = histc(abs(err_z(fix_inds)),err_magz_b);
% errz_cdist = cumsum(errz_cdist)/sum(~isnan(err_z(fix_inds)));
errz_cdist = histc(err_z(fix_inds),err_magz_b);
errz_cdist = errz_cdist/sum(errz_cdist);

% temp = randn(1e5,1);
% normz_cdf = histc(abs(temp),err_magz_b);
% normz_cdf = cumsum(normz_cdf)/length(temp);
% norm_cdf = normcdf(err_magz_b,0,1);
normz_cdf = normpdf(err_magz_b,0,1);
normz_cdf = normz_cdf/sum(normz_cdf);

h1 = figure;
subplot(2,1,1)
hold on
plot(err_mag_b,err_cdist,'linewidth',1.5);
% plot(err_mag_b,norm_cdf,'r','linewidth',1.5);
legend('Error','Normal');
xlabel('Magnitude (deg)');

% subplot(3,1,2)
% plot(err_mag_b,unc_cdist,'k','linewidth',1.5);
% xlabel('Normalized error magnitude (z)');
% 
subplot(2,1,2)
plot(err_magz_b,errz_cdist,'linewidth',1.5);
hold on
plot(err_magz_b,normz_cdf,'r','linewidth',1.5)
xlabel('Normalized error magnitude (z)');
legend('Error','Normal');

fig_width = 3.27;
rel_height = 1.6;
figufy(h1);
fname = [fig_dir 'Sim_error_analysis.pdf'];
% exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h1);


%%
print_on = true;
fig_width = 4.86; %3.27 4.86 6.83
rel_height = 0.6;
yl = [-0.4 0.4];

offset = 0.05;
close all
n_trials = length(unique(all_trialvec));
H = figure();
% for tt = 1:n_trials
% for tt = [30 31 74 113 121 179 223 239 436 522 556 563 630 658 688 736 746 757 759 913 930 949] %G086
for tt = [31] %G086
    % for tt = [96 137 154 179 376 409]
    uu = find(all_trialvec(used_inds) == tt);
    if ~isempty(uu)
        bt = all_t_axis(used_inds(uu(1)));
        et = all_t_axis(used_inds(uu(end)));
        dur = et-bt;
        if dur > 3.5
         
            hold on
%             h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,sim_tot_inf(1,uu)+offset,sim_tot_std(1,uu),{'color','b'});
            h2=shadedErrorBar(all_t_axis(used_inds(uu))-bt,sim_tot_inf(3,uu),sim_tot_std(3,uu),{'color','r'});
%             h3=shadedErrorBar(all_t_axis(used_inds(uu))-bt,sim_tot_inf(5,uu)-offset,sim_tot_std(5,uu),{'color','g'});
            
%             plot(all_t_axis(used_inds(uu)) - bt,sim_eyepos(used_inds(uu))-offset,'k','linewidth',1)
            plot(all_t_axis(used_inds(uu)) - bt,sim_eyepos(used_inds(uu)),'k','linewidth',1)
%             plot(all_t_axis(used_inds(uu)) - bt,sim_eyepos(used_inds(uu))+offset,'k','linewidth',1)
            
            xlim([0 dur]);
            ylim(yl);
            xlabel('Time (s)','fontsize',10);
            ylabel('Orthoganol position (deg)','fontsize',10);
            title(sprintf('Trial %d',tt));
                        
            if print_on
                delete(h2.edge);
                figufy(H);
                fname = [fig_dir sprintf('SIM_example_T%d',tt)];
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
fig_width = 4.86; %3.27 4.86 6.83
rel_height = 0.6;

delete(h1.edge);
delete(h2.edge);
delete(h3.edge);
figufy(h);
fname = [fig_dir sprintf('sim_RFsize_examptrace2_T%d.pdf',tt)];
exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close;


%% CREATE PLOT OF EXAMPLE neuron at gain4
cd(anal_dir);
true_moddata_name = 'monoc_eyecorr_hbar_mods';
true_eyedata_name = 'monoc_eyecorr_hbar2';
load(true_moddata_name,'all_mod*');
load(true_eyedata_name,'et_tr_set','dit_mods');
su_inds = find(all_mod_SU(et_tr_set) > 0);
cur_su = 7;

true_mod = dit_mods{end}(et_tr_set(su_inds(cur_su)));

sim_gain = 4;
if sim_gain >= 1
    gname = num2str(sim_gain);
elseif sim_gain == 0.5
    gname = 'p5';
elseif sim_gain == 0.25
    gname = 'p25';
else
    error('unsupported');
end
data_name = ['monoc_eyecorr_hbar_sim4_g' gname];
fprintf('Loading %s\n',data_name);
load(data_name,'drift_*','it_fix*','it_R2','dit_R2','it_mods','dit_mods');

pre_mod = it_mods{1}(su_inds(cur_su));
post_mod = dit_mods{end}(su_inds(cur_su));

flen = 12;
dt = 0.01;
use_nPix_us = use_nPix*spatial_usfac;
sp_dx = 0.0565/spatial_usfac/sim_gain;

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
true_filts = true_filts(:,:,[1 3 2]);

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
    set(gca,'xtick',[-0.1:0.05:0.1]);
    set(gca,'ydir','normal');
    %     xlim([-0.5 0.5]);

    subplot(n_stimfilts,3,(ii-1)*3+2);
    imagesc(pix_axis(disp_pix),lag_axis,squeeze(uncor_filts(:,disp_pix,ii))); 
%     caxis([-true_max_vals(ii) true_max_vals(ii)]);
    caxis([-uc_max_vals(ii) uc_max_vals(ii)]);
%     line(pix_axis([1 end]),lag_axis(best_lags_true([ii ii])),'color','r');
    xlabel('Position (deg)');
    ylabel('Lag (ms)');
    set(gca,'xtick',[-0.1:0.05:0.1]);
    set(gca,'ydir','normal');
%     xlim([-0.5 0.5]);

subplot(n_stimfilts,3,(ii-1)*3+3);
    imagesc(pix_axis(disp_pix),lag_axis,squeeze(cor_filts(:,disp_pix,ii))); 
%     caxis([-true_max_vals(ii) true_max_vals(ii)]);
    caxis([-max_vals(ii) max_vals(ii)]);
%     line(pix_axis([1 end]),lag_axis(best_lags_true([ii ii])),'color','r');
    xlabel('Position (deg)');
    ylabel('Lag (ms)');
    set(gca,'xtick',[-0.1:0.05:0.1]);
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
fname = [fig_dir 'RFsim_exampmod2.pdf'];
exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close;

%% 
% all_R2_imp = all_fin_R2./all_init_R2;
% ravg_R2_imp = squeeze(mean(all_R2_imp,3));
% mean_R2_imp = mean(ravg_R2_imp,2);
% std_R2_imp = std(ravg_R2_imp,[],2);
% prc_R2_imp = prctile(all_R2_imp',[25 50 75])';
% 
% fig_width = 3.27;
% rel_height = 0.8;
% 
% uset = find(poss_sim_gains <= 4);
% 
% h2 = figure(); hold on
% % eh = errorbarxy(avg_rf_width./poss_sim_gains(uset),mean_R2_imp(uset),std_rf_width./poss_sim_gains(uset),std_R2_imp(uset),{'k','k','k'});
% eh = errorbarxy(prc_rf_widths_deg(:,2),prc_R2_imp(uset,2),prc_rf_widths_deg(:,2)-prc_rf_widths_deg(:,1),prc_rf_widths_deg(:,3)-prc_rf_widths_deg(:,2),...
%     prc_R2_imp(uset,2)-prc_R2_imp(uset,1),prc_R2_imp(uset,3)-prc_R2_imp(uset,2),{'k','k','k'});
% set(eh(1),'linewidth',1.5);
% set(eh(2:end),'linewidth',0.75);
% % plot(avg_rf_width./poss_sim_gains(uset),mean_R2_imp(uset),'ko','linewidth',1,'markersize',4)
% plot(prc_rf_widths_deg(:,2),prc_R2_imp(uset,2),'ko','linewidth',1,'markersize',4)
% xlabel('RF width');
% ylabel('Model improvement');
% xlim([0 0.6]);
% set(gca,'xscale','log');
% xlim([.025 0.6]);
% 
% figufy(h2);
% fname = [fig_dir 'modimp_vs_RFwidth.pdf'];
% exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close;
