clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');
% fig_dir = '/Volumes/james/Analysis/bruce/ET_final/';
fig_dir = '/home/james/Analysis/bruce/ET_final/';

%run on 86 [0, 90];
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

bar_ori = 0;

base_Nrpts = 1;
poss_Nunits = [90 80 70 60 50 40 30 20 10 5];
poss_Nrpts = round(linspace(3,5,length(poss_Nunits)));

% poss_Nrpts = base_Nrpts*ones(size(poss_Nunits));

%add MU only data point
poss_Nunits = [poss_Nunits nan];
poss_Nrpts = [poss_Nrpts 1];

if bar_ori == 0
    true_mod_data = 'monoc_eyecorr_hbar_mods';
    true_data = 'monoc_eyecorr_hbar2';
else
    true_mod_data = 'monoc_eyecorr_vbar_mods';
    true_data = 'monoc_eyecorr_vbar2';
end

use_measured_pos = 0;
use_sac_kerns = 1;
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

n_fix_inf_it = 4; %4
n_drift_inf_it = 1; %2

fix_prior_sigma = 0.15;
fix_noise_sigma = 0.1;
drift_noise_sigma = 0.003;
drift_prior_sigma = 0.004; %.004 may be best here
drift_jump_sigma = 0.075; %0.05 start
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
            if ~isempty(all_stim_times)
                if any(cur_stim_times < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
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
tr_inds = full_inds;


%% LOAD BEST EYE POSITION DATA AND MODELS
%LOAD DATA HERE
cd(anal_dir)
load(true_mod_data);
load(true_data,'et_tr_set','drift_post*','*_R2','it_fix*');

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
% for ii = 1:n_fixs
%     cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
%     if length(cur_inds) > sac_shift
%         fin_drift_corr(cur_inds(1:end-sac_shift+1)) = fin_drift_corr(cur_inds(sac_shift:end));
%         fin_drift_std(cur_inds(1:end-sac_shift+1)) = fin_drift_std(cur_inds(sac_shift:end));
%     end
% end
for ii = 1:length(trial_start_inds)
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
    fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
end
fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT);
fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);

best_inf_ep = fin_fix_corr(:) + fin_drift_corr(:);
best_std_ep = sqrt(fin_fix_std(:).^2 + fin_drift_std(:).^2);
best_init_R2 = it_R2(1,:)';
best_fin_R2 = dit_R2(end,:)';

%%
clear it_inf_ep it_std_ep it_fin_R2
cd(anal_dir)
max_nrpts = max(poss_Nrpts);
it_inf_ep = nan(length(fin_drift_corr),length(poss_Nunits),max_nrpts);
it_std_ep = nan(length(fin_drift_corr),length(poss_Nunits),max_nrpts);
it_fin_R2 = nan(length(et_tr_set),length(poss_Nunits),max_nrpts);

for cc = 1:length(poss_Nunits)
    
    Nunits = poss_Nunits(cc);
    Nrpts = poss_Nrpts(cc);
    
    for rr = 1:Nrpts

%         cur_data_name = sprintf('popsize_%d_it%d',Nunits,rr);
        cur_data_name = sprintf('popsize_v2_%d_it%d',Nunits,rr);
        fprintf('Loading data %s\n',cur_data_name);
        load(cur_data_name);
        
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
        %         for ii = 1:n_fixs
        %             cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
        %             if length(cur_inds) > sac_shift
        %                 fin_drift_corr(cur_inds(1:end-sac_shift+1)) = fin_drift_corr(cur_inds(sac_shift:end));
        %                 fin_drift_std(cur_inds(1:end-sac_shift+1)) = fin_drift_std(cur_inds(sac_shift:end));
        %             end
        %         end
        for ii = 1:length(trial_start_inds)
            cur_inds = trial_start_inds(ii):trial_end_inds(ii);
            fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
            fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
        end
        fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT);
        fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);
        
        it_inf_ep(:,cc,rr) = fin_fix_corr(:) + fin_drift_corr(:);
        it_std_ep(:,cc,rr) = sqrt(fin_fix_std(:).^2 + fin_drift_std(:).^2);
        
        if size(dit_R2,2) > length(et_tr_set)
            it_fin_R2(:,cc,rr) = dit_R2(end,et_tr_set);
        else
            it_fin_R2(:,cc,rr) = dit_R2(end,:);
        end
    end
end

%%
% clear Cit_inf_ep Cit_std_ep Cit_fin_R2
% cd(anal_dir)
% max_nrpts = max(poss_Nrpts);
% Cit_inf_ep = nan(length(fin_drift_corr),length(poss_Nunits),max_nrpts);
% Cit_std_ep = nan(length(fin_drift_corr),length(poss_Nunits),max_nrpts);
% Cit_fin_R2 = nan(length(et_tr_set),length(poss_Nunits),max_nrpts);
% 
% for cc = 1:length(poss_Nunits)
%     
%     Nunits = poss_Nunits(cc);
%     Nrpts = poss_Nrpts(cc);
%     
%     for rr = 1:Nrpts
% 
%         cur_data_name = sprintf('popsize_Cprior_%d_it%d',Nunits,rr);
%         fprintf('Loading data %s\n',cur_data_name);
%         load(cur_data_name);
%         
%         fin_fix_corr = nan(NT,1);
%         fin_fix_std = nan(NT,1);
%         fin_fix_corr(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
%         fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
%         fin_fix_std(~isnan(fix_ids)) = it_fix_post_std(end,fix_ids(~isnan(fix_ids)));
%         fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);
%         
%         fin_fix_corr = fin_fix_corr*sp_dx;
%         fin_fix_std = fin_fix_std*sp_dx;
%         
%         fin_drift_corr = drift_post_mean(end,:)*sp_dx;
%         fin_drift_std = drift_post_std(end,:)*sp_dx;
%         for ii = 1:n_fixs
%             cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
%             if length(cur_inds) > sac_shift
%                 fin_drift_corr(cur_inds(1:end-sac_shift+1)) = fin_drift_corr(cur_inds(sac_shift:end));
%                 fin_drift_std(cur_inds(1:end-sac_shift+1)) = fin_drift_std(cur_inds(sac_shift:end));
%             end
%         end
%         
%         Cit_inf_ep(:,cc,rr) = fin_fix_corr(:) + fin_drift_corr(:);
%         Cit_std_ep(:,cc,rr) = sqrt(fin_fix_std(:).^2 + fin_drift_std(:).^2);
%         
%         if size(dit_R2,2) > length(et_tr_set)
%             Cit_fin_R2(:,cc,rr) = dit_R2(end,et_tr_set);
%         else
%             Cit_fin_R2(:,cc,rr) = dit_R2(end,:);
%         end
%     end
% end

%%
% best_imps = best_fin_R2(et_tr_set)./best_init_R2(et_tr_set);
% it_imps = bsxfun(@rdivide,it_fin_R2,best_init_R2(et_tr_set));
% avg_it_imps = squeeze(mean(it_imps));
% 
% % Cit_imps = bsxfun(@rdivide,Cit_fin_R2,best_init_R2(et_tr_set));
% % Cavg_it_imps = squeeze(mean(Cit_imps));
% 
% h1=figure;
% errorbar(poss_Nunits,nanmean(avg_it_imps,2),nanstd(avg_it_imps,[],2),'color','k');
% hold on
% % errorbar(poss_Nunits,nanmean(Cavg_it_imps,2),nanstd(Cavg_it_imps,[],2),'color','r');
% xl = xlim();
% line(xl,[avg_it_imps(end,1) avg_it_imps(end,1)],'color','r','linestyle','--');
% line(xl,[mean(best_imps) mean(best_imps)],'color','g','linestyle','--');
% ylim([1 2.5]);
% xlabel('Population size');
% ylabel('Average model improvement');


%%
% fig_width = 3.27; 
% rel_height = 0.8;
% 
% figufy(h1);
% fname = [fig_dir 'popsize_modimp.pdf'];
% exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h1);


%% FIG WITH ROBUST ERROR
ep_diffs = bsxfun(@minus,it_inf_ep,best_inf_ep);
ep_errs = squeeze(robust_std_dev(ep_diffs));
ep_unc = squeeze(nanmedian(it_std_ep));
ep_tot_err = robust_std_dev(best_inf_ep);
ep_best_unc = nanmedian(best_std_ep);
ep_R2 = squeeze(1-bsxfun(@rdivide,ep_errs.^2,ep_tot_err.^2));
% ep_err_prc = prctile(abs(ep_diffs),[25 50 75]);
% ep_err_meds = squeeze(ep_err_prc(2,:,:));
% ep_err_LE = squeeze(ep_err_prc(2,:,:)-ep_err_prc(1,:,:));
% ep_err_RE = squeeze(ep_err_prc(3,:,:)-ep_err_prc(2,:,:));

% Cep_diffs = bsxfun(@minus,Cit_inf_ep,best_inf_ep);
% Cep_errs = squeeze(robust_std_dev(Cep_diffs));
% Cep_tot_err = robust_std_dev(best_inf_ep);
% Cep_R2 = squeeze(1-bsxfun(@rdivide,Cep_errs.^2,Cep_tot_err.^2));

h2=figure; hold on
errorbar(poss_Nunits,nanmean(ep_errs,2),nanstd(ep_errs,[],2),'color','k');
% errorbar(poss_Nunits,nanmean(Cep_errs,2),nanstd(Cep_errs,[],2),'color','r');
% errorbar(poss_Nunits,nanmean(ep_err_meds,2),nanmean(ep_err_LE,2),nanmean(ep_err_RE,2),'color','k');
xl = xlim();
line(xl,[ep_errs(end,1) ep_errs(end,1)],'color','r','linestyle','--');
line(xl,[ep_tot_err ep_tot_err],'color','k','linestyle','--');
xlabel('Population size');
ylabel('SD of difference (deg)');

% cmap = jet(10);
% for ii = 1:10
%     plot(poss_Nunits(ii),nanmean(ep_errs(ii,:),2),'o','color',cmap(ii,:),'linewidth',2);
% end
% 

h3 = figure;
errorbar(poss_Nunits,nanmean(ep_unc,2),nanstd(ep_unc,[],2),'color','k');
xl = xlim();
line(xl,[ep_unc(end,1) ep_unc(end,1)],'color','r','linestyle','--');
line(xl,[ep_tot_err ep_tot_err],'color','k','linestyle','--');
xlabel('Population size');
ylabel('Uncertainty (deg)');

%% FIG WITH REGULAR SD

ep_diffs = bsxfun(@minus,it_inf_ep,best_inf_ep);
ep_errs = squeeze(std(ep_diffs));
ep_tot_err = std(best_inf_ep);
ep_R2 = squeeze(1-bsxfun(@rdivide,ep_errs.^2,ep_tot_err.^2));

h2=figure; hold on
errorbar(poss_Nunits,nanmean(ep_errs,2),nanstd(ep_errs,[],2),'color','k');
% errorbar(poss_Nunits,nanmean(Cep_errs,2),nanstd(Cep_errs,[],2),'color','r');
xl = xlim();
line(xl,[ep_errs(end,1) ep_errs(end,1)],'color','r','linestyle','--');
line(xl,[ep_tot_err ep_tot_err],'color','k','linestyle','--');
xlabel('Population size');
ylabel('Error SD (deg)');

%% FIG WITH UNC-WEIGHTED SD

min_unc = 0;
ep_diffs = bsxfun(@minus,it_inf_ep,best_inf_ep);
ep_uncs = max(cat(4,it_std_ep,repmat(best_std_ep,[1 size(it_std_ep,2) size(it_std_ep,3)])),[],4);
ep_uncs(ep_uncs < min_unc) = min_unc;
ep_uncs = 1./ep_uncs;
ep_uncs = bsxfun(@rdivide,ep_uncs,sum(ep_uncs));
ep_errs = sqrt(squeeze(sum(ep_uncs.*ep_diffs.^2)));
tot_unc = best_std_ep;
tot_unc(tot_unc < min_unc) = min_unc;
tot_unc = 1./tot_unc;
tot_unc = tot_unc/sum(tot_unc);
ep_tot_err = sqrt(sum(tot_unc.*best_inf_ep.^2));
ep_R2 = squeeze(1-bsxfun(@rdivide,ep_errs.^2,ep_tot_err.^2));


h2=figure; hold on
errorbar(poss_Nunits,nanmean(ep_errs,2),nanstd(ep_errs,[],2),'color','k');
% errorbar(poss_Nunits,nanmean(Cep_errs,2),nanstd(Cep_errs,[],2),'color','r');
xl = xlim();
line(xl,[ep_errs(end,1) ep_errs(end,1)],'color','r','linestyle','--');
line(xl,[ep_tot_err ep_tot_err],'color','k','linestyle','--');
xlabel('Population size');
ylabel('Error SD (deg)');

cmap = jet(10);
for ii = 1:10
    plot(poss_Nunits(ii),nanmean(ep_errs(ii,:),2),'o','color',cmap(ii,:),'linewidth',2);
end

%%
fig_width = 3.27; 
rel_height = 0.8;

figufy(h2);
fname = [fig_dir 'popsize_errrSD.pdf'];
exportfig(h2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h2);
%% COMPARE ROBUST STATS
ep_diffs = bsxfun(@minus,it_inf_ep,best_inf_ep);

n_pts = 200;
xx = linspace(0,0.125,n_pts);
ydist = nan(11,5,n_pts);
for cc = 1:11
    for rr = 1:5
        yobs = abs(squeeze(ep_diffs(:,cc,rr)));
        temp = histc(yobs,xx);
        temp = cumsum(temp)/sum(~isnan(yobs));
        ydist(cc,rr,:) = temp;
    end
end
yobs = squeeze(abs(ep_diffs(:,end,1)));
temp = histc(yobs,xx);
temp = cumsum(temp)/sum(~isnan(yobs));
ydist_mua = temp;

yobs = squeeze(abs(best_inf_ep));
temp = histc(yobs,xx);
temp = cumsum(temp)/sum(~isnan(yobs));
ydist_best = temp;

cmap = jet(10);
h1=figure; hold on
for ii = 1:10
    plot(xx,squeeze(nanmean(ydist(ii,:,:),2)),'color',cmap(ii,:),'linewidth',1);
%     pause
end
plot(xx,ydist_mua,'k--','linewidth',2);
plot(xx,ydist_best,'k','linewidth',2)
yl = ylim();
% ylim([1e-4 yl(2)]);
% set(gca,'yscale','log');
xlim(xx([1 end]));

%%
fig_width = 3.27; 
rel_height = 0.8;

figufy(h1);
fname = [fig_dir 'popsize_errDISTS.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h1);

%%
cc = 1; rr = 3;
yobs = squeeze(ep_diffs(:,cc,rr));
cur_unc = sqrt(best_std_ep.^2 + squeeze(it_std_ep(:,cc,rr)).^2);
cur_unc = 1./cur_unc;
cur_unc = cur_unc/sum(cur_unc);
wSD = sqrt(sum(cur_unc.*squeeze(it_std_ep(:,cc,rr)).^2));
IQRSD = diff(prctile(yobs,[25 75]))/1.35;

xx = linspace(-0.06,0.06,500);
cur_ydist = ksdensity(yobs,xx);
g1 = normpdf(xx,0,std(yobs));
g2 = normpdf(xx,0,robust_std_dev(yobs));
g3 = normpdf(xx,0,wSD);

h1 = figure; hold on
plot(xx,cur_ydist,'linewidth',1)
plot(xx,g1,'r',xx,g2,'k',xx,g3,'m','linewidth',1);
xlim(xx([1 end]));

% xx = linspace(-0.5,0.5,200);
% cur_ydist = ksdensity(yobs,xx);
% g1 = normpdf(xx,0,std(yobs));
% g2 = normpdf(xx,0,robust_std_dev(yobs));
% g3 = normpdf(xx,0,wSD);
% 
% figure; hold on
% plot(xx,cur_ydist)
% plot(xx,g1,'r',xx,g2,'k',xx,g3,'m');
% xlim(xx([1 end]));
% set(gca,'yscale','log');
% ylim([1e-4 100]);

%%
fig_width = 3.27; 
rel_height = 0.8;

figufy(h1);
fname = [fig_dir 'popsize_robust_examp.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h1);