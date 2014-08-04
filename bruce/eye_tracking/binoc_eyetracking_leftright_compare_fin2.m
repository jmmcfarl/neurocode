clear all
% close all

% dir_prefix = '~/';
dir_prefix = '/Volumes/james/';

addpath([dir_prefix 'James_scripts/bruce/eye_tracking/']);
addpath([dir_prefix 'James_scripts/bruce/processing/']);
addpath([dir_prefix 'James_scripts/general_functions']);
Expt_num = 88;
Expt_name = sprintf('G0%d',Expt_num);
data_dir = [dir_prefix 'Data/bruce/' Expt_name];
cd(data_dir);

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = [dir_prefix 'Analysis/bruce/' Expt_name];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

bar_ori = 0;

if bar_ori == 90
left_anal_name = 'binoc_eyecorr_vbar_left_Fin';
right_anal_name = 'binoc_eyecorr_vbar_right_Fin';
else
 left_anal_name = 'binoc_eyecorr_hbar_left_Fin';
right_anal_name = 'binoc_eyecorr_hbar_right_Fin';   
end
recompute_init_mods = 0;

% which_eye_stim = 'left';
use_measured_pos = 0;
use_sac_kerns = 1;
ignore_blocks = [];

%%
xv_frac = 0;

n_fix_inf_it = 3;
n_drift_inf_it = 1;
fix_prior_sigma = 0.175;
drift_sigma = 0.01; %.015
drift_jump_sigma = 0.06; %0.05 start

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
cur_block_set = find(included_type & expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori);
% cur_block_set = find(expt_binoc' & expt_Fr == 1 & expt_bar_ori == bar_ori);

if exist('./fin_aclust_data.mat','file')
    load('./fin_aclust_data.mat');
    [n_aclust_expts,n_aclust_probes] = size(autoclust);
else
    disp('No fin_aclust_data found.');
    autoclust = [];
    n_aclust_expts = 0; n_aclust_probes = 0;
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
%             if strcmp(which_eye_stim,'left')
%                 cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
%             elseif strcmp(which_eye_stim,'right')
%                 cur_stim_mat = double(right_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
%             else
%                 error('Specify right or left eye');
%             end
            
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

%% select submatrix with central pixels
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
temp = find(isnan(interp_sac_stop_inds));
interp_sac_stop_inds(temp) = interp_sac_start_inds(temp);
temp = find(isnan(interp_sac_stop_inds));
interp_sac_stop_inds(temp) = interp_sac_start_inds(temp);

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

%% DEFINE POINTS TO ALLOW RAPID CHANGES (TRIAL STARTS AFTER SACCADES)
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

%%
use_fix_n = 5;
use_drift_n = 3;

cd(anal_dir);
load(left_anal_name);
left_it_post_mean = it_fix_post_mean(use_fix_n,:);
left_it_post_std = it_fix_post_std(use_fix_n,:);
left_drift_post_mean = drift_post_mean(use_drift_n,:);
left_drift_post_std = drift_post_std(use_drift_n,:);

cd(anal_dir);
load(right_anal_name);
right_it_post_mean = it_fix_post_mean(use_fix_n,:);
right_it_post_std = it_fix_post_std(use_fix_n,:);
right_drift_post_mean = drift_post_mean(use_drift_n,:);
right_drift_post_std = drift_post_std(use_drift_n,:);

%% MAKE FINAL DRIFT AND FIXATION CORRECTIONS

finL_fix_corr = nan(NT,1);
finL_fix_std = nan(NT,1);
finR_fix_corr = nan(NT,1);
finR_fix_std = nan(NT,1);
finL_fix_corr(~isnan(fix_ids)) = left_it_post_mean(fix_ids(~isnan(fix_ids)));
finL_fix_corr = interp1(find(~isnan(fix_ids)),finL_fix_corr(~isnan(fix_ids)),1:NT);
finR_fix_corr(~isnan(fix_ids)) = right_it_post_mean(fix_ids(~isnan(fix_ids)));
finR_fix_corr = interp1(find(~isnan(fix_ids)),finR_fix_corr(~isnan(fix_ids)),1:NT);

finL_fix_std(~isnan(fix_ids)) = left_it_post_std(fix_ids(~isnan(fix_ids)));
finL_fix_std = interp1(find(~isnan(fix_ids)),finL_fix_std(~isnan(fix_ids)),1:NT);
finR_fix_std(~isnan(fix_ids)) = right_it_post_std(fix_ids(~isnan(fix_ids)));
finR_fix_std = interp1(find(~isnan(fix_ids)),finR_fix_std(~isnan(fix_ids)),1:NT);

finL_fix_corr = finL_fix_corr*sp_dx;
finL_fix_std = finL_fix_std*sp_dx;
finR_fix_corr = finR_fix_corr*sp_dx;
finR_fix_std = finR_fix_std*sp_dx;

finL_drift_corr = left_drift_post_mean*sp_dx;
finL_drift_std = left_drift_post_std*sp_dx;
finR_drift_corr = right_drift_post_mean*sp_dx;
finR_drift_std = right_drift_post_std*sp_dx;
min_fix_dur = 0.15;
fix_inds = [];
long_fix_inds = [];
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > sac_shift
        finL_drift_corr(cur_inds(1:end-sac_shift+1)) = finL_drift_corr(cur_inds(sac_shift:end));
        finL_drift_std(cur_inds(1:end-sac_shift+1)) = finL_drift_std(cur_inds(sac_shift:end));
        finR_drift_corr(cur_inds(1:end-sac_shift+1)) = finR_drift_corr(cur_inds(sac_shift:end));
        finR_drift_std(cur_inds(1:end-sac_shift+1)) = finR_drift_std(cur_inds(sac_shift:end));
    end
    fix_inds = [fix_inds cur_inds];
    if length(cur_inds)*dt >= min_fix_dur
       long_fix_inds = [long_fix_inds cur_inds]; 
    end
end

finL_tot_corr = finL_fix_corr + finL_drift_corr;
finR_tot_corr = finR_fix_corr + finR_drift_corr;
finL_tot_std = sqrt(finL_fix_std.^2 + finL_drift_std.^2);
finR_tot_std = sqrt(finR_fix_std.^2 + finR_drift_std.^2);

%%
usable_inds = find(~ismember(all_blockvec(used_inds),ignore_blocks));
usable_fix_inds = fix_inds(ismember(fix_inds,usable_inds));
usable_long_fix_inds = long_fix_inds(ismember(long_fix_inds,usable_inds));

sparse_blocks = find(expt_dds(cur_block_set) == 12);
sparse_inds = find(ismember(all_blockvec(used_inds),sparse_blocks));
sparse_usable_inds = sparse_inds(ismember(sparse_inds,usable_inds));
sparse_fix_inds = sparse_inds(ismember(sparse_inds,usable_fix_inds));
sparse_long_fix_inds = sparse_inds(ismember(sparse_inds,usable_long_fix_inds));

min_unc = 0.001;
avg_unc = sqrt(finL_tot_std.^2 + finR_tot_std.^2);
avg_unc(avg_unc < min_unc) = min_unc;
weights = 1./avg_unc;
[weighted_r,raw_r] = weighted_correlation(finR_tot_corr(sparse_inds),finL_tot_corr(sparse_inds),weights(sparse_inds));

%%
cd ~/Analysis/bruce/summary_analysis/eyetrack_figs/

max_unc = 0.02;
cert_data = sparse_inds(finR_tot_std(sparse_inds) < max_unc & finL_tot_std(sparse_inds) < max_unc);

unc_corr_fix = corr(finL_fix_corr(cert_data)',finR_fix_corr(cert_data)');
unc_corr_tot = corr(finL_tot_corr(cert_data)',finR_tot_corr(cert_data)');


bin_ax = linspace(-0.15,0.15,50);
LRdiff = finL_tot_corr(cert_data) - finR_tot_corr(cert_data);
figure;hold on
y = histc(LRdiff,bin_ax);
y = y/sum(y);
stairs(bin_ax,y,'r')
xlim(bin_ax([1 end]));
% set(gca,'xtick',[],'ytick',[]);
yl = ylim();
line([0 0],yl,'color','k');
box off
set(gca,'fontsize',10,'fontname','arial');
title('Diff position','fontsize',10);
fillPage(gcf,'papersize',[4 4]);
fname = sprintf('LR_diff_dist_sure_bar%d',bar_ori);

%%

figure
[h,dens_det] = DensityPlot_jmm(finL_tot_corr(sparse_inds),finR_tot_corr(sparse_inds),'ynormal','xrange',[-0.5 0.5],'yrange',[-0.5 0.5],'sd',[3 3],'sqrtsc');
line([-0.5 0.5],[-0.5 0.5],'color','w');
xlabel('Left eye position (deg)','fontsize',12);
ylabel('Right eye position (deg)','fontsize',12);
set(gca,'fontname','arial','fontsize',10);
axis square
fillPage(gcf,'papersize',[5 5]);
set(gca,'xtick',-0.5:0.25:0.5,'ytick',-0.5:0.25:5);
fname = sprintf('Left_right_dens_bar%d',bar_ori);

figure;
bin_ax = linspace(-0.5,0.5,50);
y = histc(finL_tot_corr(sparse_inds),bin_ax);
y = y/sum(y);
stairs(bin_ax,y);
xlim(bin_ax([1 end]));
set(gca,'xtick',[],'ytick',[]);
yl = ylim();
line([0 0],yl,'color','k');
title('Left position','fontsize',10);
fillPage(gcf,'papersize',[4 4]);
fname = sprintf('Left_dist_bar%d',bar_ori);

figure;
y = histc(finR_tot_corr(sparse_inds),bin_ax);
y = y/sum(y);
stairs(bin_ax,y);
xlim(bin_ax([1 end]));
set(gca,'xtick',[],'ytick',[]);
yl = ylim();
line([0 0],yl,'color','k');
title('Right position','fontsize',10);
fillPage(gcf,'papersize',[4 4]);
fname = sprintf('Right_dist_bar%d',bar_ori);

LRdiff = finL_tot_corr(sparse_inds) - finR_tot_corr(sparse_inds);
figure;
y = histc(LRdiff,bin_ax);
y = y/sum(y);
stairs(bin_ax,y);
xlim(bin_ax([1 end]));
set(gca,'xtick',[],'ytick',[]);
yl = ylim();
line([0 0],yl,'color','k');
title('Diff position','fontsize',10);
fillPage(gcf,'papersize',[4 4]);
fname = sprintf('LR_diff_dist_bar%d',bar_ori);

%%
close all
n_trials = length(unique(all_trialvec));
for tt = 1:n_trials
    fprintf('Trial %d of %d\n',tt,n_trials);
% for tt = [12 15 29 170 192 242]
uu = find(all_trialvec(used_inds) == tt);
    if ~isempty(uu)
        bt = all_t_axis(used_inds(uu(1)));
        et = all_t_axis(used_inds(uu(end)));
        dur = et-bt;
        if dur > 3
            hold on
            h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,finL_tot_corr(uu),finL_tot_std(uu),{'color','b'});
            h2=shadedErrorBar(all_t_axis(used_inds(uu))-bt,finR_tot_corr(uu),finR_tot_std(uu),{'color','r'});
            if bar_ori == 0
                h3=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),2),'k');
%                 h4=plot(all_t_axis(used_inds(uu))-bt,corrected_interp_eyevals(used_inds(uu),4)-median(corrected_interp_eyevals(used_inds(uu),4)),'color',[0.2 0.8 0.2])
            else
                h3=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),1),'k');
%                 h4=plot(all_t_axis(used_inds(uu))-bt,corrected_interp_eyevals(used_inds(uu),3)-median(corrected_interp_eyevals(used_inds(uu),3)),'color',[0.2 0.8 0.2])
            end
%             legend([h1.mainLine h2.mainLine h3 h4],{'Left-eye inferred','Right-eye inferred','Left-eye measured','Right-eye measured'})
            xlim([0 dur]);
            ylim([-0.5 0.5]);
            xlabel('Time (s)','fontsize',10);
            ylabel('Orthoganol position (deg)','fontsize',10);
            title(sprintf('Trial %d',tt));
            set(gca,'fontsize',8,'fontname','arial');
            box off
            fillPage(gcf,'papersize',[8 5]);
            pause
            clf
        end
    end
end
