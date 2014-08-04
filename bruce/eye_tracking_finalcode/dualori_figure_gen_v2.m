clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

Expt_num = 91;
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
fig_dir = '/home/james/Analysis/bruce/ET_final/';

cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
% cluster_dir = ['/Volumes/james/Analysis/bruce/' Expt_name '/clustering'];

mod_data_name = 'dualori_eyecorr_mods2';
old_anal_name = 'dualori_eyecorr5';
anal_name = 'dualori_eyecorr_highres';
hrmod_data_name = 'dualori_eyecorr_mods_hres';

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

n_drift_inf_it = 1; %3

drift_noise_sigma = sqrt(2*0.004^2);
drift_prior_sigma = sqrt(2*0.004^2); 
drift_jump_sigma = sqrt(2*0.075^2); 
% drift_noise_sigma = 0.003;
% drift_prior_sigma = 0.004; %.004 may be best here
% drift_jump_sigma = 0.075; %0.05 start
drift_dsf = 2;

min_trial_dur = 0.75;

spatial_usfac = 4;
old_usfac = 2;

%%
eps = -1e3;

stim_fs = 100; %in Hz
dt = 0.01;
full_nPix = 36;
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

max_Dshift = round(8*spatial_usfac);

%%
include_expts = {'rls.orRC'};
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
cur_block_set = find(included_type & ~expt_binoc' & expt_Fr == 1);


cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];

n_blocks = length(cur_block_set);

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
all_stim_ori = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
% all_spk_times = cell(n_probes,1);
% all_spk_inds = cell(n_probes,1);
% all_clust_ids = cell(n_probes,1);

trial_toffset = zeros(length(cur_block_set),1);
cur_spkind_offset = 0;
cur_toffset = 0;
for ee = 1:n_blocks;
    fprintf('Expt %d Block %d of %d;  UNMATCHED EXPT TYPE\n',Expt_num,ee,n_blocks);
    cur_block = cur_block_set(ee);
    
%     fname = [cluster_dir sprintf('/Block%d_Clusters.mat',cur_block)];
%     load(fname,'Clusters');
%     for cc = 1:n_probes
%         all_spk_times{cc} = cat(1,all_spk_times{cc},Clusters{cc}.times + cur_toffset);
%         all_spk_inds{cc} = cat(1,all_spk_inds{cc},Clusters{cc}.spk_inds + cur_spkind_offset);
%         all_clust_ids{cc} = cat(1,all_clust_ids{cc},Clusters{cc}.spike_clusts);
%     end
    
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
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start/1e4;
        cur_stim_ori = Expts{cur_block}.Trials(use_trials(tt)).or;
        n_frames = size(left_stim_mats{use_trials(tt)},1);
        if n_frames > 0
            if length(cur_stim_times) == 1
                cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
                cur_stim_times(cur_stim_times > trial_end_times(use_trials(tt))) = [];
            end
        end
        cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if ~any(isnan(left_stim_mats{use_trials(tt)}(:))) && n_frames > min_trial_dur/dt
            use_frames = min(length(cur_stim_times),n_frames);
%             cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
%             cur_stim_ori = cur_stim_ori(1:use_frames);
            %             bar_Xmat = create_time_embedding(cur_stim_mat,stim_params);
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times];
            all_t_axis = [all_t_axis; cur_t_axis];
%             all_stim_mat = [all_stim_mat; cur_stim_mat];
%             all_stim_ori = [all_stim_ori; cur_stim_ori];
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
lin_correction = false;
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,zeros(length(cur_block_set),1),used_inds,lin_correction);

[saccades,et_params] = detect_saccades(corrected_eye_vals,all_eye_speed,all_eye_ts,et_params);

par_thresh = 1.25;
orth_thresh = 1.25;
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
cd(anal_dir);

fprintf('Loading pre-computed initial models\n');
load(mod_data_name);
load(old_anal_name,'dit_mods','tr_set','it_fix_post_mean','drift_post_mean','it_fix_post_std','drift_post_std','it_LL*','dit_LL*');
old_best_mods = dit_mods{end};

best_fix_cor = squeeze(it_fix_post_mean(end,:,:))*spatial_usfac/old_usfac;
best_fix_std = squeeze(it_fix_post_std(end,:,:))*spatial_usfac/old_usfac;
best_drift_cor = squeeze(drift_post_mean(end,:,:))*spatial_usfac/old_usfac;
best_drift_std = squeeze(drift_post_std(end,:,:))*spatial_usfac/old_usfac;
best_LLimp = it_LLimp(end,:);
best_dLLimp = dit_LLimp(end,:);

clear it_fix_post_mean drift_post_mean
clear all_mod_fits

%%
load(anal_name);

%%
fin_fix_corr = nan(NT,2);
fin_fix_std = nan(NT,2);
fin_fix_corr(~isnan(fix_ids),:) = best_fix_cor(fix_ids(~isnan(fix_ids)),:);
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids),:),1:NT);
fin_fix_std(~isnan(fix_ids),:) = best_fix_std(fix_ids(~isnan(fix_ids)),:);
fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids),:),1:NT);

fin_fix_corr = fin_fix_corr*sp_dx;
fin_fix_std = fin_fix_std*sp_dx;

fin_drift_corr = squeeze(drift_post_mean(end,:,:))*sp_dx;
fin_drift_std = squeeze(drift_post_std(end,:,:))*sp_dx;

min_fix_dur = 0.15;
fix_inds = [];
long_fix_inds = [];
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > sac_shift
        fin_drift_corr(cur_inds(1:end-sac_shift+1),:) = fin_drift_corr(cur_inds(sac_shift:end),:);
        fin_drift_std(cur_inds(1:end-sac_shift+1),:) = fin_drift_std(cur_inds(sac_shift:end),:);
    end
    fix_inds = [fix_inds cur_inds];
    if length(cur_inds)*dt >= min_fix_dur
        long_fix_inds = [long_fix_inds cur_inds];
    end
end

fin_tot_corr = fin_fix_corr + fin_drift_corr;
fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2);
    
%%
for bb = 1:length(cur_block_set)
%    binds = find(all_blockvec == bb);
   binds = used_inds(all_blockvec(used_inds) == bb);
   corrected_eye_vals_interp(binds,:) = bsxfun(@minus,corrected_eye_vals_interp(binds,:),nanmedian(corrected_eye_vals_interp(binds,:)));
end


%% MICROSACCADE AMP COMPARISON
% measured_seq = smooth_eyepos(:,[2 4]);
measured_seq = corrected_eye_vals_interp(used_inds,:);

sac_buff_inds = 2;
min_fix_dur = 0.1;

saccade_prefix = fix_ids(saccade_start_inds);
saccade_postfix = fix_ids(saccade_stop_inds);
saccade_prefix_dur = nan(size(saccade_start_inds));
saccade_postfix_dur = nan(size(saccade_start_inds));
saccade_prefix_dur(~isnan(saccade_prefix)) = fix_durs(saccade_prefix(~isnan(saccade_prefix)));
saccade_postfix_dur(~isnan(saccade_postfix)) = fix_durs(saccade_postfix(~isnan(saccade_postfix)));

too_short = find(saccade_prefix_dur  < min_fix_dur | saccade_postfix_dur < min_fix_dur);
long_enough = setdiff(1:length(saccade_start_inds),too_short);

start_pts = saccade_start_inds(long_enough) - sac_buff_inds;
end_pts = saccade_stop_inds(long_enough) + sac_buff_inds;
start_pts(start_pts < 1) = 1; end_pts(end_pts > length(used_inds)) = length(used_inds);
m_pre_pos = measured_seq(start_pts,:);
m_post_pos = measured_seq(end_pts,:);
m_delta_pos = (m_post_pos - m_pre_pos);

inferred_pre_pos = fin_tot_corr(start_pts,:);
inferred_post_pos = fin_tot_corr(end_pts,:);
inferred_delta_pos = (inferred_post_pos - inferred_pre_pos);

saccade_blocks = all_blockvec(used_inds(saccade_start_inds));

use_micros = ismember(long_enough',micro_sacs) & ~isnan(m_delta_pos(:,1)) & ~isnan(inferred_delta_pos(:,1));
use_nonmicros = ~ismember(long_enough',micro_sacs) & ~isnan(m_delta_pos(:,1)) & ~isnan(inferred_delta_pos(:,1));
% use_micros = ismember(saccade_blocks(long_enough),sim_sac_expts) & ismember(long_enough',micro_sacs) & ~isnan(m_delta_pos(:,1)) & ~isnan(inferred_delta_pos);

[msac_corrs,msac_pvals] = corr(m_delta_pos(use_micros,:),inferred_delta_pos(use_micros,:),'type','spearman');

%%
% measured_seq = smooth_eyepos;
measured_seq = corrected_eye_vals_interp(used_inds,:);
min_fix_dur = 0.1;
long_fix_dur = 0.1;

inferred_drift = nan(size(fin_tot_corr));
measured_drift = nan(length(fin_tot_corr),4);
inferred_fix_avg = nan(n_fixs,2);
measured_fix_avg = nan(n_fixs,4);
fix_ind_vec = nan(length(fin_tot_corr),1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    cur_inf = fin_tot_corr(cur_inds,:);
    inferred_fix_avg(ii,:) = nanmedian(fin_tot_corr(cur_inds,:));
    inferred_drift(cur_inds,:) = bsxfun(@minus,cur_inf,inferred_fix_avg(ii,:));
    
    measured_fix_avg(ii,:) = nanmedian(measured_seq(cur_inds,:));
    measured_drift(cur_inds,:) = bsxfun(@minus,measured_seq(cur_inds,:),measured_fix_avg(ii,:));
    
    fix_ind_vec(cur_inds) = ii;
end
fix_dur_vec = nan(length(fin_tot_corr),1);
fix_dur_vec(~isnan(fix_ind_vec)) = fix_durs(fix_ind_vec(~isnan(fix_ind_vec)));

used_fixs = find(fix_durs > long_fix_dur);
used_fix_inds = find(fix_dur_vec >= long_fix_dur);

%%
close all
sqrt_scale = true;
n_bins = 40;
xr = [-0.1 0.1];
h1=figure;
[h,det] = DensityPlot_jmm(inferred_drift(used_fix_inds,2),measured_drift(used_fix_inds,2),'ynormal','xrange',xr,'yrange',xr,'sd',[3 3]);
% DensityPlot_jmm(inferred_drift(used_fix_inds),measured_drift(used_fix_inds,1),'ynormal','xrange',xr,'yrange',xr,'sd',[3 3],'sqrtsc');
eps = -0.5;
lZ = log10(det.z);
lZ(lZ < eps) = eps;
imagesc(det.x(1,:),det.y(:,1),lZ);
set(gca,'ydir','normal');
caxis([eps 3.25]);
line([-0.5 0.5],[-0.5 0.5],'color','w');
xlabel('Inferred amplitude (deg)');
ylabel('Measured amplitude (deg)');
xlim(xr); ylim(xr);
set(gca,'xtick',-0.2:0.05:0.2,'ytick',-0.2:0.05:0.2);
axis square

b_ax = linspace(xr(1),xr(2),n_bins);

INF_nx = histc(inferred_drift(used_fix_inds,2),b_ax);
INF_nx = INF_nx/sum(INF_nx);

MEAS_nx = histc(measured_drift(used_fix_inds,2),b_ax);
MEAS_nx = MEAS_nx/sum(MEAS_nx);

DIFF_nx = histc(inferred_drift(used_fix_inds,2) - measured_drift(used_fix_inds,2),b_ax);
DIFF_nx = DIFF_nx/sum(DIFF_nx);

h5=figure; hold on
stairs(b_ax,INF_nx,'b');
stairs(b_ax,MEAS_nx,'r');
stairs(b_ax,DIFF_nx,'k');
xlim(xr); 
set(gca,'xtick',-0.2:0.05:0.2,'ytick',[]);
% set(gca,'yscale','log');
% ylim([0.0001 2]);

%%
fig_width = 3.27; 
rel_height = 0.85;

figufy(h1);
fname = [fig_dir 'dualori_drift_density.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h1);

figufy(h5);
fname = [fig_dir 'dualori_drift_dists.pdf'];
exportfig(h5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h5);

%%
sqrt_scale = true;
n_bins = 40;
xr = [-0.4 0.4];

close all

h1=figure;
if sqrt_scale
    [hs,det] = DensityPlot_jmm(inferred_delta_pos(use_micros,2),m_delta_pos(use_micros,2),'ynormal','xrange',[-0.5 0.5],'yrange',[-0.5 0.5],'sd',[3 3],'sqrtsc');
else
    DensityPlot_jmm(inferred_delta_pos(use_micros,2),m_delta_pos(use_micros,2),'ynormal','xrange',[-0.5 0.5],'yrange',[-0.5 0.5],'sd',[3 3]);
end
% plot(inferred_delta_pos(use_micros),m_delta_pos(use_micros,1),'.','markersize',4);
% xlim([-0.5 0.5]); ylim([-0.5 0.5]);
line([-0.5 0.5],[-0.5 0.5],'color','w');
xlabel('Inferred amplitude (deg)');
ylabel('Measured amplitude (deg)');
xlim(xr); ylim(xr);
set(gca,'xtick',-0.5:0.1:0.5,'ytick',-0.5:0.1:0.5);
axis square


b_ax = linspace(-0.5,0.5,n_bins);

INF_nx = histc(inferred_delta_pos(use_micros,2),b_ax);
INF_nx = INF_nx/sum(INF_nx);

MEAS_nx = histc(m_delta_pos(use_micros,2),b_ax);
MEAS_nx = MEAS_nx/sum(MEAS_nx);

DIFF_nx = histc(inferred_delta_pos(use_micros,2) - m_delta_pos(use_micros,2),b_ax);
DIFF_nx = DIFF_nx/sum(DIFF_nx);

h5=figure; hold on
stairs(b_ax,INF_nx,'b');
stairs(b_ax,MEAS_nx,'r');
stairs(b_ax,DIFF_nx,'k');
xlim(b_ax([1 end])); 
xlim(xr); 
set(gca,'xtick',-0.5:0.1:0.5,'ytick',[]);

%% PRINT FILTER BEFORE/AFTER COMPARISONS
% fig_width = 3.27; 
fig_width = 3.27; 
rel_height = 0.85;

figufy(h1);
fname = [fig_dir 'dualori_microsac_density.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h1);

figufy(h5);
fname = [fig_dir 'dualori_microsac_alldiff.pdf'];
exportfig(h5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h5);

%%
n_bins = 40;
xr = [-0.4 0.4];
sqrt_scale = true;

usable_inds = 1:length(fin_tot_corr);
h1=figure;
if sqrt_scale
    DensityPlot_jmm(fin_tot_corr(usable_inds,2),measured_seq(usable_inds,2),'ynormal','xrange',[-0.5 0.5],'yrange',[-0.5 0.5],'sd',[3 3],'sqrtsc');
else
    DensityPlot_jmm(fin_tot_corr(usable_inds,2),measured_seq(usable_inds,2),'ynormal','xrange',[-0.5 0.5],'yrange',[-0.5 0.5],'sd',[3 3]);    
end
% plot(fin_tot_corr(usable_inds),measured_seq(usable_inds,1),'.','markersize',1);
line([-0.5 0.5],[-0.5 0.5],'color','w');
xlabel('Inferred position (deg)','fontsize',12);
ylabel('Measured position (deg)','fontsize',12);
set(gca,'fontsize',10,'fontname','arial');
set(gca,'xtick',-0.5:0.1:0.5,'ytick',-0.5:0.1:0.5);
axis square
xlim(xr); ylim(xr);
fillPage(gcf,'papersize',[5 5]);

b_ax = linspace(-0.5,0.5,n_bins);
INF_nx = histc(fin_tot_corr(usable_inds,2)',b_ax);
INF_nx = INF_nx/sum(INF_nx);

h3=figure; hold on
MEAS_nx = histc(measured_seq(usable_inds,2),b_ax);
MEAS_nx = MEAS_nx/sum(MEAS_nx);
stairs(b_ax,MEAS_nx,'r');
xlim([-0.4 0.4]);
yl = ylim();
line([0 0],yl,'color','k')
xlim(b_ax([1 end])); 
set(gca,'xtick',[],'ytick',[]);
title('Measured position','fontsize',10);

DIFF_nx = histc(fin_tot_corr(usable_inds,2) - measured_seq(usable_inds,2),b_ax);
DIFF_nx = DIFF_nx/sum(DIFF_nx);

h5=figure; hold on
stairs(b_ax,INF_nx,'b');
stairs(b_ax,MEAS_nx,'r');
stairs(b_ax,DIFF_nx,'k');
set(gca,'xtick',-0.5:0.1:0.5,'ytick',[]);
xlim(xr);
%%
figufy(h5);
fname = [fig_dir 'dualori_alleyepos_alldiff.pdf'];
exportfig(h5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h5);

figufy(h1);
fname = [fig_dir 'dualori_alleyepos_dens.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h1);


%% FOR SACCADES
% 
% measured_pre_Ly = [saccades(:).pre_Ly];
% measured_post_Ly = [saccades(:).post_Ly];
% measured_pre_Lx = [saccades(:).pre_Lx];
% measured_post_Lx = [saccades(:).post_Lx];
% measured_pre_Ry = [saccades(:).pre_Ry];
% measured_post_Ry = [saccades(:).post_Ry];
% measured_pre_Rx = [saccades(:).pre_Rx];
% measured_post_Rx = [saccades(:).post_Rx];
% use_ss = find(saccade_start_inds > 5 & saccade_start_inds < (NT-10));
% 
% min_fix_dur = 0.15;
% s_amp_buffer = 2;
% 
% fix_durs = (fix_stop_inds - fix_start_inds)*dt;
% sac_fix_dur = nan(length(use_ss),1);
% for ii = 1:length(use_ss)
%     cor_fix = find(fix_start_inds == saccade_stop_inds(use_ss(ii)));
%     if isempty(cor_fix)
%         cor_fix = find(fix_start_inds > saccade_stop_inds(use_ss(ii)),1);
%     end
%     sac_fix_dur(ii) = fix_durs(cor_fix);
% end
% too_short = find(sac_fix_dur < min_fix_dur);
% fprintf('Eliminating %d of %d too short sacs\n',length(too_short),length(use_ss));
% use_ss(too_short) = [];
% 
% use_micros = (ismember(use_ss,micro_sacs));
% use_nonmicros = (~ismember(use_ss,micro_sacs));
% 
% clear measured_delta_*
% measured_delta_posL(:,1) = measured_post_Lx(used_saccade_set(use_ss)) - measured_pre_Lx(used_saccade_set(use_ss));
% measured_delta_posL(:,2) = measured_post_Ly(used_saccade_set(use_ss)) - measured_pre_Ly(used_saccade_set(use_ss));
% measured_delta_posR(:,1) = measured_post_Rx(used_saccade_set(use_ss)) - measured_pre_Rx(used_saccade_set(use_ss));
% measured_delta_posR(:,2) = measured_post_Ry(used_saccade_set(use_ss)) - measured_pre_Ry(used_saccade_set(use_ss));
% measured_delta_posA = 0.5*measured_delta_posL + 0.5*measured_delta_posR;
% 
% saccade_blocks = cur_block_set(all_blockvec(used_inds(saccade_start_inds)));
% 
% m_pre_pos = smooth_eyepos(used_inds(saccade_start_inds(use_ss)- s_amp_buffer),:);
% m_post_pos = smooth_eyepos(used_inds(saccade_stop_inds(use_ss) + s_amp_buffer),:);
% m_delta_pos = (m_post_pos - m_pre_pos);
% 
% inferred_pre_pos = fin_tot_corr(saccade_start_inds(use_ss) - s_amp_buffer,:);
% inferred_post_pos = fin_tot_corr(saccade_stop_inds(use_ss) + s_amp_buffer,:);
% inferred_delta_pos = (inferred_post_pos - inferred_pre_pos);
% 
% u = find(all(~isnan(measured_delta_posA),2) & all(~isnan(inferred_delta_pos),2));
% [sac_corrs,sac_pvals] = corr([measured_delta_posL(u,:) measured_delta_posR(u,:) measured_delta_posA(u,:)],inferred_delta_pos(u,:),'type','spearman');
% [sm_sac_corrs,sm_sac_pvals] = corr(m_delta_pos(u,:),inferred_delta_pos(u,:),'type','spearman');
% 
% u_micro = u(use_micros(u));
% [msac_corrs,msac_pvals] = corr([measured_delta_posL(u_micro,:) measured_delta_posR(u_micro,:) measured_delta_posA(u_micro,:)],inferred_delta_pos(u_micro,:),'type','spearman');
% [sm_msac_corrs,sm_msac_pvals] = corr(m_delta_pos(u_micro,:),inferred_delta_pos(u_micro,:),'type','spearman');
% [msac_corrs_LR,msac_pvals_LR] = corr([measured_delta_posL(u_micro,:) measured_delta_posR(u_micro,:)],'type','spearman');
% [sm_msac_corrs_LR,sm_msac_pvals_LR] = corr(m_delta_pos(u_micro,:),'type','spearman');
% 
% fprintf('\nSmooth sac corr:\nLH:%.3f LV:%.3f \nRH:%.3f RV:%.3f\n',sm_msac_corrs(1,1),sm_msac_corrs(2,2),sm_msac_corrs(3,1),sm_msac_corrs(4,2));
% fprintf('\nDsac corr:\nLH:%.3f LV:%.3f \nRH:%.3f RV:%.3f\n',msac_corrs(1,1),msac_corrs(2,2),msac_corrs(3,1),msac_corrs(4,2));
% 
%%
sqrt_scale = true;
n_bins = 40;
xr = [-0.4 0.4];

close all

h1=figure;
if sqrt_scale
    [hs,det] = DensityPlot_jmm(inferred_delta_pos(use_micros,1),m_delta_pos(use_micros,1),'ynormal','xrange',[-0.5 0.5],'yrange',[-0.5 0.5],'sd',[3 3],'sqrtsc');
else
    DensityPlot_jmm(inferred_delta_pos(use_micros,1),m_delta_pos(use_micros,1),'ynormal','xrange',[-0.5 0.5],'yrange',[-0.5 0.5],'sd',[3 3]);
end
% plot(inferred_delta_pos(use_micros),m_delta_pos(use_micros,1),'.','markersize',4);
% xlim([-0.5 0.5]); ylim([-0.5 0.5]);
line([-0.5 0.5],[-0.5 0.5],'color','w');
xlabel('Inferred amplitude (deg)');
ylabel('Measured amplitude (deg)');
xlim(xr); ylim(xr);
set(gca,'xtick',-0.5:0.1:0.5,'ytick',-0.5:0.1:0.5);
axis square

h2=figure;
if sqrt_scale
    [hs,det] = DensityPlot_jmm(inferred_delta_pos(use_micros,2),m_delta_pos(use_micros,2),'ynormal','xrange',[-0.5 0.5],'yrange',[-0.5 0.5],'sd',[3 3],'sqrtsc');
else
    DensityPlot_jmm(inferred_delta_pos(use_micros,2),m_delta_pos(use_micros,2),'ynormal','xrange',[-0.5 0.5],'yrange',[-0.5 0.5],'sd',[3 3]);
end
% plot(inferred_delta_pos(use_micros),m_delta_pos(use_micros,1),'.','markersize',4);
% xlim([-0.5 0.5]); ylim([-0.5 0.5]);
line([-0.5 0.5],[-0.5 0.5],'color','w');
xlabel('Inferred amplitude (deg)');
ylabel('Measured amplitude (deg)');
xlim(xr); ylim(xr);
set(gca,'xtick',-0.5:0.1:0.5,'ytick',-0.5:0.1:0.5);
axis square

%% FIGURE GENERATOR
% fig_width = 6.83; %3.27 4.86 6.83
% relheight = 0.8;
% 
% print_on = false;
% 
% close all
% n_trials = length(unique(all_trialvec));
% H = figure();
% yl = [-0.4 0.4];
% interesting_trials = [6 52 148 170 200 216 217 249 254 301 386 415 418 429 460 531 602];
% print_trials = [52 148 170 249];
% % print_trials = [249];
% % for tt = 1:n_trials
% % for tt = interesting_trials
%     for tt = print_trials
%     uu = find(all_trialvec(used_inds) == tt);
%     if ~isempty(uu)
%         bt = all_t_axis(used_inds(uu(1)));
%         et = all_t_axis(used_inds(uu(end)));
%         dur = et-bt;
%         if dur > 3
%             cur_sac_inds = find(ismember(saccade_start_inds,uu));
%             rel_sac_start_times = all_t_axis(used_inds(saccade_start_inds(cur_sac_inds))) - bt;
%             rel_sac_end_times = all_t_axis(used_inds(saccade_stop_inds(cur_sac_inds))) - bt;
%             
%             plot_dim = 1;
% %             subplot(2,1,1)
%             title(sprintf('Trial %d',tt));
%             hold on
%             h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_tot_corr(uu,plot_dim),fin_tot_std(uu,plot_dim),{'color','k'});
% %             plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),plot_dim),'r','linewidth',1.5);
% %             plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),plot_dim+2),'color','b','linewidth',1.5);
%             
%             xlim([0 dur]);
%             ylim(yl);
%              xlabel('Time (s)');
%             ylabel('Horizontal position (deg)');
%             for ii = 1:length(rel_sac_start_times)
%                 line(rel_sac_start_times([ii ii]),yl,'color',[0.2 0.8 0.2],'linestyle','--');
%             end           
%             line([0 dur],[0 0],'color','k','linestyle','--');
%             
%             plot_dim = 2;
% %             subplot(2,1,2)
%             hold on
%             h2=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_tot_corr(uu,plot_dim),fin_tot_std(uu,plot_dim),{'color','r'});
% %             plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),plot_dim),'r','linewidth',1.5);
% %             plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),plot_dim+2),'color','b','linewidth',1.5);
%             xlim([0 dur]);
%             ylim(yl);
%              xlabel('Time (s)');
%             ylabel('Vertical position (deg)');
%             for ii = 1:length(rel_sac_start_times)
%                 line(rel_sac_start_times([ii ii]),yl,'color',[0.2 0.8 0.2],'linestyle','--');
%             end
%              line([0 dur],[0 0],'color','k','linestyle','--');
% %              xlim([0 2.5]);
%             if print_on
%                 delete(h1.edge);
%                 delete(h2.edge);
%                 figufy(H);
%                 fname = [fig_dir sprintf('examp_dualori_%s_T%d',Expt_name,tt)];
%                 exportfig(H,fname,'width',fig_width,'height',relheight*fig_width,'fontmode','scaled','fontsize',1);
%                 close(H);
%             else
%                 pause
%                 clf(H);
%             end
% 
%         end
%     end
% end

%%
%Example trial 671
fig_dir = ['/home/james/Analysis/bruce/ET_final/'];

yr = [-0.4 0.4];
print_on = false;
fig_width = 4.86; %3.27 4.86 6.83
rel_height = 1.2;

close all
H = figure();
n_trials = length(unique(all_trialvec));
% for tt = 1:n_trials
%     for tt = [96 137 154 179 376 409]
    for tt = [96 217 249 320 366 400 410 453 485 507 510 610 671]
% for tt = [217 671]
    uu = find(all_trialvec(used_inds) == tt);
    if ~isempty(uu)
        bt = all_t_axis(used_inds(uu(1)));
        et = all_t_axis(used_inds(uu(end)));
        dur = et-bt;
        if dur > 3
            cur_sac_inds = find(ismember(saccade_start_inds,uu));
            rel_sac_start_times = all_t_axis(used_inds(saccade_start_inds(cur_sac_inds))) - bt;
            rel_sac_end_times = all_t_axis(used_inds(saccade_stop_inds(cur_sac_inds))) - bt;
            
            plot_dim = 2;
            subplot(2,1,1)
            hold on
            h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_tot_corr(uu,plot_dim),fin_tot_std(uu,plot_dim),{'color','k'});
            l1=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),plot_dim),'r','linewidth',2);
            xlim([0 dur]);
            ylim(yr);
            
            yl = ylim();
            for ii = 1:length(rel_sac_start_times)
                line(rel_sac_start_times([ii ii]),yl,'color','g','linestyle','--');
            end            
            xlabel('Time (s)','fontsize',10);
            ylabel('Orthoganol position (deg)','fontsize',10);
            title(sprintf('Trial %d',tt));
            set(gca,'fontsize',8,'fontname','arial');

            plot_dim = 1;
            subplot(2,1,2)
            hold on
            h2=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_tot_corr(uu,plot_dim),fin_tot_std(uu,plot_dim),{'color','k'});
            
            l2=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),plot_dim),'r','linewidth',2);
            
            xlim([0 dur]);
            ylim(yr);
            yl = ylim();
            for ii = 1:length(rel_sac_start_times)
                line(rel_sac_start_times([ii ii]),yl,'color','g','linestyle','--');
            end            
            xlabel('Time (s)','fontsize',10);
            ylabel('Orthoganol position (deg)','fontsize',10);
            title(sprintf('Trial %d',tt));
            
            if print_on
                delete(h1.edge);
                delete(h2.edge);
                figufy(H);
                fname = [fig_dir sprintf('Dualori_example_trace_%s_T%d',Expt_name,tt)];
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
% close all
% f1 = figure();
% f2 = figure();
% for ss = 1:length(tr_set)
%     % sbeg = find(all_mod_SU(tr_set) > 0,1);
%     % for ss = sbeg:length(tr_set)
%     ss
%     init_mod = all_mod_fits(tr_set(ss));
%     xtargs = [init_mod.mods(:).Xtarget];
%     kmat = [init_mod.mods(xtargs == 1).filtK];
%     figure(f1); clf
%     subplot(2,2,1)
%     imagesc(reshape(kmat(:,1),flen,use_nPix_us));
%     ca = max(abs(kmat(:,1))); caxis([-ca ca]);
%     for ii = 1:(size(kmat,2)-1)
%         subplot(2,2,2+ii)
%         imagesc(reshape(kmat(:,ii+1),flen,use_nPix_us));
%         ca = max(abs(kmat(:,ii+1))); caxis([-ca ca]);
%     end
%     colormap(gray)
%     
%     fin_mod = dit_mods{end}(tr_set(ss));
%     xtargs = [fin_mod.mods(:).Xtarget];
%     kmat = [fin_mod.mods(xtargs == 1).filtK];
%     figure(f2); clf
%     subplot(2,2,1)
%     imagesc(reshape(kmat(:,1),flen,use_nPix_us));
%     ca = max(abs(kmat(:,1))); caxis([-ca ca]);
%     for ii = 1:(size(kmat,2)-1)
%         subplot(2,2,2+ii)
%         imagesc(reshape(kmat(:,ii+1),flen,use_nPix_us));
%         ca = max(abs(kmat(:,ii+1))); caxis([-ca ca]);
%     end
%     colormap(gray)
%     
%     fprintf('Cell %d of %d\n',ss,length(tr_set));
%     fprintf('Original: %.4f  Fin: %.4f\n',all_mod_LLimp(tr_set(ss)),dit_LLimp(end,tr_set(ss)));
%     pause
% end
% 
