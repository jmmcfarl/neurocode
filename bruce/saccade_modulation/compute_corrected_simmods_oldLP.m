clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

Expt_num = 232;

switch Expt_num
    case 232
        bar_ori = 50;
    case 235 
        bar_ori = 30;
    case 239
        bar_ori = 130;
end

Expt_name = sprintf('M%d',Expt_num);
data_dir = ['~/Data/bruce/' Expt_name];
cd(data_dir);

load ./random_bar_eyedata_ftime.mat bar_expts

load(['lemM' num2str(Expt_num) 'Expts.mat']);
load ./bar_params.mat

anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final/'];
fin_anal_dir = ['~/Analysis/bruce/' Expt_name '/stim_mods/'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end
if ~exist(fin_anal_dir,'dir')
    system(['mkdir ' fin_anal_dir]);
end

mod_data_name = 'monoc_eyecorr_mods';
et_anal_name = 'monoc_eyecorr';
fin_mod_name = 'corr_mods_quick';


%dont fit stim models using these blocks
if Expt_num == 235
    ignore_blocks = [51]; %G086
elseif Expt_num == 239
    ignore_blocks = [40];
else
   ignore_blocks = []; 
end

%%
xv_frac = 0;

flen = 15;

n_bar_pos = bar_params.n_bars;
bar_axis = bar_params.bar_axis;

min_trial_dur = 0.75;

spatial_usfac = 4;
spatial_tbspace = 2;

%%
eps = -1e3;

stim_fs = 100; %in Hz
dt = 0.01;
stim_params = NIMcreate_stim_params([flen n_bar_pos],dt);
Fr = 1;

beg_buffer = 0.2;
end_buffer = 0.05;
trial_dur = 2;

n_probes = 24;

use_right_eye = true;

n_use_blocks = Inf;

use_nPix_us = n_bar_pos*spatial_usfac;
use_nPix_tb = use_nPix_us/spatial_tbspace;
klen_us = use_nPix_us*flen;
klen_tb = use_nPix_tb*flen;

sac_backlag = round(0.05/dt);
sac_forlag = round(0.3/dt);
sac_bincents = -sac_backlag:sac_forlag;
n_sac_bins = length(sac_bincents);

sp_dx = 0.125/spatial_usfac;
max_shift = round(10*spatial_usfac);

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
cur_block_set = bar_expts;
cur_block_set(ismember(cur_block_set,ignore_blocks)) = [];

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
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);

cur_spkind_offset = 0;
trial_toffset = zeros(length(cur_block_set),1);
cur_toffset = 0;
for ee = 1:n_blocks;
    fprintf('Expt %d Block %d of %d;  UNMATCHED EXPT TYPE\n',Expt_num,ee,n_blocks);
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
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start/1e4;
        cur_bar_Op = [Expts{cur_block}.Trials(use_trials(tt)).Op];
        n_frames = length(cur_bar_Op);
        cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
        cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
        
        if n_frames > min_trial_dur/dt
            use_frames = min(length(cur_stim_times),n_frames);
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times+cur_toffset < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times + cur_toffset];
            all_t_axis = [all_t_axis; cur_t_axis + cur_toffset];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges + cur_toffset];
            all_stim_mat = [all_stim_mat cur_bar_Op];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
    end
    trial_cnt = trial_cnt + n_trials;
    trial_toffset(ee) = all_t_bin_edges(end);
    cur_toffset = trial_toffset(ee);
    cur_spkind_offset = cur_spkind_offset + round(32e3*(max(trial_end_times)-min(trial_start_times) + 50));
end

%%
all_bar_mat = zeros(length(all_stim_times),n_bar_pos);
for b = 1:n_bar_pos
    cur_set = find(all_stim_mat==bar_axis(b));
    all_bar_mat(cur_set,b) = 1;
end

full_nPix_us = spatial_usfac*n_bar_pos;
if spatial_usfac == 2
    all_bar_mat_up = zeros(length(all_stim_times),full_nPix_us);
    for ii = 1:n_bar_pos
        all_bar_mat_up(:,2*(ii-1)+1) = all_bar_mat(:,ii);
        all_bar_mat_up(:,2*(ii-1)+2) = all_bar_mat(:,ii);
    end
elseif spatial_usfac == 1
    all_bar_mat_up = all_bar_mat;
elseif spatial_usfac == 4
    all_bar_mat_up = zeros(length(all_stim_times),full_nPix_us);
    for ii = 1:n_bar_pos
        all_bar_mat_up(:,4*(ii-1)+(1:4)) = repmat(all_bar_mat(:,ii),1,4);
   end
else
    error('Unsupported spatial Up-sample factor!');
end

full_nPix_tb = full_nPix_us/spatial_tbspace;

stim_params_us = NMMcreate_stim_params([flen full_nPix_us],dt);
stim_params_tb = NMMcreate_stim_params([flen full_nPix_tb],dt);

all_Xmat = create_time_embedding(all_bar_mat,stim_params);
all_Xmat_us = create_time_embedding(all_bar_mat_up,stim_params_us);

if spatial_tbspace > 1
    %create a tent-basis (triangle) filter
    tent_filter = [(1:spatial_tbspace)/spatial_tbspace 1-(1:spatial_tbspace-1)/spatial_tbspace]/spatial_tbspace;
    %apply to the stimulus
    filtered_stim = zeros(size(all_bar_mat_up));
    for i = 1:length(tent_filter)
        filtered_stim = filtered_stim + shift_mat_zpad(all_bar_mat_up,i-spatial_tbspace,2)*tent_filter(i);
    end
    all_bar_mat_tb = filtered_stim(:,spatial_tbspace:spatial_tbspace:end);
    all_Xmat_tb = create_time_embedding(all_bar_mat_tb,stim_params_tb);
    
    n_tbasis = size(all_bar_mat_tb,2);
    tent_basis = zeros(n_tbasis,use_nPix_us);
    for ii = 1:n_tbasis
        cur_inds = (ii-1)*spatial_tbspace + (1:length(tent_filter));
        uu = find(cur_inds >= 1 & cur_inds <= use_nPix_us);
        tent_basis(ii,cur_inds(uu)) = tent_filter(uu);
    end
else
    all_Xmat_tb = create_time_embedding(all_bar_mat_up,stim_params_tb);
end

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

double_spike_buffer = 3; %number of samples (in either direction) to exclude double spikes from adjacent-probe SUs
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
%     cur_mua_inds = find(all_clust_ids{cc} >= 1);
    cur_mua_inds = find(all_clust_ids{cc} >= 0);
    
%     %remove spikes from isolated SUs on the same probe from the MU
%     for ss = 1:length(unique_su_nums)
%         cur_mua_inds(ismember(all_spk_inds{cc}(cur_mua_inds),all_su_spk_inds{unique_su_nums(ss)})) = [];
%     end
% 
%     nearby_probes = [cc-1 cc+1]; nearby_probes(nearby_probes < 1 | nearby_probes > n_probes) = [];
%     cur_set = find(ismember(SU_block_probes,nearby_probes));
%     if ~isempty(cur_set)
%         [cur_SS,cur_BB] = ind2sub([length(SU_numbers) length(cur_used_blocks)],cur_set);
%     else
%         cur_SS = [];
%     end
%     unique_su_nums = unique(cur_SS); %set of SU numbers picked up on adjacent probes
%     if ~isempty(unique_su_nums)
%         double_spikes = [];
%         for ss = 1:length(unique_su_nums)
%             cur_blocked_inds = bsxfun(@plus,all_su_spk_inds{unique_su_nums(ss)},-double_spike_buffer:double_spike_buffer);
%             double_spikes = [double_spikes; find(ismember(all_spk_inds{cc}(cur_mua_inds),cur_blocked_inds))];
%         end
%         fprintf('Eliminating %d of %d double spikes in MUA\n',length(double_spikes),length(cur_mua_inds));
%         cur_mua_inds(double_spikes) = [];
%     end
    [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc}(cur_mua_inds),all_t_bin_edges);
    cur_spkhist(all_bin_edge_pts) = [];
    all_binned_mua(:,cc) = cur_spkhist;
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
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

%compute corrected eye data in bar-oriented frame
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,bar_ori*ones(length(cur_block_set),1),used_inds);

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

%% INCORPORATE INFERRED EYE-POSITIONS
cd(anal_dir)
load(et_anal_name);

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

fin_tot_corr = fin_fix_corr(:) + fin_drift_corr(:);
fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2);

fin_tot_corr_rnd = round(fin_tot_corr/sp_dx);
fin_tot_corr_rnd(isnan(fin_tot_corr_rnd)) = 0;
max_shift = size(all_bar_mat_up,2);
fin_tot_corr_rnd(fin_tot_corr_rnd > max_shift) = max_shift;
fin_tot_corr_rnd(fin_tot_corr_rnd < - max_shift) = -max_shift;
%%
all_stimmat_cor = all_bar_mat_up;
for ii = 1:NT
    all_stimmat_cor(used_inds(ii),:) = shift_matrix_Nd(all_stimmat_cor(used_inds(ii),:),-fin_tot_corr_rnd(ii),2);
end
all_stimmat_cor = all_stimmat_cor*tent_basis';
    all_Xmat_cor = create_time_embedding(all_stimmat_cor,stim_params_tb);
all_Xmat_cor = all_Xmat_cor(used_inds,:);

%% Create set of TR and XV trials
use_trials = unique(all_trialvec(used_inds));
tr_inds = 1:length(used_inds);

%% SELECT USABLE UNITS AND make Robs_mat
cd(anal_dir);
load(mod_data_name,'all_mod_SU*');

tr_set = et_tr_set;
full_n_chs = length(all_mod_SU);

% make Robs_mat
poss_blocks = 1:(n_blocks-1);
full_Robs_mat = nan(length(used_inds),full_n_chs);
for ss = 1:full_n_chs
    if all_mod_SU(ss) > 0
        su_probe_ind = find(SU_numbers == all_mod_SUnum(ss));
        full_Robs_mat(:,ss) = all_binned_sua(used_inds,su_probe_ind);
    elseif ~isnan(all_mod_SU(ss))
        full_Robs_mat(:,ss) = all_binned_mua(used_inds,ss);
    end
end

%%
cd(fin_anal_dir);
if exist([fin_mod_name '.mat'],'file')
fprintf('Loading precomputed models\n');
load(fin_mod_name);
end
%%
base_lambda_d2XT = 200;
base_lambda_L1 = 10;
base_qlambda_d2XT = 50;
base_qlambda_L1 = 5;

close all
silent =1;

% for ss = (n_probes+1):full_n_chs
for ss = 1:full_n_chs
    fprintf('Unit %d/%d\n',ss,full_n_chs);
    
    cur_Robs = full_Robs_mat(:,ss);
    cur_tr_inds = tr_inds(~isnan(cur_Robs(tr_inds)));
    if ~isempty(cur_tr_inds)
    init_reg_params = NMMcreate_reg_params('lambda_d2XT',10);
    fin_stim_params = NMMcreate_stim_params([flen use_nPix_tb],dt);
    
    stim_mod_signs = [1];
    stim_NL_types = {'lin'};
    qfilts = find(strcmp('quad',stim_NL_types));
    nqfilts = find(~strcmp('quad',stim_NL_types));
    init_mod = NMMinitialize_model( fin_stim_params, stim_mod_signs, stim_NL_types, init_reg_params);
    Xtargs = [init_mod.mods(:).Xtarget];
    
    init_fitN = ceil(length(cur_tr_inds)/10);
    init_fit_subset = randperm(length(cur_tr_inds));
    init_fit_subset = init_fit_subset(1:init_fitN);
    
    init_mod = NMMfit_filters(init_mod,cur_Robs(cur_tr_inds(init_fit_subset)),all_Xmat_cor(cur_tr_inds(init_fit_subset),:),[],[],silent);
    [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(init_mod,cur_Robs(cur_tr_inds),all_Xmat_cor(cur_tr_inds,:));
    init_mod = NMMadjust_regularization(init_mod,nqfilts,'lambda_d2XT',base_lambda_d2XT./var(gint(:,nqfilts))');
    init_mod = NMMadjust_regularization(init_mod,nqfilts,'lambda_L1',base_lambda_L1./std(gint(:,nqfilts))');
%     init_mod = NMMadjust_regularization(init_mod,qfilts,'lambda_d2XT',base_qlambda_d2XT./var(gint(:,qfilts))');
%     init_mod = NMMadjust_regularization(init_mod,qfilts,'lambda_L1',base_qlambda_L1./std(gint(:,qfilts))');
%     init_mod = NMMfit_filters(init_mod,cur_Robs(cur_tr_inds(init_fit_subset)),all_Xmat_cor(cur_tr_inds(init_fit_subset),:),[],[],silent);
    cor_gqm(ss).mod_fit = NMMfit_filters(init_mod,cur_Robs(cur_tr_inds),all_Xmat_cor(cur_tr_inds,:),[],[],silent);
    [LL, penLL, pred_rate, G, gint,fgint,nullLL] = NMMmodel_eval(cor_gqm(ss).mod_fit,cur_Robs(cur_tr_inds),all_Xmat_cor(cur_tr_inds,:));

    full_LL = nansum(cur_Robs.*log(cur_Robs)-cur_Robs)/sum(cur_Robs);
    cor_gqm(ss).nullLL = nullLL;
    cor_gqm(ss).fullLL = full_LL;
    cor_gqm(ss).filt_mags = std(gint);
    cor_gqm(ss).filt_conts = std(fgint);
    cor_gqm(ss).tot_spks = sum(cur_Robs);
    cor_gqm(ss).tot_samps = length(cur_Robs);
   
    %%
    save(fin_mod_name,'cor_gqm');
    end
end
