clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

Expt_num = 277;

switch Expt_num
    case 266
        bar_ori = 80;
    case 270
        bar_ori = 60;
    case 275 
        bar_ori = 135;
    case 277
        bar_ori = 70;
end

if Expt_num == 275 || Expt_num == 277
    rpt_seed = 1001; %M275 M277
else
    rpt_seed = 1e4; %m270 and 266
end

Expt_name = sprintf('M%d',Expt_num);
data_dir = ['~/Data/bruce/' Expt_name];
cd(data_dir);

load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final/'];
fin_anal_dir = ['~/Analysis/bruce/' Expt_name '/stim_mods/'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end
if ~exist(fin_anal_dir,'dir')
    system(['mkdir ' fin_anal_dir]);
end

cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];

mod_data_name = 'monoc_eyecorr_mods';
et_anal_name = 'monoc_eyecorr';
fin_mod_name = 'corr_mods';

recompute_init_mods = 0;
use_measured_pos = 0;
use_sac_kerns = 1;
use_LOOXV = 0;
use_coils = [0 0]; %[L R]       

if any(use_coils > 0)
   anal_name = [anal_name '_Cprior'];
   old_anal_name = [old_anal_name '_Cprior'];
end
if use_measured_pos == 1
    mod_data_name = [mod_data_name '_Cinit'];
    anal_name = [anal_name '_Cinit'];
end

%dont fit stim models using these blocks
if Expt_num == 270
    ignore_blocks = [5 19];
elseif Expt_num == 275
    ignore_blocks = 15;
else
   ignore_blocks = []; 
end

if Expt_num==270
    scale_fac = 1.72;
else
    scale_fac = 1;
end

%%
xv_frac = 'rpt';

flen = 12;
use_nPix = 32;

n_fix_inf_it = 1; %3
n_drift_inf_it = 1; %3

min_trial_dur = 0.75;

spatial_usfac = 2;
%%
eps = -1e3;

stim_fs = 100; %in Hz
dt = 0.01;

if Expt_num==270
    full_nPix=32;
else
    full_nPix = 36;
end
stim_params = NIMcreate_stim_params([flen full_nPix],dt);
Fr = 1;
sp_dx = 0.0565/spatial_usfac/scale_fac;

beg_buffer = 0.2;
end_buffer = 0.05;
trial_dur = 4;

n_probes = 24;

use_right_eye = false;

n_use_blocks = Inf;

use_nPix_us = use_nPix*spatial_usfac;
klen_us = use_nPix_us*flen;

sac_backlag = round(0.2/dt);
sac_forlag = round(0.4/dt);
sac_bincents = -sac_backlag:sac_forlag;
n_sac_bins = length(sac_bincents);

%%
include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaXFaXFs','rls.AllSac','rls.imiXFa'};
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

if Expt_num==275
    expt_bar_ori(expt_bar_ori > 1e4) = bar_ori;
end

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
% all_Xmat = [];
all_stim_mat = [];
all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trial_Se = [];
all_trial_blk = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
trial_toffset = zeros(length(cur_block_set),1);
cur_toffset = 0;
cur_spkind_offset = 0;
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
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)' + cur_toffset);
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)' + cur_toffset);
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    
    trial_Se = [Expts{cur_block}.Trials(:).se];
    trial_Se = trial_Se(id_inds);
    all_trial_Se = cat(1,all_trial_Se,trial_Se(use_trials)');
    all_trial_blk = cat(1,all_trial_blk,ones(length(use_trials),1)*ee);

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
    cur_spkind_offset = cur_spkind_offset+ round(32e3*(max(trial_end_times)-min(trial_start_times) + 50));
end

%%
full_nPix_us = spatial_usfac*full_nPix;
if spatial_usfac == 2
    all_stimmat_up = zeros(size(all_stim_mat,1),full_nPix_us);
    for ii = 1:size(all_stim_mat,2)
        all_stimmat_up(:,2*(ii-1)+1) = all_stim_mat(:,ii);
        all_stimmat_up(:,2*(ii-1)+2) = all_stim_mat(:,ii);
    end
    elseif spatial_usfac == 4
    all_stimmat_up = zeros(size(all_stim_mat,1),full_nPix_us);
    for ii = 1:size(all_stim_mat,2)
        all_stimmat_up(:,4*(ii-1)+1) = all_stim_mat(:,ii);
        all_stimmat_up(:,4*(ii-1)+2) = all_stim_mat(:,ii);
        all_stimmat_up(:,4*(ii-1)+3) = all_stim_mat(:,ii);
        all_stimmat_up(:,4*(ii-1)+4) = all_stim_mat(:,ii);
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
% LOAD REFCLUSTERS
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

%for SU probes
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
        cur_blocks = [cur_blocks; find(SU_ID_mat(:,cur_clust) == SU_numbers(ss))];
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
    if ~isempty(all_su_spk_times{ss})
    cur_suahist = histc(all_su_spk_times{ss},all_t_bin_edges);
    cur_suahist(all_bin_edge_pts) = [];
    cur_id_set = ismember(all_blockvec,cur_blocks);
    all_binned_sua(cur_id_set,ss) = cur_suahist(cur_id_set);
    su_probes(ss) = mode(SU_block_probes(ss,~isnan(SU_block_probes(ss,:))));
    end
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
    cur_mua_inds = find(all_clust_ids{cc} >= 1);
    
    %remove spikes from isolated SUs on the same probe from the MU
    for ss = 1:length(unique_su_nums)
        cur_mua_inds(ismember(all_spk_inds{cc}(cur_mua_inds),all_su_spk_inds{unique_su_nums(ss)})) = [];
    end

    nearby_probes = [cc-1 cc+1]; nearby_probes(nearby_probes < 1 | nearby_probes > n_probes) = [];
    cur_set = find(ismember(SU_block_probes,nearby_probes));
    if ~isempty(cur_set)
        [cur_SS,cur_BB] = ind2sub([length(SU_numbers) length(cur_used_blocks)],cur_set);
    else
        cur_SS = [];
    end
    unique_su_nums = unique(cur_SS); %set of SU numbers picked up on adjacent probes
    if ~isempty(unique_su_nums)
        double_spikes = [];
        for ss = 1:length(unique_su_nums)
            cur_blocked_inds = bsxfun(@plus,all_su_spk_inds{unique_su_nums(ss)},-double_spike_buffer:double_spike_buffer);
            double_spikes = [double_spikes; find(ismember(all_spk_inds{cc}(cur_mua_inds),cur_blocked_inds))];
        end
        fprintf('Eliminating %d of %d double spikes in MUA\n',length(double_spikes),length(cur_mua_inds));
        cur_mua_inds(double_spikes) = [];
    end
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
    all_t_axis,all_blockvec,expt_bar_ori(cur_block_set),used_inds);

[saccades,et_params] = detect_saccades(corrected_eye_vals,all_eye_speed,all_eye_ts,et_params);

par_thresh = 5;
orth_thresh = 1.25;
target_coils = [1 0];
[out_of_range] = detect_bad_fixation(corrected_eye_vals_interp,all_trialvec,used_inds,par_thresh,orth_thresh,target_coils);
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

%%

saccade_start_inds = find(ismember(used_inds,interp_sac_start_inds));
used_saccade_set = find(ismember(interp_sac_start_inds,used_inds));
%nearest index in the used data set of the saccade stop time
saccade_stop_inds = round(interp1(used_inds,1:length(used_inds),interp_sac_stop_inds(used_saccade_set)))';

saccade_stop_inds(isnan(saccade_stop_inds)) = length(used_inds);

saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));
saccade_trial_inds_ends = all_trialvec(used_inds(saccade_stop_inds));

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

%% Measure saccade magnitudes
sac_amps_inf = abs(fin_tot_corr(saccade_stop_inds) - fin_tot_corr(saccade_start_inds));

min_sac_amp = 0.1;

sac_amps = [saccades(:).amplitude];
is_micro = sac_amps(used_saccade_set) < 1 & sac_amps(used_saccade_set >= min_sac_amp);
is_macro = sac_amps(used_saccade_set) >= 1;
big_sacs = find(is_macro);
micro_sacs = find(is_micro);

big_unc_times = find(fin_tot_std(1:end-1) < 0.05 & fin_tot_std(2:end) >= 0.05);

% big_inf_sacs = find(sac_amps_inf >= prctile(sac_amps_inf,75));
% small_inf_sacs = find(sac_amps_inf <= prctile(sac_amps_inf,25));
big_inf_sacs = big_sacs(sac_amps_inf(big_sacs) >= prctile(sac_amps_inf(big_sacs),75));
small_inf_sacs = big_sacs(sac_amps_inf(big_sacs) <= prctile(sac_amps_inf(big_sacs),25));

%MAKE BIG SACS AND MICRO SACS BASED ON ORTH MOVEMENT
big_sacs = big_inf_sacs;
micro_sacs = small_inf_sacs;

Xsac = zeros(NT,n_sac_bins);
Xmsac = zeros(NT,n_sac_bins);
% Xsac_end = zeros(NT,n_sac_bins);
% Xmsac_end = zeros(NT,n_sac_bins);
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

%     cur_sac_target = saccade_start_inds(big_inf_sacs) + sac_bincents(ii);
%     uu = find(cur_sac_target > 1 & cur_sac_target < NT);
%     cur_sac_target = cur_sac_target(uu);
%     cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_trial_inds(big_inf_sacs(uu))) = [];
%     Xsac(cur_sac_target,ii) = 1;
%     
%     cur_sac_target = saccade_start_inds(small_inf_sacs) + sac_bincents(ii);
%     uu = find(cur_sac_target > 1 & cur_sac_target < NT);
%     cur_sac_target = cur_sac_target(uu);
%     cur_sac_target(all_trialvec(used_inds(cur_sac_target)) ~= saccade_trial_inds(small_inf_sacs(uu))) = [];
%     Xmsac(cur_sac_target,ii) = 1;

%     cur_sac_target = big_unc_times + sac_bincents(ii);
%     uu = find(cur_sac_target > 1 & cur_sac_target < NT);
%     cur_sac_target = cur_sac_target(uu);
%     Xmsac(cur_sac_target,ii) = 1;

end

%%
all_stimmat_cor = all_stimmat_up;
for ii = 1:NT
    all_stimmat_cor(used_inds(ii),:) = shift_matrix_Nd(all_stimmat_cor(used_inds(ii),:),-fin_tot_corr_rnd(ii),2);
end

%% randomize stimulus during saccades as a test
% sac_rand_buff = 4;
% for ii = 1:length(saccade_start_inds)
%     cur_inds = used_inds(saccade_start_inds(ii)):(used_inds(saccade_start_inds(ii))+sac_rand_buff);
%     cur_shift = round(rand*2*max_shift - max_shift);
%     all_stimmat_cor(cur_inds,:) = shift_matrix_Nd(all_stimmat_cor(cur_inds,:),cur_shift,2);
% end
%%
all_Xmat_cor = create_time_embedding(all_stimmat_cor,stim_params_us);
all_Xmat_cor = all_Xmat_cor(used_inds,use_kInds_up);

all_Xmat_uncor = create_time_embedding(all_stimmat_up,stim_params_us);
all_Xmat_uncor = all_Xmat_uncor(used_inds,use_kInds_up);

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

%% SELECT USABLE UNITS AND make Robs_mat
load(mod_data_name,'all_mod_SU*');

tr_set = et_tr_set;
full_n_chs = length(all_mod_SU);

% make Robs_mat
poss_blocks = 1:(n_blocks-1);
full_Robs_mat = nan(length(used_inds),full_n_chs);
for ss = 1:full_n_chs
    if ~isnan(all_mod_SU(ss))
    if all_mod_SU(ss) > 0
        su_probe_ind = find(SU_numbers == all_mod_SUnum(ss));
        full_Robs_mat(:,ss) = all_binned_sua(used_inds,su_probe_ind);
    else
        full_Robs_mat(:,ss) = all_binned_mua(used_inds,ss);
    end
    end
end

%%
cd(fin_anal_dir);
load(fin_mod_name);

%%
sac_stim_params(1) = NMMcreate_stim_params([2+n_blocks 1],dt);
sac_stim_params(2) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(3) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(4) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(5) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(6) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(7) = NMMcreate_stim_params([n_sac_bins 1],dt);

base_lambda_d2T = 100;
base_lambda_L2 = 2;

gain_lambda_d2T = 100;
gain_lambda_L2 = 5;

sac_BCs = [0 0 0];
sac_reg_params = NMMcreate_reg_params('lambda_d2T',base_lambda_d2T,'lambda_L2',base_lambda_L2,'boundary_conds',sac_BCs);

%%
silent = 1;
close all
clear all_sac_filters
% for ss = 1:full_n_chs
for ss = (n_probes+1):full_n_chs
% for ss = [8 23 28 36 40 50 58]
    ss
    if ~isnan(all_mod_SU(ss))
        cur_Robs = full_Robs_mat(:,ss);
        cur_tr_inds = full_inds(~isnan(cur_Robs(full_inds)));
        
        if ~isempty(cur_tr_inds)
            cur_mod = cor_gqm(ss).mod_fit;
            Xtargs = [cur_mod.mods(:).Xtarget];
            stim_NL_types = {cur_mod.mods(:).NLtype};
            stim_mod_signs = [cur_mod.mods(:).sign];
            stim_filters = [cur_mod.mods(Xtargs == 1).filtK];
            stim_outs = all_Xmat_cor*stim_filters;
            qfilts = find(strcmp(stim_NL_types,'quad'));
            tlinfilts = find(strcmp(stim_NL_types,'threshlin'));
            stim_outs(:,qfilts) = stim_outs(:,qfilts).^2;
            for jj = 1:length(tlinfilts)
                stim_outs(stim_outs(:,tlinfilts(jj)) < 0,tlinfilts(jj)) = 0;
            end
            
            exc_stim_out = sum(stim_outs(:,stim_mod_signs == 1),2);
            sup_stim_out = sum(stim_outs(:,stim_mod_signs == -1),2);
            tot_stim_out = sum(bsxfun(@times,stim_outs,stim_mod_signs),2);
            
            Xsac_excstim = bsxfun(@times,Xsac,exc_stim_out);
            Xsac_supstim = bsxfun(@times,Xsac,sup_stim_out);
            % Xsac_totstim = bsxfun(@times,Xsac,tot_stim_out);
            Xmsac_excstim = bsxfun(@times,Xmsac,exc_stim_out);
            Xmsac_supstim = bsxfun(@times,Xmsac,sup_stim_out);
            % Xmsac_totstim = bsxfun(@times,Xmsac,tot_stim_out);
            
            tr_stim{1} = [exc_stim_out(cur_tr_inds,:) sup_stim_out(cur_tr_inds,:) Xblock(cur_tr_inds,:)];
            tr_stim{2} = Xsac(cur_tr_inds,:);
            tr_stim{3} = Xmsac(cur_tr_inds,:);
            tr_stim{4} = Xsac_excstim(cur_tr_inds,:);
            tr_stim{5} = Xsac_supstim(cur_tr_inds,:);
            tr_stim{6} = Xmsac_excstim(cur_tr_inds,:);
            tr_stim{7} = Xmsac_supstim(cur_tr_inds,:);
            
            % tr_stim{1} = [tot_stim_out(cur_tr_inds,:) Xblock(cur_tr_inds,:)];
            % tr_stim{2} = Xsac(cur_tr_inds,:);
            % tr_stim{3} = Xmsac(cur_tr_inds,:);
            % tr_stim{4} = Xsac_totstim(cur_tr_inds,:);
            % tr_stim{5} = Xsac_totstim(cur_tr_inds,:);
            % tr_stim{6} = Xmsac_totstim(cur_tr_inds,:);
            % tr_stim{7} = Xmsac_totstim(cur_tr_inds,:);
            
            mod_signs = [1 1 1];
            Xtargets = [1 2 3];
            NL_types = {'lin','lin','lin','lin'};
            sac_full_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
            sac_full_mod.mods(1).reg_params = NMMcreate_reg_params();
            sac_full_mod = NMMfit_filters(sac_full_mod,cur_Robs(cur_tr_inds),tr_stim,[],[],silent);
            
            % filt_weights = [std(exc_stim_out) std(sup_stim_out)];
            filt_weights = abs(sac_full_mod.mods(1).filtK(1:2));
            filt_weights = [filt_weights; filt_weights];
            % filt_weights = [1; 1; 1; 1];
            cur_reg_params = NMMcreate_reg_params('lambda_d2T',gain_lambda_d2T*filt_weights.^2,...
                'lambda_L2',gain_lambda_L2*filt_weights.^2,'boundary_conds',[sac_BCs;sac_BCs;sac_BCs;sac_BCs]);
            sac_full_mod = NMMadd_NLinput(sac_full_mod,{'lin'},1,4,zeros(n_sac_bins,1),cur_reg_params(1));
            sac_full_mod = NMMadd_NLinput(sac_full_mod,{'lin'},-1,5,zeros(n_sac_bins,1),cur_reg_params(2));
            sac_full_mod = NMMadd_NLinput(sac_full_mod,{'lin'},1,6,zeros(n_sac_bins,1),cur_reg_params(3));
            sac_full_mod = NMMadd_NLinput(sac_full_mod,{'lin'},-1,7,zeros(n_sac_bins,1),cur_reg_params(4));
            
            sac_full_mod = NMMfit_filters(sac_full_mod,cur_Robs(cur_tr_inds),tr_stim,[],[2 3 4 5 6 7],silent);
            % sac_full_mod = NMMfit_filters(sac_full_mod,cur_Robs(cur_tr_inds),tr_stim,[],[2 3 4 5],silent);
            
            sac_filters = [sac_full_mod.mods([sac_full_mod.mods(:).Xtarget] > 1).filtK];
            all_sac_filters(ss,:,:) = sac_filters;
            
            [LL, penLL, pred_rate] = NMMmodel_eval( sac_full_mod, cur_Robs(cur_tr_inds), tr_stim);
            bsac_rel_inds = find(ismember(cur_tr_inds,saccade_start_inds(big_sacs)));
            msac_rel_inds = find(ismember(cur_tr_inds,saccade_start_inds(micro_sacs)));
            avg_sac_rate = get_event_trig_avg(pred_rate,bsac_rel_inds,-min(sac_bincents),max(sac_bincents));
            avg_msac_rate = get_event_trig_avg(pred_rate,msac_rel_inds,-min(sac_bincents),max(sac_bincents));
            %     avg_sac_rate = get_event_trig_avg(cur_Robs(cur_tr_inds),bsac_rel_inds,-min(sac_bincents),max(sac_bincents));
            %     avg_msac_rate = get_event_trig_avg(cur_Robs(cur_tr_inds),msac_rel_inds,-min(sac_bincents),max(sac_bincents));
            
    subplot(3,1,1)
            plot(sac_bincents*dt,sac_filters(:,[1 3 4]),'.-');
            subplot(3,1,2)
            plot(sac_bincents*dt,sac_filters(:,[2 5 6]),'.-');
           subplot(3,1,3)
            plot(sac_bincents*dt,avg_sac_rate,'.-',sac_bincents*dt,avg_msac_rate,'r.-');
            % figure
            % subplot(2,1,1)
            % plot(sac_bincents*dt,sac_filters(:,[1 2]));
            % subplot(2,1,2)
            % plot(sac_bincents*dt,sac_filters(:,[3 5]));
            
            NMMdisplay_model(cur_mod)
            pause
            close all
        end
    end
end

%%
sac_stim_params(1) = NMMcreate_stim_params([1 1],dt);
sac_stim_params(2) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(3) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(4) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(5) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(6) = NMMcreate_stim_params([n_sac_bins 1],dt);
sac_stim_params(7) = NMMcreate_stim_params([n_sac_bins 1],dt);

base_lambda_d2T = 50;
base_lambda_L2 = 2;

gain_lambda_d2T = 50;
gain_lambda_L2 = 2;

sac_BCs = [0 0 0];
sac_reg_params = NMMcreate_reg_params('lambda_d2T',base_lambda_d2T,'lambda_L2',base_lambda_L2,'boundary_conds',sac_BCs);

silent = 1;
close all
for ss = 1:full_n_chs
% for ss = (n_probes+1):full_n_chs
% for ss = [8 23 28 36 40 50 58]
    ss
    if ~isnan(all_mod_SU(ss))
        cur_Robs = full_Robs_mat(:,ss);
        cur_tr_inds = full_inds(~isnan(cur_Robs(full_inds)));
        Robs = cur_Robs(cur_tr_inds);
            cur_mod = dit_mods{end}(ss);
%             cur_mod = it_mods{1}(ss);
%             cur_mod = cor_gqm(ss).mod_fit;
            cur_mod.mods(4:end) = [];
        if ~isempty(cur_tr_inds) & ~isempty(cur_mod.LL_seq)
            Xtargs = [cur_mod.mods(:).Xtarget];
            stim_NL_types = {cur_mod.mods(:).NLtype};
            stim_mod_signs = [cur_mod.mods(:).sign];
            stim_filters = [cur_mod.mods(Xtargs == 1).filtK];
            stim_outs = all_Xmat_cor*stim_filters;
%             stim_outs = all_Xmat_uncor*stim_filters;
%             stim_outs(:,1) = 0;
            qfilts = find(strcmp(stim_NL_types,'quad'));
            tlinfilts = find(strcmp(stim_NL_types,'threshlin'));
            stim_outs(:,qfilts) = stim_outs(:,qfilts).^2;
            for jj = 1:length(tlinfilts)
                stim_outs(stim_outs(:,tlinfilts(jj)) < 0,tlinfilts(jj)) = 0;
            end
            
%             exc_stim_out = sum(stim_outs(:,stim_mod_signs == 1),2);
%             sup_stim_out = sum(stim_outs(:,stim_mod_signs == -1),2);
            tot_stim_out = sum(bsxfun(@times,stim_outs,stim_mod_signs),2);
            
            Xsac_totstim = bsxfun(@times,Xsac,tot_stim_out);
            Xmsac_totstim = bsxfun(@times,Xmsac,tot_stim_out);
            
            tr_stim{1} = [tot_stim_out(cur_tr_inds,:)];
            tr_stim{2} = Xsac(cur_tr_inds,:);
            tr_stim{3} = Xmsac(cur_tr_inds,:);
            tr_stim{4} = Xsac_totstim(cur_tr_inds,:);
            tr_stim{5} = Xmsac_totstim(cur_tr_inds,:);
 
            mod_signs = [1 1];
            Xtargets = [2 3];
            NL_types = {'lin','lin','lin','lin'};
            null_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
            null_mod = NMMfit_filters(null_mod,Robs,tr_stim,[],[],silent);
            [nullLL, penLL, null_pred_rate] = NMMmodel_eval( null_mod, Robs, tr_stim);
            
            mod_signs = [1 1 1];
            Xtargets = [1 2 3];
            NL_types = {'lin','lin','lin','lin'};
            sac_full_mod = NMMinitialize_model(sac_stim_params,mod_signs,NL_types,sac_reg_params,Xtargets);
            sac_full_mod.mods(1).reg_params = NMMcreate_reg_params();
            sac_full_mod = NMMfit_filters(sac_full_mod,Robs,tr_stim,[],[],silent);
            
            filt_weights = abs(sac_full_mod.mods(1).filtK(1));
            filt_weights = [filt_weights; filt_weights];
            cur_reg_params = NMMcreate_reg_params('lambda_d2T',gain_lambda_d2T*filt_weights.^2,...
                'lambda_L2',gain_lambda_L2*filt_weights.^2,'boundary_conds',[sac_BCs;sac_BCs]);
            sac_full_mod = NMMadd_NLinput(sac_full_mod,{'lin'},1,4,0.1*randn(n_sac_bins,1),cur_reg_params(1));
            sac_full_mod = NMMadd_NLinput(sac_full_mod,{'lin'},1,5,0.1*randn(n_sac_bins,1),cur_reg_params(2));
            
            dopt.optTol = 1e-10;
            dopt.progTol = 1e-10;
            sac_full_mod = NMMfit_filters(sac_full_mod,Robs,tr_stim,[],[2 3 4 5],silent);
            
            sac_filters = [sac_full_mod.mods([sac_full_mod.mods(:).Xtarget] > 1).filtK];
            
            SacMod(ss).sac_Toffset = sac_filters(:,1);
            SacMod(ss).msac_Toffset = sac_filters(:,2);
            SacMod(ss).sac_Tgain = sac_filters(:,3);
            SacMod(ss).msac_Tgain = sac_filters(:,4);
            
            [LL, penLL, pred_rate] = NMMmodel_eval( sac_full_mod,Robs, tr_stim);
            bsac_rel_inds = find(ismember(cur_tr_inds,saccade_start_inds(big_sacs)));
            msac_rel_inds = find(ismember(cur_tr_inds,saccade_start_inds(micro_sacs)));
            avg_sac_prate = get_event_trig_avg(pred_rate,bsac_rel_inds,-min(sac_bincents),max(sac_bincents));
            avg_msac_prate = get_event_trig_avg(pred_rate,msac_rel_inds,-min(sac_bincents),max(sac_bincents));
            avg_sac_rate = get_event_trig_avg(Robs,bsac_rel_inds,-min(sac_bincents),max(sac_bincents));
            avg_msac_rate = get_event_trig_avg(Robs,msac_rel_inds,-min(sac_bincents),max(sac_bincents));

            null_LLvec = Robs.*log2(null_pred_rate)-null_pred_rate;
            sacmod_LLvec = Robs.*log2(pred_rate) - pred_rate;
            avg_sac_info_rate = get_event_trig_avg(sacmod_LLvec-null_LLvec,bsac_rel_inds,-min(sac_bincents),max(sac_bincents));
            avg_sac_info_rate = jmm_smooth_1d_cor(avg_sac_info_rate,1)';
            avg_sac_rate = jmm_smooth_1d_cor(avg_sac_rate,1)';
            avg_sac_info = avg_sac_info_rate./avg_sac_rate;
            avg_msac_info_rate = get_event_trig_avg(sacmod_LLvec-null_LLvec,msac_rel_inds,-min(sac_bincents),max(sac_bincents));
            avg_msac_info_rate = jmm_smooth_1d_cor(avg_msac_info_rate,1)';
            avg_msac_rate = jmm_smooth_1d_cor(avg_msac_rate,1)';
            avg_msac_info = avg_msac_info_rate./avg_msac_rate;
            
            avg_rate = mean(Robs);
            ov_info = (LL-nullLL)/log(2);
            ov_info_rate = ov_info*avg_rate;
%             SacMod(ss).ov_info = ov_info_cont;
%             SacMod(ss).ov_info_ps = ov_info_cont_ps;
%             SacMod(ss).sac_trig_info = info_cont;
%             SacMod(ss).sac_trig_info_ps = info_cont_ps;
%             

NMMdisplay_model(cur_mod,[],[],1);

figure
subplot(2,2,1)
plot(sac_bincents*dt,sac_filters(:,[1 3]))
subplot(2,2,2)
plot(sac_bincents*dt,avg_sac_info,sac_bincents*dt,avg_msac_info,'r');
line(sac_bincents([1 end])*dt,[ov_info ov_info],'color','k')
subplot(2,2,3)
plot(sac_bincents*dt,avg_sac_prate,sac_bincents*dt,avg_msac_prate,'r');
            pause
            close all
        end
    end
end

