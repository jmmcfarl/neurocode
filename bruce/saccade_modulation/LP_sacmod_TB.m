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

% %MAKE BIG SACS AND MICRO SACS BASED ON ORTH MOVEMENT
% big_sacs = big_inf_sacs;
% micro_sacs = small_inf_sacs;

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
time_to_gsac = nan(length(all_t_axis),1);
time_to_msac = nan(length(all_t_axis),1);
max_dist = round(0.4/dt);
rel_tax = (-max_dist:max_dist)';
for ii = 1:length(big_sacs)
    cur_inds = (used_inds(saccade_start_inds(big_sacs(ii))) - max_dist):(used_inds(saccade_start_inds(big_sacs(ii))) + max_dist);
    cur_uset = find(cur_inds > 1 & cur_inds < length(all_t_axis));
    cur_uset(abs(time_to_gsac(cur_inds(cur_uset))) < abs(rel_tax(cur_uset))) = [];
    time_to_gsac(cur_inds(cur_uset)) = rel_tax(cur_uset);
end
for ii = 1:length(micro_sacs)
    cur_inds = (used_inds(saccade_start_inds(micro_sacs(ii))) - max_dist):(used_inds(saccade_start_inds(micro_sacs(ii))) + max_dist);
    cur_uset = find(cur_inds > 1 & cur_inds < length(all_t_axis));
    cur_uset(abs(time_to_gsac(cur_inds(cur_uset))) < abs(rel_tax(cur_uset))) = [];
    time_to_msac(cur_inds(cur_uset)) = rel_tax(cur_uset);
end
time_to_gsac = time_to_gsac*dt;
time_to_msac = time_to_msac*dt;

%%
back_t = 0.2;
for_t = 0.4;
n_Gbins = 35;
Xtick = -(back_t-dt/2):(dt):(for_t+dt/2);
n_sbins = length(Xtick);
addpath('~/James_scripts/TentBasis2D/');
%%
silent = 1;
close all
clear all_sac_filters
% for ss = 1:full_n_chs
for ss = (n_probes+1):full_n_chs
% % for ss = [8 23 28 36 40 50 58]
    ss
    if ~isnan(all_mod_SU(ss))
        cur_Robs = full_Robs_mat(:,ss);
        cur_tr_inds = full_inds(~isnan(cur_Robs(full_inds)));
        
        if ~isempty(cur_tr_inds)
            cur_mod = dit_mods{end}(ss);
%             cur_mod = it_mods{end}(ss);
            Xtargs = [cur_mod.mods(:).Xtarget];
            stim_NL_types = {cur_mod.mods(:).NLtype};
            stim_mod_signs = [cur_mod.mods(Xtargs==1).sign];
            stim_filters = [cur_mod.mods(Xtargs == 1).filtK];
            stim_outs = all_Xmat_cor*stim_filters;
%             stim_outs = all_Xmat_uncor*stim_filters;
            qfilts = find(strcmp(stim_NL_types,'quad'));
            tlinfilts = find(strcmp(stim_NL_types,'threshlin'));
            stim_outs(:,qfilts) = stim_outs(:,qfilts).^2;
            for jj = 1:length(tlinfilts)
                stim_outs(stim_outs(:,tlinfilts(jj)) < 0,tlinfilts(jj)) = 0;
            end
            G = sum(bsxfun(@times,stim_outs,stim_mod_signs),2);
            
            TB_stim = [time_to_gsac(used_inds(cur_tr_inds)) G];
            Ytick = linspace(my_prctile(TB_stim(:,2),0.1),my_prctile(TB_stim(:,2),100-1),n_Gbins);
            TB = TentBasis2D(Xtick, Ytick);
            
            used_data = find(TB_stim(:,1) >= Xtick(1) & TB_stim(:,1) <= Xtick(end) & ...
                TB_stim(:,2) >= Ytick(1) & TB_stim(:,2) <= Ytick(end));
            [TB_Xmat,TB_counts] = TB.InputNL2D(TB_stim(used_data,:));
            L2_params = create_L2_params([],[1 n_sbins*n_Gbins],[n_sbins n_Gbins],2,3,[Inf Inf],[0.1 1]);
            TB_fitmod = regGLM_fit(TB_Xmat,cur_Robs(used_data),L2_params,10,[],[],1);
            [LL, penLL, pred_rate, G] = regGLM_eval(TB_fitmod,cur_Robs(used_data),TB_Xmat);
            TB_K = reshape(TB_fitmod.K,n_sbins,n_Gbins)';
            bin_areas = TB.GetBinAreas();
            gsac_TB_dist = TB_counts./bin_areas;
            gsac_TB_dist = gsac_TB_dist'/sum(gsac_TB_dist(:));
            gsac_TB_rate = log(1 + exp(TB_K + TB_fitmod.theta));
            
            %INFO CALS
            cur_avg_rate = mean(cur_Robs(used_data));
            marg_gdist = sum(gsac_TB_dist,2);
            marg_sdist = sum(gsac_TB_dist);
            marg_gsacrate = sum(gsac_TB_dist.*gsac_TB_rate)./marg_sdist;
            marg_grate = sum(gsac_TB_dist.*gsac_TB_rate,2)./marg_gdist;
            ov_info_gsac = sum(marg_gdist.*marg_grate.*log2(marg_grate/cur_avg_rate));
        ov_info_gsac_ps = ov_info_gsac/cur_avg_rate;

        gsacdep_info = nan(n_sbins,1);
            for tt = 1:n_sbins
                gsacdep_info(tt) = sum(gsac_TB_dist(:,tt).*gsac_TB_rate(:,tt).*log2(gsac_TB_rate(:,tt)/marg_gsacrate(tt)))/sum(gsac_TB_dist(:,tt));
            end
        gsacdep_info_ps = gsacdep_info./marg_gsacrate';
          
                figure
        imagesc(Xtick,Ytick,gsac_TB_rate);
        set(gca,'ydir','normal');

        figure
        subplot(2,1,1);
        plot(Xtick,marg_gsacrate)
        line(Xtick([1 end]),[cur_avg_rate cur_avg_rate],'color','k');
        subplot(2,1,2);
        plot(Xtick,gsacdep_info);
        hold on
        plot(Xtick,gsacdep_info_ps,'r');
        line(Xtick([1 end]),[ov_info_gsac ov_info_gsac],'color','k');
        line(Xtick([1 end]),[ov_info_gsac_ps ov_info_gsac_ps],'color','g');
        end
    end
end

