clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');
fig_dir = '/home/james/Analysis/bruce/ET_final/';

Expt_num = 270;

switch Expt_num
    case 266
        bar_ori = 80;
    case 270
        bar_ori = 60;
    case 275
        bar_ori = 135;
    case 277
        bar_ori = 70;
    case 281
        bar_ori = 140;
    case 287
        bar_ori = 90;
    case 289
        bar_ori = 160;
    case 294
        bar_ori = 40;
end

if Expt_num >= 275
    rpt_seed = 1001; %M275 M277 M281
else
    rpt_seed = 1e4; %m270 and 266
end


Expt_name = sprintf('M%d',Expt_num);
if Expt_num > 280 & Expt_num < 289
    data_dir = ['/media/NTlab_data2/Data/bruce/' Expt_name];
elseif Expt_num >= 289
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
else
    data_dir = ['~/Data/bruce/' Expt_name];
end
cd(data_dir);

load(sprintf('lem%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final/'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

mod_data_name = 'monoc_eyecorr_mods';
anal_name = 'monoc_eyecorr2';
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];

recompute_init_mods = 0;
use_measured_pos = 0;
use_sac_kerns = 1;
use_LOOXV = 1;
use_coils = [1 1]; %[L R]

if any(use_coils > 0)
    anal_name = [anal_name '_Cprior'];
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
elseif Expt_num == 289
    ignore_blocks = [27]; %27 is off somehow
elseif Expt_num == 294
    ignore_blocks = [37 38 39]; %37-39 have slightly different dw used in these blocks
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

spatial_usfac = 2;
%these recs have larger bar widths
if ismember(Expt_num,[287 289 294])
    use_nPix = 15;
    spatial_usfac = 4;
end

min_trial_dur = 0.75;


%%
eps = -1e3;

stim_fs = 100; %in Hz
dt = 0.01;

if Expt_num==270
    full_nPix=32;
elseif Expt_num == 287 || Expt_num == 289
    full_nPix = 22;
elseif Expt_num == 294
    full_nPix = 20;
else
    full_nPix=36;
end
stim_params = NIMcreate_stim_params([flen full_nPix],dt);
Fr = 1;

beg_buffer = 0.2;
end_buffer = 0.05;
trial_dur = 4;

n_probes = 24;

use_right_eye = false;

n_use_blocks = Inf;

use_nPix_us = use_nPix*spatial_usfac;
klen_us = use_nPix_us*flen;

sac_backlag = round(0.05/dt);
sac_forlag = round(0.3/dt);
sac_bincents = -sac_backlag:sac_forlag;
n_sac_bins = length(sac_bincents);

%%
include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaXFaXFs','rls.AllSac','rls.imiXFa'};
expt_names = cell(1,length(Expts));
expt_dds = nan(1,length(Expts));
expt_bar_ori = nan(1,length(Expts));
expt_sac_dir = nan(1,length(Expts));
expt_Fr = nan(1,length(Expts));
expt_imback = nan(1,length(Expts));
included_type = false(1,length(Expts));
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

%%
base_sp_dx = Expts{cur_block_set(2)}.Stimvals.dw; %recorded bar width (deg)
sp_dx = base_sp_dx/spatial_usfac/scale_fac; %model dx in deg

max_shift = round(0.8475/sp_dx);
dshift = 1;
shifts = -max_shift:dshift:max_shift;
n_shifts = length(shifts);
zero_frame = find(shifts==0);

max_Dshift = round(0.452/sp_dx);
Dshifts = -max_Dshift:dshift:max_Dshift;
n_Dshifts = length(Dshifts);
zero_Dframe = find(Dshifts==0);

max_Tshift = max_shift + max_Dshift;
Tshifts = -max_Tshift:dshift:max_Tshift;


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
    if buffer_pix == -1
        for ii = 1:length(left_stim_mats)
            left_stim_mats{ii} = [zeros(size(left_stim_mats{ii},1),1) left_stim_mats{ii} zeros(size(left_stim_mats{ii},1),1)];
        end
    buffer_pix = 0;
    end
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
if spatial_usfac > 1
    all_stimmat_up = zeros(size(all_stim_mat,1),full_nPix_us);
    for ii = 1:size(all_stim_mat,2)
        for jj = 1:spatial_usfac
            all_stimmat_up(:,spatial_usfac*(ii-1)+jj) = all_stim_mat(:,ii);
        end
    end
elseif spatial_usfac == 1
    all_stimmat_up = all_stim_mat;
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
    use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1) & Xinds_up(:) <= cur_use_pix(end));
    cnt = 0;
    while length(use_kInds_up)/flen < use_nPix_us
        use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1)-cnt/spatial_usfac & Xinds_up(:) <= cur_use_pix(end) + cnt/spatial_usfac);
        if length(use_kInds_up)/flen < use_nPix_us
            use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1)-(cnt+1)/spatial_usfac & Xinds_up(:) <= cur_use_pix(end) + cnt/spatial_usfac);
        end
        cnt = cnt + 1;
    end
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
        double_spikes = unique(double_spikes);
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
fix_start_inds(fix_durs==0) = [];
fix_stop_inds(fix_durs == 0) = [];
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

%%
%shift matrices for images (non-time-embedded)
n_Tshifts = length(Tshifts);
Tshift_mat = cell(n_Tshifts,1);
for xx = 1:n_Tshifts
    Tshift_mat{xx} = spdiags( ones(full_nPix_us,1), -Tshifts(xx), full_nPix_us, full_nPix_us);
end


%%
cd(anal_dir)
load(anal_name);
load(mod_data_name,'all_mod_*');
tr_set = et_tr_set;
% su_inds = find(all_mod_SU(tr_set) > 0);
su_inds = 1:length(tr_set);

%%
n_tr_chs = length(tr_set);
% make Robs_mat
poss_blocks = 1:(n_blocks-1);
Robs_mat = nan(length(used_inds),n_tr_chs);
for ss = 1:n_tr_chs
    if all_mod_SU(tr_set(ss)) > 0
        %         su_probe_ind = find(su_probes == all_mod_SU(tr_set(ss)));
        su_probe_ind = find(SU_numbers == all_mod_SUnum(tr_set(ss)));
        Robs_mat(:,ss) = all_binned_sua(used_inds,su_probe_ind);
    else
        Robs_mat(:,ss) = all_binned_mua(used_inds,tr_set(ss));
    end
end

%%
n_pts = 500;
%generate shift matrices. Must be applied to the stimulus (not the filters)
It = speye(flen);
shift_mat = cell(n_shifts,1);
for xx = 1:n_shifts
    temp = spdiags( ones(full_nPix_us,1), -shifts(xx), full_nPix_us, full_nPix_us);
    temp = kron(temp,It);
    shift_mat{xx} = temp(:,use_kInds_up);
end

use_eye_range = [-0.4 0.4];
eye_xx = linspace(use_eye_range(1),use_eye_range(2),n_pts);
srange = find(shifts*sp_dx >= use_eye_range(1) & shifts*sp_dx <= use_eye_range(2));

%%
rpt_taxis = (1:round(trial_dur)/dt)*dt-dt/2;
rpt_taxis(rpt_taxis < beg_buffer) = [];
rpt_taxis(trial_dur - rpt_taxis < end_buffer) = [];
all_trial_dur = all_trial_end_times-all_trial_start_times;

rpt_trials = find(all_trial_Se == rpt_seed);
n_rpts = length(rpt_trials);
all_rpt_inds = find(ismember(all_trialvec(used_inds),rpt_trials));

% tr_SUs = find(all_mod_SUnum(tr_set) > 0);
tr_SUs = 1:length(tr_set);
% tr_SUs = [25:27];
full_psth = nan(n_rpts,length(rpt_taxis),length(tr_SUs));
for ss = 1:length(tr_SUs)
    fprintf('Unit %d of %d\n',ss,length(tr_SUs));
%     su_ind = find(SU_numbers == all_mod_SUnum(tr_set(tr_SUs(ss))));
%     
%     cur_unit_ind = find(all_mod_SUnum==SU_numbers(su_ind));
%     cur_runit_ind = find(all_mod_SUnum(tr_set)==SU_numbers(su_ind));
%     cur_SU_ind = find(su_inds == cur_runit_ind);
    
    cur_unit_ind = tr_set(tr_SUs(ss));
    cur_runit_ind = find(tr_set == cur_unit_ind);

    fin_fix_corr = nan(NT,1);
    fin_fix_std = nan(NT,1);
    fin_fix_corr(~isnan(fix_ids)) = it_fix_post_mean_LOO(cur_runit_ind,end,fix_ids(~isnan(fix_ids)));
    fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
    fin_fix_std(~isnan(fix_ids)) = it_fix_post_std_LOO(cur_runit_ind,end,fix_ids(~isnan(fix_ids)));
    fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);
    fin_fix_corr = fin_fix_corr*sp_dx;
    fin_fix_std = fin_fix_std*sp_dx;
    fin_drift_corr = squeeze(drift_post_mean_LOO(cur_runit_ind,end,:)*sp_dx);
    fin_drift_std = squeeze(drift_post_std_LOO(cur_runit_ind,end,:)*sp_dx);
    fin_drift_corr_neural = fin_drift_corr;
    
    for ii = 1:length(trial_start_inds)
        cur_inds = trial_start_inds(ii):trial_end_inds(ii);
        fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
        fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
    end
    fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT);
    fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);
    
    fin_tot_corr = fin_fix_corr + fin_drift_corr;
    fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2);
    
%     for ii = 1:length(fix_start_inds)
%         cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
%         if length(cur_inds) > sac_shift
%             fin_drift_corr_neural(cur_inds(1:sac_shift)) = nan;
%         else
%             fin_drift_corr_neural(cur_inds) = nan;
%         end
%     end
    for ii = 1:length(fix_start_inds)-1
        cur_inds = fix_stop_inds(ii):fix_start_inds(ii+1);
        fin_drift_corr_neural(cur_inds) = nan;
    end
    fin_tot_corr_neural = fin_fix_corr + fin_drift_corr_neural';
    
%     cur_used_blocks = find(~isnan(SU_block_probes(su_ind,:)));
%     cur_rpt_trials = rpt_trials(ismember(all_trial_blk(rpt_trials),cur_used_blocks));
    all_used_rpt_inds = all_rpt_inds(~isnan(Robs_mat(all_rpt_inds,cur_runit_ind)));
    used_rpt_inds = find(ismember(all_rpt_inds,all_used_rpt_inds));
    
    for ii = 1:length(rpt_trials)
        cur_inds = find(all_trialvec(used_inds) == rpt_trials(ii));
        cur_inds(size(full_psth,2)+1:end) = [];
%         full_psth(ii,1:length(cur_inds),ss) = all_binned_sua(used_inds(cur_inds),su_ind);
        full_psth(ii,1:length(cur_inds),ss) = Robs_mat(cur_inds,cur_runit_ind);
    end
        
    all_post_cor = fin_tot_corr/sp_dx;
    all_post_cor(isnan(all_post_cor)) = 0;
    all_post_cor = round(all_post_cor) + max_Tshift + 1;
    
    %RECOMPUTE XMAT
    all_shift_stimmat_up = all_stimmat_up;
    for i=1:NT
        all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_post_cor(i)};
    end
    drift_Xmat = create_time_embedding(all_shift_stimmat_up,stim_params_us);
    drift_Xmat = drift_Xmat(used_inds,use_kInds_up);
    
    clear cor_X uncor_X
    cor_X = drift_Xmat(all_rpt_inds,:);
    uncor_X = all_Xmat_us(used_inds(all_rpt_inds),use_kInds_up);
    
    %CORRECTED AND UNCORRECTED MODEL FITS
    cur_cormod = dit_mods_spkNL_LOO{cur_runit_ind,end}(cur_unit_ind);
    cur_Xtargs = [cur_cormod.mods(:).Xtarget];
    cur_cormod.mods(cur_Xtargs > 1) = [];
%     cur_cormod = NMMfit_logexp_spkNL(cur_cormod,all_binned_sua(used_inds(all_used_rpt_inds),su_ind),cor_X(used_rpt_inds,:),[],[],[1 2]);
    cur_cormod = NMMfit_logexp_spkNL(cur_cormod,Robs_mat(all_used_rpt_inds,cur_runit_ind),cor_X(used_rpt_inds,:),[],[],[1 2]);
    
    cur_uncormod = it_mods_spkNL{1}(cur_unit_ind);
    cur_uncormod.mods(cur_Xtargs > 1) = [];
%     cur_uncormod = NMMfit_logexp_spkNL(cur_uncormod,all_binned_sua(used_inds(all_used_rpt_inds),su_ind),uncor_X(used_rpt_inds,:),[],[],[1 2]);
    cur_uncormod = NMMfit_logexp_spkNL(cur_uncormod,Robs_mat(all_used_rpt_inds,cur_runit_ind),uncor_X(used_rpt_inds,:),[],[],[1 2]);
    
%     [corrLL, penLL, cor_pred_rate(:,ss),~,~,~,nullLL] = NMMmodel_eval(cur_cormod,all_binned_sua(used_inds(all_rpt_inds),su_ind),cor_X);
%     [uncorrLL, penLL, uncor_pred_rate(:,ss),~,~,~,nullLL] = NMMmodel_eval(cur_uncormod,all_binned_sua(used_inds(all_rpt_inds),su_ind),uncor_X);    
    [corrLL, penLL, cor_pred_rate(:,ss),~,~,~,nullLL] = NMMmodel_eval(cur_cormod,Robs_mat(all_rpt_inds,cur_runit_ind),cor_X);
    [uncorrLL, penLL, uncor_pred_rate(:,ss),~,~,~,nullLL] = NMMmodel_eval(cur_uncormod,Robs_mat(all_rpt_inds,cur_runit_ind),uncor_X);    
    
    neur_eye_pos(:,:,ss) = nan(length(rpt_taxis),length(rpt_trials));
    for ii = 1:length(rpt_trials)
        test_rpt_inds = find(all_trialvec(used_inds(all_rpt_inds)) == rpt_trials(ii));
        test_rpt_inds(length(rpt_taxis)+1:end) = [];
        neur_eye_pos(1:length(test_rpt_inds),ii,ss) = fin_tot_corr_neural(all_rpt_inds(test_rpt_inds));
    end
    
    
    upts = sum(~isnan(neur_eye_pos(:,:,ss)),1);
    [~,ut] = max(upts);
    test_rpt_inds = find(all_trialvec(used_inds(all_rpt_inds)) == rpt_trials(ut));
    eyefun_pred_out = nan(length(srange),length(rpt_taxis));
    for ii = 1:length(srange)
%         fprintf('%d of %d\n',ii,length(srange));
        shift_data = uncor_X(test_rpt_inds,:)*shift_mat{srange(ii)}';
        [~, ~, temp] = NMMmodel_eval(cur_cormod,[],shift_data(:,use_kInds_up));
        eyefun_pred_out(ii,:) = temp;
    end
    eyefun_pred_out_int(:,:,ss) = interp1(-shifts(srange)*sp_dx,eyefun_pred_out,eye_xx');

end

%% Compute overall eye pos estimate
fin_fix_corr = nan(NT,1);
fin_fix_std = nan(NT,1);
fin_fix_corr(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
fin_fix_std(~isnan(fix_ids)) = it_fix_post_std(end,fix_ids(~isnan(fix_ids)));
fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);
fin_fix_corr = fin_fix_corr*sp_dx;
fin_fix_std = fin_fix_std*sp_dx;
fin_drift_corr = squeeze(drift_post_mean(end,:)*sp_dx);
fin_drift_std = squeeze(drift_post_std(end,:)*sp_dx);
fin_drift_corr_neural = fin_drift_corr;

for ii = 1:length(trial_start_inds)
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
    fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
end
fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT);
fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);

fin_tot_corr = fin_fix_corr + fin_drift_corr;
fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2);

for ii = 1:length(fix_start_inds)
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > sac_shift
        fin_drift_corr_neural(cur_inds(1:sac_shift)) = nan;
    else
        fin_drift_corr_neural(cur_inds) = nan;
    end
end

for ii = 1:length(fix_start_inds)-1
    cur_inds = fix_stop_inds(ii):fix_start_inds(ii+1);
    fin_drift_corr_neural(cur_inds) = nan;
end
fin_tot_corr_neural = fin_fix_corr + fin_drift_corr_neural;

tot_neur_eye_pos = nan(length(rpt_taxis),length(rpt_trials));
for ii = 1:length(rpt_trials)
    test_rpt_inds = find(all_trialvec(used_inds(all_rpt_inds)) == rpt_trials(ii));
    test_rpt_inds(length(rpt_taxis)+1:end) = [];
    tot_neur_eye_pos(1:length(test_rpt_inds),ii) = fin_tot_corr_neural(all_rpt_inds(test_rpt_inds));
end


%%
%COMPUTE TRIAL-BY-TRIAL PREDICTED RATES
tbt_cor_pred_rate = nan(length(rpt_taxis),length(rpt_trials),length(tr_SUs));
tbt_sac_pred_rate = nan(length(rpt_taxis),length(rpt_trials),length(tr_SUs));
tbt_uncor_pred_rate = nan(length(rpt_taxis),length(rpt_trials),length(tr_SUs));
inf_eye_pos_t = nan(length(rpt_taxis),length(rpt_trials));
row_ind = nan(length(rpt_trials),1);
for ii = 1:length(rpt_trials)
    row_ind(ii) = ii;
    test_rpt_inds = find(all_trialvec(used_inds(all_rpt_inds)) == rpt_trials(ii));
    test_rpt_inds(length(rpt_taxis)+1:end) = [];
    tbt_cor_pred_rate(1:length(test_rpt_inds),ii,:) = cor_pred_rate(test_rpt_inds,1:length(tr_SUs));
    tbt_sac_pred_rate(1:length(test_rpt_inds),ii,:) = sac_pred_rate(test_rpt_inds,1:length(tr_SUs));
    tbt_uncor_pred_rate(1:length(test_rpt_inds),ii,:) =uncor_pred_rate(test_rpt_inds,1:length(tr_SUs));    
    inf_eye_pos_t(1:length(test_rpt_inds),ii) = fin_tot_corr(all_rpt_inds(test_rpt_inds));
end

[sorted_inf_eyepos,inf_eyepos_ord] = sort(inf_eye_pos_t);


%%
use_eye_range = [-0.4 0.4];
n_pts = 500;

n_knots = 50;
knot_pts = linspace(-0.4,0.4,n_knots);

all_trials = 1:n_rpts;
xv_frac = 0.2;
xv_trials = randperm(n_rpts);
xv_trials(round(xv_frac*n_rpts)+1:end) = [];
tr_trials = setdiff(all_trials,xv_trials);
[RR,TT] = meshgrid(1:length(rpt_taxis),1:n_rpts);

%%
% eye_xx = linspace(use_eye_range(1),use_eye_range(2),n_pts);
% sp = NMMcreate_stim_params([1 n_knots]);
% spr = NMMcreate_reg_params('lambda_d2X',10,'boundary_conds',[Inf Inf Inf]);
% close all
% smoothed_eye_obs_fun = nan(length(eye_xx),length(rpt_taxis));
% np_prate = nan(length(rpt_trials),length(rpt_taxis),length(tr_SUs));
% for ss = 1:length(tr_SUs)
%     fprintf('Unit %d of %d\n',ss,length(tr_SUs));
%     for ii = 1:length(rpt_taxis)
% %         fprintf('Time %d of %d\n',ii,length(rpt_taxis));
%         xx = squeeze(neur_eye_pos(ii,:,ss))';
%         xx(xx > use_eye_range(2) | xx < use_eye_range(1)) = nan;
%         
%         yy_obs = squeeze(full_psth(:,ii,ss));
%         
%         uu = find(~isnan(xx) & ~isnan(yy_obs));
%         uu2 = uu(ismember(TT(uu,ii),tr_trials));
%         if sum(yy_obs(uu2)>0) < 1
% %             smoothed_eye_obs_fun(:,ii) = 0;
% %             cur_eyepos_p = zeros(size(uu));
%         np_prate(uu,ii,ss) = zeros(size(uu));
%         else
%             temp = fnval(csape(knot_pts,eye(length(knot_pts)),'var'),xx(uu2));
%             mod = NMMinitialize_model(sp,1,{'lin'},spr,[],[],'logexp');
%             mod = NMMfit_filters(mod,yy_obs(uu2),temp');
%             beta = mod.mods(1).filtK;
%             pfun = csape(knot_pts,beta,xx(uu),'var');
%             smoothed_eye_obs_fun(:,ii) = log(1+exp(mod.spk_NL_params(1) + fnval(pfun,eye_xx)));
%             out_bounds = find(eye_xx > max(xx(uu)) | eye_xx < min(xx(uu)));
%             smoothed_eye_obs_fun(out_bounds,ii) = nan;
%             
%             %         uu = find(~isnan(xx) & ~isnan(yy_obs));
%             cur_eyepos_p = log(1+exp(mod.spk_NL_params(1) + fnval(pfun,xx(uu))));
%             np_prate(uu,ii,ss) = cur_eyepos_p;
%             
% %             [xs,ord] = sort(xx(uu));
% %             plot(eye_xx,smoothed_eye_obs_fun(:,ii),'r'); hold on
% %             plot(xs,smooth(yy_obs(uu(ord)),10,'loess'),'b.-')
% %             plot(xs,yy_obs(uu(ord)),'k.-')
% % %             plot(eye_xx,eyefun_pred_out_int(:,ii),'g','linewidth',2);
% %             %         [xs,ord] = sort(xx(uu));
% %             %         plot(xs,yy_rate(uu(ord)),'b.-')
% %             xlim(use_eye_range);
% %             pause
% %             clf
%             
%         end
%     end
% end
% 
% ov_avg_rates = nanmean(full_psth(all_trials,:,:));
% ov_tbt_psth_mod_rate = permute(repmat(ov_avg_rates,[length(rpt_trials),1,1]),[2 1 3]);
% psth_mod = permute(ov_tbt_psth_mod_rate,[2 1 3]);
% np_prate(np_prate < psth_mod) = psth_mod(np_prate < psth_mod);

%%
[II,JJ] = meshgrid(1:n_rpts);
Tinds = [II(:) JJ(:)];

for ss = 1:length(tr_SUs);
    eps_radius = sp_dx/2;
    eye_pos_bin_edges = -0.4:sp_dx:0.4;
    in_ball_var = nan(length(rpt_taxis),1);
    in_ball_meanvar = in_ball_var;
    % in_ball_means = [];
    min_occ = 5;
    cond_mean_var = nan(length(rpt_taxis),1);
    cond_mean_invar = nan(length(rpt_taxis),1);
    cond_mean_cor = nan(length(rpt_taxis),1);
    tsub_psth_mean = nan(length(rpt_taxis),1);
    tsub_psth_var = nan(length(rpt_taxis),1);
    n_pts = nan(length(rpt_taxis),1);
    for tt = 1:length(rpt_taxis)
        cur_eye_pos = squeeze(neur_eye_pos(tt,:,ss))';
        cur_spks = squeeze(full_psth(:,tt,ss));
        
        [n,bin] = histc(cur_eye_pos,eye_pos_bin_edges);
        uset = find(n > min_occ);
        cond_mean = nan(length(eye_pos_bin_edges)-1,1);
        cond_var = nan(length(eye_pos_bin_edges)-1,1);
        for ii = 1:length(uset)
            cond_mean(uset(ii)) = nanmean(cur_spks(bin==uset(ii)));
            cond_var(uset(ii)) = nanvar(cur_spks(bin==uset(ii)));
        end
        if length(cur_eye_pos) > 10
            w = n(1:end-1)/sum(n);
            cond_mean = cond_mean-nanmean(cond_mean);
            cond_mean_var(tt) = nansum(cond_mean.^2.*w);
            cond_mean_cor(tt) = nansum(cond_mean.^2.*w./n(1:end-1));
            %    cond_mean_invar(tt) = nansum(cond_var.*n);
            %    tsub = randperm(length(cur_spks));
            %    tsub_psth_mean(tt) = nanmean(cur_spks(tsub(1:min_occ)));
            %    tsub_psth_var(tt) = nanvar(cur_spks(tsub(1:min_occ)));
        end
        n_pts(tt) = sum(~isnan(cur_eye_pos));
        %    eye_pos_dists = squareform(pdist(cur_eye_pos));
        % %    eye_pos_dists(logical(eye(n_rpts))) = nan;
        %    eye_pos_dists(II>=JJ) = nan;
        %    in_ball = find(eye_pos_dists <= eps_radius);
        %    if ~isempty(in_ball)
        %    in_ball_var(tt) = nanmean(diff(cur_spks(Tinds(in_ball,:)),[],2).^2/2);
        %    between_ball_var(tt) = nanvar(nanmean(cur_spks(Tinds(in_ball,:)),2));
        %    end
    end
    u_tpts = find(n_pts > 50);
    
    % tot_var = nanvar(reshape(full_psth(:,:,ss),[],1));
    % tot_mean = nanmean(reshape(full_psth(:,:,ss),[],1));
    psth_var = nanvar(squeeze(nanmean(full_psth(:,u_tpts,ss))));
    % rate_var = tot_var - nanmean(in_ball_var) - tot_mean;
    % em_var = rate_var - psth_var ;
    % em_var_frac = em_var/rate_var;
    
    em_var = nanmean(cond_mean_var-cond_mean_cor);
    em_var_frac(ss) = em_var/(em_var+psth_var);
end
%%
avg_rates = nanmean(full_psth(all_trials,:,:));
tr_sub_psth = bsxfun(@minus,full_psth,nanmean(full_psth,2));
tbt_psth_mod_rate = permute(repmat(avg_rates,[length(rpt_trials),1,1]),[2 1 3]);

[II,JJ] = meshgrid(1:n_rpts);
Tinds = [II(:) JJ(:)];
[a,b] = sort(tot_neur_eye_pos,2);
% b(isnan(a)) = nan;
b = b';a = a';
adiff = diff(a);
maxadiff = sp_dx;
not_usable = isnan(a);
full_psth_eyeperm1 = tr_sub_psth;
full_psth_eyeperm2 = tr_sub_psth;
full_psth_eyeperm3 = tr_sub_psth;
for ii = 1:length(rpt_taxis)
    temp = squeeze(tr_sub_psth(:,ii,:));
    temp = temp(b(:,ii),:);
    temp(not_usable(:,ii),:) = nan;
   full_psth_eyeperm1(:,ii,:) = temp;
    temp = squeeze(tr_sub_psth(:,ii,:));
    temp(1:end-1,:) = temp(b(2:end,ii),:);
    temp(not_usable(2:end,ii),:) = nan;
%     temp(adiff(:,ii) >= maxadiff,:) = nan;
   full_psth_eyeperm2(:,ii,:) = temp;
    temp = squeeze(tr_sub_psth(:,ii,:));
    temp(2:end,:) = temp(b(1:end-1,ii),:);
    temp(not_usable(1:end-1,ii),:) = nan;
   full_psth_eyeperm3(:,ii,:) = temp;
end

vec_obs_psth = reshape(permute(tr_sub_psth,[2 1 3]),[],length(tr_SUs));
vec_obs_eyeperm1 = reshape(permute(full_psth_eyeperm1,[2 1 3]),[],length(tr_SUs));
vec_obs_eyeperm2 = reshape(permute(full_psth_eyeperm2,[2 1 3]),[],length(tr_SUs));
vec_obs_eyeperm3 = reshape(permute(full_psth_eyeperm3,[2 1 3]),[],length(tr_SUs));
vec_psth_mod_rate = reshape(tbt_psth_mod_rate,[],length(tr_SUs));
vec_cor_pred_rate = reshape(tbt_cor_pred_rate,[],length(tr_SUs));
% vec_cor_pred_rate = bsxfun(@minus,vec_cor_pred_rate,nanmean(vec_cor_pred_rate));
vec_sac_pred_rate = reshape(tbt_sac_pred_rate,[],length(tr_SUs));
% vec_sac_pred_rate = bsxfun(@minus,vec_sac_pred_rate,nanmean(vec_sac_pred_rate));

maxlag = round(0.15/dt);
lags = (-maxlag:maxlag);
Cmat_obs = nan(length(tr_SUs),length(tr_SUs),2*maxlag+1);
Cmat_mod = nan(length(tr_SUs),length(tr_SUs),2*maxlag+1);
Cmat_sac = nan(length(tr_SUs),length(tr_SUs),2*maxlag+1);
Cmat_mod_res = nan(length(tr_SUs),length(tr_SUs),2*maxlag+1);
Cmat_psth_res = Cmat_obs;
Cmat_eye = Cmat_obs;
Cmat_psthmod = Cmat_obs;
Cmat_psthmod2 = Cmat_obs;
for ii = 1:length(tr_SUs)
    ii
    for jj = 1:length(tr_SUs)
        if jj > ii
            uset = ~isnan(vec_obs_psth(:,ii)) & ~isnan(vec_obs_psth(:,jj));
            uset2 = ~isnan(vec_obs_eyeperm1(:,ii)) & ~isnan(vec_obs_eyeperm2(:,jj));
            uset3 = ~isnan(vec_obs_eyeperm1(:,ii)) & ~isnan(vec_obs_eyeperm3(:,jj));
            if sum(uset) >0
                Cmat_obs(ii,jj,:) = xcov(vec_obs_psth(uset,ii),vec_obs_psth(uset,jj),maxlag,'biased');
                Cmat_mod(ii,jj,:) = xcov(vec_cor_pred_rate(uset,ii),vec_cor_pred_rate(uset,jj),maxlag,'biased');
                Cmat_sac(ii,jj,:) = xcov(vec_sac_pred_rate(uset,ii),vec_sac_pred_rate(uset,jj),maxlag,'biased');
                Cmat_psth_res(ii,jj,:) = xcov(vec_obs_psth(uset,ii) - vec_psth_mod_rate(uset,ii),...
                    vec_obs_psth(uset,jj) - vec_psth_mod_rate(uset,jj),maxlag,'biased');
                Cmat_mod_res(ii,jj,:) = xcov(vec_obs_psth(uset,ii) - vec_cor_pred_rate(uset,ii),...
                    vec_obs_psth(uset,jj) - vec_cor_pred_rate(uset,jj),maxlag,'biased');
                
                Cmat_psthmod(ii,jj,:) = xcov(vec_psth_mod_rate(uset,ii) - vec_sac_pred_rate(uset,ii),vec_psth_mod_rate(uset,jj)-vec_sac_pred_rate(uset,jj),maxlag,'biased');
                Cmat_psthmod2(ii,jj,:) = xcov(vec_psth_mod_rate(uset,ii),vec_psth_mod_rate(uset,jj),maxlag,'biased');
            end
            if sum(uset2) > 0
            temp1 = xcov(vec_obs_eyeperm1(uset2,ii),vec_obs_eyeperm2(uset2,jj),maxlag,'biased');
            else
                temp1 = nan(2*maxlag+1,1);
            end
            if sum(uset3) > 0
            temp2 = xcov(vec_obs_eyeperm1(uset3,ii),vec_obs_eyeperm3(uset3,jj),maxlag,'biased');
            else
                temp2 = nan(2*maxlag+1,1);
            end
            Cmat_eye(ii,jj,:) = nanmean([temp1 temp2],2);
        elseif jj < ii
            Cmat_obs(ii,jj,:) = Cmat_obs(jj,ii,:);
            Cmat_mod(ii,jj,:) = Cmat_mod(jj,ii,:);
            Cmat_sac(ii,jj,:) = Cmat_sac(jj,ii,:);
            Cmat_mod_res(ii,jj,:) = Cmat_mod_res(jj,ii,:);
            Cmat_eye(ii,jj,:) = Cmat_eye(jj,ii,:);
            Cmat_psth_res(ii,jj,:) = Cmat_psth_res(jj,ii,:);
            Cmat_psthmod(ii,jj,:) = Cmat_psthmod(jj,ii,:);
             Cmat_psthmod2(ii,jj,:) = Cmat_psthmod(jj,ii,:);
       end
    end
end

Cmat_eye_res = Cmat_obs - Cmat_eye;

tot_avg_rates = nanmean(vec_obs_psth);
tot_var_rates = nanvar(vec_obs_psth);
tot_sc_mat = sqrt(tot_var_rates'*tot_var_rates);
Cmat_obs = bsxfun(@rdivide,Cmat_obs,tot_sc_mat);
Cmat_mod = bsxfun(@rdivide,Cmat_mod,tot_sc_mat);
Cmat_sac = bsxfun(@rdivide,Cmat_sac,tot_sc_mat);
Cmat_mod_res = bsxfun(@rdivide,Cmat_mod_res,tot_sc_mat);
Cmat_eye = bsxfun(@rdivide,Cmat_eye,tot_sc_mat);
Cmat_eye_res = bsxfun(@rdivide,Cmat_eye_res,tot_sc_mat);
Cmat_psth_res = bsxfun(@rdivide,Cmat_psth_res,tot_sc_mat);
Cmat_psthmod = bsxfun(@rdivide,Cmat_psthmod,tot_sc_mat);
Cmat_psthmod2 = bsxfun(@rdivide,Cmat_psthmod2,tot_sc_mat);

for ii = 1:length(tr_SUs)
    neighbs = [ii-1 ii+1];
    neighbs(neighbs > length(tr_SUs) | neighbs < 1) = [];
    Cmat_obs(ii,neighbs,maxlag+1) = nan;
    Cmat_sac(ii,neighbs,maxlag+1) = nan;
    Cmat_eye_res(ii,neighbs,maxlag+1) = nan;
    Cmat_psth_res(ii,neighbs,maxlag+1) = nan;
    Cmat_mod_res(ii,neighbs,maxlag+1) = nan;
end

for ii = 1:length(tr_SUs)
    temp = squeeze(Cmat_obs(ii,:,:));
    temp2 = squeeze(Cmat_psth_res(ii,:,:));
    temp3 = squeeze(Cmat_eye_res(ii,:,:));
    beta_psth(ii) = regress(temp2(:,maxlag+1),temp(:,maxlag+1));
    beta_eye(ii) = regress(temp3(:,maxlag+1),temp(:,maxlag+1));
    
    beta_psth_sur(ii) = regress(reshape(temp2(:,[maxlag maxlag+2]),2*length(tr_SUs),[]),reshape(temp(:,[maxlag maxlag+2]),2*length(tr_SUs),[]));
    beta_eye_sur(ii) = regress(reshape(temp3(:,[maxlag maxlag+2]),2*length(tr_SUs),[]),reshape(temp(:,[maxlag maxlag+2]),2*length(tr_SUs),[]));
end

%%
% close all
figure
for ii = 1:length(tr_SUs)
   subplot(2,2,1)
   imagescnan(lags*dt,1:length(tr_SUs),squeeze(Cmat_obs(ii,:,:)));
   ca = caxis(); cam = 0.75*max(abs(ca)); caxis([-cam cam]);
   yl = ylim();
   line(lags([maxlag+1 maxlag+1])*dt,yl,'color','w');
   subplot(2,2,2)
   imagescnan(lags*dt,1:length(tr_SUs),squeeze(Cmat_psthmod(ii,:,:)));
   caxis([-cam cam]);
   line(lags([maxlag+1 maxlag+1])*dt,yl,'color','w');
   subplot(2,2,3)
   imagescnan(lags*dt,1:length(tr_SUs),squeeze(Cmat_psthmod2(ii,:,:)));
   line(lags([maxlag+1 maxlag+1])*dt,yl,'color','w');
   caxis([-cam cam]);
   subplot(2,2,4)
   imagescnan(lags*dt,1:length(tr_SUs),squeeze(Cmat_eye_res(ii,:,:)));
   caxis([-cam cam]);
   line(lags([maxlag+1 maxlag+1])*dt,yl,'color','w');
   
   pause
   clf
end
%%
udata = find(ismember(TT',xv_trials));
% udata = find(ismember(TT',all_trials));
% udata = find(ismember(TT',tr_trials));

ov_avg_rates = nanmean(full_psth(all_trials,:,:));
ov_tbt_psth_mod_rate = permute(repmat(ov_avg_rates,[length(rpt_trials),1,1]),[2 1 3]);

avg_rates = nanmean(full_psth(tr_trials,:,:));
avg_rates(avg_rates < ov_avg_rates) = ov_avg_rates(avg_rates < ov_avg_rates);
tbt_psth_mod_rate = permute(repmat(avg_rates,[length(rpt_trials),1,1]),[2 1 3]);
tbt_psth_mod_rate(tbt_psth_mod_rate < ov_tbt_psth_mod_rate) = ov_tbt_psth_mod_rate(tbt_psth_mod_rate<ov_tbt_psth_mod_rate);


tr_sub_psth = bsxfun(@rdivide,full_psth,nanmean(full_psth,2));
% vec_obs_psth = reshape(permute(tr_sub_psth,[2 1 3]),[],length(tr_SUs));
vec_obs_psth = reshape(permute(full_psth,[2 1 3]),[],length(tr_SUs));
nov_avg_rates = nanmean(vec_obs_psth);

vec_cor_pred_rate = reshape(tbt_cor_pred_rate,[],length(tr_SUs));
vec_uncor_pred_rate = reshape(tbt_uncor_pred_rate,[],length(tr_SUs));
vec_psth_mod_rate = reshape(tbt_psth_mod_rate,[],length(tr_SUs));
% vec_np_prate = reshape(permute(np_prate,[2 1 3]),[],length(tr_SUs));
vec_null_prate = repmat(nov_avg_rates,size(vec_obs_psth,1),1);

% vec_cor_pred_rate = bsxfun(@rdivide,vec_cor_pred_rate,nanmean(vec_cor_pred_rate));
% vec_uncor_pred_rate = bsxfun(@rdivide,vec_uncor_pred_rate,nanmean(vec_uncor_pred_rate));

% % np_LL = nansum(vec_obs_psth(udata,:).*log(vec_np_prate(udata,:)) - vec_np_prate(udata,:))./nansum(vec_obs_psth(udata,:));
% full_LL = nansum(vec_obs_psth(udata,:).*log(vec_obs_psth(udata,:)) - vec_obs_psth(udata,:))./nansum(vec_obs_psth(udata,:));
% cormod_LL = nansum(vec_obs_psth(udata,:).*log(vec_cor_pred_rate(udata,:)) - vec_cor_pred_rate(udata,:))./nansum(vec_obs_psth(udata,:));
% uncormod_LL = nansum(vec_obs_psth(udata,:).*log(vec_uncor_pred_rate(udata,:)) - vec_uncor_pred_rate(udata,:))./nansum(vec_obs_psth(udata,:));
% psthmod_LL = nansum(vec_obs_psth(udata,:).*log(vec_psth_mod_rate(udata,:)) - vec_psth_mod_rate(udata,:))./nansum(vec_obs_psth(udata,:));
% null_LL = nansum(vec_obs_psth(udata,:).*log(vec_null_prate(udata,:)) - vec_null_prate(udata,:))./nansum(vec_obs_psth(udata,:));
% 
% uncor_r2 = 1-(full_LL-uncormod_LL)./(full_LL-null_LL);
% cor_r2 = 1-(full_LL-cormod_LL)./(full_LL-null_LL);
% psth_r2 = 1-(full_LL-psthmod_LL)./(full_LL-null_LL);
% % fprintf('NP: %.3f\n',1-(full_LL-np_LL)./(full_LL-null_LL));
% % fprintf('Uncorr: %.3f\n',1-(full_LL-uncormod_LL)./(full_LL-null_LL));
% % fprintf('Corr: %.3f\n',1-(full_LL-cormod_LL)./(full_LL-null_LL));
% % fprintf('PSTH: %.3f\n',1-(full_LL-psthmod_LL)./(full_LL-null_LL));
% % fprintf('\n');
%%
maxlag = round(0.25/dt);
lags = -(maxlag:maxlag);
Cmat_obs = nan(length(tr_SUs),length(tr_SUs),2*maxlag+1);
Cmat_mod = Cmat_obs;
Cmat_np = Cmat_obs;
Cmat_psthmod = Cmat_obs;
Cmat_mod_res = Cmat_obs;
Cmat_np_res = Cmat_obs;
Cmat_psth_res = Cmat_obs;
for ii = 1:length(tr_SUs)
    ii
    for jj = 1:length(tr_SUs)
        if jj > ii
        uset = ~isnan(vec_obs_psth(:,ii)) & ~isnan(vec_obs_psth(:,jj));
        if sum(uset) >0
        Cmat_obs(ii,jj,:) = xcov(vec_obs_psth(uset,ii),vec_obs_psth(uset,jj),maxlag,'biased');
        Cmat_mod(ii,jj,:) = xcov(vec_cor_pred_rate(uset,ii),vec_cor_pred_rate(uset,jj),maxlag,'biased');
        Cmat_np(ii,jj,:) = xcov(vec_np_prate(uset,ii),vec_np_prate(uset,jj),maxlag,'biased');
        
        Cmat_mod_res(ii,jj,:) = xcov(vec_obs_psth(uset,ii) - vec_cor_pred_rate(uset,ii),...
            vec_obs_psth(uset,jj) - vec_cor_pred_rate(uset,jj),maxlag,'biased');
        Cmat_np_res(ii,jj,:) = xcov(vec_obs_psth(uset,ii) - vec_np_prate(uset,ii),...
            vec_obs_psth(uset,jj,:) - vec_np_prate(uset,jj),maxlag,'biased');
        Cmat_psth_res(ii,jj,:) = xcov(vec_obs_psth(uset,ii) - vec_psth_mod_rate(uset,ii),...
            vec_obs_psth(uset,jj) - vec_psth_mod_rate(uset,jj),maxlag,'biased');
        
        Cmat_psthmod(ii,jj,:) = xcov(vec_psth_mod_rate(uset,ii),vec_psth_mod_rate(uset,jj),maxlag,'biased');
        
        end
        elseif jj < ii
            Cmat_obs(ii,jj,:) = Cmat_obs(jj,ii,:);
            Cmat_mod(ii,jj,:) = Cmat_mod(jj,ii,:);
            Cmat_mod_res(ii,jj,:) = Cmat_mod_res(jj,ii,:);
            Cmat_np(ii,jj,:) = Cmat_mod(jj,ii,:);
            Cmat_np_res(ii,jj,:) = Cmat_mod_res(jj,ii,:);
            Cmat_psth_res(ii,jj,:) = Cmat_psth_res(jj,ii,:);
            Cmat_psthmod(ii,jj,:) = Cmat_psthmod(jj,ii,:);
        end
    end
end

tot_avg_rates = nanmean(vec_obs_psth);
tot_var_rates = nanvar(vec_obs_psth);
% tot_sc_mat = sqrt(tot_avg_rates'*tot_avg_rates);
tot_sc_mat = sqrt(tot_var_rates'*tot_var_rates);
Cmat_obs = bsxfun(@rdivide,Cmat_obs,tot_sc_mat);
Cmat_mod = bsxfun(@rdivide,Cmat_mod,tot_sc_mat);
Cmat_np = bsxfun(@rdivide,Cmat_np,tot_sc_mat);
Cmat_psth_res = bsxfun(@rdivide,Cmat_psth_res,tot_sc_mat);
Cmat_mod_res = bsxfun(@rdivide,Cmat_mod_res,tot_sc_mat);
Cmat_psthmod = bsxfun(@rdivide,Cmat_psthmod,tot_sc_mat);
Cmat_np_res = bsxfun(@rdivide,Cmat_np_res,tot_sc_mat);
%%
% tbt_rate_ms = bsxfun(@minus,permute(np_prate,[2 1 3]),reshape(nanmean(vec_np_prate),[1 1 length(tr_SUs)]));
% tbt_pred_diff = bsxfun(@minus,tbt_rate_ms,nanmean(tbt_rate_ms,2));
% 
% np_stimlock_var = squeeze(nanvar(nanmean(tbt_rate_ms,2)));
% np_tot_var = nanvar(reshape(tbt_rate_ms,[],length(tr_SUs)));
% np_EM_var = nanvar(reshape(tbt_pred_diff,[],length(tr_SUs)));
% np_EM_var_frac = np_EM_var./np_tot_var;

tbt_rate_ms = bsxfun(@minus,tbt_cor_pred_rate,reshape(nanmean(vec_cor_pred_rate),[1 1 length(tr_SUs)]));
tbt_pred_diff = bsxfun(@minus,tbt_rate_ms,nanmean(tbt_rate_ms,2));

cor_stimlock_var = squeeze(nanvar(nanmean(tbt_rate_ms,2)));
cor_tot_var = nanvar(reshape(tbt_rate_ms,[],length(tr_SUs)));
cor_EM_var = nanvar(reshape(tbt_pred_diff,[],length(tr_SUs)));
cor_EM_var_frac = cor_EM_var./cor_tot_var;

%%
clear res
for ss = 1:length(tr_SUs)
[res(:,ss)] = reshape(compute_poisreg_residuals(full_psth(:,:,ss),tbt_cor_pred_rate(:,:,ss)),1,[]);
[respsth(:,ss)] = reshape(compute_poisreg_residuals(full_psth(:,:,ss),tbt_psth_mod_rate(:,:,ss)),1,[]);
end
obs = reshape(full_psth,[],length(tr_SUs));
mod = reshape(tbt_cor_pred_rate,[],length(tr_SUs));
modPSTH = reshape(tbt_psth_mod_rate,[],length(tr_SUs));
%%
Cmat_res = nan(length(tr_SUs));
Cmat_resPSTH = nan(length(tr_SUs));
Cmat_obs = nan(length(tr_SUs));
Cmat_mod = nan(length(tr_SUs));
Cmat_psthmod = nan(length(tr_SUs));
Cmat_corfac = nan(length(tr_SUs));
for ii = 1:length(tr_SUs)
    for jj = 1:length(tr_SUs)
        uset = ~isnan(res(:,ii)) & ~isnan(res(:,jj));
        if sum(uset) >0
        Cmat_res(ii,jj) = corr(res(uset,ii),res(uset,jj));
        Cmat_resPSTH(ii,jj) = corr(respsth(uset,ii),respsth(uset,jj));
        Cmat_obs(ii,jj) = corr(obs(uset,ii),obs(uset,jj));
        Cmat_mod(ii,jj) = corr(mod(uset,ii),mod(uset,jj));
        
        Cmat_psthmod(ii,jj) = corr(modPSTH(uset,ii),modPSTH(uset,jj));
        
        c1 = (1+nanmean(mod(uset,ii))/nanvar(mod(uset,ii)));
        c2 = (1+nanmean(mod(uset,jj))/nanvar(mod(uset,jj)));
        Cmat_corfac(ii,jj) = 1/sqrt(c1*c2);
        end
    end
end

%%
filt_mat = [cur_cormod.mods(:).filtK];
filt_mat = reshape(filt_mat,[flen use_nPix_us 3]);
temp_kerns = squeeze(std(filt_mat,[],2));
temp_kern = mean(temp_kerns,2);
temp_kern = flipud(temp_kern) - min(temp_kern);
temp_kern = [temp_kern; zeros(length(temp_kern),1)];
temp_kern = temp_kern/sum(temp_kern);
temp_kern = flipud(temp_kern);
filtered_eyepos = inf_eye_pos_t_neu;
filtered_eyepos_vis = inf_eye_pos_t;
for ii = 1:length(cur_rpt_trials)
    cur_inds = find(~isnan(inf_eye_pos_t(ii,:)));
    %     filtered_eyepos_vis(ii,cur_inds) = conv(inf_eye_pos_t(ii,cur_inds),temp_kern,'same');
    filtered_eyepos(ii,cur_inds(sac_shift+1:end)) = inf_eye_pos_t(ii,cur_inds(1:end-sac_shift));
end

%%
n_pts = 500;
%generate shift matrices. Must be applied to the stimulus (not the filters)
It = speye(flen);
shift_mat = cell(n_shifts,1);
for xx = 1:n_shifts
    temp = spdiags( ones(full_nPix_us,1), -shifts(xx), full_nPix_us, full_nPix_us);
    temp = kron(temp,It);
    shift_mat{xx} = temp(:,use_kInds_up);
end

use_eye_range = [-0.4 0.4];
eye_xx = linspace(use_eye_range(1),use_eye_range(2),n_pts);
srange = find(shifts*sp_dx >= use_eye_range(1) & shifts*sp_dx <= use_eye_range(2));
upts = sum(~isnan(inf_eye_pos_t),2);
[~,ut] = max(upts);
test_rpt_inds = find(all_trialvec(used_inds(all_used_rpt_inds)) == cur_rpt_trials(ut));
eyefun_pred_out = nan(length(srange),length(rpt_taxis));
for ss = 1:length(srange)
    fprintf('%d of %d\n',ss,length(srange));
    shift_data = uncor_X*shift_mat{srange(ss)}';
    [~, ~, cor_pred_rate] = NMMmodel_eval(cur_cormod,[],shift_data(test_rpt_inds,use_kInds_up));
    eyefun_pred_out(ss,:) = cor_pred_rate;
end

eyefun_pred_out_int = interp1(-shifts(srange)*sp_dx,eyefun_pred_out,eye_xx');

%%
all_trials = 1:size(filtered_eyepos,1);
xv_frac = 0;
xv_trials = randperm(length(all_trials));
xv_trials(round(xv_frac*length(all_trials))+1:end) = [];
tr_trials = setdiff(all_trials,xv_trials);

use_eye_range = [-0.4 0.4];
n_pts = 500;

n_knots = 20;
knot_pts = linspace(-0.4,0.4,n_knots);

eye_xx = linspace(use_eye_range(1),use_eye_range(2),n_pts);
sp = NMMcreate_stim_params([1 n_knots]);
spr = NMMcreate_reg_params('lambda_d2X',0.1,'boundary_conds',[Inf Inf Inf]);
close all
smoothed_eye_obs_fun = nan(length(eye_xx),length(rpt_taxis));
np_prate = nan(length(all_trials),length(rpt_taxis));
for ii = 1:length(rpt_taxis)
    %     fprintf('Time %d of %d\n',ii,length(rpt_taxis));
    xx = filtered_eyepos(tr_trials,ii);
    xx(xx > use_eye_range(2) | xx < use_eye_range(1)) = nan;
    yy_rate = corr_pred_rate(tr_trials,ii);
    yy_obs = full_psth(tr_trials,ii);
    
    %     knot_pts2 = prctile(xx,linspace(1,99,n_knots));
    
    
    %     uu = find(~isnan(xx) & ~isnan(yy_obs));
    uu = find(~isnan(xx) & ~isnan(yy_rate));
    %     fprintf('%d used trials\n',length(uu));
    %     temp = fnval(fnxtr(csape(knot_pts,eye(length(knot_pts)),'var')),xx(uu));
    temp = fnval(csape(knot_pts,eye(length(knot_pts)),'var'),xx(uu));
    %     temp2 = fnval(csape(knot_pts2,eye(length(knot_pts)),'var'),xx(uu));
    if sum(yy_obs(uu)>0) < 1
        smoothed_eye_obs_fun(:,ii) = 0;
        cur_eyepos_p = zeros(size(uu));
    else
        mod = NMMinitialize_model(sp,1,{'lin'},spr,[],[],'logexp');
        mod = NMMfit_filters(mod,yy_obs(uu),temp');
        beta = mod.mods(1).filtK;
        %                 pfun = fnxtr(csape(knot_pts,beta,xx(uu),'var'));
        pfun = csape(knot_pts,beta,xx(uu),'var');
        smoothed_eye_obs_fun(:,ii) = log(1+exp(mod.spk_NL_params(1) + fnval(pfun,eye_xx)));
        
        out_bounds = find(eye_xx > max(xx(uu)) | eye_xx < min(xx(uu)));
        smoothed_eye_obs_fun(out_bounds,ii) = nan;
        
        uu = find(~isnan(xx) & ~isnan(yy_obs));
        cur_eyepos_p = log(1+exp(mod.spk_NL_params(1) + fnval(pfun,xx(uu))));
        np_prate(tr_trials(uu),ii) = cur_eyepos_p;
        %         mod = NMMinitialize_model(sp,1,{'lin'},spr,[],[],'logexp');
        %         mod = NMMfit_filters(mod,yy_obs(uu),temp2');
        %         beta2 = mod.mods(1).filtK;
        % %         pfun = fnxtr(csape(knot_pts,beta,xx(uu),'var'));
        %         pfun2 = csape(knot_pts2,beta2,xx(uu),'var');
        %         smoothed_eye_obs_fun2(:,ii) = log(1+exp(mod.spk_NL_params(1) + fnval(pfun2,eye_xx)));
        
        
        %         [xs,ord] = sort(xx(uu));
        %         plot(eye_xx,smoothed_eye_obs_fun(:,ii),'r'); hold on
        %         plot(xs,smooth(yy_obs(uu(ord)),10,'loess'),'k.-')
        %     %     plot(xs,yy_obs(uu(ord)),'k.-')
        % %         plot(eye_xx,eyefun_pred_out_int(:,ii),'g','linewidth',2);
        % %         [xs,ord] = sort(xx(uu));
        % %         plot(xs,yy_rate(uu(ord)),'b.-')
        %         xlim(use_eye_range);
        %         pause
        %         clf
        %
    end
    
    
    
end

%%
% min_rate = 1e-2;
np_pred_rate = nan(size(corr_pred_rate));
np_pred_rate_vis = nan(size(corr_pred_rate));
mod_pred_rate_vis = nan(size(corr_pred_rate));
mod_pred_rate = nan(size(corr_pred_rate));
for ii = 1:length(rpt_taxis)
    np_pred_rate(:,ii) = interp1(eye_xx,smoothed_eye_obs_fun(:,ii),filtered_eyepos(:,ii));
    np_pred_rate_vis(:,ii) = interp1(eye_xx,smoothed_eye_obs_fun(:,ii),filtered_eyepos_vis(:,ii));
    mod_pred_rate(:,ii) = interp1(eye_xx,eyefun_pred_out_int(:,ii),filtered_eyepos(:,ii));
    mod_pred_rate_vis(:,ii) = interp1(eye_xx,eyefun_pred_out_int(:,ii),filtered_eyepos_vis(:,ii));
end
% np_pred_rate(np_pred_rate < min_rate) = min_rate;
% uu = find(~isnan(np_pred_rate) & filtered_eyepos >= use_eye_range(1) & filtered_eyepos <= use_eye_range(2));
% np_LL = nansum(full_psth(uu).*log(np_pred_rate(uu)) - np_pred_rate(uu))/nansum(full_psth(uu))
uu = find(~isnan(np_prate) & filtered_eyepos >= use_eye_range(1) & filtered_eyepos <= use_eye_range(2));
np_LL = nansum(full_psth(uu).*log(np_prate(uu)) - np_prate(uu))/nansum(full_psth(uu))

fprintf('Uncorr: %.4f\n',1-(np_LL-uncorrLL)/(np_LL-nullLL));
fprintf('Corr: %.6f\n',1-(np_LL-corrLL)/(np_LL-nullLL));
fprintf('PSTH: %.4f\n',1-(np_LL-psth_LL)/(np_LL-nullLL));

%%
uset = find(~isnan(np_prate));
Robs = full_psth(uset);
np_rpred = np_prate(uset);
q_rpred = mod_pred_rate(uset);
psth_pred = repmat(nanmean(full_psth),[size(full_psth,1) 1]);
psth_pred = psth_pred(uset);

other_Robs = all_binned_mua(used_inds(full_inds(uset)),:);
other_Robs(:,su_probes(su_num)) = [];

X = [np_rpred other_Robs];
[B,dev,stats] = glmfit(X,Robs,'poisson');
X = [q_rpred other_Robs];
[B2,dev,stats2] = glmfit(X,Robs,'poisson');
X = [psth_pred other_Robs];
[B3,dev,stats2] = glmfit(X,Robs,'poisson');

%%
corr_pred_rate_ms = corr_pred_rate - nanmean(corr_pred_rate(:));
corr_pred_diff = bsxfun(@minus,corr_pred_rate_ms,nanmean(corr_pred_rate_ms));

psth_var = nanvar(nanmean(corr_pred_rate_ms));
tot_var = nanvar(corr_pred_rate_ms(:));
EM_var = nanvar(corr_pred_diff(:));

EM_var_frac = EM_var/tot_var;
