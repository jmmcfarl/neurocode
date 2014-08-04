clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');
fig_dir = '/home/james/Analysis/bruce/ET_final/';

Expt_num = 287;

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
    used_clust_set = unique(CC(SU_ID_mat==ss)); %set of clusters used to capture this SU
    SU_block_probes(ss,:) = nan(1,length(cur_block_set));
    cur_su_spk_times = [];
    cur_su_spk_inds = [];
    cur_blocks = [];
    for cc = 1:length(used_clust_set)
        cur_clust = used_clust_set(cc);
        cur_probe = SU_clust_data(cur_clust).probe_num;
        cur_clust_label = SU_clust_data(cur_clust).cluster_label;
        cur_blocks = [cur_blocks find(SU_ID_mat(:,cur_clust) == ss)];
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
su_inds = find(all_mod_SU(tr_set) > 0);

%%
rpt_taxis = (1:round(trial_dur)/dt)*dt-dt/2;
rpt_taxis(rpt_taxis < beg_buffer) = [];
rpt_taxis(trial_dur - rpt_taxis < end_buffer) = [];
all_trial_dur = all_trial_end_times-all_trial_start_times;

close all
% su_num =5; %6 for M270
% su_num = 3; %3 for M275
su_num = 7; %3 for M287

cur_unit_ind = find(all_mod_SUnum==SU_numbers(su_num));
cur_runit_ind = find(all_mod_SUnum(tr_set)==SU_numbers(su_num));
cur_SU_ind = find(su_inds == cur_runit_ind);

fin_fix_corr = nan(NT,1);
fin_fix_std = nan(NT,1);
fin_fix_corr(~isnan(fix_ids)) = it_fix_post_mean_LOO(cur_SU_ind,end,fix_ids(~isnan(fix_ids)));
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
fin_fix_std(~isnan(fix_ids)) = it_fix_post_std_LOO(cur_SU_ind,end,fix_ids(~isnan(fix_ids)));
fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids)),1:NT);
fin_fix_corr = fin_fix_corr*sp_dx;
fin_fix_std = fin_fix_std*sp_dx;
fin_drift_corr = squeeze(drift_post_mean_LOO(cur_SU_ind,end,:)*sp_dx);
fin_drift_std = squeeze(drift_post_std_LOO(cur_SU_ind,end,:)*sp_dx);
fin_drift_corr_neural = fin_drift_corr;

for ii = 1:length(trial_start_inds)
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
    fin_drift_std(cur_inds(1:end-sac_shift)) = fin_drift_std(cur_inds(sac_shift+1:end));
end
fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT);
fin_drift_std = interp1(find(~isnan(fix_ids)),fin_drift_std(~isnan(fix_ids)),1:NT);
% fin_drift_corr_neural = interp1(find(~isnan(fix_ids)),fin_drift_corr_neural(~isnan(fix_ids)),1:NT);

fin_tot_corr = fin_fix_corr + fin_drift_corr;
fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2);
% fin_tot_corr = interp1(find(~isnan(fin_tot_corr)),fin_tot_corr(~isnan(fin_tot_corr)),1:NT);
% fin_tot_std = interp1(find(~isnan(fin_tot_std)),fin_tot_std(~isnan(fin_tot_std)),1:NT);

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
fin_tot_corr_neural = fin_fix_corr + fin_drift_corr_neural';

cur_used_blocks = find(~isnan(SU_block_probes(su_num,:)));
cur_rpt_trials = find(all_trial_Se == rpt_seed & ismember(all_trial_blk,cur_used_blocks));
all_used_rpt_inds = find(ismember(all_trialvec(used_inds),cur_rpt_trials));
all_used_inds = find(ismember(all_blockvec(used_inds),cur_used_blocks));

cur_n_rpts = length(cur_rpt_trials);
full_psth = nan(cur_n_rpts,length(rpt_taxis));
all_spk_inds = [];
all_rel_times = [];
all_spk_trials = [];
all_abs_times = [];
full_inds = nan(cur_n_rpts,length(rpt_taxis));
for ii = 1:cur_n_rpts
    cur_inds = find(all_trialvec(used_inds) == cur_rpt_trials(ii));
    cur_inds(size(full_psth,2)+1:end) = [];
    full_psth(ii,1:length(cur_inds)) = all_binned_sua(used_inds(cur_inds),su_num);
    
    cur_spk_set = find(all_su_spk_times{su_num} > all_trial_start_times(cur_rpt_trials(ii)) & ...
        all_su_spk_times{su_num} < all_trial_end_times(cur_rpt_trials(ii)));
    spk_times = all_su_spk_times{su_num}(cur_spk_set);
    if ~isempty(cur_inds)
        bad = find(spk_times > all_t_axis(used_inds(cur_inds(end))) | spk_times < all_t_axis(used_inds(cur_inds(1))));
        cur_spk_set(bad) = [];
        spk_times(bad) = [];
    end
    all_spk_inds = [all_spk_inds; cur_spk_set];
    all_spk_trials = [all_spk_trials; ones(size(cur_spk_set))*ii];
    all_rel_times = [all_rel_times; spk_times-all_trial_start_times(cur_rpt_trials(ii))];
    all_abs_times = [all_abs_times; spk_times];
    full_inds(ii,1:length(cur_inds)) = cur_inds;
end

all_spk_eyepos = interp1(all_t_axis(used_inds),fin_tot_corr,all_abs_times);

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
cor_X = drift_Xmat(all_used_rpt_inds,:);
uncor_X = all_Xmat_us(used_inds(all_used_rpt_inds),use_kInds_up);

%CORRECTED AND UNCORRECTED MODEL FITS
cur_cormod = dit_mods_spkNL_LOO{cur_SU_ind,end}(cur_unit_ind);
cur_Xtargs = [cur_cormod.mods(:).Xtarget];
cur_cormod.mods(cur_Xtargs > 1) = [];
cur_cormod = NMMfit_logexp_spkNL(cur_cormod,all_binned_sua(used_inds(all_used_rpt_inds),su_num),cor_X,[],[],[1 2]);

cur_uncormod = it_mods_spkNL{1}(cur_unit_ind);
cur_uncormod.mods(cur_Xtargs > 1) = [];
cur_uncormod = NMMfit_logexp_spkNL(cur_uncormod,all_binned_sua(used_inds(all_used_rpt_inds),su_num),uncor_X,[],[],[1 2]);

[corrLL, penLL, cor_pred_rate,~,~,~,nullLL] = NMMmodel_eval(cur_cormod,all_binned_sua(used_inds(all_used_rpt_inds),su_num),cor_X);
[uncorrLL, penLL, uncor_pred_rate,~,~,~,nullLL] = NMMmodel_eval(cur_uncormod,all_binned_sua(used_inds(all_used_rpt_inds),su_num),uncor_X);
psth_pred = repmat(nanmean(full_psth),cur_n_rpts,1);
psth_LL = nansum(full_psth(:).*log(psth_pred(:)) - psth_pred(:))/nansum(full_psth(:));
full_LL = nansum(full_psth(:).*log(full_psth(:)) - full_psth(:))/nansum(full_psth(:));
fprintf('Uncorr: %.4f\n',1-(full_LL-uncorrLL)/(full_LL-nullLL));
fprintf('Corr: %.6f\n',1-(full_LL-corrLL)/(full_LL-nullLL));
fprintf('PSTH: %.4f\n',1-(full_LL-psth_LL)/(full_LL-nullLL));

%COMPUTE TRIAL-BY-TRIAL PREDICTED RATES
corr_pred_rate = nan(length(cur_rpt_trials),length(rpt_taxis));
uncorr_pred_rate = nan(length(cur_rpt_trials),length(rpt_taxis));
inf_eye_pos_t = nan(length(cur_rpt_trials),length(rpt_taxis));
inf_eye_pos_t_neu = nan(length(cur_rpt_trials),length(rpt_taxis));
row_ind = nan(length(cur_rpt_trials),1);
for ii = 1:length(cur_rpt_trials)
    row_ind(ii) = ii;
    test_rpt_inds = find(all_trialvec(used_inds(all_used_rpt_inds)) == cur_rpt_trials(ii));
    test_rpt_inds(length(rpt_taxis)+1:end) = [];
    [~, ~, corr_pred_rate(ii,1:length(test_rpt_inds))] = NMMmodel_eval(cur_cormod,all_binned_sua(used_inds(all_used_rpt_inds(test_rpt_inds)),su_num),cor_X(test_rpt_inds,:));
    [~, ~, uncorr_pred_rate(ii,1:length(test_rpt_inds))] = NMMmodel_eval(cur_uncormod,all_binned_sua(used_inds(all_used_rpt_inds(test_rpt_inds)),su_num),uncor_X(test_rpt_inds,:));
    
    inf_eye_pos_t(ii,1:length(test_rpt_inds)) = fin_tot_corr(all_used_rpt_inds(test_rpt_inds));
    inf_eye_pos_t_neu(ii,1:length(test_rpt_inds)) = fin_tot_corr_neural(all_used_rpt_inds(test_rpt_inds));
end

[sorted_inf_eyepos,inf_eyepos_ord] = sort(inf_eye_pos_t);

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

%%
np_pred_rate_ms = np_pred_rate - nanmean(np_pred_rate(:));
np_pred_diff = bsxfun(@minus,np_pred_rate_ms,nanmean(np_pred_rate_ms));

psth_var = nanvar(nanmean(np_pred_rate_ms));
tot_var = nanvar(np_pred_rate_ms(:));
EM_var = nanvar(np_pred_diff(:));

EM_var_frac_np = EM_var/tot_var;

%%
uu = 1:length(rpt_taxis);
cur_utrials = find(~any(isnan(np_pred_rate_vis(:,uu)),2));

mod_prediction = full_psth(cur_utrials,uu);
% mod_prediction = np_pred_rate_vis(cur_utrials,uu);
% mod_prediction = mod_pred_rate(utrials,:);
% mod_prediction = mod_pred_rate_vis(cur_utrials,uu);

Trials = full_psth(cur_utrials,uu);
Trials(isnan(mod_prediction)) = nan;

pred = nanmean(mod_prediction);
nt = sum(~isnan(mod_prediction),2);

Resp = nanmean(Trials);
NTrials = size(Trials,1);

mTrials = nanmean(Trials,2);
mPred = nanmean(pred);
VarResp = nanvar(Resp);

VarError = nanvar(Resp - pred);

varRespTrial = zeros(NTrials,1);
varErrorTrial = zeros(NTrials,1);
for ii = 1:NTrials
    varErrorTrial(ii) = nansum((Trials(ii,:) - pred - mTrials(ii) + mPred).^2)./(nt(ii)-1);
    varRespTrial(ii) = nansum((Trials(ii,:) - mTrials(ii)).^2)./(nt(ii)-1);
end
VarErrorInd = nanmean(varErrorTrial);
VarRespInd = nanmean(varRespTrial);
VarSignal = 1/(NTrials-1) * (NTrials*VarResp - VarRespInd);
VarNoise = VarRespInd - VarSignal;
Beta = (VarResp - VarError)/VarSignal;

%%
fig_dir = '/home/james/Analysis/bruce/WIP/';

close all
utrials = find(~any(isnan(full_psth),2));
xr = [1 3];
ep = find(rpt_taxis >= xr(2),1);
sp = find(rpt_taxis >= xr(1),1);

bad_trials = 194;
utrials(ismember(utrials,bad_trials)) = [];

h = figure;
hold on;
plot(rpt_taxis(sp:ep),nanmean(full_psth(utrials,sp:ep)/dt),'k','linewidth',2);
plot(rpt_taxis(sp:ep),nanmean(np_pred_rate(utrials,sp:ep)/dt),'r','linewidth',2);
plot(rpt_taxis(sp:ep),nanmean(mod_pred_rate(utrials,sp:ep)/dt),'b','linewidth',2);
xlim(xr);
xlabel('Time (s)');
ylabel('Firing rate (Hz)');
yl = ylim();
ylim([0 yl(2)]);

fig_width = 4.6; %3.27 4.86 6.83
rel_height = 0.7;
figufy(h);
% fname = [fig_dir 'examp_psth2.pdf'];
% exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h);

%%

xr = [1 3];

n_utrials = length(utrials);
lw = 0.5;
theight = 0.75;
% theight = 0.95;

h1 = figure(); hold on
for nn = 1:length(utrials)
    cur_spikes = find(all_spk_trials == utrials(nn));
    for ii = 1:length(cur_spikes)
        line(all_rel_times(cur_spikes([ii ii])),nn + [0 theight],'color','k','linewidth',lw);
    end
end
xlim(xr);
ylim([1 n_utrials+1]);
xlabel('Time (s)');
ylabel('Trial');

fig_width = 4.6; %3.27 4.86 6.83
rel_height = 0.8;
figufy(h1);
% fname = [fig_dir 'examp_rpt_raster.pdf'];
% exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(h1);

%%

h2 = figure();
imagesc(rpt_taxis(sp:ep),1:length(utrials),np_pred_rate_vis(utrials,sp:ep)/dt)
xlabel('Time (s)');
ylabel('Trial');
% colorbar
caxis([0 1]/dt);
colormap(gray);

fig_width = 4.6; %3.27 4.86 6.83
rel_height = 0.8;
figufy(h2);
fname = [fig_dir 'examp_rpt_prates3.pdf'];
exportfig(h2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h2);

%%
psth_pred_rate = repmat(nanmean(full_psth(utrials,:)),length(utrials),1);
h = figure();
imagesc(rpt_taxis(sp:ep),1:length(utrials),psth_pred_rate(:,sp:ep)/dt)
xlabel('Time (s)');
ylabel('Trial');
% colorbar
caxis([0 1]/dt);
colormap(gray);


fig_width = 4.6; %3.27 4.86 6.83
rel_height = 0.8;
figufy(h);
fname = [fig_dir 'psth_fullmat.pdf'];
exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h);

%%
h = figure();
imagesc(rpt_taxis(sp:ep),1:length(utrials),mod_pred_rate_vis(utrials,sp:ep)/dt)
xlabel('Time (s)');
ylabel('Trial');
% colorbar
caxis([0 1]/dt);
colormap(gray);


fig_width = 4.6; %3.27 4.86 6.83
rel_height = 0.8;
figufy(h);
fname = [fig_dir 'modpred_fullmat.pdf'];
exportfig(h,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h);

%%
% examp_t = 180;
% examp_t = 192;
examp_t = 205;
h3 = figure();
plot(eye_xx,smoothed_eye_obs_fun(:,examp_t))
hold on
k = ksdensity(filtered_eyepos_vis(:),eye_xx);
k = k/max(k);
plot(eye_xx,k,'k')
xlim([-0.3 0.3]);
%%
xr = [1.5 2.5];
ep = find(rpt_taxis >= xr(2),1);
sp = find(rpt_taxis >= xr(1),1);

fig_width = 3.5; %3.27 4.86 6.83
rel_height = 0.8;

t1 = 12;
t2 = 30;
% for t2 = 13:100
h5 = figure();
%     subplot(2,1,1)
plot(rpt_taxis(sp:ep),inf_eye_pos_t(utrials(t1),sp:ep));
hold on
plot(rpt_taxis(sp:ep),inf_eye_pos_t(utrials(t2),sp:ep),'r');
xlim(rpt_taxis([sp ep]));
ylim([-0.25 0.25])
xlabel('Time (s)');
ylabel('Eye position (deg)');

h6 = figure();
%     subplot(2,1,2)
plot(rpt_taxis(sp:ep),np_pred_rate_vis(utrials(t1),sp:ep)/dt);
hold on
plot(rpt_taxis(sp:ep),np_pred_rate_vis(utrials(t2),sp:ep)/dt,'r');
xlim(rpt_taxis([sp ep]));
xlabel('Time (s)');
ylabel('Predicted rate (Hz)');
%     t2
%     pause
%     clf
% end

figufy(h5);
fname = [fig_dir 'examp_trial_EM.pdf'];
exportfig(h5,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h5);

figufy(h6);
fname = [fig_dir 'examp_trial_rates2.pdf'];
exportfig(h6,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h6);

%%
hh = figure();
imagesc(rpt_taxis(sp:ep),eye_xx,smoothed_eye_obs_fun(:,sp:ep));
hold on
ylim([-0.3 0.3]);
caxis([0 1]);
colormap(gray);
% plot(rpt_taxis(sp:ep),inf_eye_pos_t(utrials(t1),sp:ep));
% hold on
% plot(rpt_taxis(sp:ep),inf_eye_pos_t(utrials(t2),sp:ep),'r');
% plot(rpt_taxis(sp:ep),inf_eye_pos_t_neu(utrials(t1),sp:ep));
% hold on
% plot(rpt_taxis(sp:ep),inf_eye_pos_t_neu(utrials(t2),sp:ep),'r');
plot(rpt_taxis(sp:ep),filtered_eyepos_vis(utrials(t1),sp:ep));
hold on
plot(rpt_taxis(sp:ep),filtered_eyepos_vis(utrials(t2),sp:ep),'r');
xlabel('Time (s)');
ylabel('Eye position (deg)');
figufy(hh);


hh2 = figure();
ksdensity(inf_eye_pos_t(:),eye_xx);
xlim([-0.3 0.3])
xlabel('Eye position (deg)');
ylabel('Probability');
figufy(hh2);

figufy(hh);
fname = [fig_dir 'spline_pred.pdf'];
exportfig(hh,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(hh);

figufy(hh2);
fname = [fig_dir 'EM_dist.pdf'];
exportfig(hh2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(hh2);

%%
[I,J] = meshgrid(1:length(utrials));
temp = corr(corr_pred_rate(utrials,:)','type','spearman');
temp(I <= J) = nan;
avg_spear_corr = nanmean(temp(:));


corr_pred_rate_ms = corr_pred_rate - nanmean(corr_pred_rate(:));
corr_pred_diff = bsxfun(@minus,corr_pred_rate_ms,nanmean(corr_pred_rate_ms));

np_rate_ms = np_pred_rate - nanmean(np_pred_rate(:));
np_diff = bsxfun(@minus,np_pred_rate,nanmean(np_pred_rate));


cor_psth_var = nanvar(nanmean(corr_pred_rate_ms));
cor_tot_var = nanvar(corr_pred_rate_ms(:));
cor_EM_var = nanvar(corr_pred_diff(:));
cor_EM_var_frac = cor_EM_var/cor_tot_var;

np_psth_var = nanvar(nanmean(np_rate_ms));
np_tot_var = nanvar(np_rate_ms(:));
np_EM_var = nanvar(np_diff(:));
np_EM_var_frac = np_EM_var/np_tot_var;

%%
base_Xmat = all_Xmat_us(used_inds,:);
disc_fin_tot_corr = round(fin_tot_corr/sp_dx);
rate_out = nan(n_shifts,NT);
for ii = 1:n_shifts
    fprintf('Shift %d of %d\n',ii,n_shifts);
    shiftX = base_Xmat*shift_mat{ii};
    [~,~,rate_out(ii,:)] = NMMmodel_eval(cur_cormod,[],shiftX);
end
uu = find(disc_fin_tot_corr >= -max_shift & disc_fin_tot_corr <= max_shift);
ep_dist = hist(disc_fin_tot_corr(uu),shifts);
ep_dist = ep_dist/sum(ep_dist);

cond_mrate = sum(bsxfun(@times,rate_out,ep_dist'));
rate_out_ms = bsxfun(@minus,rate_out,cond_mrate);
cond_var = sum(bsxfun(@times,rate_out_ms.^2,ep_dist'));

[~,~,base_rate_out] = NMMmodel_eval(cur_cormod,[],base_Xmat(:,use_kInds_up));

%%
n_tr_set = length(tr_set);
n_squared_filts = 1;
all_mods = dit_mods{end}(tr_set);
Xtargets = [all_mods(1).mods(:).Xtarget];
filt_bank = zeros(n_tr_set,klen_us,n_squared_filts+1);
lin_kerns = nan(n_tr_set,n_blocks);
mod_spkNL_params = nan(n_tr_set,3);
for ss = 1:n_tr_set
    cur_Xtargs = [all_mods(ss).mods(:).Xtarget];
    cur_k = [all_mods(ss).mods(cur_Xtargs == 1).filtK];
    n_used_filts = size(cur_k,2);
    filt_bank(ss,:,1:n_used_filts) = cur_k;
    mod_spkNL_params(ss,:) = all_mods(ss).spk_NL_params;
    lin_kerns(ss,:) = all_mods(ss).mods(cur_Xtargs == 2).filtK;
end
filt_bank = permute(filt_bank,[2 1 3]);

mod_spkNL_params(:,1) = mod_spkNL_params(:,1) + mean(lin_kerns,2);

use_NT = round(NT/5);
base_Xmat = all_Xmat_us(used_inds(1:use_NT),:);

disc_fin_tot_corr = round(fin_tot_corr/sp_dx);
full_rate_out = nan(n_shifts,use_NT,n_tr_set);
for ii = 1:n_shifts
    fprintf('Shift %d of %d\n',ii,n_shifts);
    shiftX = base_Xmat*shift_mat{ii};
    
    %outputs of stimulus models at current X-matrix shift
    gfuns = ones(use_NT,n_tr_set);
    gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,1)');
    gfuns = gfuns + shiftX*squeeze(filt_bank(:,:,1));
    for ff = 2:(n_squared_filts+1)
        gfuns = gfuns + (shiftX*squeeze(filt_bank(:,:,ff))).^2;
    end
    
    %incorporate beta
    gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,2)');
    
    %handle numerical overflow with log(1+exp)
    too_large = gfuns > 50;
    pred_rate = log(1+exp(gfuns));
    pred_rate(too_large) = gfuns(too_large);
    
    %incorporate alpha
    pred_rate = bsxfun(@times,pred_rate,mod_spkNL_params(:,3)');
    
    %enforce min predicted rate
    pred_rate(pred_rate < 1e-50) = 1e-50;
    full_rate_out(ii,:,:) = pred_rate;
end

cond_mrate = (sum(bsxfun(@times,full_rate_out,ep_dist')));
rate_out_ms = bsxfun(@minus,full_rate_out,cond_mrate);
cond_mrate = squeeze(cond_mrate);
full_cond_var = nan(n_tr_set,n_tr_set,use_NT);
for nn = 1:n_tr_set
    nn
    for jj = (nn):n_tr_set
        cond_var = sum(bsxfun(@times,rate_out_ms(:,:,nn).*rate_out_ms(:,:,jj),ep_dist'));
        full_cond_var(nn,jj,:) = cond_var;
    end
end
full_EM_var = mean(full_cond_var,3);
aa = diag(full_EM_var);
full_EM_corr = full_EM_var./sqrt(aa*aa');

base_rate_out = squeeze(full_rate_out(max_shift+1,:,:));
full_base_cvar = cov(base_rate_out);
full_base_mean = mean(base_rate_out);
full_base_var = var(base_rate_out);

full_EM_corr2 = full_EM_var./sqrt(full_base_var'*full_base_var);


%%
full_base_mean_r = full_base_mean/dt;
full_base_var_r = full_base_var/dt^2;

T = dt;
sc_facs = T + full_base_mean_r./full_base_var_r;
sc_fac_mat = sc_facs'*sc_facs;
sc_fac_mat = T./sqrt(sc_fac_mat);

pred_corr_mat = full_EM_corr2.*sc_fac_mat;

%%
all_trials = 1:size(filtered_eyepos,1);
[RR,TT] = meshgrid(1:size(corr_pred_rate,2),1:size(corr_pred_rate,1));
xv_frac = 0.5;
xv_trials = randperm(length(all_trials));
xv_trials(round(xv_frac*length(all_trials))+1:end) = [];
tr_trials = setdiff(all_trials,xv_trials);

use_eye_range = [-0.4 0.4];
n_pts = 500;

n_knots = 15;
knot_pts = linspace(-0.4,0.4,n_knots);

eye_xx = linspace(use_eye_range(1),use_eye_range(2),n_pts);
sp = NMMcreate_stim_params([1 n_knots]);
spr = NMMcreate_reg_params('lambda_d2X',0.5,'boundary_conds',[Inf Inf Inf]);
close all
smoothed_eye_obs_fun = nan(length(eye_xx),length(rpt_taxis));
smoothed_eye_obs_fun2 = nan(length(eye_xx),length(rpt_taxis));
for ii = 1:length(rpt_taxis)
    %     fprintf('Time %d of %d\n',ii,length(rpt_taxis));
    xx = filtered_eyepos(all_trials,ii);
    xx(xx > use_eye_range(2) | xx < use_eye_range(1)) = nan;
    yy_rate = corr_pred_rate(all_trials,ii);
    yy_obs = full_psth(all_trials,ii);
    
    uu = find(~isnan(xx) & ~isnan(yy_rate));
    uu2 = uu(ismember(TT(uu),tr_trials));
    
    temp = fnval(csape(knot_pts,eye(length(knot_pts)),'var'),xx(uu2));
    if sum(yy_obs(uu2)>0) < 1
        smoothed_eye_obs_fun(:,ii) = 0;
        smoothed_eye_obs_fun2(:,ii) = 0;
        cur_eyepos_p = zeros(size(uu2));
    else
        mod = NMMinitialize_model(sp,1,{'lin'},spr,[],[],'logexp');
        mod = NMMfit_filters(mod,yy_obs(uu2),temp');
        beta = mod.mods(1).filtK;
        %         pfun = fnxtr(csape(knot_pts,beta,xx(uu),'var'));
        pfun = csape(knot_pts,beta,xx(uu),'var');
        smoothed_eye_obs_fun(:,ii) = log(1+exp(mod.spk_NL_params(1) + fnval(pfun,eye_xx)));
        
        temp = fnval(csape(knot_pts,eye(length(knot_pts)),'var'),xx(uu));
        mod = NMMinitialize_model(sp,1,{'lin'},spr,[],[],'logexp');
        mod = NMMfit_filters(mod,yy_obs(uu),temp');
        beta = mod.mods(1).filtK;
        %         pfun = fnxtr(csape(knot_pts,beta,xx(uu),'var'));
        pfun = csape(knot_pts,beta,xx(uu),'var');
        smoothed_eye_obs_fun2(:,ii) = log(1+exp(mod.spk_NL_params(1) + fnval(pfun,eye_xx)));
        
        cur_eyepos_p = log(1+exp(mod.spk_NL_params(1) + fnval(pfun,xx(uu))));
    end
    
    if ~isempty(uu)
        [xs,ord] = sort(xx(uu));
        plot(eye_xx,smoothed_eye_obs_fun(:,ii),'r'); hold on
        plot(eye_xx,smoothed_eye_obs_fun2(:,ii),'g'); hold on
        plot(xs,smooth(yy_obs(uu(ord)),10,'loess'),'k.-')
        %     plot(xs,yy_obs(uu(ord)),'k.-')
        %     plot(eye_xx,eyefun_pred_out_int(:,ii),'g','linewidth',2);
        %     [xs,ord] = sort(xx(uu));
        %     plot(xs,yy_rate(uu(ord)),'b.-')
        xlim(use_eye_range);
        pause
        clf
    end
end

%%
[RR,TT] = meshgrid(1:size(corr_pred_rate,2),1:size(corr_pred_rate,1));

% min_rate = 1e-2;
np_pred_rate = nan(size(corr_pred_rate));
np_pred_rate2 = nan(size(corr_pred_rate));
mod_pred_rate = nan(size(corr_pred_rate));
for ii = 1:length(rpt_taxis)
    np_pred_rate(:,ii) = interp1(eye_xx,smoothed_eye_obs_fun(:,ii),filtered_eyepos(:,ii));
    np_pred_rate2(:,ii) = interp1(eye_xx,smoothed_eye_obs_fun2(:,ii),filtered_eyepos(:,ii));
    mod_pred_rate(:,ii) = interp1(eye_xx,eyefun_pred_out_int(:,ii),filtered_eyepos(:,ii));
end
% np_pred_rate(np_pred_rate < min_rate) = min_rate;
% uu = find(~isnan(np_pred_rate) & filtered_eyepos >= use_eye_range(1) & filtered_eyepos <= use_eye_range(2));
uu = find(ismember(TT,tr_trials) & ~isnan(np_pred_rate) & filtered_eyepos >= use_eye_range(1) & filtered_eyepos <= use_eye_range(2));
np_LL = nansum(full_psth(uu).*log(np_pred_rate(uu)) - np_pred_rate(uu))/nansum(full_psth(uu))
np_LL2 = nansum(full_psth(uu).*log(np_pred_rate2(uu)) - np_pred_rate2(uu))/nansum(full_psth(uu))
corrLL = nansum(full_psth(uu).*log(corr_pred_rate(uu)) - corr_pred_rate(uu))/nansum(full_psth(uu))
uncorrLL = nansum(full_psth(uu).*log(uncorr_pred_rate(uu)) - uncorr_pred_rate(uu))/nansum(full_psth(uu))

fprintf('Uncorr: %.4f\n',1-(np_LL-uncorrLL)/(np_LL-nullLL));
fprintf('Corr: %.6f\n',1-(np_LL-corrLL)/(np_LL-nullLL));
fprintf('PSTH: %.4f\n',1-(np_LL-psth_LL)/(np_LL-nullLL));
fprintf('Sub NP: %.4f\n',1-(np_LL-np_LL2)/(np_LL-nullLL));

%%
n_rpts = 20;
poss_xv_fracs = [0.1 0.25 0.5 0.75 0.9];
for pp = 1:length(poss_xv_fracs)
    xv_frac = poss_xv_fracs(pp);
    for rr = 1:n_rpts
        rset = randperm(length(rpt_taxis));
        rset = rset(1:round(length(rpt_taxis)*xv_frac));
        uu = find(ismember(RR,rset) & ~isnan(np_pred_rate) & filtered_eyepos >= use_eye_range(1) & filtered_eyepos <= use_eye_range(2));
        np_LL = nansum(full_psth(uu).*log(np_pred_rate(uu)) - np_pred_rate(uu))/nansum(full_psth(uu))
        corrLL = nansum(full_psth(uu).*log(corr_pred_rate(uu)) - corr_pred_rate(uu))/nansum(full_psth(uu))
        uncorrLL = nansum(full_psth(uu).*log(uncorr_pred_rate(uu)) - uncorr_pred_rate(uu))/nansum(full_psth(uu))
        
        uncorr_(rr,pp) = 1-(np_LL-uncorrLL)/(np_LL-nullLL);
        corr_(rr,pp) = 1-(np_LL-corrLL)/(np_LL-nullLL);
        psth_(rr,pp) = 1-(np_LL-psth_LL)/(np_LL-nullLL);
    end
end


%%
use_eye_range = [-0.4 0.4];
n_pts = 500;

poss_smooth = [0.01 0.1 0.5 1 5 10 50];
        eye_xx = linspace(use_eye_range(1),use_eye_range(2),n_pts);

all_trials = 1:size(filtered_eyepos,1);
[RR,TT] = meshgrid(1:size(corr_pred_rate,2),1:size(corr_pred_rate,1));
xv_frac = 0.25;

min_prate = 1e-5;

n_rpts = 5;
for rr = 1:n_rpts
    xv_trials = randperm(length(all_trials));
    xv_trials(round(xv_frac*length(all_trials))+1:end) = [];
    tr_trials = setdiff(all_trials,xv_trials);
    
    n_knots = 15;
    knot_pts = linspace(-0.4,0.4,n_knots);
    
    for ss = 1:length(poss_smooth)
        smooth_lambda = poss_smooth(ss);
        
        sp = NMMcreate_stim_params([1 n_knots]);
        spr = NMMcreate_reg_params('lambda_d2X',smooth_lambda,'boundary_conds',[Inf Inf Inf]);
        
        close all
        smoothed_eye_obs_fun = nan(length(eye_xx),length(rpt_taxis));
        eval_eye_obs_fun = nan(length(all_trials),length(rpt_taxis));
        for ii = 1:length(rpt_taxis)
            %     fprintf('Time %d of %d\n',ii,length(rpt_taxis));
            xx = filtered_eyepos(tr_trials,ii);
            xx(xx > use_eye_range(2) | xx < use_eye_range(1)) = nan;
            yy_rate = corr_pred_rate(tr_trials,ii);
            yy_obs = full_psth(tr_trials,ii);
            
            uu = find(~isnan(xx) & ~isnan(yy_rate));
            
            temp = fnval(csape(knot_pts,eye(length(knot_pts)),'var'),xx(uu));
            if sum(yy_obs(uu)>0) < 1
                smoothed_eye_obs_fun(:,ii) = 0;
                cur_eyepos_p = zeros(size(uu));
            else
                mod = NMMinitialize_model(sp,1,{'lin'},spr,[],[],'logexp');
                mod = NMMfit_filters(mod,yy_obs(uu),temp');
                beta = mod.mods(1).filtK;
                %         pfun = fnxtr(csape(knot_pts,beta,xx(uu),'var'));
                pfun = csape(knot_pts,beta,xx(uu),'var');
                smoothed_eye_obs_fun(:,ii) = log(1+exp(mod.spk_NL_params(1) + fnval(pfun,eye_xx)));
                
                xx = filtered_eyepos(xv_trials,ii);
                uu = find(xx >= use_eye_range(1) & xx <= use_eye_range(2));
                cur_eyepos_p = log(1+exp(mod.spk_NL_params(1) + fnval(pfun,xx(uu))));
                eval_eye_obs_fun(xv_trials(uu),ii) = cur_eyepos_p;
            end
            
        end
        eval_eye_obs_fun(eval_eye_obs_fun < min_prate) = min_prate;
        uu = find(~isnan(eval_eye_obs_fun));
        np_xvLL(rr,ss) = nansum(full_psth(uu).*log(eval_eye_obs_fun(uu)) - eval_eye_obs_fun(uu))/nansum(full_psth(uu))
        
    end
end