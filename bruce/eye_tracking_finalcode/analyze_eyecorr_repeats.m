clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');
fig_dir = '/home/james/Analysis/bruce/ET_final/';

Expt_num = 275;

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
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

mod_data_name = 'monoc_eyecorr_mods';
anal_name = 'monoc_eyecorr';

recompute_init_mods = 0;
use_measured_pos = 0;
use_sac_kerns = 1;
use_LOOXV = 1;
use_coils = [1 1]; %[L R]

if any(use_coils > 0)
   anal_name = [anal_name '_Cpriortight'];
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

sp_dx = 0.0565/spatial_usfac/scale_fac;
max_shift = round(15*spatial_usfac/scale_fac);
dshift = 1;
shifts = -max_shift:dshift:max_shift;
n_shifts = length(shifts);
zero_frame = find(shifts==0);
temp_dx = [-2*max_shift:dshift:2*max_shift];
shift_dx_edges = [temp_dx*sp_dx-dshift*sp_dx/2 temp_dx(end)*sp_dx+dshift*sp_dx/2];
ovshift_dx_edges = [shifts*sp_dx-dshift*sp_dx/2 shifts(end)*sp_dx+dshift*sp_dx/2];

max_Dshift = round(8*spatial_usfac/scale_fac);
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
% include_expts = {'rls.Fa', 'rls.FaXimi','rls.FaXFaXFs','rls.AllSac'};
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
%shift matrices for images (non-time-embedded)    
n_Tshifts = length(Tshifts);
Tshift_mat = cell(n_Tshifts,1);
for xx = 1:n_Tshifts
    Tshift_mat{xx} = spdiags( ones(full_nPix_us,1), -Tshifts(xx), full_nPix_us, full_nPix_us);
end

%%
cd(anal_dir)
anal_name = 'monoc_eyecorr_Cprior';
% anal_name = 'monoc_eyecorr';
load(anal_name);
mod_data_name = 'monoc_eyecorr_mods';
load(mod_data_name,'all_mod_*');
tr_set = et_tr_set;
su_inds = find(all_mod_SU(tr_set) > 0);

%%
rpt_taxis = (1:round(trial_dur)/dt)*dt-dt/2;
rpt_taxis(rpt_taxis < beg_buffer) = [];
rpt_taxis(trial_dur - rpt_taxis < end_buffer) = [];
all_trial_dur = all_trial_end_times-all_trial_start_times;

close all
% su_num =6; %for M270
su_num = 3;

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
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > sac_shift
        fin_drift_corr(cur_inds(1:end-sac_shift+1)) = fin_drift_corr(cur_inds(sac_shift:end));
        fin_drift_std(cur_inds(1:end-sac_shift+1)) = fin_drift_std(cur_inds(sac_shift:end));
    end
end

fin_tot_corr = fin_fix_corr' + fin_drift_corr;
fin_tot_std = sqrt(fin_fix_std'.^2 + fin_drift_std.^2);
fin_tot_corr = interp1(find(~isnan(fin_tot_corr)),fin_tot_corr(~isnan(fin_tot_corr)),1:NT);
fin_tot_std = interp1(find(~isnan(fin_tot_std)),fin_tot_std(~isnan(fin_tot_std)),1:NT);

cur_used_blocks = find(~isnan(SU_block_probes(su_num,:)));
cur_rpt_trials = find(all_trial_Se == rpt_seed & ismember(all_trial_blk,cur_used_blocks));
all_used_rpt_inds = find(ismember(all_trialvec(used_inds),cur_rpt_trials));
all_used_inds = find(ismember(all_blockvec(used_inds),cur_used_blocks));

cur_n_rpts = length(cur_rpt_trials);
full_psth = nan(cur_n_rpts,length(rpt_taxis));
y_cnt = 0; y_space = 1;
for ii = 1:cur_n_rpts
    cur_inds = find(all_trialvec(used_inds) == cur_rpt_trials(ii));
    cur_inds(size(full_psth,2)+1:end) = [];
    full_psth(ii,1:length(cur_inds)) = all_binned_sua(used_inds(cur_inds),su_num);
    
    cur_spk_set = find(all_su_spk_times{su_num} > all_trial_start_times(cur_rpt_trials(ii)) & ...
        all_su_spk_times{su_num} < all_trial_end_times(cur_rpt_trials(ii)));
    spk_times = all_su_spk_times{su_num}(cur_spk_set) - all_trial_start_times(cur_rpt_trials(ii));
    
    plot(spk_times,ones(size(spk_times))*y_cnt,'k.')
    hold on
    y_cnt = y_cnt + y_space;
end

%NOW CONSTRUCT CORRECTED STIMULUS
fix_post_cor = nan(NT,1);
fix_post_cor(~isnan(fix_ids)) = it_fix_post_mean_LOO(cur_SU_ind,end,fix_ids(~isnan(fix_ids)));
fix_post_cor = interp1(find(~isnan(fix_ids)),fix_post_cor(~isnan(fix_ids)),1:NT);
fix_post_cor(isnan(fix_post_cor)) = 0;
%back project drift (within-fixation) by sac_shift
drift_post_cor = squeeze(drift_post_mean_LOO(cur_SU_ind,end,:));
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > sac_shift
        drift_post_cor(cur_inds(1:end-sac_shift+1)) = drift_post_cor(cur_inds(sac_shift:end));
    end
end
drift_post_cor(isnan(drift_post_cor)) = 0;

all_post_cor = round((fix_post_cor'+drift_post_cor)) + max_Tshift + 1;

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
row_ind = nan(length(cur_rpt_trials),1);
for ii = 1:length(cur_rpt_trials)
    row_ind(ii) = ii;
    test_rpt_inds = find(all_trialvec(used_inds(all_used_rpt_inds)) == cur_rpt_trials(ii));
    %     if length(test_rpt_inds) >= length(rpt_taxis)
    test_rpt_inds(length(rpt_taxis)+1:end) = [];
    [~, ~, corr_pred_rate(ii,1:length(test_rpt_inds))] = NMMmodel_eval(cur_cormod,all_binned_sua(used_inds(all_used_rpt_inds(test_rpt_inds)),su_num),cor_X(test_rpt_inds,:));
    [~, ~, uncorr_pred_rate(ii,1:length(test_rpt_inds))] = NMMmodel_eval(cur_uncormod,all_binned_sua(used_inds(all_used_rpt_inds(test_rpt_inds)),su_num),uncor_X(test_rpt_inds,:));
    %     end
    
    inf_eye_pos_t(ii,1:length(test_rpt_inds)) = fin_tot_corr(all_used_rpt_inds(test_rpt_inds));
end

bad_rows = find(any(isnan(corr_pred_rate),2));
% bad_rows = find(any(inf_eye_pos_t <= -0.25,2));
% uu_inds = find(fin_tot_corr(all_used_rpt_inds) > -0.25);
psth = full_psth;
avg_psth = nanmean(psth);

corr_pred_rate(bad_rows,:) = [];
uncorr_pred_rate(bad_rows,:) = [];
psth(bad_rows,:) = [];
inf_eye_pos_t(bad_rows,:) = [];
row_ind(bad_rows) = [];

psth_sm_sigma = 1;
smooth_psth = psth;
for jj = 1:size(psth,1)
    smooth_psth(jj,:) = jmm_smooth_1d_cor(psth(jj,:),psth_sm_sigma);
end

n_full_rpts = size(psth,1);

rnull_rate = mean(psth);
rnull_pred_rate = repmat(rnull_rate,size(psth,1),1);
[pseudo_R2,mod_dev,null_dev] = pseudo_r2(psth,corr_pred_rate,rnull_pred_rate);

%%
figure
plot(rpt_taxis,pseudo_R2)

%% PLOT STACK OF SINGLE-TRIAL OBSERVED AND PREDICTED RATES
xr = [0.6 1.6]; %1-2 2-24
t1 = 2;
t2 = 24;
fig_width = 4.86; %3.27 4.86 6.83
rel_height = 2;

h=figure;
% subplot(2,1,1)
% imagesc(rpt_taxis,1:n_full_rpts,psth);
% %c=colorbar; ylabel(c,'Spike count');
% % caxis([0 5]);
% xlim(xr);
% xlabel('Time (s)');
% ylabel('Trial number');
% title('Observed spiking');

% subplot(2,1,2)
imagesc(rpt_taxis,1:n_full_rpts,corr_pred_rate/dt);
caxis([0 200]);
%c=colorbar; ylabel(c,'Spike rate');
xlim(xr);
xlabel('Time (s)');
ylabel('Trial number');
title('Predicted with eye correction');

% subplot(3,1,3)
% % imagesc(rpt_taxis,1:n_full_rpts,inf_eye_pos_t);
% plot(rpt_taxis,inf_eye_pos_t,'k','linewidth',0.5);
% hold on
% plot(rpt_taxis,inf_eye_pos_t(t1,:),'linewidth',2);
% plot(rpt_taxis,inf_eye_pos_t(t2,:),'r','linewidth',2);
% ylim([-0.41 0.41]);
% % c=colorbar; 
% % ylabel(c,'Eye position (deg)');
% xlim(xr);
% xlabel('Time (s)');
% ylabel('Trial number');
% title('Inferred eye position');

figufy(h);

% fname = [fig_dir 'example_repeat_rates.pdf'];
% exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close;


%%
% figure
% imagesc(rpt_taxis,1:n_full_rpts,inf_eye_pos_t);

%% PLOT EXAMPLE PREDICTED RATES FOR DIFFERENT TRIALS, ALONG WITH INFERRED EYE POSITION
% close all

fig_width = 4.86; %3.27 4.86 6.83
rel_height = 1.75;

% xr = [1 2];
xr = [0.6 1.6]; %1-2 2-24
% xr = [1 2.3];
yr = [0 250];
t1 = 2;
t2 = 24;
h=figure;

% for t2 = 14:29;


subplot(3,1,1)
plot(rpt_taxis,inf_eye_pos_t(t1,:),'linewidth',1)
hold on
plot(rpt_taxis,inf_eye_pos_t(t2,:),'r','linewidth',1)
xlim(xr);
ylim([-0.3 0.3]);
line(xr,[0 0],'color','k','linestyle','--');
xlabel('Time (s)');
ylabel('Inferred eye position (deg)')
% pause
% clf
% end

subplot(3,1,2)
plot(rpt_taxis,corr_pred_rate(t1,:)/dt,'linewidth',1);
hold on
plot(rpt_taxis,corr_pred_rate(t2,:)/dt,'r','linewidth',1)
xlim(xr);
ylim(yr);
xlabel('Time (s)');
ylabel('Predicted rate (Hz)')

subplot(3,1,3); hold on
plot(rpt_taxis,avg_psth/dt,'k','linewidth',1)
plot(rpt_taxis,mean(corr_pred_rate)/dt,'m','linewidth',1)
legend('Observed trial-avg','Predicted trial-avg');
xlim(xr);
ylim(yr);
xlabel('Time (s)');
ylabel('Predicted rate (Hz)')

% pause
% clf
% end

% figufy(h);
% 
% fname = [fig_dir 'example_repeat_trials.pdf'];
% exportfig(fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close;
% 

%%
all_rate_dists = [];
all_eye_dists = [];
for tt = 1:length(rpt_taxis);
   rate_dists = pdist(corr_pred_rate(:,tt));
   eye_dists = pdist(inf_eye_pos_t(:,tt));
   dist_corr(tt) = corr(rate_dists',eye_dists','type','spearman');
   all_rate_dists = [all_rate_dists rate_dists/avg_psth(tt)];
   all_eye_dists = [all_eye_dists eye_dists];
end


%% COMPARE MODELS BEFORE AND AFTER EYE CORRS
close all
fig_width = 4.86; %3.27 4.86 6.83

n_stim_filts = 3;
cor_filts = reshape([cur_cormod.mods(1:n_stim_filts).filtK],[flen use_nPix_us n_stim_filts]);
uncor_filts = reshape([cur_uncormod.mods(1:n_stim_filts).filtK],[flen use_nPix_us n_stim_filts]);

[~, ~, ~, ~, cor_filt_out] = NMMmodel_eval(cur_cormod,[],drift_Xmat(all_used_inds,:));
[~, ~, ~, ~, uncor_filt_out] = NMMmodel_eval(cur_uncormod,[],all_Xmat_us(used_inds(all_used_inds),use_kInds_up));


rf_pos = Expts{cur_block_set(1)}.Stimvals.rf(1:2)/scale_fac;
rf_proj = rf_pos(:,1)*cos((bar_ori-90)*pi/180) + rf_pos(:,2)*sin((bar_ori-90)*pi/180);
pix_ax = ((1:use_nPix_us)-use_nPix_us/2-0.5)*sp_dx + rf_proj;

lag_ax = (0:(flen-1))*dt*1e3;

h = figure();
n_cols = 3;
for ii = 1:n_stim_filts
   subplot(n_stim_filts,n_cols,(ii-1)*n_cols+1);
   imagesc(pix_ax,lag_ax,squeeze(cor_filts(:,:,ii)));
   set(gca,'ydir','normal');
%   xlim([3.05 3.8]);
  xlabel('Position (deg)'); ylabel('Lag (ms)');
  
   subplot(n_stim_filts,n_cols,(ii-1)*n_cols+2);
   imagesc(pix_ax,lag_ax,squeeze(uncor_filts(:,:,ii)));
   set(gca,'ydir','normal');
%    xlim([3.05 3.8]);
  xlabel('Position (deg)'); ylabel('Lag (ms)');

  subplot(n_stim_filts,n_cols,(ii-1)*n_cols+3);
   xr = prctile(cor_filt_out(:,ii),[0.1 99.9]);
   gendist_x = linspace(xr(1),xr(2),100);
   [gendist_y] = ksdensity(cor_filt_out(:,ii),gendist_x);
   [gendist_y2] = ksdensity(uncor_filt_out(:,ii),gendist_x);
   plot(gendist_x,gendist_y,gendist_x,gendist_y2,'r');
   xlim(xr);

end
colormap(gray);
subplot(n_stim_filts,n_cols,1);
title('Corrected');
subplot(n_stim_filts,n_cols,2);
title('Uncorrected');
% subplot(n_stim_filts,n_cols,3);
% legend('Corrected','Uncorrected');

figufy(h);

fname = [fig_dir 'example_repeat_cormod.pdf'];
exportfig(fname,'width',fig_width,'height',0.8*fig_width,'fontmode','scaled','fontsize',1);
close;

%%
rpt_taxis = (1:round(trial_dur)/dt)*dt-dt/2;
rpt_taxis(rpt_taxis < beg_buffer) = [];
rpt_taxis(trial_dur - rpt_taxis < end_buffer) = [];
all_trial_dur = all_trial_end_times-all_trial_start_times;

close all
% su_num =6; %6 for M270
su_num = 3; %3 for M275

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
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > sac_shift
        fin_drift_corr(cur_inds(1:end-sac_shift+1)) = fin_drift_corr(cur_inds(sac_shift:end));
        fin_drift_std(cur_inds(1:end-sac_shift+1)) = fin_drift_std(cur_inds(sac_shift:end));
    end
end

fin_tot_corr = fin_fix_corr' + fin_drift_corr;
fin_tot_std = sqrt(fin_fix_std'.^2 + fin_drift_std.^2);
fin_tot_corr = interp1(find(~isnan(fin_tot_corr)),fin_tot_corr(~isnan(fin_tot_corr)),1:NT);
fin_tot_std = interp1(find(~isnan(fin_tot_std)),fin_tot_std(~isnan(fin_tot_std)),1:NT);

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
    
end

all_spk_eyepos = interp1(all_t_axis(used_inds),fin_tot_corr,all_abs_times);


%NOW CONSTRUCT CORRECTED STIMULUS
fix_post_cor = nan(NT,1);
fix_post_cor(~isnan(fix_ids)) = it_fix_post_mean_LOO(cur_SU_ind,end,fix_ids(~isnan(fix_ids)));
fix_post_cor = interp1(find(~isnan(fix_ids)),fix_post_cor(~isnan(fix_ids)),1:NT);
fix_post_cor(isnan(fix_post_cor)) = 0;
%back project drift (within-fixation) by sac_shift
drift_post_cor = squeeze(drift_post_mean_LOO(cur_SU_ind,end,:));
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > sac_shift
        drift_post_cor(cur_inds(1:end-sac_shift+1)) = drift_post_cor(cur_inds(sac_shift:end));
    end
end
drift_post_cor(isnan(drift_post_cor)) = 0;

all_post_cor = round((fix_post_cor'+drift_post_cor)) + max_Tshift + 1;

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
row_ind = nan(length(cur_rpt_trials),1);
for ii = 1:length(cur_rpt_trials)
    row_ind(ii) = ii;
    test_rpt_inds = find(all_trialvec(used_inds(all_used_rpt_inds)) == cur_rpt_trials(ii));
    test_rpt_inds(length(rpt_taxis)+1:end) = [];
    [~, ~, corr_pred_rate(ii,1:length(test_rpt_inds))] = NMMmodel_eval(cur_cormod,all_binned_sua(used_inds(all_used_rpt_inds(test_rpt_inds)),su_num),cor_X(test_rpt_inds,:));
    [~, ~, uncorr_pred_rate(ii,1:length(test_rpt_inds))] = NMMmodel_eval(cur_uncormod,all_binned_sua(used_inds(all_used_rpt_inds(test_rpt_inds)),su_num),uncor_X(test_rpt_inds,:));
    
    inf_eye_pos_t(ii,1:length(test_rpt_inds)) = fin_tot_corr(all_used_rpt_inds(test_rpt_inds));
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
filtered_eyepos = inf_eye_pos_t;
for ii = 1:length(cur_rpt_trials)
    cur_inds = find(~isnan(inf_eye_pos_t(ii,:)));
    filtered_eyepos(ii,cur_inds) = conv(inf_eye_pos_t(ii,cur_inds),temp_kern,'same');
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
test_rpt_inds = find(all_trialvec(used_inds(all_used_rpt_inds)) == cur_rpt_trials(1));
eyefun_pred_out = nan(length(srange),length(rpt_taxis));
for ss = 1:length(srange)
    fprintf('%d of %d\n',ss,length(srange));
   shift_data = uncor_X*shift_mat{srange(ss)}';
   [~, ~, cor_pred_rate] = NMMmodel_eval(cur_cormod,[],shift_data(test_rpt_inds,use_kInds_up));
   eyefun_pred_out(ss,:) = cor_pred_rate;
end

eyefun_pred_out_int = interp1(-shifts(srange)*sp_dx,eyefun_pred_out,eye_xx');

%%
n_use_trials = size(corr_pred_rate,1);

use_eye_range = [-0.4 0.4];
n_knots = 15;
spline_order = 4;

eye_xx = linspace(use_eye_range(1),use_eye_range(2),n_pts);
sp = NMMcreate_stim_params([1 n_knots]);
spr = NMMcreate_reg_params('lambda_d2X',0.1);

knot_pts = linspace(-0.4,0.4,n_knots);
% Bspline_mat = bspline_basismatrix(spline_order,knot_pts,eye_xx);

close all
smoothed_eye_rate_fun = nan(n_pts,length(rpt_taxis));
smoothed_eye_obs_fun = nan(n_pts,length(rpt_taxis));

eye_density = ksdensity(filtered_eyepos(:),eye_xx);

for ii = 1:length(rpt_taxis)
    ii
    xx = filtered_eyepos(:,ii);
    yy_rate = corr_pred_rate(:,ii);
    yy_obs = full_psth(:,ii);
    
    % knot_pts = [-0.4*ones(1,2) knot_pts 0.4*ones(1,2)];
    uu = find(~isnan(xx) & ~isnan(yy_obs));
    
    temp = fnval(fnxtr(csape(knot_pts,eye(length(knot_pts)),'var')),xx(uu));
    if sum(yy_obs(uu)>0) < 5
        smoothed_eye_obs_fun(:,ii) = 0;
    else        
        mod = NMMinitialize_model(sp,1,{'lin'},spr,[],[],'logexp');
        mod = NMMfit_filters(mod,yy_obs(uu),temp');
        beta = mod.mods(1).filtK;
        pfun = fnxtr(csape(knot_pts,beta,xx(uu),'var'));
        smoothed_eye_obs_fun(:,ii) = log(1+exp(mod.spk_NL_params(1) + fnval(pfun,eye_xx)));
    end
    
    mod = NMMinitialize_model(sp,1,{'lin'},spr,[],[],'linear');
    mod = NMMfit_filters(mod,yy_rate(uu),temp');
    beta = mod.mods(1).filtK;
    pfun = fnxtr(csape(knot_pts,beta,xx(uu),'var'));
    smoothed_eye_rate_fun(:,ii) = mod.spk_NL_params(1) + fnval(pfun,eye_xx);
        
%     plot(xx,yy_obs,'.'); hold on
%     plot(eye_xx,smoothed_eye_obs_fun(:,ii),'r');
%     [xs,ord] = sort(xx(uu));
%     plot(xs,smooth(yy_obs(uu(ord)),50,'loess'),'k')
%     pause
%     clf
    
end

%%
% bm_pred_rate = nan(size(corr_pred_rate));
% for ii = 1:length(rpt_taxis)
%    bm_pred_rate(:,ii) = interp1(eye_xx,smoothed_eye_obs_fun(:,ii),filtered_eyepos(:,ii)); 
% end
% bm_pred_rate(bm_pred_rate < 1e-3) = 1e-3;
% bm_LL = nansum(full_psth(:).*log(bm_pred_rate(:)) - bm_pred_rate(:))/nansum(full_psth(:));
% 
% fprintf('Uncorr: %.4f\n',1-(bm_LL-uncorrLL)/(bm_LL-nullLL));
% fprintf('Corr: %.6f\n',1-(bm_LL-corrLL)/(bm_LL-nullLL));
% fprintf('PSTH: %.4f\n',1-(bm_LL-psth_LL)/(bm_LL-nullLL));

%%
close all
xr = [1 2];

ep = find(rpt_taxis > xr(2),1);
utrials = find(~any(isnan(full_psth(:,1:ep)),2));
n_utrials = length(utrials);


lw = 0.5;
theight = 0.95;

f1 = figure(); hold on
for nn = 1:length(utrials)
    cur_spikes = find(all_spk_trials == utrials(nn));
    for ii = 1:length(cur_spikes)
        line(all_rel_times(cur_spikes([ii ii])),nn + [0 theight],'color','k','linewidth',lw);
    end    
end
xlim(xr);
ylim([0 n_utrials]);
xlabel('Time (s)');
ylabel('Trial');

%%
fig_width = 3.27; 
rel_height = 0.8;

figufy(f1);
fname = [fig_dir 'examp_rpt_raster.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

%%


h1 = figure();
imagescnan(rpt_taxis,1:n_utrials,corr_pred_rate(utrials,:)/dt);
caxis([0 0.8]/dt);
%c=colorbar; ylabel(c,'Spike rate');
xlim(xr);
xlabel('Time (s)');
ylabel('Trial number');
title('Predicted with eye correction');
ylim([85 150])
%%
% close all
%for M275SU3 [50 136]
t1 = 94; %52
t2 = 136; %92
h2 = figure(); 
% for t2 = 100:n_utrials;
plot(rpt_taxis,corr_pred_rate(utrials(t1),:)/dt,'r');hold on
plot(rpt_taxis,corr_pred_rate(utrials(t2),:)/dt,'b');
xlim(xr);
xlabel('Time (s)');
ylabel('Predicted rate (Hz)');

% pause
% clf
% end


h3 = figure(); 
% plot(rpt_taxis,filtered_eyepos(utrials(t1),:),'r');hold on
% plot(rpt_taxis,filtered_eyepos(utrials(t2),:),'b');
plot(rpt_taxis,inf_eye_pos_t(utrials(t1),:),'r');hold on
plot(rpt_taxis,inf_eye_pos_t(utrials(t2),:),'b');
xlim(xr);
xlabel('Time (s)');
ylabel('Eye position (deg)');
yl = ylim(); ym = max(abs(yl)); ylim([-ym ym]);
line(xr,[0 0],'color','k','linestyle','--');

%%
fig_width = 4.6; %3.27 4.86 6.83
rel_height = 0.8;

figufy(h1);
fname = [fig_dir 'examp_rpt_predstack.pdf'];
exportfig(h1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h1);

figufy(h2);
fname = [fig_dir 'examp_rpt_predtraces.pdf'];
exportfig(h2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h2);

figufy(h3);
fname = [fig_dir 'examp_rpt_eyetraces.pdf'];
exportfig(h3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(h3);

%%

f1 = figure();
subplot(2,1,1);
imagesc(rpt_taxis,eye_xx,smoothed_eye_obs_fun/dt);
set(gca,'ydir','normal');
caxis([0 1]/dt);
xlim(xr);
xlabel('Time (s)');
ylabel('Eye position (deg)');
title('Nonparametric');
subplot(2,1,2);
imagesc(rpt_taxis,eye_xx,eyefun_pred_out_int/dt);
set(gca,'ydir','normal');
caxis([0 1]/dt);
xlim(xr);
xlabel('Time (s)');
ylabel('Eye position (deg)');
title('Model');

f2 = figure;
plot(eye_xx,eye_density);
xlabel('Eye position (deg)');
ylabel('Probability density');

f3 = figure;
plot(rpt_taxis,nanmean(full_psth)/dt,'k');
hold on
plot(rpt_taxis,nanmean(corr_pred_rate)/dt,'r');
xlabel('Time (s)');
ylabel('Firing rate (Hz)');
xlim(xr);

%%
fig_width = 4.6; 
rel_height = 1.5;

figufy(f1);
fname = [fig_dir 'examp_rpt_nonparam.pdf'];
exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f1);

fig_width = 3.27; 
rel_height = 0.8;

figufy(f2);
fname = [fig_dir 'examp_rpt_eyedist.pdf'];
exportfig(f2,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f2);

fig_width = 4.6; 
rel_height = 0.8;

figufy(f3);
fname = [fig_dir 'examp_rpt_psth.pdf'];
exportfig(f3,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
close(f3);

%%
    beta = glmfit(temp',yy_obs(uu),'poisson'); %beta(1) = [];

    sp = NMMcreate_stim_params([1 6]);
    spr = NMMcreate_reg_params('lambda_L2',1);
    
    mod = NMMinitialize_model(sp,1,{'lin'},spr,[],[],'exp');
    mod = NMMfit_filters(mod,yy_obs(uu),temp');
    
%%
n_use_trials = size(corr_pred_rate,1);

use_eye_range = [-0.35 0.35];
n_pts = 50;
n_knots = 10;
spline_order = 4;

eye_xx = linspace(use_eye_range(1),use_eye_range(2),n_pts);
knot_pts = [ones(1,spline_order-1)*use_eye_range(1) ...
    linspace(use_eye_range(1),use_eye_range(2),n_knots) ones(1,spline_order-1)*use_eye_range(2)];

Bspline_mat = bspline_basismatrix(spline_order,knot_pts,eye_xx);
pbins = [1 99];

smoothed_eye_rate_fun = nan(n_pts,length(rpt_taxis));
for ii = 1:length(rpt_taxis);
    ii
    xx = inf_eye_pos_t(:,ii);
%     yy = corr_pred_rate(:,ii);
    yy = full_psth(:,ii);
     knot_pts = prctile(xx(uu),linspace(pbins(1),pbins(end),n_knots));
    knot_pts = [prctile(xx(uu),pbins(1))*ones(1,spline_order-1) knot_pts prctile(xx(uu),pbins(2))*ones(1,spline_order-1)];
   
    uu = ~isnan(xx(:)) & ~isnan(yy(:)) & xx >= knot_pts(1) & xx <= knot_pts(end);
    temp = fnval(fnxtr(csape(knot_pts,eye(length(knot_pts)),'var')),xx(uu));
%     beta = yy(uu)'/temp;
    beta = glmfit(temp',yy(uu),'poisson'); %beta(1) = [];
%     pfun = fnxtr(csape(knot_pts,beta,xx(uu),'var'));
%     pfun = fnxtr(csape(knot_pts,beta(2:end),xx(uu),'var'));
%     smoothed_eye_rate_fun(:,ii) = fnval(pfun,eye_xx);
    pfun = fnxtr(csape(knot_pts,beta(2:end),xx(uu),'var'));
    smoothed_eye_rate_fun(:,ii) = exp(beta(1) + fnval(pfun,eye_xx));
end

%%

sm_win = 25;
% smoothed_eye_rate_fun = nan(size(corr_pred_rate));
smoothed_eye_robs_fun = nan(size(corr_pred_rate));
eye_rate_fun = nan(size(corr_pred_rate));
for ii = 1:length(rpt_taxis)
    ii
    xx = sorted_inf_eyepos(:,ii);
    yy = corr_pred_rate(inf_eyepos_ord(:,ii),ii);
    uu = find(~isnan(xx) & ~isnan(yy));
%     YY = smooth(xx(uu),yy(uu),sm_win,'rloess');
%     smoothed_eye_rate_fun(uu,ii) = YY;
    
    yy = full_psth(inf_eyepos_ord(:,ii),ii);
    YY = smooth(xx(uu),yy(uu),sm_win,'lowess');
    smoothed_eye_robs_fun(uu,ii) = YY;
%     smoothed_eye_robs_fun(uu,ii) = yy(uu);

    eye_rate_fun(:,ii) = corr_pred_rate(inf_eyepos_ord(:,ii),ii);
end

%%
xr = [0.25 1.25]; %1-2 2-24
yr = [-0.35 0.35];


h1 = figure();
pcolor(rpt_taxis,sorted_inf_eyepos,smoothed_eye_rate_fun);shading flat
xlim(xr);
ylim(yr);
caxis(cr);


