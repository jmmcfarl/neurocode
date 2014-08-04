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

cur_unit_inds = tr_set(ismember(all_mod_SUnum(tr_set),SU_numbers));
cur_runit_inds = find(ismember(SU_numbers,all_mod_SUnum(cur_unit_inds)));

if Expt_num == 277
    cur_unit_inds(2) = [];
    cur_runit_inds(2) = [];
end

fin_fix_corr = nan(NT,1);
fin_fix_corr(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
fin_fix_corr = fin_fix_corr*sp_dx;
fin_drift_corr = squeeze(drift_post_mean(end,:)*sp_dx);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > sac_shift
        fin_drift_corr(cur_inds(1:end-sac_shift+1)) = fin_drift_corr(cur_inds(sac_shift:end));
    end
end

fin_tot_corr = fin_fix_corr + fin_drift_corr;
fin_tot_corr = interp1(find(~isnan(fin_tot_corr)),fin_tot_corr(~isnan(fin_tot_corr)),1:NT);

cur_rpt_trials = find(all_trial_Se == rpt_seed);
all_used_rpt_inds = find(ismember(all_trialvec(used_inds),cur_rpt_trials));
all_used_inds = find(ismember(all_blockvec(used_inds),cur_used_blocks));

cur_n_rpts = length(cur_rpt_trials);
full_psth = nan(cur_n_rpts,length(rpt_taxis),length(cur_unit_inds));
y_cnt = 0; y_space = 1;
for ii = 1:cur_n_rpts
    cur_inds = find(all_trialvec(used_inds) == cur_rpt_trials(ii));
    cur_inds(size(full_psth,2)+1:end) = [];
    full_psth(ii,1:length(cur_inds),:) = all_binned_sua(used_inds(cur_inds),cur_runit_inds);    
end

%NOW CONSTRUCT CORRECTED STIMULUS
fix_post_cor = nan(NT,1);
fix_post_cor(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
fix_post_cor = interp1(find(~isnan(fix_ids)),fix_post_cor(~isnan(fix_ids)),1:NT);
fix_post_cor(isnan(fix_post_cor)) = 0;
%back project drift (within-fixation) by sac_shift
drift_post_cor = squeeze(drift_post_mean(end,:));
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > sac_shift
        drift_post_cor(cur_inds(1:end-sac_shift+1)) = drift_post_cor(cur_inds(sac_shift:end));
    end
end
drift_post_cor(isnan(drift_post_cor)) = 0;

all_post_cor = round((fix_post_cor+drift_post_cor)) + max_Tshift + 1;

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

for jj = 1:length(cur_runit_inds)
    %CORRECTED AND UNCORRECTED MODEL FITS
    cur_cormod = dit_mods_spkNL{end}(cur_unit_inds(jj));
    cur_Xtargs = [cur_cormod.mods(:).Xtarget];
    cur_cormod.mods(cur_Xtargs > 1) = [];
    cur_uinds = find(~isnan(all_binned_sua(used_inds(all_used_rpt_inds),cur_runit_inds(jj))));
    cur_cormod = NMMfit_logexp_spkNL(cur_cormod,all_binned_sua(used_inds(all_used_rpt_inds(cur_uinds)),cur_runit_inds(jj)),cor_X(cur_uinds,:),[],[],[1 2]);
    
    cur_uncormod = it_mods_spkNL{1}(cur_unit_inds(jj));
    cur_uncormod.mods(cur_Xtargs > 1) = [];
    cur_uncormod = NMMfit_logexp_spkNL(cur_uncormod,all_binned_sua(used_inds(all_used_rpt_inds(cur_uinds)),cur_runit_inds(jj)),uncor_X(cur_uinds,:),[],[],[1 2]);
    
    %COMPUTE TRIAL-BY-TRIAL PREDICTED RATES
    corr_pred_rate = nan(length(cur_rpt_trials),length(rpt_taxis));
    uncorr_pred_rate = nan(length(cur_rpt_trials),length(rpt_taxis));
    row_ind = nan(length(cur_rpt_trials),1);
    for ii = 1:length(cur_rpt_trials)
        row_ind(ii) = ii;
        test_rpt_inds = find(all_trialvec(used_inds(all_used_rpt_inds)) == cur_rpt_trials(ii));
        test_rpt_inds(length(rpt_taxis)+1:end) = [];
        [~, ~, corr_pred_rate(ii,1:length(test_rpt_inds))] = NMMmodel_eval(cur_cormod,[],cor_X(test_rpt_inds,:));
        [~, ~, uncorr_pred_rate(ii,1:length(test_rpt_inds))] = NMMmodel_eval(cur_uncormod,[],uncor_X(test_rpt_inds,:));
    end
    all_corr_pred_rate(:,:,jj) = corr_pred_rate;
    all_uncorr_pred_rate(:,:,jj) = uncorr_pred_rate;
end

inf_eye_pos_t = nan(length(cur_rpt_trials),length(rpt_taxis));
row_ind = nan(length(cur_rpt_trials),1);
for ii = 1:length(cur_rpt_trials)
    row_ind(ii) = ii;
    test_rpt_inds = find(all_trialvec(used_inds(all_used_rpt_inds)) == cur_rpt_trials(ii));
    test_rpt_inds(length(rpt_taxis)+1:end) = [];
    
    inf_eye_pos_t(ii,1:length(test_rpt_inds)) = fin_tot_corr(all_used_rpt_inds(test_rpt_inds));
end

all_corr_pred_rate(isnan(full_psth)) = nan;
all_uncorr_pred_rate(isnan(full_psth)) = nan;
%%
all_corr_pred_rate_reshaped = reshape(all_corr_pred_rate,length(cur_rpt_trials)*length(rpt_taxis),size(all_corr_pred_rate,3));
all_uncorr_pred_rate_reshaped = reshape(all_uncorr_pred_rate,length(cur_rpt_trials)*length(rpt_taxis),size(all_corr_pred_rate,3));
all_psth_reshaped = reshape(full_psth,length(cur_rpt_trials)*length(rpt_taxis),size(all_corr_pred_rate,3));

%%
clc
pair_id = [4 5];
n_perms = 25;
used_trials = find(~any(isnan(full_psth(:,:,pair_id(1))),2) & ~any(isnan(full_psth(:,:,pair_id(2))),2));
tlen = length(used_trials)*length(rpt_taxis);
shuff_corr = nan(n_perms,1);
shuff_pcorr = nan(n_perms,1);
cur_strains = nan(tlen,2);
for ii = 1:n_perms
   rperm = randperm(length(used_trials));
   cur_strains(:,1) = reshape(full_psth(used_trials,:,pair_id(1)),tlen,1);
   cur_strains(:,2) = reshape(full_psth(used_trials(rperm),:,pair_id(2)),tlen,1);
   shuff_corr(ii) = corr(cur_strains(:,1),cur_strains(:,2));

   cur_ptrains(:,1) = reshape(all_corr_pred_rate(used_trials,:,pair_id(1)),tlen,1);
   cur_ptrains(:,2) = reshape(all_corr_pred_rate(used_trials(rperm),:,pair_id(2)),tlen,1);
   shuff_pcorr(ii) = corr(cur_ptrains(:,1),cur_ptrains(:,2));
end
sp=signrank(shuff_corr);
fprintf('Shuffled corr: %.4f  p: %.4f\n',mean(shuff_corr),sp);
sp=signrank(shuff_pcorr);
fprintf('Shuffled pcorr: %.4f  p: %.4f\n',mean(shuff_pcorr),sp);

[a,b] = nancorr(all_psth_reshaped(:,pair_id));
fprintf('Measured corr: %.4f  p: %.4f\n',a(2,1),b(2,1));

[a,b] = nancorr(all_corr_pred_rate_reshaped(:,pair_id));
fprintf('Model pred corr: %.4f  p: %.4f\n',a(2,1),b(2,1));

[a,b] = nancorr(all_uncorr_pred_rate_reshaped(:,pair_id));
fprintf('Uncorr model predicted corr: %.4f  p: %.4f\n',a(2,1),b(2,1));


%%
%CORRECTED AND UNCORRECTED MODEL FITS
unit_1 = 1;
unit_2 = 2;

cur_cormod = dit_mods{end}(cur_unit_inds(unit_1));
cur_Xtargs = [cur_cormod.mods(:).Xtarget];
cur_cormod.mods(cur_Xtargs > 1) = [];
cur_uncormod = it_mods{1}(cur_unit_inds(unit_1));
cur_uncormod.mods(cur_Xtargs > 1) = [];

cur_uinds = find(~isnan(all_binned_sua(used_inds(all_used_rpt_inds),cur_runit_inds(unit_1))) & ...
    ~isnan(all_binned_sua(used_inds(all_used_rpt_inds),cur_runit_inds(unit_2))));
Robs = all_binned_sua(used_inds(all_used_rpt_inds(cur_uinds)),cur_runit_inds(unit_1));
Rpred = all_binned_sua(used_inds(all_used_rpt_inds(cur_uinds)),cur_runit_inds(unit_2));

cur_cormod = NMMfit_logexp_spkNL(cur_cormod,Robs,cor_X(cur_uinds,:),[],[],[1 2]);
cur_uncormod = NMMfit_logexp_spkNL(cur_uncormod,Robs,uncor_X(cur_uinds,:),[],[],[1 2]);

[LL, penLL, pred_rate, G, gint] = NMMmodel_eval(cur_cormod, Robs, cor_X(cur_uinds,:));
[LLu, penLL, pred_rate, unG, gint] = NMMmodel_eval(cur_uncormod, Robs, uncor_X(cur_uinds,:));
fullLL = nansum(Robs.*log(Robs) - Robs);

% [B,dev,stats] = glmfit([G Rpred],Robs,'poisson');
% [Bu,devu,statsu] = glmfit([unG Rpred],Robs,'poisson');
[B,dev,stats] = glmfit(G(:),Robs(:),'poisson');
[Bu,devu,statsu] = glmfit(unG(:),Robs(:),'poisson');

cur_stim_params = NMMcreate_stim_params([1 1],dt);
testmod = NMMinitialize_model(cur_stim_params,1,{'lin'});
testmod = NMMfit_filters(testmod,Robs,G);
[LL] = NMMmodel_eval(testmod, Robs, G);


[B2,dev2,stats2] = glmfit([Rpred],Robs,'poisson');

avgrate = nanmean(Robs);
nullpred = avgrate*ones(size(Robs));
nulldev = 2*(nansum(Robs.*log(Robs) - Robs) - nansum(Robs.*log(nullpred) - nullpred));

full_R2 = 1-dev/nulldev;
uR2 = 1-devu/nulldev;
spike_R2 = 1-dev2/nulldev;

%%
%%
rpt_taxis = (1:round(trial_dur)/dt)*dt-dt/2;
rpt_taxis(rpt_taxis < beg_buffer) = [];
rpt_taxis(trial_dur - rpt_taxis < end_buffer) = [];
all_trial_dur = all_trial_end_times-all_trial_start_times;

% cur_unit_inds = tr_set(ismember(all_mod_SUnum(tr_set),SU_numbers));
% cur_runit_inds = find(ismember(SU_numbers,all_mod_SUnum(cur_unit_inds)));
cur_unit_inds = tr_set;

fin_fix_corr = nan(NT,1);
fin_fix_corr(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);
fin_fix_corr = fin_fix_corr*sp_dx;
fin_drift_corr = squeeze(drift_post_mean(end,:)*sp_dx);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > sac_shift
        fin_drift_corr(cur_inds(1:end-sac_shift+1)) = fin_drift_corr(cur_inds(sac_shift:end));
    end
end

fin_tot_corr = fin_fix_corr + fin_drift_corr;
fin_tot_corr = interp1(find(~isnan(fin_tot_corr)),fin_tot_corr(~isnan(fin_tot_corr)),1:NT);

cur_rpt_trials = find(all_trial_Se == rpt_seed);
all_used_rpt_inds = find(ismember(all_trialvec(used_inds),cur_rpt_trials));
all_used_inds = find(ismember(all_blockvec(used_inds),cur_used_blocks));

all_binned_ua = [all_binned_mua all_binned_sua];

cur_n_rpts = length(cur_rpt_trials);
full_psth = nan(cur_n_rpts,length(rpt_taxis),length(cur_unit_inds));
y_cnt = 0; y_space = 1;
for ii = 1:cur_n_rpts
    cur_inds = find(all_trialvec(used_inds) == cur_rpt_trials(ii));
    cur_inds(size(full_psth,2)+1:end) = [];
    full_psth(ii,1:length(cur_inds),:) = all_binned_ua(used_inds(cur_inds),cur_unit_inds);    
end

%NOW CONSTRUCT CORRECTED STIMULUS
fix_post_cor = nan(NT,1);
fix_post_cor(~isnan(fix_ids)) = it_fix_post_mean(end,fix_ids(~isnan(fix_ids)));
fix_post_cor = interp1(find(~isnan(fix_ids)),fix_post_cor(~isnan(fix_ids)),1:NT);
fix_post_cor(isnan(fix_post_cor)) = 0;
%back project drift (within-fixation) by sac_shift
drift_post_cor = squeeze(drift_post_mean(end,:));
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
    if length(cur_inds) > sac_shift
        drift_post_cor(cur_inds(1:end-sac_shift+1)) = drift_post_cor(cur_inds(sac_shift:end));
    end
end
drift_post_cor(isnan(drift_post_cor)) = 0;

all_post_cor = round((fix_post_cor+drift_post_cor)) + max_Tshift + 1;

%RECOMPUTE XMAT
all_shift_stimmat_up = all_stimmat_up;
for i=1:NT
    all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_post_cor(i)};
end
drift_Xmat = create_time_embedding(all_shift_stimmat_up,stim_params_us);
drift_Xmat = drift_Xmat(used_inds,use_kInds_up);

clear cor_X uncor_X
cor_X{1} = drift_Xmat(all_used_rpt_inds,:);
uncor_X{1} = all_Xmat_us(used_inds(all_used_rpt_inds),use_kInds_up);
cor_X{2} = Xblock;
cor_X{3} = Xsac;
cor_X{4} = Xmsac;
uncor_X{2} = Xblock;
uncor_X{3} = Xsac;
uncor_X{4} = Xmsac;

all_corr_pred_rate = nan(length(cur_rpt_trials),length(rpt_taxis),length(cur_unit_inds));
all_uncorr_pred_rate = nan(length(cur_rpt_trials),length(rpt_taxis),length(cur_unit_inds));
for jj = 1:length(cur_unit_inds)
    jj
    %CORRECTED AND UNCORRECTED MODEL FITS
    cur_cormod = dit_mods_spkNL{end}(cur_unit_inds(jj));
    cur_Xtargs = [cur_cormod.mods(:).Xtarget];
%     cur_cormod.mods(cur_Xtargs > 1) = [];
    cur_cormod.mods(cur_Xtargs == 2) = [];
    cur_uinds = find(~isnan(all_binned_ua(used_inds(all_used_rpt_inds),cur_unit_inds(jj))));
    cur_cormod = NMMfit_logexp_spkNL(cur_cormod,all_binned_ua(used_inds(all_used_rpt_inds(cur_uinds)),cur_unit_inds(jj)),get_Xcell_tInds(cor_X,cur_uinds),[],[],[1 2]);
    
    cur_uncormod = it_mods_spkNL{1}(cur_unit_inds(jj));
%     cur_uncormod.mods(cur_Xtargs > 1) = [];
    cur_uncormod.mods(cur_Xtargs == 2) = [];
    cur_uncormod = NMMfit_logexp_spkNL(cur_uncormod,all_binned_ua(used_inds(all_used_rpt_inds(cur_uinds)),cur_unit_inds(jj)),get_Xcell_tInds(uncor_X,cur_uinds),[],[],[1 2]);
    
    %COMPUTE TRIAL-BY-TRIAL PREDICTED RATES
    corr_pred_rate = nan(length(cur_rpt_trials),length(rpt_taxis));
    uncorr_pred_rate = nan(length(cur_rpt_trials),length(rpt_taxis));
    row_ind = nan(length(cur_rpt_trials),1);
    for ii = 1:length(cur_rpt_trials)
        row_ind(ii) = ii;
        test_rpt_inds = find(all_trialvec(used_inds(all_used_rpt_inds)) == cur_rpt_trials(ii));
        test_rpt_inds(length(rpt_taxis)+1:end) = [];
        [~, ~, corr_pred_rate(ii,1:length(test_rpt_inds))] = NMMmodel_eval(cur_cormod,[],get_Xcell_tInds(cor_X,test_rpt_inds));
        [~, ~, uncorr_pred_rate(ii,1:length(test_rpt_inds))] = NMMmodel_eval(cur_uncormod,[],get_Xcell_tInds(uncor_X,test_rpt_inds));
    end
    all_corr_pred_rate(:,:,jj) = corr_pred_rate;
    all_uncorr_pred_rate(:,:,jj) = uncorr_pred_rate;
end

inf_eye_pos_t = nan(length(cur_rpt_trials),length(rpt_taxis));
row_ind = nan(length(cur_rpt_trials),1);
for ii = 1:length(cur_rpt_trials)
    row_ind(ii) = ii;
    test_rpt_inds = find(all_trialvec(used_inds(all_used_rpt_inds)) == cur_rpt_trials(ii));
    test_rpt_inds(length(rpt_taxis)+1:end) = [];
    
    inf_eye_pos_t(ii,1:length(test_rpt_inds)) = fin_tot_corr(all_used_rpt_inds(test_rpt_inds));
end

all_corr_pred_rate(isnan(full_psth)) = nan;
all_uncorr_pred_rate(isnan(full_psth)) = nan;
%%
all_corr_pred_rate_reshaped = reshape(all_corr_pred_rate,length(cur_rpt_trials)*length(rpt_taxis),size(all_corr_pred_rate,3));
all_uncorr_pred_rate_reshaped = reshape(all_uncorr_pred_rate,length(cur_rpt_trials)*length(rpt_taxis),size(all_corr_pred_rate,3));
all_psth_reshaped = reshape(full_psth,length(cur_rpt_trials)*length(rpt_taxis),size(all_corr_pred_rate,3));

%%
n_shuffs = 10;
shuff_corrs = zeros(length(cur_unit_inds));
shuff_pcorrs = zeros(length(cur_unit_inds));
for ii = 1:n_shuffs
    ii
    cur_cpred_rate = all_corr_pred_rate;
    cur_obs_rate = full_psth;
    for jj = 1:length(cur_unit_inds)
        cur_tord = randperm(size(all_corr_pred_rate,1));
        cur_cpred_rate(:,:,jj) = cur_cpred_rate(cur_tord,:,jj);
        cur_obs_rate(:,:,jj) = cur_obs_rate(cur_tord,:,jj);
    end
   shuff_corr_pred_rate_reshaped = reshape(cur_cpred_rate,length(cur_rpt_trials)*length(rpt_taxis),size(all_corr_pred_rate,3));
   cur_corr = nancorr(shuff_corr_pred_rate_reshaped);
   shuff_corrs = shuff_corrs + cur_corr;
  
   shuff_psth_reshaped = reshape(cur_obs_rate,length(cur_rpt_trials)*length(rpt_taxis),size(all_corr_pred_rate,3));
   cur_corr = nancorr(shuff_psth_reshaped);
   shuff_pcorrs = shuff_pcorrs + cur_corr;
end
shuff_corrs = shuff_corrs/n_shuffs;
shuff_pcorrs = shuff_pcorrs/n_shuffs;

corrected_rate_corrs = nancorr(all_corr_pred_rate_reshaped) - shuff_corrs;
corrected_cnt_corrs = nancorr(all_psth_reshaped)-shuff_pcorrs;
