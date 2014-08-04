clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

%run on 86 [0, 90]; 
Expt_num = 86;    
Expt_name = sprintf('G%.3d',Expt_num);
data_dir = ['~/Data/bruce/' Expt_name];
cd(data_dir);
fig_dir = '/home/james/Analysis/bruce/ET_final/';

load(sprintf('jbe%sExpts.mat',Expt_name)); %load in Expts struct
load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final/'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

bar_ori = 0;

true_data_name = 'monoc_eyecorr_hbar_highres3';
true_moddata_name = 'monoc_eyecorr_hbar_mods_hres2.mat';

sim_gain = 0.5;
if sim_gain >= 1
    gname = num2str(sim_gain);
elseif sim_gain == 0.5
    gname = 'p5';
elseif sim_gain == 0.25
    gname = 'p25';
else
    error('unsupported');
end
sim_data_name = ['monoc_eyecorr_hbar_simdata_g' gname]';

sim_spatial_usfac = 8;
spatial_usfac = 4;
old_usfac = 2;

use_measured_pos = 0;
use_smooth_eyepos = false;
use_sac_kerns = 1;
use_coils = [0 0]; %[L R]

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
xv_frac = 0;

flen = 12;
use_nPix = 16;

n_fix_inf_it = 1; %4
n_drift_inf_it = 1; %2

fix_prior_sigma = 0.15;
fix_noise_sigma = 0.1;
drift_noise_sigma = 0.003;
drift_prior_sigma = 0.004; %.004 may be best here
drift_jump_sigma = 0.1; %0.05 start
drift_dsf = 2;

min_trial_dur = 0.75;


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
sim_sp_dx = 0.0565/sim_spatial_usfac;

max_shift_deg = 1*sim_gain;
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

n_blocks = length(cur_block_set);

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
all_trial_wi = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);

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
            cur_stim_mat = double(left_stim_mats{use_trials(tt)}(1:use_frames,cur_use_pix));
            %             bar_Xmat = create_time_embedding(cur_stim_mat,stim_params);
            
            if ~isempty(all_stim_times)
                if any(cur_stim_times < all_stim_times(end))
                    fprintf('Warn trial %d\n',tt);
                end
            end
            all_stim_times = [all_stim_times; cur_stim_times];
            all_t_axis = [all_t_axis; cur_t_axis];
            %             all_Xmat = [all_Xmat; bar_Xmat];
            all_stim_mat = [all_stim_mat; cur_stim_mat];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(size(cur_t_axis))*ee];
            all_trialvec = [all_trialvec; ones(size(cur_t_axis))*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
        end
    end
    trial_cnt = trial_cnt + n_trials;
end

%%
if sim_gain > 1
    full_nPix = full_nPix + sim_gain*8;
end
stim_params = NIMcreate_stim_params([flen full_nPix],dt);
if sim_gain > 1
    extra_pix = full_nPix - 36;

    extra_stim = zeros(length(all_t_axis),extra_pix);
    rand_mat = rand(length(all_t_axis),extra_pix);
    extra_stim(rand_mat <= 0.06) = -1;
    extra_stim(rand_mat >= 0.94) = 1;
else
    extra_stim = [];
end
all_stim_mat = cat(2,all_stim_mat,extra_stim);

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
    if ~isempty(used_clust_set)
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
end

%for only-MU probes
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
    
    %remove spikes from isolated SUs from the MU
    for ss = 1:length(unique_su_nums)
       cur_mua_inds(ismember(all_spk_inds{cc}(cur_mua_inds),all_su_spk_inds{unique_su_nums(ss)})) = []; 
    end
    
    [cur_spkhist,cur_bin_inds] = histc(all_spk_times{cc}(cur_mua_inds),all_t_bin_edges);
    cur_spkhist(all_bin_edge_pts) = [];
    all_binned_mua(:,cc) = cur_spkhist;
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

%% select submatrix with central pixels
[Xinds,Tinds] = meshgrid(1:full_nPix,1:flen);
buffer_pix = floor((full_nPix - use_nPix)/2);
cur_use_pix = (1:use_nPix) + buffer_pix;
use_kInds = find(ismember(Xinds(:),cur_use_pix));

xax =1/spatial_usfac:1/spatial_usfac:full_nPix;
[Xinds_up,Tinds] = meshgrid(xax,1:flen);
if spatial_usfac > 1
%     use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1)-1/spatial_usfac & Xinds_up(:) <= cur_use_pix(end));
    use_kInds_up = find(Xinds_up(:) >= cur_use_pix(1)-2/spatial_usfac & Xinds_up(:) <= cur_use_pix(end) + 1/spatial_usfac);
else
    use_kInds_up = use_kInds;
end
use_kInds_back = find(ismember(Xinds_up(use_kInds_up),cur_use_pix));

use_xax_up = find(xax >= cur_use_pix(1) - 2/spatial_usfac & xax <= cur_use_pix(end) + 1/spatial_usfac);

full_nPix_us = spatial_usfac*full_nPix;



%%
sim_nPix_us = full_nPix*sim_spatial_usfac;

all_stimmat_up = zeros(size(all_stim_mat,1),sim_nPix_us);
for ii = 1:size(all_stim_mat,2)
    for jj = 1:sim_spatial_usfac
    all_stimmat_up(:,sim_spatial_usfac*(ii-1)+jj) = all_stim_mat(:,ii);
    end
end

sim_nPix_us = sim_spatial_usfac*use_nPix;
xax =1/sim_spatial_usfac:1/sim_spatial_usfac:full_nPix;
[Xinds_sim,Tinds] = meshgrid(1/sim_spatial_usfac:1/sim_spatial_usfac:full_nPix,1:flen);
cnt = 1;
use_xax_sim = [];
while length(use_xax_sim) < sim_nPix_us
% use_kInds_sim = find(Xinds_sim(:) >= cur_use_pix(1)-4/sim_spatial_usfac & Xinds_sim(:) <= cur_use_pix(end) + 3/sim_spatial_usfac);
use_xax_sim = find(xax >= cur_use_pix(1) - cnt/sim_spatial_usfac & xax <= cur_use_pix(end) + cnt/sim_spatial_usfac);
if length(use_xax_sim) < sim_nPix_us
use_xax_sim = find(xax >= cur_use_pix(1) - (cnt+1)/sim_spatial_usfac & xax <= cur_use_pix(end) + cnt/sim_spatial_usfac);
end
cnt = cnt + 1;
end

%%
fprintf('Incorporating inferred eye-positions\n');
cd(anal_dir)

load(true_data_name,'dit_mods','et_tr_set','best_fix_cor','drift_post_*');
true_mods = dit_mods{1};
tr_set = et_tr_set;
n_tr_chs = length(tr_set);

load(true_moddata_name,'all_mod_SU*');

fin_fix_corr = nan(NT,1);
fin_fix_std = nan(NT,1);
fin_fix_corr(~isnan(fix_ids)) = best_fix_cor(fix_ids(~isnan(fix_ids)));
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids)),1:NT);

fin_fix_corr = fin_fix_corr*sp_dx;

fin_drift_corr = drift_post_mean(end,:)*sp_dx;

n_trials = length(trial_start_inds);
for ii = 1:n_trials
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    fin_drift_corr(cur_inds(1:end-sac_shift)) = fin_drift_corr(cur_inds(sac_shift+1:end));
end
%within-saccade estimates use linear interpolation
fin_drift_corr = interp1(find(~isnan(fix_ids)),fin_drift_corr(~isnan(fix_ids)),1:NT);

inferred_eyepos = (fin_fix_corr + fin_drift_corr);


%%
if use_smooth_eyepos
    
    sim_eyepos = corrected_eye_vals_interp(:,2);
    sim_eyepos = sim_eyepos*sim_gain;
    sim_eyepos(isnan(sim_eyepos)) = 0;
    
    %smooth out fast transients in eye signal
    eye_smooth_sig = round(0.025/dt);
    interp_inds = [];
    for ii = 1:n_fixs
        cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
        if length(cur_inds) > eye_smooth_sig*5;
            sim_eyepos(used_inds(cur_inds)) = jmm_smooth_1d_cor(sim_eyepos(used_inds(cur_inds)),eye_smooth_sig,2);
        end
        interp_inds = [interp_inds; cur_inds'];
    end
    interp_inds = unique(interp_inds);
    sim_eyepos(used_inds) = interp1(used_inds(interp_inds),sim_eyepos(used_inds(interp_inds)),used_inds);
    
    for ii = 1:length(cur_block_set)
        cur_inds = find(all_blockvec == ii);
        sim_eyepos(cur_inds,:) = bsxfun(@minus,sim_eyepos(cur_inds,:),median(sim_eyepos(cur_inds,:)));
    end
    
    if use_med_sub
        %subtract off within-trial median
        for ii = 1:length(trial_start_inds)
            cur_inds = trial_start_inds(ii):trial_end_inds(ii);
            sim_eyepos(used_inds(cur_inds),:) = bsxfun(@minus,sim_eyepos(used_inds(cur_inds),:),median(sim_eyepos(used_inds(cur_inds),:)));
        end
    end
    
    %scale to match SD of inferred
    inf_std = std(inferred_eyepos);
    meas_std = std(smooth_eyepos(used_inds));
    sim_eyepos = sim_eyepos*inf_std/meas_std;

    %maximum initial corrections
    max_sim_pos = max_shift*sp_dx;
    sim_eyepos(sim_eyepos > max_sim_pos) = max_sim_pos;
    sim_eyepos(sim_eyepos < - max_sim_pos) = -max_sim_pos;
    sim_eyepos_us_rnd = round(sim_eyepos/sp_dx);
    
    %%
else
    sim_eyepos = zeros(length(all_t_axis),1);
    sim_eyepos(used_inds) = inferred_eyepos;
    sim_eyepos = sim_eyepos*sim_gain;
    
    %maximum initial corrections
    sim_eyepos(sim_eyepos > max_shift_deg) = max_shift_deg;
    sim_eyepos(sim_eyepos < - max_shift_deg) = -max_shift_deg;
    sim_eyepos_us_rnd = round(sim_eyepos/(sim_sp_dx));

end


%%
sim_eyepos_us_rnd(sim_eyepos_us_rnd > size(all_stimmat_up,2)) = size(all_stimmat_up,2);
sim_eyepos_us_rnd(sim_eyepos_us_rnd < -size(all_stimmat_up,2)) = -size(all_stimmat_up,2);
all_stimmat_true = all_stimmat_up;
for ii = 1:length(sim_eyepos)
    all_stimmat_true(ii,:) = shift_matrix_Nd(all_stimmat_up(ii,:),-sim_eyepos_us_rnd(ii),2);
end
all_stimmat_true = all_stimmat_true(:,use_xax_sim);

use_nPix_sim = use_nPix*sim_spatial_usfac;
sim_stim_params = NMMcreate_stim_params([flen use_nPix_sim],dt);
all_Xmat_true = create_time_embedding(all_stimmat_true,sim_stim_params);
all_Xmat_true = all_Xmat_true(used_inds,:);

%%
if sim_spatial_usfac > spatial_usfac
    new_ax = 1:0.5:(use_nPix_us+0.5);
    for ii = 1:n_tr_chs
        true_mods(tr_set(ii)).stim_params(1).stim_dims(2) = use_nPix*sim_spatial_usfac;
        cur_stim_filts = [true_mods(tr_set(ii)).mods(1:3).filtK];
        cur_stim_filts = reshape(cur_stim_filts,[flen use_nPix_us 3]);
        cur_stim_filts = permute(cur_stim_filts,[2 1 3]);
        cur_stim_filts = interp1((1:use_nPix_us)',cur_stim_filts,new_ax','spline');
        cur_stim_filts = permute(cur_stim_filts,[2 1 3]);
        cur_stim_filts = cur_stim_filts*spatial_usfac/sim_spatial_usfac;
        cur_stim_filts = reshape(cur_stim_filts,[flen*use_nPix*sim_spatial_usfac 3]);
        true_mods(tr_set(ii)).mods(1).filtK = cur_stim_filts(:,1);
        true_mods(tr_set(ii)).mods(2).filtK = cur_stim_filts(:,2);
        true_mods(tr_set(ii)).mods(3).filtK = cur_stim_filts(:,3);
    end
end
%% SIMULATE SPIKING RESPONSES
fprintf('Creating simulated spike trains\n');

tr_X{1} = all_Xmat_true;
tr_X{2} = Xblock(used_inds,:);
tr_X{3} = Xsac;
tr_X{4} = Xmsac;

% make Robs_mat
sim_prate_mat = nan(length(used_inds),n_tr_chs);
for ss = 1:n_tr_chs
    [~, ~, pred_rate] = NMMmodel_eval(true_mods(tr_set(ss)),[], tr_X);
    if tr_set(ss) <= n_probes
        cur_use = find(~isnan(all_binned_mua(used_inds,tr_set(ss))));
    else
        su_ind = find(SU_numbers == all_mod_SUnum(tr_set(ss)));
        cur_use = find(~isnan(all_binned_sua(used_inds,su_ind)));
    end
    sim_prate_mat(cur_use,ss) = pred_rate(cur_use);
end
Robs_mat = poissrnd(sim_prate_mat);
Robs_xv = poissrnd(sim_prate_mat);

%%
cd(anal_dir)
sim_tr_set = tr_set;
save(sim_data_name,'sim_eyepos','sim_sp_dx','sim_prate_mat','Robs_mat','Robs_xv','sim_tr_set','true_mods','all_mod_SU*','extra_stim')
