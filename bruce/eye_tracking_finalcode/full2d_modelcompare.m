clear all
% close all
addpath('~/James_scripts/bruce/eye_tracking/');

Expt_num = 103;
Expt_name = sprintf('G%d',Expt_num);
% data_dir = ['/media/NTlab_data1/Data/bruce/' Expt_name];
data_dir = ['/home/james/Data/bruce/' Expt_name];
if ~exist(data_dir,'dir');
    system(['mkdir ' data_dir]);
end
cd(data_dir);

load(sprintf('%sExpts.mat',Expt_name)); %load in Expts struct

anal_dir = ['~/Analysis/bruce/' Expt_name];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end
fig_dir = '/home/james/Analysis/bruce/ET_final/';

mod_data_name = 'Twod_eyecorr_mods_flen1_v2';
anal_name = 'Twod_eyecorr_flen1_v3';

recompute_init_mods = 0;
use_measured_pos = 2;
use_sac_kerns = 1;
use_LOOXV = 0;

if use_measured_pos == 1
    mod_data_name = [mod_data_name '_Cinit'];
    anal_name = [anal_name '_Cinit'];
elseif use_measured_pos == 2
    mod_data_name = [mod_data_name '_Rinit'];
    anal_name = [anal_name '_Rinit'];
end

use_coils = [0 0]; %[L R]

if any(use_coils > 0)
   anal_name = [anal_name '_Cprior'];
end


%%
xv_frac = 0;
rpt_seed = 1450;

dsfrac = 2;
flen = 5;
if flen > 1
    fixed_delay = 0;
else
    fixed_delay = 2;
end

% n_fix_inf_it = 4; %3
% n_drift_inf_it = 1; %3

% fix_noise_sigma = 0.1;
% fix_prior_sigma = 0.2; 
% drift_noise_sigma = 0.01;
% drift_prior_sigma = 0.01;
% drift_jump_sigma = 0.075; 
% drift_dsf = 2;

min_trial_dur = 1.75;

%%
eps = -1e3;

dt = 2*100e-4;

beg_buffer = 0.15;
end_buffer = 0.05;
trial_dur = 4;

n_probes = 96;

n_use_blocks = Inf;

sac_backlag = round(0.05/dt);
sac_forlag = round(0.3/dt);
sac_bincents = -sac_backlag:sac_forlag;
n_sac_bins = length(sac_bincents);

Pix2Deg = 0.018837;
sp_dx = Pix2Deg*dsfrac;

max_shift_deg = 0.7;
max_shift = round(max_shift_deg/sp_dx);
dshift = 1;

max_Dshift_deg = 0.3;
max_Dshift = round(max_Dshift_deg/sp_dx);

%% load overall su data
% LOAD REFCLUSTERS
% cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering2'];
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
temp = cellfun(@(x) x.Header.expname,Expts,'uniformoutput',false);

cur_block_set = find(strcmp(temp,'checker.imiRC') | strcmp(temp,'image.imiRC'));
cur_block_set(cur_block_set == 6) = []; %model fits are bad for this block for some reason
n_blocks = length(cur_block_set);

%% COMPUTE TRIAL DATA
cd(data_dir);

fprintf('Computing prep data\n');
trial_cnt = 0;

all_t_axis = [];
all_t_bin_edges = [];
all_tsince_start = [];
all_blockvec = [];
all_trialvec = [];
all_trialSe = [];
all_Seeds = [];
all_trial_blocknums = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_bin_edge_pts = [];
all_spk_times = cell(n_probes,1);
all_spk_inds = cell(n_probes,1);
all_clust_ids = cell(n_probes,1);

cur_spkind_offset = 0;
cur_toffset = 0;
for ee = 1:n_blocks;
    fprintf('Block %d of %d\n',ee,n_blocks);
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
    trial_Se = [Expts{cur_block}.Trials(:).se];
    
    n_rpt_frames = nan(length(trial_start_times),1);
    if isfield(Expts{cur_block}.Trials,'rptframes')
        for tt = 1:length(trial_start_times)
            n_rpt_frames(tt) = length(Expts{cur_block}.Trials(tt).rptframes);
        end
    end
    
    use_trials = find(trial_durs >= min_trial_dur);
    
    dset = find(n_rpt_frames(use_trials) > 0);
    fprintf('Removing %d of %d trials with rpt frames\n',length(dset),length(use_trials));
    use_trials(dset) = [];
    
    all_trial_start_times = cat(1,all_trial_start_times,trial_start_times(use_trials)');
    all_trial_end_times = cat(1,all_trial_end_times,trial_end_times(use_trials)');
    all_trial_blocknums = cat(1,all_trial_blocknums,ones(length(use_trials),1)*ee);
    all_trialSe = cat(1,all_trialSe,trial_Se(use_trials)');
    n_trials = length(use_trials);
    for tt = 1:n_trials
        cur_Seeds = downsample(Expts{cur_block}.Trials(use_trials(tt)).Seedseq,2);
        if ~isempty(cur_Seeds)
            cur_stim_times = Expts{cur_block}.Trials(use_trials(tt)).Start/1e4;
            n_frames = length(cur_stim_times);
            cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt];
            if n_frames > 0
                if length(cur_stim_times) == 1
                    cur_stim_times = (cur_stim_times:dt*Fr:(cur_stim_times + (n_frames-1)*dt*Fr))';
                    cur_stim_times(cur_stim_times > trial_end_times(use_trials(tt))) = [];
                    cur_t_edges = [cur_stim_times; cur_stim_times(end) + dt*Fr];
                end
            end
            cur_t_axis = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
            
            if isempty(cur_Seeds)
                cur_Seeds = nan(n_frames,1);
            end
            cur_Seeds(n_frames+1:end) = [];
            
            cur_Seeds((end-fixed_delay+1):end) = [];
            cur_Seeds = [nan(fixed_delay,1); cur_Seeds];
            
            cur_tsince_start = cur_t_axis - trial_start_times(use_trials(tt));
            
            all_Seeds = [all_Seeds; cur_Seeds(:)];
            
            all_t_axis = [all_t_axis; cur_t_axis];
            all_t_bin_edges = [all_t_bin_edges; cur_t_edges];
            all_tsince_start = [all_tsince_start; cur_tsince_start];
            all_blockvec = [all_blockvec; ones(length(cur_t_axis),1)*ee];
            all_trialvec = [all_trialvec; ones(length(cur_t_axis),1)*(tt + trial_cnt)];
            all_bin_edge_pts = [all_bin_edge_pts; length(all_t_bin_edges)];
            
            if length(all_t_axis) ~= length(all_Seeds)
                error('t')
            end
        end
    end
    trial_cnt = trial_cnt + n_trials;
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
used_inds = used_inds(~isnan(all_Seeds(used_inds)));
NT = length(used_inds);

%% CREATE EVENT PREDICTORS FOR REAL AND SIM SACCADES (POOLING MICRO AND MACRO SACS)
Xblock = zeros(length(all_t_axis),n_blocks);
for i = 1:n_blocks
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end

%% PROCESS EYE TRACKING DATA
trial_toffset = zeros(length(cur_block_set),1);
[all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

expt_bar_ori = zeros(size(cur_block_set));
%compute corrected eye data in bar-oriented frame
lin_correction = false;
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,expt_bar_ori,used_inds,lin_correction);

[saccades,et_params] = detect_saccades(corrected_eye_vals,all_eye_speed,all_eye_ts,et_params);

par_thresh = 1;
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

bad_sacs = find(isnan(saccade_stop_inds));
used_saccade_set(bad_sacs) = [];
saccade_start_inds(bad_sacs) = [];
saccade_stop_inds(bad_sacs) = [];

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
fix_start_inds(fix_durs <= 0) = [];
fix_stop_inds(fix_durs <= 0) = [];
n_fixs = length(fix_start_inds);
fix_durs = fix_stop_inds-fix_start_inds;

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
pfix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
    cur_inds = pfix_start_inds(ii):(pfix_stop_inds(ii));
    pfix_ids(cur_inds) = ii;
end

%% Create set of TR and XV trials
use_trials = unique(all_trialvec(used_inds));
nuse_trials = length(use_trials);

%select random subset of image frames for xval
poss_frames = unique(all_Seeds(used_inds));
n_xv_frames = round(length(poss_frames)*xv_frac);
xv_frames = randperm(length(poss_frames));
xv_frames(n_xv_frames+1:end) = [];
tr_frames = setdiff(poss_frames,xv_frames);
% save('~/Data/bruce/G103/xv_frames.mat','xv_frames','tr_frames'); 
% load('~/Data/bruce/G103/xv_frames.mat');

xv_inds = find(ismember(all_Seeds(used_inds),xv_frames));
tr_inds = find(ismember(all_Seeds(used_inds),tr_frames));

full_inds = sort([tr_inds; xv_inds]);
% full_inds = tr_inds;

%% COMPILE STIMULUS DATA FROM IMAGE FILES
Nyp = 320;
Nxp = 320;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
xax = linspace(-Nx/2,Nx/2,Nx)*sp_dx; yax = linspace(-Ny/2,Ny/2,Ny)*sp_dx;

rfx = 0.42;
rfy = 0;

im_patch = [rfx-1 rfx+1; rfy-1 rfy+1];
im_patch_pix = im_patch/sp_dx;
RF_patch = [-0.1 0.9; -0.45 0.45];

xpatch_inds = find(xax >= im_patch(1,1) & xax <= im_patch(1,2));
ypatch_inds = find(yax >= im_patch(2,1) & yax <= im_patch(2,2));
cur_len = min(length(xpatch_inds),length(ypatch_inds));
xpatch_inds(cur_len+1:end) = [];
ypatch_inds(cur_len+1:end) = [];

RF_xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
RF_ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));
cur_len = min(length(RF_xpatch_inds),length(RF_ypatch_inds));
RF_xpatch_inds(cur_len+1:end) = [];
RF_ypatch_inds(cur_len+1:end) = [];

NT = length(all_Seeds);
cd ~/James_scripts/bruce/stim_generation/gabor_noise/

all_im_patches = nan(NT,length(ypatch_inds),length(xpatch_inds));
unique_images = unique(all_Seeds(used_inds));
for ii = 1:length(unique_images)
    
    filename = sprintf('Gabor_noise_%.3d.pgm',unique_images(ii));
    IMAGEorg = imread(filename);
    IMAGEorg = double(IMAGEorg); % convert to double format
    IMAGEorg = IMAGEorg - 127.5; %shift range to [-127.5:127.5]
    IMAGE = flipud(IMAGEorg); %flip y
    
    IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
    
    cur_samp_set = find(all_Seeds == unique_images(ii));
    fprintf('Analyzing image %d, %d samps\n',ii,length(cur_samp_set));
    
    for j = 1:length(cur_samp_set)
        cur_patch = IMAGE(ypatch_inds,xpatch_inds);
        all_im_patches(cur_samp_set(j),:,:) = cur_patch;
    end
end

%% NORMALIZE STIMULUS
NT = length(used_inds);

all_im_patches = bsxfun(@minus,all_im_patches,nanmean(all_im_patches));
all_im_patches = bsxfun(@rdivide,all_im_patches,nanstd(all_im_patches));

%% INCORPORATE STIMULUS CORRECTIONS FOR MODEL INITIALIZATION
% if use_measured_pos == 1
%     fprintf('Using eye-coil initialization\n');
%     
%     init_eyepos = corrected_eye_vals_interp(:,1:2);
%     
%     smooth out fast transients in eye signal
%     eye_smooth_sig = round(0.025/dt);
%     interp_inds = [];
%     for ii = 1:n_fixs
%         cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
%         if length(cur_inds) > eye_smooth_sig*5;
%             init_eyepos(used_inds(cur_inds),1) = jmm_smooth_1d_cor(init_eyepos(used_inds(cur_inds),1),eye_smooth_sig,2);
%             init_eyepos(used_inds(cur_inds),2) = jmm_smooth_1d_cor(init_eyepos(used_inds(cur_inds),2),eye_smooth_sig,2);
%         end
%         interp_inds = [interp_inds; cur_inds'];
%     end
%     
%     for ii = 1:n_trials
%         cur_inds = trial_start_inds(ii):trial_end_inds(ii);
%         init_eyepos(used_inds(cur_inds),:) = bsxfun(@minus,init_eyepos(used_inds(cur_inds),:),median(init_eyepos(used_inds(cur_inds),:)));
%     end
%     
%     interp_inds = unique(interp_inds);
%     init_eyepos(used_inds) = interp1(used_inds(interp_inds),init_eyepos(used_inds(interp_inds)),used_inds);
%     
%     max_eyepos = max_shift*sp_dx;
%     init_eyepos(init_eyepos > max_eyepos) = max_eyepos;
%     init_eyepos(init_eyepos < -max_eyepos) = -max_eyepos;
%     
%     shift backwards in time
%     if fixed_delay > 0
%         init_eyepos((end-fixed_delay+1):end,:) = [];
%         init_eyepos = [zeros(fixed_delay,2); init_eyepos];
%     end
%     
% elseif use_measured_pos == 2
%     
%     fprintf('Incorporating random initialization \n');
%     rand_fix_std = 0.2;
%     rand_fix_pos = randn(n_fixs,2)*rand_fix_std;
%     init_eyepos = nan(length(all_t_axis),2);
%     for ii = 1:n_fixs
%         cur_inds = used_inds(fix_start_inds(ii):fix_stop_inds(ii));
%         init_eyepos(cur_inds,:) = repmat(rand_fix_pos(ii,:),length(cur_inds),1);
%     end
%     
%     init_eyepos(isnan(init_eyepos)) = 0;
%     
%     init_eyepos_rnd = round(init_eyepos/sp_dx);
%     init_im_patches = all_im_patches;
%     for ii = 1:NT
%         temp = shift_matrix_Nd(init_im_patches(used_inds(ii),:,:),-init_eyepos_rnd(used_inds(ii),1),3);
%         temp = shift_matrix_Nd(temp,-init_eyepos_rnd(used_inds(ii),2),2);
%         init_im_patches(used_inds(ii),:,:) = temp;
%     end
% end

%% NORMALIZE STIMULUS AND CREATE X-MATRIX

[XXcoords,YYcoords] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

RF_subset = find(XXcoords >= RF_patch(1,1) & XXcoords <= RF_patch(1,2) & ...
    YYcoords >= RF_patch(2,1) & YYcoords <= RF_patch(2,2));

SDIM = length(RF_xpatch_inds);
full_SDIM = length(xpatch_inds);

[Yinds,Tinds,Xinds] = meshgrid(yax(ypatch_inds),1:flen,xax(xpatch_inds));
use_kInds = find(ismember(Xinds,xax(RF_xpatch_inds)) & ismember(Yinds,yax(RF_ypatch_inds)));

full_stim_params = NMMcreate_stim_params([flen full_SDIM full_SDIM],dt);
% if use_measured_pos > 0
%     full_X = create_time_embedding(init_im_patches,full_stim_params);
% else
    full_X = create_time_embedding(all_im_patches,full_stim_params);
% end
full_X = full_X(used_inds,:);

tr_inds = full_inds;

%% DEFINE POINTS TO ALLOW RAPID CHANGES (TRIAL STARTS AFTER SACCADES)
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

trial_ids = nan(NT,1);
for ii = 1:n_trials
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    trial_ids(cur_inds) = ii;
end

%%
x_shifts = -max_shift:dshift:max_shift;
y_shifts = -max_shift:dshift:max_shift;

[Xsh,Ysh] = meshgrid(x_shifts,y_shifts);
[Xsh_inds,Ysh_inds] = meshgrid(1:length(x_shifts),1:length(y_shifts));
SH = [Xsh(:) Ysh(:)];
SH_inds = [Xsh_inds(:) Ysh_inds(:)];
n_shifts = size(SH,1);
zero_frame = find(SH(:,1) == 0 & SH(:,2) == 0);

It = speye(flen);
Ix = speye(full_SDIM);
Iy = speye(full_SDIM);
yshift_mat = cell(length(x_shifts),1);
for xx = 1:length(x_shifts)
    tempx = spdiags( ones(full_SDIM,1), -y_shifts(xx), full_SDIM, full_SDIM);
    yshift_mat{xx} = kron(Iy, kron(tempx, It));
end
xshift_mat = cell(length(y_shifts),1);
for yy = 1:length(x_shifts)
    tempy = spdiags( ones(full_SDIM,1), -x_shifts(yy), full_SDIM, full_SDIM);
    xshift_mat{yy} = kron(tempy, kron(Ix, It));
end
shift_mat = cell(n_shifts,1);
for xx = 1:n_shifts
    temp = xshift_mat{SH_inds(xx,1)}*yshift_mat{SH_inds(xx,2)};
    shift_mat{xx} = temp(:,use_kInds);
end

x_Dshifts = -max_Dshift:dshift:max_Dshift;
y_Dshifts = -max_Dshift:dshift:max_Dshift;

[DXsh,DYsh] = meshgrid(x_Dshifts,y_Dshifts);
[DXsh_inds,DYsh_inds] = meshgrid(1:length(x_Dshifts),1:length(y_Dshifts));
DSH = [DXsh(:) DYsh(:)];
DSH_inds = [DXsh_inds(:) DYsh_inds(:)];
n_Dshifts = size(DSH,1);
Dzero_frame = find(DSH(:,1) == 0 & DSH(:,2) == 0);

Dyshift_mat = cell(length(x_Dshifts),1);
for xx = 1:length(x_Dshifts)
    tempx = spdiags( ones(full_SDIM,1), -y_Dshifts(xx), full_SDIM, full_SDIM);
    Dyshift_mat{xx} = kron(Iy, kron(tempx, It));
end
Dxshift_mat = cell(length(y_Dshifts),1);
for yy = 1:length(x_Dshifts)
    tempy = spdiags( ones(full_SDIM,1), -x_Dshifts(yy), full_SDIM, full_SDIM);
    Dxshift_mat{yy} = kron(tempy, kron(Ix, It));
end
Dshift_mat = cell(n_Dshifts,1);
for xx = 1:n_Dshifts
    temp = Dxshift_mat{DSH_inds(xx,1)}*Dyshift_mat{DSH_inds(xx,2)};
    Dshift_mat{xx} = temp(:,use_kInds);
end
% 
% %shift matrices for images (non-time-embedded)
% max_Tshift = max_shift + max_Dshift;
% x_Tshifts = -max_Tshift:dshift:max_Tshift;
% y_Tshifts = -max_Tshift:dshift:max_Tshift;
% [TXsh,TYsh] = meshgrid(x_Tshifts,y_Tshifts);
% TSH = [TXsh(:) TYsh(:)];
% n_Tshifts = size(TSH,1);

Tshift_mat = cell(n_shifts,1);
for xx = 1:n_shifts
    tempx = spdiags(ones(full_SDIM,1), -SH(xx,2), full_SDIM, full_SDIM);
    xshift = kron(Iy,tempx);
    tempy = spdiags(ones(full_SDIM,1),-SH(xx,1),full_SDIM,full_SDIM);
    yshift = kron(tempy,Ix);
    Tshift_mat{xx} = xshift*yshift;
end

%%
cd(anal_dir)
load(anal_name);

%% compute corrected mat

all_fix_post_mean_cor = nan(NT,2);
% all_fix_post_mean_cor(~isnan(fix_ids),:) = fix_post_mean(fix_ids(~isnan(fix_ids)),:);
all_fix_post_mean_cor(~isnan(fix_ids),:) = squeeze(it_fix_post_mean(end,fix_ids(~isnan(fix_ids)),:));
all_fix_post_mean_cor = interp1(find(~isnan(fix_ids)),all_fix_post_mean_cor(~isnan(fix_ids),:),1:NT);

drift_post_cor = squeeze(drift_post_mean(end,:,:));
for ii = 1:length(trial_start_inds)
    cur_inds = trial_start_inds(ii):trial_end_inds(ii);
    drift_post_cor(cur_inds(1:end-sac_shift),:) = drift_post_cor(cur_inds(sac_shift+1:end),:);
end
drift_post_cor = interp1(find(~isnan(fix_ids)),drift_post_cor(~isnan(fix_ids),:),1:NT);
drift_post_cor(isnan(drift_post_cor)) = 0;

fin_tot_corr = all_fix_post_mean_cor + drift_post_cor;

% all_fix_post_map = round(all_fix_map);
all_shift_stimmat = all_im_patches;
for i=1:NT
%     all_shift_stimmat(used_inds(i),:) = all_shift_stimmat(used_inds(i),:)*Tshift_mat{all_fix_post_map(i)};
      all_shift_stimmat(used_inds(i),:,:) = shift_matrix_Nd(all_shift_stimmat(used_inds(i),:,:),-round(fin_tot_corr(i,1)),3);
      all_shift_stimmat(used_inds(i),:,:) = shift_matrix_Nd(all_shift_stimmat(used_inds(i),:,:),-round(fin_tot_corr(i,2)),2);
end
% cur_X{1} = create_time_embedding(all_shift_stimmat,full_stim_params);
% cur_X{1} = cur_X{1}(used_inds,use_kInds);

% all_dshift_stimmat = all_shift_stimmat;
% for i=1:NT
%     all_dshift_stimmat(used_inds(i),:,:) = shift_matrix_Nd(all_dshift_stimmat(used_inds(i),:,:),-squeeze(round(drift_post_mean(end,i,1))),3);
%     all_dshift_stimmat(used_inds(i),:,:) = shift_matrix_Nd(all_dshift_stimmat(used_inds(i),:,:),-round(drift_post_mean(end,i,2)),2);
% end

% cur_X{1} = create_time_embedding(all_dshift_stimmat,full_stim_params);
cur_X{1} = create_time_embedding(all_shift_stimmat,full_stim_params);
cur_X{1} = cur_X{1}(used_inds,use_kInds);

%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITS
init_stim_params = NMMcreate_stim_params([flen SDIM SDIM],dt);
init_stim_params(2) = NMMcreate_stim_params([n_blocks 1],dt);
if use_sac_kerns
    init_stim_params(3) = NMMcreate_stim_params(n_sac_bins,dt);
    init_stim_params(4) = NMMcreate_stim_params(n_sac_bins,dt);
end

null_stim_params = init_stim_params(2:end);

block_L2 = 1;
silent = 0;
sac_d2t = 100;

base_lambda_d2XT = 300; %300 
base_lambda_L1 = 50; %30

% init_optim_p.optTol = 1e-5; init_optim_p.progTol = 1e-8;

sac_reg_params = NMMcreate_reg_params('lambda_d2T',sac_d2t);
if use_sac_kerns
    null_reg_params = NMMcreate_reg_params('lambda_d2T',[0; sac_d2t; sac_d2t],'lambda_L2',[block_L2; 0; 0]);
else
    null_reg_params = NMMcreate_reg_params('lambda_L2',block_L2);
end

n_squared_filts = 2;
mod_signs = ones(1,n_squared_filts+2);
%if using more than 2 quadratic filters, take the additional ones as
%suppressive
if n_squared_filts > 2
    mod_signs(n_squared_filts+1) = -1;
end
NL_types = [{'lin'} repmat({'quad'},1,n_squared_filts) {'lin'}];
init_d2XT = [20*ones(n_squared_filts+1,1); 0;];
init_L2 = [zeros(n_squared_filts+1,1); block_L2];
init_BCs = repmat(zeros(1,3),n_squared_filts+2,1);
init_mix_prop = repmat([0.5 1 1],n_squared_filts+2,1);
init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_d2XT,'lambda_L2',init_L2,'boundary_conds',init_BCs,'mixing_prop',init_mix_prop);
init_Xtargs = [ones(n_squared_filts+1,1); 2];

init_filts = cell(length(mod_signs),1);
cd(anal_dir);

ss=2;

tic
fprintf('Computing base LLs for SU %d of %d\n',ss,length(su_probes));
cur_tr_inds = tr_inds(~isnan(all_binned_sua(used_inds(tr_inds),ss)));
cur_xv_inds = xv_inds(~isnan(all_binned_sua(used_inds(xv_inds),ss)));
tr_NT = length(cur_tr_inds);
Robs = all_binned_sua(used_inds(cur_tr_inds),ss);

tr_X{1} = cur_X{1}(cur_tr_inds,:);
tr_X{2} = Xblock(cur_tr_inds,:);
if use_sac_kerns
    tr_X{3} = Xsac(cur_tr_inds,:);
    tr_X{4} = Xmsac(cur_tr_inds,:);
    null_mod = NMMinitialize_model(null_stim_params,[1 1 1],repmat({'lin'},1,3),null_reg_params,[1 2 3]);
else
    null_mod = NMMinitialize_model(null_stim_params,[1],{'lin'},null_reg_params,[1]);
end
null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));

gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent);

[LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,Robs,tr_X);
gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
if use_sac_kerns
    gqm1 = NMMadd_NLinput(gqm1,{'lin'},1,3,zeros(n_sac_bins,1),sac_reg_params); %sac term
    gqm1 = NMMadd_NLinput(gqm1,{'lin'},1,4,zeros(n_sac_bins,1),sac_reg_params); %nsac term
end
gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent);

gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent);
examp_cor_mod(ss) = gqm1;
took = toc;
fprintf('Fit took %.2f sec\n',took);
NMMdisplay_model(examp_cor_mod(ss),tr_X,all_binned_sua(used_inds(tr_inds),ss),1);

tic
fprintf('Computing base LLs for SU %d of %d\n',ss,length(su_probes));
cur_tr_inds = tr_inds(~isnan(all_binned_sua(used_inds(tr_inds),ss)));
cur_xv_inds = xv_inds(~isnan(all_binned_sua(used_inds(xv_inds),ss)));
tr_NT = length(cur_tr_inds);
Robs = all_binned_sua(used_inds(cur_tr_inds),ss);

tr_X{1} = full_X(cur_tr_inds,use_kInds);
tr_X{2} = Xblock(cur_tr_inds,:);
if use_sac_kerns
    tr_X{3} = Xsac(cur_tr_inds,:);
    tr_X{4} = Xmsac(cur_tr_inds,:);
    null_mod = NMMinitialize_model(null_stim_params,[1 1 1],repmat({'lin'},1,3),null_reg_params,[1 2 3]);
else
    null_mod = NMMinitialize_model(null_stim_params,[1],{'lin'},null_reg_params,[1]);
end
null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));

gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent);

[LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,Robs,tr_X);
gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
if use_sac_kerns
    gqm1 = NMMadd_NLinput(gqm1,{'lin'},1,3,zeros(n_sac_bins,1),sac_reg_params); %sac term
    gqm1 = NMMadd_NLinput(gqm1,{'lin'},1,4,zeros(n_sac_bins,1),sac_reg_params); %nsac term
end
gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent);

gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent);
examp_uncor_mod(ss) = gqm1;
took = toc;
fprintf('Fit took %.2f sec\n',took);
NMMdisplay_model(examp_uncor_mod(ss),tr_X,all_binned_sua(used_inds(tr_inds),ss),1);

cd(anal_dir)
sname = sprintf('Examp_2d_mods_SU%d_v4_nocoil',ss);
save(sname,'examp_uncor_mod','examp_cor_mod');

%%
cd(anal_dir)
ss=2;
% sname = sprintf('Examp_2d_mods_SU%d_v4',ss);
% load(sname);

flen = 5;
klen = length(examp_cor_mod(ss).mods(1).filtK);
SDIM = sqrt(klen/flen);
nfilts = sum([examp_cor_mod(ss).mods(:).Xtarget] == 1);

fig_width = 4.86; %3.27 4.86 6.83
relheight = 0.5;

int_fac = 8;
[xb,yb] = meshgrid(1:24);
[xi,yi] = meshgrid(1:1/int_fac:24);
xvec = 1:1/int_fac:24;

h1 = figure();
for ii = 1:nfilts
    cur_filt = reshape([examp_cor_mod(ss).mods(ii).filtK],[flen SDIM SDIM]);
    cur_ma = 0.85*max(abs(cur_filt(:)));
    for jj = 1:flen
    subplot(nfilts,flen,(ii-1)*flen + jj);
    INT = interp2(xb,yb,squeeze(cur_filt(jj,:,:)),xi,yi,'cubic');
    imagesc(xvec,xvec,INT);
%     imagesc(squeeze(cur_filt(jj,:,:)));
%     pcolor(squeeze(cur_filt(jj,:,:))); shading interp
    set(gca,'ydir','normal');
    caxis([-cur_ma cur_ma]);
    set(gca,'xtick',[],'ytick',[]);
    xlim([1 18]); ylim([1 18])
    if ii == 1
    title(sprintf('%d ms',dt*1e3*jj-dt/2*1e3));
    end
    end
end
colormap(gray);
figufy(h1);
fname = [fig_dir 'exampSUmod_2d_cor_nocoil'];
exportfig(fname,'width',fig_width,'height',relheight*fig_width,'fontmode','scaled','fontsize',1,'format','pdf');

h2 = figure();
for ii = 1:nfilts
    cur_filt = reshape([examp_uncor_mod(ss).mods(ii).filtK],[flen SDIM SDIM]);
    cur_ma = 0.85*max(abs(cur_filt(:)));
    for jj = 1:flen
    subplot(nfilts,flen,(ii-1)*flen + jj);
    INT = interp2(xb,yb,squeeze(cur_filt(jj,:,:)),xi,yi,'cubic');
    imagesc(xvec,xvec,INT);
%     imagesc(squeeze(cur_filt(jj,:,:)));
%     pcolor(squeeze(cur_filt(jj,:,:)));shading interp
    set(gca,'ydir','normal');
    caxis([-cur_ma cur_ma]);
    set(gca,'xtick',[],'ytick',[]);
    xlim([1 18]); ylim([1 18])
    if ii == 1
    title(sprintf('%d ms',dt*1e3*jj-dt/2*1e3));
    end
    end
end
colormap(gray);
figufy(h2);
fname = [fig_dir 'exampSUmod_2d_uncor_nocoil'];
exportfig(fname,'width',fig_width,'height',relheight*fig_width,'fontmode','scaled','fontsize',1,'format','pdf');

%%
fig_width = 3.27; %3.27 4.86 6.83
relheight = 1.5;
ulag = 4;
az = -30; %azimuthal angle
el = 15; %elevation

rf_pos = Expts{cur_block_set(1)}.Stimvals.rf(1:2);
xpix_ax = ((1:SDIM)-SDIM/2-0.5)*sp_dx + rf_pos(1);
ypix_ax = ((1:SDIM)-SDIM/2-0.5)*sp_dx + rf_pos(2);

int_fac = 4;
xpix_i = ((1:1/int_fac:SDIM)-SDIM/2-0.5)*sp_dx + rf_pos(1);
ypix_i = ((1:1/int_fac:SDIM)-SDIM/2-0.5)*sp_dx + rf_pos(2);

ux = find(xpix_ax >= -0.05 & xpix_ax <= 0.45);
uy = find(ypix_ax >= -0.75 & ypix_ax <= -0.25);

ux_i = find(xpix_i >= -0.05 & xpix_i <= 0.45);
uy_i = find(ypix_i >= -0.75 & ypix_i <= -0.25);

[XX,YY] = meshgrid(xpix_ax,ypix_ax);
[XXi,YYi] = meshgrid(xpix_i,ypix_i);
UU = find(XXi >= -0.05 & XXi <= 0.45 & YYi >= -0.75 & YYi <= -0.25);

cor_filt = reshape([examp_cor_mod(ss).mods(1).filtK],[flen SDIM SDIM]);
cor_filt = squeeze(cor_filt(ulag,:,:));
uncor_filt = reshape([examp_uncor_mod(ss).mods(1).filtK],[flen SDIM SDIM]);
uncor_filt = squeeze(uncor_filt(ulag,:,:));

cur_ma = 1*max(abs(cor_filt(:)));
h3 = figure();
subplot(2,1,1);
INT = interp2(XX,YY,cor_filt,XXi,YYi,'cubic');
h=surf(xpix_i(ux_i),ypix_i(uy_i),INT(uy_i,ux_i));
set(h,'EdgeAlpha',1,'LineWidth',0.25);
% surf(xpix_ax(ux),ypix_ax(uy),cor_filt(uy,ux));
axis tight
zlim([-cur_ma cur_ma]);
xlim([0 0.4]); ylim([-0.7 -0.3]);
xlabel('Horizontal position (deg)');
ylabel('Vertical position (deg)');
zlabel('Filter amp');
title('Corrected');
view(az,el);
subplot(2,1,2);
INT = interp2(XX,YY,uncor_filt,XXi,YYi,'cubic');
h=surf(xpix_i(ux_i),ypix_i(uy_i),INT(uy_i,ux_i));
set(h,'EdgeAlpha',1,'LineWidth',0.25);
% surf(xpix_ax(ux),ypix_ax(uy),uncor_filt(uy,ux)); 
axis tight
zlim([-cur_ma cur_ma]);
xlim([0 0.4]); ylim([-0.7 -0.3]);
view(az,el);
xlabel('Horizontal position (deg)');
ylabel('Vertical position (deg)');
zlabel('Filter amp');
title('Uncorrected');
figufy(h3);
fname = [fig_dir 'exampSUmod_2d_filts_nocoil.pdf'];
exportfig(fname,'width',fig_width,'height',relheight*fig_width,'fontmode','scaled','fontsize',1,'renderer','painters');

% cur_ma = 1*max(abs(cor_filt(:)));
% h3 = figure();
% subplot(2,1,1);
% imagesc(xpix_ax,ypix_ax,cor_filt);
% caxis([-cur_ma cur_ma]);
% set(gca,'ydir','normal');
% xlabel('Horizontal position (deg)');
% ylabel('Vertical position (deg)');
% title('Corrected');
% subplot(2,1,2);
% imagesc(xpix_ax,ypix_ax,uncor_filt); 
% caxis([-cur_ma cur_ma]);
% xlabel('Horizontal position (deg)');
% ylabel('Vertical position (deg)');
% caxis([-cur_ma cur_ma]);
% title('Uncorrected');
% colormap(gray);
% figufy(h3);
% fname = [fig_dir 'exampSUmod_2d_filts_v2.pdf'];
% exportfig(fname,'width',fig_width,'height',relheight*fig_width,'fontmode','scaled','fontsize',1);

%% FOR MUA
% ss=8;
% 
% tic
% fprintf('Computing base LLs for MU %d of %d\n',ss,n_probes);
% cur_tr_inds = tr_inds(~isnan(all_binned_mua(used_inds(tr_inds),ss)));
% cur_xv_inds = xv_inds(~isnan(all_binned_mua(used_inds(xv_inds),ss)));
% tr_NT = length(cur_tr_inds);
% Robs = all_binned_mua(used_inds(cur_tr_inds),ss);
% 
% tr_X{1} = cur_X{1}(cur_tr_inds,:);
% tr_X{2} = Xblock(cur_tr_inds,:);
% if use_sac_kerns
%     tr_X{3} = Xsac(cur_tr_inds,:);
%     tr_X{4} = Xmsac(cur_tr_inds,:);
%     null_mod = NMMinitialize_model(null_stim_params,[1 1 1],repmat({'lin'},1,3),null_reg_params,[1 2 3]);
% else
%     null_mod = NMMinitialize_model(null_stim_params,[1],{'lin'},null_reg_params,[1]);
% end
% null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));
% 
% gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
% gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent);
% 
% [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,Robs,tr_X);
% gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
% gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
% if use_sac_kerns
%     gqm1 = NMMadd_NLinput(gqm1,{'lin'},1,3,zeros(n_sac_bins,1),sac_reg_params); %sac term
%     gqm1 = NMMadd_NLinput(gqm1,{'lin'},1,4,zeros(n_sac_bins,1),sac_reg_params); %nsac term
% end
% gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent);
% 
% gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent);
% examp_cor_MUmod(ss) = gqm1;
% took = toc;
% fprintf('Fit took %.2f sec\n',took);
% NMMdisplay_model(examp_cor_MUmod(ss),tr_X,all_binned_mua(used_inds(tr_inds),ss),1);
% 
% tic
% fprintf('Computing base LLs for MU %d of %d\n',ss,n_probes);
% cur_tr_inds = tr_inds(~isnan(all_binned_mua(used_inds(tr_inds),ss)));
% cur_xv_inds = xv_inds(~isnan(all_binned_mua(used_inds(xv_inds),ss)));
% tr_NT = length(cur_tr_inds);
% Robs = all_binned_mua(used_inds(cur_tr_inds),ss);
% 
% tr_X{1} = full_X(cur_tr_inds,use_kInds);
% tr_X{2} = Xblock(cur_tr_inds,:);
% if use_sac_kerns
%     tr_X{3} = Xsac(cur_tr_inds,:);
%     tr_X{4} = Xmsac(cur_tr_inds,:);
%     null_mod = NMMinitialize_model(null_stim_params,[1 1 1],repmat({'lin'},1,3),null_reg_params,[1 2 3]);
% else
%     null_mod = NMMinitialize_model(null_stim_params,[1],{'lin'},null_reg_params,[1]);
% end
% null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));
% 
% gqm1 = NMMinitialize_model(init_stim_params,mod_signs,NL_types,init_reg_params,init_Xtargs);
% gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent);
% 
% [LL, penLL, pred_rate, G, gint] = NMMmodel_eval(gqm1,Robs,tr_X);
% gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_d2XT',base_lambda_d2XT./var(gint)');
% gqm1 = NMMadjust_regularization(gqm1,find(init_Xtargs==1),'lambda_L1',base_lambda_L1./std(gint)');
% if use_sac_kerns
%     gqm1 = NMMadd_NLinput(gqm1,{'lin'},1,3,zeros(n_sac_bins,1),sac_reg_params); %sac term
%     gqm1 = NMMadd_NLinput(gqm1,{'lin'},1,4,zeros(n_sac_bins,1),sac_reg_params); %nsac term
% end
% gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent);
% 
% gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent);
% examp_uncor_MUmod(ss) = gqm1;
% took = toc;
% fprintf('Fit took %.2f sec\n',took);
% NMMdisplay_model(examp_uncor_MUmod(ss),tr_X,all_binned_mua(used_inds(tr_inds),ss),1);
% 
% 
% 
