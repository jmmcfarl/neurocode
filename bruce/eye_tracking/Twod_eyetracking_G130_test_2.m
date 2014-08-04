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
cluster_dir = [anal_dir '/clustering2'];

mod_data_name = 'Twod_eyecorr_mods3';
anal_name = 'Twod_eyecorr_noCprior3';

recompute_init_mods = 0;
use_measured_pos = 0;
use_sac_kerns = 1;
use_LOOXV = 0;

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

use_coils = [0 0]; %[L R]

%%
xv_frac = 0.15;

flen = 5;
rpt_seed = 1450;

n_fix_inf_it = 5; %3
n_drift_inf_it = 0; %3

fix_noise_sigma = 0.1;
fix_prior_sigma = 0.2; %0.125 start

% fix_prior_sigma = 0.15;
% fix_noise_sigma = 0.1;
% drift_noise_sigma = 0.003;
% drift_prior_sigma = 0.003;
% drift_jump_sigma = 0.05; %0.05 start
% drift_dsf = 3;

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
dsfrac = 3;
sp_dx = Pix2Deg*dsfrac;

max_shift = 12;
dshift = 1;
shifts = -max_shift:dshift:max_shift;
n_shifts = length(shifts);
zero_frame = find(shifts==0);
temp_dx = [-2*max_shift:dshift:2*max_shift];
shift_dx_edges = [temp_dx*sp_dx-dshift*sp_dx/2 temp_dx(end)*sp_dx+dshift*sp_dx/2];
ovshift_dx_edges = [shifts*sp_dx-dshift*sp_dx/2 shifts(end)*sp_dx+dshift*sp_dx/2];

max_Dshift = 8;
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
temp = cellfun(@(x) x.Header.expname,Expts,'uniformoutput',false);

cur_block_set = find(strcmp(temp,'checker.imiRC') | strcmp(temp,'image.imiRC'));
% cur_block_set(cur_block_set > 11 | cur_block_set == 2) = [];
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

trial_toffset = zeros(length(cur_block_set),1);
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
SU_block_probes = nan(length(SU_numbers),length(cur_block_set));
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
    
    cur_suahist = histc(all_su_spk_times{ss},all_t_bin_edges);
    cur_suahist(all_bin_edge_pts) = [];
    cur_id_set = ismember(all_blockvec,cur_blocks);
    all_binned_sua(cur_id_set,ss) = cur_suahist(cur_id_set);
    su_probes(ss) = mode(SU_block_probes(ss,~isnan(SU_block_probes(ss,:))));
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
    cur_mua_inds = all_clust_ids{cc} >= 1;
    
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
[corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
    all_t_axis,all_blockvec,expt_bar_ori,used_inds);

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
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
end

%% Create set of TR and XV trials
use_trials = unique(all_trialvec(used_inds));
nuse_trials = length(use_trials);

% n_xv_trials = round(xv_frac*nuse_trials);
% xv_trials = randperm(nuse_trials);
% xv_trials(n_xv_trials+1:end) = [];
% xv_trials = use_trials(xv_trials);
xv_trials = find(all_trialSe==rpt_seed);
tr_trials = setdiff(use_trials,xv_trials);

tr_inds = find(ismember(all_trialvec(used_inds),tr_trials));
xv_inds = find(ismember(all_trialvec(used_inds),xv_trials));

% full_inds = sort([tr_inds; xv_inds]);
full_inds = tr_inds;

%% COMPILE STIMULUS DATA FROM IMAGE FILES
Nyp = 320;
Nxp = 320;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
xax = linspace(-Nx/2,Nx/2,Nx)*sp_dx; yax = linspace(-Ny/2,Ny/2,Ny)*sp_dx;

rfx = 0.34; 
% rfy = -0.43;
rfy = 0;

im_patch = [rfx-0.82 rfx+0.82; rfy-0.8 rfy+0.8];
im_patch_pix = im_patch/sp_dx;
RF_patch = [rfx-0.6 rfx+0.6; rfy-0.6 rfy+0.6];

xpatch_inds = find(xax >= im_patch(1,1) & xax <= im_patch(1,2));
ypatch_inds = find(yax >= im_patch(2,1) & yax <= im_patch(2,2));

RF_xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
RF_ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

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

%% NORMALIZE STIMULUS AND CREATE X-MATRIX
NT = length(used_inds);

all_im_patches = bsxfun(@minus,all_im_patches,nanmean(all_im_patches));
all_im_patches = bsxfun(@rdivide,all_im_patches,nanstd(all_im_patches));

[XXcoords,YYcoords] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

RF_subset = find(XXcoords >= RF_patch(1,1) & XXcoords <= RF_patch(1,2) & ...
    YYcoords >= RF_patch(2,1) & YYcoords <= RF_patch(2,2));

SDIM = length(RF_xpatch_inds);
% stim_params = NMMcreate_stim_params([flen SDIM SDIM],dt);
% X = create_time_embedding(norm_stim(:,RF_subset),stim_params);
% X = X(used_inds,:);

full_SDIM = length(xpatch_inds);

[Xinds,Tinds,Yinds] = meshgrid(xax(xpatch_inds),1:flen,yax(ypatch_inds));
use_kInds = find(Xinds >= RF_patch(1,1) & Xinds <= RF_patch(1,2) & ...
    Yinds >= RF_patch(2,1) & Yinds <= RF_patch(2,2) & Tinds == 3);

full_stim_params = NMMcreate_stim_params([flen full_SDIM full_SDIM],dt);
full_X = create_time_embedding(all_im_patches,full_stim_params);
flen = 1;

%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITS
init_stim_params = NMMcreate_stim_params([flen SDIM SDIM],dt);
init_stim_params(2) = NMMcreate_stim_params([n_blocks 1],dt);
if use_sac_kerns
    init_stim_params(3) = NMMcreate_stim_params(n_sac_bins,dt);
    init_stim_params(4) = NMMcreate_stim_params(n_sac_bins,dt);
end

null_stim_params = init_stim_params(2:end);

block_L2 = 1;
silent = 1;
sac_d2t = 100;

base_lambda_d2XT = 500;
base_lambda_L1 = 30;

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
init_d2XT = [50*ones(n_squared_filts+1,1); 0;];
init_L2 = [zeros(n_squared_filts+1,1); block_L2];
init_BCs = repmat(zeros(1,3),n_squared_filts+2,1);
init_mix_prop = repmat([0.5 1 1],n_squared_filts+2,1);
init_reg_params = NMMcreate_reg_params('lambda_d2XT',init_d2XT,'lambda_L2',init_L2,'boundary_conds',init_BCs,'mixing_prop',init_mix_prop);
init_Xtargs = [ones(n_squared_filts+1,1); 2];

init_filts = cell(length(mod_signs),1);
cd(anal_dir);

if ~exist(['./' mod_data_name '.mat'],'file') || recompute_init_mods == 1
    tot_nUnits = length(su_probes) + n_probes;
    all_mod_SU = zeros(tot_nUnits,1);
    all_mod_SUnum = zeros(tot_nUnits,1);
    for ss = 1:n_probes
        fprintf('Computing base LLs for MU %d of %d\n',ss,n_probes);
        cur_tr_inds = tr_inds(~isnan(all_binned_mua(used_inds(tr_inds),ss)));
        cur_xv_inds = xv_inds(~isnan(all_binned_mua(used_inds(xv_inds),ss)));
        tr_NT = length(cur_tr_inds);
        xv_NT = length(cur_xv_inds);
        if ~isempty(cur_tr_inds)
            Robs = all_binned_mua(used_inds(cur_tr_inds),ss);
            Robsxv = all_binned_mua(used_inds(cur_xv_inds),ss);
            
            tr_X{1} = full_X(used_inds(cur_tr_inds),use_kInds);
            tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
            xv_X{1} = full_X(used_inds(cur_xv_inds),use_kInds);
            xv_X{2} = Xblock(used_inds(cur_xv_inds),:);
            
            if use_sac_kerns
                tr_X{3} = Xsac(cur_tr_inds,:);
                tr_X{4} = Xmsac(cur_tr_inds,:);
                xv_X{3} = Xsac(cur_xv_inds,:);
                xv_X{4} = Xmsac(cur_xv_inds,:);
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
 
            xvLL = NMMmodel_eval(gqm1,Robsxv,xv_X);
            null_xvLL(ss) = NMMmodel_eval(null_mod,Robsxv,xv_X(2:end));
            all_mod_xvLLimp(ss) = (xvLL-null_xvLL(ss))/log(2);
            
            fprintf('xvLL imp: %.4f\n',all_mod_xvLLimp(ss));
            
            %now refit model using all (usable) data
            cur_tr_inds = sort([cur_tr_inds; cur_xv_inds]);
            Robs = all_binned_mua(used_inds(cur_tr_inds),ss);
            tr_X{1} = full_X(used_inds(cur_tr_inds),use_kInds);
            tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
            if use_sac_kerns
                tr_X{3} = Xsac(cur_tr_inds,:);
                tr_X{4} = Xmsac(cur_tr_inds,:);
            end
            null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));
            all_nullmod(ss) = null_mod;

            gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent);
            all_mod_fits(ss) = gqm1;
            all_mod_fits_withspkNL(ss) = NMMfit_logexp_spkNL(gqm1,Robs,tr_X);
            
            [LL,~,pred_rate] = NMMmodel_eval(all_mod_fits_withspkNL(ss),Robs,tr_X);
            [null_LL(ss),~,null_prate] = NMMmodel_eval(null_mod,Robs,tr_X(2:end));
            [all_mod_R2(ss),all_mod_dev(ss),all_null_dev(ss)] = pseudo_r2(Robs,pred_rate,null_prate);
            all_mod_LLimp(ss) = (LL-null_LL(ss))/log(2);            
        end
    end
    
    for ss = 1:length(su_probes);
        fprintf('Computing base LLs for SU %d of %d\n',ss,length(su_probes));
        all_mod_SU(ss+n_probes) = su_probes(ss);
        all_mod_SUnum(ss+n_probes) = SU_numbers(ss);
        cur_tr_inds = tr_inds(~isnan(all_binned_sua(used_inds(tr_inds),ss)));
        cur_xv_inds = xv_inds(~isnan(all_binned_sua(used_inds(xv_inds),ss)));
        tr_NT = length(cur_tr_inds);
        xv_NT = length(cur_xv_inds);
        if ~isempty(cur_tr_inds)
            Robs = all_binned_sua(used_inds(cur_tr_inds),ss);
            Robsxv = all_binned_sua(used_inds(cur_xv_inds),ss);
            
            tr_X{1} = full_X(used_inds(cur_tr_inds),use_kInds);
            tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
            xv_X{1} = full_X(used_inds(cur_xv_inds),use_kInds);
            xv_X{2} = Xblock(used_inds(cur_xv_inds),:);
            if use_sac_kerns
                tr_X{3} = Xsac(cur_tr_inds,:);
                tr_X{4} = Xmsac(cur_tr_inds,:);
                xv_X{3} = Xsac(cur_xv_inds,:);
                xv_X{4} = Xmsac(cur_xv_inds,:);
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
            
            xvLL = NMMmodel_eval(gqm1,Robsxv,xv_X);
            null_xvLL(ss+n_probes) = NMMmodel_eval(null_mod,Robsxv,xv_X(2:end));
            all_mod_xvLLimp(ss+n_probes) = (xvLL-null_xvLL(ss+n_probes))/log(2);
            
            %now refit model using all (usable) data
            cur_tr_inds = sort([cur_tr_inds; cur_xv_inds]);
            Robs = all_binned_sua(used_inds(cur_tr_inds),ss);
            tr_X{1} = full_X(used_inds(cur_tr_inds),use_kInds);
            tr_X{2} = Xblock(used_inds(cur_tr_inds),:);
            if use_sac_kerns
                tr_X{3} = Xsac(cur_tr_inds,:);
                tr_X{4} = Xmsac(cur_tr_inds,:);
            end
            null_mod = NMMfit_filters(null_mod,Robs,tr_X(2:end));
            all_nullmod(ss+n_probes) = null_mod;

            gqm1 = NMMfit_filters(gqm1,Robs,tr_X,[],[],silent);
            all_mod_fits(ss+n_probes) = gqm1;
            all_mod_fits_withspkNL(ss+n_probes) = NMMfit_logexp_spkNL(gqm1,Robs,tr_X);
            
           [LL,~,pred_rate] = NMMmodel_eval(all_mod_fits_withspkNL(ss+n_probes),Robs,tr_X);
            [null_LL(ss+n_probes),~,null_prate] = NMMmodel_eval(null_mod,Robs,tr_X(2:end));
            [all_mod_R2(ss+n_probes),all_mod_dev(ss+n_probes),all_null_dev(ss+n_probes)] = pseudo_r2(Robs,pred_rate,null_prate);
            all_mod_LLimp(ss+n_probes) = (LL-null_LL(ss+n_probes))/log(2);            
        end
    end
    save(mod_data_name,'all_mod*','all_nullmod','su_probes','null_xvLL','null_LL','*_trials');
    
else
    fprintf('Loading pre-computed initial models\n');
    load(mod_data_name);
end

%% SELECT USABLE UNITS AND make Robs_mat
if xv_frac == 0
    LL_imp_thresh = 5e-3;
    usable_units = find(all_mod_LLimp >= LL_imp_thresh);
else
    LL_imp_thresh = 0;
    usable_units = find(all_mod_xvLLimp > LL_imp_thresh);
end

n_used_sus = sum(all_mod_SU(usable_units) ~= 0);
n_used_mus = sum(all_mod_SU(usable_units) == 0);
fprintf('Using %d SUs and %d MUs for analysis\n',n_used_sus,n_used_mus);
tr_set = usable_units;
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

%don't use separate xv set for eye-tracking
tr_inds = full_inds;

%% DIAGNOSTICS
cur_X{1} = full_X(used_inds,use_kInds);
cur_X{2} = Xblock(used_inds,:);
if use_sac_kerns
    cur_X{3} = Xsac;
    cur_X{4} = Xmsac;
end
full_prates = nan(NT,n_tr_chs);
% full_xvprates = nan(NT,n_tr_chs);
full_nullrates = nan(NT,n_tr_chs);
for cc = 1:n_tr_chs
    [~,~,full_prates(:,cc)] = NMMmodel_eval(all_mod_fits_withspkNL(tr_set(cc)),Robs_mat(:,cc),cur_X);
    %     [~,~,full_xvprates(:,cc)] = NMMmodel_eval(all_xvmod_fits(tr_set(cc)),Robs_mat(:,cc),cur_X);
    [~,~,full_nullrates(:,cc)] = NMMmodel_eval(all_nullmod(tr_set(cc)),Robs_mat(:,cc),cur_X(2:end));
end
full_modLL = Robs_mat.*log(full_prates) - full_prates;
full_nullLL = Robs_mat.*log(full_nullrates) - full_nullrates;
full_LLimp = full_modLL-full_nullLL;

trial_blocknums = nan(nuse_trials,1);
trial_LLimp = nan(nuse_trials,n_tr_chs);
trial_meanrate = nan(nuse_trials,n_tr_chs);
trial_nspks = nan(nuse_trials,n_tr_chs);
trial_durs = nan(nuse_trials,1);
trial_eyeerr = nan(nuse_trials,1);
for tt = 1:nuse_trials
    cur_used_inds = find(all_trialvec(used_inds) == use_trials(tt));
    trial_LLimp(tt,:) = nanmean(full_LLimp(cur_used_inds,:));
    trial_meanrate(tt,:) = nanmean(Robs_mat(cur_used_inds,:));
    trial_nspks(tt,:) = nansum(Robs_mat(cur_used_inds,:));
    trial_blocknums(tt) = unique(all_blockvec(used_inds(cur_used_inds)));
    trial_durs(tt) = length(cur_used_inds);
    trial_eyeerr(tt) = mean(corrected_eye_vals_interp(used_inds(cur_used_inds),2));
end

block_LLimp = nan(n_blocks,n_tr_chs);
block_meanrate = nan(n_blocks,n_tr_chs);
block_nspks = nan(n_blocks,n_tr_chs);
for tt = 1:n_blocks
    cur_used_inds = find(all_blockvec(used_inds) == tt);
    block_LLimp(tt,:) = nanmean(full_LLimp(cur_used_inds,:));
    block_meanrate(tt,:) = nanmean(Robs_mat(cur_used_inds,:));
    block_nspks(tt,:) = nansum(Robs_mat(cur_used_inds,:));
end

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
% fix_prior = diff(erf(ovshift_dx_edges/(sqrt(2)*fix_prior_sigma)));
% fix_prior = log(fix_prior/sum(fix_prior));
% eps = -1e3;
% fix_prior(fix_prior < eps) = eps;

x_shifts = -max_shift:dshift:max_shift;
y_shifts = -max_shift:dshift:max_shift;

[Xsh,Ysh] = meshgrid(x_shifts,y_shifts);
[Xsh_inds,Ysh_inds] = meshgrid(1:length(x_shifts),1:length(y_shifts));
SH = [Xsh(:) Ysh(:)];
SH_inds = [Xsh_inds(:) Ysh_inds(:)];
n_shifts = size(SH,1);
zero_frame = find(SH(:,1) == 0 & SH(:,2) == 0);

% temp_dx = [-2*max_shift:dshift:2*max_shift];
% shift_dx_edges = [temp_dx*sp_dx-dshift*sp_dx/2 temp_dx(end)*sp_dx+dshift*sp_dx/2];
% ovshift_dx_edges = [shifts*sp_dx-dshift*sp_dx/2 shifts(end)*sp_dx+dshift*sp_dx/2];

%generate shift matrices. Must be applied to the stimulus (not the filters)
flen = 5;
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

% Dshift_mat = cell(n_Dshifts,1);
% for xx = 1:n_Dshifts
%     temp = spdiags( ones(full_nPix_us,1), -Dshifts(xx), full_nPix_us, full_nPix_us);
%     temp = kron(temp,It);
%     Dshift_mat{xx} = temp(:,use_kInds_up);
% end

%shift matrices for images (non-time-embedded)    
n_Tshifts = length(n_shifts);
Tshift_mat = cell(n_shifts,1);
for xx = 1:n_shifts
    tempx = spdiags(ones(full_SDIM,1), -SH(xx,2), full_SDIM, full_SDIM);
    xshift = kron(Iy,tempx);
    tempy = spdiags(ones(full_SDIM,1),-SH(xx,1),full_SDIM,full_SDIM);
    yshift = kron(tempy,Ix);
    Tshift_mat{xx} = xshift*yshift;
end

% Tshift_mat = cell(n_shifts,1);
% for xx = 1:n_shifts
%     temp = xshift_mat{find(x_shifts == SH(xx,1))}*yshift_mat{find(y_shifts == SH(xx,2))};
% Tshift_mat{xx} = temp(:,use_kInds);
% end

%%
% shifts = -max_shift:dshift:max_shift;
% n_shifts = length(shifts);
% 
% fix_prior = diff(erf(ovshift_dx_edges/(sqrt(2)*fix_prior_sigma)));
% fix_prior = log(fix_prior/sum(fix_prior));
% eps = -1e3;
% fix_prior(fix_prior < eps) = eps;
% 
% %generate shift matrices. Must be applied to the stimulus (not the filters)
% It = speye(flen);
% Iy = speye(full_SDIM);
% shift_mat = cell(n_shifts,1);
% for xx = 1:n_shifts
%     temp = spdiags( ones(full_SDIM,1), -shifts(xx), full_SDIM, full_SDIM);
%     temp = kron(Iy, kron(temp,It));
%     shift_mat{xx} = temp(:,use_kInds);
% end
% 
% %shift matrices for images (non-time-embedded)    
% Tshift_mat = cell(n_shifts,1);
% for xx = 1:n_shifts
%     temp = kron(Iy,spdiags( ones(full_SDIM,1), -Tshifts(xx), full_SDIM, full_SDIM));
%     temp = kron(Iy,spdiags( ones(full_SDIM,1), -Tshifts(xx), full_SDIM, full_SDIM));
%     
%     Tshift_mat{xx} = temp*temp1;
% end

%%
fix_prior = -sum((SH*sp_dx).^2,2)/(2*fix_prior_sigma^2);
fix_prior = bsxfun(@minus,fix_prior,logsumexp(fix_prior)); %normalize



%%
flen = 1;
su_inds = find(all_mod_SU(tr_set) > 0);
klen = flen*SDIM^2;

%% ITERATE FIXATION-BASED CORRECTIONS

it_mods{1} = all_mod_fits;
it_mods_spkNL{1} = all_mod_fits_withspkNL;
it_LLimp(1,:) = all_mod_LLimp;
it_R2(1,:) = all_mod_R2;
it_dev(1,:) = all_mod_dev;
it_fix_sigma(1) = fix_prior_sigma;
if use_LOOXV
    it_LLimp_LOO(length(su_inds),1,:) = all_mod_LLimp;
    it_R2_LOO(length(su_inds),1,:) = all_mod_R2;
    it_dev_LOO(length(su_inds),:) = all_mod_dev;
    for xv = 1:length(su_inds)
       it_mods_LOO{xv,1} = all_mod_fits;
       it_mods_spkNL_LOO{xv,1} = all_mod_fits_withspkNL;
    end
end
for nn = 1:n_fix_inf_it
    fprintf('Inferring fixation corrections, iter %d of %d\n',nn,n_fix_inf_it);
    
    %% PREPROCESS MODEL COMPONENTS
    filt_bank = zeros(n_tr_chs,klen,n_squared_filts+1);
    lin_kerns = nan(n_tr_chs,n_blocks);
    if use_sac_kerns
        sac_kerns = nan(n_tr_chs,n_sac_bins);
        msac_kerns = nan(n_tr_chs,n_sac_bins);
    end
    mod_spkNL_params = nan(n_tr_chs,3);
    for ss = 1:n_tr_chs
        cur_Xtargs = [it_mods{nn}(tr_set(ss)).mods(:).Xtarget];
        cur_k = [it_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 1).filtK];
        n_used_filts = size(cur_k,2);
        filt_bank(ss,:,1:n_used_filts) = cur_k;
        mod_spkNL_params(ss,:) = it_mods_spkNL{nn}(tr_set(ss)).spk_NL_params;
        lin_kerns(ss,:) = it_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 2).filtK;
        if use_sac_kerns
            sac_kerns(ss,:) = it_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 3).filtK;
            msac_kerns(ss,:) = it_mods{nn}(tr_set(ss)).mods(cur_Xtargs == 4).filtK;
        end
    end
    filt_bank = permute(filt_bank,[2 1 3]);
    
    %indicator predictions
    block_out = Xblock(used_inds,:)*lin_kerns';
    if use_sac_kerns
        sac_out = Xsac*sac_kerns';
        msac_out = Xmsac*msac_kerns';
    end
    
    %% ESTIMATE LL for each shift in each stimulus frame
    cur_Xmat = full_X(used_inds,:);
    %precompute LL at all shifts for all units
    frame_LLs = nan(length(used_inds),n_shifts);
    for xx = 1:n_shifts
        fprintf('Shift %d of %d\n',xx,n_shifts);
        cur_stim_shift = cur_Xmat*shift_mat{xx};
        
        %outputs of stimulus models at current X-matrix shift
        gfuns = ones(length(used_inds),n_tr_chs);
        gfuns = bsxfun(@times,gfuns,mod_spkNL_params(:,1)');
        gfuns = gfuns + cur_stim_shift*squeeze(filt_bank(:,:,1));
        for ff = 2:(n_squared_filts+1)
            gfuns = gfuns + mod_signs(ff)*(cur_stim_shift*squeeze(filt_bank(:,:,ff))).^2;
        end
        
        %add contributions from extra lin kernels
        gfuns = gfuns + block_out;
        if use_sac_kerns
            gfuns = gfuns + sac_out + msac_out;
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
        
        frame_LLs(:,xx) = squeeze(nansum(Robs_mat.*log(pred_rate) - pred_rate,2));
    end
    
    %% INFER MICRO-SAC SEQUENCE
    fix_LLs = nan(n_fixs,n_shifts);
    temp_rate = nan(n_fixs,length(tr_set));
    temp_nspk = nan(n_fixs,length(tr_set));
    for ii = 1:n_fixs
        cur_inds = pfix_start_inds(ii):pfix_stop_inds(ii);
        fix_LLs(ii,:) = sum(frame_LLs(cur_inds,:));
%         temp_rate(ii,:) = mean(Robs_mat(cur_inds,:));
%         temp_nspk(ii,:) = sum(Robs_mat(cur_inds,:));
    end
    
    if all(use_coils==0)
%                 lgamma = bsxfun(@plus,fix_LLs,fix_prior');
        lgamma = fix_LLs;
        lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    else
        lalpha=zeros(n_fixs,n_shifts);
        lbeta = zeros(n_fixs,n_shifts);
        %compute forward messages
        lalpha(1,:) = fix_prior' + fix_LLs(1,:);
        for t=2:n_fixs
            %             gprobs = diff(erf((shift_dx_edges + measured_fix_deltas(t))/(fix_noise_sigma*sqrt(2))));
            %             cur_lA = log(toeplitz(gprobs(n_shifts:end)',gprobs(n_shifts:-1:1)'));
            %             cur_lA(isinf(cur_lA)) = eps;
            %             cur_lA = bsxfun(@plus,cur_lA,fix_prior);
            %             cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2));
            
            cdist = pdist2(SH*sp_dx,[SH(:,1)*sp_dx - measured_fix_deltas(t,1) SH(:,2)*sp_dx - measured_fix_deltas(t,2)]);
            cur_lA = -cdist.^2/(2*fix_noise_sigma^2);
            cur_lA = bsxfun(@plus,cur_lA,fix_prior);
            cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2));
                         
            lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + fix_LLs(t,:);
        end
        
        %compute backward messages
        lbeta(n_fixs,:)=log(ones(1,n_shifts));
        for t=n_fixs-1:-1:1
%             gprobs = diff(erf((shift_dx_edges + measured_fix_deltas(t+1))/(fix_noise_sigma*sqrt(2))));
%             cur_lA = log(toeplitz(gprobs(n_shifts:end)',gprobs(n_shifts:-1:1)'));
%             cur_lA(isinf(cur_lA)) = eps;
%             cur_lA = bsxfun(@plus,cur_lA,fix_prior);            
%             cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2));

            cdist = pdist2(SH*sp_dx,[SH(:,1)*sp_dx - measured_fix_deltas(t+1,1) SH(:,2)*sp_dx - measured_fix_deltas(t+1,2)]);
            cur_lA = -cdist.^2/(2*fix_noise_sigma^2);
            cur_lA = bsxfun(@plus,cur_lA,fix_prior);
            cur_lA = bsxfun(@minus,cur_lA,logsumexp(cur_lA,2));

            lf1 = lbeta(t+1,:) + fix_LLs(t+1,:);
            lbeta(t,:) = logmulexp(lf1,cur_lA');
        end
        lgamma= lalpha + lbeta;
        lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    end
    
   gamma = exp(lgamma);
    
   fix_post_mean(:,1) = sum(bsxfun(@times,gamma,SH(:,1)'),2);
   fix_post_mean(:,2) = sum(bsxfun(@times,gamma,SH(:,2)'),2);

   cur_xdiff = bsxfun(@minus,fix_post_mean(:,1)',SH(:,1)).^2';
   fix_post_std(:,1) = sqrt(sum(cur_xdiff.*gamma,2));
   cur_ydiff = bsxfun(@minus,fix_post_mean(:,2)',SH(:,2)).^2';
   fix_post_std(:,2) = sqrt(sum(cur_ydiff.*gamma,2));
    
   [~,fix_post_max] = max(gamma,[],2);
   
    %back-project saccade-times
    all_fix_post_mean_cor = nan(NT,2);
    all_fix_post_mean_cor(~isnan(fix_ids),:) = fix_post_mean(fix_ids(~isnan(fix_ids)),:);
    all_fix_post_mean_cor = interp1(find(~isnan(fix_ids)),all_fix_post_mean_cor(~isnan(fix_ids),:),1:NT);
    
    all_fix_map = nan(NT,1);
    all_fix_map(~isnan(fix_ids)) = fix_post_max(fix_ids(~isnan(fix_ids)));
    all_fix_map = interp1(find(~isnan(fix_ids)),all_fix_map(~isnan(fix_ids)),1:NT);
    all_fix_map(isnan(all_fix_map)) = zero_frame;
    %%
    all_fix_post_map = round(all_fix_map);
    all_shift_stimmat = all_im_patches;
    for i=1:NT
        all_shift_stimmat(used_inds(i),:) = all_im_patches(used_inds(i),:)*Tshift_mat{all_fix_post_map(i)};
    end
    cur_X{1} = create_time_embedding(all_shift_stimmat,full_stim_params);
    
    %% REFIT ALL CELLS
    cur_X{1} = cur_X{1}(used_inds,use_kInds);
    cur_X{2} = Xblock(used_inds,:);
    if use_sac_kerns
        cur_X{3} = Xsac;
        cur_X{4} = Xmsac;
    end
    
    silent = 1;
    for ss = 1:length(tr_set)
        cur_cell = tr_set(ss);
        fprintf('Refitting model for tr cell %d of %d\n',ss,length(tr_set));
        cur_unit_ind = find(tr_set == cur_cell);
        cur_tr_uset = tr_inds(~isnan(Robs_mat(tr_inds,cur_unit_ind)));
        
        tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
        
        it_mods{nn+1}(cur_cell) = it_mods{nn}(cur_cell);
        it_mods{nn+1}(cur_cell) = NMMfit_filters(it_mods{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),...
            tr_X,[],[],silent); %fit stimulus filters
        
        %refit spk NL
        it_mods_spkNL{nn+1}(cur_cell) = NMMfit_logexp_spkNL(it_mods{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
        
        [newLL,~,new_prate] = NMMmodel_eval(it_mods_spkNL{nn+1}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
        [~,~,null_prate] = NMMmodel_eval(all_nullmod(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X(2:end));
        [it_R2(nn+1,cur_cell),it_dev(nn+1,cur_cell)] = pseudo_r2(Robs_mat(cur_tr_uset,cur_unit_ind),new_prate,null_prate);
         
        it_LLimp(nn+1,cur_cell) = (newLL - null_LL(cur_cell))/log(2);
        if nn > 1
            fprintf('Original: %.4f  Prev: %.4f  New: %.4f\n',all_mod_LLimp(cur_cell),it_LLimp(nn,cur_cell),it_LLimp(nn+1,cur_cell));
        else
            fprintf('Original: %.4f  New: %.4f\n',all_mod_LLimp(cur_cell),it_LLimp(nn+1,cur_cell));
        end
    end
    
end

%%
cd(anal_dir)
save(anal_name,'it_*','fix_*','tr_set')

%% SAVE EYE-TRACKING RESULTS
% % et_params = struct('beg_buffer',beg_buffer,'end_buffer',end_buffer,'min_trial_dur',min_trial_dur,'bar_ori',bar_ori,...
% %     'use_nPix',use_nPix,'flen',flen,'dt',dt,'drift_jump_sigma',drift_jump_sigma,'drift_prior_sigma',drift_prior_sigma,...
% %     'fix_prior_sigma',fix_prior_sigma,'fix_noise_sigma',fix_noise_sigma,'drift_noise_sigma',drift_noise_sigma,...
% %     'drift_dsf',drift_dsf,'n_fix_inf_it',n_fix_inf_it,'n_drift_inf_it',n_drift_inf_it,'use_sac_kerns',use_sac_kerns,'shifts',shifts,...
% %     'use_measured_pos',use_measured_pos,'sac_bincents',sac_bincents,'spatial_usfac',spatial_usfac,'sac_shift',sac_shift);
% % 
% % et_used_inds = used_inds;
% % et_tr_set = tr_set;
% % et_tr_trials = tr_trials;
% % et_xv_trials = xv_trials;
% % et_saccade_inds = saccade_start_inds;
% % cd(anal_dir);
% % save(anal_name,'it_*','drift_post_*','fix_ids','dit_*','et_used_inds','et_tr_set','et_saccade_inds','et_params','et_tr_trials','et_xv_trials');
% 
% %% CHECK LL IMPROVEMENTS SPECIFICALLY FOR SIM SAC BLOCKS
% sim_sac_tr_inds = tr_inds(ismember(all_trialvec(tr_inds),sim_sac_expts));
% for xv = 1:length(su_inds)
%     
%     fix_post_cor = nan(NT,1);
%     fix_post_cor(~isnan(fix_ids)) = squeeze(it_fix_post_mean_LOO(xv,end,fix_ids(~isnan(fix_ids))));
%     fix_post_cor = interp1(find(~isnan(fix_ids)),fix_post_cor(~isnan(fix_ids)),1:NT);
%     
%     %back project drift (within-fixation) by sac_shift
%     drift_post_cor = squeeze(drift_post_mean_LOO(xv,end,:));
%     for ii = 1:n_fixs
%         cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
%         if length(cur_inds) > sac_shift
%             drift_post_cor(cur_inds(1:end-sac_shift+1)) = drift_post_cor(cur_inds(sac_shift:end));
%         end
%     end
%     
%     all_post_cor = round(fix_post_cor+drift_post_cor') + max_Tshift + 1;
%     
%     %RECOMPUTE XMAT
%     all_shift_stimmat_up = all_stimmat_up;
%     for i=1:NT
%         all_shift_stimmat_up(used_inds(i),:) = all_stimmat_up(used_inds(i),:)*Tshift_mat{all_post_cor(i)};
%     end
%     cur_X{1} = create_time_embedding(all_shift_stimmat_up,stim_params_us);
%     cur_X{1} = cur_X{1}(used_inds,use_kInds_up);
%     
%     cur_cell = tr_set(su_inds(xv));
%     cur_unit_ind = find(tr_set == cur_cell);
%     cur_tr_uset = sim_sac_tr_inds(~isnan(Robs_mat(sim_sac_tr_inds,cur_unit_ind)));
%     
%     tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
%     LL_cor_simsac(xv) = NMMmodel_eval(dit_mods_spkNL_LOO{xv,end}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
%     LL_cor_simsac(xv) = LL_cor_simsac(xv) - null_LL(tr_set(su_inds(xv)));
%     
%     cur_X{1} = all_Xmat_us(used_inds,use_kInds_up);
%     tr_X = get_Xcell_tInds(cur_X,cur_tr_uset);
%     LL_uncor_simsac(xv) = NMMmodel_eval(dit_mods_spkNL_LOO{xv,end}(cur_cell),Robs_mat(cur_tr_uset,cur_unit_ind),tr_X);
%     LL_uncor_simsac(xv) = LL_uncor_simsac(xv) - null_LL(tr_set(su_inds(xv)));
%         
% end
% 
%%
fin_fix_corr = nan(NT,2);
fin_fix_std = nan(NT,2);
fin_fix_corr(~isnan(fix_ids),:) = fix_post_mean(fix_ids(~isnan(fix_ids)),:);
fin_fix_corr = interp1(find(~isnan(fix_ids)),fin_fix_corr(~isnan(fix_ids),:),1:NT);
fin_fix_std(~isnan(fix_ids),:) = fix_post_std(fix_ids(~isnan(fix_ids)),:);
fin_fix_std = interp1(find(~isnan(fix_ids)),fin_fix_std(~isnan(fix_ids),:),1:NT);

fin_fix_corr = fin_fix_corr*sp_dx;
fin_fix_std = fin_fix_std*sp_dx;

% fin_drift_corr = drift_post_mean(end,:)*sp_dx;
% fin_drift_std = drift_post_std(end,:)*sp_dx;
% min_fix_dur = 0.15;
% fix_inds = [];
% long_fix_inds = [];
% for ii = 1:n_fixs
%     cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
%     if length(cur_inds) > sac_shift
%         fin_drift_corr(cur_inds(1:end-sac_shift+1)) = fin_drift_corr(cur_inds(sac_shift:end));
%         fin_drift_std(cur_inds(1:end-sac_shift+1)) = fin_drift_std(cur_inds(sac_shift:end));
%     end
%     fix_inds = [fix_inds cur_inds];
%     if length(cur_inds)*dt >= min_fix_dur
%         long_fix_inds = [long_fix_inds cur_inds];
%     end
% end
% 
% fin_tot_corr = fin_fix_corr + fin_drift_corr;
% fin_tot_std = sqrt(fin_fix_std.^2 + fin_drift_std.^2);
% % fin_tot_corr = interp1(find(~isnan(fix_ids)),fin_tot_corr(~isnan(fix_ids)),1:NT);
% % fin_tot_std = interp1(find(~isnan(fix_ids)),fin_tot_std(~isnan(fix_ids)),1:NT);
% 
% if use_measured_pos==1
%     fin_tot_corr = fin_tot_corr + init_eyepos(used_inds)';
% end
% %%
% 
% measured_seqL = corrected_eye_vals_interp(used_inds,2);
% measured_seqR = corrected_eye_vals_interp(used_inds,4);
% 
% min_fix_dur = 0;
% % inferred_drift = nan(size(fin_tot_corr));
% % measured_driftL = nan(size(fin_tot_corr));
% % measured_driftR = nan(size(fin_tot_corr));
% inferred_fix_avg = nan(n_fixs,1);
% measured_fix_avgL = nan(n_fixs,1);
% measured_fix_avgR = nan(n_fixs,1);
% for ii = 1:n_fixs
%     cur_inds = fix_start_inds(ii):fix_stop_inds(ii);
%     if length(cur_inds)*dt >= min_fix_dur 
% %         cur_inf = fin_tot_corr(cur_inds);
% %         inferred_fix_avg(ii) = median(fin_tot_corr(cur_inds));
% %         inferred_drift(cur_inds) = cur_inf - inferred_fix_avg(ii);
%         
%         measured_fix_avgL(ii) = median(measured_seqL(cur_inds));
%         measured_fix_avgR(ii) = median(measured_seqR(cur_inds));
% %         measured_driftL(cur_inds) = measured_seqL(cur_inds) - measured_fix_avgL(ii);
% %         measured_driftR(cur_inds) = measured_seqR(cur_inds) - measured_fix_avgR(ii);
%     end
% end
% 
% % u = find(~isnan(measured_driftL) & ~isnan(inferred_drift));
% % [drift_corrs,drif_pvals] = corr([measured_driftL(u)' measured_driftR(u)'],inferred_drift(u)','type','spearman');
% % u = find(~isnan(measured_fix_avgL) & ~isnan(inferred_fix_avg));
% % [fix_corrs,fix_pvals] = corr([measured_fix_avgL(u) measured_fix_avgR(u)],inferred_fix_avg(u),'type','spearman');
% % [tot_corrs,tot_pvals] = corr([measured_seqL(usable_inds) measured_seqR(usable_inds)],fin_tot_corr(usable_inds)','type','spearman');
% % [tot_corrs_pear,tot_pvals_pear] = corr([measured_seqL(usable_inds) measured_seqR(usable_inds)],fin_tot_corr(usable_inds)','type','pearson');
% 
% %% QUANTIFY LL IMPROVEMENTS
% sus = tr_set(all_mod_SU(tr_set) > 0);
% mus = tr_set(all_mod_SU(tr_set) == 0);
% 
% it_LLimp_LOO(1,:) = it_LLimp(1,:);
% 
% fix_LL_imp = bsxfun(@rdivide,it_LLimp,it_LLimp(1,:));
% full_LL_imp = bsxfun(@rdivide,dit_LLimp,it_LLimp(1,:));
% 
% fix_LL_imp_LOO = bsxfun(@rdivide,it_LLimp_LOO(:,sus),it_LLimp_LOO(1,sus));
% full_LL_imp_LOO = bsxfun(@rdivide,dit_LLimp_LOO(:,sus),it_LLimp_LOO(1,sus));
% 
% % figure
% % plot(it_xvLLimp(1,mus),it_xvLLimp(end,mus),'k.','markersize',8)
% % hold on
% % plot(it_xvLLimp(1,sus),it_xvLLimp(end,sus),'r.','markersize',8)
% % xlim([0 0.6]); ylim([0 0.6])
% % line([0 0.6],[0 0.6],'color','k')
% % xlabel('Initial xvLL (bits/spk)','fontsize',12);
% % ylabel('Final xvLL (bits/spk)','fontsize',12);
% % legend('MU','SU','Location','Southeast');
% % box off
% % set(gca,'fontname','arial','fontsize',10);
% % fillPage(gcf,'papersize',[4 4]);
% %
% % figure
% % plot(dit_xvLLimp(1,mus),dit_xvLLimp(end,mus),'k.','markersize',8)
% % hold on
% % plot(dit_xvLLimp(1,sus),dit_xvLLimp(end,sus),'r.','markersize',8)
% % xlim([0 0.75]); ylim([0 0.75])
% % line([0 0.75],[0 0.75],'color','k')
% % xlabel('Initial xvLL (bits/spk)','fontsize',12);
% % ylabel('Final xvLL (bits/spk)','fontsize',12);
% % legend('MU','SU','Location','Southeast');
% % box off
% % set(gca,'fontname','arial','fontsize',10);
% % fillPage(gcf,'papersize',[4 4]);
% 
% figure;hold on
% errorbar(1:(n_fix_inf_it+1),nanmean(fix_LL_imp(:,mus),2),nanstd(fix_LL_imp(:,mus),[],2)/sqrt(length(mus)));
% errorbar(1:(n_fix_inf_it+1),nanmean(fix_LL_imp(:,sus),2),nanstd(fix_LL_imp(:,sus),[],2)/sqrt(length(sus)),'r');
% 
% errorbar((n_fix_inf_it+2):(n_fix_inf_it+n_drift_inf_it+2),nanmean(full_LL_imp(:,mus),2),nanstd(full_LL_imp(:,mus),[],2)/sqrt(length(mus)));
% errorbar((n_fix_inf_it+2):(n_fix_inf_it+n_drift_inf_it+2),nanmean(full_LL_imp(:,sus),2),nanstd(full_LL_imp(:,sus),[],2)/sqrt(length(sus)),'r');
% xlabel('Iterations','fontsize',12);
% ylabel('LL ratio','fontsize',12);
% box off
% set(gca,'fontname','arial','fontsize',10);
% fillPage(gcf,'papersize',[4 4]);
% 
% figure;hold on
% errorbar(1:(n_fix_inf_it+1),nanmean(fix_LL_imp(:,sus),2),nanstd(fix_LL_imp(:,sus),[],2)/sqrt(length(sus)),'k');
% errorbar((n_fix_inf_it+2):(n_fix_inf_it+n_drift_inf_it+2),nanmean(full_LL_imp(:,sus),2),nanstd(full_LL_imp(:,sus),[],2)/sqrt(length(sus)),'k');
% errorbar(1:(n_fix_inf_it+1),nanmean(fix_LL_imp_LOO,2),nanstd(fix_LL_imp_LOO,[],2)/sqrt(length(sus)),'r');
% errorbar((n_fix_inf_it+2):(n_fix_inf_it+n_drift_inf_it+2),nanmean(full_LL_imp_LOO,2),nanstd(full_LL_imp_LOO,[],2)/sqrt(length(sus)),'r');
% xlabel('Iterations','fontsize',12);
% ylabel('LL ratio','fontsize',12);
% box off
% set(gca,'fontname','arial','fontsize',10);
% fillPage(gcf,'papersize',[4 4]);
% 
%%
close all
n_trials = length(unique(all_trialvec));
for tt = 1:n_trials
    % for tt = [96 137 154 179 376 409]
    uu = find(all_trialvec(used_inds) == tt);
    if ~isempty(uu)
        bt = all_t_axis(used_inds(uu(1)));
        et = all_t_axis(used_inds(uu(end)));
        dur = et-bt;
        if dur > 3.5
            hold on
            h1=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_fix_corr(uu,1),fin_fix_std(uu,1),{'color','m'});
            %             h2=shadedErrorBar(all_t_axis(used_inds(uu))-bt,fin_tot_corr(uu),fin_tot_std(uu),{'color','k'});
            h3=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),1),'r','linewidth',2);
            h4=plot(all_t_axis(used_inds(uu))-bt,corrected_eye_vals_interp(used_inds(uu),3)-median(corrected_eye_vals_interp(used_inds(uu),4)),'color',[0.2 0.8 0.2],'linewidth',2);
            if use_measured_pos==1
                plot(all_t_axis(used_inds(uu))-bt,init_eyepos(used_inds(uu)),'c','linewidth',2)
            end
            %             plot(all_t_axis(used_inds(uu))-bt,nanmean(Robs_mat(uu,:),2)/5,'k');
            
            %             legend([h1.mainLine h2.mainLine h3 h4],{'Fixation corrections','Drift corrections','Left eye','Right eye'})
            xlim([0 dur]);
            ylim([-0.5 0.5]);
            xlabel('Time (s)','fontsize',10);
            ylabel('Orthoganol position (deg)','fontsize',10);
            title(sprintf('Trial %d',tt));
            set(gca,'fontsize',8,'fontname','arial');
            fillPage(gcf,'papersize',[8 5]);
            pause
            clf
        end
    end
end
% 
% %%
% close all
% f1 = figure();
% f2 = figure();
% for ss = 1:length(tr_set)
% % sbeg = find(all_mod_SU(tr_set) > 0,1);
% % for ss = sbeg:length(tr_set)
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
