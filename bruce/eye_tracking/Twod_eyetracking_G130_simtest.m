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

mod_data_name = 'Twod_eyecorr_mods7';
anal_name = 'Twod_eyecorr_noCprior7';

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

flen = 1;
rpt_seed = 1450;

n_fix_inf_it = 3; %3
n_drift_inf_it = 0; %3

fix_noise_sigma = 0.1;
fix_prior_sigma = 0.2; %0.125 start

min_trial_dur = 1.75;

fixed_delay = 2;
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
% [corrected_eye_vals,corrected_eye_vals_interp]  = get_corrected_ET_data(Expts(cur_block_set),all_eye_vals,all_eye_ts,...
%     all_t_axis,all_blockvec,expt_bar_ori,used_inds);
corrected_eye_vals = bsxfun(@minus,all_eye_vals,nanmedian(all_eye_vals));
corrected_eye_vals_interp = interp1(all_eye_ts,corrected_eye_vals,all_t_axis);
corrected_eye_vals_interp(isnan(corrected_eye_vals_interp)) = 0;

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

% rfx = 0.34; 
rfx = 0.42; 
% rfy = -0.43;
rfy = 0;

im_patch = [rfx-1 rfx+1; rfy-1 rfy+1];
im_patch_pix = im_patch/sp_dx;
% RF_patch = [rfx-0.3 rfx+0.9; rfy-0.6 rfy+0.6];
RF_patch = [rfx-0.7 rfx+0.7; rfy-0.7 rfy+0.7];

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

%% NORMALIZE STIMULUS AND CREATE X-MATRIX
NT = length(used_inds);

all_im_patches = bsxfun(@minus,all_im_patches,nanmean(all_im_patches));
all_im_patches = bsxfun(@rdivide,all_im_patches,nanstd(all_im_patches));

[XXcoords,YYcoords] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

RF_subset = find(XXcoords >= RF_patch(1,1) & XXcoords <= RF_patch(1,2) & ...
    YYcoords >= RF_patch(2,1) & YYcoords <= RF_patch(2,2));

SDIM = length(RF_xpatch_inds);
full_SDIM = length(xpatch_inds);

[Yinds,Tinds,Xinds] = meshgrid(yax(ypatch_inds),1:flen,xax(xpatch_inds));
% use_kInds = find(Xinds >= RF_patch(1,1) & Xinds <= RF_patch(1,2) & ...
%     Yinds >= RF_patch(2,1) & Yinds <= RF_patch(2,2));
use_kInds = find(ismember(Xinds,xax(RF_xpatch_inds)) & ismember(Yinds,yax(RF_ypatch_inds)));

full_stim_params = NMMcreate_stim_params([flen full_SDIM full_SDIM],dt);
full_X = create_time_embedding(all_im_patches,full_stim_params);

%%
% % clear sta sta2
% for ss = 1:96
%     ss
%     Robs = convert_to_spikebins(all_binned_mua(used_inds,ss));
% %     sta(ss,:) = mean(full_X(used_inds(Robs),:)) - mean(full_X(used_inds,:));
%     sta2(ss,:) = mean(abs(full_X(used_inds(Robs),:))) - mean(abs(full_X(used_inds,:)));
% end
% %%
% close all
% for ss = 1:96
%     ss
%     
% %     temp2 = reshape(sta2(ss,:),[flen SDIM SDIM]);
%     temp = reshape(sta2(ss,use_kInds),[flen SDIM SDIM]);
%     cm = max(abs(temp(:)))*0.85;
%     for ii = 1:flen
%     subplot(1,flen,ii)
%     imagesc(squeeze(temp(ii,:,:)));
%     caxis([-cm cm]);
%     end
%     
% %     subplot(2,1,2)
% %     imagesc(squeeze(temp2(3,:,:)));colorbar
%     pause
%     clf
% end

%% EVALUATE LLS, FIT SPK NLS, EVALUATE NULL_LLS AND THRESHOLD LL IMP VALUES TO GET USABLE UNITS
init_stim_params = NMMcreate_stim_params([flen SDIM SDIM],dt);
init_stim_params(2) = NMMcreate_stim_params([n_blocks 1],dt);
if use_sac_kerns
    init_stim_params(3) = NMMcreate_stim_params(n_sac_bins,dt);
    init_stim_params(4) = NMMcreate_stim_params(n_sac_bins,dt);
end

null_stim_params = init_stim_params(2:end);

n_squared_filts = 2;
mod_signs = ones(1,n_squared_filts+2);

cd(anal_dir);

mod_data_name = 'Twod_eyecorr_mods5';
anal_name = 'Twod_eyecorr_noCprior5';

fprintf('Loading pre-computed initial models\n');
load(mod_data_name);
load(anal_name,'it_mods*','tr_set');
true_mods = it_mods_spkNL{end};
n_tr_chs = length(tr_set);

%% SELECT USABLE UNITS AND make Robs_mat

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
max_shift = 12;
max_shift = 12;
dshift = 1;
x_shifts = -max_shift:dshift:max_shift;
y_shifts = -max_shift:dshift:max_shift;

[Xsh,Ysh] = meshgrid(x_shifts,y_shifts);
[Xsh_inds,Ysh_inds] = meshgrid(1:length(x_shifts),1:length(y_shifts));
SH = [Xsh(:) Ysh(:)];
SH_inds = [Xsh_inds(:) Ysh_inds(:)];
n_shifts = size(SH,1);
zero_frame = find(SH(:,1) == 0 & SH(:,2) == 0);

fix_prior = -sum((SH*sp_dx).^2,2)/(2*fix_prior_sigma^2);
fix_prior = bsxfun(@minus,fix_prior,logsumexp(fix_prior)); %normalize

Ix = speye(full_SDIM);
Iy = speye(full_SDIM);
yshift_mat = cell(length(x_shifts),1);
for xx = 1:length(x_shifts)
    tempx = spdiags( ones(full_SDIM,1), y_shifts(xx), full_SDIM, full_SDIM);
    yshift_mat{xx} = kron(Iy, tempx);
end
xshift_mat = cell(length(y_shifts),1);
for yy = 1:length(x_shifts)
    tempy = spdiags( ones(full_SDIM,1), x_shifts(yy), full_SDIM, full_SDIM);
    xshift_mat{yy} = kron(tempy, Ix);
end
shift_mat = cell(n_shifts,1);
Tshift_mat = cell(n_shifts,1);
for xx = 1:n_shifts
    temp = xshift_mat{SH_inds(xx,1)}*yshift_mat{SH_inds(xx,2)};
    shift_mat{xx} = temp(:,use_kInds);
    Tshift_mat{xx} = temp;
end

%%
sim_range = 2*max_shift;

rand_pos = rand(n_fixs,2)*sim_range - max_shift;
rand_pos = round(rand_pos);

all_shift_stimmat = all_im_patches;
for ii=1:n_fixs
    cur_inds = used_inds(fix_ids == ii);
    temp = shift_matrix_Nd(all_im_patches(cur_inds,:,:),-rand_pos(ii,1),3);
    temp = shift_matrix_Nd(temp,-rand_pos(ii,2),2);
    all_shift_stimmat(cur_inds,:,:) = temp;
end
true_X = create_time_embedding(all_shift_stimmat,full_stim_params);

%%
cur_X{1} = true_X(used_inds,use_kInds);
cur_X{2} = Xblock(used_inds,:);
cur_X{3} = Xsac;
cur_X{4} = Xmsac;

pred_rate = nan(length(tr_set),length(used_inds));
for ii = 1:length(tr_set)
ii
[LL, penLL, pred_rate(ii,:)] = NMMmodel_eval(true_mods(tr_set(ii)),Robs_mat(:,ii),cur_X);
end

sim_Robs_mat = poissrnd(pred_rate)';

%%

su_inds = find(all_mod_SU(tr_set) > 0);
klen = flen*SDIM^2;

%% ITERATE FIXATION-BASED CORRECTIONS

it_mods{1} = true_mods;
it_mods_spkNL{1} = true_mods;
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
        
        frame_LLs(:,xx) = squeeze(nansum(sim_Robs_mat.*log(pred_rate) - pred_rate,2));
    end
    
    %% INFER MICRO-SAC SEQUENCE
    fix_LLs = nan(n_fixs,n_shifts);
    temp_rate = nan(n_fixs,length(tr_set));
    temp_nspk = nan(n_fixs,length(tr_set));
    for ii = 1:n_fixs
        cur_inds = pfix_start_inds(ii):pfix_stop_inds(ii);
        fix_LLs(ii,:) = sum(frame_LLs(cur_inds,:));
        temp_rate(ii,:) = mean(Robs_mat(cur_inds,:));
        temp_nspk(ii,:) = sum(Robs_mat(cur_inds,:));
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
            
            cdist = pdist2(SH*sp_dx,[SH(:,1)*sp_dx + measured_fix_deltas(t,1) SH(:,2)*sp_dx + measured_fix_deltas(t,2)]);
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

            cdist = pdist2(SH*sp_dx,[SH(:,1)*sp_dx + measured_fix_deltas(t+1,1) SH(:,2)*sp_dx + measured_fix_deltas(t+1,2)]);
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
