% clear all
% close all

addpath('~/other_code/fastBSpline/');

global Expt_name bar_ori monk_name rec_type rec_number

% Expt_name = 'M011';
% monk_name = 'jbe';
% bar_ori = 160; %bar orientation to use (only for UA recs)
% rec_number = 1;

% [266-80 270-60 275-135 277-70 281-140 287-90 289-160 294-40 296-45 297-0/90 5-50 9-0 10-60 11-160 12-0 13-100 14-40 320-100]

sname = 'rpt_variability_compact_FIN6_noxc';

% et_mod_data_name = 'full_eyetrack_initmods_Rinit';
% et_anal_name = 'full_eyetrack_Rinit';
et_mod_data_name = 'full_eyetrack_initmods_FIN2_Rinit';
et_anal_name = 'full_eyetrack_FIN2_Rinit';
mod_name = 'corrected_models_comp_FIN2';

use_MUA = false; EP_params.use_MUA = use_MUA; %use MUA in model-fitting
use_hres_ET = true; EP_params.use_hres_ET = use_hres_ET; %use high-res eye-tracking?
exclude_sacs = false; EP_params.exclude_sacs = exclude_sacs; %exclude data surrounding microsaccades?
sub_trialavgs = false; EP_params.sub_trialavgs = sub_trialavgs; %subtract out trial avg spike counts?
do_xcorrs = false; EP_params.do_xcorrs = do_xcorrs; %compute pairwise stats
compute_sims = false; EP_params.compute_sims = compute_sims; %do simulated calcs for alphas
compute_PF_rate = false;

poss_bin_dts = [0.005 0.01 0.02 0.04 0.08 0.16 0.32]; EP_params.poss_bin_dts = poss_bin_dts; %possible time bins to test
direct_bin_dts = [0.005 0.01 0.02 0.04 0.08 0.16 0.32]; EP_params.direct_bin_dts = direct_bin_dts; %time bins to use for direct estimates
mod_bin_dts = [0.005 0.01 0.02 0.04 0.08 0.16 0.32]; EP_params.mod_bin_dts = mod_bin_dts; %possible time bins for model-based analysis
% poss_bin_dts = [0.01 0.02 0.05 0.1]; EP_params.poss_bin_dts = poss_bin_dts; %possible time bins to test
% direct_bin_dts = [0.01 0.02 0.05 0.1]; EP_params.direct_bin_dts = direct_bin_dts; %time bins to use for direct estimates
% mod_bin_dts = [0.01 0.02 0.05 0.1]; EP_params.mod_bin_dts = mod_bin_dts; %possible time bins for model-based analysis

max_tlag = 10; EP_params.max_tlag = max_tlag; %max time lag for computing xcorrs (units of dt bins)

sim_n_rpts = 500; EP_params.sim_n_rpts = sim_n_rpts; %number of repeats for simulation calcs

maxD_prc = 100; %maximum delta_X percentile to model with spline fit
n_EP_bins = 100; EP_params.n_EP_bins = n_EP_bins; %number of quantiles of delta_X for binned estimates
n_spline_knots = 4; EP_params.n_spline_knots = 4;  %number of spline knot pts
poss_eps_sizes = [.005 .01 .02 .04]; EP_params.poss_eps_sizes = poss_eps_sizes;   %possible epsilon balls to test
n_eval_pts = 100; EP_params.n_eval_pts = n_eval_pts; %number of points to evaluate spline fit
ball_nboots = 100; EP_params.ball_nboots = ball_nboots; %number of bootstrap samples for computing eps ball vars

use_LOOXV = 1; %[0 is no LOO; 1 is SUs only; 2 is SU + MU]

%%
data_dir = ['~/Data/bruce/' Expt_name];
if ~exist(data_dir,'dir')
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
end
if ~exist(data_dir,'dir');
    error('Couldnt find data directory');
end
anal_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
if ~exist(anal_dir,'dir')
    mkdir(anal_dir)
end

Edata_file = strcat(data_dir,'/',monk_name,Expt_name,'Expts');
load(Edata_file);

%is this a laminar probe or utah array rec?
ff = find(cellfun(@(x) ~isempty(x),Expts),1);
if strcmp(Expts{ff}.Header.DataType,'GridData 96')
    rec_type = 'UA';
    n_probes = 96;
elseif strcmp(Expts{ff}.Header.DataType,'Spike2')
    rec_type = 'LP';
    n_probes = 24;
end

%load in packaged data
data_name = sprintf('%s/packaged_data_ori%d',data_dir,bar_ori);
if rec_number > 1
    data_name = strcat(data_name,sprintf('_r%d',rec_number));
end
fprintf('Loading %s\n',data_name);
load(data_name);

ov_RF_pos = Expts{expt_data.used_blocks(1)}.Stimvals.rf(1:2)/params.scale_fac;

if sub_trialavgs
    sname = strcat(sname,'_subTrial');
end

%directories
et_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
model_dir = ['~/Analysis/bruce/' Expt_name '/models'];

if rec_number > 1
    cluster_dir = [cluster_dir sprintf('/rec%d',rec_number)];
end

%if using coil info
if any(params.use_coils > 0)
    et_anal_name = [et_anal_name '_Cprior'];
end
et_hres_anal_name = strcat(et_anal_name,'_hres');

et_mod_data_name = [et_mod_data_name sprintf('_ori%d',bar_ori)];
et_anal_name = [et_anal_name sprintf('_ori%d',bar_ori)];
et_hres_anal_name = [et_hres_anal_name sprintf('_ori%d',bar_ori)];
mod_name = [mod_name sprintf('_ori%d',bar_ori)];

if rec_number > 1
    mod_name = strcat(mod_name,sprintf('_r%d',rec_number));
    et_mod_data_name = strcat(et_mod_data_name,sprintf('r%d',rec_number));
    et_hres_anal_name = strcat(et_hres_anal_name,sprintf('r%d',rec_number));
    et_anal_name = strcat(et_anal_name,sprintf('r%d',rec_number));
end

% et_hres_anal_name = strcat(et_hres_anal_name,'_fullLOO');

%% LOAD EYE-TRACKING DATA
cd(et_dir)
load(et_mod_data_name,'all_mod*');
if use_hres_ET
    fprintf('Loading ET data %s\n',et_hres_anal_name);
    load(et_hres_anal_name)
else
    fprintf('Loading ET data %s\n',et_anal_name);
    load(et_anal_name);
end
tr_set = et_tr_set; %set of units used in ET

%for some recs the saccades for ET_data were detected using a slightly
%different algo. Eliminate the saccades causing the difference
if length(ET_data.saccades) ~= length(et_saccades)
    sac_start_times = [ET_data.saccades(:).start_time];
    old_sac_start_times = [et_saccades(:).start_time];
    if length(sac_start_times) > length(old_sac_start_times)
        extra_sacs = find(~ismember(sac_start_times,old_sac_start_times));
        ET_data.saccades(extra_sacs) = [];
        ET_data.is_blink(extra_sacs) = [];
        fprintf('Difference in saccade detection from eye-tracking data, eliminating %d/%d saccades\n',length(extra_sacs),length(sac_start_times));
    else
        error('Fewer saccades than in ETdata');
    end
end

%% LOAD IN MODEL FITS
load([model_dir '/' mod_name]);

%% BIN SPIKES FOR MU AND SU
all_binned_mua = spikes_int82double(spike_data.binned_mua);
all_binned_sua = spikes_int82double(spike_data.binned_sua);
Clust_data = spike_data.Clust_data;
su_probes = Clust_data.SU_probes;
SU_numbers = Clust_data.SU_numbers;

n_units = size(all_binned_mua,2) + size(all_binned_sua,2);

%% extract time series to keep track of trial number and block number
orig_dt = params.dt;

NT = length(used_inds);
fullNT = size(all_binned_mua,1);
n_trials = length(time_data.trial_flip_ids);
n_blocks = length(time_data.block_flip_ids);

all_t_axis = time_data.t_axis;
trial_start_inds = [1; 1+time_data.trial_flip_inds(2:end)];
trial_end_inds = [time_data.trial_flip_inds(2:end); fullNT];
all_trialvec = nan(fullNT,1); %trial index vector
for ii = 1:n_trials
    all_trialvec(trial_start_inds(ii):trial_end_inds(ii)) = time_data.trial_flip_ids(ii);
end

block_start_inds = [1; 1+time_data.block_flip_inds(2:end)];
block_end_inds = [time_data.block_flip_inds(2:end); fullNT];
all_blockvec = nan(fullNT,1); %block index vector
for ii = 1:n_blocks
    all_blockvec(block_start_inds(ii):block_end_inds(ii)) = time_data.block_flip_ids(ii);
end
Xblock = zeros(fullNT,n_blocks);
for i = 1:n_blocks
    cur_set = find(all_blockvec==i);
    Xblock(cur_set,i) = 1;
end

%% EXTRACT RELEVANT SACCADE DATA
corrected_eye_vals_interp = ET_data.interp_eye_pos;
sac_start_times = [ET_data.saccades(:).start_time];
sac_stop_times = [ET_data.saccades(:).stop_time];

interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_stop_times));

used_saccade_set = find(ismember(interp_sac_start_inds,used_inds)); %saccades occuring during used data

saccade_start_inds = find(ismember(used_inds,interp_sac_start_inds));

%nearest index in the used data set of the saccade stop time
saccade_stop_inds = round(interp1(used_inds,1:length(used_inds),interp_sac_stop_inds(used_saccade_set)))';
saccade_stop_inds(isnan(saccade_stop_inds)) = length(used_inds);

saccades = ET_data.saccades(used_saccade_set);
saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));

used_is_blink = ET_data.is_blink(used_saccade_set); %which used events are blinks

sac_durs = [saccades(:).duration];
sac_prepos = reshape([saccades(:).pre_pos],[],length(saccades));
sac_postpos = reshape([saccades(:).post_pos],[],length(saccades));
sac_deltaX = sac_postpos(1,:) - sac_prepos(1,:);

sacburst_set = find([saccades(:).isi] < params.sac_burst_isi | [saccades(:).next_isi] < params.sac_burst_isi);
micro_sacs = find([saccades(:).amplitude] < params.micro_thresh & ~used_is_blink');

msac_bursts = micro_sacs(ismember(micro_sacs,sacburst_set));
% micro_sacs(ismember(micro_sacs,sacburst_set)) = []; %eliminate microsacs that are part of a 'burst'

%guided saccades are those whose parallel component is large enough and
%that aren't blinks (and whose duration is not too long to be suspicious
big_sacs = find(abs(sac_deltaX) > params.gsac_thresh & ~used_is_blink' & sac_durs <= params.max_gsac_dur);
used_is_blink(sac_durs >= params.max_gsac_dur) = true;

%% DEFINE FIXATIONS
trial_start_inds = [1; find(diff(all_trialvec(used_inds)) ~= 0) + 1];
trial_end_inds = [find(diff(all_trialvec(used_inds)) ~= 0); NT];

fix_start_inds = sort([trial_start_inds; saccade_stop_inds]);
fix_stop_inds = sort([trial_end_inds; saccade_start_inds]);
fix_durs = fix_stop_inds-fix_start_inds;
fix_start_inds(fix_durs<=0) = [];
fix_stop_inds(fix_durs <= 0) = [];
fix_post_blink = ismember(fix_start_inds,saccade_stop_inds(used_is_blink)); %is this fixation following a blink
n_fixs = length(fix_start_inds);

%index values numbering the different fixations
fix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
end

%% Recon retinal stim for non LOO data
%set of units where we have LOOXV on eye-tracking
if use_LOOXV == 2
    loo_set = tr_set;
elseif use_LOOXV == 1
    loo_set = tr_set(all_mod_SU(tr_set) > 0);
else
    loo_set = [];
end

%for model fitting
if use_MUA
    targs = 1:n_units; %SU and MU
else
    targs = setdiff(1:n_units,1:params.n_probes); %SU only
end
targs(targs > length(ModData)) = [];

sac_shift = et_params.sac_shift; %forward projection of saccade start times
if use_hres_ET %if using high-res ET
    [post_mean_EP,post_std_EP] = construct_eye_position(best_fix_cor,best_fix_std,...
        drift_post_mean,drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);
    sp_dx = et_params.sp_dx;
    post_mean_EP = post_mean_EP*sp_dx; %in deg
    
    %now get separate eye-tracking sequences for each LOO set
    if exist('drift_post_mean_LOO','var')
        post_mean_EP_LOO = nan(length(loo_set),NT);
        for ss = 1:length(loo_set)
            [post_mean_EP_LOO(ss,:)] = construct_eye_position(best_fix_cor,best_fix_std,...
                drift_post_mean_LOO(ss,:),drift_post_std_LOO(ss,:),fix_ids,trial_start_inds,trial_end_inds,sac_shift);
        end
        post_mean_EP_LOO = post_mean_EP_LOO*sp_dx;
    else
        error('No LOO variables detected');
    end
else %if using base resolution ET
    [post_mean_EP,post_std_EP] = construct_eye_position(it_fix_post_mean(end,:),it_fix_post_std(end,:),...
        drift_post_mean(end,:),drift_post_std(end,:),fix_ids,trial_start_inds,trial_end_inds,sac_shift);
    sp_dx = et_params.sp_dx;
    post_mean_EP = post_mean_EP*sp_dx;
    
    post_mean_EP_LOO = nan(length(loo_set),NT);
    for ss = 1:length(loo_set)
        [post_mean_EP_LOO(ss,:)] = construct_eye_position(squeeze(it_fix_post_mean_LOO(ss,end,:)),squeeze(it_fix_post_std_LOO(ss,end,:)),...
            squeeze(drift_post_mean_LOO(ss,end,:)),squeeze(drift_post_std_LOO(ss,end,:)),fix_ids,trial_start_inds,trial_end_inds,sac_shift);
    end
    post_mean_EP_LOO = post_mean_EP_LOO*sp_dx;
end

%% IDENTIFY REPEAT TRIALS
all_trial_Se = [trial_data(:).se];
n_rpt_seeds = length(params.rpt_seeds); EP_params.n_rpt_seeds = n_rpt_seeds; %seed values for repeat sequences

all_rpt_trials = find(ismember(all_trial_Se,params.rpt_seeds)); %trial indices of repeat trials
all_nrpt_trials = find(~ismember(all_trial_Se,params.rpt_seeds)); %trial indices of non-repeat trials
rptframe_trials = find([trial_data(all_rpt_trials).nrpt_frames] > 0); %identify repeat trials where there were repeated frames
fprintf('Detected %d/%d trials with rpt frames\n',length(rptframe_trials),length(all_rpt_trials));

%get number of frames for each repeat trial
T = tabulate(all_trialvec);
rpt_tdurs = T(all_rpt_trials,2);

too_short = find(rpt_tdurs < 390); %dont use repeat trials that dont have ~400 frames
if ~isempty(too_short)
    fprintf('Eliminating %d/%d repeat trials without enough frames\n',length(too_short),length(rpt_tdurs));
    all_rpt_trials(too_short) = [];
    rptframe_trials = find([trial_data(all_rpt_trials).nrpt_frames] > 0); %recompute which repeat trials have repeat frames
end

if any(arrayfun(@(x) any(x.rpt_frames == 0),trial_data(all_rpt_trials))) %if there are rpt frames listed as frame 0 (these are special cases)
    %find which repeat trials have a mix of 0 and non-0 rpt_frames
    mixed_trials = false(length(all_rpt_trials),1);
    for ii = 1:length(all_rpt_trials)
        if any(trial_data(all_rpt_trials(ii)).rpt_frames == 0) & any(trial_data(all_rpt_trials(ii)).rpt_frames > 0)
            mixed_trials(ii) = true;
        end
    end
    %get rid of any mixed trials (these are too difficult to handle)
    if any(mixed_trials)
        fprintf('Eliminating %d/%d repeat trials with mixed rpt frame types\n',sum(mixed_trials),length(mixed_trials));
        all_rpt_trials(mixed_trials) = [];
        rptframe_trials = find([trial_data(all_rpt_trials).nrpt_frames] > 0); %recompute which rpt trials have rpt frames
    end
end

all_stim_mat = decompressTernNoise(stimComp);

%find any trials where the stimulus isn't displayed
temp_rpt_inds = used_inds(ismember(all_trialvec(used_inds),all_rpt_trials));
test_mat = reshape(all_stim_mat(temp_rpt_inds,:),[],length(all_rpt_trials),size(all_stim_mat,2));
bad_stim_trials = find(all(all(test_mat == 0),3)); %these trials have all 0 values (only once have i seen this, not sure what happened)
if ~isempty(bad_stim_trials)
    fprintf('Eliminating %d/%d bad stim trials\n',length(bad_stim_trials),length(all_rpt_trials));
    all_rpt_trials(bad_stim_trials) = [];
    rptframe_trials = find([trial_data(all_rpt_trials).nrpt_frames] > 0); %recompute which rpt trials have rpt frames
end

%keep track of which rpt sequence is displayed in each repeat trial
all_rpt_seqnum = nan(size(all_rpt_trials));
for ii = 1:n_rpt_seeds
    cur_trials = find(all_trial_Se(all_rpt_trials) == params.rpt_seeds(ii));
    all_rpt_seqnum(cur_trials) = ii;
end

tot_nrpts = length(all_rpt_seqnum); %number of used repeat trials
fprintf('Using %d repeat trials, %d sequences\n',tot_nrpts,length(params.rpt_seeds));

rpt_trial_block = [trial_data(all_rpt_trials).block_nums]; %which block was each repeat trial in

%% IDENTIFY TIMES WITHIN SACCADES AND BLINKS
sac_buff = round(0.05/params.dt); EP_params.sac_buff = sac_buff; %window of data to exclude during saccades
sac_delay = round(0.03/params.dt); EP_params.sac_delay = sac_delay; %shift exclusion window to account for neural delay
blink_buff = round(0.1/params.dt); EP_params.blink_buff = blink_buff; %window for excluding blinks

%these are the non-blink saccades (real saccades)
nblink_start_inds = saccade_start_inds(~used_is_blink);
nblink_stop_inds = saccade_stop_inds(~used_is_blink);
in_sac_inds = false(NT,1); %store indices that count as 'affected' by a saccade
for ii = 1:length(nblink_start_inds)
    cur_inds = (nblink_start_inds(ii):(nblink_stop_inds(ii) + sac_buff)) + sac_delay;
    cur_inds(cur_inds > NT) = [];
    in_sac_inds(cur_inds) = true;
end

%find indices that count as affected by a blink
blink_start_inds = saccade_start_inds(used_is_blink);
blink_stop_inds = saccade_stop_inds(used_is_blink);
in_blink_inds = false(NT,1);
for ii = 1:length(blink_start_inds)
    cur_inds = (blink_start_inds(ii):(blink_stop_inds(ii) + blink_buff)) + sac_delay;
    cur_inds(cur_inds > NT) = [];
    in_blink_inds(cur_inds) = true;
end

used_rpt_inds = find(ismember(all_trialvec(used_inds),all_rpt_trials)); %indices of repeat trials within used_inds vector

nf = 400; %number of frames per trial
used_nf = nf-(params.beg_buffer + params.end_buffer)/params.dt; %number of frames per trial used in analysis

%reshape these vectors into trial-by-trial matrices
in_sac_inds = reshape(in_sac_inds(used_rpt_inds),used_nf,tot_nrpts);
in_blink_inds = reshape(in_blink_inds(used_rpt_inds),used_nf,tot_nrpts);

%% PROCESS MODEL FITS

has_stim_mod = false(length(targs),1);
for cc = 1:length(targs) %loop over units used in analysis
    if ~isempty(ModData(targs(cc)).bestGQM) %if we have a model for this unit
        cur_mod = ModData(targs(cc)).bestGQM; %use the optimized GQM
        stim_mod_withblock(cc) = cur_mod;
                
        %remove the block-by-block variability in the model by
        %incorporating the avg output of the block-filter
        cur_block_filt = cur_mod.mods(1).filtK; %filter applied to the block index
        cur_used_blocks = ModData(targs(cc)).unit_data.used_blocks; %which blocks was this neuron isolated during
        poss_used_blocks = ModData(targs(cc)).unit_data.poss_used_blocks; %total set of used blocks
        cur_used_blocks = find(ismember(poss_used_blocks,cur_used_blocks)); %indices of blocks where this neuron was isolated
        cur_mod.spk_NL_params(1) = cur_mod.spk_NL_params(1) + mean(cur_block_filt(cur_used_blocks)); %add the average output of the block-filter to the spkNL offset
        cur_mod.mods(1) = []; %eliminate block filter
        
        %store model data
        stim_mod(cc) = cur_mod;
        EP_data(cc,1).unit_data = ModData(targs(cc)).unit_data;
        EP_data(cc,1).tune_props = ModData(targs(cc)).tune_props;
        EP_data(cc,1).bestGQM = ModData(targs(cc)).bestGQM;
        EP_data(cc,1).nullMod = ModData(targs(cc)).nullMod;
        EP_data(cc,1).modFitParams = modFitParams;
        EP_data(cc,1).useMod = cur_mod;
        has_stim_mod(cc) = true;
    end
end

%% determine which rpt indices need to be blanked out because of repeat frames
used_frame_inds = (params.beg_buffer/params.dt+1):(params.trial_dur - params.end_buffer)/params.dt;

rpt_blanked = false(nf,tot_nrpts);
N_wrong_start_trials = 0;
N_rptframe_trials = 0;
if ~isempty(rptframe_trials)
    post_rpt_buffer = round(0.1/params.dt); %exclude data for this duration following each rpt frame
    
    for ii = 1:length(rptframe_trials)
        cur_trial = all_rpt_trials(rptframe_trials(ii));
        cur_rpt_frames = trial_data(cur_trial).rpt_frames;
        if all(cur_rpt_frames == 0) %for all zero frames the entire stimulus was shifted by an amount = num_zero_frames
            shift_amount = length(cur_rpt_frames);
            to_blank_inds = false(nf,1);
            to_blank_inds(1:(shift_amount + params.beg_buffer/params.dt)) = true;
            rpt_blanked(to_blank_inds,rptframe_trials(ii)) = true;
            N_wrong_start_trials = N_wrong_start_trials + 1;
        elseif ~any(cur_rpt_frames == 0)
            new_frame_ids = 1:nf;
            for jj = 1:length(cur_rpt_frames)
                target_inds = (cur_rpt_frames(jj) + 1):nf;
                map_to = target_inds + 1; map_to(map_to > nf) = nf;
                new_frame_ids(target_inds) = new_frame_ids(map_to);
                
                rpt_blanked(cur_rpt_frames(jj):(cur_rpt_frames(jj)+post_rpt_buffer),rptframe_trials(ii)) = true;
                rpt_blanked((nf-length(cur_rpt_frames)*1):nf,rptframe_trials(ii)) = true;
            end
            N_rptframe_trials = N_rptframe_trials + 1;
        else
            error('mix of 0 and nonzero rpt frames');
        end
        
    end
end
usedrpt_blanked = rpt_blanked(used_frame_inds,:);

%% get stimulus xmat during repeat trials

%spatial up-sampling of the stimulus
full_nPix_us = modFitParams.spatial_usfac*params.full_nPix;
if modFitParams.spatial_usfac > 1
    all_stimmat_up = zeros(size(all_stim_mat,1),full_nPix_us);
    for ii = 1:size(all_stim_mat,2)
        for jj = 1:modFitParams.spatial_usfac
            all_stimmat_up(:,modFitParams.spatial_usfac*(ii-1)+jj) = all_stim_mat(:,ii);
        end
    end
elseif modFitParams.spatial_usfac == 1
    all_stimmat_up = all_stim_mat;
end

%select subset of pixels used for model fitting
buffer_pix = floor((params.full_nPix - modFitParams.use_nPix_us/modFitParams.spatial_usfac)/2);
[Xinds_up,~] = meshgrid(1/modFitParams.spatial_usfac:1/modFitParams.spatial_usfac:params.full_nPix,1:modFitParams.flen);
cur_use_pix = (1/modFitParams.spatial_usfac:1/modFitParams.spatial_usfac:(modFitParams.use_nPix_us/modFitParams.spatial_usfac)) + buffer_pix;
use_kInds_up = find(ismember(Xinds_up(:),cur_use_pix));

stim_params_full = NMMcreate_stim_params([modFitParams.flen full_nPix_us]);

%indices of repeat trials
full_rpt_inds = find(ismember(all_trialvec,all_rpt_trials));

if compute_PF_rate
    all_Xmat_PF = create_time_embedding(all_stimmat_up(full_rpt_inds,:),stim_params_full); %incorporate time embedding for whole-trials on repeats
    all_Xmat_PF = all_Xmat_PF(ismember(full_rpt_inds,used_inds),use_kInds_up); %take only usable repeat-trial data
    %if using tent-basis spatial-upsampling
    if modFitParams.add_usfac > 1
        all_Xmat_PF = tb_proc_stim(all_Xmat_PF,modFitParams.add_usfac,modFitParams.flen);
    end
end
%% get model-predicted trial-by-trial firing rates
% make Robs_mat
tot_sus = size(all_binned_sua,2);
tot_nUnits = length(su_probes) + n_probes;
Robs_mat = nan(length(used_inds),n_probes + tot_sus);
for ss = 1:size(Robs_mat,2)
    if ss > n_probes
        Robs_mat(:,ss) = all_binned_sua(used_inds,ss-n_probes);
    else
        Robs_mat(:,ss) = all_binned_mua(used_inds,ss);
    end
end

X{2} = Xblock(used_rpt_inds,:);

all_mod_emp_prates = nan(length(used_rpt_inds),length(targs));
for cc = 1:length(targs)
    if has_stim_mod(cc)
        
        loo_ind = find(loo_set == targs(cc));
        if ~isempty(loo_ind)
            cur_Robs = Robs_mat(used_rpt_inds,targs(cc));
            cur_uinds = find(~isnan(cur_Robs));
            cur_rpt_trialset = unique(all_trialvec(used_inds(used_rpt_inds(cur_uinds))));
            EP_data(cc,1).rpt_trialset = cur_rpt_trialset;
            
            cur_nullMod = ModData(targs(cc)).nullMod;
            
            post_mean_EP_rpt = post_mean_EP_LOO(loo_ind,used_rpt_inds);
            fin_shift_cor = round(post_mean_EP_rpt/modFitParams.sp_dx); %use overall EP estimate
            fin_shift_cor(isnan(fin_shift_cor)) = 0;
            
            %RECOMPUTE XMAT
            best_shift_stimmat_up = all_stimmat_up;
            for i=1:length(used_rpt_inds) %correct repeat trial data
                best_shift_stimmat_up(used_inds(used_rpt_inds(i)),:) = shift_matrix_Nd(all_stimmat_up(used_inds(used_rpt_inds(i)),:),-fin_shift_cor(i),2);
            end
            X{1} = create_time_embedding(best_shift_stimmat_up(full_rpt_inds,:),stim_params_full); %incorporate time embedding for whole-trials on repeats
            X{1} = X{1}(ismember(full_rpt_inds,used_inds),use_kInds_up); %take only usable repeat-trial data
            
            %if using tent-basis spatial-upsampling
            if modFitParams.add_usfac > 1
                X{1} = tb_proc_stim(X{1},modFitParams.add_usfac,modFitParams.flen);
            end
            
            rpt_LL = NMMeval_model(stim_mod_withblock(cc),cur_Robs,X,[],cur_uinds);
            rpt_nullLL = NMMeval_model(cur_nullMod,cur_Robs,X,[],cur_uinds);
            EP_data(cc,1).rpt_LL = rpt_LL;
            EP_data(cc,1).rpt_nullLL = rpt_nullLL;
            [~,~,all_mod_emp_prates(:,cc)] = NMMmodel_eval(stim_mod(cc),[],X);
        end
    end
end
all_mod_emp_prates_noEM = all_mod_emp_prates; %make a copy of the model-predicted rates that will be nan whenever theres a blink (or sac), as with the spk data
all_mod_emp_prates_noEM(in_blink_inds(:),:) = nan;
if exclude_sacs
    all_mod_emp_prates_noEM(in_sac_inds(:),:) = nan;
end
all_mod_emp_prates = reshape(all_mod_emp_prates,used_nf,tot_nrpts,length(targs)); %convert to tr-x-tr
all_mod_emp_prates_noEM = reshape(all_mod_emp_prates_noEM,used_nf,tot_nrpts,length(targs)); %convert to tr-x-tr

%if computing model-predicted rates for perfect-fixation case
if compute_PF_rate
    all_mod_PF_prates = nan(length(used_rpt_inds),length(targs));
    for cc = 1:length(targs)
        if has_stim_mod(cc)
            [~,~,all_mod_PF_prates(:,cc)] = NMMmodel_eval(stim_mod(cc),[],all_Xmat_PF);
        end
    end
    all_mod_PF_prates = reshape(all_mod_PF_prates,used_nf,tot_nrpts,length(targs));
end

%% shift model-predicted rates to realign repeats on trials with rpt frames
if ~isempty(rptframe_trials)
    shifted_mod_prates = nan(nf,tot_nrpts,length(targs)); %buffer with nans to trial-length
    shifted_mod_prates_noEM = nan(nf,tot_nrpts,length(targs));
    shifted_mod_prates(used_frame_inds,:,:) = all_mod_emp_prates;
    shifted_mod_prates_noEM(used_frame_inds,:,:) = all_mod_emp_prates_noEM;
    
    for ii = 1:length(rptframe_trials)
        cur_trial = all_rpt_trials(rptframe_trials(ii));
        cur_rpt_frames = trial_data(cur_trial).rpt_frames;
        new_spike_frame_ids = 1:nf;
        if all(cur_rpt_frames == 0)
            shift_amount = length(cur_rpt_frames);
            shifted_mod_prates(:,rptframe_trials(ii),:) = shift_matrix_Nd(shifted_mod_prates(:,rptframe_trials(ii),:),shift_amount,1);
            shifted_mod_prates_noEM(:,rptframe_trials(ii),:) = shift_matrix_Nd(shifted_mod_prates_noEM(:,rptframe_trials(ii),:),shift_amount,1);
        else
            for jj = 1:length(cur_rpt_frames)
                target_inds = (cur_rpt_frames(jj) + 1):nf;
                map_to = target_inds + 1; map_to(map_to > nf) = nf;
                new_spike_frame_ids(target_inds) = new_spike_frame_ids(map_to);
            end
            shifted_mod_prates(:,rptframe_trials(ii),:) = shifted_mod_prates(new_spike_frame_ids,rptframe_trials(ii),:);
            shifted_mod_prates_noEM(:,rptframe_trials(ii),:) = shifted_mod_prates_noEM(new_spike_frame_ids,rptframe_trials(ii),:);
        end
    end
    shifted_mod_prates = shifted_mod_prates(used_frame_inds,:,:);
    shifted_mod_prates = reshape(shifted_mod_prates,[],length(targs));
    shifted_mod_prates(usedrpt_blanked(:),:) = nan;
    shifted_mod_prates = reshape(shifted_mod_prates,length(used_frame_inds),tot_nrpts,length(targs));
    all_mod_emp_prates = shifted_mod_prates;
    
    shifted_mod_prates_noEM = shifted_mod_prates_noEM(used_frame_inds,:,:);
    shifted_mod_prates_noEM = reshape(shifted_mod_prates_noEM,[],length(targs));
    shifted_mod_prates_noEM(usedrpt_blanked(:),:) = nan;
    shifted_mod_prates_noEM = reshape(shifted_mod_prates_noEM,length(used_frame_inds),tot_nrpts,length(targs));
    all_mod_emp_prates_noEM = shifted_mod_prates_noEM;
end
if compute_PF_rate
    nrptframe_trials = setdiff(1:length(all_rpt_trials),rptframe_trials);
    all_mod_PF_prates = squeeze(all_mod_PF_prates(:,nrptframe_trials(1),:));
end

rpt_EP_SD = robust_std_dev(post_mean_EP_rpt); %SD of EP during repeat trials

rpt_data.tot_nrpts = tot_nrpts;
rpt_data.N_rptframe_trials = N_rptframe_trials;
rpt_data.N_wrong_start_trials = N_wrong_start_trials;
%% loop over possible time windows
for bbb = 1:length(poss_bin_dts)
    bin_dt = poss_bin_dts(bbb);
    fprintf('Running analysis at dt = %.3f sec\n',bin_dt);
    
    %% Construct Time embedded eye position sequences for repeat trials
    emb_win = (bin_dt + 0.05); EP_params.emb_win = emb_win; %look back this many time steps to parse EP trajectories
    emb_shift = 0.03; EP_params.emb_shift = emb_shift; %account for neural delay (relevant EP data is this far in the past relative to spiking data)
    emb_win = round(emb_win/orig_dt); %in units of original time bins
    emb_shift = round(emb_shift/orig_dt);
    
    tbt_EP = reshape(post_mean_EP_rpt,used_nf,tot_nrpts);
    
    %initialize a time-embedded version of the EPs
    sp = NMMcreate_stim_params(emb_win + emb_shift);
    tbt_EP_emb = create_time_embedding(tbt_EP(:),sp);
    tbt_EP_emb(in_blink_inds(:),:) = nan; %exclude in-blink data
    if exclude_sacs %if excluding in-sac data
        tbt_EP_emb(in_sac_inds(:),:) = nan;
    end
    tbt_EP_emb = reshape(tbt_EP_emb(:,(emb_shift+1):end),used_nf,tot_nrpts,[]); %drop the first emb_shift time lags, and reshape into a trial-by-trial array
    
    % now compile all time-embedded LOO EP sequences
    ms = size(tbt_EP_emb);
    if ~isempty(loo_set)
        loo_tbt_EP_emb = nan(ms(1),ms(2),ms(3),length(loo_set));
        for cc = 1:length(targs)
            loo_ind = find(loo_set == targs(cc));
            if ~isempty(loo_ind)
                post_mean_rpt = post_mean_EP_LOO(loo_ind,used_rpt_inds);
                cur_tbt_EP = reshape(post_mean_rpt,used_nf,tot_nrpts); %make trial-by-trial EP matrix
                
                %initialize a time-embedded version of the EPs
                cur_tbt_EP_emb = create_time_embedding(cur_tbt_EP(:),sp);
                cur_tbt_EP_emb(in_blink_inds(:),:) = nan;
                if exclude_sacs
                    cur_tbt_EP_emb(in_sac_inds(:),:) = nan;
                end
                loo_tbt_EP_emb(:,:,:,loo_ind) = reshape(cur_tbt_EP_emb(:,(emb_shift+1):end),used_nf,tot_nrpts,[]);
            end
        end
    end
    %% shift EP data to realign repeats on trials with rpt frames
    if ~isempty(rptframe_trials)
        post_rpt_buffer = round(0.1/params.dt); %exclude data for this duration following each rpt frame
        
        %add a buffer of nans to the beginning
        tbt_EP_emb = cat(1,nan(params.beg_buffer/params.dt,tot_nrpts,emb_win),tbt_EP_emb,nan(params.end_buffer/params.dt+1,tot_nrpts,emb_win));
        tbt_EP_emb((nf+1):end,:,:) = [];
        if ~isempty(loo_set)
            loo_tbt_EP_emb = cat(1,nan(params.beg_buffer/params.dt,tot_nrpts,emb_win,length(loo_set)),loo_tbt_EP_emb,nan(params.end_buffer/params.dt+1,tot_nrpts,emb_win,length(loo_set)));
            loo_tbt_EP_emb((nf+1):end,:,:,:) = [];
        end
        tbt_EP = cat(1,nan(params.beg_buffer/params.dt,tot_nrpts),tbt_EP,nan(params.end_buffer/params.dt+1,tot_nrpts));
        tbt_EP((nf+1):end,:) = [];
        for ii = 1:length(rptframe_trials) %loop over all rpt trials that have rpt frames
            cur_trial = all_rpt_trials(rptframe_trials(ii));
            cur_rpt_frames = trial_data(cur_trial).rpt_frames;
            if all(cur_rpt_frames == 0) %if it's a trial with 0 rpt frames then it starts on the Nth frame
                shift_amount = length(cur_rpt_frames); %number of 0s indicates which frame the trial really started on
                %shift the EP data forward in time this many frames
                tbt_EP_emb(:,rptframe_trials(ii),:) = shift_matrix_Nd(tbt_EP_emb(:,rptframe_trials(ii),:),shift_amount,1);
                if ~isempty(loo_set)
                    loo_tbt_EP_emb(:,rptframe_trials(ii),:,:) = shift_matrix_Nd(loo_tbt_EP_emb(:,rptframe_trials(ii),:,:),shift_amount,1);
                end
                tbt_EP(:,rptframe_trials(ii)) = shift_matrix_Nd(tbt_EP(:,rptframe_trials(ii)),shift_amount,1);
            elseif ~any(cur_rpt_frames == 0) %otherwise, this variable stores the index of frames that were repeated
                new_frame_ids = 1:nf; %index values to map onto
                for jj = 1:length(cur_rpt_frames)
                    target_inds = (cur_rpt_frames(jj) + 1):nf; %everything after the repeat frame
                    map_to = target_inds + 1; %needs to be shifted forward in time by 1
                    map_to(map_to > nf) = nf; %cap at nf
                    new_frame_ids(target_inds) = new_frame_ids(map_to);
                end
                tbt_EP_emb(:,rptframe_trials(ii),:) = tbt_EP_emb(new_frame_ids,rptframe_trials(ii),:);
                if ~isempty(loo_set)
                    loo_tbt_EP_emb(:,rptframe_trials(ii),:,:) = loo_tbt_EP_emb(new_frame_ids,rptframe_trials(ii),:,:);
                end
                tbt_EP(:,rptframe_trials(ii)) = tbt_EP(new_frame_ids,rptframe_trials(ii));
            else
                error('mix of 0 and nonzero rpt frames'); %shouldnt be any of these
            end
        end
        tbt_EP_emb = tbt_EP_emb(used_frame_inds,:,:);
        tbt_EP = tbt_EP(used_frame_inds,:);
        if ~isempty(loo_set); loo_tbt_EP_emb = loo_tbt_EP_emb(used_frame_inds,:,:,:); end
    end
    %convert to Tx1 arrays
    tbt_EP_emb = reshape(tbt_EP_emb,[],emb_win);
    %make sure the data that needs to be blanked out is
    tbt_EP_emb(usedrpt_blanked(:),:) = nan;
    %flip back to trxtr arrays
    tbt_EP_emb = reshape(tbt_EP_emb,used_nf,tot_nrpts,emb_win);
    
    %     tbt_EP = reshape(tbt_EP,[],1);
    %     tbt_EP(usedrpt_blanked(:),:) = nan;
    %     tbt_EP = reshape(tbt_EP,used_nf,tot_nrpts);
    
    if ~isempty(loo_set)
        loo_tbt_EP_emb = reshape(loo_tbt_EP_emb,[],emb_win,length(loo_set));
        loo_tbt_EP_emb(usedrpt_blanked(:),:,:) = nan;
        loo_tbt_EP_emb = reshape(loo_tbt_EP_emb,used_nf,tot_nrpts,emb_win,length(loo_set));
    end
    %% get binned spike data
    %if using a coarser time-binning, first compute spikes at native dt
    %resolution
    if bin_dt > params.dt
        cur_bin_dt = params.dt;
    else
        cur_bin_dt = bin_dt;
    end
    
    %centers of spike count time bins
    full_bin_taxis = (1:round(params.trial_dur)/cur_bin_dt)*cur_bin_dt-cur_bin_dt/2;
    bin_taxis = full_bin_taxis;
    bin_taxis(bin_taxis < params.beg_buffer) = [];
    bin_taxis(params.trial_dur - bin_taxis < params.end_buffer) = [];
    
    n_Tbins = length(full_bin_taxis); %total number of time bins
    n_used_Tbins = length(bin_taxis); %number of used time bins
    used_Tinds = find(ismember(full_bin_taxis,bin_taxis)); %index values of used time bins
    
    tbt_binned_spikes = nan(n_Tbins,tot_nrpts,length(SU_numbers));
    tbt_t_axis = nan(n_Tbins,tot_nrpts); %absolute times
    for ii = 1:tot_nrpts
        cur_bin_edges = [trial_data(all_rpt_trials(ii)).start_times:cur_bin_dt:(trial_data(all_rpt_trials(ii)).start_times + cur_bin_dt*(n_Tbins))];
        cur_bin_cents = 0.5*cur_bin_edges(1:end-1) + 0.5*cur_bin_edges(2:end); %bin centers
        for cc = 1:length(SU_numbers) %count spike for each unit in this trial
            if ~isnan(Clust_data.SU_block_probes(cc,rpt_trial_block(ii))) %was this unit isolated during this block?
                cur_hist = histc(spike_data.SU_spk_times{cc},cur_bin_edges);
                tbt_binned_spikes(:,ii,cc) = cur_hist(1:end-1);
            end
        end
        tbt_t_axis(:,ii) = cur_bin_cents;
    end
    
    %get the closest index of the original time axis for each of the new time
    %bins
    orig_t_ind = round(interp1(time_data.t_axis,1:fullNT,tbt_t_axis(:)));
    orig_t_ind(isnan(orig_t_ind)) = 1;
    orig_t_ind = reshape(orig_t_ind,n_Tbins,tot_nrpts);
    
    %determine which of the spkcnt time bins is in the set of original used
    %indices
    up_used_inds = false(size(orig_t_ind));
    for ii = 1:tot_nrpts
        if any(ismember(orig_t_ind(:,ii),used_inds))
            up_used_inds(used_Tinds,ii) = true;
        end
    end
    
    %for blocks where we didn't have an SU clustered, set the binned spkcnts to
    %nan (should already be nan, given the above checks) . leaving this in as
    %an additional check.
    tbt_binned_spikes = reshape(tbt_binned_spikes,[],length(SU_numbers));
    for ii = 1:length(SU_numbers)
        if any(~isnan(tbt_binned_spikes(isnan(all_binned_sua(orig_t_ind,ii)),ii)))
            error('This should already be nan');
        end
        %         tbt_binned_spikes(isnan(all_binned_sua(orig_t_ind,ii)),ii) = nan;
    end
    tbt_binned_spikes = reshape(tbt_binned_spikes,n_Tbins,tot_nrpts,length(SU_numbers));
    tbt_binned_spikes(:,:,length(targs)+1:end) = [];
    
    %% handle in-blink and in-sac data in the new spike count binning
    spk_in_blink_inds = false(n_Tbins,tot_nrpts);
    spk_in_sac_inds = false(n_Tbins,tot_nrpts);
    if bin_dt < params.dt %if we're binning at a finer time resolution
        %interpolate the indicator vectors onto the finer time axis
        temp_spk_in_blink_inds = round(interp1(all_t_axis(used_inds(used_rpt_inds)),double(in_blink_inds(:)),reshape(tbt_t_axis(used_Tinds,:),[],1)));
        temp_spk_in_blink_inds(isnan(temp_spk_in_blink_inds)) = 0; %dont worry about out-of-bounds values
        temp_spk_in_blink_inds = logical(reshape(temp_spk_in_blink_inds,length(used_Tinds),[]));
        temp_spk_in_sac_inds = round(interp1(all_t_axis(used_inds(used_rpt_inds)),double(in_sac_inds(:)),reshape(tbt_t_axis(used_Tinds,:),[],1)));
        temp_spk_in_sac_inds(isnan(temp_spk_in_sac_inds)) = 0;
        temp_spk_in_sac_inds = logical(reshape(temp_spk_in_sac_inds,length(used_Tinds),[]));
    else %otherwise these are unchanged
        temp_spk_in_blink_inds = in_blink_inds;
        temp_spk_in_sac_inds = in_sac_inds;
    end
    spk_in_blink_inds(used_Tinds,:) = temp_spk_in_blink_inds;
    spk_in_sac_inds(used_Tinds,:) = temp_spk_in_sac_inds;
    
    tbt_BS_ms = tbt_binned_spikes; %this will be the binned spike count array for analysis
    
    tbt_BS_ms = reshape(tbt_BS_ms,[],length(targs)); %flip to Tx1
    tbt_BS_ms(spk_in_blink_inds(:),:) = nan;
    if exclude_sacs
        tbt_BS_ms(spk_in_sac_inds(:),:) = nan;
    end
    tbt_BS_ms = reshape(tbt_BS_ms,n_Tbins,tot_nrpts,length(targs)); %flip back to trxtr
    
    %% shift spike data to align repeats on trials with rpt frames
    if ~isempty(rptframe_trials)
        dt_uf = params.dt/cur_bin_dt; %temporal up-sampling factor (floored at 1)
        if bin_dt < params.dt %for temporal up-sampling
            %interpolate rpt blank indices
            spk_rptblank_inds = round(interp1(all_t_axis(used_inds(used_rpt_inds)),double(usedrpt_blanked(:)),reshape(tbt_t_axis(used_Tinds,:),[],1)));
            spk_rptblank_inds(isnan(spk_rptblank_inds)) = 0; %dont worry about out of bounds
            spk_rptblank_inds = logical(spk_rptblank_inds);
        else
            spk_rptblank_inds = usedrpt_blanked;
        end
        
        for ii = 1:length(rptframe_trials)
            cur_trial = all_rpt_trials(rptframe_trials(ii));
            cur_rpt_frames = trial_data(cur_trial).rpt_frames;
            if all(cur_rpt_frames == 0) %zero values indicate which frame the trial really started on
                spk_shift_amount = length(cur_rpt_frames)*dt_uf; %number of up-sampled time steps to shift by
                tbt_BS_ms(:,rptframe_trials(ii),:) = shift_matrix_Nd(tbt_BS_ms(:,rptframe_trials(ii),:),spk_shift_amount,1);
                tbt_binned_spikes(:,rptframe_trials(ii),:) = shift_matrix_Nd(tbt_binned_spikes(:,rptframe_trials(ii),:),spk_shift_amount,1);
                
            elseif ~any(cur_rpt_frames == 0) %non-zero values indicate which frame was repeated
                new_spike_frame_ids = 1:nf*dt_uf;
                for jj = 1:length(cur_rpt_frames)
                    target_inds = (cur_rpt_frames(jj) + 1)*dt_uf:nf*dt_uf; %these time bins need to be shifted
                    map_to = target_inds + dt_uf; %amount of shift
                    map_to(map_to > nf*dt_uf) = nf*dt_uf; %cap
                    new_spike_frame_ids(target_inds) = new_spike_frame_ids(map_to);
                end
                tbt_BS_ms(:,rptframe_trials(ii),:) = tbt_BS_ms(new_spike_frame_ids,rptframe_trials(ii),:);
                tbt_binned_spikes(:,rptframe_trials(ii),:) = tbt_binned_spikes(new_spike_frame_ids,rptframe_trials(ii),:);
            else
                warning('mix of 0 and nonzero rpt frames');
            end
        end
        tbt_binned_spikes = reshape(tbt_binned_spikes(used_Tinds,:,:),[],length(targs)); %flip to Tx1
        tbt_binned_spikes(spk_rptblank_inds(:),:) = nan;%indices that need to be blanked because of repeat frames
        tbt_binned_spikes = reshape(tbt_binned_spikes,length(used_Tinds),tot_nrpts,length(targs)); %flip back to trxtr
        
        tbt_BS_ms = reshape(tbt_BS_ms(used_Tinds,:,:),[],length(targs)); %flip to Tx1
        tbt_BS_ms(spk_rptblank_inds(:),:) = nan;%indices that need to be blanked because of repeat frames
        tbt_BS_ms = reshape(tbt_BS_ms,length(used_Tinds),tot_nrpts,length(targs)); %flip back to trxtr
    else
        tbt_BS_ms = tbt_BS_ms(used_Tinds,:,:);
    end
    
    %% align relevant data to same time-binning
    if bin_dt > params.dt %if using coarser time binning
        bin_dsfac = bin_dt/params.dt; %temporal down-sampling factor
        if mod(bin_dsfac,1) ~= 0
            error('have to use integer multiple of dt for time binning');
        end
        
        n_Tbins = floor(used_nf/bin_dsfac); %new number of time bins
        
        %collect all spike counts within each of the coarser time bins
        new_BS_ms = nan(n_Tbins,tot_nrpts,length(targs),bin_dsfac);
        for ii = 1:bin_dsfac
            new_BS_ms(:,:,:,ii) = tbt_BS_ms(ii:bin_dsfac:(ii+bin_dsfac*(n_Tbins-1)),:,:,:);
        end
        new_BS_ms = squeeze(sum(new_BS_ms,4)); %now sum over the coarser time bins. This makes the whole bin a NAN if any component bins are NAN
        %         new_BS_ms = squeeze(nansum(new_BS_ms,4));
        
        %take the sub-sampled eye position history leading up to the last time bin
        ep_bin_ids = bin_dsfac:bin_dsfac:(bin_dsfac + bin_dsfac*(n_Tbins-1));
        new_EP_emb = tbt_EP_emb(ep_bin_ids,:,:);
        if ~isempty(loo_set); new_loo_EP_emb = loo_tbt_EP_emb(ep_bin_ids,:,:,:); end
        
    elseif bin_dt < params.dt %f using finer-than-dt binning
        bin_usfac = bin_dt/params.dt; %temporal up-sampling factor
        n_Tbins = length(used_Tinds);
        new_BS_ms = tbt_BS_ms;
        
        %use rounded indices to up-sample eye-position data
        ep_bin_ids = ceil(bin_usfac:bin_usfac:used_nf);
        new_EP_emb = tbt_EP_emb(ep_bin_ids,:,:);
        if ~isempty(loo_set); new_loo_EP_emb = loo_tbt_EP_emb(ep_bin_ids,:,:,:); end
        
    else %if using native binning, no changes needed
        n_Tbins = length(used_Tinds);
        new_BS_ms = tbt_BS_ms;
        new_EP_emb = tbt_EP_emb;
        if ~isempty(loo_set); new_loo_EP_emb = loo_tbt_EP_emb; end
    end
    
    if bin_dt > params.dt %if using coarser time binning
        tlags = 0; %just use zero-lag bin
    else
        tlags = -max_tlag:max_tlag; %range of time lags
    end
    
    %% BASIC STATS
    %first compute avg spike rates
    ov_avg_BS = nanmean(reshape(new_BS_ms,[],length(targs))); %overall avg rate
    trial_avg_BS = nanmean(new_BS_ms); %within-trial average rates
    for cc = 1:length(targs) %store data
        EP_data(cc,bbb).ov_avg_BS = ov_avg_BS(cc);
        EP_data(cc,bbb).trial_avg_BS = squeeze(trial_avg_BS(:,:,cc));
        EP_data(cc,1).EP_SD = rpt_EP_SD; %overall eye position SD
        EP_data(cc,bbb).tlags = tlags;
    end
    for rr = 1:n_rpt_seeds %for each repeat sequence
        cur_trial_set = find(all_rpt_seqnum == rr); %find the set of repeat trials
        cur_nrpts = length(cur_trial_set); %number of repeats
        
        n_utrials = squeeze(mean(sum(~isnan(new_BS_ms(:,cur_trial_set,:)),2))); %across-time avg of the number of used repeat trials for each unit
        n_spikes = squeeze(nansum(reshape(new_BS_ms,[],length(SU_numbers)))); %total number of spikes for this repeat sequence
        for cc = 1:length(targs) %store data
            EP_data(cc,bbb).n_utrials(rr) = n_utrials(cc);
            EP_data(cc,bbb).n_spikes(rr) = n_spikes(cc);
        end
    end
    
    %subtract out trial-avg spk counts if desired
    if sub_trialavgs
        new_BS_ms = bsxfun(@minus,new_BS_ms,trial_avg_BS);
    else
        %if not subtracting trial avgs, subtract out within-block avgs
        for ii = 1:length(expt_data.used_blocks)
            cur_block_trials = find(rpt_trial_block == ii); %find all repeat trials in this block
            if ~isempty(cur_block_trials)
                cur_block_avgs = nanmean(reshape(new_BS_ms(:,cur_block_trials,:),[],length(targs))); %within block avg rates
                new_BS_ms(:,cur_block_trials,:) = bsxfun(@minus,new_BS_ms(:,cur_block_trials,:),reshape(cur_block_avgs,1,1,length(targs)));
            end
        end
    end
    trial_avg_BS = squeeze(trial_avg_BS);
    
    %now subtract out overall avg spike count to ensure that this is
    %exactly zero
    new_BS_ms = bsxfun(@minus,new_BS_ms,reshape(nanmean(reshape(new_BS_ms,[],length(targs))),[1 1 length(targs)]));
    
    %% COMPUTE BASIC ACROSS-TRIAL STATS
    for rr = 1:n_rpt_seeds %separately for each repeat sequence
        cur_trial_set = find(all_rpt_seqnum == rr);
        cur_nrpts = length(cur_trial_set);
        
        psths = squeeze(nanmean(new_BS_ms(:,cur_trial_set,:),2)); %PSTHs for this sequence
        psth_var = nanvar(psths); %variance of PSTHs
        tot_resp_var = nanvar(reshape(new_BS_ms(:,cur_trial_set,:),[],length(SU_numbers))); %total spk cnt variance in these trials
        
        avg_temp_var = squeeze(nanmean(nanvar(new_BS_ms(:,cur_trial_set,:)))); %avg (across trials) of across-time variance
        psth_var_cor = psth_var.*(n_utrials'./(n_utrials'-1)) - avg_temp_var'./(n_utrials-1)'; %sahani linden correction for PSTH sampling noise
        
        avg_acrossTrial_var = squeeze(nanmean(nanvar(new_BS_ms(:,cur_trial_set,:),[],2))); %across-time avg of across-trial variances
        trial_avg_var = squeeze(nanvar(trial_avg_BS)); %variance of trial-avg rates
        
        for cc = 1:length(targs) %store data
            EP_data(cc,bbb).psths(rr,:) = psths(:,cc);
            EP_data(cc,bbb).psth_var(rr) = psth_var(cc);
            EP_data(cc,bbb).psth_var_cor(rr) = psth_var_cor(cc);
            EP_data(cc,bbb).across_trial_var(rr) = avg_acrossTrial_var(cc);
            EP_data(cc,bbb).tot_var(rr) = tot_resp_var(cc);
            EP_data(cc,bbb).trial_avg_var(rr) = trial_avg_var(cc);
        end
    end
    
    %% GET A SAMPLE OF THE BETWEEN TRIAL EYE POSITION SIMILARITY DISTRIBUTION AND DETERMINE ITS QUANTILES
    
    % compute the distribution of delta_X
    rand_deltaX = [];
    for rr = 1:n_rpt_seeds
        cur_trial_set = find(all_rpt_seqnum == rr);
        
        % estimate quantiles of the distribution of pairwise EP similarities by this metric
        rset = randi(length(used_frame_inds),100,1);
        for ii = 1:length(rset)
            cur_Dmat = abs(squareform(pdist(squeeze(tbt_EP_emb(rset(ii),cur_trial_set,:)))))/sqrt(emb_win);
            cur_Dmat(logical(eye(length(cur_trial_set)))) = nan;
            cur_Dmat = cur_Dmat(~isnan(cur_Dmat));
            rand_deltaX = cat(1,rand_deltaX,cur_Dmat);
        end
    end
    
    %bins for estimating variances as a function of delta_X (should be ~
    %equipopulated)
    EP_bin_edges = prctile(rand_deltaX,linspace(0,maxD_prc,n_EP_bins+1));
    EP_bin_centers = (EP_bin_edges(1:end-1)+EP_bin_edges(2:end))/2;
    maxD = prctile(rand_deltaX,maxD_prc);
    
    %% MAIN WITHIN-CELL ANALYSIS LOOP
    eval_xx = unique([0 prctile(rand_deltaX,linspace(0,maxD_prc,n_eval_pts))]); %x-axis for evaluating spline models
    
    if ismember(bin_dt,direct_bin_dts) %if we're doing direct calcs at this bin window
        
        for cc = 1:length(targs)
            fprintf('SU %d/%d\n',cc,length(targs));
            loo_ind = find(loo_set == targs(cc));
            if ~isempty(loo_ind)
                cur_tbt_EP_emb = squeeze(new_loo_EP_emb(:,:,:,loo_ind));
            end
            all_LOO_D = []; %LOO delta_X values
            all_base_D = []; %delta_X without LOO
            all_X = []; %Y_i*Y_j
            for rr = 1:n_rpt_seeds %loop over unique repeat seeds
                cur_trial_set = find(all_rpt_seqnum == rr);
                cur_nrpts = length(cur_trial_set);
                [II,JJ] = meshgrid(1:cur_nrpts);
                uset = JJ > II; %only need to count each unique trial pair once
                n_unique_pairs = sum(uset(:));
                
                cur_LOO_D = nan(n_unique_pairs*n_Tbins,1);
                cur_base_D = nan(n_unique_pairs*n_Tbins,1);
                cur_X = nan(n_unique_pairs*n_Tbins,1);
                for tt = 1:n_Tbins %loop over time bins
                    cur_inds = (tt-1)*n_unique_pairs + (1:n_unique_pairs); %range of index values for storage
                    Y1 = squeeze(new_BS_ms(tt,cur_trial_set,cc)); %responses at this time point
                    
                    %get matrix of delta_X using LOO EP signal
                    if ~isempty(loo_ind)
                        cur_Dmat = abs(squareform(pdist(squeeze(cur_tbt_EP_emb(tt,cur_trial_set,:)))))/sqrt(emb_win);
                        cur_Dmat(logical(eye(cur_nrpts))) = nan;
                        cur_LOO_D(cur_inds) = cur_Dmat(uset); %use only unique unequal trial pairs
                    end
                    
                    %now repeat using overall EP signal
                    cur_Dmat = abs(squareform(pdist(squeeze(new_EP_emb(tt,cur_trial_set,:)))))/sqrt(emb_win);
                    cur_Dmat(logical(eye(cur_nrpts))) = nan;
                    cur_base_D(cur_inds) = cur_Dmat(uset);
                    
                    %get matrix of products between responses
                    cur_Xmat = bsxfun(@times,Y1(II),Y1(JJ));
                    cur_X(cur_inds) = cur_Xmat(uset);
                end
                
                if ~isempty(loo_ind)
                    cur_upts = find(~isnan(cur_LOO_D) & ~isnan(cur_X)); %points we need to keep
                    all_LOO_D = cat(1,all_LOO_D,cur_LOO_D(cur_upts));
                else
                    cur_upts = find(~isnan(cur_base_D) & ~isnan(cur_X));
                end
                all_base_D = cat(1,all_base_D,cur_base_D(cur_upts));
                all_X = cat(1,all_X,cur_X(cur_upts));
            end
            n_data_points = length(all_X);
            
            if n_data_points > 0
                
                %compute histogram-based YY vs deltaX
                if ~isempty(loo_ind)
                    [bincnts,binids] = histc(all_LOO_D,EP_bin_edges);
                else
                    [bincnts,binids] = histc(all_base_D,EP_bin_edges);
                end
                var_ep_binned = nan(n_EP_bins,1);
                for bb = 1:n_EP_bins
                    var_ep_binned(bb) = mean(all_X(binids == bb));
                end
                
                %location of spline knot pts set by quantiles of deltaX dist
                if ~isempty(loo_ind)
                    spline_DS = prctile(all_LOO_D,maxD_prc/(n_spline_knots-1):maxD_prc/(n_spline_knots-1):(maxD_prc-maxD_prc/(n_spline_knots-1)));
                    knot_pts = [0 0 0 0 spline_DS maxD maxD maxD]; %handle knots at boundaries
                    upts = find(all_LOO_D <= knot_pts(end));
                else
                    spline_DS = prctile(all_base_D,maxD_prc/(n_spline_knots-1):maxD_prc/(n_spline_knots-1):(maxD_prc-maxD_prc/(n_spline_knots-1)));
                    knot_pts = [0 0 0 0 spline_DS maxD maxD maxD]; %handle knots at boundaries
                    upts = find(all_base_D <= knot_pts(end));
                end
                
                %compute <YY> for an epsilon-ball surrounding deltaX=0
                eps_ball_var = nan(length(poss_eps_sizes),1);
                eps_ball_var_noLOO = nan(length(poss_eps_sizes),1);
                eps_ball_npts = nan(length(poss_eps_sizes),1);
                eps_ball_boot = nan(length(poss_eps_sizes),2);
                bootfun = @(x) mean(all_X(x(randi(length(x),length(x),1)))); %bootstrap function which computes the eps-ball mean over a random-resampling of the subset of points within deltaE
                for bb = 1:length(poss_eps_sizes)
                    curset = find(all_base_D < poss_eps_sizes(bb));
                    eps_ball_var_noLOO(bb) = mean(all_X(curset));
                    if ~isempty(loo_ind)
                        curset = find(all_LOO_D < poss_eps_sizes(bb));
                        eps_ball_var(bb) = mean(all_X(curset));
                        eps_ball_npts(bb) = length(curset);
                        if eps_ball_npts(bb) > 1
                            bootsamps = bootstrp(ball_nboots,bootfun,curset);
                            eps_ball_boot(bb,:) = [mean(bootsamps),std(bootsamps)];
                        end
                    end
                end
                clear bootfun
                
                %fit cubic spline using LOO deltaX
                if ~isempty(loo_ind)
                    loo_sp = fastBSpline.lsqspline(knot_pts,3,all_LOO_D(upts),all_X(upts));
                end
                
                %fit cubic spline using overall deltaX
                upts = find(all_base_D <= knot_pts(end));
                base_sp = fastBSpline.lsqspline(knot_pts,3,all_base_D(upts),all_X(upts));
                
                %store data
                if ~isempty(loo_ind)
                    EP_data(cc,bbb).spline_looEP = loo_sp;
                end
                EP_data(cc,bbb).spline_baseEP = base_sp;
                EP_data(cc,bbb).eps_ball_var = eps_ball_var;
                EP_data(cc,bbb).eps_ball_var_noLOO = eps_ball_var_noLOO;
                EP_data(cc,bbb).eps_ball_npts = eps_ball_npts;
                EP_data(cc,bbb).eps_ball_boot = eps_ball_boot;
                EP_data(cc,bbb).var_ep_binned = var_ep_binned;
                EP_data(cc,bbb).pair_psth_var = mean(all_X); %estimate of PSTH variance as <Y_i*Y_j> independent of deltaX
                EP_data(cc,bbb).EP_bin_centers = EP_bin_centers;
                EP_data(cc,bbb).eval_xx = eval_xx;
            end
        end
        %         end
        
        %% COMPUTE FFs
        avg_rates = [EP_data(:,bbb).ov_avg_BS];
        %     ucells = find(~isnan(avg_rates));
        %     for cc = ucells
        for cc = 1:length(targs)
            if ~isnan(avg_rates(cc))
                tot_var = mean(EP_data(cc,bbb).tot_var);
                psth_noise_var = tot_var - EP_data(cc,bbb).pair_psth_var; %estimated noise variance using PSTH
                spline_noise_var = tot_var - EP_data(cc,bbb).spline_looEP.weights(1); %spline-based estimate of noise variance
                EP_data(cc,bbb).spline_FF = spline_noise_var/avg_rates(cc); %spline-based FF
                EP_data(cc,bbb).psth_FF = psth_noise_var/avg_rates(cc); %PSTH-based FF
                for bb = 1:length(poss_eps_sizes) %compute epsilon-ball based FF
                    cur_ball_noise_var = tot_var - EP_data(cc,bbb).eps_ball_var(bb);
                    EP_data(cc,bbb).ball_FF(bb) =  cur_ball_noise_var/avg_rates(cc);
                end
            end
        end
        
        %% compute pairwise covariances
        if do_xcorrs
            if length(targs) > 1 %if there's at least one SU pair
                for rr = 1:n_rpt_seeds; %loop over unique repeat seeds
                    cur_trial_set = find(all_rpt_seqnum == rr);
                    cur_nrpts = length(cur_trial_set);
                    [II,JJ] = meshgrid(1:cur_nrpts);
                    uset = JJ ~= II; %use all unequal trial pairs
                    n_unique_pairs = sum(uset(:));
                    
                    %make shift-embedded spike mat
                    allY1 = new_BS_ms(:,cur_trial_set,:);
                    allY2 = nan(n_Tbins,length(cur_trial_set),length(targs),length(tlags));
                    for tt = 1:length(tlags)
                        allY2(:,:,:,tt) = shift_matrix_Nd(squeeze(new_BS_ms(:,cur_trial_set,:)),tlags(tt),1,'nan'); %shift matrix, using nan padding
                    end
                    
                    %                     cur_D = nan(n_unique_pairs*n_Tbins,1);
                    %                     cur_D_LOO = nan(n_unique_pairs*n_Tbins,length(loo_set));
                    %                     cur_X = nan(n_unique_pairs*n_Tbins,length(targs),length(targs),length(tlags));
                    cur_X_sum = zeros(length(poss_eps_sizes),length(targs),length(targs),length(tlags));
                    cur_X_cnt = zeros(length(poss_eps_sizes),length(targs),length(targs),length(tlags));
                    cur_ovX_sum = zeros(length(targs),length(targs),length(tlags));
                    cur_ovX_cnt = zeros(length(targs),length(targs),length(tlags));
                    cur_X_sum_LOO = zeros(length(poss_eps_sizes),length(targs),length(targs),length(tlags),length(loo_set));
                    cur_X_cnt_LOO = zeros(length(poss_eps_sizes),length(targs),length(targs),length(tlags),length(loo_set));
                    for tt = 1:n_Tbins
                        if mod(tt,10) == 0
                            fprintf('Pair time bin %d/%d\n',tt,n_Tbins);
                        end
                        cur_inds = (tt-1)*n_unique_pairs + (1:n_unique_pairs);
                        Y1 = squeeze(allY1(tt,:,:)); %[NxCxC] mat
                        Y2 = reshape(squeeze(allY2(tt,:,:,:)),[length(cur_trial_set) 1 length(targs) length(tlags)]); %make this a [Nx1xCxL] mat
                        
                        %compute all trial-pair products Y*Y
                        cur_Xmat = bsxfun(@times,Y1(II(uset),:),Y2(JJ(uset),:,:,:));
                        cur_ovX_sum = cur_ovX_sum + squeeze(nansum(cur_Xmat));
                        cur_ovX_cnt = cur_ovX_cnt + squeeze(sum(~isnan(cur_Xmat)));
                        
                        %calculate deltaX mat using overall EP
                        cur_Dmat = abs(squareform(pdist(squeeze(new_EP_emb(tt,cur_trial_set,:)))))/sqrt(emb_win);
                        cur_Dmat(logical(eye(cur_nrpts))) = nan;
                        %                         cur_D(cur_inds) = cur_Dmat(uset);
                        
                        cur_Dmat = cur_Dmat(uset);
                        
                        %loop over possible epsilon ball sizes
                        for pp = 1:length(poss_eps_sizes)
                            cur_set = find(cur_Dmat < poss_eps_sizes(pp));
                            cur_X_sum(pp,:,:,:) = cur_X_sum(pp,:,:,:) + nansum(cur_Xmat(cur_set,:,:,:));
                            cur_X_cnt(pp,:,:,:) = cur_X_cnt(pp,:,:,:) + sum(~isnan(cur_Xmat(cur_set,:,:,:)));
                        end
                        
                        % calculate a separate deltaX mat for each LOO EP signal
                        for ll = 1:length(loo_set)
                            cur_Dmat = abs(squareform(pdist(squeeze(new_loo_EP_emb(tt,cur_trial_set,:,ll)))))/sqrt(emb_win);
                            cur_Dmat(logical(eye(cur_nrpts))) = nan;
                            
                            cur_Dmat = cur_Dmat(uset);
                            for pp = 1:length(poss_eps_sizes)
                                cur_set = find(cur_Dmat < poss_eps_sizes(pp));
                                cur_X_sum_LOO(pp,:,:,:,ll) = cur_X_sum_LOO(pp,:,:,:,ll) + nansum(cur_Xmat(cur_set,:,:,:));
                                cur_X_cnt_LOO(pp,:,:,:,ll) = cur_X_cnt_LOO(pp,:,:,:,ll) + sum(~isnan(cur_Xmat(cur_set,:,:,:)));
                            end
                            
                            %                             cur_D_LOO(cur_inds,ll) = cur_Dmat(uset);
                        end
                        
                        %                         cur_X(cur_inds,:,:,:) = cur_Xmat;
                    end
                    
                    %this gives PSTH-based estimator of signal covariances
                    %(marginalizing over deltaX)
                    %                     pair_xcovar = squeeze(nanmean(cur_X));
                    pair_xcovar = cur_ovX_sum./cur_ovX_cnt;
                    
                    %                 [bincnts,binids] = histc(cur_D,EP_bin_edges);
                    %                 var_ep_binned = nan(n_EP_bins,length(targs),length(targs),length(tlags));
                    %                 for bb = 1:n_EP_bins
                    %                     var_ep_binned(bb,:,:,:) = nanmean(cur_X(binids == bb,:,:,:));
                    %                 end
                    
                    %compute epsilon-ball based covariances
                    %                     eps_ball_var = nan(length(poss_eps_sizes),length(targs),length(targs),length(tlags));
                    %                     for bb = 1:length(poss_eps_sizes)
                    %                         curset = find(cur_D < poss_eps_sizes(bb));
                    %                         if ~isempty(curset)
                    %                             eps_ball_var(bb,:,:,:) = nanmean(cur_X(curset,:,:,:),1);
                    %                         end
                    %                     end
                    eps_ball_var = cur_X_sum./cur_X_cnt;
                    
                    eps_ball_var_LOO1 = nan(length(poss_eps_sizes),length(targs),length(targs),length(tlags));
                    eps_ball_var_LOO2 = nan(length(poss_eps_sizes),length(targs),length(targs),length(tlags));
                    for cc1 = 1:length(targs) %loop over pairs of neurons
                        for cc2 = 1:length(targs)
                            %                             for ll = 1:length(tlags)
                            
                            %find LOO indices of the units in this pair
                            loo_ind1 = find(ismember(targs(cc1),loo_set));
                            loo_ind2 = find(ismember(targs(cc2),loo_set));
                            if ~isempty(loo_ind1) && ~isempty(loo_ind2)
                                %compute epsilon-ball based covariances using these LOO EP deltaXs
                                %                                     for bb = 1:length(poss_eps_sizes)
                                %                                         curset = find(cur_D_LOO(:,loo_ind1) < poss_eps_sizes(bb));
                                %                                         eps_ball_var_LOO1(bb,cc1,cc2,:) = nanmean(cur_X(curset,cc1,cc2,:));
                                %                                         curset = find(cur_D_LOO(:,loo_ind2) < poss_eps_sizes(bb));
                                %                                         eps_ball_var_LOO2(bb,cc1,cc2,:) = nanmean(cur_X(curset,cc1,cc2,:));
                                %                                     end
                                eps_ball_var_LOO1(:,cc1,cc2,:) = cur_X_sum_LOO(:,cc1,cc2,:,loo_ind1)./cur_X_cnt_LOO(:,cc1,cc2,:,loo_ind1);
                                eps_ball_var_LOO2(:,cc1,cc2,:) = cur_X_sum_LOO(:,cc1,cc2,:,loo_ind2)./cur_X_cnt_LOO(:,cc1,cc2,:,loo_ind2);
                            end
                            %                             end
                        end
                    end
                    
                    %normalization based on average across-trial variance
                    %('psth-based' estimate of noise variance)
                    at_vars = squeeze(nanmean(nanvar(new_BS_ms,[],2)));
                    noisevar_norm = at_vars*at_vars';
                    
                    [CI,CJ] = meshgrid(1:length(targs));
                    un_pairs = CJ >= CI; %all unique pairs of neurons (including self-pairs)
                    Cpairs = [CI(un_pairs) CJ(un_pairs)]; %I and J indices for each pair
                    n_cell_pairs = size(Cpairs,1);
                    for cc = 1:n_cell_pairs
                        Y1 = reshape(squeeze(allY1(:,:,Cpairs(cc,1))),[],1); %response of neuron I
                        Y2 = reshape(squeeze(allY2(:,:,Cpairs(cc,2),:)),[],length(tlags)); %shift-embedded response of neuron J
                        
                        EP_pairs(cc,bbb).ids = Cpairs(cc,:); %store index values of the neurons in this pair
                        pair_rpt_set = intersect(EP_data(Cpairs(cc,1),1).rpt_trialset,EP_data(Cpairs(cc,2),1).rpt_trialset);
                        EP_pairs(cc,bbb).pair_rpt_set = pair_rpt_set;
                        EP_pairs(cc,bbb).tot_xcovar(rr,:) = nanmean(bsxfun(@times,Y1,Y2)); %raw spk count covariances
                        EP_pairs(cc,bbb).at_var_norm(rr) = sqrt(prod(noisevar_norm(Cpairs(cc,1),Cpairs(cc,2)))); %normalization by product of noise SDs
                        
                        %average covariance for this cell pair between using trial
                        %pairs ij in both orders
                        EP_pairs(cc,bbb).pair_xcovar(rr,:) = pair_xcovar(Cpairs(cc,1),Cpairs(cc,2),:);
                        
                        %                         %histogram based covariance vs deltaX
                        %                         EP_pairs(cc,bbb).xcovar_ep_binned(rr,:,:) = squeeze(var_ep_binned(:,Cpairs(cc,1),Cpairs(cc,2),:));
                        
                        %epsilon-ball based covariances (using total EP signal)
                        EP_pairs(cc,bbb).eps_xcovar(rr,:,:) = squeeze(eps_ball_var(:,Cpairs(cc,1),Cpairs(cc,2),:));
                        
                        %epsilon-ball based covariance using each LOO EP signal
                        %of the pair, and then avg them together.
                        eps_covar_LOO1 = squeeze(eps_ball_var_LOO1(:,Cpairs(cc,1),Cpairs(cc,2),:));
                        eps_covar_LOO2 = squeeze(eps_ball_var_LOO2(:,Cpairs(cc,1),Cpairs(cc,2),:));
                        EP_pairs(cc,bbb).eps_xcovar_LOO(rr,:,:) = 0.5*eps_covar_LOO1 + 0.5*eps_covar_LOO2;
                    end
                end
            else
                EP_pairs = [];
            end
        end
    end
    
    if ismember(bin_dt,mod_bin_dts)
        %% rebin model predicted firing rates
        if bin_dt > params.dt %if we're downsapmlign in time
            new_mod_prates = nan(n_Tbins,tot_nrpts,length(targs),bin_dsfac);
            new_mod_prates_noEM = nan(n_Tbins,tot_nrpts,length(targs),bin_dsfac);
            for ii = 1:bin_dsfac
                new_mod_prates(:,:,:,ii) = all_mod_emp_prates(ii:bin_dsfac:(ii+bin_dsfac*(n_Tbins-1)),:,:,:);
                new_mod_prates_noEM(:,:,:,ii) = all_mod_emp_prates_noEM(ii:bin_dsfac:(ii+bin_dsfac*(n_Tbins-1)),:,:,:);
            end
            new_mod_prates = squeeze(mean(new_mod_prates,4))*bin_dsfac; %if any in-bin values are NAN, the whole bin is NAN
            new_mod_prates_noEM = squeeze(mean(new_mod_prates_noEM,4))*bin_dsfac; %if any in-bin values are NAN, the whole bin is NAN
            
        elseif bin_dt < params.dt
%             error('Havent incorporated upsampling for model fits');
            new_mod_prates = interp1(1:used_nf,all_mod_emp_prates,bin_usfac:bin_usfac:used_nf)*bin_usfac;
            new_mod_prates(1:(1/bin_usfac-1),:,:) = new_mod_prates(1/bin_usfac,:,:); %handle the nans that arise from the first upsampled bin centers being out of range
            new_mod_prates_noEM = interp1(1:used_nf,all_mod_emp_prates_noEM,bin_usfac:bin_usfac:used_nf)*bin_usfac;
            new_mod_prates_noEM(1:(1/bin_usfac-1),:,:) = new_mod_prates_noEM(1/bin_usfac,:,:); %handle the nans that arise from the first upsampled bin centers being out of range
        else
            n_Tbins = used_nf;
            new_mod_prates = all_mod_emp_prates;
            new_mod_prates_noEM = all_mod_emp_prates_noEM;
        end
        
        new_mod_prates_noEM(isnan(new_BS_ms)) = nan; %for model data that is matched to spk data, use only where unit was isolated
        %% analyze model predicted firing rates
        for rr = 1:n_rpt_seeds
            cur_trial_set = find(all_rpt_seqnum == rr);
            n_utrials = squeeze(mean(sum(~isnan(new_mod_prates_noEM(:,cur_trial_set,:)),2))); %avg (across time) number of usable trials
            
            mod_psths_noEM = squeeze(nanmean(new_mod_prates_noEM(:,cur_trial_set,:),2)); %across-trial avg rate
            mod_cond_vars_noEM = squeeze(nanvar(new_mod_prates_noEM(:,cur_trial_set,:),[],2)); %across-trial variances
            mod_tot_vars_noEM = squeeze(nanvar(reshape(new_mod_prates_noEM(:,cur_trial_set,:),[],length(targs)))); %total rate variance
            mod_psth_vars_noEM = nanvar(mod_psths_noEM); %raw PSTH variance
            avg_temp_var = squeeze(nanmean(nanvar(new_mod_prates(:,cur_trial_set,:)))); %avg (across trials) of across-time variance
            mod_psth_vars_cor_noEM = mod_psth_vars_noEM.*(n_utrials'./(n_utrials'-1)) - avg_temp_var'./(n_utrials-1)'; %sahani linden correction for PSTH sampling noise
            
            mod_psths = squeeze(nanmean(new_mod_prates(:,cur_trial_set,:),2)); %across-trial avg rate
            mod_cond_vars = squeeze(nanvar(new_mod_prates(:,cur_trial_set,:),[],2)); %across-trial variances
            mod_tot_vars = squeeze(nanvar(reshape(new_mod_prates(:,cur_trial_set,:),[],length(targs)))); %total rate variance
            mod_psth_vars = nanvar(mod_psths); %raw PSTH variance
            
            %save data
            for cc = 1:length(targs)
                if has_stim_mod(cc)
                    EP_data(cc,bbb).mod_psths_noEM(rr,:) = mod_psths_noEM(:,cc);
                    EP_data(cc,bbb).mod_cond_vars_noEM(rr,:) = mod_cond_vars_noEM(:,cc);
                    EP_data(cc,bbb).mod_tot_vars_noEM(rr) = mod_tot_vars_noEM(cc);
                    EP_data(cc,bbb).mod_psth_vars_noEM(rr) = mod_psth_vars_noEM(cc);
                    EP_data(cc,bbb).mod_psth_vars_cor_noEM(rr) = mod_psth_vars_cor_noEM(cc);
                    EP_data(cc,bbb).mod_ep_vars_noEM(rr) = nanmean(mod_cond_vars_noEM(:,cc));
                    EP_data(cc,bbb).mod_alphas_noEM(rr) = EP_data(cc,bbb).mod_ep_vars_noEM(rr)/EP_data(cc,bbb).mod_tot_vars_noEM(rr);
                    
                    EP_data(cc,bbb).mod_psths(rr,:) = mod_psths(:,cc);
                    EP_data(cc,bbb).mod_cond_vars(rr,:) = mod_cond_vars(:,cc);
                    EP_data(cc,bbb).mod_tot_vars(rr) = mod_tot_vars(cc);
                    EP_data(cc,bbb).mod_psth_vars(rr) = mod_psth_vars(cc);
                    EP_data(cc,bbb).mod_ep_vars(rr) = nanmean(mod_cond_vars(:,cc));
                    EP_data(cc,bbb).mod_alphas(rr) = EP_data(cc,bbb).mod_ep_vars(rr)/EP_data(cc,bbb).mod_tot_vars(rr);
                end
            end
        end
        
        %% model-predicted rate covariances
        if do_xcorrs
            if length(targs) > 1
                allY1 = bsxfun(@minus,new_mod_prates,reshape(nanmean(reshape(new_mod_prates,[],length(targs))),[1 1 length(targs)])); %remove overall avg
                for rr = 1:n_rpt_seeds
                    cur_trial_set = find(all_rpt_seqnum == rr);
                    
                    %construct shifted rate array
                    allY2 = nan(n_Tbins,length(cur_trial_set),length(targs),length(tlags));
                    for tt = 1:length(tlags)
                        allY2(:,:,:,tt) = shift_matrix_Nd(squeeze(allY1(:,cur_trial_set,:)),tlags(tt),1,'nan');
                    end
                    
                    %loop over cell pairs
                    for cc = 1:n_cell_pairs
                        EP_pairs(cc,bbb).mod_tot_covar(rr,:) = nan(1,length(tlags));
                        EP_pairs(cc,bbb).mod_psth_covar(rr,:) = nan(1,length(tlags));
                        Y1 = squeeze(allY1(:,cur_trial_set,Cpairs(cc,1)));
                        for ll = 1:length(tlags)
                            Y2 = squeeze(allY2(:,:,Cpairs(cc,2),ll));
                            
                            EP_pairs(cc,bbb).mod_tot_covar(rr,ll) = nanmean(Y1(:).*Y2(:)); %raw rate covariance
                            EP_pairs(cc,bbb).mod_psth_covar(rr,ll) = nanmean(nanmean(Y1,2).*nanmean(Y2,2)); %PSTH covariance
                        end
                    end
                end
            end
        end
        %% use direct estimate to compute stats for simulated data
        if compute_sims
            
            %first get all pairwise deltaX vals
            all_delta_X = [];
            for rr = 1:n_rpt_seeds %loop over unique repeat seeds
                cur_trial_set = find(all_rpt_seqnum == rr);
                cur_nrpts = length(cur_trial_set);
                [II,JJ] = meshgrid(1:cur_nrpts);
                uset = JJ > II; %only need to count each unique trial pair once
                n_unique_pairs = sum(uset(:));
                
                cur_D = nan(n_unique_pairs*n_Tbins,1);
                for tt = 1:n_Tbins
                    cur_inds = (tt-1)*n_unique_pairs + (1:n_unique_pairs);
                    
                    cur_Dmat = abs(squareform(pdist(squeeze(new_EP_emb(tt,cur_trial_set,:)))))/sqrt(emb_win);
                    cur_Dmat(logical(eye(cur_nrpts))) = nan;
                    cur_D(cur_inds) = cur_Dmat(uset);
                end
                all_delta_X = cat(1,all_delta_X,cur_D);
            end
            
            %now run a series of spiking simulations
            clear sim_stats
            reverseStr = '';
            for sr = 1:sim_n_rpts
                msg = sprintf('Simulating rpt: %d/%d',sr,sim_n_rpts);
                fprintf([reverseStr msg]);
                reverseStr = repmat(sprintf('\b'),1,length(msg));
                
                sim_spikes = poissrnd(new_mod_prates_noEM); %generate rate-mod poisson spikes
                sim_trial_avg = nanmean(sim_spikes); %simulated trial-avgs
                %subtract out trial-avg spk counts if desired (to better
                %match real calcs)
                if sub_trialavgs
                    sim_spikes = bsxfun(@minus,sim_spikes,sim_trial_avg);
                end
                %now subtract out overall avg spike count
                sim_spikes = bsxfun(@minus,sim_spikes,reshape(nanmean(reshape(sim_spikes,[],length(targs))),[1 1 length(targs)]));
                
                %compute all pairwise products Y_I*Y_J
                all_X = [];
                for rr = 1:n_rpt_seeds %loop over unique repeat seeds
                    cur_trial_set = find(all_rpt_seqnum == rr);
                    cur_nrpts = length(cur_trial_set);
                    [II,JJ] = meshgrid(1:cur_nrpts);
                    uset = JJ > II; %only need to count each unique trial pair once
                    n_unique_pairs = sum(uset(:));
                    
                    cur_X = nan(n_unique_pairs*n_Tbins,length(targs));
                    for tt = 1:n_Tbins
                        cur_inds = (tt-1)*n_unique_pairs + (1:n_unique_pairs);
                        Y1 = squeeze(sim_spikes(tt,cur_trial_set,:));
                        if length(targs) == 1
                            Y1 = Y1';
                        end
                        cur_Xmat = bsxfun(@times,Y1(II,:),Y1(JJ,:));
                        cur_X(cur_inds,:) = cur_Xmat(uset,:);
                    end
                    all_X = cat(1,all_X,cur_X);
                end
                
                n_data_points = length(all_delta_X);
                
                %                 %get spline knot points
                %                 spline_DS = prctile(all_delta_X,maxD_prc/(n_spline_knots-1):maxD_prc/(n_spline_knots-1):(maxD_prc-maxD_prc/(n_spline_knots-1)));
                %                 knot_pts = [0 0 0 0 spline_DS maxD maxD maxD];
                
                %                 %for each cell compute the total rate variance from a
                %                 %spline regression
                %                 spline_tot_var = nan(1,length(targs));
                %                 for cc = 1:length(targs)
                %                     upts = find(all_delta_X <= knot_pts(end) & ~isnan(all_X(:,cc)));
                %                     sp = fastBSpline.lsqspline(knot_pts,3,all_delta_X(upts),all_X(upts,cc));
                %                     spline_tot_var(cc) = sp.evalAt(0);
                %                 end
                %                 all_psth_vars = nanmean(all_X); %marginal avg Y*Y gives estimate of "PSTH variance"
                %                 spline_alpha_ests = all_psth_vars./spline_tot_var; %spline-based estimates of alpha
                
                %compute epsilon-ball estimates
                eps_ball_vars = nan(length(poss_eps_sizes),length(targs));
                for bb = 1:length(poss_eps_sizes)
                    curset = find(all_delta_X < poss_eps_sizes(bb));
                    eps_ball_vars(bb,:) = nanmean(all_X(curset,:));
                end
                eps_alpha_ests = bsxfun(@rdivide,all_psth_vars,eps_ball_vars); %epsilon-based estimates of alpha
                
                for cc = 1:length(targs)
                    sim_stats(sr,cc).psth_vars = all_psth_vars(cc);
                    %                     sim_stats(sr,cc).spline_vars = spline_tot_var(cc);
                    sim_stats(sr,cc).eps_vars = eps_ball_vars(:,cc);
                    sim_stats(sr,cc).eps_alphas = all_psth_vars(cc)./eps_ball_vars(:,cc);
                    %                     sim_stats(sr,cc).spline_alpha = all_psth_vars(cc)/spline_tot_var(cc);
                end
            end
            
            %store the model-based stats for each neuron
            for cc = 1:length(targs)
                EP_data(cc,bbb).mod_sim_stats = sim_stats(:,cc);
            end
            fprintf('\n');
        end
    end
    
    %% export trial-by-trial data for making figures
    %     fig_dname = [anal_dir 'tbt_fig_data'];
    %     fig_data.tbt_EP = tbt_EP;
    %     fig_data.tbt_EP_emb = tbt_EP_emb;
    %     fig_data.pred_rates = new_mod_prates;
    %     fig_data.binned_spks = tbt_binned_spikes;
    %     fig_data.binned_spks_nan = new_BS_ms;
    %     fig_data.EP_data = EP_data;
    %     fig_data.EP_params = EP_params;
    %     fig_data.modFitParams = modFitParams;
    %
    %     if compute_PF_rate
    %         fig_data.PF_prates = all_mod_PF_prates;
    %     end
    %
    %     n_probes = 24;
    %     trial_dur = 4;
    %     poss_targs = (n_probes + 1):(n_probes + length(SU_numbers));
    %     for ss = 1:length(targs)
    %         cur_SU_ind = find(targs(ss) == poss_targs);
    %         cur_spk_times = spike_data.SU_spk_times{cur_SU_ind};
    %         for ii = 1:tot_nrpts
    %             trial_start_time = trial_data(all_rpt_trials(ii)).start_times;
    %             trial_spk_times = cur_spk_times(cur_spk_times >= trial_start_time & cur_spk_times <= (trial_start_time + trial_dur));
    %             tbt_spk_times{ss,ii} = trial_spk_times - trial_start_time;
    %         end
    %     end
    %     fig_data.tbt_spk_times = tbt_spk_times;
    %
    %     %handle any trials with repeat frames
    %     to_eliminate = [];
    %     for ii = 1:length(rptframe_trials)
    %         cur_trial = all_rpt_trials(rptframe_trials(ii));
    %         cur_rpt_frames = trial_data(cur_trial).rpt_frames;
    %         if all(cur_rpt_frames == 0) %zero values indicate which frame the trial really started on
    %             spk_shift_amount = length(cur_rpt_frames)*dt_uf; %number of up-sampled time steps to shift by
    %             for ss = 1:length(targs)
    %                 tbt_spk_times{ss,rptframe_trials(ii)} = tbt_spk_times{ss,rptframe_trials(ii)} + spk_shift_amount*bin_dt;
    %             end
    %         elseif ~any(cur_rpt_frames == 0) %in this case if there's a repeat frame in the middle of the trial, just remove the trial (for plotting purposes this would be a mess)
    %             to_eliminate = [to_eliminate rptframe_trials(ii)];
    %         end
    %     end
    %
    %     fig_data.tbt_spk_times(:,to_eliminate) = [];
    %     fig_data.tbt_EP(:,to_eliminate) = [];
    %     fig_data.tbt_EP_emb(:,to_eliminate) = [];
    %     fig_data.pred_rates(:,to_eliminate,:) = [];
    %     fig_data.binned_spks(:,to_eliminate,:) = [];
    %
    %     fprintf('Saving %s\n',fig_dname);
    %     save(fig_dname,'fig_data');
end

%%
cd(anal_dir);

sname = [sname sprintf('_ori%d',bar_ori)];
if rec_number > 1
    sname = strcat(sname,sprintf('_r%d',rec_number));
end
if do_xcorrs
    save(sname,'targs','EP_data','EP_pairs','EP_params','rpt_data');
else
    save(sname,'targs','EP_data','EP_params','rpt_data');
end