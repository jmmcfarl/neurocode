clear all
close all

addpath('~/other_code/fastBSpline/');

global Expt_name bar_ori monk_name rec_type rec_number

Expt_name = 'M296';
monk_name = 'lem';
bar_ori = 45; %bar orientation to use (only for UA recs)
rec_number = 1;
% %
% [266-80 270-60 275-135 277-70 281-140 287-90 289-160 294-40 296-45 297-0/90 5-50 9-0 10-60 11-160 12-0 13-100 14-40]

use_MUA = false; %use MUA in model-fitting
use_hres_ET = true; EP_params.use_hres_ET = use_hres_ET; %use high-res eye-tracking?
exclude_sacs = true; EP_params.exclude_sacs = exclude_sacs;
sub_trialavgs = true; EP_params.sub_trialavgs = sub_trialavgs; %subtract out trial avg spike counts

poss_bin_dts = [0.01 0.05 0.1 0.2]; EP_params.poss_bin_dts = poss_bin_dts;
max_tlag = 0;
tlags = -max_tlag:max_tlag;

use_LOOXV = 1; %[0 is no LOO; 1 is SUs only; 2 is SU + MU]

data_dir = ['~/Data/bruce/' Expt_name];
if ~exist(data_dir,'dir')
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
end
if ~exist(data_dir,'dir');
    error('Couldnt find data directory');
end
Edata_file = strcat(data_dir,'/',monk_name,Expt_name,'Expts');
load(Edata_file);

%is this a laminar probe or utah array rec?
ff = find(cellfun(@(x) ~isempty(x),Expts),1);
if strcmp(Expts{ff}.Header.DataType,'GridData 96')
    rec_type = 'UA';
elseif strcmp(Expts{ff}.Header.DataType,'Spike2')
    rec_type = 'LP';
end

%load in packaged data
data_name = sprintf('%s/packaged_data_ori%d',data_dir,bar_ori);
if rec_number > 1
    data_name = strcat(data_name,sprintf('_r%d',rec_number));
end
fprintf('Loading %s\n',data_name);
load(data_name);

ov_RF_pos = Expts{expt_data.used_blocks(1)}.Stimvals.rf(1:2)/params.scale_fac;

sname = 'rpt_variability_compact_multDT';
if sub_trialavgs
    sname = strcat(sname,'_subTrial');
end
%%
Expt_num = str2num(Expt_name(2:end));

cd(data_dir);
% load([data_dir '/stims/expt_data.mat']); %load in stim-alignment meta-data

et_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
cluster_dir = ['~/Analysis/bruce/' Expt_name '/clustering'];
model_dir = ['~/Analysis/bruce/' Expt_name '/models'];
et_mod_data_name = 'full_eyetrack_initmods_Rinit';
et_anal_name = 'full_eyetrack_Rinit';
mod_name = 'corrected_models_comp';

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
tr_set = et_tr_set;

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
% n_blocks = length(expt_data.used_blocks);
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

%% CREATE SACCADE PREDICTOR MATS
corrected_eye_vals_interp = ET_data.interp_eye_pos;
sac_start_times = [ET_data.saccades(:).start_time];
sac_stop_times = [ET_data.saccades(:).stop_time];

interp_sac_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_start_times));
interp_sac_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),sac_stop_times));

saccade_start_inds = find(ismember(used_inds,interp_sac_start_inds));
used_saccade_set = find(ismember(interp_sac_start_inds,used_inds));
%nearest index in the used data set of the saccade stop time
saccade_stop_inds = round(interp1(used_inds,1:length(used_inds),interp_sac_stop_inds(used_saccade_set)))';
saccade_stop_inds(isnan(saccade_stop_inds)) = length(used_inds);

used_is_blink = ET_data.is_blink(used_saccade_set);

saccade_trial_inds = all_trialvec(used_inds(saccade_start_inds));

saccades = ET_data.saccades(used_saccade_set);

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
fix_post_blink = ismember(fix_start_inds,saccade_stop_inds(used_is_blink));
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
    post_mean_EP = post_mean_EP*sp_dx;
    
    if exist('drift_post_mean_LOO','var')
        post_mean_EP_LOO = nan(length(loo_set),NT);
        for ss = 1:length(loo_set)
            [post_mean_EP_LOO(ss,:)] = construct_eye_position(best_fix_cor,best_fix_std,...
                drift_post_mean_LOO(ss,:),drift_post_std_LOO(ss,:),fix_ids,trial_start_inds,trial_end_inds,sac_shift);
        end
        post_mean_EP_LOO = post_mean_EP_LOO*sp_dx;
    else
        warning('No LOO variables detected');
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
n_rpt_seeds = length(params.rpt_seeds); EP_params.n_rpt_seeds = n_rpt_seeds;

all_rpt_trials = find(ismember(all_trial_Se,params.rpt_seeds));
rptframe_trials = find([trial_data(all_rpt_trials).nrpt_frames] > 0); %identify repeat trials where there were repeat frames
fprintf('Detected %d/%d trials with rpt frames\n',length(rptframe_trials),length(all_rpt_trials));

% all_rpt_trials(rptframe_trials) = [];
% rptframe_trials = [];

T = tabulate(all_trialvec);
rpt_tdurs = T(all_rpt_trials,2);

too_short = find(rpt_tdurs < 390);
if ~isempty(too_short)
    fprintf('Eliminating %d/%d repeat trials without enough frames\n',length(too_short),length(rpt_tdurs));
    all_rpt_trials(too_short) = [];
    rptframe_trials = find([trial_data(all_rpt_trials).nrpt_frames] > 0); %get rid of any repeat trials where there were repeat frames
end

if any(arrayfun(@(x) any(x.rpt_frames == 0),trial_data(all_rpt_trials)))
    mixed_trials = false(length(all_rpt_trials),1);
    for ii = 1:length(all_rpt_trials)
        if any(trial_data(all_rpt_trials(ii)).rpt_frames == 0) & any(trial_data(all_rpt_trials(ii)).rpt_frames > 0)
            mixed_trials(ii) = true;
        end
    end
    if any(mixed_trials)
        fprintf('Eliminating %d/%d repeat trials with mixed rpt frame types\n',sum(mixed_trials),length(mixed_trials));
        all_rpt_trials(mixed_trials) = [];
        rptframe_trials = find([trial_data(all_rpt_trials).nrpt_frames] > 0); %get rid of any repeat trials where there were repeat frames
    end
end

all_rpt_seqnum = nan(size(all_rpt_trials));
for ii = 1:n_rpt_seeds
    cur_trials = find(all_trial_Se(all_rpt_trials) == params.rpt_seeds(ii));
    all_rpt_seqnum(cur_trials) = ii;
end

tot_nrpts = length(all_rpt_seqnum);
fprintf('Using %d repeat trials, %d sequences\n',tot_nrpts,length(params.rpt_seeds));

rpt_trial_block = [trial_data(all_rpt_trials).block_nums];

%% IDENTIFY TIMES WITHIN SACCADES AND BLINKS
sac_buff = round(0.05/params.dt); EP_params.sac_buff = sac_buff; %window of data to exclude during saccades
sac_delay = round(0.03/params.dt); EP_params.sac_delay = sac_delay; %shift exclusion window to account for neural delay
blink_buff = round(0.1/params.dt); EP_params.blink_buff = blink_buff; %window for excluding blinks

% in_sac_inds = zeros(NT,1);
% nblink_start_inds = saccade_start_inds(~used_is_blink);
% for ii = 1:(sac_buff+1)
%     cur_inds = nblink_start_inds + sac_delay + (ii-1);
%     uu = find(cur_inds <= NT);
%     uu(all_trialvec(used_inds(cur_inds(uu))) ~= all_trialvec(used_inds(nblink_start_inds(uu)))) = [];
%     in_sac_inds(cur_inds(uu)) = 1;
% end
%
% blink_start_inds = saccade_start_inds(used_is_blink);
% in_blink_inds = zeros(NT,1);
% for ii = 1:(blink_buff+1)
%     cur_inds = blink_start_inds + sac_delay + (ii-1);
%     uu = find(cur_inds <= NT);
%     uu(all_trialvec(used_inds(cur_inds(uu))) ~= all_trialvec(used_inds(blink_start_inds(uu)))) = [];
%     in_blink_inds(cur_inds(uu)) = 1;
% end

in_sac_inds = false(NT,1);
nblink_start_inds = saccade_start_inds(~used_is_blink);
nblink_stop_inds = saccade_stop_inds(~used_is_blink);
for ii = 1:length(nblink_start_inds)
    cur_inds = (nblink_start_inds(ii):(nblink_stop_inds(ii) + sac_buff)) + sac_delay;
    cur_inds(cur_inds > NT) = [];
    in_sac_inds(cur_inds) = true;
end

in_blink_inds = false(NT,1);
blink_start_inds = saccade_start_inds(used_is_blink);
blink_stop_inds = saccade_stop_inds(used_is_blink);
for ii = 1:length(blink_start_inds)
    cur_inds = (blink_start_inds(ii):(blink_stop_inds(ii) + blink_buff)) + sac_delay;
    cur_inds(cur_inds > NT) = [];
    in_blink_inds(cur_inds) = true;
end

used_rpt_inds = find(ismember(all_trialvec(used_inds),all_rpt_trials)); %indices of repeat trials within used_inds vector
nf = 400;
used_nf = nf-(params.beg_buffer + params.end_buffer)/params.dt;


% in_blink_inds = reshape(ceil(interp1(all_t_axis(used_inds),in_blink_inds,tbt_t_axis(up_used_inds))),used_up_nf,tot_nrpts);
in_sac_inds = reshape(in_sac_inds(used_rpt_inds),used_nf,tot_nrpts);
in_blink_inds = reshape(in_blink_inds(used_rpt_inds),used_nf,tot_nrpts);
in_sac_inds(isnan(in_sac_inds)) = 0; in_blink_inds(isnan(in_blink_inds)) = 0;
in_sac_inds = logical(in_sac_inds); in_blink_inds = logical(in_blink_inds);

%% PROCESS MODEL FITS
has_stim_mod = false(length(targs),1);
for cc = 1:length(targs)
    if ~isempty(ModData(targs(cc)).bestGQM)
        cur_mod = ModData(targs(cc)).bestGQM;
        cur_block_filt = cur_mod.mods(1).filtK;
        cur_used_blocks = ModData(targs(cc)).unit_data.used_blocks;
        poss_used_blocks = ModData(targs(cc)).unit_data.poss_used_blocks;
        cur_used_blocks = find(ismember(cur_used_blocks,poss_used_blocks));
        cur_mod.spk_NL_params(1) = cur_mod.spk_NL_params(1) + mean(cur_block_filt(cur_used_blocks));
        cur_mod.mods(1) = []; %eliminate block filter
        stim_mod(cc) = cur_mod;
        EP_data(cc,1).unit_data = ModData(targs(cc)).unit_data;
        EP_data(cc,1).tune_props = ModData(targs(cc)).tune_props;
        EP_data(cc,1).bestGQM = ModData(targs(cc)).bestGQM;
        EP_data(cc,1).nullMod = ModData(targs(cc)).nullMod;
        has_stim_mod(cc) = true;
    end
end

%% shift original EP data to realign repeats on trials with rpt frames
post_mean_EP_rpt = post_mean_EP(used_rpt_inds);
tbt_EP = reshape(post_mean_EP_rpt,used_nf,tot_nrpts);
used_NF = (params.beg_buffer/params.dt+1):(params.trial_dur - params.end_buffer)/params.dt;

rpt_blanked = false(nf,tot_nrpts);
shifted_EP = tbt_EP;
if ~isempty(rptframe_trials)
    post_rpt_buffer = round(0.1/params.dt); %exclude data for this duration following each rpt frame
    
    shifted_EP = cat(1,nan(params.beg_buffer/params.dt,tot_nrpts),tbt_EP,nan(params.end_buffer/params.dt+1,tot_nrpts));
    shifted_EP((nf+1):end,:) = [];
    for ii = 1:length(rptframe_trials)
        cur_trial = all_rpt_trials(rptframe_trials(ii));
        cur_rpt_frames = trial_data(cur_trial).rpt_frames;
        if all(cur_rpt_frames == 0)
            shift_amount = length(cur_rpt_frames);
            shifted_EP(:,rptframe_trials(ii)) = shift_matrix_Nd(shifted_EP(:,rptframe_trials(ii)),shift_amount,1);
            to_blank_inds = false(nf,1);
            to_blank_inds(1:(shift_amount + params.beg_buffer/params.dt)) = true;
            rpt_blanked(to_blank_inds,rptframe_trials(ii)) = true;
            
        elseif ~any(cur_rpt_frames == 0)
            new_frame_ids = 1:nf;
            for jj = 1:length(cur_rpt_frames)
                target_inds = (cur_rpt_frames(jj) + 1):nf;
                map_to = target_inds + 1; map_to(map_to > nf) = nf;
                new_frame_ids(target_inds) = new_frame_ids(map_to);
                
                rpt_blanked(cur_rpt_frames(jj):(cur_rpt_frames(jj)+post_rpt_buffer),rptframe_trials(ii)) = true;
                rpt_blanked((nf-length(cur_rpt_frames)*1):nf,rptframe_trials(ii)) = true;
            end
            shifted_EP(:,rptframe_trials(ii)) = shifted_EP(new_frame_ids,rptframe_trials(ii));
            
        else
            warning('mix of 0 and nonzero rpt frames');
        end
        
    end
    shifted_EP = shifted_EP(used_NF,:);
end
usedrpt_blanked = rpt_blanked(used_NF,:);

%% get stimulus xmat during repeat trials

all_stim_mat = decompressTernNoise(stimComp);

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
% used_rpt_inds = find(ismember(all_trialvec(used_inds),all_rpt_trials));

% fin_shift_cor = round(post_mean_EP/modFitParams.sp_dx); %use overall EP estimate
fin_shift_cor = round(reshape(shifted_EP,[],1)/modFitParams.sp_dx); %use overall EP estimate
fin_shift_cor(isnan(fin_shift_cor)) = 0;

%RECOMPUTE XMAT
best_shift_stimmat_up = all_stimmat_up;
for i=1:length(used_rpt_inds)
    %     best_shift_stimmat_up(used_inds(used_rpt_inds(i)),:) = shift_matrix_Nd(all_stimmat_up(used_inds(used_rpt_inds(i)),:),-fin_shift_cor(used_rpt_inds(i)),2);
    best_shift_stimmat_up(used_inds(used_rpt_inds(i)),:) = shift_matrix_Nd(all_stimmat_up(used_inds(used_rpt_inds(i)),:),-fin_shift_cor(i),2);
end
all_Xmat_shift = create_time_embedding(best_shift_stimmat_up(full_rpt_inds,:),stim_params_full);
all_Xmat_shift = all_Xmat_shift(ismember(full_rpt_inds,used_inds),use_kInds_up);

%%
all_mod_emp_prates = nan(length(used_rpt_inds),length(targs));
for cc = 1:length(targs)
    if has_stim_mod(cc)
        [~,~,all_mod_emp_prates(:,cc)] = NMMmodel_eval(stim_mod(cc),[],all_Xmat_shift);
    end
end
all_mod_emp_prates(in_blink_inds(:),:) = nan;
if exclude_sacs
    all_mod_emp_prates(in_sac_inds(:),:) = nan;
end
all_mod_emp_prates = reshape(all_mod_emp_prates,used_nf,tot_nrpts,length(targs));

%% shift model-predicted rates to realign repeats on trials with rpt frames
if ~isempty(rptframe_trials)
    shifted_mod_prates = nan(nf,tot_nrpts,length(targs));
    shifted_mod_prates(used_NF,:,:) = all_mod_emp_prates;
    for ii = 1:length(rptframe_trials)
        cur_trial = all_rpt_trials(rptframe_trials(ii));
        cur_rpt_frames = trial_data(cur_trial).rpt_frames;
        new_spike_frame_ids = 1:nf;
        if all(cur_rpt_frames == 0)
            shift_amount = length(cur_rpt_frames);
            shifted_mod_prates(:,rptframe_trials(ii),:) = shift_matrix_Nd(shifted_mod_prates(:,rptframe_trials(ii),:),shift_amount,1);
        else
            for jj = 1:length(cur_rpt_frames)
                target_inds = (cur_rpt_frames(jj) + 1):nf;
                map_to = target_inds + 1; map_to(map_to > nf) = nf;
                new_spike_frame_ids(target_inds) = new_spike_frame_ids(map_to);
            end
            shifted_mod_prates(:,rptframe_trials(ii),:) = shifted_mod_prates(new_spike_frame_ids,rptframe_trials(ii),:);
        end
    end
    shifted_mod_prates = shifted_mod_prates(used_NF,:,:);
    shifted_mod_prates = reshape(shifted_mod_prates,[],length(targs));
    shifted_mod_prates(usedrpt_blanked(:),:) = nan;
    shifted_mod_prates = reshape(shifted_mod_prates,length(used_NF),tot_nrpts,length(targs));
    all_mod_emp_prates = shifted_mod_prates;
end

%%
for bbb = 1:length(poss_bin_dts)
    bin_dt = poss_bin_dts(bbb);
    fprintf('Running analysis at dt = %.3f sec\n',bin_dt);
    
    %% Construct Time embedded eye position sequences for repeat trials
    emb_win = (bin_dt + 0.05); EP_params.emb_win = emb_win; %look back this many time steps to parse EP trajectories
    emb_shift = 0.03; EP_params.emb_shift = emb_shift;
    emb_win = round(emb_win/orig_dt);
    emb_shift = round(emb_shift/orig_dt);
    
    %initialize a time-embedded version of the EPs
    sp = NMMcreate_stim_params(emb_win + emb_shift);
    tbt_EP_emb = create_time_embedding(tbt_EP(:),sp);
    tbt_EP_emb(in_blink_inds(:),:) = nan;
    if exclude_sacs
        tbt_EP_emb(in_sac_inds(:),:) = nan;
    end
    tbt_EP_emb = reshape(tbt_EP_emb(:,(emb_shift+1):end),used_nf,tot_nrpts,[]);
    
    %% compile all time-embedded LOO EP sequences
    ms = size(tbt_EP_emb);
    loo_tbt_EP_emb = nan(ms(1),ms(2),ms(3),length(loo_set));
    
    for cc = 1:length(targs)
        loo_ind = find(loo_set == targs(cc));
        if ~isempty(loo_ind)
            %         interp_post_mean_EP = interp1(time_data.t_axis(used_inds),post_mean_EP_LOO(loo_ind,:),tbt_t_axis(up_used_inds));
            post_mean_rpt = post_mean_EP_LOO(loo_ind,used_rpt_inds);
            cur_tbt_EP = reshape(post_mean_rpt,used_nf,tot_nrpts);
            
            %initialize a time-embedded version of the EPs
            sp = NMMcreate_stim_params(emb_win+emb_shift);
            cur_tbt_EP_emb = create_time_embedding(cur_tbt_EP(:),sp);
            cur_tbt_EP_emb(in_blink_inds(:),:) = nan;
            if exclude_sacs
                cur_tbt_EP_emb(in_sac_inds(:),:) = nan;
            end
            loo_tbt_EP_emb(:,:,:,loo_ind) = reshape(cur_tbt_EP_emb(:,(emb_shift+1):end),used_nf,tot_nrpts,[]);
        end
    end
    
    %% shift EP data to realign repeats on trials with rpt frames
    if ~isempty(rptframe_trials)
        post_rpt_buffer = round(0.1/params.dt); %exclude data for this duration following each rpt frame
        
        tbt_EP_emb = cat(1,nan(params.beg_buffer/params.dt,tot_nrpts,emb_win),tbt_EP_emb,nan(params.end_buffer/params.dt+1,tot_nrpts,emb_win));
        loo_tbt_EP_emb = cat(1,nan(params.beg_buffer/params.dt,tot_nrpts,emb_win,length(loo_set)),loo_tbt_EP_emb,nan(params.end_buffer/params.dt+1,tot_nrpts,emb_win,length(loo_set)));
        tbt_EP_emb((nf+1):end,:,:) = [];
        loo_tbt_EP_emb((nf+1):end,:,:,:) = [];
        for ii = 1:length(rptframe_trials)
            cur_trial = all_rpt_trials(rptframe_trials(ii));
            cur_rpt_frames = trial_data(cur_trial).rpt_frames;
            if all(cur_rpt_frames == 0)
                shift_amount = length(cur_rpt_frames);
                tbt_EP_emb(:,rptframe_trials(ii),:) = shift_matrix_Nd(tbt_EP_emb(:,rptframe_trials(ii),:),shift_amount,1);
                loo_tbt_EP_emb(:,rptframe_trials(ii),:,:) = shift_matrix_Nd(loo_tbt_EP_emb(:,rptframe_trials(ii),:,:),shift_amount,1);
            elseif ~any(cur_rpt_frames == 0)
                new_frame_ids = 1:nf;
                for jj = 1:length(cur_rpt_frames)
                    target_inds = (cur_rpt_frames(jj) + 1):nf;
                    map_to = target_inds + 1; map_to(map_to > nf) = nf;
                    new_frame_ids(target_inds) = new_frame_ids(map_to);
                end
                tbt_EP_emb(:,rptframe_trials(ii),:) = tbt_EP_emb(new_frame_ids,rptframe_trials(ii),:);
                loo_tbt_EP_emb(:,rptframe_trials(ii),:,:) = loo_tbt_EP_emb(new_frame_ids,rptframe_trials(ii),:,:);
                
            else
                warning('mix of 0 and nonzero rpt frames');
            end
            
        end
        tbt_EP_emb = tbt_EP_emb(used_NF,:,:);
        loo_tbt_EP_emb = loo_tbt_EP_emb(used_NF,:,:,:);
    end
    tbt_EP_emb = reshape(tbt_EP_emb,[],emb_win);
    loo_tbt_EP_emb = reshape(loo_tbt_EP_emb,[],emb_win,length(loo_set));
    tbt_EP_emb(usedrpt_blanked(:),:) = nan;
    loo_tbt_EP_emb(usedrpt_blanked(:),:,:) = nan;
    tbt_EP_emb = reshape(tbt_EP_emb,used_nf,tot_nrpts,emb_win);
    loo_tbt_EP_emb = reshape(loo_tbt_EP_emb,used_nf,tot_nrpts,emb_win,length(loo_set));
    
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
    
    n_Tbins = length(full_bin_taxis);
    n_used_Tbins = length(bin_taxis);
    used_Tinds = find(ismember(full_bin_taxis,bin_taxis));
    
    tbt_binned_spikes = nan(n_Tbins,tot_nrpts,length(SU_numbers));
    tbt_t_axis = nan(n_Tbins,tot_nrpts);
    for ii = 1:tot_nrpts
        cur_bin_edges = [trial_data(all_rpt_trials(ii)).start_times:cur_bin_dt:(trial_data(all_rpt_trials(ii)).start_times + cur_bin_dt*(n_Tbins))];
        cur_bin_cents = 0.5*cur_bin_edges(1:end-1) + 0.5*cur_bin_edges(2:end);
        for cc = 1:length(SU_numbers)
            if ~isnan(Clust_data.SU_block_probes(cc,rpt_trial_block(ii)))
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
    %nan
    tbt_binned_spikes = reshape(tbt_binned_spikes,[],length(SU_numbers));
    for ii = 1:length(SU_numbers)
        tbt_binned_spikes(isnan(all_binned_sua(orig_t_ind,ii)),ii) = nan;
    end
    tbt_binned_spikes = reshape(tbt_binned_spikes,n_Tbins,tot_nrpts,length(SU_numbers));
    tbt_binned_spikes(:,:,length(targs)+1:end) = [];
    
    %%
    if bin_dt < params.dt
        spk_in_blink_inds = round(interp1(all_t_axis(used_inds(used_rpt_inds)),double(in_blink_inds(:)),reshape(tbt_t_axis(used_Tinds,:),[],1)));
        spk_in_blink_inds(isnan(spk_in_blink_inds)) = 0;
        spk_in_blink_inds = logical(spk_in_blink_inds);
        spk_in_sac_inds = round(interp1(all_t_axis(used_inds(used_rpt_inds)),double(in_sac_inds(:)),reshape(tbt_t_axis(used_Tinds,:),[],1)));
        spk_in_sac_inds(isnan(spk_in_sac_inds)) = 0;
        spk_in_sac_inds = logical(spk_in_sac_inds);
    else
        spk_in_blink_inds = in_blink_inds;
        spk_in_sac_inds = in_sac_inds;
    end
    temp = reshape(tbt_binned_spikes(used_Tinds,:,:),[],length(targs));
    temp(spk_in_blink_inds(:),:) = nan;
    if exclude_sacs
       temp(spk_in_sac_inds(:),:) = nan; 
    end
    tbt_binned_spikes(used_Tinds,:,:) = reshape(temp,length(used_Tinds),tot_nrpts,length(targs));
    
    %% shift spike data to realign repeats on trials with rpt frames
    if ~isempty(rptframe_trials)
        dt_uf = params.dt/cur_bin_dt;
        
        for ii = 1:length(rptframe_trials)
            cur_trial = all_rpt_trials(rptframe_trials(ii));
            cur_rpt_frames = trial_data(cur_trial).rpt_frames;
            if all(cur_rpt_frames == 0)
                spk_shift_amount = length(cur_rpt_frames)*dt_uf;
                tbt_binned_spikes(:,rptframe_trials(ii),:) = shift_matrix_Nd(tbt_binned_spikes(:,rptframe_trials(ii),:),spk_shift_amount,1);
                
            elseif ~any(cur_rpt_frames == 0)
                new_frame_ids = 1:nf;
                new_spike_frame_ids = 1:nf*dt_uf;
                for jj = 1:length(cur_rpt_frames)
                    
                    target_inds = (cur_rpt_frames(jj) + 1)*dt_uf:nf*dt_uf;
                    map_to = target_inds + dt_uf; map_to(map_to > nf*dt_uf) = nf*dt_uf;
                    new_spike_frame_ids(target_inds) = new_spike_frame_ids(map_to);
                end
                tbt_binned_spikes(:,rptframe_trials(ii),:) = tbt_binned_spikes(new_spike_frame_ids,rptframe_trials(ii),:);
                
            else
                warning('mix of 0 and nonzero rpt frames');
            end
        end
    end
    
    %% process trial-by-trial binned spike data (exclude blinks, sacs, subtract trial avgs)
    tbt_BS_ms = reshape(tbt_binned_spikes(used_Tinds,:,:),[],length(targs));
    if bin_dt < params.dt
        spk_rptblank_inds = round(interp1(all_t_axis(used_inds(used_rpt_inds)),double(usedrpt_blanked(:)),reshape(tbt_t_axis(used_Tinds,:),[],1)));
        spk_rptblank_inds(isnan(spk_rptblank_inds)) = 0;
        spk_rptblank_inds = logical(spk_rptblank_inds);
    else
        spk_rptblank_inds = usedrpt_blanked;
    end
    
    tbt_BS_ms(spk_rptblank_inds(:),:) = nan;
    tbt_BS_ms = reshape(tbt_BS_ms,length(used_Tinds),tot_nrpts,length(targs));
    
    %% if using coarser than dt time binning, do rebinning here
    if bin_dt > params.dt
        bin_usfac = bin_dt/params.dt;
        if mod(bin_usfac,1) ~= 0
            error('have to use integer multiple of dt for time binning');
        end
        
        n_Tbins = floor(used_nf/bin_usfac);
        
        new_BS_ms = nan(n_Tbins,tot_nrpts,length(targs),bin_usfac);
        for ii = 1:bin_usfac
            new_BS_ms(:,:,:,ii) = tbt_BS_ms(ii:bin_usfac:(ii+bin_usfac*(n_Tbins-1)),:,:,:);
        end
        new_BS_ms = squeeze(sum(new_BS_ms,4));
%         new_BS_ms = squeeze(nansum(new_BS_ms,4));
        
        %take the eye position history leading up to the last time bin
        ep_bin_ids = bin_usfac:bin_usfac:(bin_usfac +bin_usfac*(n_Tbins-1));
        new_EP_emb = tbt_EP_emb(ep_bin_ids,:,:);
        new_loo_EP_emb = loo_tbt_EP_emb(ep_bin_ids,:,:,:);
        
    elseif bin_dt < params.dt %f using finer-than-dt binning
        bin_usfac = bin_dt/params.dt;
        n_Tbins = length(used_Tinds);
        new_BS_ms = tbt_BS_ms;
        ep_bin_ids = ceil(bin_usfac:bin_usfac:used_nf);
        new_EP_emb = tbt_EP_emb(ep_bin_ids,:,:);
        new_loo_EP_emb = loo_tbt_EP_emb(ep_bin_ids,:,:,:);
    else %if using native binning
        n_Tbins = length(used_Tinds);
        new_BS_ms = tbt_BS_ms;
        new_EP_emb = tbt_EP_emb;
        new_loo_EP_emb = loo_tbt_EP_emb;
    end
    
    %% basic avg spk count calculations
    %first compute avg spike rates
    ov_avg_BS = nanmean(reshape(new_BS_ms,[],length(targs)));
    trial_avg_BS = nanmean(new_BS_ms);
    for cc = 1:length(targs)
        EP_data(cc,bbb).ov_avg_BS = ov_avg_BS(cc);
        EP_data(cc,bbb).trial_avg_BS = squeeze(trial_avg_BS(:,:,cc));
    end
    for rr = 1:n_rpt_seeds
        cur_trial_set = find(all_rpt_seqnum == rr);
        cur_nrpts = length(cur_trial_set);
        
        n_utrials = squeeze(mean(sum(~isnan(new_BS_ms(:,cur_trial_set,:)),2)));
        n_spikes = squeeze(nansum(reshape(new_BS_ms,[],length(SU_numbers))));
        
        for cc = 1:length(targs)
            EP_data(cc,bbb).n_utrials(rr) = n_utrials(cc);
            EP_data(cc,bbb).n_spikes(rr) = n_spikes(cc);
        end
    end
    
    %subtract out trial-avg spk counts if desired
    if sub_trialavgs
        new_BS_ms = bsxfun(@minus,new_BS_ms,trial_avg_BS);
    else
        %if not subtracting trial avgs subtract out within-block avgs
        for ii = 1:length(expt_data.used_blocks)
            cur_block_trials = find(rpt_trial_block == ii);
            if ~isempty(cur_block_trials)
                cur_block_avgs = nanmean(reshape(new_BS_ms(:,cur_block_trials,:),[],length(targs)));
                new_BS_ms(:,cur_block_trials,:) = bsxfun(@minus,new_BS_ms(:,cur_block_trials,:),reshape(cur_block_avgs,1,1,length(targs)));
            end
        end
    end
    trial_avg_BS = squeeze(trial_avg_BS);
    
    %now subtract out overall avg spike count
    new_BS_ms = bsxfun(@minus,new_BS_ms,reshape(nanmean(reshape(new_BS_ms,[],length(targs))),[1 1 length(targs)]));
    
    %% COMPUTE BASIC STATS FOR REPEAT TRIALS
    for rr = 1:n_rpt_seeds
        cur_trial_set = find(all_rpt_seqnum == rr);
        cur_nrpts = length(cur_trial_set);
        
        psths = squeeze(nanmean(new_BS_ms(:,cur_trial_set,:),2));
        psth_var = nanvar(psths);
        tot_resp_var = nanvar(reshape(new_BS_ms(:,cur_trial_set,:),[],length(SU_numbers)));
        
        avg_temp_var = squeeze(nanmean(nanvar(new_BS_ms(:,cur_trial_set,:)))); %avg (across trials) of across-time variance
        psth_var_cor = psth_var.*(n_utrials'./(n_utrials'-1)) - avg_temp_var'./n_utrials'; %sahani linden correction for PSTH sampling noise
        
        avg_acrossTrial_var = squeeze(nanmean(nanvar(new_BS_ms(:,cur_trial_set,:),[],2)));
        trial_avg_var = squeeze(nanvar(trial_avg_BS)); %variance of trial-avg rates
        
        for cc = 1:length(targs)
            EP_data(cc,bbb).psths(rr,:) = psths(:,cc);
            EP_data(cc,bbb).psth_var(rr) = psth_var(cc);
            EP_data(cc,bbb).psth_var_cor(rr) = psth_var_cor(cc);
            EP_data(cc,bbb).across_trial_var(rr) = avg_acrossTrial_var(cc);
            EP_data(cc,bbb).tot_var(rr) = tot_resp_var(cc);
            EP_data(cc,bbb).trial_avg_var(rr) = trial_avg_var(cc);
        end
    end
    
    %% ESTIMATE OVERALL DISTRIBUTION OF DELTA_X to determine quantiles
    
    % compute the distribution of delta_X
    rand_T = [];
    for rr = 1:n_rpt_seeds
        cur_trial_set = find(all_rpt_seqnum == rr);
        
        % estimate quantiles of the distribution of pairwise EP similarities by this metric
        rset = randi(length(used_NF),100,1);
        for ii = 1:length(rset)
            cur_Dmat = abs(squareform(pdist(squeeze(tbt_EP_emb(rset(ii),cur_trial_set,:)))))/sqrt(emb_win);
            cur_Dmat(logical(eye(length(cur_trial_set)))) = nan;
            cur_Dmat = cur_Dmat(~isnan(cur_Dmat));
            rand_T = cat(1,rand_T,cur_Dmat);
        end
    end
    
    %% MAIN WITHIN-CELL ANALYSIS LOOP
    maxD_prc = 50;
    
    % xvfold = 10; EP_params.xvfold = xvfold; %cross-val fold for estimating optimal number of splines
    poss_n_splines = [3:8]; EP_params.poss_N_splines = poss_n_splines; %range of possible values for number of splines.
    best_n_knots = 4; EP_params.best_n_knots = 4;
    % n_boot_samps = 2; EP_params.n_boot_samps = n_boot_samps; %number of bootstrap samples for estimating spline uncertainty
    n_EP_bins = 50; EP_params.n_EP_bins = n_EP_bins;
    EP_bin_edges = prctile(rand_T,linspace(0,maxD_prc,n_EP_bins+1));
    EP_bin_centers = (EP_bin_edges(1:end-1)+EP_bin_edges(2:end))/2;  EP_params.EP_bin_centers = EP_bin_centers;
    maxD = prctile(rand_T,maxD_prc);
    
    poss_eps_sizes = [.005 .01 .02]; EP_params.poss_eps_sizes = poss_eps_sizes;
    
    n_eval_pts = 100; EP_params.n_eval_pts = n_eval_pts; %number of points to evaluate spline fit
    eval_xx = unique([0 prctile(rand_T,linspace(0,maxD_prc,n_eval_pts))]); %x-axis for evaluating spline models
    EP_params.eval_xx = eval_xx;
    
    for cc = 1:length(targs)
        fprintf('SU %d/%d\n',cc,length(targs));
        loo_ind = find(loo_set == targs(cc));
        if ~isempty(loo_ind)
            cur_tbt_EP_emb = squeeze(new_loo_EP_emb(:,:,:,loo_ind));
            
            all_D = [];
            all_base_D = [];
            %         all_sing_D = [];
            all_X = [];
            for rr = 1:n_rpt_seeds %loop over unique repeat seeds
                cur_trial_set = find(all_rpt_seqnum == rr);
                cur_nrpts = length(cur_trial_set);
                [II,JJ] = meshgrid(1:cur_nrpts);
                uset = JJ > II; %only need to count each unique trial pair once
                n_unique_pairs = sum(uset(:));
                
                cur_D = nan(n_unique_pairs*n_Tbins,1);
                cur_base_D = nan(n_unique_pairs*n_Tbins,1);
                cur_sing_D = nan(n_unique_pairs*n_Tbins,1);
                cur_X = nan(n_unique_pairs*n_Tbins,1);
                for tt = 1:n_Tbins
                    cur_inds = (tt-1)*n_unique_pairs + (1:n_unique_pairs);
                    Y1 = squeeze(new_BS_ms(tt,cur_trial_set,cc));
                    
                    cur_Dmat = abs(squareform(pdist(squeeze(cur_tbt_EP_emb(tt,cur_trial_set,:)))))/sqrt(emb_win);
                    cur_Dmat(logical(eye(cur_nrpts))) = nan;
                    cur_D(cur_inds) = cur_Dmat(uset);
                    %                 cur_Dmat = abs(squareform(pdist(squeeze(tbt_EP_emb(tt,cur_trial_set,:)))))/sqrt(emb_win);
                    cur_Dmat = abs(squareform(pdist(squeeze(new_EP_emb(tt,cur_trial_set,:)))))/sqrt(emb_win);
                    cur_Dmat(logical(eye(cur_nrpts))) = nan;
                    cur_base_D(cur_inds) = cur_Dmat(uset);
                    
                    cur_Xmat = bsxfun(@times,Y1(II),Y1(JJ));
                    cur_X(cur_inds) = cur_Xmat(uset);
                end
                cur_upts = find(~isnan(cur_D) & ~isnan(cur_X));
                all_D = cat(1,all_D,cur_D(cur_upts));
                all_base_D = cat(1,all_base_D,cur_base_D(cur_upts));
                %             all_sing_D = cat(1,all_sing_D,cur_sing_D(cur_upts));
                all_X = cat(1,all_X,cur_X(cur_upts));
            end
            n_data_points = length(all_X);
            
            [bincnts,binids] = histc(all_D,EP_bin_edges);
            var_ep_binned = nan(n_EP_bins,1);
            for bb = 1:n_EP_bins
                var_ep_binned(bb) = mean(all_X(binids == bb));
            end
            
            spline_DS = prctile(all_D,maxD_prc/(best_n_knots-1):maxD_prc/(best_n_knots-1):(maxD_prc-maxD_prc/(best_n_knots-1)));
            knot_pts = [0 0 0 0 spline_DS maxD maxD maxD];
            upts = find(all_D <= knot_pts(end));
            
            %             eps_ball_boot = nan(eps_boots,length(poss_eps_sizes));
            eps_ball_var = nan(length(poss_eps_sizes),1);
            eps_ball_cnts = nan(length(poss_eps_sizes),1);
            for bb = 1:length(poss_eps_sizes)
                curset = find(all_D < poss_eps_sizes(bb));
                eps_ball_var(bb) = mean(all_X(curset));
                eps_ball_cnts(bb) = length(curset);
                
                %                 for dd = 1:eps_boots
                %                     newset = curset(randi(eps_ball_cnts(bb),eps_ball_cnts(bb),1));
                %                     eps_ball_boot(dd,bb) = mean(all_X(newset));
                %                 end
            end
            EP_data(cc,bbb).eps_ball_var = eps_ball_var;
            %             EP_data(cc).eps_ball_sd = std(eps_ball_boot);
            
            sp = fastBSpline.lsqspline(knot_pts,3,all_D(upts),all_X(upts));
            EP_data(cc,bbb).spline_pred_looEP = sp.evalAt(eval_xx);
            EP_data(cc,bbb).spline_looEP = sp;
            
            upts = find(all_base_D <= knot_pts(end));
            sp = fastBSpline.lsqspline(knot_pts,3,all_base_D(upts),all_X(upts));
            EP_data(cc,bbb).spline_pred_baseEP = sp.evalAt(eval_xx);
            
            %         upts = find(all_sing_D <= knot_pts(end));
            %         sp = fastBSpline.lsqspline(knot_pts,3,all_sing_D(upts),all_X(upts));
            %         EP_data(cc).spline_pred_singEP = sp.evalAt(eval_xx);
            
            EP_data(cc,bbb).pair_psth_var = mean(all_X);
            EP_data(cc,bbb).n_knots = best_n_knots;
            EP_data(cc,bbb).var_ep_binned = var_ep_binned;
        end
    end
    
    %% COMPUTE FFs
    avg_rates = [EP_data(:,bbb).ov_avg_BS];
    ucells = find(~isnan(avg_rates));
    for cc = ucells
        tot_var = mean(EP_data(cc,bbb).tot_var);
        psth_noise_var = tot_var - EP_data(cc,bbb).pair_psth_var;
        spline_noise_var = tot_var - EP_data(cc,bbb).spline_pred_looEP(1);
        EP_data(cc,bbb).spline_FF = spline_noise_var/avg_rates(cc);
        EP_data(cc,bbb).psth_FF = psth_noise_var/avg_rates(cc);
        for bb = 1:length(poss_eps_sizes)
            cur_ball_noise_var = tot_var - EP_data(cc,bbb).eps_ball_var(bb);
            EP_data(cc,bbb).ball_FF(bb) =  cur_ball_noise_var/avg_rates(cc);
        end
    end
    
    %%
    %     covar_epsilon = 0.01;
    if length(targs) > 1 %if there's at least one SU pair
        for rr = 1:n_rpt_seeds;
            cur_trial_set = find(all_rpt_seqnum == rr);
            cur_nrpts = length(cur_trial_set);
            [II,JJ] = meshgrid(1:cur_nrpts);
            uset = JJ > II; %only need to count each unique trial pair once
            n_unique_pairs = sum(uset(:));
            
            %randomly partition trials into two non-overlapping sets
            tset1 = randperm(length(cur_trial_set));
            tset1(round(length(cur_trial_set)/2)+1:end) = [];
            tset2 = setdiff(1:length(cur_trial_set),tset1);
            
            %indices for pairs of unequal trials in each group
            uset_inds = find(uset(:));
            uset1 = find(ismember(JJ(uset_inds),tset1) & ismember(II(uset_inds),tset1)); %using only trials in set 1
            uset2 = find(ismember(JJ(uset_inds),tset2) & ismember(II(uset_inds),tset2)); %using only trials in set 2
            
            %make shift-embedded spike mat
            %         allY1 = tbt_BS_ms(:,cur_trial_set,:);
            allY1 = new_BS_ms(:,cur_trial_set,:);
            allY2 = nan(n_Tbins,length(cur_trial_set),length(targs),length(tlags));
            for tt = 1:length(tlags)
                %             allY2(:,:,:,tt) = shift_matrix_Nd(squeeze(tbt_BS_ms(:,cur_trial_set,:)),tlags(tt),1);
                allY2(:,:,:,tt) = shift_matrix_Nd(squeeze(new_BS_ms(:,cur_trial_set,:)),tlags(tt),1);
            end
            
            %             %the 3-element dimension stores the total averages, as well as
            %             %the averages restricted to subset1 and subset2 of trials
            %             cur_sumX = zeros(length(targs),length(targs),length(tlags),3,length(poss_eps_sizes));
            %             cur_cnts = zeros(length(targs),length(targs),length(tlags),3,length(poss_eps_sizes));
            %             cur_sumX_LOO = zeros(length(targs),length(targs),length(tlags),length(loo_set),3,length(poss_eps_sizes));
            %             cur_cnts_LOO = zeros(length(targs),length(targs),length(tlags),length(loo_set),3,length(poss_eps_sizes));
            %             cur_sumX_rand = zeros(length(targs),length(targs),length(tlags),3);
            %             cur_cnts_rand = zeros(length(targs),length(targs),length(tlags),3);
            %             for tt = 1:n_Tbins
            % %                 tt
            %                 Y1 = squeeze(allY1(tt,:,:));
            %                 Y2 = reshape(squeeze(allY2(tt,:,:,:)),[length(cur_trial_set) 1 length(targs) length(tlags)]);
            %                 cur_Xmat = bsxfun(@times,Y1(II(uset),:),Y2(JJ(uset),:,:,:));
            %
            %                 cur_Dmat = abs(squareform(pdist(squeeze(new_EP_emb(tt,cur_trial_set,:)))))/sqrt(emb_win);
            %
            %                 %get Dmatrices for each LOO eye-tracking vector
            %                 loo_Dmats = nan(length(loo_set),length(cur_trial_set),length(cur_trial_set));
            %                 for ll = 1:length(loo_set)
            %                     loo_Dmats(ll,:,:) = abs(squareform(pdist(squeeze(new_loo_EP_emb(tt,cur_trial_set,:,ll)))))/sqrt(emb_win);
            %                 end
            %
            %                 for pp = 1:length(poss_eps_sizes)
            %                     cur_epsilon = poss_eps_sizes(pp);
            %                     upairs = cur_Dmat(uset) < cur_epsilon;
            %                     cur_sumX(:,:,:,1,pp) = cur_sumX(:,:,:,1,pp) + squeeze(nansum(cur_Xmat(upairs,:,:,:),1));
            %                     cur_cnts(:,:,:,1,pp) = cur_cnts(:,:,:,1,pp) + squeeze(sum(~isnan(cur_Xmat(upairs,:,:,:)),1));
            %                     upairs = uset1(cur_Dmat(uset_inds(uset1)) < cur_epsilon);
            %                     cur_sumX(:,:,:,2,pp) = cur_sumX(:,:,:,2,pp) + squeeze(nansum(cur_Xmat(upairs,:,:,:),1));
            %                     cur_cnts(:,:,:,2,pp) = cur_cnts(:,:,:,2,pp) + squeeze(sum(~isnan(cur_Xmat(upairs,:,:,:)),1));
            %                     upairs = uset2(cur_Dmat(uset_inds(uset2)) < cur_epsilon);
            %                     cur_sumX(:,:,:,3,pp) = cur_sumX(:,:,:,3,pp) + squeeze(nansum(cur_Xmat(upairs,:,:,:),1));
            %                     cur_cnts(:,:,:,3,pp) = cur_cnts(:,:,:,3,pp) + squeeze(sum(~isnan(cur_Xmat(upairs,:,:,:)),1));
            %
            %                     for ll = 1:length(loo_set)
            %                         temp_Dmat = squeeze(loo_Dmats(ll,:,:));
            %                         upairs = temp_Dmat(uset) < cur_epsilon;
            %                         cur_sumX_LOO(:,:,:,ll,1,pp) = cur_sumX_LOO(:,:,:,ll,1,pp) + squeeze(nansum(cur_Xmat(upairs,:,:,:),1));
            %                         cur_cnts_LOO(:,:,:,ll,1,pp) = cur_cnts_LOO(:,:,:,ll,1,pp) + squeeze(sum(~isnan(cur_Xmat(upairs,:,:,:)),1));
            %                         upairs = uset1(temp_Dmat(uset_inds(uset1)) < cur_epsilon);
            %                         cur_sumX_LOO(:,:,:,ll,2,pp) = cur_sumX_LOO(:,:,:,ll,2,pp) + squeeze(nansum(cur_Xmat(upairs,:,:,:),1));
            %                         cur_cnts_LOO(:,:,:,ll,2,pp) = cur_cnts_LOO(:,:,:,ll,2,pp) + squeeze(sum(~isnan(cur_Xmat(upairs,:,:,:)),1));
            %                         upairs = uset2(temp_Dmat(uset_inds(uset2)) < cur_epsilon);
            %                         cur_sumX_LOO(:,:,:,ll,3,pp) = cur_sumX_LOO(:,:,:,ll,3,pp) + squeeze(nansum(cur_Xmat(upairs,:,:,:),1));
            %                         cur_cnts_LOO(:,:,:,ll,3,pp) = cur_cnts_LOO(:,:,:,ll,3,pp) + squeeze(sum(~isnan(cur_Xmat(upairs,:,:,:)),1));
            %                     end
            %                 end
            %                 cur_sumX_rand(:,:,:,1) = cur_sumX_rand(:,:,:,1) + squeeze(nansum(cur_Xmat));
            %                 cur_cnts_rand(:,:,:,1) = cur_cnts_rand(:,:,:,1) + squeeze(sum(~isnan(cur_Xmat),1));
            %                 cur_sumX_rand(:,:,:,2) = cur_sumX_rand(:,:,:,2) + squeeze(nansum(cur_Xmat(uset1,:,:,:)));
            %                 cur_cnts_rand(:,:,:,2) = cur_cnts_rand(:,:,:,2) + squeeze(sum(~isnan(cur_Xmat(uset1,:,:,:)),1));
            %                 cur_sumX_rand(:,:,:,3) = cur_sumX_rand(:,:,:,3) + squeeze(nansum(cur_Xmat(uset2,:,:,:)));
            %                 cur_cnts_rand(:,:,:,3) = cur_cnts_rand(:,:,:,3) + squeeze(sum(~isnan(cur_Xmat(uset2,:,:,:)),1));
            %
            %                 %                 [bincnts,binids] = histc(cur_Dmat(uset),EP_bin_edges);
            % %                 for bb = 1:n_EP_bins
            % %                     var_ep_binned(tt,bb,:,:) = nanmean(cur_Xmat(binids == bb,:,:),1);
            % %                 end
            %             end
            %             EP_xcovar = cur_sumX./cur_cnts;
            %             pair_xcovar = cur_sumX_rand./cur_cnts_rand;
            %             EP_xcovar_LOO = cur_sumX_LOO./cur_cnts_LOO;
            
            %             %make tlags the last dimension
%             EP_xcovar = permute(EP_xcovar,[1 2 4 5 3]);
%             EP_xcovar_LOO = permute(EP_xcovar_LOO,[1 2 4 5 6 3]);
%             pair_xcovar = permute(pair_xcovar,[1 2 4 3]);
%             
%             %subtract off product of mean rates within each set of data to obtain
%             %covariance estimates
%             mean_rates = nanmean(reshape(allY1,[],length(targs)));
%             mean_rates_t1 = nanmean(reshape(allY1(:,tset1,:),[],length(targs)));
%             mean_rates_t2 = nanmean(reshape(allY1(:,tset2,:),[],length(targs)));
%             pair_xcovar(:,:,1,:) = bsxfun(@minus,pair_xcovar(:,:,1,:),mean_rates'*mean_rates);
%             pair_xcovar(:,:,2,:) = bsxfun(@minus,pair_xcovar(:,:,2,:),mean_rates_t1'*mean_rates_t1);
%             pair_xcovar(:,:,3,:) = bsxfun(@minus,pair_xcovar(:,:,3,:),mean_rates_t2'*mean_rates_t2);
%             EP_xcovar(:,:,1,:,:) = bsxfun(@minus,EP_xcovar(:,:,1,:,:),mean_rates'*mean_rates);
%             EP_xcovar(:,:,2,:,:) = bsxfun(@minus,EP_xcovar(:,:,2,:,:),mean_rates_t1'*mean_rates_t1);
%             EP_xcovar(:,:,3,:,:) = bsxfun(@minus,EP_xcovar(:,:,3,:,:),mean_rates_t2'*mean_rates_t2);
%             EP_xcovar_LOO(:,:,:,1,:,:) = bsxfun(@minus,EP_xcovar_LOO(:,:,:,1,:,:),mean_rates'*mean_rates);
%             EP_xcovar_LOO(:,:,:,2,:,:) = bsxfun(@minus,EP_xcovar_LOO(:,:,:,2,:,:),mean_rates_t1'*mean_rates_t1);
%             EP_xcovar_LOO(:,:,:,3,:,:) = bsxfun(@minus,EP_xcovar_LOO(:,:,:,3,:,:),mean_rates_t2'*mean_rates_t2);
            

            cur_D = nan(n_unique_pairs*n_Tbins,1);
            cur_D_LOO = nan(n_unique_pairs*n_Tbins,length(loo_set));
            cur_X = nan(n_unique_pairs*n_Tbins,length(targs),length(targs),length(tlags));
            for tt = 1:n_Tbins
                cur_inds = (tt-1)*n_unique_pairs + (1:n_unique_pairs);
                Y1 = squeeze(allY1(tt,:,:));
                Y2 = reshape(squeeze(allY2(tt,:,:,:)),[length(cur_trial_set) 1 length(targs) length(tlags)]);
                
                cur_Dmat = abs(squareform(pdist(squeeze(new_EP_emb(tt,cur_trial_set,:)))))/sqrt(emb_win);
                cur_Dmat(logical(eye(cur_nrpts))) = nan;
                cur_D(cur_inds) = cur_Dmat(uset);

                for ll = 1:length(loo_set)
                    cur_Dmat = abs(squareform(pdist(squeeze(new_loo_EP_emb(tt,cur_trial_set,:,ll)))))/sqrt(emb_win);
                    cur_Dmat(logical(eye(cur_nrpts))) = nan;
                    cur_D_LOO(cur_inds,ll) = cur_Dmat(uset);
                end
                
                cur_Xmat = bsxfun(@times,Y1(II(uset),:),Y2(JJ(uset),:,:,:));
                cur_X(cur_inds,:,:,:) = cur_Xmat;
            end

            pair_xcovar = squeeze(nanmean(cur_X));
            
            [bincnts,binids] = histc(cur_D,EP_bin_edges);
            var_ep_binned = nan(n_EP_bins,length(targs),length(targs),length(tlags));
            for bb = 1:n_EP_bins
                var_ep_binned(bb,:,:,:) = nanmean(cur_X(binids == bb,:,:,:));
            end
            eps_ball_var = nan(length(poss_eps_sizes),length(targs),length(targs),length(tlags));
            eps_ball_var_LOO1 = nan(length(poss_eps_sizes),length(targs),length(targs),length(tlags));
            eps_ball_var_LOO2 = nan(length(poss_eps_sizes),length(targs),length(targs),length(tlags));
            for bb = 1:length(poss_eps_sizes)
                curset = find(cur_D < poss_eps_sizes(bb));
                if ~isempty(curset)
                eps_ball_var(bb,:,:,:) = nanmean(cur_X(curset,:,:,:),1);
                end
            end
            
            spline_DS = prctile(cur_D,maxD_prc/(best_n_knots-1):maxD_prc/(best_n_knots-1):(maxD_prc-maxD_prc/(best_n_knots-1)));
            knot_pts = [0 0 0 0 spline_DS maxD maxD maxD];
            upts = find(cur_D <= knot_pts(end));
            all_spline_pred = nan(length(eval_xx),length(targs),length(targs),length(tlags));
            all_spline_pred_LOO1 = nan(length(eval_xx),length(targs),length(targs),length(tlags));
            all_spline_pred_LOO2 = nan(length(eval_xx),length(targs),length(targs),length(tlags));
            for cc1 = 1:length(targs)
                for cc2 = 1:length(targs)
                    for ll = 1:length(tlags)
                        cur_upts = upts(~isnan(cur_X(upts,cc1,cc2,ll)));
                        sp = fastBSpline.lsqspline(knot_pts,3,cur_D(cur_upts),cur_X(cur_upts,cc1,cc2,ll));
                        all_spline_pred(:,cc1,cc2,ll) = sp.evalAt(eval_xx);
                        
                        loo_ind1 = find(ismember(targs(cc1),loo_set));
                        loo_ind2 = find(ismember(targs(cc2),loo_set));
                        if ~isempty(loo_ind1) && ~isempty(loo_ind2)
                        cur_upts = find(cur_D_LOO(:,loo_ind1) <= knot_pts(end) & ~isnan(cur_X(:,cc1,cc2,ll)));
                        sp = fastBSpline.lsqspline(knot_pts,3,cur_D_LOO(cur_upts,loo_ind1),cur_X(cur_upts,cc1,cc2,ll));
                        all_spline_pred_LOO1(:,cc1,cc2,ll) = sp.evalAt(eval_xx);
                        cur_upts = find(cur_D_LOO(:,loo_ind2) <= knot_pts(end) & ~isnan(cur_X(:,cc1,cc2,ll)));
                        sp = fastBSpline.lsqspline(knot_pts,3,cur_D_LOO(cur_upts,loo_ind2),cur_X(cur_upts,cc1,cc2,ll));
                        all_spline_pred_LOO2(:,cc1,cc2,ll) = sp.evalAt(eval_xx);
                        
                        for bb = 1:length(poss_eps_sizes)
                            curset = find(cur_D_LOO(:,loo_ind1) < poss_eps_sizes(bb));
                            eps_ball_var_LOO1(bb,cc1,cc2,:) = nanmean(cur_X(curset,cc1,cc2,:));
                            curset = find(cur_D_LOO(:,loo_ind2) < poss_eps_sizes(bb));
                            eps_ball_var_LOO2(bb,cc1,cc2,:) = nanmean(cur_X(curset,cc1,cc2,:));
                        end
                        end
                    end
                end
            end
             
            %normalization based on average across-trial variance
            at_vars = squeeze(nanmean(nanvar(new_BS_ms,[],2)));
            noisevar_norm = at_vars*at_vars';
            
            [CI,CJ] = meshgrid(1:length(targs));
            un_pairs = CJ >= CI;
            Cpairs = [CI(un_pairs) CJ(un_pairs)];
            n_cell_pairs = size(Cpairs,1);
            tot_xcovar = nan(n_cell_pairs,length(tlags));
            for cc = 1:n_cell_pairs
                Y1 = reshape(squeeze(allY1(:,:,Cpairs(cc,1))),[],1);
                Y2 = reshape(squeeze(allY2(:,:,Cpairs(cc,2),:)),[],length(tlags));
                
                tot_xcovar(cc,:) = nanmean(bsxfun(@times,Y1,Y2));
                
                EP_pairs(cc,bbb).ids = Cpairs(cc,:);
                EP_pairs(cc,bbb).tot_xcovar(rr,:) = tot_xcovar(cc,:);
                EP_pairs(cc,bbb).at_var_norm(rr) = sqrt(prod(noisevar_norm(Cpairs(cc,1),Cpairs(cc,2))));
                
                %average covariance for this cell pair between using trial
                %pairs ij in both orders
                EP_pairs(cc,bbb).pair_xcovar(rr,:) = 0.5*pair_xcovar(Cpairs(cc,1),Cpairs(cc,2),:) + ...
                    0.5*pair_xcovar(Cpairs(cc,2),Cpairs(cc,1),:);
                
%                 EP_pairs(cc,bbb).EP_xcovar(rr,:,:,:) = 0.5*squeeze(EP_xcovar(Cpairs(cc,1),Cpairs(cc,2),:,:,:)) + ...
%                     0.5*squeeze(EP_xcovar(Cpairs(cc,2),Cpairs(cc,1),:,:,:));
                EP_pairs(cc,bbb).spline_xcovar(rr,:,:) = 0.5*squeeze(all_spline_pred(1,Cpairs(cc,1),Cpairs(cc,2),:)) + ...
                    0.5*squeeze(all_spline_pred(1,Cpairs(cc,2),Cpairs(cc,1),:));
                spline_xcovar_LOO1 = 0.5*squeeze(all_spline_pred_LOO1(1,Cpairs(cc,1),Cpairs(cc,2),:)) + ...
                    0.5*squeeze(all_spline_pred_LOO1(1,Cpairs(cc,2),Cpairs(cc,1),:));
                spline_xcovar_LOO2 = 0.5*squeeze(all_spline_pred_LOO2(1,Cpairs(cc,1),Cpairs(cc,2),:)) + ...
                    0.5*squeeze(all_spline_pred_LOO2(1,Cpairs(cc,2),Cpairs(cc,1),:));
                EP_pairs(cc,bbb).spline_xcovar_LOO(rr,:,:) = 0.5*spline_xcovar_LOO1 + 0.5*spline_xcovar_LOO2;
                
                EP_pairs(cc,bbb).splinefun_xcovar(rr,:,:) = 0.5*squeeze(all_spline_pred(:,Cpairs(cc,1),Cpairs(cc,2),:)) + ...
                    0.5*squeeze(all_spline_pred(:,Cpairs(cc,2),Cpairs(cc,1),:));
                
                EP_pairs(cc,bbb).xcovar_ep_binned(rr,:,:) = 0.5*squeeze(var_ep_binned(:,Cpairs(cc,1),Cpairs(cc,2),:)) + ...
                    0.5*squeeze(var_ep_binned(:,Cpairs(cc,2),Cpairs(cc,1),:));
                
                EP_pairs(cc,bbb).eps_xcovar(rr,:,:) = 0.5*squeeze(eps_ball_var(:,Cpairs(cc,1),Cpairs(cc,2),:)) + ...
                    0.5*squeeze(eps_ball_var(:,Cpairs(cc,2),Cpairs(cc,1),:));
                eps_covar_LOO1 = 0.5*squeeze(eps_ball_var_LOO1(:,Cpairs(cc,1),Cpairs(cc,2),:)) + ...
                    0.5*squeeze(eps_ball_var_LOO1(:,Cpairs(cc,2),Cpairs(cc,1),:));
                eps_covar_LOO2 = 0.5*squeeze(eps_ball_var_LOO2(:,Cpairs(cc,1),Cpairs(cc,2),:)) + ...
                    0.5*squeeze(eps_ball_var_LOO2(:,Cpairs(cc,2),Cpairs(cc,1),:));
                EP_pairs(cc,bbb).eps_xcovar_LOO(rr,:,:) = 0.5*eps_covar_LOO1 + 0.5*eps_covar_LOO2;
                
%                 %compute EP_corrected covariance using LOO EPs for each
%                 %unit of the pair, averaged together
%                 loo_ind1 = find(ismember(targs(Cpairs(cc,1)),loo_set));
%                 loo_ind2 = find(ismember(targs(Cpairs(cc,2)),loo_set));
%                 if ~isempty(loo_ind1) && ~isempty(loo_ind2)
%                     loo_EP_xcov1 = 0.5*squeeze(EP_xcovar_LOO(Cpairs(cc,1),Cpairs(cc,2),loo_ind1,:,:,:)) + 0.5*squeeze(EP_xcovar_LOO(Cpairs(cc,2),Cpairs(cc,1),loo_ind1,:,:,:));
%                     loo_EP_xcov2 = 0.5*squeeze(EP_xcovar_LOO(Cpairs(cc,1),Cpairs(cc,2),loo_ind2,:,:,:)) + 0.5*squeeze(EP_xcovar_LOO(Cpairs(cc,2),Cpairs(cc,1),loo_ind2,:,:,:));
%                     EP_pairs(cc,bbb).EP_xcovar_LOO(rr,:,:,:) = 0.5*loo_EP_xcov1 + 0.5*loo_EP_xcov2;
%                 else
%                     EP_pairs(cc,bbb).EP_xcovar_LOO(rr,:,:,:) = nan(3,length(poss_eps_sizes),length(tlags));
%                 end
                
%                 %store noisecov estimates using [1=all trials, 2=subset2,
%                 %and 3=subset1]. These noisecov estimates can be compared
%                 %against signal cov estimates based on these same sets of
%                 %trials, allowing comparisons of signal and noise covs
%                 %based on non-overlapping trial subsets
%                 EP_pairs(cc,bbb).psth_noisecov_ests(rr,1,:) = reshape(tot_xcovar(cc,:),[1,1,length(tlags)]) - EP_pairs(cc,bbb).pair_xcovar(rr,1,:);
%                 EP_pairs(cc,bbb).psth_noisecov_ests(rr,2,:) = reshape(tot_xcovar(cc,:),[1 1 length(tlags)]) - EP_pairs(cc,bbb).pair_xcovar(rr,3,:);
%                 EP_pairs(cc,bbb).psth_noisecov_ests(rr,3,:) = reshape(tot_xcovar(cc,:),[1 1 length(tlags)]) - EP_pairs(cc,bbb).pair_xcovar(rr,2,:);
%                 
%                 EP_pairs(cc,bbb).EP_noisecov_ests(rr,1,:,:) = bsxfun(@plus,-EP_pairs(cc,bbb).EP_xcovar(rr,1,:,:),reshape(tot_xcovar(cc,:),[1 1 1 length(tlags)]));
%                 EP_pairs(cc,bbb).EP_noisecov_ests(rr,2,:,:) = bsxfun(@plus,-EP_pairs(cc,bbb).EP_xcovar(rr,3,:,:),reshape(tot_xcovar(cc,:),[1 1 1 length(tlags)]));
%                 EP_pairs(cc,bbb).EP_noisecov_ests(rr,3,:,:) = bsxfun(@plus,-EP_pairs(cc,bbb).EP_xcovar(rr,2,:,:),reshape(tot_xcovar(cc,:),[1 1 1 length(tlags)]));
%                 EP_pairs(cc,bbb).EP_noisecov_LOO_ests(rr,1,:,:) = bsxfun(@plus,-EP_pairs(cc,bbb).EP_xcovar_LOO(rr,1,:,:),reshape(tot_xcovar(cc,:),[1 1 1 length(tlags)]));
%                 EP_pairs(cc,bbb).EP_noisecov_LOO_ests(rr,2,:,:) = bsxfun(@plus,-EP_pairs(cc,bbb).EP_xcovar_LOO(rr,3,:,:),reshape(tot_xcovar(cc,:),[1 1 1 length(tlags)]));
%                 EP_pairs(cc,bbb).EP_noisecov_LOO_ests(rr,3,:,:) = bsxfun(@plus,-EP_pairs(cc,bbb).EP_xcovar_LOO(rr,2,:,:),reshape(tot_xcovar(cc,:),[1 1 1 length(tlags)]));
            end
        end
    else
        EP_pairs = [];
    end
    
    %% rebin model predicted firing rates
    if bin_dt > params.dt
        new_mod_prates = nan(n_Tbins,tot_nrpts,length(targs),bin_usfac);
        for ii = 1:bin_usfac
            new_mod_prates(:,:,:,ii) = all_mod_emp_prates(ii:bin_usfac:(ii+bin_usfac*(n_Tbins-1)),:,:,:);
        end
        new_mod_prates = squeeze(nanmean(new_mod_prates,4))*bin_usfac;
        
    else
        n_Tbins = used_nf;
        new_mod_prates = all_mod_emp_prates;
    end
    %% analyze model predcted firing rates
    for rr = 1:n_rpt_seeds
        cur_trial_set = find(all_rpt_seqnum == rr);
        mod_psths = squeeze(nanmean(new_mod_prates(:,cur_trial_set,:),2));
        mod_cond_vars = squeeze(nanvar(new_mod_prates(:,cur_trial_set,:),[],2));
        mod_tot_vars = squeeze(nanvar(reshape(new_mod_prates(:,cur_trial_set,:),[],length(targs))));
        mod_psth_vars = nanvar(mod_psths);
        
        n_utrials = squeeze(mean(sum(~isnan(new_mod_prates(:,cur_trial_set,:)),2)));
        avg_temp_var = squeeze(nanmean(nanvar(new_mod_prates(:,cur_trial_set,:)))); %avg (across trials) of across-time variance
        mod_psth_vars_cor = mod_psth_vars.*(n_utrials'./(n_utrials'-1)) - avg_temp_var'./n_utrials'; %sahani linden correction for PSTH sampling noise
        
        
        for cc = 1:length(targs)
            if has_stim_mod(cc)
                EP_data(cc,bbb).mod_psths(rr,:) = mod_psths(:,cc);
                EP_data(cc,bbb).mod_cond_vars(rr,:) = mod_cond_vars(:,cc);
                EP_data(cc,bbb).mod_tot_vars(rr) = mod_tot_vars(cc);
                EP_data(cc,bbb).mod_psth_vars(rr) = mod_psth_vars(cc);
                EP_data(cc,bbb).mod_psth_vars_cor(rr) = mod_psth_vars_cor(cc);
                
                EP_data(cc,bbb).mod_ep_vars(rr) = mean(mod_cond_vars(:,cc));
                EP_data(cc,bbb).mod_alphas(rr) = EP_data(cc,bbb).mod_psth_vars_cor(rr)/EP_data(cc,bbb).mod_tot_vars(rr);
            end
        end
    end
    
    %% model-predicted rate covariances
    if length(targs) > 1
        allY1 = bsxfun(@minus,new_mod_prates,reshape(nanmean(reshape(new_mod_prates,[],length(targs))),[1 1 length(targs)]));
        for rr = 1:n_rpt_seeds
            cur_trial_set = find(all_rpt_seqnum == rr);
            
            mod_psths = squeeze(nanmean(allY1(:,cur_trial_set,:),2));
            mod_cond_vars = squeeze(nanvar(allY1(:,cur_trial_set,:),[],2));
            mod_tot_vars = squeeze(nanvar(reshape(allY1(:,cur_trial_set,:),[],length(targs))));
            mod_psth_vars = nanvar(mod_psths);
            
            n_utrials = squeeze(mean(sum(~isnan(allY1(:,cur_trial_set,:)),2)));
            avg_temp_var = squeeze(nanmean(nanvar(allY1(:,cur_trial_set,:)))); %avg (across trials) of across-time variance
            mod_psth_vars_cor = mod_psth_vars.*(n_utrials'./(n_utrials'-1)) - avg_temp_var'./n_utrials'; %sahani linden correction for PSTH sampling noise
            
            allY2 = nan(n_Tbins,length(cur_trial_set),length(targs),length(tlags));
            for tt = 1:length(tlags)
                allY2(:,:,:,tt) = shift_matrix_Nd(squeeze(allY1(:,cur_trial_set,:)),tlags(tt),1);
            end
            
            for cc = 1:n_cell_pairs
                EP_pairs(cc,bbb).mod_tot_covar(rr,:) = nan(1,length(tlags));
                EP_pairs(cc,bbb).mod_psth_covar(rr,:) = nan(1,length(tlags));
                Y1 = squeeze(allY1(:,cur_trial_set,Cpairs(cc,1)));
                for ll = 1:length(tlags)
                    Y2 = squeeze(allY2(:,:,Cpairs(cc,2),ll));
                    
                    EP_pairs(cc,bbb).mod_tot_covar(rr,ll) = nanmean(Y1(:).*Y2(:));
                    EP_pairs(cc,bbb).mod_psth_covar(rr,ll) = nanmean(nanmean(Y1,2).*nanmean(Y2,2));
                end
            end
        end
    end
    
end

%%
anal_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
if ~exist(anal_dir,'dir')
    mkdir(anal_dir)
end
cd(anal_dir);

sname = [sname sprintf('_ori%d',bar_ori)];
if rec_number > 1
    sname = strcat(sname,sprintf('_r%d',rec_number));
end

save(sname,'targs','EP_data','EP_pairs','EP_params','use_MUA','tlags');