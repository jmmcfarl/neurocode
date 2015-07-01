clear all

Expt_name = 'M012';
monk_name = 'jbe';
bar_ori = 0; %bar orientation to use (only for UA recs)
rec_number = 1;

old_anal_name = 'full_eyetrack_FIN';
mod_name = 'corrected_models_comp';

%%
Expt_num = str2num(Expt_name(2:end));

data_dir = ['~/Data/bruce/' Expt_name];
if ~exist(data_dir,'dir')
    data_dir = ['/media/NTlab_data3/Data/bruce/' Expt_name];
end
if ~exist(data_dir,'dir');
    error('Couldnt find data directory');
end
Edata_file = strcat(data_dir,'/',monk_name,Expt_name,'Expts');
load(Edata_file);

data_name = sprintf('%s/packaged_data_ori%d',data_dir,bar_ori);
if rec_number > 1
    data_name = strcat(data_name,sprintf('_r%d',rec_number));
end
fprintf('Loading %s\n',data_name);
load(data_name);

anal_dir = ['~/Analysis/bruce/' Expt_name '/ET_final_imp/'];
model_dir = ['~/Analysis/bruce/' Expt_name '/models'];
if ~exist(anal_dir,'dir')
    system(['mkdir ' anal_dir]);
end

old_anal_name = [old_anal_name '_Rinit'];

%if using coil info
if any(params.use_coils > 0)
    old_anal_name = [old_anal_name '_Cprior'];
end

hr_anal_name = [old_anal_name '_hres'];

hr_anal_name = [hr_anal_name sprintf('_ori%d',bar_ori)];
mod_name = [mod_name sprintf('_ori%d',bar_ori)];

if rec_number > 1
    hr_anal_name = strcat(hr_anal_name,sprintf('r%d',rec_number));
    mod_name = strcat(mod_name,sprintf('_r%d',rec_number));
end

%% get trial and block data
dt = params.dt;
flen = modFitParams.flen;
NT = length(used_inds);
fullNT = size(spike_data.binned_mua,1);
n_trials = length(time_data.trial_flip_ids);
n_blocks = length(expt_data.used_blocks);

all_t_axis = time_data.t_axis;
trial_start_inds = [1; 1+time_data.trial_flip_inds(2:end)];
trial_end_inds = [time_data.trial_flip_inds(2:end); fullNT];
all_trialvec = nan(fullNT,1);
for ii = 1:n_trials
    all_trialvec(trial_start_inds(ii):trial_end_inds(ii)) = time_data.trial_flip_ids(ii);
end

%% get needed saccade data
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

%% DEFINE FIXATION POINTS
trial_start_inds = [1; find(diff(all_trialvec(used_inds)) ~= 0) + 1];
trial_end_inds = [find(diff(all_trialvec(used_inds)) ~= 0); NT];

fix_start_inds = sort([trial_start_inds; saccade_stop_inds]);
fix_stop_inds = sort([trial_end_inds; saccade_start_inds]);
fix_durs = fix_stop_inds-fix_start_inds;
fix_start_inds(fix_durs<=0) = [];
fix_stop_inds(fix_durs <= 0) = [];
fix_post_blink = ismember(fix_start_inds,saccade_stop_inds(used_is_blink));
n_fixs = length(fix_start_inds);

%push the effects of saccades forward in time
sac_shift = round(0.05/dt);
fix_ids = nan(NT,1);
for ii = 1:n_fixs
    cur_inds = fix_start_inds(ii):(fix_stop_inds(ii));
    fix_ids(cur_inds) = ii;
end

%% load in and process HRES eye tracking data
cd(anal_dir)
load(hr_anal_name,'best_fix_cor','best_fix_std','drift_post_mean','drift_post_std','et_params');

[best_tot_corr,best_tot_std] = construct_eye_position(best_fix_cor,best_fix_std,...
    drift_post_mean,drift_post_std,fix_ids,trial_start_inds,trial_end_inds,sac_shift);

ep_in_deg = best_tot_corr*et_params.sp_dx;

%% identify stimulus repeat trials
if ~isnan(params.rpt_seeds)
    xv_type = 'rpt';
else
    xv_type = 'uni';
end

if strcmp(xv_type,'rpt')
    rpt_trials = find(ismember([trial_data(:).se],params.rpt_seeds));
else
    rpt_trials = [];
end

n_rpts = length(rpt_trials);
used_rpt_inds = find(ismember(all_trialvec(used_inds),rpt_trials));
nf = 400 - round((params.end_buffer + params.beg_buffer)/dt); %number of used time samples per trial

%% LOAD IN MODEL FITS
load([model_dir '/' mod_name]);

%% process stimulus
all_stim_mat = decompressTernNoise(stimComp);
full_nPix_us = et_params.spatial_usfac*params.full_nPix;
if et_params.spatial_usfac > 1
    all_stimmat_up = zeros(size(all_stim_mat,1),full_nPix_us);
    for ii = 1:size(all_stim_mat,2)
        for jj = 1:et_params.spatial_usfac
            all_stimmat_up(:,et_params.spatial_usfac*(ii-1)+jj) = all_stim_mat(:,ii);
        end
    end
elseif et_params.spatial_usfac == 1
    all_stimmat_up = all_stim_mat;
end
stim_params_us = NMMcreate_stim_params([flen full_nPix_us],dt);

%repeat for up-sampled versions of the Xmatrix
use_nPix = modFitParams.use_nPix_us/modFitParams.spatial_usfac;
buffer_pix = floor((params.full_nPix - use_nPix)/2);
[Xinds_up,~] = meshgrid(1/et_params.spatial_usfac:1/et_params.spatial_usfac:params.full_nPix,1:flen);
cur_use_pix = (1/et_params.spatial_usfac:1/et_params.spatial_usfac:use_nPix) + buffer_pix;
use_kInds_up = find(ismember(Xinds_up(:),cur_use_pix));

%%


