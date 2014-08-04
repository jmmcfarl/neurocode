clear all
% close all
addpath(genpath('~/James_scripts'));
cd ~/Data/bruce/2_27_12/
load Blocks.mat
accept_window = [-5 5;-5 5];
sac_eyespeed = 10;
thresh_eyespeed = 2.5;
sac_buffer = 0;
min_fix_dur = 0.2;
blink_thresh = 5;
min_blink_dur = 0.05;
% nlags = 4;

Pix2Deg = 0.018837;
dsfrac = 2;
Fsd = 1/Pix2Deg/dsfrac;

% desired temporal resolution
stimres = 0.025; %in s


spk_cnts = [];
mua_cnts = [];
all_sac_amps = [];
all_sac_vels = [];
all_sac_durs = [];
all_fix_start_times = [];
all_fix_stop_times = [];
all_stim_filtered = [];
all_sac_locs_left = [];
all_sac_locs_right = [];
all_stim_nums = [];
blockids = [];

%%
cd /Users/James/Data/bruce/2_27_12/stimrecon

load fixation_data
load fixation_image_patches_d2 used*
  
cd /Users/James/Data/bruce/2_27_12/stimrecon
close all
load left_eye_pgabortrack_d2v1
left_X_mean = cur_X_seq(end,:)/Fsd;
left_X_std = post_stdX(end,:)/Fsd;
left_Y_mean = cur_Y_seq(end,:)/Fsd;
left_Y_std = post_stdY(end,:)/Fsd;
left_x_pos = gabor_params{4}(:,1);
left_y_pos = gabor_params{4}(:,2);

load right_eye_pgabortrack_d2v1
n_iter = size(cur_X_seq,1);
right_X_mean = cur_X_seq(end,:)/Fsd;
right_X_std = post_stdX(end,:)/Fsd;
right_Y_mean = cur_Y_seq(end,:)/Fsd;
right_Y_std = post_stdY(end,:)/Fsd;
right_x_pos = gabor_params{n_iter}(:,1);
right_y_pos = gabor_params{n_iter}(:,2);

avg_x_pos = 0.5*left_x_pos + 0.5*right_x_pos;
left_xoffset = mean(left_x_pos - avg_x_pos);
right_xoffset = mean(right_x_pos - avg_x_pos);
left_X_mean = left_X_mean + left_xoffset/Fsd;
right_X_mean = right_X_mean + right_xoffset/Fsd;

avg_y_pos = 0.5*left_y_pos + 0.5*right_y_pos;
left_yoffset = mean(left_y_pos - avg_y_pos);
right_yoffset = mean(right_y_pos - avg_y_pos);
left_Y_mean = left_Y_mean + left_yoffset/Fsd;
right_Y_mean = right_Y_mean + right_yoffset/Fsd;

cor_left(:,1) = all_sac_locs_left(used_fixs,1)+left_X_mean';
cor_left(:,2) = all_sac_locs_left(used_fixs,2)+left_Y_mean';
cor_right(:,1) = all_sac_locs_right(used_fixs,1)+right_X_mean';
cor_right(:,2) = all_sac_locs_right(used_fixs,2)+right_Y_mean';

%%
cd ~/Data/bruce/2_27_12/saccades/
load(sprintf('lemM232.5%d.em.sac.mat',1))
Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
EyeFs = 1/Eyedt;

maxlag = round(EyeFs*0.5);
ov_avg_hdisp = [];
ov_avg_vdisp = [];
for blockid = 1:4
    
    cd ~/Data/bruce/2_27_12/saccades/
    load(sprintf('lemM232.5%d.em.sac.mat',blockid))
    
    % identify saccade start and stop times
    EyeStartT = Expt.Trials.Start/10000; % time of first eye sample
    EyeEndT = Expt.Trials.End/10000; % time of last eye sample
    Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
    EyeFs = 1/Eyedt;
    eyets = EyeStartT:Eyedt:EyeEndT; %eye tracking time axis (sec)
    sac_buffer_inds = round(sac_buffer/Eyedt);
    
    reye_pos = [Expt.Trials.Eyevals.rh Expt.Trials.Eyevals.rv];
    leye_pos = [Expt.Trials.Eyevals.lh Expt.Trials.Eyevals.lv];    
        
    disparity = reye_pos - leye_pos;
    
    cur_fix_set = find(blockids(used_fixs) == blockid);
    cur_fix_start_times = all_fix_start_times(used_fixs(cur_fix_set));
    cur_fix_stop_times = all_fix_stop_times(used_fixs(cur_fix_set));
    cur_fix_lefterror = [left_X_mean(cur_fix_set)' left_Y_mean(cur_fix_set)'];
    cur_fix_righterror = [right_X_mean(cur_fix_set)' right_Y_mean(cur_fix_set)'];
    
    cur_fix_start_times2 = [EyeStartT; cur_fix_start_times; EyeEndT];
    cur_fix_lefterror = [cur_fix_lefterror(1,:); cur_fix_lefterror; cur_fix_lefterror(end,:)];
    cur_fix_righterror = [cur_fix_righterror(1,:); cur_fix_righterror; cur_fix_righterror(end,:)];
    
    interp_lefterror = interp1(cur_fix_start_times2,cur_fix_lefterror,eyets);
    interp_righterror = interp1(cur_fix_start_times2,cur_fix_righterror,eyets);
    
    inferred_leftpos = leye_pos + interp_lefterror(1:size(leye_pos,1),:);
    inferred_rightpos = reye_pos + interp_righterror(1:size(reye_pos,1),:);
    
    inferred_disp = inferred_rightpos - inferred_leftpos;
    
    cur_fix_start_inds = round(interp1(eyets,1:length(eyets),cur_fix_start_times));
    cur_fix_stop_inds = round(interp1(eyets,1:length(eyets),cur_fix_stop_times));
    bad_fixs = find(isnan(cur_fix_start_inds) | isnan(cur_fix_stop_inds));
    cur_fix_start_inds(bad_fixs) = [];
    cur_fix_stop_inds(bad_fixs) = [];
    avg_hdisp = nan(length(cur_fix_start_inds),maxlag+1);
    avg_vdisp = nan(length(cur_fix_start_inds),maxlag+1);
    for i = 1:length(cur_fix_start_inds)
        cur_set = cur_fix_start_inds(i):cur_fix_stop_inds(i);
        if length(cur_set) > maxlag+1
            cur_set((maxlag+1):end) = [];
        end
        cur_len = length(cur_set);
        avg_hdisp(i,1:cur_len) = inferred_disp(cur_set,1);      
        avg_vdisp(i,1:cur_len) = inferred_disp(cur_set,2);      
    end
     
    ov_avg_hdisp = [ov_avg_hdisp; avg_hdisp];
    ov_avg_vdisp = [ov_avg_vdisp; avg_vdisp];
    
end

lags = (0:maxlag)/EyeFs;
