clear all
close all
cd
addpath('~/Data/bruce/7_15_12/')

cd ~/Data/bruce/7_15_12

cd G029/
% load ./jbeG029.em.mat
% em_data = Expt; clear Expt
% load ./CellList.mat
% load ./G029Expts.mat
load ./eye_calibration_data
cd ..

cd G034
load ./jbeG034.em.mat
em_data = Expt; clear Expt
load ./CellList.mat
load ./G034Expts.mat
% load ./eye_calibration_data

% cd G035
% load ./jbeG035.em.mat
% em_data = Expt; clear Expt
% load ./CellList.mat
% load ./G035Expts.mat
% load ./eye_calibration_data


dt = 118*2/1e4;
fst = 1/dt;

Pix2Deg = 0.018837;
dsfrac = 1.25;
Fsd = 1/Pix2Deg/dsfrac;
Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);

min_trial_dur = 2;
max_sac_amp = 0.5;

eye_fs = 588.24;
hcf = 100;
[bb,aa] = butter(2,hcf/(eye_fs/2),'low');

%%
% Expt_nu = [1 6 16 17 20 25 28]; %expt 1
% Expt_nu = [3 8 15 19 26 29]; %expt 2
% Expt_nu = [4 5 9 23 32]; %expt 3
% Expt_nu = [7:12]; %expt 3 34
% Expt_nu = [1:5 8]; %expt 4 35
Expt_nu = [21 22]; %expt 4 34
% Expt_nu = [14 22 27]; %expt 4
% Expt_nu = [1 2 17 18 19 23 24]; %expt 1 34


n_allunits = 96;
all_expt_vec = [];
all_trial_vec = [];
all_trial_durs = [];
all_eyepos = [];
all_inblink = [];
% all_inhf = [];
all_insac = [];
all_t = [];
all_sac_start_inds = [];
all_sac_stop_inds = [];
all_sac_amps = [];
% all_lsac_stop_inds = [];
% all_rsac_stop_inds = [];
% all_lsac_camps = [];
% all_rsac_camps = [];
% all_lsac_amps = [];
% all_rsac_amps = [];
% all_hf_start_inds = [];
% all_hf_stop_inds = [];
all_blink_start_inds = [];
all_blink_stop_inds = [];
all_eyespeed = [];

for ee = 1:length(Expt_nu)
    
    fprintf('Analyzing Expt %d of %d\n',ee,length(Expt_nu));
    single_units{ee} = find(CellList(Expt_nu(ee),:,1) > 0);
    n_sus = length(single_units{ee});
    multi_units{ee} = setdiff(1:n_allunits,single_units{ee});
    n_mus = 96;
    load(sprintf('Expt%dClusterTimes.mat',Expt_nu(ee)));
    
    Trial_starts = [Expts{Expt_nu(ee)}.Trials(:).Start]/1e4;
    Trial_ends = [Expts{Expt_nu(ee)}.Trials(:).End]/1e4;
    endEvents = [Expts{Expt_nu(ee)}.Trials(:).endevent];
    Trial_durs = (Trial_ends-Trial_starts);
    %     used_trials = 1:length(Trial_durs); %use all trials
    used_trials = find(Trial_durs > min_trial_dur);
    
    [eye_vals,eye_ts,eye_insample] = get_eye_tracking_data(em_data,[Trial_starts(1) Trial_ends(end)]);
    bad_samps = find(eye_ts > Trial_ends(end) | eye_ts < Trial_starts(1));
    eye_vals(bad_samps,:) = [];
    eye_ts(bad_samps) = [];
    eye_insample(bad_samps) = [];
    
    clear lcorrected_* corrected*
    %correct eye positions
    corrected_left = bsxfun(@plus,eye_vals(:,1:2)*left_gain,left_offset);
    corrected_right = bsxfun(@plus,eye_vals(:,3:4)*right_gain,right_offset);
    eye_vals = [corrected_left corrected_right];
    eye_fvals = filtfilt(bb,aa,eye_vals);
    avg_eyepos = 0.5*eye_fvals(:,1:2) + 0.5*eye_fvals(:,3:4);
    
    leye_vel = [0 0; diff(eye_fvals(:,1:2))];
    reye_vel = [0 0; diff(eye_fvals(:,3:4))];
    leye_speed = sqrt(sum(leye_vel.^2,2));
    reye_speed = sqrt(sum(reye_vel.^2,2));
    eye_speed = [leye_speed reye_speed];
    
    [blink_data,in_blink,tot_disp_f] = get_blinks_v2(eye_vals(:,3:4),eye_vals(:,1:2),eye_ts);
    
%     [hf_start_inds,hf_stop_inds,sac_data] = get_saccades_dualdet(eye_fvals(:,3:4),eye_fvals(:,1:2),eye_ts,in_blink);
    [sac_data,in_sac,eye_speed] = get_saccades(eye_fvals(:,3:4),eye_fvals(:,1:2),eye_ts,in_blink);
    sac_amps = [sac_data(:).amplitude];
    sac_start_inds = sac_data.start_inds;
    sac_stop_inds = sac_data.stop_inds;
    all_sac_start_inds = [all_sac_start_inds; length(all_t) + sac_start_inds];
    all_sac_stop_inds = [all_sac_stop_inds; length(all_t) + sac_stop_inds];
%     all_lsac_stop_inds = [all_lsac_stop_inds; length(all_t) + sac_data.left_stop_inds];
%     all_rsac_stop_inds = [all_lsac_stop_inds; length(all_t) + sac_data.right_stop_inds];
%     all_lsac_camps = [all_lsac_camps; sac_data.left_camps];
%     all_rsac_camps = [all_rsac_camps; sac_data.right_camps];
%     all_lsac_amps = [all_lsac_amps; sac_data.left_amps];
%     all_rsac_amps = [all_rsac_amps; sac_data.right_amps];
    
%     all_hf_start_inds = [all_hf_start_inds; length(all_t) + hf_start_inds];
%     all_hf_stop_inds = [all_hf_stop_inds; length(all_t) + hf_stop_inds];
    
    blink_start_inds = 1+find(in_blink(1:end-1)==0 & in_blink(2:end)==1);
    blink_stop_inds = 1+find(in_blink(1:end-1)==1 & in_blink(2:end)==0);
    all_blink_start_inds = [all_blink_start_inds length(all_t) + blink_start_inds];
    all_blink_stop_inds = [all_blink_stop_inds length(all_t) + blink_stop_inds];
    
%     cur_in_hf = zeros(size(in_blink));
%     for i = 1:length(hf_start_inds)
%         cur_in_hf(hf_start_inds(i):hf_stop_inds(i)) = 1;
%     end
    cur_in_sac = zeros(size(in_blink));
    for i = 1:length(sac_start_inds)
        cur_in_sac(sac_start_inds(i):sac_stop_inds(i)) = 1;
    end
    
    all_t = [all_t; eye_ts'];
    all_sac_amps = [all_sac_amps; sac_amps'];
    all_inblink = [all_inblink; in_blink'];
%     all_inhf = [all_inhf; cur_in_hf'];
    all_insac = [all_insac; cur_in_sac'];
    all_eyepos = [all_eyepos; eye_fvals];
    all_eyespeed = [all_eyespeed; eye_speed];
    
end

%%
cd ~/Data/bruce/7_15_12/G034
save all_eyedata_expt4_34 all*

%%
