clear all
% close all
cd ~/Data/bruce/G075/
load jbeG075Expts.mat
addpath('~/Data/bruce/7_15_12/')


load ./jbeG075.em.mat
em_data = Expt; clear Expt
load ./CellList.mat

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
hcf = 80;
[bb,aa] = butter(2,hcf/(eye_fs/2),'low');

%%
% sim_sac_blocks = [10 14 18 22 24 28 32 37 40 42 50] - 6;
sim_sac_blocks = [14 18 24 28 37 40 50] - 6;
% sim_sac_blocks = [18 28 40] - 6; %sim sac
% sim_sac_blocks = [10 22 32 42] - 6; %sim sac a
% sim_sac_blocks = [14 24 37 50] - 6; %sim sac b

n_allunits = 96;
all_expt_vec = [];
all_trial_vec = [];
all_trial_durs = [];
all_eyepos = [];
all_inblink = [];
all_inhf = [];
all_insac = [];
all_t = [];
all_sac_start_inds = [];
all_sac_stop_inds = [];
all_hf_start_inds = [];
all_hf_stop_inds = [];
all_blink_start_inds = [];
all_blink_stop_inds = [];
all_eyespeed = [];
all_sac_amps = [];
for ee = 1:length(sim_sac_blocks)
    
    fprintf('Analyzing Expt %d of %d\n',ee,length(sim_sac_blocks));
    single_units{ee} = find(CellList(sim_sac_blocks(ee),:,1) > 0);
    n_sus = length(single_units{ee});
    multi_units{ee} = setdiff(1:n_allunits,single_units{ee});
    n_mus = 96;
    load(sprintf('Expt%dClusterTimes.mat',sim_sac_blocks(ee)));
    
        cur_n_trials = length(Expts{sim_sac_blocks(ee)}.Trials);
    bad_t = [];
    for i = 1:cur_n_trials
        if length(Expts{sim_sac_blocks(ee)}.Trials(i).Start) ~= 1
            bad_t = [bad_t i];
        end
    end
    Expts{sim_sac_blocks(ee)}.Trials(bad_t) = [];

    Trial_starts = [Expts{sim_sac_blocks(ee)}.Trials(:).Start]/1e4;
    Trial_ends = [Expts{sim_sac_blocks(ee)}.Trials(:).End]/1e4;
%     endEvents = [Expts{sim_sac_blocks(ee)}.Trials(:).endevent];
    Trial_durs = (Trial_ends-Trial_starts);
    %     used_trials = 1:length(Trial_durs); %use all trials
    used_trials = find(Trial_durs > min_trial_dur);
    
    [eye_vals,eye_ts,eye_insample] = get_eye_tracking_data_new(em_data,[Trial_starts(1) Trial_ends(end)]);
    use_samps = find(eye_ts >= Trial_starts(1) & eye_ts <= Trial_ends(end));
    eye_vals = eye_vals(use_samps,:);
    eye_insample = eye_insample(use_samps);
    eye_ts = eye_ts(use_samps);
    
    
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
    sac_start_inds = sac_data.start_inds;
    sac_stop_inds = sac_data.stop_inds;
    all_sac_start_inds = [all_sac_start_inds; length(all_t) + sac_start_inds];
    all_sac_stop_inds = [all_sac_stop_inds; length(all_t) + sac_stop_inds];
    all_sac_amps = [all_sac_amps sac_data(:).amplitude];
    
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
    all_inblink = [all_inblink; in_blink'];
%     all_inhf = [all_inhf; cur_in_hf'];
    all_insac = [all_insac; cur_in_sac'];
    all_eyepos = [all_eyepos; eye_fvals];
    all_eyespeed = [all_eyespeed; eye_speed];
    
end

%%
cd ~/Data/bruce/G075/

save expt3_eyedata_NS_sensitive all*






