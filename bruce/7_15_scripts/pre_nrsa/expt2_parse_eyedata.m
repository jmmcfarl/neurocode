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

dt = 118*2/1e4;
fst = 1/dt;

Pix2Deg = 0.018837;
dsfrac = 1.25;
Fsd = 1/Pix2Deg/dsfrac;
Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
use_win = [-5 5;-5 5]; 
min_trial_dur = 0.25;

%%
% Expt_nu = [3 8 15 19 26 29];
Expt_nu = [13 14 15 16 25 28 29];
n_allunits = 96;

all_expt_vec = [];
all_trial_vec = [];
all_trial_durs = [];
all_eyepos = [];
all_insac = [];
all_inblink = [];
all_inout = [];
all_t = [];
all_eyespeed = [];
all_imnum = [];

all_blink_times = [];
all_sac_times = [];
all_out_times = [];

all_sacdx = [];
all_sacdy = [];
all_sacamp = [];
for ee = 1:length(Expt_nu)
    
    fprintf('Analyzing Expt %d of %d\n',ee,length(Expt_nu));
    
    Trial_starts = [Expts{Expt_nu(ee)}.Trials(:).Start]/1e4;
    Trial_ends = [Expts{Expt_nu(ee)}.Trials(:).End]/1e4;
    endEvents = [Expts{Expt_nu(ee)}.Trials(:).endevent];
    Imnums = [Expts{Expt_nu(ee)}.Trials(:).IB];
    Trial_durs = (Trial_ends-Trial_starts);
    used_trials = find(Trial_durs > min_trial_dur);
    
    [eye_vals,eye_ts,eye_insample] = get_eye_tracking_data(em_data,[Trial_starts(1) Trial_ends(end)]);
    
    clear lcorrected_* corrected*
%     %correct eye positions
%     lcorrected_left = bsxfun(@plus,eye_vals(:,1:2)*left_gain,left_offset);
%     lcorrected_right = bsxfun(@plus,eye_vals(:,3:4)*right_gain,right_offset);
    %correct eye positions (quadratic)
    X = [eye_vals(:,1:2) eye_vals(:,1:2).^2 eye_vals(:,1).*eye_vals(:,2)];
    corrected_left(:,1) = bsxfun(@plus,X*bx_left(2:end),bx_left(1));
    corrected_left(:,2) = bsxfun(@plus,X*by_left(2:end),by_left(1));
    X = [eye_vals(:,3:4) eye_vals(:,3:4).^2 eye_vals(:,3).*eye_vals(:,4)];
    corrected_right(:,1) = bsxfun(@plus,X*bx_right(2:end),bx_right(1));
    corrected_right(:,2) = bsxfun(@plus,X*by_right(2:end),by_right(1));
    
%     avg_eyepos = 0.5*eye_vals(:,1:2) + 0.5*eye_vals(:,3:4);
%     avg_eyepos = 0.5*corrected_left + 0.5*corrected_right;
    avg_eyepos = corrected_left;
    out_window = find(avg_eyepos(:,1) < use_win(1,1) | avg_eyepos(:,1) > use_win(1,2) | ...
        avg_eyepos(:,2) < use_win(2,1) | avg_eyepos(:,2) > use_win(2,2));
    out_log = zeros(size(eye_ts));
    out_log(out_window) = 1;
    out_starts = find(out_log(1:end-1) == 0 & out_log(2:end) == 1)+1;
    out_ends = find(out_log(1:end-1)==1 & out_log(2:end)==0)+1;
    if out_log(1)==1
        out_starts = [1 out_starts];
    end
    if out_log(end) == 1
        out_ends = [out_ends length(eye_ts)];
    end
    
    [blink_data,in_blink,tot_disp_f] = get_blinks(eye_vals(:,3:4),eye_vals(:,1:2),eye_ts);
    temp_use = in_blink;
    temp_use(out_window) = 1;
    
%     [sac_data,in_sac,eye_speed] = get_saccades(eye_vals(:,3:4),eye_vals(:,1:2),eye_ts,in_blink);
%     [lsac_data,lin_sac,leye_speed] = get_saccades(lcorrected_right,lcorrected_left,eye_ts,in_blink);
    [sac_data,in_sac,eye_speed] = get_saccades(corrected_right,corrected_left,eye_ts,temp_use);
%     [sac_data,in_sac,eye_speed] = get_saccades_dualdet(corrected_right,corrected_left,eye_ts,in_blink);
    all_sacdx = [all_sacdx; [sac_data(:).dx_left]' [sac_data(:).dx_right]'];
    all_sacdy = [all_sacdy; [sac_data(:).dy_left]' [sac_data(:).dy_right]'];
    all_sacamp = [all_sacamp; [sac_data(:).amplitude]'];
    
    
    blink_start_times = [blink_data(:).start_times];
    blink_end_times = [blink_data(:).stop_times];
    sac_amps = [sac_data(:).amplitude];
    sac_start_times = [sac_data(:).start_time];
    sac_end_times = [sac_data(:).stop_time];
    out_start_times = eye_ts(out_starts);
    out_end_times = eye_ts(out_ends);
    
    
    cur_imnum = nan(size(eye_ts));
    cur_exptnum = ones(size(eye_ts))*ee;
    cur_trialnum = nan(size(eye_ts));
    use_trials = find(Trial_durs > min_trial_dur);
    for i = 1:length(use_trials)
       cur_set = find(eye_ts >= Trial_starts(use_trials(i)) & eye_ts <= Trial_ends(use_trials(i)));
       cur_imnum(cur_set) = Imnums(use_trials(i));
       cur_trialnum(cur_set) = use_trials(i);
    end
    
    all_t = [all_t; eye_ts'];
    all_eyepos = [all_eyepos; eye_vals];
    all_eyespeed = [all_eyespeed; eye_speed];
    all_inblink = [all_inblink; in_blink'];
    all_inout = [all_inout; out_log'];
    all_insac = [all_insac; in_sac'];
    all_imnum = [all_imnum; cur_imnum'];
    all_expt_vec = [all_expt_vec; cur_exptnum'];
    all_trial_vec = [all_trial_vec; cur_trialnum'];
    all_trial_durs = [all_trial_durs; Trial_durs(use_trials)'];
    
    all_blink_times = [all_blink_times; blink_start_times' blink_end_times'];
    all_sac_times = [all_sac_times; sac_start_times' sac_end_times'];
    all_out_times = [all_out_times; out_start_times' out_end_times'];
    
end

%%
save expt2_parsed_g034 all_*