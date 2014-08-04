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
% Expt_nu = [1 6 16 17 20 25 28]; %expt 1
% Expt_nu = [3 8 15 19 26 29]; %expt 2
% Expt_nu = [4 5 9 23 32]; %expt 3
% Expt_nu = [14 22 27]; %expt 4
Expt_nu = [1 2 17 18 19 23 24]; %expt 1 34


n_allunits = 96;
all_expt_vec = [];
all_trial_vec = [];
all_trial_durs = [];
all_eyepos = [];
all_eyespeed = [];
% all_ceyepos = [];
% all_qeyepos = [];
all_eyed = [];
all_leyed = [];
all_qeyed = [];
all_insac = [];
all_inblink = [];
all_t = [];
all_sac_start_times = [];
all_sac_stop_times = [];
all_sac_amps = [];

pool_dx = [];
pool_dy = [];
pool_ldx = [];
pool_ldy = [];
pool_qdx = [];
pool_qdy = [];
pool_inblink = [];
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
    outrange = find(eye_ts > Trial_ends(end) | eye_ts < Trial_starts(1));
    eye_vals(outrange,:) = [];
    eye_ts(outrange) = [];
    eye_insample(outrange) = [];
    
    clear lcorrected_* corrected*
    %correct eye positions
    lcorrected_left = bsxfun(@plus,eye_vals(:,1:2)*left_gain,left_offset);
    lcorrected_right = bsxfun(@plus,eye_vals(:,3:4)*right_gain,right_offset);
    %correct eye positions (quadratic)
    X = [eye_vals(:,1:2) eye_vals(:,1:2).^2 eye_vals(:,1).*eye_vals(:,2)];
    corrected_left(:,1) = bsxfun(@plus,X*bx_left(2:end),bx_left(1));
    corrected_left(:,2) = bsxfun(@plus,X*by_left(2:end),by_left(1));
    X = [eye_vals(:,3:4) eye_vals(:,3:4).^2 eye_vals(:,3).*eye_vals(:,4)];
    corrected_right(:,1) = bsxfun(@plus,X*bx_right(2:end),bx_right(1));
    corrected_right(:,2) = bsxfun(@plus,X*by_right(2:end),by_right(1));
    
    leye_vals = [lcorrected_left lcorrected_right];
    qeye_vals = [corrected_left corrected_right];
    
    eye_fvals = filtfilt(bb,aa,eye_vals);
    leye_fvals = filtfilt(bb,aa,leye_vals);
    qeye_fvals = filtfilt(bb,aa,qeye_vals);
    
    avg_eyepos = 0.5*eye_vals(:,1:2) + 0.5*eye_vals(:,3:4);
    
    [blink_data,in_blink,tot_disp_f] = get_blinks(eye_vals(:,3:4),eye_vals(:,1:2),eye_ts);
    
    [sac_data,in_sac,eye_speed] = get_saccades(eye_vals(:,3:4),eye_vals(:,1:2),eye_ts,in_blink);
%     [lsac_data,lin_sac,leye_speed] = get_saccades(lcorrected_right,lcorrected_left,eye_ts,in_blink);
%     [qsac_data,qin_sac,qeye_speed] = get_saccades(corrected_right,corrected_left,eye_ts,in_blink);    
    
    blink_start_times = [blink_data(:).start_times];
    blink_stop_times = [blink_data(:).stop_times];
    sac_start_times = [sac_data(:).start_time];
    sac_end_times = [sac_data(:).stop_time];
%     lsac_start_times = [lsac_data(:).start_time];
%     lsac_end_times = [lsac_data(:).stop_time];
%     qsac_start_times = [qsac_data(:).start_time];
%     qsac_end_times = [qsac_data(:).stop_time];
    
    all_sac_start_times = [all_sac_start_times sac_start_times];
    all_sac_stop_times = [all_sac_stop_times sac_end_times];
    all_sac_amps = [all_sac_amps sac_data(:).amplitude];
% %     check for blinks or big sacs within each trial
%     bad_trials = [];
%     for i = 1:length(used_trials)
%         if ~isempty(find(big_sac_starts > Trial_starts(used_trials(i)) & ...
%                 big_sac_starts < Trial_ends(used_trials(i)), 1)) ...
%                  | ~isempty(find(big_sac_ends >= Trial_starts(used_trials(i)) & ...
%                 big_sac_ends <= Trial_ends(used_trials(i)), 1))
%            bad_trials = [bad_trials i];
%         end
%     end
%     fprintf('Eliminating %d of %d sac trials\n',length(bad_trials),length(used_trials));
%     used_trials(bad_trials) = [];
    
        bad_trials = [];
    for i = 1:length(used_trials)
        if ~isempty(find(blink_start_times >= Trial_starts(used_trials(i)) & ...
                blink_start_times <= Trial_ends(used_trials(i)), 1)) ...
                | ~isempty(find(blink_stop_times >= Trial_starts(used_trials(i)) & ...
                blink_stop_times <= Trial_ends(used_trials(i)), 1))
            bad_trials = [bad_trials i];
        end
    end
    fprintf('Eliminating %d of %d blink trials\n',length(bad_trials),length(used_trials));
    used_trials(bad_trials) = [];
    
    eye_d = [zeros(1,4); diff(eye_fvals)];
%     leye_d = [zeros(1,4); diff(leye_fvals)];
%     qeye_d = [zeros(1,4); diff(qeye_fvals)];
    
    used_inds = [];
    used_sacs = [];
    used_lsacs = [];
    used_qsacs = [];
    for i = 1:length(used_trials)
        cur_set = find(eye_ts > Trial_starts(used_trials(i)) & eye_ts < Trial_ends(used_trials(i)));
        used_inds = [used_inds cur_set];
        
        all_eyed = [all_eyed; eye_d(cur_set,:)];
%         all_leyed = [all_leyed; leye_d(cur_set,:)];
%         all_qeyed = [all_qeyed; qeye_d(cur_set,:)];
        
        cur_sac_set = find(sac_start_times > Trial_starts(used_trials(i)) & sac_end_times < Trial_ends(used_trials(i)));
        used_sacs = [used_sacs cur_sac_set];
%         cur_sac_set = find(lsac_start_times > Trial_starts(used_trials(i)) & lsac_end_times < Trial_ends(used_trials(i)));
%         used_lsacs = [used_lsacs cur_sac_set];
%         cur_sac_set = find(qsac_start_times > Trial_starts(used_trials(i)) & qsac_end_times < Trial_ends(used_trials(i)));
%         used_qsacs = [used_qsacs cur_sac_set];
    end
    
%     used_inds(in_blink(used_inds)==1 | in_sac(used_inds)==1) = [];
    used_inds = 1:length(eye_ts);
    
    pool_inblink = [pool_inblink in_blink];
    pool_dx = [pool_dx; [sac_data(used_sacs).dx_left]' [sac_data(used_sacs).dx_right]'];
    pool_dy = [pool_dy; [sac_data(used_sacs).dy_left]' [sac_data(used_sacs).dy_right]'];
%     pool_ldx = [pool_ldx; [lsac_data(used_lsacs).dx_left]' [lsac_data(used_lsacs).dx_right]'];
%     pool_ldy = [pool_ldy; [lsac_data(used_lsacs).dy_left]' [lsac_data(used_lsacs).dy_right]'];
%     pool_qdx = [pool_qdx; [qsac_data(used_qsacs).dx_left]' [qsac_data(used_qsacs).dx_right]'];
%     pool_qdy = [pool_qdy; [qsac_data(used_qsacs).dy_left]' [qsac_data(used_qsacs).dy_right]'];
    
    
    all_expt_vec = [all_expt_vec ee*ones(size(used_inds))];
    
    %     used_trials(used_trials==length(Trial_durs)) = []; %to prevent errors when testing misalignment
    %      used_trials(used_trials==1) = []; %to prevent errors when testing misalignment
%     use = find(in_blink==0 & in_sac==0);
    all_t = [all_t; eye_ts(used_inds)'];
    all_eyepos = [all_eyepos; eye_vals(used_inds,:)];
    all_eyespeed = [all_eyespeed; eye_speed(used_inds,:)];
%     all_ceyepos = [all_ceyepos; lcorrected_left(used_inds,:) lcorrected_right(used_inds,:)];
%     all_qeyepos = [all_qeyepos; corrected_left(used_inds,:) corrected_right(used_inds,:)];
    
%     figure
%     plot(eye_ts-eye_ts(1),leye_fvals(:,2)-mean(leye_fvals(:,2)))
%     hold on
%     plot(eye_ts-eye_ts(1),leye_fvals(:,4)-mean(leye_fvals(:,4)),'r')
%     plot(eye_ts-eye_ts(1),in_blink/3-0.4,'k','linewidth',2)
%     plot(eye_ts-eye_ts(1),in_sac/3-0.4,'g')
%     ylim([-0.7 0.7])
% % 
%         figure
%     plot(eye_ts-eye_ts(1),leye_fvals(:,1)-mean(leye_fvals(:,1)))
%     hold on
%     plot(eye_ts-eye_ts(1),leye_fvals(:,3)-mean(leye_fvals(:,3)),'r')
%     plot(eye_ts-eye_ts(1),in_blink/3-0.4,'k','linewidth',2)
%     plot(eye_ts-eye_ts(1),in_sac/3-0.4,'g')
%     ylim([-0.7 0.7])
% 
%     pause
    
end

%%
save expt1_eyedata_34_new all*






