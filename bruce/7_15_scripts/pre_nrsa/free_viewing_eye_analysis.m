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
use_win = [-5 5;-5 5]; 
min_trial_dur = 2;
% max_sac_amp = 0.5;
max_sac_amp = 10;

sac_dt = 0.01;
maxlag = round(1/sac_dt);
%%
% Expt_nu = [3 8 15 19 26 29];
Expt_nu = [13 14 15 16 25 28 29];
n_allunits = 96;
all_expt_vec = [];
all_trial_vec = [];
all_trial_durs = [];
all_eyepos = [];
all_ceyepos = [];
all_qeyepos = [];
all_insac = [];
all_inblink = [];
all_t = [];
all_eyespeed = [];

all_sac_start_times = [];
all_sac_end_times = [];
all_sac_amps = [];
all_sac_expt_vec = [];

pool_dx = [];
pool_dy = [];
pool_ldx = [];
pool_ldy = [];
pool_qdx = [];
pool_qdy = [];
pool_inblink = [];
for ee = 1:length(Expt_nu)
    
    fprintf('Analyzing Expt %d of %d\n',ee,length(Expt_nu));
%     single_units{ee} = find(CellList(Expt_nu(ee),:,1) > 0);
%     n_sus = length(single_units{ee});
%     multi_units{ee} = setdiff(1:n_allunits,single_units{ee});
%     n_mus = 96;
%     load(sprintf('Expt%dClusterTimes.mat',Expt_nu(ee)));
    
    Trial_starts = [Expts{Expt_nu(ee)}.Trials(:).Start]/1e4;
    Trial_ends = [Expts{Expt_nu(ee)}.Trials(:).End]/1e4;
    endEvents = [Expts{Expt_nu(ee)}.Trials(:).endevent];
    Trial_durs = (Trial_ends-Trial_starts);
    %     used_trials = 1:length(Trial_durs); %use all trials
    used_trials = find(Trial_durs > min_trial_dur);
    
    [eye_vals,eye_ts,eye_insample] = get_eye_tracking_data(em_data,[Trial_starts(1) Trial_ends(end)]);
    
    clear lcorrected_* corrected*
    %correct eye positions
    corrected_left = bsxfun(@plus,eye_vals(:,1:2)*left_gain,left_offset);
    corrected_right = bsxfun(@plus,eye_vals(:,3:4)*right_gain,right_offset);
%     %correct eye positions (quadratic)
%     X = [eye_vals(:,1:2) eye_vals(:,1:2).^2 eye_vals(:,1).*eye_vals(:,2)];
%     corrected_left(:,1) = bsxfun(@plus,X*bx_left(2:end),bx_left(1));
%     corrected_left(:,2) = bsxfun(@plus,X*by_left(2:end),by_left(1));
%     X = [eye_vals(:,3:4) eye_vals(:,3:4).^2 eye_vals(:,3).*eye_vals(:,4)];
%     corrected_right(:,1) = bsxfun(@plus,X*bx_right(2:end),bx_right(1));
%     corrected_right(:,2) = bsxfun(@plus,X*by_right(2:end),by_right(1));
    
%     avg_eyepos = 0.5*eye_vals(:,1:2) + 0.5*eye_vals(:,3:4);
    avg_eyepos = 0.5*corrected_left + 0.5*corrected_right;
    out_window = find(avg_eyepos(:,1) < use_win(1,1) | avg_eyepos(:,1) > use_win(1,2) | ...
        avg_eyepos(:,2) < use_win(2,1) | avg_eyepos(:,2) > use_win(2,2));
    
    [blink_data,in_blink,tot_disp_f] = get_blinks(eye_vals(:,3:4),eye_vals(:,1:2),eye_ts);
    in_blink(out_window) = 1;
    
%     [sac_data,in_sac,eye_speed] = get_saccades(eye_vals(:,3:4),eye_vals(:,1:2),eye_ts,in_blink);
%     [lsac_data,lin_sac,leye_speed] = get_saccades(lcorrected_right,lcorrected_left,eye_ts,in_blink);
    [sac_data,in_sac,eye_speed] = get_saccades(corrected_right,corrected_left,eye_ts,in_blink);
%     [sac_data,in_sac,eye_speed] = get_saccades_dualdet(corrected_right,corrected_left,eye_ts,in_blink);
    
    pool_inblink = [pool_inblink in_blink];
    pool_dx = [pool_dx; [sac_data(:).dx_left]' [sac_data(:).dx_right]'];
    pool_dy = [pool_dy; [sac_data(:).dy_left]' [sac_data(:).dy_right]'];
%     pool_ldx = [pool_ldx; [lsac_data(:).dx_left]' [lsac_data(:).dx_right]'];
%     pool_ldy = [pool_ldy; [lsac_data(:).dy_left]' [lsac_data(:).dy_right]'];
%     pool_qdx = [pool_qdx; [qsac_data(:).dx_left]' [qsac_data(:).dx_right]'];
%     pool_qdy = [pool_qdy; [qsac_data(:).dy_left]' [qsac_data(:).dy_right]'];
    
    
    blink_start_times = [blink_data(:).start_times];
    blink_stop_times = [blink_data(:).stop_times];
    sac_amps = [sac_data(:).amplitude];
    sac_start_times = [sac_data(:).start_time];
    sac_end_times = [sac_data(:).stop_time];
    big_sacs = find(sac_amps > max_sac_amp);
    big_sac_starts = sac_start_times(big_sacs);
    big_sac_ends = sac_end_times(big_sacs);
    
    all_sac_start_times = [all_sac_start_times; sac_start_times(:)];
    all_sac_end_times = [all_sac_end_times; sac_end_times(:)];
    all_sac_amps = [all_sac_amps; sac_amps(:)];
    all_sac_expt_vec = [all_sac_expt_vec; ee*ones(length(sac_amps),1)];
    
    
    macro_sacs = find(sac_amps > 1);
    micro_sacs = find(sac_amps < 1);
    sac_t_axis = Trial_starts(1):sac_dt:Trial_ends(end);
    sac_times_binned = histc(sac_start_times,sac_t_axis);
    micro_times_binned = histc(sac_start_times(micro_sacs),sac_t_axis);
    macro_times_binned = histc(sac_start_times(macro_sacs),sac_t_axis);
    
    all_acorr(ee,:) = xcov(sac_times_binned,maxlag,'coeff');
    macro_acorr(ee,:) = xcov(macro_times_binned,maxlag,'coeff');
     micro_acorr(ee,:) = xcov(micro_times_binned,maxlag,'coeff');
   
    
%     
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
%     
%     bad_trials = [];
%     for i = 1:length(used_trials)
%         if ~isempty(find(blink_start_times >= Trial_starts(used_trials(i)) & ...
%                 blink_start_times <= Trial_ends(used_trials(i)), 1)) ...
%                 | ~isempty(find(blink_stop_times >= Trial_starts(used_trials(i)) & ...
%                 blink_stop_times <= Trial_ends(used_trials(i)), 1))
%             bad_trials = [bad_trials i];
%         end
%     end
%     fprintf('Eliminating %d of %d blink trials\n',length(bad_trials),length(used_trials));
%     used_trials(bad_trials) = [];
%     
%     bad_trials = [];
%     for i = 1:length(used_trials)
%         if isfield(Expts{Expt_nu(ee)}.Trials(used_trials(i)),'rptframes')
%             if ~isempty(Expts{Expt_nu(ee)}.Trials(used_trials(i)).rptframes)
%                 bad_trials = [bad_trials; i];
%             end
%         end
%     end
%     fprintf('Eliminating %d of %d framerep trials\n',length(bad_trials),length(used_trials));
%     used_trials(bad_trials) = [];
    
    %     used_trials(used_trials==length(Trial_durs)) = []; %to prevent errors when testing misalignment
    %      used_trials(used_trials==1) = []; %to prevent errors when testing misalignment
%     use = find(in_blink==0 & in_sac==0);
    use = find(in_blink==0);
    all_t = [all_t; eye_ts(use)'];
    all_eyepos = [all_eyepos; eye_vals(use,:)];
    all_eyespeed = [all_eyespeed; eye_speed(use,:)];
%     all_ceyepos = [all_ceyepos; lcorrected_left(use,:) lcorrected_right(use,:)];    
%     all_qeyepos = [all_qeyepos; corrected_left(use,:) corrected_right(use,:)];  
    
end


%%
lags = -maxlag:maxlag;
lags = lags*sac_dt;
zl = find(lags == 0);
all_acorr(:,zl) = 0;
macro_acorr(:,zl) = 0;
micro_acorr(:,zl) = 0;
shadedErrorBar(lags,mean(all_acorr),std(all_acorr)/sqrt(7))
hold on
shadedErrorBar(lags,mean(macro_acorr),std(macro_acorr)/sqrt(7),{'color','b'})
shadedErrorBar(lags,mean(micro_acorr),std(micro_acorr)/sqrt(7),{'color','r'})


%%
% save expt2_eyeanalysis_g034 all_*