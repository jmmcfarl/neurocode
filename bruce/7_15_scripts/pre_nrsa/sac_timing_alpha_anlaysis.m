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

sac_dt = 0.005;
maxlag = round(1/sac_dt);
sac_dp = 0.3;
maxlag_p = round(10*2*pi/sac_dp);



Fs = 3e4;
dsf = 100;Fsd = Fs/dsf;
use_lfps = [6];
[b_alpha,a_alpha] = butter(2,[5 12]/(Fsd/2));
alpha_ch = 1;
alpha_dp = 0.15;

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
all_t = [];
all_eyespeed = [];

all_sac_start_times = [];
all_sac_end_times = [];
all_sac_amps = [];
all_sac_expt_vec = [];

for ee = 1:length(Expt_nu)
    
    fprintf('Analyzing Expt %d of %d\n',ee,length(Expt_nu));
    
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

    avg_eyepos = 0.5*corrected_left + 0.5*corrected_right;
    out_window = find(avg_eyepos(:,1) < use_win(1,1) | avg_eyepos(:,1) > use_win(1,2) | ...
        avg_eyepos(:,2) < use_win(2,1) | avg_eyepos(:,2) > use_win(2,2));
    
    [blink_data,in_blink,tot_disp_f] = get_blinks(eye_vals(:,3:4),eye_vals(:,1:2),eye_ts);
    in_blink(out_window) = 1;
    
    [sac_data,in_sac,eye_speed] = get_saccades(corrected_right,corrected_left,eye_ts,in_blink);    
    
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
    
    
    clear all_Valpha
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d\n',ll);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = decimate(V,dsf);
        V = V(:);
        V_alpha = filtfilt(b_alpha,a_alpha,V);
%         V_alpha = filter(b_alpha,a_alpha,V);
        all_Valpha(:,ll) = V_alpha;
        alpha_phase = angle(hilbert(V_alpha));
        alpha_phase = unwrap_phase_monotonic(alpha_phase);
    end
%         avg_alpha = mean(all_Valpha,2);
    avg_alpha = all_Valpha(:,alpha_ch);
    avg_alpha_phase = angle(hilbert(avg_alpha));
    avg_alpha_phase = unwrap_phase_monotonic(avg_alpha_phase);
    avg_alpha_phase = smooth(avg_alpha_phase,15,'lowess');
   
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);

    uniform_phase_grid = avg_alpha_phase(1):sac_dp:avg_alpha_phase(end);
    t_ax_interp = interp1(avg_alpha_phase,t_ax,uniform_phase_grid);

    sac_start_phases = interp1(t_ax,avg_alpha_phase,sac_start_times);
    
%     sac_start_phase_bins = round(interp1(t_ax_interp,1:length(t_ax_interp),sac_start_times));
     
    macro_sacs = find(sac_amps > 1);
    micro_sacs = find(sac_amps < 1);
    sac_t_axis = Trial_starts(1):sac_dt:Trial_ends(end);
    sac_times_binned = histc(sac_start_times,sac_t_axis);
    micro_times_binned = histc(sac_start_times(micro_sacs),sac_t_axis);
    macro_times_binned = histc(sac_start_times(macro_sacs),sac_t_axis);

             all_acorr(ee,:) = xcov(sac_times_binned,maxlag,'coeff');
    macro_acorr(ee,:) = xcov(macro_times_binned,maxlag,'coeff');
     micro_acorr(ee,:) = xcov(micro_times_binned,maxlag,'coeff');
 
%     usable_sacs = find(~isnan(sac_start_phase_bins));
%     macro_sacs = macro_sacs(ismember(macro_sacs,usable_sacs));
%     micro_sacs = micro_sacs(ismember(micro_sacs,usable_sacs));
%     sac_phases_binned = zeros(size(t_ax_interp));
%     sac_phases_binned(sac_start_phase_bins(usable_sacs)) = 1;
%     
%     micro_phases_binned = zeros(size(t_ax_interp));
%     micro_phases_binned(sac_start_phase_bins(micro_sacs)) = 1;
%     macro_phases_binned = zeros(size(t_ax_interp));
%     macro_phases_binned(sac_start_phase_bins(macro_sacs)) = 1;

     usable_sacs = find(~isnan(sac_start_phases));
    macro_sacs = macro_sacs(ismember(macro_sacs,usable_sacs));
    micro_sacs = micro_sacs(ismember(micro_sacs,usable_sacs));
    sac_phases_binned = histc(sac_start_phases,uniform_phase_grid);
    micro_phases_binned = histc(sac_start_phases(micro_sacs),uniform_phase_grid);
    macro_phases_binned = histc(sac_start_phases(macro_sacs),uniform_phase_grid);
    
   
    all_acorr_p(ee,:) = xcov(sac_phases_binned,maxlag_p,'coeff');
    macro_acorr_p(ee,:) = xcov(macro_phases_binned,maxlag_p,'coeff');
     micro_acorr_p(ee,:) = xcov(micro_phases_binned,maxlag_p,'coeff');
   
        
end


%%

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
% shadedErrorBar(lags,mean(micro_acorr),std(micro_acorr)/sqrt(7),{'color','r'})

figure
lags_p = -maxlag_p:maxlag_p;
lags_p = lags_p*sac_dp/(2*pi);
zl = find(lags_p == 0);
all_acorr_p(:,zl) = 0;
macro_acorr_p(:,zl) = 0;
micro_acorr_p(:,zl) = 0;
shadedErrorBar(lags_p,mean(all_acorr_p),std(all_acorr_p)/sqrt(7))
hold on
shadedErrorBar(lags_p,mean(macro_acorr_p),std(macro_acorr_p)/sqrt(7),{'color','b'})
% shadedErrorBar(lags_p,mean(micro_acorr_p),std(micro_acorr_p)/sqrt(7),{'color','r'})

%%
sm_fac = 5;
for ee = 1:7
    all_sm_acorr(ee,:) = smooth(all_acorr(ee,:),sm_fac,'lowess');
    macro_sm_acorr(ee,:) = smooth(macro_acorr(ee,:),sm_fac,'lowess');
    micro_sm_acorr(ee,:) = smooth(micro_acorr(ee,:),sm_fac,'lowess');
    all_sm_acorr_p(ee,:) = smooth(all_acorr_p(ee,:),sm_fac,'lowess');
    macro_sm_acorr_p(ee,:) = smooth(macro_acorr_p(ee,:),sm_fac,'lowess');
    micro_sm_acorr_p(ee,:) = smooth(micro_acorr_p(ee,:),sm_fac,'lowess');
end
figure
% shadedErrorBar(lags,mean(all_sm_acorr),std(all_sm_acorr)/sqrt(7))
hold on
shadedErrorBar(lags,mean(macro_sm_acorr),std(macro_sm_acorr)/sqrt(7),{'color','b'})
shadedErrorBar(lags,mean(micro_sm_acorr),std(micro_sm_acorr)/sqrt(7),{'color','r'})

figure
% shadedErrorBar(lags_p,mean(all_sm_acorr_p),std(all_sm_acorr_p)/sqrt(7))
hold on
shadedErrorBar(lags_p,mean(macro_sm_acorr_p),std(macro_sm_acorr_p)/sqrt(7),{'color','b'})
% shadedErrorBar(lags_p,mean(micro_sm_acorr_p),std(micro_sm_acorr_p)/sqrt(7),{'color','r'})

%%
% save expt2_eyeanalysis_g034 all_*