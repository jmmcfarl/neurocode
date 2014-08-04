%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12
cd G029/
load ./eye_calibration_data
cd ~/Data/bruce/7_15_12/G034/
load ./G034Expts.mat
load ./jbeG034.em.mat
em_data = Expt; clear Expt

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt3_fixbased_data 
trial_expt_num = full_expt_vec(trial_start_inds);
trial_stop_times = full_t(trial_stop_inds);
trial_start_times = full_t(trial_start_inds);
trial_durs = trial_stop_times - trial_start_times;
min_trial_dur = 0.5;
u_trials = find(trial_durs > min_trial_dur);
trial_start_times = trial_start_times(u_trials);
trial_stop_times = trial_stop_times(u_trials);
trial_imnum = trial_imnum(u_trials);
trial_expt_num = trial_expt_num(u_trials);

Pix2Deg = 0.018837;
Nyp = 1024;
Nxp = 1280;
use_win = [-7 7;-7 7];

min_fix_dur = 0.4;
use_lfps = [1:2:96];

Fs = 3e4;
dsf = 100; Fsd = Fs/dsf;
[b2,a2] = butter(2,[30 60]/(Fsd/2));
[b,a] = butter(2,[0.5]/(Fsd/2),'high');

%%
% Expt_nu = [3 8 15 19 26 29];
% Expt_nu = [13 14 15 16 25 28 29];
Expt_nu = [7:12]; %expt 3 34
n_allunits = 96;
pooled_Vbb_early = [];
pooled_Vbb_late = [];
pooled_Vgamma = [];
for ee = 1:length(Expt_nu)
    
    fprintf('Analyzing Expt %d of %d\n',ee,length(Expt_nu));
    load(sprintf('Expt%dClusterTimes.mat',Expt_nu(ee)));
    
    Trial_starts = [Expts{Expt_nu(ee)}.Trials(:).Start]/1e4;
    Trial_ends = [Expts{Expt_nu(ee)}.Trials(:).End]/1e4;
    endEvents = [Expts{Expt_nu(ee)}.Trials(:).endevent];
    Trial_durs = (Trial_ends-Trial_starts);
    
    [eye_vals,eye_ts,eye_insample] = get_eye_tracking_data(em_data,[Trial_starts(1) Trial_ends(end)]);
    
    clear lcorrected_* corrected*
    %correct eye positions
    corrected_left = bsxfun(@plus,eye_vals(:,1:2)*left_gain,left_offset);
    corrected_right = bsxfun(@plus,eye_vals(:,3:4)*right_gain,right_offset);
    
    avg_eyepos = 0.5*corrected_left + 0.5*corrected_right;
    out_window = find(avg_eyepos(:,1) < use_win(1,1) | avg_eyepos(:,1) > use_win(1,2) | ...
        avg_eyepos(:,2) < use_win(2,1) | avg_eyepos(:,2) > use_win(2,2));        
    
    %%
    clear ampgram all_V*
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d\n',ll);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = decimate(V,dsf);
        V_bb = filtfilt(b,a,V);
        V_gamma = filtfilt(b2,a2,V);
        all_Vbb(:,ll) = V_bb;
        all_Vgamma(:,ll) = V_gamma;
    end
    
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
    
    all_Vbb = zscore(all_Vbb);
    all_gphase = angle(hilbert(all_Vgamma));
    
    
    %%
    cur_trials = find(trial_expt_num == Expt_nu(ee));
    cur_trial_start_inds = round(interp1(t_ax,1:length(t_ax),trial_start_times(cur_trials)));
    cur_trial_stop_inds = round(interp1(t_ax,1:length(t_ax),trial_stop_times(cur_trials)));
    
    n_trials = length(cur_trial_start_inds);
    
    early_inds = [cur_trial_start_inds(:) cur_trial_start_inds(:)+round(0.15*Fsd)];
    late_inds = [cur_trial_start_inds(:)+round(0.15*Fsd) cur_trial_start_inds(:)+round(0.5*Fsd)];
    all_inds = [cur_trial_start_inds(:) cur_trial_start_inds(:)+round(0.5*Fsd)];
    
    %%
    trg_Vbb_early = zeros(round(0.15*Fsd)+1,n_trials,length(use_lfps));
    trg_Vbb_late = zeros(round(0.35*Fsd)+1,n_trials,length(use_lfps));
    trg_Vgamma = zeros(round(0.5*Fsd)+1,n_trials,length(use_lfps));
    for i = 1:n_trials
        cur_inds = early_inds(i,1):early_inds(i,2);
        trg_Vbb_early(:,i,:) = all_Vbb(cur_inds,:);
        cur_inds = late_inds(i,1):late_inds(i,2);
        trg_Vbb_late(:,i,:) = all_Vbb(cur_inds,:);
        cur_inds = all_inds(i,1):all_inds(i,2);
        trg_Vgamma(:,i,:) = all_gphase(cur_inds,:);
    end
    %%
    pooled_Vbb_early = cat(2,pooled_Vbb_early,trg_Vbb_early);
    pooled_Vbb_late = cat(2,pooled_Vbb_late,trg_Vbb_late);
    pooled_Vgamma = cat(2,pooled_Vgamma,trg_Vgamma);
end

%%
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;

n_bins = round(Fsd*0.5)+1;
n_pairs = length(use_lfps)*(length(use_lfps)-1)/2;
phase_cor = nan(length(use_lfps)*(length(use_lfps)-1),n_bins);
cur_dist = nan(length(use_lfps)*(length(use_lfps)-1),1);
cnt = 1;
for ch1 = 1:length(use_lfps)
    ch1
    for ch2 = (ch1 + 1):length(use_lfps)
        cur_dist(cnt) = sqrt((X_pos(use_lfps(ch1))-X_pos(use_lfps(ch2)))^2 + (Y_pos(use_lfps(ch1))-Y_pos(use_lfps(ch2)))^2);
        cur_phase_sig1 = squeeze(pooled_Vgamma(:,:,ch1))';
        cur_phase_sig2 = squeeze(pooled_Vgamma(:,:,ch2))';
        for ll = 1:n_bins
            [phase_cor(cnt,ll)] = circ_corrcc(cur_phase_sig1(:,ll),cur_phase_sig2(:,ll));
%             [phase_cor(ch1,ch2,ll)] = circ_kappa(cur_phase_sig1(:,ll)-cur_phase_sig2(:,ll));
        end
        cnt = cnt + 1;
    end
end

cnt = 1;
phase_cor_sm = nan(size(phase_cor));
for ch1 = 1:length(use_lfps)
    for ch2 = (ch1 + 1):length(use_lfps)
        phase_cor_sm(cnt,:) = smooth(phase_cor(cnt,:),5);
        cnt = cnt + 1;
    end
end
[~,dist_ord] = sort(cur_dist);

%%
params.Fs = Fsd;
params.tapers = [1 1];
params.trialave = 1;
for ch1 = 1:length(use_lfps);
    ch1
    for ch2 = ch1 + 1:length(use_lfps);
        [C_early(ch1,ch2,:),phi,S12,S1,S2,f_early]=coherencyc(squeeze(pooled_Vbb_early(:,:,ch1)),squeeze(pooled_Vbb_early(:,:,ch2)),params);
        [C_late(ch1,ch2,:),phi,S12,S1,S2,f_late]=coherencyc(squeeze(pooled_Vbb_late(:,:,ch1)),squeeze(pooled_Vbb_late(:,:,ch2)),params);
    end
end

%%
pooled_adj_early = bsxfun(@minus,pooled_Vbb_early,mean(pooled_Vbb_early));
pooled_adj_late = bsxfun(@minus,pooled_Vbb_late,mean(pooled_Vbb_late));
%%
[N,nt,nch] = size(pooled_adj_late);
% [N,nt,nch] = size(pooled_adj_early);
nfft = 2^nextpow2(N);
[f,findx]=getfgrid(Fsd,nfft,[0 Fsd/2]);
for ch1 = 1:length(use_lfps);
    ch1
    for ch2 = ch1 + 1:length(use_lfps);
        data1 = squeeze(pooled_adj_late(:,:,ch1));
        data2 = squeeze(pooled_adj_late(:,:,ch2));
%         data1 = squeeze(pooled_adj_early(:,:,ch1));
%         data2 = squeeze(pooled_adj_early(:,:,ch2));
        j1 = fft(data1,nfft);
        j2 = fft(data2,nfft);
        s12 = conj(j1).*j2;
        s1 = conj(j1).*j1;
        s2 = conj(j2).*j2;
        s12 = mean(s12,2);
        s1 = mean(s1,2);
        s2 = mean(s2,2);
        C12=s2./sqrt(s1.*s2);
        late_C(ch1,ch2,:) = C12(findx);
%         early_C(ch1,ch2,:) = C12(findx);
    end
end

%%
 figure
 plot(f,C12(findx))