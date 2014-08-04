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

Pix2Deg = 0.018837;
Nyp = 1024;
Nxp = 1280;
use_win = [-7 7;-7 7];

min_fix_dur = 0.15;

Fs = 3e4;
dsf = 30;Fsd = Fs/dsf;
[b,a] = butter(2,[0.2 400]/(Fsd/2));
[b2,a2] = butter(2,[4 15]/(Fsd/2));
[b3,a3] = butter(2,[25 60]/(Fsd/2));
scales = logspace(log10(2.5),log10(40),25);
% scales = [5:10 12 14 16 20 25 30 35 40 45 50 55 60];
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

%%
% Expt_nu = [3 8 15 19 26 29];
Expt_nu = [13 14 15 16 25 28 29];
n_allunits = 96;
all_eyespeed = [];
all_sac_start_times = [];
all_sac_end_times = [];
all_sac_amps = [];
all_sac_expt_vec = [];
all_intrial = [];
all_t = [];
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
    
    [blink_data,in_blink,tot_disp_f] = get_blinks(eye_vals(:,3:4),eye_vals(:,1:2),eye_ts);
    in_blink(out_window) = 1;
    
    [sac_data,in_sac,eye_speed] = get_saccades(corrected_right,corrected_left,eye_ts,in_blink);
    
    blink_start_times = [blink_data(:).start_times];
    blink_stop_times = [blink_data(:).stop_times];
    sac_amps = [sac_data(:).amplitude];
    sac_start_times = [sac_data(:).start_time];
    sac_end_times = [sac_data(:).stop_time];
    
    
    fix_start_times = sac_end_times(1:end-1);
    fix_stop_times = sac_start_times(2:end);
    fix_start_inds = round(interp1(eye_ts,1:length(eye_ts),fix_start_times));
    fix_stop_inds = round(interp1(eye_ts,1:length(eye_ts),fix_stop_times));
    fix_sac_amps = sac_amps(1:end-1);
    blink_fix = zeros(length(fix_start_times),1);
    for i = 1:length(fix_start_times)
       if any(in_blink(fix_start_inds(i):fix_stop_inds(i)))
           blink_fix(i) = 1;
       end
    end
    fix_durs = fix_stop_times-fix_start_times;
    use_fixs = find(fix_durs >= min_fix_dur);
    fix_start_times = fix_start_times(use_fixs);
    fix_stop_times = fix_stop_times(use_fixs);
    fix_sac_amps = fix_sac_amps(use_fixs);
    blink_fix = blink_fix(use_fixs);
    
    %%
    clear ampgram all_V*
    use_lfps = [1:2:96];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d\n',ll);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = decimate(V,dsf);
        V = V(:);
        V_alpha = filtfilt(b2,a2,V);
        V_gamma = filtfilt(b3,a3,V);
        V_bb = filtfilt(b,a,V);
        all_Vbb(:,ll) = V_bb;
        all_Valpha(:,ll) = V_alpha;
        all_Vgamma(:,ll) = V_gamma;
        ampgram(:,:,ll) = abs(cwt(V,scales,'cmor1-1'))';
    end
    
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
    
    ampgram = zscore(ampgram);
    all_Vbb = zscore(all_Vbb);
    all_Valpha = zscore(all_Valpha);
    all_Vgamma = zscore(all_Vgamma);
    
%     sm_win = 5;
%     cur_binned_spks = nan(96,length(t_ax));
%     for j = 1:96
%         cur_binned_spks(j,:) = smooth(histc(Clusters{j}.times,t_ax),sm_win);
%     end

    %%
    forwardlag = round(Fsd*0.75);
    backlag = round(Fsd*0.1);
    lags = -backlag:forwardlag;
    
    fix_start_inds = round(interp1(t_ax,1:length(t_ax),fix_start_times));
    fix_stop_inds = round(interp1(t_ax,1:length(t_ax),fix_stop_times));
%     use_fixs = find(fix_durs >= min_fix_dur & blink_fix' == 0);
    use_fixs = find(fix_start_inds > forwardlag & fix_stop_inds < (length(t_ax)-forwardlag) & blink_fix' == 0);
    fix_start_inds = fix_start_inds(use_fixs);
    fix_stop_inds = fix_stop_inds(use_fixs);
    fix_sac_amps = fix_sac_amps(use_fixs);
    n_fixs = length(fix_stop_inds);
    
%     micro_sac_inds = find(fix_sac_amps < 1);
%     macro_sac_inds = find(fix_sac_amps > 1);
    
    %%
    trg_ampgrams = zeros(length(lags),length(scales),length(use_lfps));
    trg_Valpha = zeros(length(lags),length(use_lfps));
    trg_Vbb = zeros(length(lags),length(use_lfps));
    n_cnts = zeros(length(lags),1);
    
%     trg_ampgrams_mic = zeros(length(use_lfps),length(scales),length(lags));
%     trg_V_mic = zeros(length(use_lfps),length(lags));
%     n_cnts_mic = zeros(1,length(lags));
%     
%     trg_ampgrams_mac = zeros(length(use_lfps),length(scales),length(lags));
%     trg_V_mac = zeros(length(use_lfps),length(lags));
%     n_cnts_mac = zeros(1,length(lags));
    for i = 1:n_fixs
        cur_inds = (fix_start_inds(i)-backlag):(fix_start_inds(i)+forwardlag);
        cur_inds(cur_inds > fix_stop_inds(i)) = [];
        cl = length(cur_inds);
        trg_ampgrams(1:cl,:,:) = trg_ampgrams(1:cl,:,:) + ampgram(cur_inds,:,:);
        trg_Vbb(1:cl,:) = trg_Vbb(1:cl,:) + all_Vbb(cur_inds,:);
        trg_Valpha(1:cl,:) = trg_Valpha(1:cl,:) + all_Valpha(cur_inds,:);
        n_cnts(1:cl) = n_cnts(1:cl) + 1;
        
%         if ismember(i,micro_sac_inds)
%             trg_ampgrams_mic(:,:,1:cl) = trg_ampgrams_mic(:,:,1:cl) + ampgram(:,:,cur_inds);
%             trg_V_mic(:,1:cl) = trg_V_mic(:,1:cl) + all_V(:,cur_inds);
%             n_cnts_mic(1:cl) = n_cnts_mic(1:cl) + 1;
%         end
%         if ismember(i,macro_sac_inds)
%             trg_ampgrams_mac(:,:,1:cl) = trg_ampgrams_mac(:,:,1:cl) + ampgram(:,:,cur_inds);
%             trg_V_mac(:,1:cl) = trg_V_mac(:,1:cl) + all_V(:,cur_inds);
%             n_cnts_mac(1:cl) = n_cnts_mac(1:cl) + 1;
%         end
    end
    trg_Vbb = bsxfun(@rdivide,trg_Vbb,n_cnts);
    trg_Valpha = bsxfun(@rdivide,trg_Valpha,n_cnts);
    trg_ampgrams = bsxfun(@rdivide,trg_ampgrams,n_cnts);
%     trg_V_mic = bsxfun(@rdivide,trg_V_mic,n_cnts_mic);
%     trg_ampgrams_mic = bsxfun(@rdivide,trg_ampgrams_mic,reshape(n_cnts_mic,1,1,length(lags)));
%     trg_V_mac = bsxfun(@rdivide,trg_V_mac,n_cnts_mac);
%     trg_ampgrams_mac = bsxfun(@rdivide,trg_ampgrams_mac,reshape(n_cnts_mac,1,1,length(lags)));
    
    all_trg_Vbb(ee,:,:) = trg_Vbb;
    all_trg_Valpha(ee,:,:) = trg_Valpha;
    all_trg_ampgrams(ee,:,:,:) = trg_ampgrams;
%     all_trg_V_mic(ee,:,:) = trg_V_mic;
%     all_trg_ampgrams_mic(ee,:,:,:) = trg_ampgrams_mic;
%     all_trg_V_mac(ee,:,:) = trg_V_mac;
%     all_trg_ampgrams_mac(ee,:,:,:) = trg_ampgrams_mac;
    %%
    
    all_sac_start_times = [all_sac_start_times; sac_start_times(:)];
    all_sac_end_times = [all_sac_end_times; sac_end_times(:)];
    all_sac_amps = [all_sac_amps; sac_amps(:)];
    all_sac_expt_vec = [all_sac_expt_vec; ee*ones(length(sac_amps),1)];
    all_eyespeed = [all_eyespeed; eye_speed(:)];
    all_t = [all_t; eye_ts(:)];
    
end

%%
avg_trg_Vbb = squeeze(nanmean(all_trg_Vbb));
sem_trg_Vbb = squeeze(nanstd(all_trg_Vbb))/sqrt(7);
avg_trg_Valpha = squeeze(nanmean(all_trg_Valpha));
sem_trg_Valpha = squeeze(nanstd(all_trg_Valpha))/sqrt(7);
avg_trg_Vampgrams = squeeze(nanmean(all_trg_ampgrams));

