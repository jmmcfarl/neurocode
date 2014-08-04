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
use_lfps = [1:5:96];

Fs = 3e4;
dsf = 60;Fsd = Fs/dsf;
[b,a] = butter(2,[1 100]/(Fsd/2));
[b2,a2] = butter(2,[4 12]/(Fsd/2));
[b3,a3] = butter(2,[30 70]/(Fsd/2));

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
    sac_start_times = sac_start_times(use_fixs);
    blink_fix = blink_fix(use_fixs);
    
    %%
    clear ampgram all_V*
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
        all_Vgamma(:,ll) = abs(hilbert(V_gamma));
    end
    
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
    
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
    backlag = round(Fsd*0.25);
    lags = -backlag:forwardlag;
    
    sac_start_inds = round(interp1(t_ax,1:length(t_ax),sac_start_times));
    fix_start_inds = round(interp1(t_ax,1:length(t_ax),fix_start_times));
    fix_stop_inds = round(interp1(t_ax,1:length(t_ax),fix_stop_times));
    use_fixs = find(fix_start_inds > forwardlag & fix_stop_inds < (length(t_ax)-forwardlag) & blink_fix' == 0 & ~isnan(sac_start_inds));
    fix_start_inds = fix_start_inds(use_fixs);
    fix_stop_inds = fix_stop_inds(use_fixs);
    fix_sac_amps = fix_sac_amps(use_fixs);
    sac_start_inds = sac_start_inds(use_fixs);
    fix_start_times = fix_start_times(use_fixs);
    fix_stop_times = fix_stop_times(use_fixs);
    n_fixs = length(fix_stop_inds);
    
    big_sacs = find(fix_sac_amps > 1);
    micro_sacs = find(fix_sac_amps < 1);
    
    %%
    trg_Valpha = zeros(length(lags),length(use_lfps));
    trg_Vgamma = zeros(length(lags),length(use_lfps));
    trg_Vbb = zeros(length(lags),length(use_lfps));
    n_cnts = zeros(length(lags),1);
    strg_Valpha = zeros(length(lags),length(use_lfps));
    strg_Vgamma = zeros(length(lags),length(use_lfps));
    strg_Vbb = zeros(length(lags),length(use_lfps));
    n_scnts = zeros(length(lags),1);
    mactrg_Valpha = zeros(length(lags),length(use_lfps));
    mactrg_Vgamma = zeros(length(lags),length(use_lfps));
    mactrg_Vbb = zeros(length(lags),length(use_lfps));
    n_maccnts = zeros(length(lags),1);
    mictrg_Valpha = zeros(length(lags),length(use_lfps));
    mictrg_Vgamma = zeros(length(lags),length(use_lfps));
    mictrg_Vbb = zeros(length(lags),length(use_lfps));
    n_miccnts = zeros(length(lags),1);
    for i = 1:n_fixs
        cur_inds = (fix_start_inds(i)-backlag):(fix_start_inds(i)+forwardlag);
        cur_inds(cur_inds > fix_stop_inds(i)) = [];
        cl = length(cur_inds);
        trg_Vbb(1:cl,:) = trg_Vbb(1:cl,:) + all_Vbb(cur_inds,:);
        trg_Valpha(1:cl,:) = trg_Valpha(1:cl,:) + all_Valpha(cur_inds,:);
        trg_Vgamma(1:cl,:) = trg_Vgamma(1:cl,:) + all_Vgamma(cur_inds,:);
        n_cnts(1:cl) = n_cnts(1:cl) + 1;
        
        if ismember(i,big_sacs)
            mactrg_Vbb(1:cl,:) = mactrg_Vbb(1:cl,:) + all_Vbb(cur_inds,:);
            mactrg_Valpha(1:cl,:) = mactrg_Valpha(1:cl,:) + all_Valpha(cur_inds,:);
            mactrg_Vgamma(1:cl,:) = mactrg_Vgamma(1:cl,:) + all_Vgamma(cur_inds,:);
            n_maccnts(1:cl) = n_maccnts(1:cl) + 1;          
        end
        if ismember(i,micro_sacs)
            mictrg_Vbb(1:cl,:) = mictrg_Vbb(1:cl,:) + all_Vbb(cur_inds,:);
            mictrg_Valpha(1:cl,:) = mictrg_Valpha(1:cl,:) + all_Valpha(cur_inds,:);
            mictrg_Vgamma(1:cl,:) = mictrg_Vgamma(1:cl,:) + all_Vgamma(cur_inds,:);
            n_miccnts(1:cl) = n_miccnts(1:cl) + 1;          
        end
        
        cur_inds = (sac_start_inds(i)-backlag):(sac_start_inds(i)+forwardlag);
        cur_inds(cur_inds > fix_stop_inds(i)) = [];
        cl = length(cur_inds);
        strg_Vbb(1:cl,:) = strg_Vbb(1:cl,:) + all_Vbb(cur_inds,:);
        strg_Valpha(1:cl,:) = strg_Valpha(1:cl,:) + all_Valpha(cur_inds,:);
        strg_Vgamma(1:cl,:) = strg_Vgamma(1:cl,:) + all_Vgamma(cur_inds,:);
        n_scnts(1:cl) = n_scnts(1:cl) + 1;
        
    end
    trg_Vbb = bsxfun(@rdivide,trg_Vbb,n_cnts);
    trg_Valpha = bsxfun(@rdivide,trg_Valpha,n_cnts);
    trg_Vgamma = bsxfun(@rdivide,trg_Vgamma,n_cnts);
    strg_Vbb = bsxfun(@rdivide,strg_Vbb,n_scnts);
    strg_Valpha = bsxfun(@rdivide,strg_Valpha,n_scnts);
    strg_Vgamma = bsxfun(@rdivide,strg_Vgamma,n_scnts);
    mactrg_Vbb = bsxfun(@rdivide,mactrg_Vbb,n_maccnts);
    mactrg_Valpha = bsxfun(@rdivide,mactrg_Valpha,n_maccnts);
    mactrg_Vgamma = bsxfun(@rdivide,mactrg_Vgamma,n_maccnts);
    mictrg_Vbb = bsxfun(@rdivide,mictrg_Vbb,n_miccnts);
    mictrg_Valpha = bsxfun(@rdivide,mictrg_Valpha,n_miccnts);
    mictrg_Vgamma = bsxfun(@rdivide,mictrg_Vgamma,n_miccnts);
    
    all_trg_Vbb(ee,:,:) = trg_Vbb;
    all_trg_Valpha(ee,:,:) = trg_Valpha;
    all_trg_Vgamma(ee,:,:) = trg_Vgamma;
    all_strg_Vbb(ee,:,:) = strg_Vbb;
    all_strg_Valpha(ee,:,:) = strg_Valpha;
    all_strg_Vgamma(ee,:,:) = strg_Vgamma;
    all_mactrg_Vbb(ee,:,:) = mactrg_Vbb;
    all_mactrg_Valpha(ee,:,:) = mactrg_Valpha;
    all_mactrg_Vgamma(ee,:,:) = mactrg_Vgamma;
    all_mictrg_Vbb(ee,:,:) = mictrg_Vbb;
    all_mictrg_Valpha(ee,:,:) = mictrg_Valpha;
    all_mictrg_Vgamma(ee,:,:) = mictrg_Vgamma;
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
avg_trg_Vgamma = squeeze(nanmean(all_trg_Vgamma));
sem_trg_Vgamma = squeeze(nanstd(all_trg_Vgamma))/sqrt(7);

avg_strg_Vbb = squeeze(nanmean(all_strg_Vbb));
sem_strg_Vbb = squeeze(nanstd(all_strg_Vbb))/sqrt(7);
avg_strg_Valpha = squeeze(nanmean(all_strg_Valpha));
sem_strg_Valpha = squeeze(nanstd(all_strg_Valpha))/sqrt(7);
avg_strg_Vgamma = squeeze(nanmean(all_strg_Vgamma));
sem_strg_Vgamma = squeeze(nanstd(all_strg_Vgamma))/sqrt(7);

avg_mactrg_Vbb = squeeze(nanmean(all_mactrg_Vbb));
sem_mactrg_Vbb = squeeze(nanstd(all_mactrg_Vbb))/sqrt(7);
avg_mactrg_Valpha = squeeze(nanmean(all_mactrg_Valpha));
sem_mactrg_Valpha = squeeze(nanstd(all_mactrg_Valpha))/sqrt(7);
avg_mactrg_Vgamma = squeeze(nanmean(all_mactrg_Vgamma));
sem_mactrg_Vgamma = squeeze(nanstd(all_mactrg_Vgamma))/sqrt(7);

avg_mictrg_Vbb = squeeze(nanmean(all_mictrg_Vbb));
sem_mictrg_Vbb = squeeze(nanstd(all_mictrg_Vbb))/sqrt(7);
avg_mictrg_Valpha = squeeze(nanmean(all_mictrg_Valpha));
sem_mictrg_Valpha = squeeze(nanstd(all_mictrg_Valpha))/sqrt(7);
avg_mictrg_Vgamma = squeeze(nanmean(all_mictrg_Vgamma));
sem_mictrg_Vgamma = squeeze(nanstd(all_mictrg_Vgamma))/sqrt(7);

%%
save sac_trg_avgs_v2 lags Fsd avg_* sem_* 

%%
for cc = 1:length(use_lfps)
    shadedErrorBar(lags/Fsd,avg_trg_Vbb(:,cc),sem_trg_Vbb(:,cc));
    hold on
    shadedErrorBar(lags/Fsd,avg_trg_Vgamma(:,cc),sem_trg_Vgamma(:,cc),{'color','r'});
    shadedErrorBar(lags/Fsd,avg_trg_Valpha(:,cc),sem_trg_Valpha(:,cc),{'color','b'});
    xlim([-0.1 0.55])
    pause
    clf
end

%%
for cc = 1:length(use_lfps)
    shadedErrorBar(lags/Fsd,avg_strg_Vbb(:,cc),sem_strg_Vbb(:,cc));
    hold on
    shadedErrorBar(lags/Fsd,avg_trg_Vbb(:,cc),sem_strg_Vbb(:,cc),{'color','b'});
    xlim([-0.1 0.55])
    pause
    clf
end

%%
for cc = 1:length(use_lfps)
    plot(lags/Fsd,avg_trg_Vbb(:,cc),'k','linewidth',2);
    hold on
    shadedErrorBar(lags/Fsd,avg_mactrg_Vbb(:,cc),sem_mactrg_Vbb(:,cc),{'color','r'});
    shadedErrorBar(lags/Fsd,avg_mictrg_Vbb(:,cc),sem_mictrg_Vbb(:,cc),{'color','b'});
    xlim([-0.1 0.55])
    pause
    clf
end

%%
for cc = 1:length(use_lfps)
    plot(lags/Fsd,avg_trg_Vgamma(:,cc),'k','linewidth',2);
    hold on
    shadedErrorBar(lags/Fsd,avg_mactrg_Vgamma(:,cc),sem_mactrg_Vgamma(:,cc),{'color','r'});
    shadedErrorBar(lags/Fsd,avg_mictrg_Vgamma(:,cc),sem_mictrg_Vgamma(:,cc),{'color','b'});
    xlim([-0.1 0.55])
    pause
    clf
end