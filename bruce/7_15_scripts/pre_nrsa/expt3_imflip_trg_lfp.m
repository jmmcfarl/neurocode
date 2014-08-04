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

load ./Expt3_newcompiled_data_d1p25_htres full_t full_stim_ids full_image_vec full_expt_vec
flip_ids = 1 + find(diff(full_image_vec)~=0);
shift_ids = 1 + find(diff(full_stim_ids)~=0);
shift_times = full_t(shift_ids);
is_flip = ismember(shift_ids,flip_ids);
% shift_times = shift_times(is_flip);
% shift_ids = shift_ids(is_flip);

Pix2Deg = 0.018837;
Nyp = 1024;
Nxp = 1280;
use_win = [-7 7;-7 7];

min_fix_dur = 0.15;
use_lfps = [1];

Fs = 3e4;
dsf = 30;Fsd = Fs/dsf;
[b,a] = butter(2,[1 100]/(Fsd/2));
[b2,a2] = butter(2,[4 12]/(Fsd/2));
[b3,a3] = butter(2,[30 70]/(Fsd/2));

forwardlag = round(Fsd*0.7);
backlag = round(Fsd*0.4);
lags = -backlag:forwardlag;

%%
% Expt_nu = [3 8 15 19 26 29];
% Expt_nu = [13 14 15 16 25 28 29];
Expt_nu = [7:12]; %expt 3 34
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
    
    uset = find(full_expt_vec==Expt_nu(ee));
    cur_im_inds = ~isnan(full_image_vec(uset));
    
    sac_start_inds = round(interp1(t_ax,1:length(t_ax),sac_start_times));
    fix_start_inds = round(interp1(t_ax,1:length(t_ax),fix_start_times));
    fix_stop_inds = round(interp1(t_ax,1:length(t_ax),fix_stop_times));
    use_fixs = find(fix_start_inds > forwardlag & fix_stop_inds < (length(t_ax)-forwardlag) & blink_fix' == 0 & ~isnan(sac_start_inds));
    fix_start_inds = fix_start_inds(use_fixs);
    fix_stop_inds = fix_stop_inds(use_fixs);
    fix_sac_amps = fix_sac_amps(use_fixs);
    sac_start_inds = sac_start_inds(use_fixs);
    n_fixs = length(fix_stop_inds);
    big_sacs = find(fix_sac_amps > 1);
%     old_fix_start_ids
%     image_sacs = find(cur_im_inds(fix_start_inds));
    
    %%
    trg_Valpha = zeros(length(lags),length(use_lfps));
    trg_Vgamma = zeros(length(lags),length(use_lfps));
    trg_Vbb = zeros(length(lags),length(use_lfps));
    n_cnts = zeros(length(lags),1);
    strg_Valpha = zeros(length(lags),length(use_lfps));
    strg_Vgamma = zeros(length(lags),length(use_lfps));
    strg_Vbb = zeros(length(lags),length(use_lfps));
    n_scnts = zeros(length(lags),1);
    for i = 1:n_fixs
        cur_inds = (fix_start_inds(i)-backlag):(fix_start_inds(i)+forwardlag);
        cur_inds(cur_inds > fix_stop_inds(i)) = [];
        cl = length(cur_inds);
        trg_Vbb(1:cl,:) = trg_Vbb(1:cl,:) + all_Vbb(cur_inds,:);
        trg_Valpha(1:cl,:) = trg_Valpha(1:cl,:) + all_Valpha(cur_inds,:);
        trg_Vgamma(1:cl,:) = trg_Vgamma(1:cl,:) + all_Vgamma(cur_inds,:);
        n_cnts(1:cl) = n_cnts(1:cl) + 1;
        
        if ismember(i,big_sacs)
        cur_inds = (fix_start_inds(i)-backlag):(fix_start_inds(i)+forwardlag);
        cur_inds(cur_inds > fix_stop_inds(i)) = [];
        cl = length(cur_inds);
        strg_Vbb(1:cl,:) = strg_Vbb(1:cl,:) + all_Vbb(cur_inds,:);
        strg_Valpha(1:cl,:) = strg_Valpha(1:cl,:) + all_Valpha(cur_inds,:);
        strg_Vgamma(1:cl,:) = strg_Vgamma(1:cl,:) + all_Vgamma(cur_inds,:);
        n_scnts(1:cl) = n_scnts(1:cl) + 1;
        end
    end
    trg_Vbb = bsxfun(@rdivide,trg_Vbb,n_cnts);
    trg_Valpha = bsxfun(@rdivide,trg_Valpha,n_cnts);
    trg_Vgamma = bsxfun(@rdivide,trg_Vgamma,n_cnts);
    strg_Vbb = bsxfun(@rdivide,strg_Vbb,n_scnts);
    strg_Valpha = bsxfun(@rdivide,strg_Valpha,n_scnts);
    strg_Vgamma = bsxfun(@rdivide,strg_Vgamma,n_scnts);
    
    all_trg_Vbb(ee,:,:) = trg_Vbb;
    all_trg_Valpha(ee,:,:) = trg_Valpha;
    all_trg_Vgamma(ee,:,:) = trg_Vgamma;
    all_strg_Vbb(ee,:,:) = strg_Vbb;
    all_strg_Valpha(ee,:,:) = strg_Valpha;
    all_strg_Vgamma(ee,:,:) = strg_Vgamma;
 
    %%
    cur_shifts = find(full_expt_vec(shift_ids) == Expt_nu(ee));
    cur_shifts(shift_ids(cur_shifts) < backlag) = [];
    cur_shift_inds = round(interp1(t_ax,1:length(t_ax),full_t(shift_ids(cur_shifts))));
%     cur_shift_inds(cur_shift_inds < backlag) = [];
    
    ftrg_Valpha = zeros(length(lags),length(use_lfps));
    ftrg_Vgamma = zeros(length(lags),length(use_lfps));
    ftrg_Vbb = zeros(length(lags),length(use_lfps));
    n_fcnts = zeros(length(lags),1);
    for i = 1:length(cur_shift_inds)
        cur_inds = (cur_shift_inds(i)-backlag):(cur_shift_inds(i)+forwardlag);
        cl = length(cur_inds);
        ftrg_Vbb(1:cl,:) = ftrg_Vbb(1:cl,:) + all_Vbb(cur_inds,:);
        ftrg_Valpha(1:cl,:) = ftrg_Valpha(1:cl,:) + all_Valpha(cur_inds,:);
        ftrg_Vgamma(1:cl,:) = ftrg_Vgamma(1:cl,:) + all_Vgamma(cur_inds,:);
        n_fcnts(1:cl) = n_fcnts(1:cl) + 1;
    end
    ftrg_Vbb = bsxfun(@rdivide,ftrg_Vbb,n_fcnts);
    ftrg_Valpha = bsxfun(@rdivide,ftrg_Valpha,n_fcnts);
    ftrg_Vgamma = bsxfun(@rdivide,ftrg_Vgamma,n_fcnts);
    
    all_ftrg_Vbb(ee,:,:) = ftrg_Vbb;
    all_ftrg_Valpha(ee,:,:) = ftrg_Valpha;
    all_ftrg_Vgamma(ee,:,:) = ftrg_Vgamma;

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

avg_ftrg_Vbb = squeeze(nanmean(all_ftrg_Vbb));
sem_ftrg_Vbb = squeeze(nanstd(all_ftrg_Vbb))/sqrt(7);
avg_ftrg_Valpha = squeeze(nanmean(all_ftrg_Valpha));
sem_ftrg_Valpha = squeeze(nanstd(all_ftrg_Valpha))/sqrt(7);
avg_ftrg_Vgamma = squeeze(nanmean(all_ftrg_Vgamma));
sem_ftrg_Vgamma = squeeze(nanstd(all_ftrg_Vgamma))/sqrt(7);

%%
save expt3_flip_sac_trg_avgs avg_trg* avg_strg* avg_ftrg* lags Fsd sem_* 

%%
cd ~/Data/bruce/7_15_12/G034/
load ./expt3_flip_sac_trg_avgs
lags_3 = lags/Fsd;
avg_3_ftrg_Vbb = avg_ftrg_Vbb;
sem_3_ftrg_Vbb = sem_ftrg_Vbb;
load ./sac_trg_avgs_v3.mat
lags_2 = lags/Fsd;
avg_2_strg_Vbb = avg_strg_Vbb;
sem_2_strg_Vbb = sem_strg_Vbb;
% cd ~/Data/bruce/7_15_12/G035/
% load ./expt4_flip_sac_trg_avgs
% lags_4 = lags/Fsd;
close all
for cc = 1:length(use_lfps)
    shadedErrorBar(lags_3,avg_3_ftrg_Valpha(:,cc),sem_3_ftrg_Valpha(:,cc));
    hold on
%     shadedErrorBar(lags_2,avg_2_strg_Vbb(:,cc),sem_2_strg_Vbb(:,cc),{'color','r'});
%     shadedErrorBar(lags_4,avg_ftrg_Vbb(:,cc),sem_ftrg_Vbb(:,cc),{'color','b'});
    xlim([-0.1 0.5])
    yl = ylim();
    line([0 0],yl,'color','k')
    pause
    clf
end

%%
cd ~/Data/bruce/7_15_12/G034/
load ./expt3_flip_sac_trg_avgs
lags_3 = lags/Fsd;
avg_3_ftrg_Valpha = avg_ftrg_Valpha;
sem_3_ftrg_Valpha = sem_ftrg_Valpha;
load ./sac_trg_avgs_v3
lags_2 = lags/Fsd;
avg_2_strg_Valpha = avg_strg_Valpha;
sem_2_strg_Valpha = sem_strg_Valpha;
cd ~/Data/bruce/7_15_12/G035/
load ./expt4_flip_sac_trg_avgs
lags_4 = lags/Fsd;
close all
for cc = 1:length(use_lfps)
    shadedErrorBar(lags_3,avg_3_ftrg_Valpha(:,cc),sem_3_ftrg_Valpha(:,cc));
    hold on
    shadedErrorBar(lags_2,avg_2_strg_Valpha(:,cc),sem_2_strg_Valpha(:,cc),{'color','r'});
%     shadedErrorBar(lags_4,avg_ftrg_Valpha(:,cc),sem_ftrg_Valpha(:,cc),{'color','b'});
    xlim([-0.2 0.5])
    yl = ylim();
    line([0 0],yl,'color','k')
    pause
    clf
end

%%
cd ~/Data/bruce/7_15_12/G034/
load ./expt3_flip_sac_trg_avgs
lags_3 = lags/Fsd;
avg_3_ftrg_Vgamma = avg_ftrg_Vgamma;
sem_3_ftrg_Vgamma = sem_ftrg_Vgamma;
load ./sac_trg_avgs_v3
lags_2 = lags/Fsd;
avg_2_strg_Vgamma = avg_strg_Vgamma;
sem_2_strg_Vgamma = sem_strg_Vgamma;
cd ~/Data/bruce/7_15_12/G035/
load ./expt4_flip_sac_trg_avgs
lags_4 = lags/Fsd;
close all
for cc = 1:length(use_lfps)
    shadedErrorBar(lags_3,avg_3_ftrg_Vgamma(:,cc),sem_3_ftrg_Vgamma(:,cc));
    hold on
    shadedErrorBar(lags_2,avg_2_strg_Vgamma(:,cc),sem_2_strg_Vgamma(:,cc),{'color','r'});
    shadedErrorBar(lags_4,avg_ftrg_Vgamma(:,cc),sem_ftrg_Vgamma(:,cc),{'color','b'});
    xlim([-0.2 0.5])
    yl = ylim();
    line([0 0],yl,'color','k')
    pause
    clf
end

%%
load ./expt3_flip_sac_trg_avgs
lags_3 = lags/Fsd;
load ./sac_trg_avgs_v3
close all
for cc = 1:length(use_lfps)
    shadedErrorBar(lags_3,avg_ftrg_Valpha(:,cc),sem_ftrg_Valpha(:,cc));
    hold on
    shadedErrorBar(lags/Fsd,avg_strg_Valpha(:,cc),sem_strg_Valpha(:,cc),{'color','r'});
    xlim([-0.4 0.7])
    pause
    clf
end

%%
for cc = 1:length(use_lfps)
    shadedErrorBar(lags/Fsd,avg_trg_Vbb(:,cc),sem_trg_Vbb(:,cc));
    hold on
    shadedErrorBar(lags/Fsd,avg_ftrg_Vbb(:,cc),sem_ftrg_Vbb(:,cc),{'color','b'});
    shadedErrorBar(lags/Fsd,avg_strg_Vbb(:,cc),sem_strg_Vbb(:,cc),{'color','r'});
    xlim([-0.1 0.55])
    pause
    clf
end
