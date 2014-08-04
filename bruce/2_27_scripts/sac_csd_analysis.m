clear all
% close all
addpath(genpath('~/James_scripts'));

cd ~/Data/bruce/2_27_12/
load Blocks.mat
accept_window = [-5 5;-5 5];
sac_eyespeed = 10;
thresh_eyespeed = 2.5;
sac_buffer = 0;
min_fix_dur = 0.2;
blink_thresh = 5;
min_blink_dur = 0.05;
% nlags = 4;

Pix2Deg = 0.018837;

cd ~/Data/bruce/2_27_12/stimrecon
load fixation_data_v4

cd ~/Data/bruce/2_27_12/stimrecon
load sing_eye_pgabortrack_fin_est_d2
cd ~/Data/bruce/2_27_12/ExptA

%%
Fs = 1000;
dsf = 4;
Fsd = Fs/dsf;
niqf = Fs/2;
niqfd = Fsd/2;
[b,a] = butter(2,[1 50]/niqf);
% [b2,a2] = butter(2,10/niqfd,'low');

backlag = round(0.2*Fsd);
forwardlag = round(0.4*Fsd);
lags = -backlag:forwardlag;

%%
all_lfp_mat = [];
all_fix_csd_mat = [];
all_sac_csd_mat = [];
all_used_fixs = [];
vall_sac_amps = [];
for blockid = 1:3;
    fprintf('Block %d of %d\n',blockid,3);
    
    %%
    cd ~/Data/bruce/2_27_12
    load(sprintf('lemM232A.5%d.lfp.mat',blockid));
    %get start times of each LFP trial
    n_lfp_trials = length(LFP.Trials);
    lfp_trial_start = nan(n_lfp_trials,1);
    lfp_trial_stop = nan(n_lfp_trials,1);
    for i = 1:n_lfp_trials
        lfp_trial_start(i) = LFP.Trials(i).Start/1e4; %originally in tenths of ms
        lfp_trial_stop(i) = LFP.Trials(i).End/1e4;
        lfp_dur(i) = size(LFP.Trials(i).LFP,1)/1000;
    end
    
    lfp_time = [];
    lfp_samps = [];
    for i = 1:n_lfp_trials
        lfp_time = [lfp_time linspace(lfp_trial_start(i),lfp_trial_start(i)+lfp_dur(i),size(LFP.Trials(i).LFP,1))];
        lfp_samps = [lfp_samps; LFP.Trials(i).LFP];
    end
    
    lfp_sampsd = filtfilt(b,a,lfp_samps);
    lfp_sampsd = downsample(lfp_sampsd,dsf);
    lfp_timed = downsample(lfp_time,dsf);
%     lfp_sampsd_sm = lfp_sampsd;
%     lfp_sampsd_sm(:,2:end-1) = 0.23*lfp_sampsd(:,1:end-2) + 0.54*lfp_sampsd(:,2:end-1) + 0.23*lfp_sampsd(:,3:end);
%     lfp_sampsd_sm(:,1) = 0.54*lfp_sampsd(:,1) + 0.46*lfp_sampsd(:,2);
%     lfp_sampsd_sm(:,end) = 0.54*lfp_sampsd(:,end) + 0.46*lfp_sampsd(:,end-1);
%     cscd = 2*lfp_sampsd_sm(:,2:end-1) - lfp_sampsd_sm(:,1:end-2) - lfp_sampsd_sm(:,3:end);
%     cscd = filtfilt(b2,a2,cscd);
%     cscd = zscore(cscd);
    %%
    all_sac_start_times = all_fix_start_times - all_sac_durs;
    
    cur_fix_set = find(blockids == blockid);
    cur_fix_start_inds = round(interp1(lfp_timed,1:length(lfp_timed),all_fix_start_times(cur_fix_set)));
    cur_fix_stop_inds = round(interp1(lfp_timed,1:length(lfp_timed),all_fix_stop_times(cur_fix_set)));
    cur_fix_stop_inds(isnan(cur_fix_stop_inds)) = length(lfp_timed);
    cur_sac_start_inds = round(interp1(lfp_timed,1:length(lfp_timed),all_sac_start_times(cur_fix_set)));
   
    bad = find(cur_sac_start_inds < backlag | cur_fix_stop_inds > length(lfp_timed) - forwardlag);
    cur_fix_start_inds(bad) = [];
    cur_fix_stop_inds(bad) = [];
    cur_fix_set(bad) = [];
    cur_sac_start_inds(bad) = [];
    
    cur_fix_durs = cur_fix_stop_inds - cur_fix_start_inds;
    cur_n_fixs = length(cur_fix_start_inds);
        
    
    fix_trig_lfp_mat = nan(24,cur_n_fixs,length(lags));
    sac_trig_lfp_mat = nan(24,cur_n_fixs,length(lags));
    for i = 1:cur_n_fixs
        cur_set = (cur_fix_start_inds(i)-backlag):(cur_fix_start_inds(i)+forwardlag);
        if cur_fix_durs(i) < forwardlag
            cur_set(backlag+cur_fix_durs(i):end) = [];
        end
        fix_trig_lfp_mat(:,i,1:length(cur_set)) = lfp_sampsd(cur_set,:)';        
         
        cur_set = (cur_sac_start_inds(i)-backlag):(cur_sac_start_inds(i)+forwardlag);
        if cur_fix_durs(i) < forwardlag
            cur_set(backlag+cur_fix_durs(i):end) = [];
        end
        sac_trig_lfp_mat(:,i,1:length(cur_set)) = lfp_sampsd(cur_set,:)';        
   end
    
    %%
    vars.Fs = Fsd;
    vars.BrainBound = 1;
    vars.ChanSep = 0.05;
    vars.useVaknin = 'false';
    vars.useHamming = 'true';
    vars.diam = 2;
    Data = permute(fix_trig_lfp_mat,[1 3 2]);
    fix_CSD = PettersenCSD(Data,'spline',vars);
    fix_CSD = permute(fix_CSD,[1 3 2]);
    
    Data = permute(sac_trig_lfp_mat,[1 3 2]);
    sac_CSD = PettersenCSD(Data,'spline',vars);
    sac_CSD = permute(sac_CSD,[1 3 2]);

    %%
    
    cur_sac_amps = all_sac_amps(cur_fix_set);
    cur_sac_durs = all_sac_durs(cur_fix_set);
    cur_sac_dirs = atan2(all_sac_dy(cur_fix_set),all_sac_dx(cur_fix_set));
    
    %%
    vall_sac_amps = [vall_sac_amps; cur_sac_amps];
    
    all_lfp_mat = cat(2,all_lfp_mat,sac_trig_lfp_mat);
    all_fix_csd_mat = cat(2,all_fix_csd_mat,fix_CSD);
    all_sac_csd_mat = cat(2,all_sac_csd_mat,sac_CSD);
    
    all_used_fixs = [all_used_fixs; cur_fix_set];
    
end

%%
small_sacs = find(vall_sac_amps < 1);
big_sacs = find(vall_sac_amps > 1);

msac_CSD_avg = squeeze(nanmean(all_sac_csd_mat(:,small_sacs,:),2));
bsac_CSD_avg = squeeze(nanmean(all_sac_csd_mat(:,big_sacs,:),2));
sac_CSD_avg = squeeze(nanmean(all_sac_csd_mat,2));

fix_CSD_avg = squeeze(nanmean(all_fix_csd_mat,2));
mfix_CSD_avg = squeeze(nanmean(all_fix_csd_mat(:,small_sacs,:),2));
bfix_CSD_avg = squeeze(nanmean(all_fix_csd_mat(:,big_sacs,:),2));

%%
figure
imagesc(lags/Fsd,(1:24)*0.05,sac_CSD_avg);
xlim([-0.2 0.4])
xlabel('Time lag (s)','fontsize',16)
ylabel('Depth (um)','fontsize',16)
colorbar
colormap(colormap_redblackblue);
ca = max(abs(caxis()));
caxis([-ca ca]*0.95);
yl = ylim();
line([0 0],yl,'color','w');

figure
imagesc(lags/Fsd,(1:24)*0.05,fix_CSD_avg);
xlim([-0.2 0.4])
xlabel('Time lag (s)','fontsize',16)
ylabel('Depth (um)','fontsize',16)
colorbar
colormap(colormap_redblackblue);
ca = max(abs(caxis()));
caxis([-ca ca]*0.95);
yl = ylim();
line([0 0],yl,'color','w');

%%
imagesc(lags/Fsd,(1:24)*0.05,msac_CSD_avg);
xlim([-0.2 0.4])
xlabel('Time lag (s)','fontsize',16)
ylabel('Depth (um)','fontsize',16)
colorbar
colormap(colormap_redblackblue);
ca = max(abs(caxis()));
caxis([-ca ca]*0.95);
yl = ylim();
line([0 0],yl,'color','w');

figure
imagesc(lags/Fsd,(1:24)*0.05,bsac_CSD_avg);
xlim([-0.2 0.4])
xlabel('Time lag (s)','fontsize',16)
ylabel('Depth (um)','fontsize',16)
colorbar
colormap(colormap_redblackblue);
% ca = max(abs(caxis()));
caxis([-ca ca]*0.95);
yl = ylim();
line([0 0],yl,'color','w');



