clear all
% close all

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

% down-sampling fraction for image
dsfrac = 2;
Fsd = 1/Pix2Deg/dsfrac;

Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
% [XX,YY] = meshgrid(xax,yax);

% RF_patch = [3.5 6; -3.5 -1]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [2.5 6.5; -4 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [2.75 6.25; -3.5 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

% patch_inds = find(XX >= RF_patch(1,1) & XX <= RF_patch(1,2) & YY >= RF_patch(2,1) & YY <= RF_patch(2,2));
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_data_v3

cd /Users/James/Data/bruce/2_27_12/stimrecon
load sing_eye_pgabortrack_d2
load sing_eye_pgabortrack_fin_est_d2
cd ~/Data/bruce/2_27_12/ExptA

SDIM = length(xpatch_inds);
cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_image_patches_corrected
cd /Users/James/James_scripts/bruce/modelfits/
% load fixation_gabor_models
load gabor_tempmodfits_en_lin

cd /Users/James/Data/bruce/2_27_12/stimrecon
load sing_eye_pgabortrack_d2
load sing_eye_pgabortrack_fin_est_d2


NT = length(used_fixs);
X_resh = reshape(X,NT,SDIM^2);

%%
clear gmask_out total_out all_loc*
for c = 1:18
cur_mask1 = get_pgabor_mask(gabor_params_fin(c,1:6),0,[SDIM SDIM]);
cur_mask2 = get_pgabor_mask(gabor_params_fin(c,1:6),pi/2,[SDIM SDIM]);
cur_gmask = get_pgauss_mask(gabor_params_fin(c,1:6),[SDIM SDIM]);
cur_gmask = cur_gmask(:);
mask1_out = X_resh*cur_mask1(:);
mask2_out = X_resh*cur_mask2(:);
gmask_out = bsxfun(@times,X_resh,cur_gmask');
loc_mean = mean(gmask_out,2);
loc_cont = std(gmask_out,[],2);

% mask1_out = (mask1_out - loc_mean)./loc_cont;
% mask2_out = (mask2_out - loc_mean)./loc_cont;

all_loc_mean(c,:) = zscore(loc_mean(all_model_fixids));
all_loc_cont(c,:) = zscore(loc_cont(all_model_fixids));

lin_out = gabor_params_fin(c,8)*mask1_out(all_model_fixids) + gabor_params_fin(c,9)*mask2_out(all_model_fixids);
energy_out = gabor_params_fin(c,7)*sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
total_out(c,:) = lin_out + energy_out;
total_out(c,:) = zscore(total_out(c,:));
end
%%
cd ~/Data/bruce/2_27_12

for blockid = 1:4;
    
    fprintf('Processing block %d of %d...\n',blockid,4);
    
    block_times = Blocks{blockid}.blocktimes;
    stim_times = Blocks{blockid}.stimtime;
    stimID = Blocks{blockid}.stimids;
    
    n_stims = length(stim_times);
    
    cd ~/Data/bruce/2_27_12/saccades/
    load(sprintf('lemM232.5%d.em.sac.mat',blockid))
    %     load(sprintf('lemM232.5%d.em.hor.mat',blockid))
    
    % identify saccade start and stop times
    EyeStartT = Expt.Trials.Start/10000; % time of first eye sample
    EyeEndT = Expt.Trials.End/10000; % time of last eye sample
    Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
    eyets = EyeStartT:Eyedt:EyeEndT; %eye tracking time axis (sec)
    sac_buffer_inds = round(sac_buffer/Eyedt);
    
    reye_pos = [Expt.Trials.Eyevals.rh Expt.Trials.Eyevals.rv];
    leye_pos = [Expt.Trials.Eyevals.lh Expt.Trials.Eyevals.lv];
    
    avg_eyepos = (reye_pos + leye_pos)/2;
    clear sm_avg_eyepos eye_vel
    sm_avg_eyepos(:,1) = smooth(avg_eyepos(:,1),3);
    sm_avg_eyepos(:,2) = smooth(avg_eyepos(:,2),3);
    eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/Eyedt;
    eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/Eyedt;
    eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
    vel_thresh = 5;
    eye_vel(eye_vel(:,1) > vel_thresh,1) = vel_thresh;
    eye_vel(eye_vel(:,1) < -vel_thresh,1) = -vel_thresh;
    eye_vel(eye_vel(:,2) > vel_thresh,2) = vel_thresh;
    eye_vel(eye_vel(:,2) < -vel_thresh,2) = -vel_thresh;
    
    %%
    cd /Users/James/Data/bruce/2_27_12
    load(sprintf('lemM232A.5%d.lfp.mat',blockid));
        %get start times of each LFP trial
    n_lfp_trials = length(LFP.Trials);
    lfp_trial_start = nan(n_lfp_trials,1);
    for i = 1:n_lfp_trials
        lfp_trial_start(i) = LFP.Trials(i).Start/1e4; %originally in tenths of ms
        lfp_trial_stop(i) = LFP.Trials(i).End/1e4;
        lfp_dur(i) = size(LFP.Trials(i).LFP,1)/1000;
    end
    % lfp_trial_stop = lfp_trial_start+lfp_dur';
    block_trial_times = Blocks{blockid}.blocktimes;
    
    lfp_time = [];
    lfp_samps = [];
    for i = 1:n_lfp_trials
        lfp_time = [lfp_time linspace(lfp_trial_start(i),lfp_trial_start(i)+lfp_dur(i),size(LFP.Trials(i).LFP,1))];
        lfp_samps = [lfp_samps; LFP.Trials(i).LFP];
    end
    % lfp_samps = zscore(lfp_samps);
    
    Fs = 1000;
    niqf = Fs/2;
    [b,a] = butter(2,[1 200]/niqf);
    lfp_samps = filtfilt(b,a,lfp_samps);
    [b,a] = butter(2,[4 15]/niqf);
    lfp_samps_alpha = filtfilt(b,a,lfp_samps);
    [b,a] = butter(2,[20 45]/niqf);
    lfp_samps_hf = filtfilt(b,a,lfp_samps);
    
    dsf = 1;
    lfp_sampsd = downsample(lfp_samps,dsf);
    lfp_sampsd = zscore(lfp_sampsd);
    lfp_sampsd_hf = downsample(lfp_samps_hf,dsf);
    lfp_sampsd_hf = zscore(lfp_sampsd_hf);
    lfp_sampsd_alpha = downsample(lfp_samps_alpha,dsf);
    lfp_sampsd_alpha = zscore(lfp_sampsd_alpha);
    lfp_timed = downsample(lfp_time,dsf);
    Fsd = Fs/dsf;
    
    lfp_sampsd_alpha_h = hilbert(lfp_sampsd_alpha);
    lfp_sampsd_alpha_phase = angle(lfp_sampsd_alpha_h);
    lfp_sampsd_alpha_phase = mod(lfp_sampsd_alpha_phase,2*pi);
    phase_revs = 1+find(diff(lfp_sampsd_alpha_phase(:,2)) < -pi);
    %%
    spikes_binned = zeros(24,length(lfp_timed));
    for c = 1:10
        c
        temp = round(interp1(lfp_timed,1:length(lfp_timed),Blocks{blockid}.spktimes{c}));
        temp(isnan(temp)) = [];
        all_spike_lfp_inds{c} = temp;
        
        cur_hist = hist(Blocks{blockid}.spktimes{c},lfp_timed);
        cur_hist(end) = 0;
        spikes_binned(c,:) = smooth(cur_hist,10);
    end
     for c = 1:14
         c
        temp = round(interp1(lfp_timed,1:length(lfp_timed),Blocks{blockid}.mutimes{c}));
        temp(isnan(temp)) = [];
        all_mua_lfp_inds{c} = temp;
        
        cur_hist = hist(Blocks{blockid}.mutimes{c},lfp_timed);
        cur_hist(end) = 0;
        spikes_binned(c+10,:) = smooth(cur_hist,10);
     end
     
    spike_rates = spikes_binned*Fsd;
    spike_rates = zscore(spike_rates')';
    net_rate = mean(spike_rates);
    
    suprobes = Blocks{1}.suprobes;
    muprobes = Blocks{1}.muprobes;
    allprobes = [suprobes muprobes];
%     set1 = find(allprobes < 10);
%     set2 = find(allprobes >= 10);
%     net_rate_1 = mean(spike_rates(set1,:));
%     net_rate_2 = mean(spike_rates(set2,:));
    
    %%
    all_sac_dirs = atan2(all_sac_dy,all_sac_dx)*180/pi;
%     blockid = 1;
    close all
    cur_lfp1 = 2;
    cur_lfp2 = 15;
    cur_lfp3 = 10;
    all_sus = [1 2 3 4 5 6 7 8 9 10];
    backwin = 0.25;
    forwin = 1;
    yl = [-6 4];
    cur_sus = 6;
    figure
    subplot(3,1,[1 2])
    plot(lfp_timed,lfp_sampsd(:,cur_lfp1)+2,'linewidth',1);
    hold on
    plot(lfp_timed,lfp_sampsd(:,cur_lfp2)+2,'r','linewidth',1);
%     plot(lfp_timed,lfp_sampsd(:,cur_lfp3)+2,'k','linewidth',1);
    plot(lfp_timed,lfp_sampsd_alpha(:,cur_lfp1)/2+2,'k')
    plot(lfp_timed,lfp_sampsd_alpha_phase(:,cur_lfp1)/4,'g')
    plot(lfp_timed(phase_revs),lfp_sampsd_alpha_phase(phase_revs,cur_lfp1)/4,'ro')
    
    plot(lfp_timed,lfp_sampsd_hf(:,cur_lfp2)/2-2,'c')
    
    cc = find(all_model_blockids==blockid);
     plot(lfp_timed,net_rate-4.5,'r','linewidth',1.5)
     plot(lfp_timed,spike_rates(6,:)-4.5,'k')
   ylim(yl);
%     for c = 1:length(lfp_sus)
%         plot(Blocks{blockid}.spktimes{lfp_sus(c)},-ones(size(Blocks{blockid}.spktimes{lfp_sus(c)}))*0.25*(c-1)-6.5,'.','color',cmap(c,:));
% %         plot(lfp_timed(all_spike_lfp_inds{lfp_sus(c)}),lfp_sampsd_hf(all_spike_lfp_inds{lfp_sus(c)},cur_lfp1)-3,'.','color',cmap(c,:));
%     end
%         plot(Blocks{blockid}.spktimes{cur_sus},-ones(size(Blocks{blockid}.spktimes{cur_sus}))*0.25*(cur_sus-1),'k.');

plot(all_model_time_axis(cc),total_out(1:3,cc),'.-')
hold on
plot(all_model_time_axis(cc),mean(all_loc_cont(:,cc)),'k.-','linewidth',2)
subplot(3,1,3)
    hold on
    plot(eyets(1:end-1),log(eye_speed),'k')
    plot(eyets(1:end-1),eye_vel(:,1),'r')
    plot(eyets(1:end-1),eye_vel(:,2),'b')

     temp_fix_set = find(blockids == blockid & all_fix_start_times > lfp_timed(1));
   for i = 1:length(temp_fix_set)
         start_T = all_fix_start_times(temp_fix_set(i));
        end_T = all_fix_stop_times(temp_fix_set(i));
        sacend_T = start_T-all_sac_durs(temp_fix_set(i));
      subplot(3,1,[1 2])
        line([start_T start_T],yl,'color','r')
        line([sacend_T sacend_T],yl,'color','g')
    subplot(3,1,3)
        line([start_T start_T],yl,'color','r')
        line([sacend_T sacend_T],yl,'color','g')
    end
    
    cur_fix_set = find(blockids(used_fixs) == blockid & all_fix_start_times(used_fixs) > lfp_timed(1));
    for i = 1:length(cur_fix_set)
        start_T = all_fix_start_times(used_fixs(cur_fix_set(i)));
        end_T = all_fix_stop_times(used_fixs(cur_fix_set(i)));
        sacend_T = start_T-all_sac_durs(used_fixs(cur_fix_set(i)));
    subplot(3,1,[1 2])
        title(sprintf('Fixation Number %d',cur_fix_set(i)));
        xlim([start_T-backwin start_T+forwin])
%         line([start_T start_T],yl,'color','r')
%         line([sacend_T sacend_T],yl,'color','g')
    subplot(3,1,3)
        title(sprintf('Sac Amp: %.3f   Sac Dir: %.3f',all_sac_amps(used_fixs(cur_fix_set(i))),all_sac_dirs(used_fixs(cur_fix_set(i)))));
        xlim([start_T-backwin start_T+forwin])
%         line([start_T start_T],yl,'color','r')
%         line([sacend_T sacend_T],yl,'color','g')
        input('');
    end

    %%
    
end








