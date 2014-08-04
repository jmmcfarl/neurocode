clear all
% close all
cd ~/Data/bruce/G081/
load jbeG081Expts.mat
addpath('~/James_scripts/bruce/G081/');
%%
stim_fs = 100; %in Hz
Fs = 3e4;
dsf = 80;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[1 100]/niqf);
use_lfps = [1:8:96];

backlag = round(Fsd*0.1);
forwardlag = round(Fsd*0.4);
lags = -backlag:forwardlag;

%%
for i = 1:length(Expts)
    if strcmp(Expts{i}.Header.expname,'grating.OpXseRC') | strcmp(Expts{i}.Header.expname,'grating.OpRC')
        is_bar_expt(i) = 1;
    else
        is_bar_expt(i) = 0;
    end
    
    if strcmp(Expts{i}.Stimvals.Bs,'image')
        expt_image_back(i) = 1;
    else
        expt_image_back(i) = 0;
    end
    
    expt_sim_sacs(i) = Expts{i}.Stimvals.ijump;
    expt_bar_ori(i) = Expts{i}.Stimvals.or;
    
end
expt_bar_ori(expt_bar_ori == -45) = 135;

%%
load ./jbeG081.em.mat

%%
% bar_oris = [0 45 90 135];
bar_oris = [45];
fprintf('Analyzing %d ori expts\n',bar_oris);

%         cur_expt_set = find(is_bar_expt == 1 & expt_image_back == 1 & expt_bar_ori == bar_oris & expt_sim_sacs == 0);
cur_expt_set = find(is_bar_expt == 1 & expt_image_back == 0 & expt_bar_ori == bar_oris);

cur_expt_set(cur_expt_set<= 8) = []; %this expt had continuous bar luminance
%
all_stim_times = [];
all_rel_stime = [];
all_rel_etime = [];
all_Op = [];
all_phase = [];
all_trial_start_times = [];
all_trial_end_times = [];
all_trial_durs = [];
all_trial_expt = [];
for ee = 1:length(cur_expt_set);
    fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
    cur_expt = cur_expt_set(ee);
    
    n_trials = length(Expts{cur_expt}.Trials);
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_expt}.Trials(tt).Start/1e4;
        cur_Op = Expts{cur_expt}.Trials(tt).Op;
        cur_phase = Expts{cur_expt}.Trials(tt).ph;
        
        all_stim_times = [all_stim_times; cur_stim_times];
        all_rel_stime = [all_rel_stime; cur_stim_times - trial_start_times(tt)];
        all_rel_etime = [all_rel_etime; trial_end_times(tt) - cur_stim_times];
        all_Op = [all_Op; cur_Op];
        all_phase = [all_phase; cur_phase'];
    end
    
    all_trial_start_times = [all_trial_start_times; trial_start_times'];
    all_trial_end_times = [all_trial_end_times; trial_end_times'];
    all_trial_durs = [all_trial_durs; trial_durs'];
    all_trial_expt = [all_trial_expt; ee*ones(length(trial_durs),1)];
end
%%
all_Vmat = [];
all_tax = [];
all_outgoing_sacs = [];
all_returning_sacs = [];
all_stim_on = [];
all_stim_off = [];
all_msacs = [];
for ee = 1:length(cur_expt_set);
    cur_expt = cur_expt_set(ee);
    Vmat = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d of %d\n',ll,length(use_lfps));
        filename = sprintf('Expt%d.p%dFullV.mat',cur_expt,use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        %     V = V + FullV.sumscale*sumv;
        V = V*FullV.intscale(1)/FullV.intscale(2);
        nparts = length(FullV.blklen);
        dV = [];
        %splice together multiple blocks
        cur_pt = 1;
        for pp = 1:nparts
            cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
            cur_range(cur_range > length(V)) = [];
            curV = decimate(V(cur_range),dsf);
            curV = filtfilt(filt_b,filt_a,curV);
            dV = [dV curV];
            cur_pt = cur_pt + FullV.blklen(pp);
        end
        Vmat(:,ll) = dV;
    end
    
    t_ax = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        t_ax = [t_ax downsample(cur_t_ax,dsf)];
    end
    t_ax(size(Vmat,1)+1:end) = [];
    
    fprintf('LFP len: %d\n',range(t_ax));
    
    %% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)
    [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data_new(Expt,t_ax(([1 end])));
    
    eye_dt = median(diff(eye_ts_interp));
    eye_fs = 1/eye_dt;
    lEyeXY = eye_vals_interp(:,1:2);
    clear sm_avg_eyepos eye_vel
    sm_avg_eyepos(:,1) = smooth(lEyeXY(:,1),3);
    sm_avg_eyepos(:,2) = smooth(lEyeXY(:,2),3);
    eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/eye_dt;
    eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/eye_dt;
    
    eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
    
    interp_eye_vals = interp1(eye_ts_interp,lEyeXY,t_ax);
    interp_eye_speed = interp1(eye_ts_interp,eye_speed,t_ax);
    
    orth_eye_proj = (cos(bar_oris/180*pi)*interp_eye_vals(:,1) - ...
        sin(bar_oris/180*pi)*interp_eye_vals(:,2));
    
    %find threshold crossings
    sac_thresh = 25;
    saccade_inds = 1 + find(interp_eye_speed((1:end-1)) < sac_thresh & interp_eye_speed((2:end)) > sac_thresh);
    saccade_inds = unique(saccade_inds);
    
    isi = [0 diff(saccade_inds)]/Fsd;
    min_isi = 0.05;
    bad_sacs = find(isi < min_isi);
    saccade_inds(bad_sacs) = [];
    
    peri_thresh = 3;
    sac_starts = nan(size(saccade_inds));
    sac_stops = nan(size(saccade_inds));
    for i = 1:length(saccade_inds)
        cur_up = find(interp_eye_speed(1:saccade_inds(i)) < peri_thresh,1,'last');
        cur_down = find(interp_eye_speed(saccade_inds(i):end) < peri_thresh,1,'first');
        if ~isempty(cur_up) & ~isempty(cur_down)
            sac_starts(i) = cur_up;
            sac_stops(i) = saccade_inds(i) + cur_down;
        end
    end
    bad = find(isnan(sac_starts));
    saccade_inds(bad) = []; sac_starts(bad) = []; sac_stops(bad) = [];
    
    sac_starts = sac_starts - round(Fsd*0.01);
    sac_stops = sac_stops + round(Fsd*0.01);
    pre_pos = orth_eye_proj(sac_starts);
    post_pos = orth_eye_proj(sac_stops);
    sac_dist = abs(post_pos-pre_pos);
    
    cur_trial_set = find(all_trial_expt==ee);
    n_trials = length(cur_trial_set);
    in_trial = zeros(size(t_ax));
    for tt = 1:n_trials
        cur_set = find(t_ax >= all_trial_start_times(cur_trial_set(tt)) & t_ax < all_trial_end_times(cur_trial_set(tt)));
        in_trial(cur_set) = 1;
    end
    use_saccade_inds = find(ismember(saccade_inds,find(in_trial==1)));
    
    interp_rel_stime = interp1(all_stim_times,all_rel_stime,t_ax);
    %         interp_rel_etime = interp1(all_stim_times,all_rel_etime,t_ax);
    
    outgoing_sacs = use_saccade_inds(pre_pos(use_saccade_inds) > -0.6 & post_pos(use_saccade_inds) < -1);
    returning_sacs = use_saccade_inds(pre_pos(use_saccade_inds) < -1 & post_pos(use_saccade_inds) > -0.6);
    outgoing_sacs = saccade_inds(outgoing_sacs);
    returning_sacs = saccade_inds(returning_sacs);
    msacs = saccade_inds(use_saccade_inds(sac_dist(use_saccade_inds)' < 1 & ~ismember(use_saccade_inds,outgoing_sacs) & ...
        ~ismember(use_saccade_inds,returning_sacs)));
    
    outgoing_sacs(outgoing_sacs < backlag | outgoing_sacs > size(Vmat,1) - forwardlag) = [];
    returning_sacs(returning_sacs < backlag | returning_sacs > size(Vmat,1) - forwardlag) = [];
    msacs(msacs < backlag | msacs > size(Vmat,1) - forwardlag) = [];
    
    stim_on = 1+find(interp_rel_stime(1:end-1) < 0.7 & interp_rel_stime(2:end) > 0.7);
    stim_off = 1+find(interp_rel_stime(1:end-1) < 1.4 & interp_rel_stime(2:end) > 1.4);
    
    all_stim_on = [all_stim_on; stim_on' + length(all_tax)];
    all_stim_off = [all_stim_off; stim_off' + length(all_tax)];
    all_outgoing_sacs = [all_outgoing_sacs; outgoing_sacs'+length(all_tax)];
    all_returning_sacs = [all_returning_sacs; returning_sacs' + length(all_tax)];
    all_msacs = [all_msacs; msacs' + length(all_tax)];
    
    all_Vmat = [all_Vmat; Vmat];
    all_tax = [all_tax; t_ax'];
    
end

%%
interp_rel_stime = interp1(all_stim_times,all_rel_stime,all_tax);
interp_rel_etime = interp1(all_stim_times,all_rel_etime,all_tax);

trial_lfp_inds = round(interp1(all_tax,1:length(all_tax),all_trial_start_times));
%     stim_lfp_inds = round(interp1(all_tax,1:length(all_tax),all_stim_times));
n_trials = length(trial_lfp_inds);

%     use_stims = find(all_rel_stime > backlag & all_rel_etime > forwardlag/Fsd);
use_outgoing_inds = all_outgoing_sacs(interp_rel_stime(all_outgoing_sacs) > 0.25 & interp_rel_etime(all_outgoing_sacs) > 0.15);
use_returning_inds = all_returning_sacs(interp_rel_stime(all_returning_sacs) > 0.25 & interp_rel_etime(all_returning_sacs) > 0.15);
use_msac_inds = all_msacs(interp_rel_stime(all_msacs) > 0.25 & interp_rel_etime(all_msacs) > forwardlag/Fsd);

use_stimon_inds = all_stim_on(interp_rel_stime(all_stim_on) > 0.25 & interp_rel_etime(all_stim_on) > 0.15);
use_stimoff_inds = all_stim_off(interp_rel_stime(all_stim_off) > 0.25 & interp_rel_etime(all_stim_off) > 0.15);

on_trig_mat = zeros(length(use_stimon_inds),length(lags),length(use_lfps));
for ii = 1:length(use_stimon_inds)
    cur_inds = (use_stimon_inds(ii) - backlag):(use_stimon_inds(ii) + forwardlag);
    on_trig_mat(ii,:,:) = all_Vmat(cur_inds,:);
end
on_trig_avgs = squeeze(mean(on_trig_mat));
on_trig_sem = squeeze(std(on_trig_mat))/sqrt(length(use_stimon_inds));

off_trig_mat = zeros(length(use_stimoff_inds),length(lags),length(use_lfps));
for ii = 1:length(use_stimoff_inds)
    cur_inds = (use_stimoff_inds(ii) - backlag):(use_stimoff_inds(ii) + forwardlag);
    off_trig_mat(ii,:,:) = all_Vmat(cur_inds,:);
end
off_trig_avgs = squeeze(mean(off_trig_mat));
off_trig_sem = squeeze(std(off_trig_mat))/sqrt(length(use_stimoff_inds));

outsac_trig_mat = zeros(length(use_outgoing_inds),length(lags),length(use_lfps));
for ii = 1:length(use_outgoing_inds)
    cur_inds = (use_outgoing_inds(ii) - backlag):(use_outgoing_inds(ii) + forwardlag);
    outsac_trig_mat(ii,:,:) = all_Vmat(cur_inds,:);
end
outsac_trig_avgs = squeeze(mean(outsac_trig_mat));
outsac_trig_sem = squeeze(std(outsac_trig_mat))/sqrt(length(use_outgoing_inds));

retsac_trig_mat = zeros(length(use_returning_inds),length(lags),length(use_lfps));
for ii = 1:length(use_returning_inds)
    cur_inds = (use_returning_inds(ii) - backlag):(use_returning_inds(ii) + forwardlag);
    retsac_trig_mat(ii,:,:) = all_Vmat(cur_inds,:);
end
retsac_trig_avgs = squeeze(mean(retsac_trig_mat));
retsac_trig_sem = squeeze(std(retsac_trig_mat))/sqrt(length(use_returning_inds));

msac_trig_mat = zeros(length(use_msac_inds),length(lags),length(use_lfps));
for ii = 1:length(use_msac_inds)
    cur_inds = (use_msac_inds(ii) - backlag):(use_msac_inds(ii) + forwardlag);
    msac_trig_mat(ii,:,:) = all_Vmat(cur_inds,:);
end
msac_trig_avgs = squeeze(mean(msac_trig_mat));
msac_trig_sem = squeeze(std(msac_trig_mat))/sqrt(length(use_msac_inds));

fmaxlag = round(Fsd*2);
flags = (0:round(2*Fsd));
full_trig_mat = zeros(n_trials,length(flags),length(use_lfps));
for ii = 1:n_trials
    cur_inds = (trial_lfp_inds(ii)):(trial_lfp_inds(ii) + fmaxlag);
    full_trig_mat(ii,:,:) = all_Vmat(cur_inds,:);
end
full_trig_avgs = squeeze(mean(full_trig_mat));
full_trig_sem = squeeze(std(full_trig_mat))/sqrt(n_trials);


sac_t_ax = linspace(0,2,200);
all_sac_inds = [all_msacs; all_outgoing_sacs; all_returning_sacs];
sac_dist = hist(interp_rel_stime(all_sac_inds),sac_t_ax);

out_dist = hist(interp_rel_stime(use_outgoing_inds),sac_t_ax);
ret_dist = hist(interp_rel_stime(use_returning_inds),sac_t_ax);
msac_dist = hist(interp_rel_stime(use_msac_inds),sac_t_ax);

%%
% NT = size(all_Vmat,1);
% on_inds = zeros(NT,1);
% off_inds = zeros(NT,1);
% outsac_inds = zeros(NT,1);
% retsac_inds = zeros(NT,1);
% msac_inds = zeros(NT,1);
% on_inds(use_stimon_inds) = 1;
% off_inds(use_stimoff_inds) = 1;
% outsac_inds(use_outgoing_inds) =1;
% retsac_inds(use_returning_inds) = 1;
% msac_inds(use_msac_inds) = 1;
% 
% flen = round(Fsd*0.4);
% Xmat_on = makeStimRows(on_inds,flen);
% Xmat_off = makeStimRows(off_inds,flen);
% Xmat_out = makeStimRows(outsac_inds,flen);
% Xmat_ret = makeStimRows(retsac_inds,flen);
% Xmat_msac = makeStimRows(msac_inds,flen);
% clear *_kern
% for ll = 1:length(use_lfps)
%     ll
%     temp = regress(all_Vmat(:,ll),[Xmat_on Xmat_out Xmat_off Xmat_ret Xmat_msac]);
%     on_kern(ll,:) = flipud(temp(1:flen));
%     out_kern(ll,:) = flipud(temp(flen+1:2*flen));
%     off_kern(ll,:) = flipud(temp(2*flen+1:3*flen));
%     ret_kern(ll,:) = flipud(temp(3*flen+1:4*flen));
%     msac_kern(ll,:) = flipud(temp(4*flen+1:5*flen));
% end
%%
% save event_trg_avg_lfps_45deg_gray lags Fsd flags flen sac_t_ax *_dist *_trig_avgs* *_kern 
save event_trg_avg_lfps_45deg_im lags Fsd flags sac_t_ax *_dist *_trig_avgs* 
%%
close all
load ./event_trg_avg_lfps_0deg_gray
gray0_out_avgs = outsac_trig_avgs;
gray0_ret_avgs = retsac_trig_avgs;
gray0_msac_avgs = msac_trig_avgs;
gray0_on_avgs = on_trig_avgs;
gray0_off_avgs = off_trig_avgs;

load ./event_trg_avg_lfps_0deg_im.mat
im0_out_avgs = outsac_trig_avgs;
im0_ret_avgs = retsac_trig_avgs;
im0_msac_avgs = msac_trig_avgs;
im0_on_avgs = on_trig_avgs;
im0_off_avgs = off_trig_avgs;

load ./event_trg_avg_lfps_90deg_gray
gray90_out_avgs = outsac_trig_avgs;
gray90_ret_avgs = retsac_trig_avgs;
gray90_msac_avgs = msac_trig_avgs;
gray90_on_avgs = on_trig_avgs;
gray90_off_avgs = off_trig_avgs;

load ./event_trg_avg_lfps_90deg_im.mat
im90_out_avgs = outsac_trig_avgs;
im90_ret_avgs = retsac_trig_avgs;
im90_msac_avgs = msac_trig_avgs;
im90_on_avgs = on_trig_avgs;
im90_off_avgs = off_trig_avgs;

for ll = 1:length(use_lfps)
    subplot(2,2,1)
    plot(lags/Fsd,gray0_out_avgs(:,ll),'b')
    hold on
    plot(lags/Fsd,gray0_ret_avgs(:,ll),'r')
    plot(lags/Fsd,gray0_msac_avgs(:,ll),'g')
    axis tight
    subplot(2,2,3)
    plot(lags/Fsd,gray90_out_avgs(:,ll),'b')
    hold on
    plot(lags/Fsd,gray90_ret_avgs(:,ll),'r')
    plot(lags/Fsd,gray90_msac_avgs(:,ll),'g')
    axis tight

        subplot(2,2,2)
    plot(lags/Fsd,im0_out_avgs(:,ll),'b')
    hold on
    plot(lags/Fsd,im0_ret_avgs(:,ll),'r')
    plot(lags/Fsd,im0_msac_avgs(:,ll),'g')
    axis tight
    subplot(2,2,4)
    plot(lags/Fsd,im90_out_avgs(:,ll),'b')
    hold on
    plot(lags/Fsd,im90_ret_avgs(:,ll),'r')
    plot(lags/Fsd,im90_msac_avgs(:,ll),'g')
    axis tight

%     subplot(2,2,1)
%     plot(lags/Fsd,gray0_on_avgs(:,ll),'b')
%     hold on
%     plot(lags/Fsd,gray0_off_avgs(:,ll),'r')
%     axis tight
%     subplot(2,2,3)
%     plot(lags/Fsd,gray90_on_avgs(:,ll),'b')
%     hold on
%     plot(lags/Fsd,gray90_off_avgs(:,ll),'r')
%     axis tight
% 
%         subplot(2,2,2)
%     plot(lags/Fsd,im0_on_avgs(:,ll),'b')
%     hold on
%     plot(lags/Fsd,im0_off_avgs(:,ll),'r')
%     axis tight
%     subplot(2,2,4)
%     plot(lags/Fsd,im90_on_avgs(:,ll),'b')
%     hold on
%     plot(lags/Fsd,im90_off_avgs(:,ll),'r')
%     axis tight

    pause
    clf
end

%%
close all
load ./event_trg_avg_lfps_0deg_gray
gray_out_avgs = outsac_trig_avgs;
gray_ret_avgs = retsac_trig_avgs;
gray_msac_avgs = msac_trig_avgs;
gray_on_avgs = on_trig_avgs;
gray_off_avgs = off_trig_avgs;

load ./event_trg_avg_lfps_0deg_im.mat
im_out_avgs = outsac_trig_avgs;
im_ret_avgs = retsac_trig_avgs;
im_msac_avgs = msac_trig_avgs;
im_on_avgs = on_trig_avgs;
im_off_avgs = off_trig_avgs;

for ll = 1:length(use_lfps)
    subplot(2,2,1)
    plot(lags/Fsd,gray_out_avgs(:,ll),'b')
    hold on
    plot(lags/Fsd,gray_ret_avgs(:,ll),'r')
    plot(lags/Fsd,gray_msac_avgs(:,ll),'g')
    axis tight
    subplot(2,2,3)
    plot(lags/Fsd,gray_on_avgs(:,ll),'b')
    hold on
    plot(lags/Fsd,gray_off_avgs(:,ll),'r')
    axis tight

        subplot(2,2,2)
    plot(lags/Fsd,im_out_avgs(:,ll),'b')
    hold on
    plot(lags/Fsd,im_ret_avgs(:,ll),'r')
    plot(lags/Fsd,im_msac_avgs(:,ll),'g')
    axis tight
    subplot(2,2,4)
    plot(lags/Fsd,im_on_avgs(:,ll),'b')
    hold on
    plot(lags/Fsd,im_off_avgs(:,ll),'r')
    axis tight

    pause
    clf
end
%%
for ll = 1:length(use_lfps)
    %         subplot(2,2,1)
    %         plot(lags/Fsd,outsac_trig_avgs(:,ll))
    %         hold on
    %         plot(lags/Fsd,retsac_trig_avgs(:,ll),'r')
    %         plot(lags/Fsd,msac_trig_avgs(:,ll),'k'); axis tight
    %         subplot(2,2,3)
    %         plot(lags/Fsd,on_trig_avgs(:,ll))
    %         hold on
    %         plot(lags/Fsd,off_trig_avgs(:,ll),'r')
    %         axis tight
    %         subplot(2,2,2)
    %         plot(flags/Fsd,full_trig_avgs(:,ll))
    %         subplot(2,2,4)
    %         plot(sac_t_ax,sac_dist,'k')
    %         hold on
    %         plot(sac_t_ax,out_dist,'b',sac_t_ax,ret_dist,'r')
    %
    
    subplot(2,3,1)
    errorbar(lags/Fsd,outsac_trig_avgs(:,ll),outsac_trig_sem(:,ll))
    hold on
    errorbar(lags/Fsd,retsac_trig_avgs(:,ll),retsac_trig_sem(:,ll),'r')
    errorbar(lags/Fsd,msac_trig_avgs(:,ll),msac_trig_sem(:,ll),'g'); axis tight
    xlim([0 0.4])
    subplot(2,3,4)
    errorbar(lags/Fsd,on_trig_avgs(:,ll),on_trig_sem(:,ll))
    hold on
    errorbar(lags/Fsd,off_trig_avgs(:,ll),off_trig_sem(:,ll),'r')
    axis tight
    xlim([0 0.4])
    
    subplot(2,3,2)
    plot((1:flen)/Fsd,out_kern(ll,:),'b')
    hold on
    plot((1:flen)/Fsd,ret_kern(ll,:),'r')
    plot((1:flen)/Fsd,msac_kern(ll,:),'g')
    axis tight
    subplot(2,3,5)
    plot((1:flen)/Fsd,on_kern(ll,:),'b')
    hold on
    plot((1:flen)/Fsd,off_kern(ll,:),'r')
    axis tight
    
    subplot(2,3,3)
    plot(flags/Fsd,full_trig_avgs(:,ll))
    subplot(2,3,6)
    plot(sac_t_ax,sac_dist,'k')
    hold on
    plot(sac_t_ax,out_dist,'b',sac_t_ax,ret_dist,'r',sac_t_ax,msac_dist,'g')
    
    pause
    clf
end

