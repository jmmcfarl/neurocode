clear all
% close all
cd ~/Data/bruce/G081/
load jbeG081Expts.mat
load ./CellList.mat
good_sus = find(all(CellList(:,:,1) > 0));

stim_fs = 100; %in Hz
Fs = 3e4;
dsf = 100;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[1]/niqf,'high');
use_lfps = [1:4:96];

use_units = 1:96;

% scales = logspace(log10(0.2),log10(4.5),25);
scales = logspace(log10(5),log10(200),40);
scales = scales*60/dsf;
% scales = scales*2;
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

beg_buffer = round(Fsd*0.2); %don't use the first X data after start of a trial.
end_buffer = round(Fsd*0.2);

%% PARSE TRIAL DATA STRUCTURES
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

load ./all_un_bar_pos
n_bar_pos = size(all_un_bar_pos,1);


%% DETERMINE SET OF EXPERIMENTS TO USE
flen = 12;
bar_oris = [0];
un_bar_pos = all_un_bar_pos(:,1);

fprintf('Analyzing %d ori expts\n',bar_oris);

%     cur_expt_set = find(is_bar_expt == 1 & expt_image_back == 1 & expt_bar_ori == bar_oris(EE));

% %expts with X deg bars and gray back (sim sacs)
% cur_expt_set = find(is_bar_expt==1 & expt_image_back == 0 & expt_bar_ori == bar_oris);

%expts with X deg bars and any back (including sim sacs)
% cur_expt_set = find(is_bar_expt==1 & expt_bar_ori == bar_oris);

% cur_expt_set = find(is_bar_expt==1 & expt_image_back == 0);
% cur_expt_set = find(is_bar_expt==1 & expt_image_back == 1 & expt_sim_sacs == 0);
cur_expt_set = find(is_bar_expt==1 & expt_image_back == 1 & expt_sim_sacs > 0);

cur_expt_set(cur_expt_set<=8) = []; %this expt had continuous bar luminance
cur_expt_set(cur_expt_set == 13) = []; %problem with LFP data?
% cur_expt_set(cur_expt_set >= 46 & cur_expt_set <= 51) = []; %problem with image background
cur_expt_set(ismember(cur_expt_set,[46 48 49 51])) = []; %problem with image background
cur_expt_set(cur_expt_set > 60) = []; %no rect

%% COMPUTE TRIAL DATA
all_stim_times = [];

all_t_axis = [];
all_rel_stimes = [];
all_rel_etimes = [];
all_binned_spks = [];
all_used_inds = [];
all_exptvec = [];
all_trialvec = [];
all_trial_start_times = [];
all_trial_end_times = [];
for ee = 1:length(cur_expt_set);
    fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
    cur_expt = cur_expt_set(ee);
    fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
    load(fname);
    
    n_trials = length(Expts{cur_expt}.Trials);
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    
    all_trial_start_times = [all_trial_start_times; trial_start_times(:)];
    all_trial_end_times = [all_trial_end_times; trial_end_times(:)];
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_expt}.Trials(tt).Start/1e4;
        
        cur_t_edges = [(Expts{cur_expt}.Trials(tt).Start(1)/1e4):1/Fsd:(Expts{cur_expt}.Trials(tt).End(end)/1e4)];
        cur_t_cents = cur_t_edges(1:end-1) + 2/Fsd;
        cur_binned_spks = nan(length(cur_t_cents),length(use_units));
        for cc = 1:length(use_units)
            cur_hist = histc(Clusters{use_units(cc)}.times,cur_t_edges);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
        end
        
        cur_used_inds = ones(length(cur_t_cents),1);
        cur_used_inds(1:beg_buffer) = 0; %get rid of first piece in each trial
%         cur_used_inds(end-end_buffer+1:end) = 0; %get rid of last piece in each trial
        
        all_stim_times = [all_stim_times; cur_stim_times];
        all_t_axis = [all_t_axis; cur_t_cents(:)];
        all_rel_stimes = [all_rel_stimes; cur_t_cents(:)- trial_start_times(tt)];
        all_rel_etimes = [all_rel_etimes; trial_end_times(tt) - cur_t_cents(:)];
        all_used_inds = [all_used_inds; cur_used_inds];
        all_binned_spks = [all_binned_spks; cur_binned_spks];
        all_exptvec = [all_exptvec; ones(length(cur_t_cents),1)*ee];
        all_trialvec = [all_trialvec; ones(length(cur_t_cents),1)*tt];
    end

end

%%
all_phasegrams = [];
all_ampgrams = [];
all_Vmat = [];
all_lfp_taxis = [];
all_lfp_expt = [];
for ee = 1:length(cur_expt_set);
    
    %%
    fprintf('Analyzing block %d of %d\n',ee,length(cur_expt_set));
    Vmat = [];
    phasegrams = [];
    ampgrams = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d of %d\n',ll,length(use_lfps));
        filename = sprintf('Expt%d.p%dFullV.mat',cur_expt_set(ee),use_lfps(ll));
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
        
        temp = cwt(dV,scales,'cmor1-1');
        phasegrams = cat(3,phasegrams,angle(temp)');
        ampgrams = cat(3,ampgrams,abs(temp)');
    end
    
    t_ax = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        t_ax = [t_ax downsample(cur_t_ax,dsf)];
    end
    t_ax(size(Vmat,1)+1:end) = [];
    
    cur_expt_inds = find(all_exptvec == ee);
    
    
    all_lfp_taxis = [all_lfp_taxis; t_ax(:)];
    
    all_lfp_expt = [all_lfp_expt; ones(length(t_ax),1)*ee];
    all_phasegrams = cat(1,all_phasegrams,phasegrams);
    all_ampgrams = cat(1,all_ampgrams,ampgrams);
    all_Vmat = cat(1,all_Vmat,Vmat);
    
%     unwr_phasegram = unwrap(phasegrams);
%     interp_phasegrams = interp1(t_ax,unwr_phasegram,all_t_axis(cur_expt_inds));
%     interp_phasegrams = mod(interp_phasegrams+pi,2*pi)-pi;
%     interp_ampgrams = interp1(t_ax,ampgrams,all_t_axis(cur_expt_inds));
%     
%     interp_Vmat = interp1(t_ax,Vmat,all_t_axis(cur_expt_inds));
%     
%     all_interp_Vmat = cat(1,all_interp_Vmat,interp_Vmat);
%     all_interp_phasegrams = cat(1,all_interp_phasegrams,interp_phasegrams);
%     all_interp_ampgrams = cat(1,all_interp_ampgrams,interp_ampgrams);
end

%%
corresp_expt_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_lfp_taxis));
corresp_used_inds = zeros(size(corresp_expt_inds));
uset = find(~isnan(corresp_expt_inds));
corresp_used_inds(uset) = all_used_inds(corresp_expt_inds(uset));
%%
full_ampgrams_norm = all_ampgrams;
full_ampgrams_norm(corresp_used_inds==0,:,:) = nan;
full_ampgrams_norm = bsxfun(@minus,full_ampgrams_norm,nanmean(full_ampgrams_norm));
full_ampgrams_norm = bsxfun(@rdivide,full_ampgrams_norm,nanstd(full_ampgrams_norm));

full_phasegrams_used = all_phasegrams;
full_phasegrams_used(corresp_used_inds==0,:,:) = nan;

full_lfps_norm = all_Vmat;
full_lfps_norm(corresp_used_inds==0,:) = nan;
full_lfps_norm = bsxfun(@minus,full_lfps_norm,nanmean(full_lfps_norm));
full_lfps_norm = bsxfun(@rdivide,full_lfps_norm,nanstd(full_lfps_norm));

%%

clear all_ampgrams all_phasegrams all_Vmat

trial_set = unique(all_trialvec);
n_trials = length(trial_set);


%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)
load ./jbeG081.em.mat
all_eye_vals = [];
all_eye_speed = [];
all_eye_ts = [];
for ee = 1:length(cur_expt_set);
    cur_set = find(all_exptvec==ee);
    [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data_new(Expt,all_t_axis(cur_set([1 end])));

    eye_dt = median(diff(eye_ts_interp));
    eye_fs = 1/eye_dt;
    lEyeXY = eye_vals_interp(:,1:2);
    rEyeXY = eye_vals_interp(:,3:4);
    clear sm_avg_eyepos eye_vel
    sm_avg_eyepos(:,1) = smooth(lEyeXY(:,1),3);
    sm_avg_eyepos(:,2) = smooth(lEyeXY(:,2),3);
    eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/eye_dt;
    eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/eye_dt;

    eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);

    all_eye_vals = [all_eye_vals; lEyeXY rEyeXY];
    all_eye_speed = [all_eye_speed; eye_speed];
    all_eye_ts = [all_eye_ts; eye_ts_interp'];
end
back_pts = 1 + find(diff(all_eye_ts) <= 0);
double_samples = [];
for i = 1:length(back_pts)
    next_forward = find(all_eye_ts > all_eye_ts(back_pts(i)-1),1,'first');
    double_samples = [double_samples back_pts(i):next_forward];
end
all_eye_ts(double_samples) = [];
all_eye_speed(double_samples) = [];
all_eye_vals(double_samples,:) = [];
eye_fs = 1/median(diff(all_eye_ts));

interp_eye_vals = interp1(all_eye_ts,all_eye_vals,all_lfp_taxis);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_lfp_taxis);

orth_eye_pos = interp_eye_vals(:,2);
par_eye_pos = interp_eye_vals(:,1);
bad_pts = find(abs(orth_eye_pos) > 1); %cheap way of detecting blinks because the right eye signal was unreliable even for that at times

all_eye_intrial = zeros(size(all_eye_ts));
for i = 1:length(all_trial_start_times)
    curset = find(all_eye_ts >= all_trial_start_times(i) & all_eye_ts <= all_trial_end_times(i));
    all_eye_intrial(curset) = 1;
end

%compute times saccades
sac_thresh = 10;
min_isi = 0.05;
orig_saccade_inds = 1 + find(all_eye_speed(1:end-1) < sac_thresh & all_eye_speed(2:end) > sac_thresh);
orig_saccade_inds = unique(orig_saccade_inds);
isis = [Inf; diff(orig_saccade_inds)]/eye_fs;
double_sacs = find(isis < min_isi);
orig_saccade_inds(double_sacs) = []; isis(double_sacs) = [];
not_intrial = find(all_eye_intrial(orig_saccade_inds) == 0);
orig_saccade_inds(not_intrial) = []; isis(not_intrial) = [];
interp_saccade_inds = round(interp1(all_lfp_taxis,1:length(all_lfp_taxis),all_eye_ts(orig_saccade_inds)));


peri_thresh = 3;
sac_starts = nan(size(interp_saccade_inds));
sac_stops = nan(size(interp_saccade_inds));
for i = 1:length(interp_saccade_inds)
    cur_up = find(interp_eye_speed(1:interp_saccade_inds(i)) < peri_thresh,1,'last');
    cur_down = find(interp_eye_speed(interp_saccade_inds(i):end) < peri_thresh,1,'first');
    if ~isempty(cur_up) & ~isempty(cur_down)
        sac_starts(i) = cur_up;
        sac_stops(i) = interp_saccade_inds(i) + cur_down;
    end
end
bad = find(isnan(sac_starts));
interp_saccade_inds(bad) = []; sac_starts(bad) = []; sac_stops(bad) = []; isis(bad) = [];


sac_starts = sac_starts - round(Fsd*0.01);
sac_stops = sac_stops + round(Fsd*0.01);
pre_pos = par_eye_pos(sac_starts);
post_pos = par_eye_pos(sac_stops);

sac_dist = abs(post_pos-pre_pos);

outgoing_sacs = [];
returning_sacs = [];
for ee = 1:length(cur_expt_set);
    cur_sac_set = find(all_lfp_expt(interp_saccade_inds) == ee);
    
    if expt_bar_ori(cur_expt_set(ee)) == 0
        cur_out = find(pre_pos(cur_sac_set) > -0.9 & post_pos(cur_sac_set) < -0.9);
        cur_ret = find(pre_pos(cur_sac_set) < -1 & post_pos(cur_sac_set) > -1);
    elseif expt_bar_ori(cur_expt_set(ee)) == 135
        cur_out = find(pre_pos(cur_sac_set) > -1 & post_pos(cur_sac_set) < -1);
        cur_ret = find(pre_pos(cur_sac_set) < -1.1 & post_pos(cur_sac_set) > -1.2);
    elseif expt_bar_ori(cur_expt_set(ee)) == 45
        cur_out = find(pre_pos(cur_sac_set) > -1 & post_pos(cur_sac_set) < -1);
        cur_ret = find(pre_pos(cur_sac_set) < -1.1 & post_pos(cur_sac_set) > -1);
    elseif expt_bar_ori(cur_expt_set(ee)) == 90
        cur_out = find(pre_pos(cur_sac_set) > -1 & post_pos(cur_sac_set) < -1);
        cur_ret = find(pre_pos(cur_sac_set) < -1 & post_pos(cur_sac_set) > -1);
    end
    outgoing_sacs = [outgoing_sacs; interp_saccade_inds(cur_sac_set(cur_out))];
    returning_sacs = [returning_sacs; interp_saccade_inds(cur_sac_set(cur_ret))];    
end

msacs = interp_saccade_inds(sac_dist' < 1); msacs(ismember(msacs,outgoing_sacs)) = []; msacs(ismember(msacs,returning_sacs)) = [];

uset = find(~isnan(corresp_expt_inds));
corresp_rel_stimes = nan(size(corresp_expt_inds));
corresp_rel_stimes(uset) = all_rel_stimes(corresp_expt_inds(uset));
corresp_rel_etimes = nan(size(corresp_expt_inds));
corresp_rel_etimes(uset) = all_rel_etimes(corresp_expt_inds(uset));

stim_on = 1+find(corresp_rel_stimes(1:end-1) < 0.7 & corresp_rel_stimes(2:end) > 0.7);
stim_off = 1+find(corresp_rel_stimes(1:end-1) < 1.4 & corresp_rel_stimes(2:end) > 1.4);


%%
use_out = find(corresp_rel_stimes(outgoing_sacs) > 0.8 & corresp_rel_stimes(outgoing_sacs) < 1.3);
use_ret = find(corresp_rel_stimes(returning_sacs) > 1.4 & corresp_rel_stimes(returning_sacs) < 1.9);

use_outgoing_sacs = outgoing_sacs(use_out);
use_returning_sacs = returning_sacs(use_ret);

msac_inds = find(ismember(interp_saccade_inds,msacs));
use_msacs_nsi = msacs(corresp_rel_stimes(msacs) >= 0.2 & isis(msac_inds) >= 0.2);
use_msacs = msacs(corresp_rel_stimes(msacs) >= 0.2);

%%
trial_lfp_inds = round(interp1(all_lfp_taxis,1:length(all_lfp_taxis),all_trial_start_times));

%%
backlag = round(Fsd*0.5);
forwardlag = round(Fsd*0.5);

[micro_trig_avgs,lags] = get_event_trig_avg(full_lfps_norm,use_msacs,backlag,forwardlag);
% [first_trig_avgs,lags] = get_event_trig_avg(full_lfps_norm,use_outgoing_sacs,backlag,forwardlag);
% [second_trig_avgs,lags] = get_event_trig_avg(full_lfps_norm,use_returning_sacs,backlag,forwardlag);

[on_trig_avgs,lags] = get_event_trig_avg(full_lfps_norm,stim_on,backlag,forwardlag);
[off_trig_avgs,lags] = get_event_trig_avg(full_lfps_norm,stim_off,backlag,forwardlag);

[micro_nsi_trig_avgs,lags] = get_event_trig_avg(full_lfps_norm,use_msacs_nsi,backlag,forwardlag);

%%
% [first_trig_specgram,lags] = get_event_trig_avg(full_ampgrams_norm,use_outgoing_sacs,backlag,forwardlag);
% [second_trig_specgram,lags] = get_event_trig_avg(full_ampgrams_norm,use_returning_sacs,backlag,forwardlag);
[on_trig_specgram,lags] = get_event_trig_avg(full_ampgrams_norm,stim_on,backlag,forwardlag);
[off_trig_specgram,lags] = get_event_trig_avg(full_ampgrams_norm,stim_off,backlag,forwardlag);
[micro_trig_specgram,lags] = get_event_trig_avg(full_ampgrams_norm,use_msacs,backlag,forwardlag);


% [second_trig_phaselock,lags] = get_event_trig_avg(exp(full_phasegrams_used*(1i)),use_returning_sacs,backlag,forwardlag);
% second_trig_phaselock = abs(second_trig_phaselock);
% [first_trig_phaselock,lags] = get_event_trig_avg(exp(full_phasegrams_used*(1i)),use_outgoing_sacs,backlag,forwardlag);
% first_trig_phaselock = abs(first_trig_phaselock);
[on_trig_phaselock,lags] = get_event_trig_avg(exp(full_phasegrams_used*(1i)),stim_on,backlag,forwardlag);
on_trig_phaselock = abs(on_trig_phaselock);
[off_trig_phaselock,lags] = get_event_trig_avg(exp(full_phasegrams_used*(1i)),stim_off,backlag,forwardlag);
off_trig_phaselock = abs(off_trig_phaselock);
[micro_trig_phaselock,lags] = get_event_trig_avg(exp(full_phasegrams_used*(1i)),use_msacs,backlag,forwardlag);
micro_trig_phaselock = abs(micro_trig_phaselock);

[micro_nsi_trig_specgram,lags] = get_event_trig_avg(full_ampgrams_norm,use_msacs_nsi,backlag,forwardlag);
[micro_nsi_trig_phaselock,lags] = get_event_trig_avg(exp(full_phasegrams_used*(1i)),use_msacs_nsi,backlag,forwardlag);
micro_nsi_trig_phaselock = abs(micro_nsi_trig_phaselock);

%%


%%
f_backlag = round(0);
f_forwardlag = round(Fsd*2);
[trial_trig_specgram,f_lags] = get_event_trig_avg(full_ampgrams_norm,trial_lfp_inds,f_backlag,f_forwardlag);
% 
% for i = 1:length(use_lfps)
%     fprintf('LFP %d of %d\n',i,length(use_lfps));
%     
%     pcolor(f_lags/Fsd,wfreqs,squeeze(trial_trig_specgram(:,:,i))');shading flat
%     set(gca,'yscale','log');  colorbar
% 
%     pause
%     clf
% end
%%
cd ~/Data/bruce/G081
save G081_lfp_analysis_imback use_lfps wfreqs lags Fsd f_lags *trig_avgs *trig_phaselock *trig_specgram use_* *lag cur_expt_set

%%
yscale = 'log';
spec_cax = [-0.2 0.2];
phase_cax = [0 0.3];
close all
for i = 1:length(use_lfps)
    fprintf('LFP %d of %d\n',i,length(use_lfps));
    
    subplot(3,3,1)
%     plot(lags/Fsd,first_trig_avgs(:,i))
    plot(lags/Fsd,on_trig_avgs(:,i))
    ylim([-0.4 0.4])
    yl = ylim();
   line([0 0],yl,'color','k')
   
    subplot(3,3,4)
%     pcolor(lags/Fsd,wfreqs,squeeze(first_trig_specgram(:,:,i))');shading flat
    pcolor(lags/Fsd,wfreqs,squeeze(on_trig_specgram(:,:,i))');shading flat
    set(gca,'yscale',yscale); %colorbar
    yl = ylim();
   line([0 0],yl,'color','k')
    caxis(spec_cax)
   
   subplot(3,3,7)
%     pcolor(lags/Fsd,wfreqs,squeeze(first_trig_phaselock(:,:,i))');shading flat
    pcolor(lags/Fsd,wfreqs,squeeze(on_trig_phaselock(:,:,i))');shading flat
    set(gca,'yscale',yscale); % colorbar
    caxis(phase_cax)
    yl = ylim();
   line([0 0],yl,'color','k')
  
    subplot(3,3,2)
%     plot(lags/Fsd,second_trig_avgs(:,i))
    plot(lags/Fsd,off_trig_avgs(:,i))
    ylim([-0.4 0.4])
    xlim([-0.5 0.5])
   yl = ylim();
   line([0 0],yl,'color','k')
    
    subplot(3,3,5)
%     pcolor(lags/Fsd,wfreqs,squeeze(second_trig_specgram(:,:,i))');shading flat
    pcolor(lags/Fsd,wfreqs,squeeze(off_trig_specgram(:,:,i))');shading flat
    set(gca,'yscale',yscale); %colorbar
    xlim([-0.5 0.5])
   yl = ylim();
   line([0 0],yl,'color','k')
    caxis(spec_cax)

    subplot(3,3,8)
%     pcolor(lags/Fsd,wfreqs,squeeze(second_trig_phaselock(:,:,i))');shading flat
    pcolor(lags/Fsd,wfreqs,squeeze(off_trig_phaselock(:,:,i))');shading flat
    set(gca,'yscale',yscale); %colorbar
    xlim([-0.5 0.5])
   caxis(phase_cax)
   yl = ylim();
   line([0 0],yl,'color','k')

       subplot(3,3,3)
    plot(lags/Fsd,micro_trig_avgs(:,i))
%     plot(lags/Fsd,micro_nsi_trig_avgs(:,i))
    ylim([-0.3 0.3])
   yl = ylim();
   line([0 0],yl,'color','k')
    
    subplot(3,3,6)
    pcolor(lags/Fsd,wfreqs,squeeze(micro_trig_specgram(:,:,i))');shading flat
%     pcolor(lags/Fsd,wfreqs,squeeze(micro_nsi_trig_specgram(:,:,i))');shading flat
    set(gca,'yscale',yscale); %colorbar
   yl = ylim();
   line([0 0],yl,'color','k')
    caxis([-0.15 0.15])

    subplot(3,3,9)
    pcolor(lags/Fsd,wfreqs,squeeze(micro_trig_phaselock(:,:,i))');shading flat
%     pcolor(lags/Fsd,wfreqs,squeeze(micro_nsi_trig_phaselock(:,:,i))');shading flat
    set(gca,'yscale',yscale); %colorbar
   caxis([0 0.2])
   yl = ylim();
   line([0 0],yl,'color','k')

   pause
    clf
end

%%
yscale = 'linear';

spec_cax = [-0.2 0.2];
phase_cax = [0 0.2];
mspec_cax = [-0.1 0.1];
mphase_cax = [0 0.1];

figure
subplot(3,1,1)
pcolor(lags/Fsd,wfreqs,squeeze(mean(first_trig_specgram,3))'); shading flat
    caxis(spec_cax)
    set(gca,'yscale',yscale);
subplot(3,1,2)
pcolor(lags/Fsd,wfreqs,squeeze(mean(second_trig_specgram,3))'); shading flat
    caxis(spec_cax)
    set(gca,'yscale',yscale);

subplot(3,1,3)
pcolor(lags/Fsd,wfreqs,squeeze(mean(micro_trig_specgram,3))'); shading flat
    caxis(mspec_cax)
    set(gca,'yscale',yscale);
    
    figure
subplot(3,1,1)
pcolor(lags/Fsd,wfreqs,squeeze(mean(first_trig_phaselock,3))'); shading flat
    caxis(phase_cax)
    set(gca,'yscale',yscale);
subplot(3,1,2)
pcolor(lags/Fsd,wfreqs,squeeze(mean(second_trig_phaselock,3))'); shading flat
    caxis(phase_cax)
    set(gca,'yscale',yscale);

subplot(3,1,3)
pcolor(lags/Fsd,wfreqs,squeeze(mean(micro_trig_phaselock,3))'); shading flat
    caxis(mphase_cax)
        set(gca,'yscale',yscale);
