clear all
close all

cd ~/Data/bruce/2_27_12

load Blocks
load lemM232A.51.lfp.mat

blockid = 1;
%%
%get start times of each LFP trial
n_lfp_trials = length(LFP.Trials);
lfp_trial_start = nan(n_lfp_trials,1);
for i = 1:n_lfp_trials
    lfp_trial_start(i) = LFP.Trials(i).Start/1e4; %originally in tenths of ms
    lfp_trial_stop(i) = LFP.Trials(i).End/1e4;
    lfp_dur(i) = size(LFP.Trials(i).LFP,1)/1000;
end
% lfp_trial_stop = lfp_trial_start+lfp_dur';

%get block start times
block_trial_times = Blocks{blockid}.blocktimes;

%% compute LFP signal
ov_t = block_trial_times(1,1):.001:block_trial_times(2,end);

lfp_time = [];
lfp_samps = [];
for i = 1:n_lfp_trials
    lfp_time = [lfp_time linspace(lfp_trial_start(i),lfp_trial_start(i)+lfp_dur(i),size(LFP.Trials(i).LFP,1))];
    lfp_samps = [lfp_samps; LFP.Trials(i).LFP];
end
% lfp_samps = zscore(lfp_samps);

Fs = 1000;
niqf = Fs/2;
dsf = 4;
Fsd = Fs/dsf;
s
[b,a] = butter(2,[0.5 100]/niqf);
lfp_samps = filtfilt(b,a,lfp_samps);
lfp_samps_d = downsample(lfp_samps,dsf);
lfp_time_d = downsample(lfp_time,dsf);

%% Compute LFP pasegram
cur_lfp_id = 20;

scales = logspace(log10(3),log10(150),50);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

coefs = cwt(lfp_samps_d(:,cur_lfp_id),scales,'cmor1-1');

phasegram = angle(coefs);
% scalgram = abs(coefs);
% ampgram = real(coefs);
%
% scalgram = bsxfun(@minus,scalgram,mean(scalgram(:,~art_inds),2));
% scalgram = bsxfun(@rdivide,scalgram,std(scalgram(:,~art_inds),[],2));
%
% ampgram = bsxfun(@minus,ampgram,mean(ampgram(:,~art_inds),2));
% ampgram = bsxfun(@rdivide,ampgram,std(ampgram(:,~art_inds),[],2));

%% compute binned spike times and overall phase-locking stats
cmean = zeros(10,length(wfreqs));
ckapp = zeros(10,length(wfreqs));
cpval = zeros(10,length(wfreqs));
for cellid = 1:10;
    cellid
    cur_spktimes = Blocks{blockid}.spktimes{cellid};
    spikes_binned = hist(cur_spktimes,lfp_time_d);
    spikes_binned(end) = 0;
    bad_inds = 1+find(diff(lfp_time_d) > 1/Fsd+.001);
    spikes_binned(bad_inds) = 0;
    un_spk_cnts = unique(spikes_binned);
    spk_bins{cellid} = [];
    for i = 1:length(un_spk_cnts)
        cur_set = find(spikes_binned == un_spk_cnts(i));
        spk_bins{cellid} = [spk_bins{cellid}; repmat(cur_set(:),un_spk_cnts(i),1)];
    end
    spk_bins{cellid} = sort(spk_bins{cellid});
    
    for wbin = 1:length(wfreqs)
        cmean(cellid,wbin) = circ_mean(phasegram(wbin,spk_bins{cellid})');
        ckapp(cellid,wbin) = circ_kappa(phasegram(wbin,spk_bins{cellid})');
        cpval(cellid,wbin) = circ_otest(phasegram(wbin,spk_bins{cellid})');
    end
    
end
%%
% cellid = 6;
% close all
% wbin = 3;
% x = linspace(-pi,pi,30);
% [tout,rout] = rose(phasegram(wbin,:),x);
% rout = rout/sum(rout);
% [tout,routspk] = rose(phasegram(wbin,spk_bins{cellid}),x);
% routspk = routspk/sum(routspk);
% % polar(tout,rout)
% % hold on
% polar(tout,routspk,'r')
%
% figure
% polar(tout,routspk-rout,'k')
%
%% FIND FIXATION TIMES
accept_window = [-8 8;-8 8];
sac_eyespeed = 10;
thresh_eyespeed = 2.5;
sac_buffer = 0;
min_fix_dur = 0.1;

cd ~/Data/bruce/2_27_12/saccades/
load(sprintf('lemM232.5%d.em.sac.mat',blockid))

% identify saccade start and stop times
EyeStartT = Expt.Trials.Start/10000; % time of first eye sample
EyeEndT = Expt.Trials.End/10000; % time of last eye sample
Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
eyets = EyeStartT:Eyedt:EyeEndT; %eye tracking time axis (sec)
sac_buffer_inds = round(sac_buffer/Eyedt);

reye_pos = [Expt.Trials.Eyevals.rh Expt.Trials.Eyevals.rv];
leye_pos = [Expt.Trials.Eyevals.lh Expt.Trials.Eyevals.lv];

fprintf('Computing saccade times\n');
avg_eyepos = (reye_pos + leye_pos)/2;
clear sm_avg_eyepos eye_vel
sm_avg_eyepos(:,1) = smooth(avg_eyepos(:,1),3);
sm_avg_eyepos(:,2) = smooth(avg_eyepos(:,2),3);
eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/Eyedt;
eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/Eyedt;
eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);

%find saccade start and stop indices
sac_inds = find(eye_speed(1:end-1) < sac_eyespeed & eye_speed(2:end) > sac_eyespeed);

sac_start_inds = nan(size(sac_inds));
sac_stop_inds = nan(size(sac_inds));
for i = 1:length(sac_inds)
    temp = find(eye_speed(1:sac_inds(i)) < thresh_eyespeed,1,'last');
    if ~isempty(temp)
        sac_start_inds(i) = temp;
    end
    temp = find(eye_speed(sac_inds(i):end) < thresh_eyespeed,1,'first');
    if ~isempty(temp)
        sac_stop_inds(i) = sac_inds(i)+temp-1;
    end
end

%identify start and stop times of unique saccades
sac_vec = zeros(size(reye_pos,1),1);
for i = 1:length(sac_start_inds)
    if ~isnan(sac_start_inds(i)) & ~isnan(sac_stop_inds(i))
        sac_vec(sac_start_inds(i):sac_stop_inds(i)+sac_buffer_inds) = 1;
    end
end
sac_vec(length(eyets)+1:end) = [];
sac_vec([1 end]) = 0;
sac_start_indsn = find(sac_vec(1:end-1) == 0 & sac_vec(2:end) == 1);
sac_stop_indsn = find(sac_vec(1:end-1) == 1 & sac_vec(2:end) == 0);
if length(sac_start_indsn) ~= length(sac_stop_indsn)
    error('saccade mis-alignment');
end

sac_start_times = eyets(sac_start_indsn);
sac_stop_times = eyets(sac_stop_indsn);

%compute saccade amplitudes, velocities, and durations
sac_dx = avg_eyepos(sac_stop_indsn,1) - avg_eyepos(sac_start_indsn,1);
sac_dy = avg_eyepos(sac_stop_indsn,2) - avg_eyepos(sac_start_indsn,2);
sac_amps = sqrt(sac_dx.^2+sac_dy.^2);
sac_peakvel = zeros(size(sac_start_indsn));
for i = 1:length(sac_start_indsn)
    sac_peakvel(i) = max(eye_speed(sac_start_indsn(i):sac_stop_indsn(i)));
end
sac_durs = sac_stop_times-sac_start_times;

% interpolate eye signal onto stimulus time axis
%     eye_interp = interp1(eyets(1:end-1),reye_pos,recon_t);
eye_interp = interp1(eyets(1:end-1),leye_pos,lfp_time_d);

%create an interpolated 0-1 vector of saccade times
sac_vec_interp = zeros(size(lfp_time_d));
for i = 1:length(sac_start_times)
    cur_set = find(lfp_time_d >= sac_start_times(i) & lfp_time_d <= sac_stop_times(i));
    sac_vec_interp(cur_set) = 1;
end
sac_vec_interp([1 end]) = 0;

%indices of the start and stops of fixations
fix_start_inds = 1+[0 find(sac_vec_interp(1:end-1) == 1 & sac_vec_interp(2:end) == 0)];
fix_stop_inds = 1+[find(sac_vec_interp(1:end-1) == 0 & sac_vec_interp(2:end) == 1) length(lfp_time_d)-1];
if length(fix_start_inds) ~= length(fix_stop_inds)
    error('Fixation mis-alignment');
end
fix_durs = (fix_stop_inds - fix_start_inds)*1/Fsd;
too_short = find(fix_durs < min_fix_dur); %only keep fixations that are minimum duration

fprintf('%d of %d fixations too short\n',length(too_short),length(fix_durs));
fix_start_inds(too_short) = []; fix_stop_inds(too_short) = []; fix_durs(too_short) = [];
fix_start_times = lfp_time_d(fix_start_inds);

%determine correspondence between fixation starts and previous saccades
fix_prev_sac_inds = zeros(size(fix_start_inds));
for i = 1:length(fix_start_inds)
    [~,fix_prev_sac_inds(i)] = min(abs(sac_stop_times-fix_start_times(i)));
end
fix_stop_times = lfp_time_d(fix_stop_inds);
%%
dt = 0.010;
backwin = .1;

cellid = 6;

recon_dt = lfp_time_d(1):dt:lfp_time_d(end);
time_bin_edges = [(recon_dt(1)-dt/2) (recon_dt+dt/2)];

phasegram_d = interp1(lfp_time_d,phasegram',recon_dt);
lfp_samps_dd = interp1(lfp_time_d,lfp_samps_d,recon_dt);

spike_times = Blocks{blockid}.spktimes{cellid};

time_since_fix = [];
fix_inds = [];
spikebins = [];
t_inds = [];

spike_bin_vecs = zeros(10,length(recon_dt));
for cellid = 1:10
    spikes_binned = histc(Blocks{blockid}.spktimes{cellid},time_bin_edges);
    spikes_binned(end) = [];
    spike_bin_vecs(cellid,:) = spikes_binned;
end
for i = 1:length(fix_start_times)
    cur_fix_inds = find(recon_dt > fix_start_times(i)-backwin & recon_dt <= fix_stop_times(i));
    t_inds = [t_inds; cur_fix_inds'];
    time_since_fix = [time_since_fix; (1:length(cur_fix_inds))'];
    fix_inds = [fix_inds; ones(length(cur_fix_inds),1)*i];
    spikebins = [spikebins; spike_bin_vecs(:,cur_fix_inds)'];
end

%%
for_win = 0.5;
cur_taxis = (round(-backwin/dt):round(for_win/dt))*dt;
nspikes = zeros(10,length(cur_taxis),1);
nsamps = zeros(length(cur_taxis),1);
wkappa = zeros(length(cur_taxis),length(wfreqs));
wmean = zeros(length(cur_taxis),length(wfreqs));
spk_kappa = zeros(10,length(cur_taxis),length(wfreqs));
spk_mean = zeros(10,length(cur_taxis),length(wfreqs));
lfp_trg_avg = zeros(24,length(cur_taxis));
for i = 1:length(cur_taxis)
    i
    cur_set = find(time_since_fix == i);
    nsamps(i) = length(cur_set);
    lfp_trg_avg(:,i) = mean(lfp_samps_dd(t_inds(cur_set),:));
    for w = 1:length(wfreqs)
        wkappa(i,w) = circ_kappa(phasegram_d(t_inds(cur_set),w));
        wmean(i,w) = circ_mean(phasegram_d(t_inds(cur_set),w));
    end
    for cellid = 1:10
        nspikes(cellid,i) = sum(spikebins(cur_set,cellid));
        for w = 1:length(wfreqs)
            spk_kappa(cellid,i,w) = circ_kappa(phasegram_d(t_inds(cur_set),w),spikebins(cur_set,cellid));
            spk_mean(cellid,i,w) = circ_mean(phasegram_d(t_inds(cur_set),w),spikebins(cur_set,cellid));
        end
    end
end
