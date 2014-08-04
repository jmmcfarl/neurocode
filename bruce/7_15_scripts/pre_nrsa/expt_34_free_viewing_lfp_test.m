clear all
% close all
addpath('~/Data/bruce/7_15_12/')

cd ~/Data/bruce/7_15_12

cd G034/
load ./CellList.mat
load ./G034Expts.mat

samp_fac = 0.5;
dt = samp_fac*118/1e4;
fst = 1/dt;
frames_per_jump = 44/samp_fac;

Pix2Deg = 0.018837;
dsfrac = 1.25;
Fsd = 1/Pix2Deg/dsfrac;
Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);

rf_cent = [0.34 -0.43];

min_trial_dur = 1;
max_sac_amp = 0.5;
single_units = find(CellList(1,:,1) > 0);
n_sus = length(single_units);

%%
% Expt_nu = [3 8 15 19 26 29];
Expt_nu = [13 14 15 16 25 28 29];
% Expt_nu = [13];
% Expt_nu = [7]; %expt3
n_allunits = 96;

full_expt_vec = [];
full_trial_vec = [];
full_binned_spks = [];
full_trial_durs = [];
full_t = [];
full_im_flips = [];
for ee = 1:length(Expt_nu)
    
    fprintf('Analyzing Expt %d of %d\n',ee,length(Expt_nu));
    n_mus = 96;
    load(sprintf('Expt%dClusterTimes.mat',Expt_nu(ee)));
    
    Trial_starts = [Expts{Expt_nu(ee)}.Trials(:).Start]/1e4;
    Trial_ends = [Expts{Expt_nu(ee)}.Trials(:).End]/1e4;
            
    endEvents = [Expts{Expt_nu(ee)}.Trials(:).endevent];
    Trial_durs = (Trial_ends-Trial_starts);
    used_trials = find(Trial_durs > min_trial_dur);
    
    tt = nan(length(used_trials),1);
    for i = 1:length(used_trials)
        cur_n_tbins = floor(Trial_durs(used_trials(i))/dt);
        cur_n_shifts = ceil(cur_n_tbins/frames_per_jump);
        tt(i) = cur_n_shifts;
        
        cur_t_edges = Trial_starts(used_trials(i)):dt:(Trial_starts(used_trials(i)) + dt*cur_n_tbins);
        cur_t_cents = 0.5*cur_t_edges(1:end-1) + 0.5*cur_t_edges(2:end);
        
        cur_im_flips = Trial_starts(used_trials(i)):dt*frames_per_jump:Trial_ends(used_trials(i));
        
        cur_binned_spks = nan(n_mus,length(cur_t_cents));
        for j = 1:n_mus
            temp = histc(Clusters{j}.times,cur_t_edges);
            cur_binned_spks(j,:) = temp(1:end-1);
        end
        
        full_expt_vec = [full_expt_vec; ones(length(cur_t_cents),1)*Expt_nu(ee)];
        full_trial_vec = [full_trial_vec; ones(length(cur_t_cents),1)*used_trials(i)];
        full_binned_spks = [full_binned_spks; cur_binned_spks'];
        full_trial_durs = [full_trial_durs; ones(length(cur_t_cents),1)*Trial_durs(used_trials(i))];
        full_t = [full_t cur_t_cents];
        full_im_flips = [full_im_flips cur_im_flips];
    end
end

%%
sm_fac = 1;
clear smooth_binned_spks
for t = 1:96
smooth_binned_spks(:,t) = smooth(full_binned_spks(:,t),sm_fac);
end
smooth_binned_spks = zscore(smooth_binned_spks);
full_binned_spks = zscore(full_binned_spks);
load ./expt2_eyeanalysis_g034.mat
% cd ~/Data/bruce/7_15_12/G029
% load ./all_eyedata_expt3_34
% cd ~/Data/bruce/7_15_12/G034
net_spk = mean(smooth_binned_spks,2);

%%
load Expt15.p1FullV.mat
% load Expt7.p70FullV.mat
Fs = 3e4;
V = double(FullV.V);
t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
dsf = 60;Fsd = Fs/dsf;
V = decimate(V,dsf);
t_ax = downsample(t_ax,dsf);
[b,a] = butter(2,[4 14]/(Fsd/2));
V_f2 = zscore(filtfilt(b,a,V));
%%
beg_pt = find(full_expt_vec == 15,1,'first');
beg_t = full_t(beg_pt);
figure
plot(all_t-beg_t,all_eyespeed/10,'b')
hold on
plot(full_t-beg_t,full_binned_spks(:,4),'g')
% plot(full_t,smooth_binned_spks(:,23),'b')
% xlim([7754 7755.5])
plot(t_ax-beg_t,zscore(V_f2)-3,'k')
plot(t_ax-beg_t,zscore(V)-3,'r')
% plot(full_t,smooth_binned_spks(:,4),'k')
plot(full_t-beg_t,net_spk*5,'k')


%%
params.Fs = Fsd;
params.tapers = [5 9];
win = 20;
[S,f] = mtspectrumsegc(V,win,params,1);
