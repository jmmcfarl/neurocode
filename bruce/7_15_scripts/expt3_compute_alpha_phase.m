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

Fs = 3e4;
dsf = 120;Fsd = Fs/dsf;
use_lfps = [1:96];

scales = logspace(log10(7),log10(75),25);
scales = [scales 85 100 115 130];
scales = scales*60/dsf;
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

load ./Expt3_fixbased_data.mat

%%
cur_dt = median(diff(full_t));
desired_dt = cur_dt/2;
trial_start_inds = 1+find(diff(full_trial_vec) ~= 0);
trial_stop_inds = trial_start_inds - 1;
trial_start_inds = [1; trial_start_inds];
trial_stop_inds = [trial_stop_inds; length(full_trial_vec)];
new_t_axis = [];
new_expt_vec = [];
new_trial_vec = [];
for i = 1:length(trial_stop_inds)
    temp = full_t(trial_start_inds(i)):desired_dt:full_t(trial_stop_inds(i));
    new_t_axis = [new_t_axis temp];
    new_expt_vec = [new_expt_vec repmat(full_expt_vec(trial_start_inds(i)),1,length(temp))];
    toffset = (full_expt_vec(trial_start_inds(i))-7)*100;
    new_trial_vec = [new_trial_vec repmat(toffset+full_trial_vec(trial_start_inds(i)),1,length(temp))];
end

%%
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;
nearest_lfps = nan(length(use_lfps),1);
for ll = 1:96
   all_dists = sqrt((X_pos-X_pos(ll)).^2 + (Y_pos-Y_pos(ll)).^2); 
%    all_dists(ll) = inf;
   [~,best_loc] = min(all_dists(use_lfps));
   nearest_lfps(ll) = best_loc;
end
%%
% Expt_nu = [3 8 15 19 26 29];
% Expt_nu = [13 14 15 16 25 28 29];
Expt_nu = [7:12]; %expt 3 34
% Expt_nu = [13];
% Expt_nu = [7]; %expt3
n_allunits = 96;

all_interp_phases = [];
all_interp_amps = [];
for ee = 1:length(Expt_nu)
    
all_phasegram = [];
all_ampgram = [];
for ll = 1:length(use_lfps)
        fprintf('Electrode %d\n',ll);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = decimate(V(:),dsf);
        temp = cwt(V,scales,'cmor1-1');
        all_phasegram(:,:,ll) = angle(temp)';
        all_ampgram(:,:,ll) = abs(temp)';
    end
    
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
    
    cur_inds = find(new_expt_vec==Expt_nu(ee));
    interp_phases = interp1(t_ax,all_phasegram,new_t_axis(cur_inds));
    interp_amps = interp1(t_ax,all_ampgram,new_t_axis(cur_inds));
   
    all_interp_phases = cat(1,all_interp_phases,interp_phases);
    all_interp_amps = cat(1,all_interp_amps,interp_amps);
    
end

%%
cd ~/Data/bruce/7_15_12/G034/

save('expt3_lfp_alpha_phase_delay_all.mat','-v7.3','all_interp_phases','all_interp_amps','nearest_lfps','new_*','desired_dt','wfreqs')

%%
