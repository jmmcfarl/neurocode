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
dsf = 60;Fsd = Fs/dsf;
[b_alpha,a_alpha] = butter(2,[7 15]/(Fsd/2));
use_lfps = [1:4:96];
scales = logspace(log10(10),log10(75),40);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

load ./Expt1_newcompiled_data_fixeddelay_d1p25_34.mat full_t full_expt_vec

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
Expt_nu = [1 2 17:19 23 24];
% Expt_nu = [13];
% Expt_nu = [7]; %expt3
n_allunits = 96;

all_interp_alpha_phase = [];
all_interp_alpha = [];
for ee = 1:length(Expt_nu)
    
    all_alpha_phase = [];
    all_alpha = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d\n',ll);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = decimate(V(:),dsf);
        V_alpha = filtfilt(b_alpha,a_alpha,V);
        V_alpha_phase = angle(hilbert(V_alpha));

        all_alpha_phase(:,ll) = V_alpha_phase;
        all_alpha(:,ll) = V_alpha;
    end
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
        
    all_alpha = zscore(all_alpha);
    
    %buffer for spike response delay
    t_ax = t_ax - 0.05;
    
    cur_inds = find(full_expt_vec==Expt_nu(ee));
    interp_alpha_phase = interp1(t_ax,all_alpha_phase,full_t(cur_inds));
     interp_alpha = interp1(t_ax,all_alpha,full_t(cur_inds));
   
    all_interp_alpha_phase = [all_interp_alpha_phase; interp_alpha_phase];
    all_interp_alpha = [all_interp_alpha; interp_alpha];
    
end

%%
cd ~/Data/bruce/7_15_12/G034/

save expt1_lfp_alpha_phase_7_15_delay all_interp_alpha* nearest_lfps

%%
