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

% Fs = 3e4;
% dsf = 60;Fsd = Fs/dsf;
% [b_alpha,a_alpha] = butter(2,[7 14]/(Fsd/2));
% use_lfps = [1:2:96];

% scales = logspace(log10(5),log10(75),40);
% scales = [scales 90 125 150 200 250];
% wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
% nwfreqs = length(wfreqs);

Fs = 3e4;
dsf = 120;Fsd = Fs/dsf;
use_lfps = [1:96];

scales = logspace(log10(8),log10(75),30);
scales = [scales 85 100 115 130 150];
scales = scales*60/dsf;
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

load ./Expt1_newcompiled_data_fixeddelay_d1p25_34_ht.mat full_t full_expt_vec full_trial_vec full_im_ids full_image_vec

%%
% frac = 4;
frac = 4;
cur_dt = median(diff(full_t));
desired_dt = cur_dt/frac;
trial_start_inds = 1+find(diff(full_trial_vec) ~= 0);
trial_stop_inds = trial_start_inds - 1;
trial_start_inds = [1; trial_start_inds];
trial_stop_inds = [trial_stop_inds; length(full_trial_vec)];
new_t_axis = [];
new_expt_vec = [];
new_trial_vec = [];
new_trial_originds = [];
cur_trial = 1;
for i = 1:length(trial_stop_inds)
    temp = full_t(trial_start_inds(i)):desired_dt:(full_t(trial_stop_inds(i))+cur_dt-desired_dt);
   new_t_axis = [new_t_axis temp];
   new_expt_vec = [new_expt_vec repmat(full_expt_vec(trial_start_inds(i)),1,length(temp))];
   new_trial_vec = [new_trial_vec cur_trial*ones(1,length(temp))];
   temp2 = floor(trial_start_inds(i):1/frac:trial_stop_inds(i));
   temp2 = [temp2 ones(1,length(temp)-length(temp2))*trial_stop_inds(i)];
   cur_trial = cur_trial + 1;
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
Expt_nu = [1 2 17 19 23 24];
% Expt_nu = [13];
% Expt_nu = [7]; %expt3
n_allunits = 96;

% all_interp_alpha_phase = [];
% all_interp_alpha = [];
all_interp_phases = [];
all_interp_amps = [];
for ee = 1:length(Expt_nu)
    
%     all_alpha_phase = [];
%     all_alpha = [];
all_phasegram = [];
all_ampgram = [];
for ll = 1:length(use_lfps)
        fprintf('Electrode %d\n',ll);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = decimate(V(:),dsf);
%         V_alpha = filtfilt(b_alpha,a_alpha,V);
        %         V_alpha_phase = angle(hilbert(V_alpha));
        %
        %         all_alpha_phase(:,ll) = V_alpha_phase;
        %         all_alpha(:,ll) = V_alpha;
        
        temp = cwt(V,scales,'cmor1-1');
        all_phasegram(:,:,ll) = angle(temp)';
        all_ampgram(:,:,ll) = abs(temp)';
        
    end
    
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
        
%     all_alpha = zscore(all_alpha);
    
    %buffer for spike response delay
%     t_ax = t_ax - 0.05;
    
    cur_inds = find(new_expt_vec==Expt_nu(ee));
%     interp_phases = interp1(t_ax,all_phasegram,new_t_axis(cur_inds));
    interp_amps = interp1(t_ax,all_ampgram,new_t_axis(cur_inds));
%     interp_alpha_phase = interp1(t_ax,all_alpha_phase,new_t_axis(cur_inds));
%      interp_alpha = interp1(t_ax,all_alpha,new_t_axis(cur_inds));
     unwr_phasegram = unwrap(all_phasegram);
    interp_phases = interp1(t_ax,unwr_phasegram,new_t_axis(cur_inds));
    interp_phases = mod(interp_phases+pi,2*pi)-pi;
  
    all_interp_phases = cat(1,all_interp_phases,interp_phases);
    all_interp_amps = cat(1,all_interp_amps,interp_amps);
%     all_interp_alpha_phase = [all_interp_alpha_phase; interp_alpha_phase];
%     all_interp_alpha = [all_interp_alpha; interp_alpha];
    
end

%%
cd ~/Data/bruce/7_15_12/G034/

% save expt1_lfp_alpha_phase_7_14_delay_v2 all_interp_alpha* nearest_lfps new_* desired_dt
save('expt1_lfp_alpha_phase_delay_all_sameelec_new.mat','-v7.3','all_interp_*','nearest_lfps','new_*','desired_dt','wfreqs','frac')

%%
