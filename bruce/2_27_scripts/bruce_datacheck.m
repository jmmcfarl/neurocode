clear all
close all
clc

% cd /Users/James/Data/bruce/2_27_12
% load ./Blocks.mat
% 
% blockid = 5;
% 
% load(sprintf('lemM232A.5%d.lfp.mat',blockid));
% 
% n_lfp_trials = length(LFP.Trials);
% lfp_start_t = zeros(n_lfp_trials,1);
% lfp_end_t = zeros(n_lfp_trials,1);
% lfp_n_samps = zeros(n_lfp_trials,1);
% for i = 1:n_lfp_trials
%     lfp_start_t(i) = LFP.Trials(i).Start/1e4;
%     lfp_end_t(i) = LFP.Trials(i).End/1e4;
%     lfp_n_samps(i) = size(LFP.Trials(i).LFP,1);
% end
% lfp_durs = lfp_end_t - lfp_start_t;
% cd ~/Data/bruce/2_27_12/saccades/
% load(sprintf('lemM232.5%d.em.sac.mat',blockid));
% 
% n_eye_trials = length(Expt.Trials);
% eye_start_t = zeros(n_eye_trials,1);
% eye_end_t = zeros(n_eye_trials,1);
% for i = 1:n_eye_trials
% eye_start_t(i) = Expt.Trials(i).Start/1e4; % time of first eye sample
% eye_end_t(i) = Expt.Trials(i).End/1e4; % time of last eye sample
% end
%     
% disp('Block times')
% Blocks{blockid}.blocktimes'
% 
% disp('Stim times')
% Blocks{blockid}.stimtime([1 end])'
% 
% disp('LFP times')
% [lfp_start_t lfp_end_t lfp_durs lfp_n_samps/1e3]
% 
% disp('Eye times')
% [eye_start_t eye_end_t]

cd /Users/James/Data/bruce/ExptB/Raw
load ./Blocks.mat

blockid = 7;
expno = Blocks{blockid}.exptno;

load(sprintf('lemM237A.%d.lfp.mat',expno));

n_lfp_trials = length(LFP.Trials);
lfp_start_t = zeros(n_lfp_trials,1);
lfp_end_t = zeros(n_lfp_trials,1);
lfp_n_samps = zeros(n_lfp_trials,1);
for i = 1:n_lfp_trials
    lfp_start_t(i) = LFP.Trials(i).Start/1e4;
    lfp_end_t(i) = LFP.Trials(i).End/1e4;
    lfp_n_samps(i) = size(LFP.Trials(i).LFP,1);
end
lfp_durs = lfp_end_t - lfp_start_t;

cd ~/Data/bruce/ExptB/Raw/
load(sprintf('lemM237.%d.em.mat',expno));

n_eye_trials = length(Expt.Trials);
eye_start_t = zeros(n_eye_trials,1);
eye_end_t = zeros(n_eye_trials,1);
for i = 1:n_eye_trials
eye_start_t(i) = Expt.Trials(i).Start/1e4; % time of first eye sample
eye_end_t(i) = Expt.Trials(i).End/1e4; % time of last eye sample
end
    
disp('Block times')
Blocks{blockid}.blocktimes'

disp('Stim times')
Blocks{blockid}.stimtime([1 end])'
% 
disp('LFP times')
[lfp_start_t lfp_end_t lfp_durs lfp_n_samps/1e3]

disp('Eye times')
[eye_start_t eye_end_t]