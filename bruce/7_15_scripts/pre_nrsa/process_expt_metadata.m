clear all
close all
clc
cd /Users/James/Data/bruce/7_15_12/

% cd G029/
% load ./G029Expts.mat
% cd G030/
% load ./G030Expts.mat
% cd G034/
% load ./G034Expts.mat
cd G035/
load ./G035Expts.mat

for i = 1:length(Expts)
    cur_dur(i) = (Expts{i}.Header.End - Expts{i}.Header.Start)/1e4;
    num_trials(i) = length(Expts{i}.Trials);
    stimvals_fs = Expts{i}.Stimvals.Fr;
    Trial_durs = [Expts{i}.Trials(:).dur]/1e4;
    total_trial_dur(i) = sum(Trial_durs);
    fprintf('Expt: %d. %s sv: %d\n\n',i,Expts{i}.Header.expname,stimvals_fs);
end