clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');

% Expt_list = {'M266','M270'};
Expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
global Expt_name bar_ori use_LOOXV
bar_ori = 90;
use_LOOXV = 1;
for ee = 6:length(Expt_list)
    Expt_name = Expt_list{ee};
    full_eyetracking_withblinks;
    clearvars -except ee Expt_list Expt_name bar_ori use_LOOXV
end
