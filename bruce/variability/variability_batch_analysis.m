%% MODEL BASED VARIABILITY ANALYSIS

clear all
addpath('~/James_scripts/bruce/eye_tracking_improvements/');
global Expt_name bar_ori use_MUA
use_MUA = false;

Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296'};
bar_ori = nan;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    model_based_variability_v4;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori use_LOOXV
end
Expt_name = 'M297'; bar_ori = 0;
model_based_variability_v4;
Expt_name = 'M297'; bar_ori = 90;
model_based_variability_v4;


Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
bar_ori = 0;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    model_based_variability_v4;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori use_LOOXV
end

Expt_list = {'G085','G086','G087','G088','G089','G091','G093'};
bar_ori = 90;
for Elist_cnt = 1:length(Expt_list)
    Expt_name = Expt_list{Elist_cnt};
    model_based_variability_v4;
    clearvars -except Elist_cnt Expt_list Expt_name bar_ori use_LOOXV
end